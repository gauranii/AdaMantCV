# -----------------------------------------------------------------------------
# File:    R/setup.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. All figures and tables.
# Purpose: Loader sourced first by every script: finds the repository root,
#          sources the library files in R/, and defines the settings of the
#          paper, the project paths and the command line parser.
# Inputs:  None.
# Outputs: None.
# Usage:   source("R/setup.R")   # from a script started in the repository root
# Source:  R_initial.R, RPC_init.R (settings of the original runs).
# -----------------------------------------------------------------------------

REPO_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))

for (file in c("mantel.R", "ridge_cv.R", "adamant_cv.R", "simulation.R", "power_curve.R",
               "data.R", "eeg_tests.R", "liver_cv.R", "plots.R")) {
  source(file.path(REPO_ROOT, "R", file))
}

# ----- Settings of the paper, passed to the functions as `settings` -----
PAPER <- list(
  seed        = 8891,  # seed of designs, permutations and replications
  alpha       = 0.05,  # significance level
  sigma       = 0.5,   # sd of the noise added to the response
  upper       = 1000,  # largest ridge penalty of every grid
  n_perms     = 1000,  # permutations per test
  figure_size = 7      # width and height of every figure, inches
)

# ----- Paths -----

#' Folders and files of the project.
#'
#' @param root Repository root.
#' @param data Folder of the data; its subfolder 04_coherence_SNPs/ holds the
#'   EEG files.
#' @param output Folder in which results/ and figures/ are written.
#' @return Named list of paths: `eeg`, `liver`, `reference`, `archive`,
#'   `results`, `figures`.
project_paths <- function(root,
                          data = Sys.getenv("ADAMANT_DATA", file.path(root, "data")),
                          output = Sys.getenv("ADAMANT_OUTPUT", root)) {
  list(eeg       = file.path(data, "04_coherence_SNPs"),
       liver     = file.path(root, "data", "liver", "mixOmicsData.csv"),
       reference = file.path(root, "reference"),
       archive   = file.path(root, "archive"),
       results   = file.path(output, "results"),
       figures   = file.path(output, "figures"))
}

PATHS <- project_paths(REPO_ROOT)

#' Path of an output file, creating its folder when needed.
#'
#' @param folder Output folder (for example `paths$results`).
#' @param name File name.
#' @return The path `folder/name`.
output_file <- function(folder, name) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  file.path(folder, name)
}

# ----- Command line -----

#' Read `--name=value` arguments and bare `--name` flags.
#'
#' Every script accepts `--rerun`: recompute from the data instead of reading
#' the archived runs.
#'
#' @param defaults Named list of options and their default values. A value
#'   given on the command line is converted to the type of its default; a
#'   flag sets a logical default to TRUE.
#' @param args Command line arguments.
#' @return `defaults` with the values given on the command line, and `rerun`.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  defaults <- c(list(rerun = FALSE), defaults)
  pairs <- regmatches(args, regexec("^--([a-z_]+)(=(.*))?$", args))
  if (any(lengths(pairs) != 4)) stop("arguments must look like --name=value or --name")
  for (pair in pairs) {
    name <- pair[2]
    if (!name %in% names(defaults)) stop("unknown argument --", name)
    is_flag <- pair[3] == ""
    if (is_flag != is.logical(defaults[[name]])) {
      stop("--", name, if (is_flag) " needs a value" else " takes no value")
    }
    defaults[[name]] <- if (is_flag) TRUE else methods::as(pair[4], class(defaults[[name]]))
  }
  defaults
}

#' Default options of a script that computes curves over a penalty grid.
#'
#' @param grid List with the grid sizes `negative` and `positive`.
#' @param settings Settings of the paper (`n_perms`).
#' @return Named list of defaults for `parse_args()`.
curve_options <- function(grid, settings) {
  list(runs = 1000L, perms = settings$n_perms, negative = grid$negative,
       positive = grid$positive, cores = 4)
}

#' Default options of a script that runs AdaMantCV tests.
#'
#' @param settings Settings of the paper (`n_perms`).
#' @param cores Default number of parallel processes.
#' @return Named list of defaults for `parse_args()`.
test_options <- function(settings, cores = 4) {
  list(perms = settings$n_perms, negative = 50, positive = 200, cores = cores)
}

#' Turn an index text into indices.
#'
#' @param text "3:7" (3 to 7), "2,5,9", or a single number N.
#' @param count If TRUE a single number N means 1 to N.
#' @return Integer vector.
parse_index_set <- function(text, count = FALSE) {
  if (grepl("^[0-9]+:[0-9]+$", text)) {
    ends <- as.integer(strsplit(text, ":")[[1]])
    return(ends[1]:ends[2])
  }
  index <- as.integer(strsplit(text, ",")[[1]])
  if (count && length(index) == 1) seq_len(index) else index
}

#' Label of the settings a rerun was computed with.
#'
#' Used in the names of the stored parts of a rerun, so that runs with
#' different settings are kept apart.
#'
#' @param opts Options with `perms`, `negative`, `positive`.
#' @param n_runs Number of replications, or NULL to leave it out.
#' @return Character label.
settings_tag <- function(opts, n_runs = NULL) {
  runs <- if (is.null(n_runs)) "" else sprintf("runs%d_", n_runs)
  sprintf("%sperms%d_grid%dx%d", runs, opts$perms, opts$negative, opts$positive)
}

#' Read a result from a csv file if it exists; otherwise compute and store it.
#'
#' @param file Path of the csv file.
#' @param compute Function without arguments that returns a data frame.
#' @return The data frame.
cached_csv <- function(file, compute) {
  if (file.exists(file)) {
    log_message("reusing %s", file)
    return(read.csv(file))
  }
  result <- compute()
  write.csv(result, file, row.names = FALSE)
  result
}

#' Print a time-stamped message.
#'
#' @param msg Format of the message, as for `sprintf()`.
#' @param ... Values for the format.
#' @return NULL, invisibly.
log_message <- function(msg, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", sprintf(msg, ...), "\n")
}
