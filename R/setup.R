# -----------------------------------------------------------------------------
# File:    R/setup.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. All figures and tables.
# Purpose: Loader sourced first by every script: finds the repository root,
#          sources the library files in R/, and defines paths, the settings
#          shared by all scripts and the command line parser.
# Inputs:  None.
# Outputs: None (creates results/ and figures/ when a path helper is used).
# Usage:   source("R/setup.R")   # from a script started in the repository root
# Source:  R_initial.R, RPC_init.R (settings of the original runs).
# -----------------------------------------------------------------------------

REPO_ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile), ".."))

for (file in c("mantel.R", "ridge_cv.R", "adamant_cv.R", "simulation.R",
               "power_curve.R", "data.R", "plots.R")) {
  source(file.path(REPO_ROOT, "R", file))
}

# ----- Settings shared by all original runs -----
SEED        <- 8891   # seed of designs, permutations and replications
ALPHA       <- 0.05   # significance level
SIGMA       <- 0.5    # sd of the noise added to the response
UPPER       <- 1000   # largest ridge penalty of every grid
N_PERMS     <- 1000   # permutations per test
FIGURE_SIZE <- 7      # width and height of every figure, inches

# Simulation design of Tables 1 to 3 and Figure 3
SIM_N         <- 350    # sample size
SIM_RHO       <- 0.025  # correlation of the compound symmetric covariance
SIM_BETA_MEAN <- 1      # value (fixed) or mean (random) of the coefficients

# ----- Paths -----

# Raw data: the folder named by ADAMANT_DATA, or data/ in the repository.
data_path <- function(...) {
  file.path(Sys.getenv("ADAMANT_DATA", file.path(REPO_ROOT, "data")), ...)
}

reference_path <- function(...) file.path(REPO_ROOT, "reference", ...)
archive_path   <- function(...) file.path(REPO_ROOT, "archive", ...)

# Outputs go to results/ and figures/ in the folder named by ADAMANT_OUTPUT,
# or in the repository.
output_path <- function(folder, name) {
  folder <- file.path(Sys.getenv("ADAMANT_OUTPUT", REPO_ROOT), folder)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  file.path(folder, name)
}

results_path <- function(name) output_path("results", name)
figures_path <- function(name) output_path("figures", name)

# ----- Command line -----

# Read `--name=value` arguments and bare `--name` flags. `defaults` is a named
# list; each value given on the command line is converted to the type of its
# default, and a flag sets a logical default to TRUE. Every script accepts
# --rerun: recompute from the data instead of reading the archived runs.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  defaults <- c(list(rerun = FALSE), defaults)
  pairs <- regmatches(args, regexec("^--([a-z_]+)(=(.*))?$", args))
  if (any(lengths(pairs) != 4)) stop("arguments must look like --name=value or --name")
  for (pair in pairs) {
    name <- pair[2]
    if (!name %in% names(defaults)) stop("unknown argument --", name)
    is_flag <- pair[3] == ""
    if (is_flag != is.logical(defaults[[name]])) stop("--", name, if (is_flag) " needs a value" else " takes no value")
    defaults[[name]] <- if (is_flag) TRUE else methods::as(pair[4], class(defaults[[name]]))
  }
  defaults
}

# Turn "3:7" into 3:7 and "2,5,9" into c(2, 5, 9). A single number N is that
# one index, or 1:N when `count = TRUE`.
parse_index_set <- function(text, count = FALSE) {
  if (grepl("^[0-9]+:[0-9]+$", text)) {
    ends <- as.integer(strsplit(text, ":")[[1]])
    return(ends[1]:ends[2])
  }
  index <- as.integer(strsplit(text, ",")[[1]])
  if (count && length(index) == 1) seq_len(index) else index
}

# Label of the settings a rerun was computed with, used in the names of its
# stored parts so that runs with different settings are kept apart.
settings_tag <- function(opts, n_runs = NULL) {
  runs <- if (is.null(n_runs)) "" else sprintf("runs%d_", n_runs)
  sprintf("%sperms%d_grid%dx%d", runs, opts$perms, opts$negative, opts$positive)
}

# Read a result from `file` if it exists; otherwise compute it and store it.
cached_csv <- function(file, compute) {
  if (file.exists(file)) {
    log_message("reusing %s", file)
    return(read.csv(file))
  }
  result <- compute()
  write.csv(result, file, row.names = FALSE)
  result
}

log_message <- function(msg, ...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", sprintf(msg, ...), "\n")
}
