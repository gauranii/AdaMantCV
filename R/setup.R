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
data_path <- function(...) {
  file.path(Sys.getenv("ADAMANT_DATA", file.path(REPO_ROOT, "data")), ...)
}

reference_path <- function(...) file.path(REPO_ROOT, "reference", ...)

output_path <- function(folder, name) {
  dir.create(file.path(REPO_ROOT, folder), showWarnings = FALSE)
  file.path(REPO_ROOT, folder, name)
}

results_path <- function(name) output_path("results", name)
# Figures of the "paper" method get the suffix _paper (see parse_args()).
figures_path <- function(name) {
  if (METHOD == "paper") name <- sub("[.]pdf$", "_paper.pdf", name)
  output_path("figures", name)
}

# ----- Command line -----

# Method of every computation (see R/adamant_cv.R): "original" (default, the
# code that made the paper) or "paper" (the method as the paper states it).
# Scripts set it with --method=paper; parse_args() stores it here.
METHOD <- "original"

# Read `--name=value` arguments. `defaults` is a named list; each value given
# on the command line is converted to the type of its default. A `method`
# entry also sets METHOD.
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  pairs <- regmatches(args, regexec("^--([a-z_]+)=(.*)$", args))
  if (any(lengths(pairs) != 3)) stop("arguments must look like --name=value")
  for (pair in pairs) {
    name <- pair[2]
    if (!name %in% names(defaults)) stop("unknown argument --", name)
    defaults[[name]] <- methods::as(pair[3], class(defaults[[name]]))
  }
  if (!is.null(defaults$method)) {
    stopifnot(defaults$method %in% c("original", "paper"))
    METHOD <<- defaults$method
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

# Label of the settings a result file was computed with, used in file names so
# that reduced test runs never overwrite or pass for full runs.
settings_tag <- function(opts, n_runs = NULL) {
  runs <- if (is.null(n_runs)) "" else sprintf("runs%d_", n_runs)
  method <- if (METHOD == "paper") "paper_" else ""
  sprintf("%s%sperms%d_grid%dx%d", method, runs, opts$perms, opts$negative, opts$positive)
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
