# -----------------------------------------------------------------------------
# File:    checks/check_equivalence.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Compare the outputs of two versions of the code, file by file:
#          every csv under results/ (tables, curves, per-replication p-values)
#          and every figure (pdf, without its creation dates).
# Inputs:  two output folders (each with results/ and figures/), written by
#          the same commands run with ADAMANT_OUTPUT set to each folder, once
#          from a checkout of the reference commit and once from this one.
# Outputs: prints one line per file and a summary; exit status 1 unless every
#          file is equal: text columns and columns of p-values exactly, other
#          numbers within `tolerance`.
# Usage:   Rscript checks/check_equivalence.R <reference output> <output>
#          [--tolerance=1e-12]
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
folders <- args[!startsWith(args, "--")]
tolerance <- sub("--tolerance=", "", args[startsWith(args, "--tolerance=")])
tolerance <- if (length(tolerance)) as.numeric(tolerance) else 1e-12
stopifnot(length(folders) == 2, dir.exists(folders))

#' Files under a folder, relative to it.
#' @param folder Path of an output folder.
#' @return Sorted character vector of relative paths of csv and pdf files.
output_files <- function(folder) {
  sort(list.files(folder, pattern = "[.](csv|pdf)$", recursive = TRUE))
}

#' Largest difference between two csv files.
#' @param reference,current Paths of the two files.
#' @return List with `exact` (TRUE when the parsed tables are identical) and
#'   `difference` (largest absolute difference of a numeric cell; Inf when
#'   shapes, text or p-value columns are not identical).
compare_csv <- function(reference, current) {
  a <- read.csv(reference)
  b <- read.csv(current)
  if (identical(a, b)) return(list(exact = TRUE, difference = 0))
  if (!identical(dim(a), dim(b)) || !identical(names(a), names(b))) {
    return(list(exact = FALSE, difference = Inf))
  }
  exact_column <- !vapply(a, is.numeric, logical(1)) |
    grepl("^(p_|gcv$|loo$|loocv$|ga_|power|b$|run$|scenario$|triple$)", names(a))
  if (!identical(a[exact_column], b[exact_column])) return(list(exact = FALSE, difference = Inf))
  list(exact = FALSE, difference = max(abs(as.matrix(a[!exact_column]) - as.matrix(b[!exact_column]))))
}

#' Whether two pdf files are equal apart from their date entries.
#' @param reference,current Paths of the two files.
#' @return List with `exact` and `difference` (0 or Inf).
compare_pdf <- function(reference, current) {
  body <- function(file) {
    bytes <- readBin(file, "raw", file.size(file))
    for (key in c("/CreationDate", "/ModDate")) {
      for (at in grepRaw(key, bytes, fixed = TRUE, all = TRUE)) {
        bytes[at + nchar(key) + 0:18] <- as.raw(0)  # the date, " (D:yyyymmddhhmmss"
      }
    }
    bytes
  }
  same <- identical(body(reference), body(current))
  list(exact = same, difference = if (same) 0 else Inf)
}

files <- output_files(folders[1])
if (!identical(files, output_files(folders[2]))) stop("the two folders hold different files")
ok <- vapply(files, function(file) {
  compare <- if (endsWith(file, ".pdf")) compare_pdf else compare_csv
  result <- compare(file.path(folders[1], file), file.path(folders[2], file))
  pass <- result$difference <= tolerance
  cat(sprintf("%-70s %s\n", file,
              if (result$exact) "identical" else if (pass) sprintf("within %.1e", result$difference) else "NOT EQUAL"))
  pass
}, logical(1))
cat(sprintf("equivalence: %d of %d files equal\n", sum(ok), length(ok)))
quit(status = if (all(ok)) 0 else 1)
