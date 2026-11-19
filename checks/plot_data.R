# -----------------------------------------------------------------------------
# File:    checks/plot_data.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Run one script of scripts/ unchanged and, beside each figure it
#          saves, write the numbers the figure plots (the layer data of the
#          built ggplot) as csv, for the comparison with the Python figures.
# Inputs:  those of the script.
# Outputs: those of the script, and <plot data dir>/<figure name>.csv with
#          one row per plotted mark: layer, x, y, xmin, xmax, count,
#          xintercept, group, colour, value (those the layer has).
# Usage:   ADAMANT_PLOT_DATA=<plot data dir> Rscript checks/plot_data.R
#          <script name without .R> [options of the script]
# -----------------------------------------------------------------------------

script_args <- commandArgs(trailingOnly = TRUE)
source("R/setup.R")

PLOT_DATA_COLUMNS <- c("x", "y", "xmin", "xmax", "count", "xintercept", "group", "colour")

#' The numbers a plot draws, one row per mark of every layer.
#' @param plot A ggplot object.
#' @return Data frame with the column `layer` and those of PLOT_DATA_COLUMNS
#'   the layers have; `value` is the fill value of a heat map.
plotted_data <- function(plot) {
  layers <- ggplot_build(plot)$data
  rows <- lapply(seq_along(layers), function(i) {
    layer <- layers[[i]]
    out <- data.frame(layer = rep(i, nrow(layer)))
    for (column in PLOT_DATA_COLUMNS) {
      out[[column]] <- if (column %in% names(layer)) as.vector(unclass(layer[[column]])) else NA
    }
    is_tiles <- i == 1 && "fill" %in% names(layer) && "value" %in% names(plot$data) &&
      nrow(layer) == nrow(plot$data) && !"count" %in% names(layer)
    out$value <- if (is_tiles) plot$data$value else NA
    out
  })
  do.call(rbind, rows)
}

# The scripts call save_figure(); this one also writes the plotted numbers.
save_figure_only <- save_figure
save_figure <- function(plot, name, figures_dir, size) {
  save_figure_only(plot, name, figures_dir, size)
  write.csv(plotted_data(plot),
            output_file(Sys.getenv("ADAMANT_PLOT_DATA"), sub("[.]pdf$", ".csv", name)),
            row.names = FALSE)
}

# The script sources R/setup.R (done above) and reads its own options.
source <- function(...) invisible(NULL)
commandArgs <- function(...) script_args[-1]
sys.source(file.path("scripts", paste0(script_args[1], ".R")), envir = globalenv())
