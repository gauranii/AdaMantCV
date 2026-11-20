# -----------------------------------------------------------------------------
# File:    R/plots.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figures 2, 3 and 5 to 8.
# Purpose: Plots shared by several figures: the heat map of a curve over the
#          ridge penalties, the scatter of a curve against the penalty, and
#          the histogram of chosen penalties.
# Inputs:  None.
# Outputs: figures/<name>.pdf (written by save_figure()).
# Usage:   source("R/setup.R"); heat_map(curves, "power", 0.4, 3, 500, 10)
# Source:  R_Heat1.R, R_Heat2.R, R_Heat3.R, R_Output_plots.R, R_Lambda.R,
#          RPC_fcns.R.
# -----------------------------------------------------------------------------

library(ggplot2)

#' Save a square figure as pdf.
#'
#' @param plot A ggplot object.
#' @param name File name.
#' @param figures_dir Folder of the figures.
#' @param size Width and height in inches.
#' @return NULL, invisibly.
save_figure <- function(plot, name, figures_dir, size) {
  ggsave(output_file(figures_dir, name), plot, width = size, height = size)
  log_message("wrote figures/%s", name)
}

#' Heat map of one column of several curve tables against the penalties.
#'
#' @param curves Named list of curve tables; one row of tiles per table, the
#'   first table at the bottom.
#' @param column Column of the tables to show.
#' @param midpoint White point of the colour scale.
#' @param digits Decimals to which the penalties are rounded on the axis.
#' @param zero_position Position on the penalty axis of the green line that
#'   marks zero.
#' @param text_size Size of the axis text.
#' @return A ggplot object.
heat_map <- function(curves, column, midpoint, digits, zero_position, text_size) {
  long <- do.call(rbind, lapply(names(curves), function(label) {
    data.frame(lambdas = round(curves[[label]]$lambdas, digits),
               variable = factor(label, levels = names(curves)),
               value = curves[[label]][[column]])
  }))
  ggplot(long, aes(x = as.factor(lambdas), y = variable, fill = value)) +
    geom_tile() +
    xlab("lambda") + ylab("") +
    scale_fill_gradient2(low = "navy", mid = "white", high = "red",
                         midpoint = midpoint, limits = range(long$value)) +
    geom_vline(xintercept = zero_position, linewidth = 0.5, colour = "green") +
    scale_x_discrete(breaks = c(-1, 0, 1000)) +
    theme(axis.text = element_text(size = text_size),
          axis.title = element_text(size = 8))
}

#' Save the two heat maps of a figure: mean -log10 p-value (panel A) and
#' power (panel B).
#'
#' The files are named <figure>A_<subject>_heatmap_mean-log10-pvalue<tail>.pdf
#' and <figure>B_<subject>_heatmap_power<tail>.pdf.
#'
#' @param figure Figure of the paper (for example "Fig5").
#' @param subject Data shown ("liver" or "EEG").
#' @param tail End of the file names ("" or for example "_top-triples").
#' @param curves Named list of curve tables.
#' @param power_column Column shown in panel B.
#' @param midpoints White points of the colour scales, `c(logp = , power = )`.
#' @param digits,zero_position,text_size As for `heat_map()`.
#' @param figures_dir,size As for `save_figure()`.
#' @return NULL, invisibly.
save_heat_maps <- function(figure, subject, tail, curves, power_column, midpoints, digits,
                           zero_position, text_size, figures_dir, size) {
  save_figure(heat_map(curves, "ave_logpval", midpoints[["logp"]], digits, zero_position,
                       text_size),
              sprintf("%sA_%s_heatmap_mean-log10-pvalue%s.pdf", figure, subject, tail), figures_dir, size)
  save_figure(heat_map(curves, power_column, midpoints[["power"]], digits, zero_position,
                       text_size),
              sprintf("%sB_%s_heatmap_power%s.pdf", figure, subject, tail), figures_dir, size)
}

#' Points of one column of a curve table against the penalty, coloured by the
#' sign of the penalty.
#'
#' @param table A curve table.
#' @param column Column to show.
#' @param label Label of the vertical axis.
#' @return A ggplot object.
sign_scatter <- function(table, column, label) {
  ggplot(table, aes(x = lambdas, y = .data[[column]])) +
    geom_point(aes(color = 2 + sign(lambdas)), size = 1) + theme_bw() +
    xlab("lambda") + ylab(label) +
    theme(legend.position = "none")
}

#' Save the power panel (all penalties) and the GCV panel (a window of
#' penalties) of several triples.
#'
#' The power panels take the first panel letters and the GCV panels the next:
#' <figure><letter>_EEG_power-<trend>_<triple name>.pdf and
#' <figure><letter>_EEG_GCV_<triple name>.pdf.
#'
#' @param figure Figure of the paper (for example "Fig7").
#' @param trend Trend of the power in the file names ("decreasing" or
#'   "increasing").
#' @param windows Data frame with columns `triple`, `first`, `last`, `y_max`:
#'   the GCV panel shows rows `first` to `last` of the curve table, without
#'   the rows with gcv above `y_max`.
#' @param curves List of curve tables, one per row of `windows`, named
#'   "<channel1>-<channel2>-<band>" as by `eeg_curves()`.
#' @param figures_dir,size As for `save_figure()`.
#' @return NULL, invisibly.
save_power_gcv_panels <- function(figure, trend, windows, curves, figures_dir, size) {
  n <- nrow(windows)
  for (i in seq_len(n)) {
    w <- windows[i, ]
    table <- curves[[i]]
    save_figure(sign_scatter(table, "power", "power"),
                sprintf("%s%s_EEG_power-%s_%s.pdf", figure, LETTERS[i], trend, names(curves)[i]), figures_dir, size)
    shown <- table[intersect(w$first:w$last, seq_len(nrow(table))), ]
    save_figure(sign_scatter(shown[shown$gcv <= w$y_max, ], "gcv",
                             "Generalized Cross Validation"),
                sprintf("%s%s_EEG_GCV_%s.pdf", figure, LETTERS[n + i], names(curves)[i]), figures_dir, size)
  }
}

#' Histograms of the penalties chosen by GCV and by LOOCV.
#'
#' @param chosen Data frame with columns `GCV` and `LOOCV`.
#' @param tau Label of the truncation in the legend.
#' @param bin_width Width of the bins.
#' @return A ggplot object.
penalty_histogram <- function(chosen, tau, bin_width) {
  long <- data.frame(variable = rep(paste(c("GCV", "LOOCV"), "tau =", tau),
                                    each = nrow(chosen)),
                     value = c(chosen$GCV, chosen$LOOCV))
  is_gcv <- startsWith(long$variable, "GCV")
  ggplot(long, aes(x = value, color = variable)) +
    geom_histogram(data = long[is_gcv, ], fill = "red", binwidth = bin_width,
                   alpha = 0.2) +
    geom_histogram(data = long[!is_gcv, ], fill = "white", binwidth = bin_width,
                   alpha = 0.2) +
    xlab("optimal ridge penalty estimates") + ylab("") + theme_bw() +
    theme(legend.position = "bottom")
}
