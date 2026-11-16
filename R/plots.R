# -----------------------------------------------------------------------------
# File:    R/plots.R
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figures 2 and 5 to 8.
# Purpose: Plots shared by several figures: the heat map of a curve over the
#          ridge penalties and the scatter of a curve against the penalty.
# Inputs:  None.
# Outputs: figures/<name>.pdf (written by save_figure()).
# Usage:   source("R/setup.R"); heat_map(curves, "power", midpoint = 0.4)
# Source:  R_Heat1.R, R_Heat2.R, R_Heat3.R, R_Output_plots.R, RPC_fcns.R.
# -----------------------------------------------------------------------------

library(ggplot2)

save_figure <- function(plot, name) {
  ggsave(figures_path(name), plot, width = FIGURE_SIZE, height = FIGURE_SIZE)
  log_message("wrote figures/%s", name)
}

# Heat map of one column of several curve tables (one row per table, first
# table at the bottom) against the penalties rounded to `digits` decimals. A
# green line marks position `zero_position` on the penalty axis (the zero).
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

# Points of one column of a curve table against the penalty, coloured by the
# sign of the penalty.
sign_scatter <- function(table, column, label) {
  ggplot(table, aes(x = lambdas, y = .data[[column]])) +
    geom_point(aes(color = sign), size = 1) + theme_bw() +
    xlab("lambda") + ylab(label) +
    theme(legend.position = "none")
}

# Power (all penalties) and GCV (a window of penalties) panels of the triples
# in `windows` (columns triple, first, last, y_max): the GCV panel shows rows
# first to last of the curve table, without the rows with gcv above y_max.
save_power_gcv_panels <- function(figure, windows, opts) {
  for (i in seq_len(nrow(windows))) {
    w <- windows[i, ]
    table <- eeg_curves(w$triple, opts)[[1]]
    save_figure(sign_scatter(table, "power", "power"),
                sprintf("%s_power_%d.pdf", figure, w$triple))
    shown <- table[intersect(w$first:w$last, seq_len(nrow(table))), ]
    save_figure(sign_scatter(shown[shown$gcv <= w$y_max, ], "gcv",
                             "Generalized Cross Validation"),
                sprintf("%s_gcv_%d.pdf", figure, w$triple))
  }
}
