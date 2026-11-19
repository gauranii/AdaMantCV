# -----------------------------------------------------------------------------
# File:    python/adamant/plots.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Figures 2, 3 and 5 to 8.
# Purpose: Plots shared by several figures: the heat map of a curve over the
#          ridge penalties, the scatter of a curve against the penalty, and
#          the histogram of chosen penalties, drawn with matplotlib in the
#          layout of the ggplot2 figures of R/plots.R.
# Inputs:  None.
# Outputs: figures/<name>.pdf (written by save_figure()); when the
#          environment variable ADAMANT_PLOT_DATA names a folder, also
#          <folder>/<name>.csv with the numbers the figure plots.
# Usage:   from adamant.plots import heat_map
#          heat_map(curves, "power", 0.4, 3, 500, 10)
# Source:  R/plots.R.
# -----------------------------------------------------------------------------
import os

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import ListedColormap, to_rgb  # noqa: E402

from .setup import log_message, output_file, write_csv  # noqa: E402

# Columns of the csv of plotted numbers, one row per mark.
PLOT_DATA_COLUMNS = ["layer", "x", "y", "xmin", "xmax", "count", "xintercept", "group",
                     "colour", "value"]
# ggplot2's default colours of a continuous scale at its ends and middle.
GRADIENT = {0.0: "#132B43", 0.5: "#336A98", 1.0: "#56B1F7"}
PANEL_GREY = "#EBEBEB"


class Plot:
    """A figure with one panel and the numbers it plots.

    Parameters: theme ("grey": grey panel with white grid lines, or "bw":
    white panel with a black frame and light grid lines).
    """

    def __init__(self, theme="grey"):
        self.figure, self.ax = plt.subplots(layout="constrained")
        self.layers = []
        ax = self.ax
        ax.set_axisbelow(True)
        if theme == "grey":
            ax.set_facecolor(PANEL_GREY)
            ax.grid(color="white", linewidth=0.8)
            for side in ax.spines.values():
                side.set_visible(False)
        else:
            ax.grid(color="#EBEBEB", linewidth=0.6)
            for side in ax.spines.values():
                side.set_color("#333333")
        ax.tick_params(colors="#4D4D4D", labelsize=8.8)

    def record(self, layer, **columns):
        """Keep the numbers of the marks of one layer.

        Parameters: layer (layer number), columns (values per mark, named
        as PLOT_DATA_COLUMNS). Returns: None.
        """
        self.layers.append(pd.DataFrame({"layer": layer, **columns}))

    def plotted_data(self):
        """The numbers the figure plots.

        Returns: DataFrame with the columns PLOT_DATA_COLUMNS.
        """
        return pd.concat(self.layers, ignore_index=True).reindex(columns=PLOT_DATA_COLUMNS)


def save_figure(plot, name, figures_dir, size):
    """Save a square figure as pdf.

    Parameters: plot (a `Plot`), name (file name), figures_dir (folder of
    the figures), size (width and height in inches). Returns: None.
    """
    plot.figure.set_size_inches(size, size)
    plot.figure.savefig(output_file(figures_dir, name), bbox_inches=None)
    plt.close(plot.figure)
    log_message("wrote figures/%s", name)
    data_dir = os.environ.get("ADAMANT_PLOT_DATA")
    if data_dir:
        write_csv(plot.plotted_data(), output_file(data_dir, name.replace(".pdf", ".csv")))


def hue_colours(n):
    """ggplot2's default colours of n groups: hues 15, 15 + 360/n, ... at
    chroma 100 and luminance 65.

    Parameters: n (number of groups). Returns: list of hex colours.
    """
    return [_hcl(15 + 360 * i / n, 100, 65) for i in range(n)]


def _hcl(h, c, l):
    """Hex colour of polar LUV coordinates (as R's `hcl()`)."""
    white_y, white_u, white_v = 100.0, 0.1978398, 0.4683363
    u, v = c * np.cos(np.radians(h)), c * np.sin(np.radians(h))
    y = white_y * (((l + 16) / 116) ** 3 if l > 7.999592 else l / 903.3)
    u, v = u / (13 * l) + white_u, v / (13 * l) + white_v
    x, z = 9.0 * y * u / (4 * v), -(3 * u * y / (4 * v)) - 5 * y + 3 * y / v
    linear = np.array([3.240479 * x - 1.537150 * y - 0.498535 * z,
                       -0.969256 * x + 1.875992 * y + 0.041556 * z,
                       0.055648 * x - 0.204043 * y + 1.057311 * z]) / white_y
    srgb = np.where(linear > 0.0031308, 1.055 * np.maximum(linear, 0) ** (1 / 2.4) - 0.055,
                    12.92 * linear)
    return "#%02X%02X%02X" % tuple(int(round(255 * min(max(value, 0), 1))) for value in srgb)


def _lab(colour):
    """CIE Lab coordinates of a colour."""
    rgb = np.array(to_rgb(colour))
    linear = np.where(rgb > 0.04045, ((rgb + 0.055) / 1.055) ** 2.4, rgb / 12.92)
    xyz = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750],
                    [0.0193339, 0.1191920, 0.9503041]]) @ linear / [0.95047, 1.0, 1.08883]
    f = np.where(xyz > 216 / 24389, np.cbrt(xyz), (24389 / 27 * xyz + 16) / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def _from_lab(lab):
    """sRGB colour of CIE Lab coordinates."""
    fy = (lab[0] + 16) / 116
    f = np.array([fy + lab[1] / 500, fy, fy - lab[2] / 200])
    xyz = np.where(f ** 3 > 216 / 24389, f ** 3, (116 * f - 16) * 27 / 24389) \
        * [0.95047, 1.0, 1.08883]
    linear = np.array([[3.2404542, -1.5371385, -0.4985314], [-0.9692660, 1.8760108, 0.0415560],
                       [0.0556434, -0.2040259, 1.0572252]]) @ xyz
    linear = np.clip(linear, 0, 1)
    return np.where(linear > 0.0031308, 1.055 * linear ** (1 / 2.4) - 0.055, 12.92 * linear)


def diverging_colours(low="navy", mid="white", high="red", steps=511):
    """Colour map from `low` through `mid` to `high`, interpolated in Lab
    space as `scale_fill_gradient2()` does.

    Parameters: low, mid, high (colours), steps (number of colours).
    Returns: a matplotlib colour map over [0, 1], `mid` at 0.5.
    """
    anchors = [_lab(low), _lab(mid), _lab(high)]
    colours = []
    for t in np.linspace(0, 1, steps):
        a, b, w = (anchors[0], anchors[1], 2 * t) if t <= 0.5 else \
            (anchors[1], anchors[2], 2 * t - 1)
        colours.append(np.clip(_from_lab(a + w * (b - a)), 0, 1))
    return ListedColormap(colours)


def r_round(x, digits):
    """`round(x, digits)` of a vector. Returns: float array."""
    return np.array([round(float(value), digits) for value in x])


def heat_map(curves, column, midpoint, digits, zero_position, text_size):
    """Heat map of one column of several curve tables against the penalties.

    Parameters: curves (dict of curve tables; one row of tiles per table,
    the first table at the bottom), column (column of the tables to show),
    midpoint (white point of the colour scale), digits (decimals to which
    the penalties are rounded on the axis), zero_position (position on the
    penalty axis of the green line that marks zero), text_size (size of the
    axis text).
    Returns: a figure (`Plot`).
    """
    labels = list(curves)
    lambdas = np.concatenate([r_round(curves[label]["lambdas"], digits) for label in labels])
    values = np.concatenate([curves[label][column].to_numpy() for label in labels])
    rows = np.concatenate([np.full(len(curves[label]), i + 1) for i, label in enumerate(labels)])
    levels = np.unique(lambdas)
    position = np.searchsorted(levels, lambdas) + 1

    # A tile per row of the tables; tiles at the same place lie over each
    # other, the later row on top.
    tiles = np.full((len(labels), len(levels)), np.nan)
    tiles[rows - 1, position - 1] = values
    low, high = np.nanmin(values), np.nanmax(values)
    extent = 2 * max(abs(low - midpoint), abs(high - midpoint))
    plot = Plot(theme="grey")
    ax = plot.ax
    image = ax.imshow((tiles - midpoint) / extent + 0.5, cmap=diverging_colours(), vmin=0, vmax=1,
                      aspect="auto", origin="lower", interpolation="nearest",
                      extent=(0.5, len(levels) + 0.5, 0.5, len(labels) + 0.5))
    ax.grid(False)
    ax.axvline(zero_position, color="#00FF00", linewidth=1.0)  # R's "green"
    ax.set_xlim(0.4, len(levels) + 0.6)
    ax.set_ylim(0.4, len(labels) + 0.6)
    breaks = [(i + 1, level) for i, level in enumerate(levels) if level in (-1, 0, 1000)]
    ax.set_xticks([at for at, _ in breaks], [f"{level:g}" for _, level in breaks])
    ax.set_yticks(range(1, len(labels) + 1), labels)
    ax.tick_params(labelsize=text_size)
    ax.set_xlabel("lambda", fontsize=8)
    ax.set_ylabel("")
    # The colour bar shows the values between the smallest and the largest.
    bar = plot.figure.colorbar(image, ax=ax, shrink=0.2, aspect=5, pad=0.02)
    bar.ax.set_ylim((low - midpoint) / extent + 0.5, (high - midpoint) / extent + 0.5)
    ticks = matplotlib.ticker.MaxNLocator(5).tick_values(low, high)
    ticks = [tick for tick in ticks if low <= tick <= high]
    bar.set_ticks([(tick - midpoint) / extent + 0.5 for tick in ticks],
                  labels=[f"{tick:g}" for tick in ticks])
    bar.ax.set_title("value", fontsize=9, loc="left")
    bar.outline.set_visible(False)
    plot.record(1, x=position, y=rows, xmin=position - 0.5, xmax=position + 0.5, value=values)
    plot.record(2, xintercept=[zero_position], group=-1, colour="green")
    return plot


def save_heat_maps(figure, suffix, curves, power_column, midpoints, digits, zero_position,
                   text_size, figures_dir, size):
    """Save the two heat maps of a figure: mean -log10 p-value (panel A) and
    power (panel B).

    Parameters: figure (prefix of the file names, for example "fig5"),
    suffix (end of the file names, for example "eeg_heat.pdf"), curves (dict
    of curve tables), power_column (column shown in panel B), midpoints
    (white points of the colour scales, dict with logp and power), digits,
    zero_position, text_size (as for `heat_map()`), figures_dir, size (as
    for `save_figure()`). Returns: None.
    """
    save_figure(heat_map(curves, "ave_logpval", midpoints["logp"], digits, zero_position,
                         text_size),
                f"{figure}a_logp_{suffix}", figures_dir, size)
    save_figure(heat_map(curves, power_column, midpoints["power"], digits, zero_position,
                         text_size),
                f"{figure}b_power_{suffix}", figures_dir, size)


def sign_scatter(table, column, label):
    """Points of one column of a curve table against the penalty, coloured
    by the sign of the penalty.

    Parameters: table (a curve table), column (column to show), label (label
    of the vertical axis).
    Returns: a figure (`Plot`).
    """
    x, y = table["lambdas"].to_numpy(), table[column].to_numpy()
    sign = 2 + np.sign(x)
    span = sign.max() - sign.min() if len(sign) else 0
    shade = (sign - sign.min()) / span if span > 0 else np.full(len(sign), 0.5)
    colours = [GRADIENT[float(value)] for value in shade]
    plot = Plot(theme="bw")
    plot.ax.scatter(x, y, c=colours, s=5, linewidths=0)
    plot.ax.set_xlabel("lambda")
    plot.ax.set_ylabel(label)
    plot.record(1, x=x, y=y, group=-1, colour=colours)
    return plot


def save_power_gcv_panels(figure, windows, curves, figures_dir, size):
    """Save the power panel (all penalties) and the GCV panel (a window of
    penalties) of several triples.

    Parameters: figure (prefix of the file names, for example "fig7"),
    windows (DataFrame with columns triple, first, last, y_max: the GCV
    panel shows rows first to last of the curve table, without the rows with
    gcv above y_max), curves (list of curve tables, one per row of
    `windows`), figures_dir, size (as for `save_figure()`). Returns: None.
    """
    for w, table in zip(windows.to_dict("records"), curves):
        save_figure(sign_scatter(table, "power", "power"),
                    f"{figure}_power_{w['triple']}.pdf", figures_dir, size)
        shown = table.iloc[w["first"] - 1:min(w["last"], len(table))]
        save_figure(sign_scatter(shown[shown["gcv"] <= w["y_max"]], "gcv",
                                 "Generalized Cross Validation"),
                    f"{figure}_gcv_{w['triple']}.pdf", figures_dir, size)


def histogram_bins(values, x_range, bin_width):
    """Counts of values in bins of a given width, as `geom_histogram()` bins
    them: bin edges at odd multiples of half the width, bins closed on the
    right.

    Parameters: values (the values of one layer), x_range (smallest and
    largest value of all layers), bin_width (width of the bins).
    Returns: DataFrame with columns x (centre), xmin, xmax, count.
    """
    boundary = bin_width / 2
    origin = boundary + np.floor((x_range[0] - boundary) / bin_width) * bin_width
    max_x = x_range[1] + (1 - 1e-08) * bin_width
    breaks = origin + np.arange(int((max_x - origin) / bin_width + 1e-10) + 1) * bin_width
    if len(breaks) == 1:
        breaks = np.array([breaks[0], breaks[0] + bin_width])
    fuzz = 1e-08 * np.median(np.diff(breaks))
    fuzzy = breaks + np.concatenate(([-fuzz], np.full(len(breaks) - 1, fuzz)))
    index = np.searchsorted(fuzzy, values, side="left") - 1
    index[np.asarray(values) == fuzzy[0]] = 0
    inside = (index >= 0) & (index < len(breaks) - 1)
    count = np.bincount(index[inside], minlength=len(breaks) - 1)
    x = (breaks[1:] + breaks[:-1]) / 2
    width = np.diff(breaks)
    return pd.DataFrame({"x": x, "xmin": x - width / 2, "xmax": x + width / 2, "count": count})


def penalty_histogram(chosen, tau, bin_width):
    """Histograms of the penalties chosen by GCV and by LOOCV.

    Parameters: chosen (DataFrame with columns GCV and LOOCV), tau (label of
    the truncation in the legend), bin_width (width of the bins).
    Returns: a figure (`Plot`).
    """
    x_range = (min(chosen["GCV"].min(), chosen["LOOCV"].min()),
               max(chosen["GCV"].max(), chosen["LOOCV"].max()))
    outline = hue_colours(2)
    plot = Plot(theme="bw")
    for layer, (name, fill) in enumerate([("GCV", "red"), ("LOOCV", "white")], start=1):
        bins = histogram_bins(chosen[name].to_numpy(), x_range, bin_width)
        plot.ax.bar(bins["x"], bins["count"], width=bin_width, facecolor=(*to_rgb(fill), 0.2),
                    edgecolor=outline[layer - 1], linewidth=1.0, label=f"{name} tau = {tau}")
        plot.record(layer, x=bins["x"], y=bins["count"].astype(float), xmin=bins["xmin"],
                    xmax=bins["xmax"], count=bins["count"].astype(float), group=1,
                    colour=outline[layer - 1])
    plot.ax.set_xlabel("optimal ridge penalty estimates")
    plot.ax.set_ylabel("")
    plot.figure.legend(loc="outside lower center", ncol=3, frameon=False, fontsize=8.8,
                       handles=[matplotlib.patches.Patch(alpha=0, label="variable"),
                                *plot.ax.containers])
    return plot
