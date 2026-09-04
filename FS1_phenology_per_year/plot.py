"""
FS1 plot: per-year phenology curves and the distribution of per-year peak
release dates.

Supporting Information companion to F2. F2 shows the 32-year climatological
mean with its interannual 5th-95th percentile envelope; that envelope
describes the spread of success *rates* at a given release date, not the
spread of *peak dates*. This figure separates the two:

  (a) all 32 individual per-year curves, with the climatological mean
      overlaid, and each year's own peak marked
  (b) the distribution of per-year peak release dates, with the median and
      inter-quartile range, and the climatological curve's peak for contrast

The two panels share an x-axis so the scatter of per-year peaks in (b) lines
up with the curves in (a).

The distinction matters because max() is nonlinear: the climatological curve
peaks later and lower than the average per-year peak, because averaging
curves whose maxima fall on different dates broadens and flattens the
composite. Both estimators are printed to stdout.

Reads the F2 aggregation rather than re-aggregating, so the two figures
cannot disagree about the same quantities. Pass --aggregated to override.

Inputs:  ../F2_phenology_curve/data/aggregated.nc
Output:  phenology_per_year.png
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.ticker import MaxNLocator


# ---------------------------------------------------------------------------
# Style (see ../FIGURE_STYLE.md)
# ---------------------------------------------------------------------------
#
# Two panels of equal height, each at the 4:3 aspect ratio of the
# single-panel figures, so the figure is 6.5 x 9.75 in. Everything else is
# left at its rcParams default: marker sizes, legend font size, histogram
# edges and axis limits are not overridden. The one exception is the
# individual-year line width, noted below.

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "figure.figsize": (6.5, 9.75),        # 6.5 in wide, two 4:3 panels
})


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# X-axis ticks, identical to F2 so the two figures read the same way.
TICK_DAYS = [16, 47, 78, 106]
TICK_LABELS = ["Dec 1", "Jan 1", "Feb 1", "Mar 1"]

# Last release date is Mar 14 (season_day 119). Season_day 120 (Mar 15) is
# allocated in the grid but never populated.
LAST_SEASON_DAY = 119

# Histogram bin width, in days.
BIN_WIDTH = 7

# One colour per panel, so no colour serves two purposes in the same panel:
# C0 throughout (a), C1 throughout (b). Within a panel, elements are
# separated by alpha and line width rather than by hue.
COLOR_A = "C0"
COLOR_B = "C1"

# The only departure from the rcParams defaults: 32 curves at the default
# line width (1.5) fill panel (a) solid, and alpha alone does not separate
# them from the climatological mean. The mean keeps the default width, so
# the default is what the eye anchors on.
YEAR_CURVE_LINEWIDTH = 1.0
YEAR_CURVE_COLOR = "black"
YEAR_CURVE_ALPHA = 0.1

# Vertical marker lines: panel (a)'s mean-curve peak and panel (b)'s median.
# One linestyle for both, at the default line width, each in its panel's
# colour. Dashed, matching F2's peak marker.
PEAK_LINESTYLE = "--"

# Legends are opaque throughout the paper.
LEGEND_FRAMEALPHA = 1.0

# Y-axis grid, matching the other figures: default grey (#b0b0b0) at 0.5.
GRID_ALPHA = 0.5

# Panel (b) only. Alpha separates the IQR fill from the bars; rwidth insets
# the bars so the bin edges are visible.
SPAN_ALPHA = 0.2
BAR_ALPHA = 0.5
BAR_RWIDTH = 0.8


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def season_day_to_label(day: int) -> str:
    """Calendar label for a season day index (season_day 0 = Nov 15)."""
    if day <= 15:
        return f"Nov {15 + day}"
    if day <= 46:
        return f"Dec {day - 15}"
    if day <= 77:
        return f"Jan {day - 46}"
    if day <= 105:
        return f"Feb {day - 77}"
    return f"Mar {day - 105}"


def resolve_legend_color(rc_key: str, axes_attr: str):
    """
    Resolve matplotlib's "inherit" sentinel for legend colors.

    legend.facecolor and legend.edgecolor default to the literal string
    "inherit", meaning "use the axes face/edge color". That sentinel is
    legend-only, so it has to be resolved explicitly when styling a bbox
    patch to match the legend. Same helper as F2.
    """
    val = mpl.rcParams[rc_key]
    if val == "inherit":
        return mpl.rcParams[axes_attr]
    return val


def framed_label(ax, text: str) -> None:
    """Panel label in a frame styled to match the legend (as in F4/F5)."""
    bbox = dict(
        boxstyle="round",
        facecolor=resolve_legend_color("legend.facecolor", "axes.facecolor"),
        edgecolor=resolve_legend_color("legend.edgecolor", "axes.edgecolor"),
        linewidth=mpl.rcParams["axes.linewidth"],
    )
    ax.text(0.985, 0.97, text, transform=ax.transAxes,
            ha="right", va="top", bbox=bbox)


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def compute_stats(ds: xr.Dataset) -> dict:
    """Per-year curves, the climatological mean, and both peak estimators."""
    sr = ds["success_rate"].values * 100.0        # (year, season_day), percent
    has_data = ds["has_data"].values
    season_day = ds["season_day"].values
    years = ds["year"].values

    sr_masked = np.where(has_data, sr, np.nan)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean_per_day = np.nanmean(sr_masked, axis=0)

    # Climatological peak: the max of the mean curve (what F2 annotates).
    clim_peak_day = int(season_day[np.nanargmax(mean_per_day)])
    clim_peak_value = float(np.nanmax(mean_per_day))

    # Per-year peaks: the mean of the maxima. Rows that are entirely NaN are
    # excluded rather than passed to nanargmax, which would raise.
    n_years = sr_masked.shape[0]
    peak_day = np.full(n_years, -1, dtype=int)
    peak_value = np.full(n_years, np.nan)
    for i in range(n_years):
        if not np.all(np.isnan(sr_masked[i])):
            peak_day[i] = int(np.nanargmax(sr_masked[i]))
            peak_value[i] = float(np.nanmax(sr_masked[i]))

    valid = peak_day >= 0
    pd_valid = peak_day[valid]
    p25, p50, p75 = np.percentile(pd_valid, [25, 50, 75])

    return {
        "season_day": season_day,
        "years": years,
        "sr": sr_masked,
        "mean_per_day": mean_per_day,
        "clim_peak_day": clim_peak_day,
        "clim_peak_value": clim_peak_value,
        "peak_day": peak_day,
        "peak_value": peak_value,
        "valid": valid,
        "peak_day_median": int(round(p50)),
        "peak_day_p25": int(round(p25)),
        "peak_day_p75": int(round(p75)),
        "peak_day_mean": float(np.mean(pd_valid)),
        "peak_value_mean": float(np.nanmean(peak_value[valid])),
    }


def bin_edges(days: np.ndarray) -> np.ndarray:
    """Histogram edges of BIN_WIDTH days spanning the observed peak dates."""
    lo = (days.min() // BIN_WIDTH) * BIN_WIDTH
    hi = ((days.max() // BIN_WIDTH) + 1) * BIN_WIDTH
    return np.arange(lo, hi + BIN_WIDTH, BIN_WIDTH)


def print_stats(st: dict) -> None:
    """Report both estimators and the per-year peaks, for the caption."""
    print()
    print("=" * 72)
    print("Climatological curve peak (max of the mean) vs per-year peaks "
          "(mean of the max)")
    print("=" * 72)
    print(f"  Climatological curve peaks on "
          f"{season_day_to_label(st['clim_peak_day'])} "
          f"(sd {st['clim_peak_day']}) at {st['clim_peak_value']:.2f}%")
    print(f"  Per-year peak date:  median "
          f"{season_day_to_label(st['peak_day_median'])}, "
          f"mean {season_day_to_label(int(round(st['peak_day_mean'])))}, "
          f"IQR {season_day_to_label(st['peak_day_p25'])} to "
          f"{season_day_to_label(st['peak_day_p75'])} "
          f"({st['peak_day_p75'] - st['peak_day_p25']} days)")
    print(f"  Per-year peak value: mean {st['peak_value_mean']:.2f}%")
    print(f"  Offset: the mean curve peaks "
          f"{st['clim_peak_day'] - st['peak_day_median']:+d} days later and "
          f"{st['clim_peak_value'] - st['peak_value_mean']:+.2f} points lower "
          f"than the median/mean per-year peak")
    print(f"  Skew check: per-year peak date mean minus median = "
          f"{st['peak_day_mean'] - st['peak_day_median']:+.1f} days "
          f"(near zero implies a symmetric distribution)")

    # Per-year listing, sorted by peak date, so any clustering or bimodality
    # in the peak-date distribution is visible without reading the histogram.
    print()
    print("=" * 72)
    print("Per-year peaks, sorted by date:")
    print("=" * 72)
    order = np.argsort(st["peak_day"][st["valid"]])
    years = st["years"][st["valid"]][order]
    days = st["peak_day"][st["valid"]][order]
    vals = st["peak_value"][st["valid"]][order]
    for y, d, v in zip(years, days, vals):
        print(f"  {int(y)}   sd {d:3d} ({season_day_to_label(int(d)):>6s})   "
              f"{v:5.2f}%")

    # Counts per bin, as a coarse check on unimodality.
    print()
    print("=" * 72)
    print(f"Per-year peak dates, {BIN_WIDTH}-day bins:")
    print("=" * 72)
    edges = bin_edges(days)
    counts, _ = np.histogram(days, bins=edges)
    for c, e0, e1 in zip(counts, edges[:-1], edges[1:]):
        bar = "#" * int(c)
        print(f"  {season_day_to_label(int(e0)):>6s} to "
              f"{season_day_to_label(int(e1 - 1)):>6s}  {c:2d}  {bar}")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(aggregated_path: Path, out_path: Path) -> None:
    ds = xr.open_dataset(aggregated_path)
    st = compute_stats(ds)

    season_day = st["season_day"]

    fig, (ax_a, ax_b) = plt.subplots(2, 1, sharex=True,
                                     constrained_layout=True)

    # ------------------------------------------------------------------
    # (a) Individual per-year curves with the climatological mean
    # ------------------------------------------------------------------
    # Panel (a) stacking, bottom to top: individual-year curves, the two
    # vertical peak lines (blue over orange), the climatological mean, then
    # the per-year peaks. Numbering starts at 2 because the default
    # axes.axisbelow ("line") draws the grid at zorder 1.5.
    h_years = None
    for i in range(st["sr"].shape[0]):
        line, = ax_a.plot(season_day, st["sr"][i],
                          color=YEAR_CURVE_COLOR,
                          linewidth=YEAR_CURVE_LINEWIDTH,
                          alpha=YEAR_CURVE_ALPHA, zorder=2,
                          label="Individual spawning years")
        if i == 0:
            h_years = line

    h_median = ax_a.axvline(
        st["peak_day_median"], color=COLOR_B, linestyle=PEAK_LINESTYLE,
        zorder=3,
        label=f"Median peak ({season_day_to_label(st['peak_day_median'])})")

    h_peak = ax_a.axvline(
        st["clim_peak_day"], color=COLOR_A, linestyle=PEAK_LINESTYLE,
        zorder=4,
        label=f"Mean-curve peak ({season_day_to_label(st['clim_peak_day'])})")

    h_mean, = ax_a.plot(season_day, st["mean_per_day"],
                        color=COLOR_A, zorder=5, label="32-year mean")

    # Each year's own peak, which is what panel (b) counts. Open markers in
    # panel (b)'s colour, opaque so they stay legible where several years
    # peak close together.
    h_marks, = ax_a.plot(
        st["peak_day"][st["valid"]], st["peak_value"][st["valid"]],
        linestyle="none", marker="o",
        markerfacecolor="white", markeredgecolor=COLOR_B,
        zorder=6, label="Per-year peak")

    ax_a.set_ylabel("Recruitment success rate (%)")
    ax_a.grid(True, axis="y", alpha=GRID_ALPHA)
    # Explicit handle order: the climatological curve and its peak first,
    # then the per-year peaks and their median.
    ax_a.legend(handles=[h_years, h_mean, h_peak, h_marks, h_median],
                loc="upper left", framealpha=LEGEND_FRAMEALPHA)
    framed_label(ax_a, "(a)")

    # ------------------------------------------------------------------
    # (b) Distribution of per-year peak release dates
    # ------------------------------------------------------------------
    days = st["peak_day"][st["valid"]]

    ax_b.axvspan(st["peak_day_p25"], st["peak_day_p75"],
                 color=COLOR_B, alpha=SPAN_ALPHA, linewidth=0, zorder=0,
                 label="Inter-quartile range")
    ax_b.hist(days, bins=bin_edges(days), rwidth=BAR_RWIDTH,
              color=COLOR_B, alpha=BAR_ALPHA, zorder=1)
    ax_b.axvline(st["peak_day_median"], color=COLOR_B,
                 linestyle=PEAK_LINESTYLE, zorder=2,
                 label=f"Median peak "
                       f"({season_day_to_label(st['peak_day_median'])})")

    ax_b.set_ylabel("Number of spawning years")
    ax_b.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax_b.grid(True, axis="y", alpha=GRID_ALPHA)
    ax_b.legend(loc="upper left", framealpha=LEGEND_FRAMEALPHA)
    framed_label(ax_b, "(b)")

    # Shared x-axis, matching F2.
    ax_b.set_xlim(season_day[0], LAST_SEASON_DAY)
    ax_b.set_xticks(TICK_DAYS)
    ax_b.set_xticklabels(TICK_LABELS)
    ax_b.set_xlabel("Release date")

    fig.savefig(out_path, dpi=500)
    plt.close(fig)
    print(f"Wrote {out_path}")

    print_stats(st)


def main():
    here = Path(__file__).resolve().parent
    default_aggregated = here.parent / "F2_phenology_curve" / "data" / "aggregated.nc"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--aggregated", type=Path, default=default_aggregated,
        help="path to the F2 aggregation (default: ../F2_phenology_curve/data/aggregated.nc)",
    )
    args = parser.parse_args()

    if not args.aggregated.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {args.aggregated}\n"
            f"Run F2_phenology_curve/aggregate.py first, or pass --aggregated."
        )

    plot(args.aggregated, here / "phenology_per_year.png")


if __name__ == "__main__":
    main()
