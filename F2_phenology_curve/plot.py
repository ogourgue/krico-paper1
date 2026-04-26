"""
F2 plot: phenology curve.

Reads data/aggregated.nc and produces phenology_curve.png:
  - 30-year mean success rate vs. release season-day, with 5th-95th ribbon
  - Climatological peak day annotated
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Style (see ../FIGURE_STYLE.md)
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "figure.figsize": (6.5, 4.875),       # 6.5 in wide, 4:3 aspect
})


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Reference dates for x-axis tick labels (season-day -> calendar label).
# season_day 0 = Nov 15. Show first of each month only.
TICK_DAYS = [16, 47, 78, 106]
TICK_LABELS = ["Dec 1", "Jan 1", "Feb 1", "Mar 1"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def season_day_to_label(day: int) -> str:
    """Calendar label for a season day index, for annotation text."""
    # November has 30 days, December 31, January 31, February 28, March up to 15.
    # season_day 0 is Nov 15.
    if day <= 15:
        return f"Nov {15 + day}"
    if day <= 46:                                  # Dec 1 .. Dec 31
        return f"Dec {day - 15}"
    if day <= 77:                                  # Jan 1 .. Jan 31
        return f"Jan {day - 46}"
    if day <= 105:                                 # Feb 1 .. Feb 28
        return f"Feb {day - 77}"
    return f"Mar {day - 105}"                      # Mar 1 .. Mar 15


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(aggregated_path: Path, out_path: Path) -> None:
    ds = xr.open_dataset(aggregated_path)
    sr = ds["success_rate"].values * 100.0         # percent
    has_data = ds["has_data"].values
    season_day = ds["season_day"].values

    sr_masked = np.where(has_data, sr, np.nan)

    # Climatological statistics across years (axis 0). Suppress all-NaN
    # warnings — for cohorts that are missing across all years (e.g. Nov 31,
    # Feb 30, or sparse synthetic test data), nanmean/nanpercentile of an
    # all-NaN column raises a benign warning.
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean = np.nanmean(sr_masked, axis=0)
        p5 = np.nanpercentile(sr_masked, 5, axis=0)
        p95 = np.nanpercentile(sr_masked, 95, axis=0)

    # Climatological peak.
    peak_day = int(season_day[np.nanargmax(mean)])
    peak_value = float(np.nanmax(mean))

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig, ax = plt.subplots(constrained_layout=True)

    ax.fill_between(
        season_day, p5, p95,
        alpha=0.2, color="C0", linewidth=0,
        label="5th-95th percentile",
    )
    ax.plot(
        season_day, mean,
        color="C0",
        label="32-year mean",
    )

    # Mark the climatological peak with a vertical line styled like the axes.
    axes_lw = mpl.rcParams["axes.linewidth"]
    ax.axvline(peak_day, color="black", linewidth=axes_lw, linestyle="--")

    # Peak label above the mean line, horizontally centered on the peak date.
    # Frame styled to match the legend visually. Note: matplotlib's defaults
    # for legend.facecolor and legend.edgecolor are the literal string
    # "inherit", which means "use the axes face/edge color" — a legend-only
    # sentinel. We resolve it explicitly for the bbox patch.
    def resolve_legend_color(rc_key, axes_attr):
        val = mpl.rcParams[rc_key]
        if val == "inherit":
            return mpl.rcParams[axes_attr]
        return val

    legend_facecolor = resolve_legend_color("legend.facecolor", "axes.facecolor")
    legend_edgecolor = resolve_legend_color("legend.edgecolor", "axes.edgecolor")

    bbox = dict(
        boxstyle="round",
        facecolor=legend_facecolor,
        edgecolor=legend_edgecolor,
        linewidth=axes_lw,
    )
    ax.annotate(
        f"Peak: {season_day_to_label(peak_day)} ({peak_value:.1f}%)",
        xy=(peak_day, peak_value),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center", va="bottom",
        color="black",
        bbox=bbox,
    )

    # Last release date is Mar 14 (season_day 119), not Mar 15 (120).
    ax.set_xlim(season_day[0], 119)
    ax.set_xticks(TICK_DAYS)
    ax.set_xticklabels(TICK_LABELS)
    ax.set_xlabel("Release date")
    ax.set_ylabel("Recruitment success rate (%)")
    ax.grid(True, axis="y", alpha=0.5)
    ax.legend(loc="upper left", framealpha=1.0)

    fig.savefig(out_path, dpi=500)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    here = Path(__file__).resolve().parent
    aggregated_path = here / "data" / "aggregated.nc"
    out_path = here / "phenology_curve.png"

    if not aggregated_path.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {aggregated_path}\n"
            f"Run aggregate.py first."
        )

    plot(aggregated_path, out_path)


if __name__ == "__main__":
    main()