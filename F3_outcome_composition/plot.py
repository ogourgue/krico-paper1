"""
F3 plot: outcome composition.

Reads data/aggregated.nc and produces outcome_composition.png:
  - Stacked-area showing the climatological fractional breakdown of all
    outcomes vs. release date.
  - 32-year mean across spawning years.
  - Censored particles are folded into killed_M6_no_advance (they share the
    same "alive at end of tracking, no advance detected" condition; censored
    is the subset where SIC was rising at the cutoff, classification
    provisional). The killed_M6 vs censored breakdown lives in the SI.
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
    "figure.figsize": (6.5, 4.875),
})

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# X-axis ticks: first of month only.
TICK_DAYS = [16, 47, 78, 106]
TICK_LABELS = ["Dec 1", "Jan 1", "Feb 1", "Mar 1"]

# Layer order, bottom to top of the stack.
# exited_domain at the bottom: it is a modeling-domain limitation rather
# than a biological outcome, so it sits as a "noise floor" below the
# biological layers. Then killed_M* in chronological order, with success
# at the very top.
LAYER_ORDER = [
    "exited_domain",
    "killed_M1",
    "killed_M4",
    "killed_M5_no_FIV",
    "killed_M5_not_on_shelf",
    "killed_M6_no_advance",
    "success",
]

# Pretty labels for the legend.
LAYER_LABELS = {
    "killed_M1": "Killed: M1 (no spawning)",
    "killed_M4": "Killed: M4 (calyptope starvation)",
    "killed_M5_no_FIV": "Killed: M5 (no FIV)",
    "killed_M5_not_on_shelf": "Killed: M5 (off-shelf)",
    "killed_M6_no_advance": "Killed: M6 (no advance)",
    "exited_domain": "Exited domain",
    "success": "Success",
}


def layer_colors(killed_layers: list[str]) -> dict[str, tuple]:
    """
    Build a color mapping for all layers.

    - success         -> matplotlib C2 (green).
    - killed_M*       -> n-step gradient sampled from YlOrRd (truncated to
                         avoid the very-light end).
    - exited_domain   -> gray.
    - others          -> default color cycle as fallback (should not occur
                         given current LAYER_ORDER).
    """
    colors: dict[str, tuple] = {}
    colors["success"] = mpl.colors.to_rgba("C2")
    colors["exited_domain"] = mpl.colors.to_rgba("0.7")

    cmap = mpl.colormaps["YlOrRd"]
    n = len(killed_layers)
    if n == 0:
        return colors
    # Sample from 0.35 to 0.85 so the lightest tone is visible against white
    # and the darkest is not fully saturated red.
    sample_points = np.linspace(0.35, 0.85, n)
    for layer, t in zip(killed_layers, sample_points):
        colors[layer] = cmap(t)
    return colors


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(aggregated_path: Path, out_path: Path) -> None:
    ds = xr.open_dataset(aggregated_path)
    counts = ds["counts"].values                           # (year, season_day, outcome)
    total = ds["total"].values                              # (year, season_day)
    has_data = ds["has_data"].values
    season_day = ds["season_day"].values
    outcome_names = list(ds["outcome"].values.astype(str))

    n_years, n_days, n_outcomes = counts.shape

    # Fold censored into killed_M6_no_advance per the methodology decision.
    # Censored particles share the same situation as killed_M6 (alive at end
    # of tracking, no detected sea-ice advance event); the only distinction
    # is that censored marks the subset whose classification is provisional
    # because SIC was rising at the cutoff. We treat them as a subcategory
    # of killed_M6 in the figure; the breakdown lives in the SI.
    censored_idx = outcome_names.index("censored")
    m6_idx = outcome_names.index("killed_M6_no_advance")
    counts[:, :, m6_idx] += counts[:, :, censored_idx]
    counts[:, :, censored_idx] = 0

    # --- Climatological mean fractions ------------------------------------
    # For each (year, season_day) cell, fraction = count / total.
    # Then average across years (only over cells with data).
    with np.errstate(invalid="ignore", divide="ignore"):
        fractions = np.where(
            has_data[:, :, None],
            counts / np.maximum(total[:, :, None], 1),
            np.nan,
        )
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mean_fraction = np.nanmean(fractions, axis=0)       # (season_day, outcome)

    # --- All layers in the stack ------------------------------------------
    grand_total = counts.sum()
    overall_fraction_pct = 100.0 * counts.sum(axis=(0, 1)) / max(grand_total, 1)
    overall_by_name = dict(zip(outcome_names, overall_fraction_pct))

    print("Overall fractions (dataset-wide):")
    for name, pct in overall_by_name.items():
        print(f"  {name:28s} {pct:6.3f}%")
    print()

    layers = list(LAYER_ORDER)

    # Killed_M* layers (chronologically ordered slice of `layers`).
    killed_layers = [n for n in layers if n.startswith("killed_M")]
    colors = layer_colors(killed_layers)

    # --- Build stack arrays ----------------------------------------------
    # In percentages (0..100) for readability.
    layer_data = {
        name: 100.0 * mean_fraction[:, outcome_names.index(name)]
        for name in layers
    }

    # --- Figure -----------------------------------------------------------
    fig, ax = plt.subplots(constrained_layout=True)

    ax.stackplot(
        season_day,
        *[layer_data[name] for name in layers],
        labels=[LAYER_LABELS[name] for name in layers],
        colors=[colors[name] for name in layers],
        edgecolor="none",
    )

    ax.set_xlim(season_day[0], 119)
    ax.set_ylim(0, 100)
    ax.set_xticks(TICK_DAYS)
    ax.set_xticklabels(TICK_LABELS)
    ax.set_xlabel("Release date")
    ax.set_ylabel("Cumulative fraction (%)")
    # Reverse the legend order so it reads top-to-bottom matching the stack
    # (success on top in the figure -> success first in the legend).
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[::-1], labels[::-1],
        loc="lower left",
        framealpha=1.0,
    )

    fig.savefig(out_path, dpi=500)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    here = Path(__file__).resolve().parent
    aggregated_path = here / "data" / "aggregated.nc"
    out_path = here / "outcome_composition.png"

    if not aggregated_path.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {aggregated_path}\n"
            f"Run aggregate.py first."
        )

    plot(aggregated_path, out_path)


if __name__ == "__main__":
    main()