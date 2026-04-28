"""
F5 plot: destination maps by fate.

Reads data/aggregated.nc and produces destination_maps.png:
  - 3 × 2 grid of polar-stereographic maps, one panel per outcome.
  - Success in top-left, then killed_M*, then exited_domain.
  - Three colormaps echoing F3: Greens (success), Reds (4 killed),
    Greys (exited_domain).
  - Three shared colorbars at the bottom.

Identical layout to F4. The only difference is the input data: the
density grids in data/aggregated.nc reflect fate positions
(final_lon / final_lat) rather than release positions.
"""

from __future__ import annotations

import os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Style (see ../FIGURE_STYLE.md, plus per-figure deviations documented here)
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    # Panel natural aspect from the PoC F1: 6.5 in wide × 3.3 in tall for
    # one panel covering this extent. With 2 cols × 3 rows of equal panels
    # and reduced left/right margins (subplots_adjust), each panel ends up
    # ~3.1 in wide × ~1.55 in tall, total ~4.65 in for panels. Plus the
    # colorbar row (~0.2 in) and bottom margin reserved for cbar tick
    # labels and axis labels (~0.55 in).
    "figure.figsize": (6.5, 5.5),
})


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Domain extent for cropping the map (matches the aggregation grid).
LON_MIN, LON_MAX = -115.0, 40.0
LAT_MIN, LAT_MAX = -78.0, -40.0

# Panels: outcome name -> (row, col) in the 3 × 2 grid.
# Success top-left; killed_M* in time-of-action order; exited_domain last.
PANEL_LAYOUT = [
    ("success",                (0, 0)),
    ("killed_M1",              (0, 1)),
    ("killed_M4",              (1, 0)),
    ("killed_M5_no_FIV",       (1, 1)),
    ("killed_M5_not_on_shelf", (2, 0)),
    ("exited_domain",          (2, 1)),
]

PANEL_LABELS = {
    "success":                "Success",
    "killed_M1":              "Killed: M1 (no spawning)",
    "killed_M4":              "Killed: M4 (calyptope starvation)",
    "killed_M5_no_FIV":       "Killed: M5 (no FIV)",
    "killed_M5_not_on_shelf": "Killed: M5 (off-shelf)",
    "exited_domain":          "Exited domain",
}

# Color group per outcome; determines colormap and shared scale.
COLOR_GROUP = {
    "success":                "success",
    "killed_M1":              "killed",
    "killed_M4":              "killed",
    "killed_M5_no_FIV":       "killed",
    "killed_M5_not_on_shelf": "killed",
    "exited_domain":          "exited",
}

CMAPS = {
    "success": "Greens",
    "killed":  "Reds",
    "exited":  "Greys",
}

GROUP_LABELS = {
    "success": "Success",
    "killed":  "Killed",
    "exited":  "Exited domain",
}

# Percentile used to set per-group vmax. Clipping at the 99th percentile
# brightens the bulk of the signal at the cost of saturating a few hotspot
# cells.
COLOR_PERCENTILE = 99.0

# Cartopy projection (matches PoC F1).
PROJECTION = ccrs.SouthPolarStereo(central_longitude=-37.5)
DATA_TRANSFORM = ccrs.PlateCarree()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def resolve_legend_color(rc_key: str, axes_attr: str) -> str:
    """Resolve matplotlib's 'inherit' sentinel for legend colors."""
    val = mpl.rcParams[rc_key]
    if val == "inherit":
        return mpl.rcParams[axes_attr]
    return val


def make_label_bbox():
    """Bbox style matching the legend frame (white fill, gray border, opaque)."""
    return dict(
        boxstyle="round",
        facecolor=resolve_legend_color("legend.facecolor", "axes.facecolor"),
        edgecolor=resolve_legend_color("legend.edgecolor", "axes.edgecolor"),
        linewidth=mpl.rcParams["axes.linewidth"],
    )


def add_panel_labels(ax, letter: str, fate_name: str):
    """
    Add a single framed label combining the panel letter and fate name
    in the upper-right of the panel.
    """
    ax.text(
        0.97, 0.92, f"({letter}) {fate_name}",
        transform=ax.transAxes,
        ha="right", va="top",
        bbox=make_label_bbox(),
        zorder=10,
    )


def setup_polar_axes(ax):
    """Configure a polar-stereographic axis with land + coastline."""
    land = cfeature.NaturalEarthFeature(
        "physical", "land", "50m", facecolor="0.85", edgecolor="none"
    )
    ax.add_feature(land, zorder=3)
    ax.coastlines(resolution="50m", linewidth=mpl.rcParams["axes.linewidth"], zorder=4)
    # Crop to the domain extent.
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=DATA_TRANSFORM)


def load_ccamlr_geometries() -> gpd.GeoDataFrame:
    """Load CCAMLR statistical area polygons for the project domain."""
    krico_root = os.environ.get("KRICO_ROOT")
    if not krico_root:
        raise RuntimeError("KRICO_ROOT environment variable not set.")
    shp_path = (
        Path(krico_root) / "Pre" / "ccamlr-data" / "geographical_data"
        / "asd" / "CCAMLR_ASD_EPSG4326.shp"
    )
    if not shp_path.exists():
        raise FileNotFoundError(f"CCAMLR shapefile not found: {shp_path}")
    ccamlr = gpd.read_file(shp_path)
    areas = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3"]
    return ccamlr[ccamlr["GAR_Long_L"].isin(areas)]


def add_overlays(ax, ccamlr: gpd.GeoDataFrame) -> None:
    """Add CCAMLR area outlines (plain) and model domain boundary (dashed).

    Both are drawn above data but below continents (zorder=2).
    """
    # CCAMLR area boundaries: plain light-gray lines.
    for _, row in ccamlr.iterrows():
        ax.add_geometries(
            [row.geometry],
            crs=DATA_TRANSFORM,
            facecolor="none",
            edgecolor="0.7",
            linewidth=mpl.rcParams["axes.linewidth"],
            zorder=2,
        )

    # Model domain boundary: dashed light-gray, densified for smooth polar curves.
    n = 100
    lon_boundary = np.concatenate([
        np.linspace(LON_MIN, LON_MAX, n),
        np.full(n, LON_MAX),
        np.linspace(LON_MAX, LON_MIN, n),
        np.full(n, LON_MIN),
    ])
    lat_boundary = np.concatenate([
        np.full(n, LAT_MIN),
        np.linspace(LAT_MIN, LAT_MAX, n),
        np.full(n, LAT_MAX),
        np.linspace(LAT_MAX, LAT_MIN, n),
    ])
    ax.plot(
        lon_boundary, lat_boundary,
        color="0.7",
        linestyle="--",
        linewidth=mpl.rcParams["axes.linewidth"],
        transform=DATA_TRANSFORM,
        zorder=2,
    )


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(aggregated_path: Path, out_path: Path) -> None:
    ds = xr.open_dataset(aggregated_path)
    density = ds["density"].values                          # (outcome, lat, lon)
    lon_centers = ds["lon"].values
    lat_centers = ds["lat"].values
    outcome_names = list(ds["outcome"].values.astype(str))

    # Cell edges for pcolormesh.
    dlon = lon_centers[1] - lon_centers[0]
    dlat = lat_centers[1] - lat_centers[0]
    lon_edges = np.concatenate([lon_centers - dlon / 2, [lon_centers[-1] + dlon / 2]])
    lat_edges = np.concatenate([lat_centers - dlat / 2, [lat_centers[-1] + dlat / 2]])

    # ------------------------------------------------------------------
    # Compute a single shared vmax across all 6 panels, using the 99th
    # percentile of all nonzero density cells pooled together.
    # All three colormaps use the same numeric scale, so a given hue
    # intensity means the same density everywhere in the figure.
    # ------------------------------------------------------------------
    panel_indices = [outcome_names.index(name) for name, _ in PANEL_LAYOUT]
    pooled = np.concatenate([
        density[idx][density[idx] > 0].ravel()
        for idx in panel_indices
    ])
    if pooled.size == 0:
        raise RuntimeError("All panels are empty.")
    shared_vmax = float(np.percentile(pooled, COLOR_PERCENTILE))
    print(f"Shared vmax (P{COLOR_PERCENTILE:.0f} of all nonzero cells): {shared_vmax:.4e}")
    group_max = {"success": shared_vmax, "killed": shared_vmax, "exited": shared_vmax}

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    # Use a GridSpec with one extra row at the bottom for colorbars.
    fig = plt.figure()
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.10)
    gs = fig.add_gridspec(
        nrows=4, ncols=2,
        height_ratios=[1, 1, 1, 0.12],
        hspace=0.02, wspace=0.02,
    )

    # ------------------------------------------------------------------
    # Load CCAMLR area geometries once (shared across panels).
    # ------------------------------------------------------------------
    ccamlr = load_ccamlr_geometries()

    # Track an example mappable per group for the colorbars.
    group_mappable = {}
    panel_letters = "abcdef"

    for i, (name, (row, col)) in enumerate(PANEL_LAYOUT):
        idx = outcome_names.index(name)
        group = COLOR_GROUP[name]
        cmap = CMAPS[group]
        vmax = group_max[group]

        ax = fig.add_subplot(gs[row, col], projection=PROJECTION)
        setup_polar_axes(ax)

        # Mask cells with zero density so the white end of the colormap
        # (≈ background) is not drawn over the land.
        d = density[idx]
        d_masked = np.where(d > 0, d, np.nan)

        mappable = ax.pcolormesh(
            lon_edges, lat_edges, d_masked,
            cmap=cmap,
            vmin=0.0, vmax=vmax,
            transform=DATA_TRANSFORM,
            shading="flat",
            zorder=1,
        )
        group_mappable.setdefault(group, mappable)

        add_overlays(ax, ccamlr)
        add_panel_labels(ax, panel_letters[i], PANEL_LABELS[name])

    # ------------------------------------------------------------------
    # Three colorbars at the bottom, one per group, side by side.
    # ------------------------------------------------------------------
    # We use a sub-GridSpec inside the bottom row so the three bars
    # share alignment with each other and with the panels above.
    cbar_gs = gs[3, :].subgridspec(
        nrows=1, ncols=3, wspace=0.4,
    )
    for j, group in enumerate(("success", "killed", "exited")):
        cax = fig.add_subplot(cbar_gs[0, j])
        cb = fig.colorbar(group_mappable[group], cax=cax, orientation="horizontal")
        cb.set_label(f"{GROUP_LABELS[group]} density")
        # Force scientific notation with the multiplier on the right of the
        # axis (offset text), and shrink it so it doesn't collide with the
        # axis label below.
        cb.ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3))
        cb.ax.xaxis.get_offset_text().set_size(mpl.rcParams["font.size"] - 1)

    fig.savefig(out_path, dpi=500)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    here = Path(__file__).resolve().parent
    aggregated_path = here / "data" / "aggregated.nc"
    out_path = here / "destination_maps.png"

    if not aggregated_path.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {aggregated_path}\n"
            f"Run aggregate.py first."
        )

    plot(aggregated_path, out_path)


if __name__ == "__main__":
    main()
