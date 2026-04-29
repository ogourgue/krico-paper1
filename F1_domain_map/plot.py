"""
F1 plot: domain map.

Method figure showing:
  - the computational domain extent (dashed gray boundary),
  - the bathymetry-defined spawning zone (1000-2000 m, blue fill),
  - the continental shelf (<1000 m, orange fill),
  - CCAMLR Subareas 48.1-48.6 and 88.3 (light gray outlines + framed labels),
  - cartopy land + coastlines.

Subarea 48.6 is split at 60°S for the spatial analyses in F4/F5
(48.6N north of 60°S, 48.6S south of 60°S). The split is visualized
here with a dashed gray line within the 48.6 polygon and two separate
centroid labels. The rationale for the split lives in the F1 caption.

Style aligned with F4 / F5: South Polar Stereographic projection,
Arial 9 pt, light-gray boundaries (`"0.7"`), axis spines pushed above
land/coastlines, panel margins via subplots_adjust.

Reads data/aggregated.nc produced by aggregate.py.
"""

from __future__ import annotations

import os
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from shapely.geometry import LineString, box
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    # PoC F1 used 6.5 × 3.3 for one panel of this domain extent — that's
    # the natural projected aspect we keep for F1.
    "figure.figsize": (6.5, 3.3),
})


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Domain extent (matches F4 / F5).
LON_MIN, LON_MAX = -115.0, 40.0
LAT_MIN, LAT_MAX = -78.0, -40.0

# CCAMLR statistical areas covered, with full geographic names for the legend.
# 48.6 is split at 60°S for spatial analyses (F4/F5); the legend uses the
# compact "Bouvet Island (North)" / "(South)" form. The 60°S split rationale
# lives in the F1 caption.
CCAMLR_AREA_NAMES = {
    "48.1": "Antarctic Peninsula",
    "48.2": "South Orkney Islands",
    "48.3": "South Georgia",
    "48.4": "South Sandwich Islands",
    "48.5": "Weddell Sea",
    "48.6N": "Northern Bouvet",
    "48.6S": "Southern Bouvet",
    "88.3": "Amundsen Sea",
}
# Underlying CCAMLR codes in the shapefile (48.6 unsplit).
CCAMLR_AREAS = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3"]

# Bathymetry depth thresholds for the two zones (meters).
SHELF_MAX = 1000.0
SPAWNING_MIN = 1000.0
SPAWNING_MAX = 2000.0

# Latitude split for 48.6 (consistent with F4/F5 ccamlr_summary.py).
SPLIT_LAT_486 = -60.0

# Number of points used to densify the 60°S divider line so it follows the
# parallel as a smooth curve when rendered in South Polar Stereographic.
N_DIVIDER_POINTS = 200

# Cartopy projection (matches F4 / F5).
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


def make_label_bbox(alpha: float = 1.0):
    """
    Bbox style matching the legend frame.

    `alpha` controls the box facecolor opacity. Defaults to fully opaque
    (matching F4 / F5 panel labels). Pass a lower value for the in-map
    CCAMLR area labels, where the underlying polygon and bathymetry
    should remain visible through the label frame.
    """
    facecolor = resolve_legend_color("legend.facecolor", "axes.facecolor")
    return dict(
        boxstyle="round",
        facecolor=mpl.colors.to_rgba(facecolor, alpha=alpha),
        edgecolor=resolve_legend_color("legend.edgecolor", "axes.edgecolor"),
        linewidth=mpl.rcParams["axes.linewidth"],
    )


def setup_polar_axes(ax):
    """Configure a polar-stereographic axis with land + coastline.

    Zorder stack (bottom to top):
      bathymetry fills (1) < CCAMLR + domain boundaries (2) < land (3) <
      coastlines (4) < axis spines (5) < CCAMLR labels (10)
    """
    land = cfeature.NaturalEarthFeature(
        "physical", "land", "50m", facecolor="0.85", edgecolor="none"
    )
    ax.add_feature(land, zorder=3)
    ax.coastlines(resolution="50m", linewidth=mpl.rcParams["axes.linewidth"], zorder=4)
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=DATA_TRANSFORM)
    for spine in ax.spines.values():
        spine.set_zorder(5)


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
    return ccamlr[ccamlr["GAR_Long_L"].isin(CCAMLR_AREAS)]


def split_486(geom_486):
    """
    Split the 48.6 polygon at 60°S into northern and southern halves.

    Returns (north_half, south_half) as shapely geometries. Either may be
    a Polygon or MultiPolygon depending on the original 48.6 shape.
    """
    minx, miny, maxx, maxy = geom_486.bounds
    pad = 1.0  # degrees
    north_box = box(minx - pad, SPLIT_LAT_486, maxx + pad, maxy + pad)
    south_box = box(minx - pad, miny - pad, maxx + pad, SPLIT_LAT_486)
    north_half = geom_486.intersection(north_box)
    south_half = geom_486.intersection(south_box)
    return north_half, south_half


def add_ccamlr_outlines(ax, ccamlr: gpd.GeoDataFrame) -> None:
    """Plot CCAMLR area boundaries (plain light-gray lines)."""
    for _, row in ccamlr.iterrows():
        ax.add_geometries(
            [row.geometry],
            crs=DATA_TRANSFORM,
            facecolor="none",
            edgecolor="0.7",
            linewidth=mpl.rcParams["axes.linewidth"],
            zorder=2,
        )


def add_486_split_line(ax, ccamlr: gpd.GeoDataFrame) -> None:
    """
    Draw a dotted gray line at 60°S inside the 48.6 polygon, marking the
    split between 48.6N and 48.6S used in F4/F5 spatial analyses.

    The line is densified (N_DIVIDER_POINTS along the lon span) so that
    when rendered on the South Polar Stereographic projection it follows
    the curved 60°S parallel rather than a chord. It is then clipped to
    the 48.6 polygon so it does not extend past the subarea boundary.
    Dotted style distinguishes this internal subdivision from the dashed
    style used for the external computational domain boundary.
    """
    rows = ccamlr[ccamlr["GAR_Long_L"] == "48.6"]
    if rows.empty:
        return
    geom_486 = unary_union(rows.geometry.values)

    # Build a densified horizontal line at 60°S spanning the full lon range
    # of 48.6, with a small pad so clipping by the polygon does not miss
    # the polygon boundary.
    minx, _, maxx, _ = geom_486.bounds
    pad = 1.0
    lons = np.linspace(minx - pad, maxx + pad, N_DIVIDER_POINTS)
    coords = list(zip(lons, np.full_like(lons, SPLIT_LAT_486)))
    full_line = LineString(coords)

    # Clip to the 48.6 polygon. The result may be a MultiLineString if the
    # polygon is multi-part or non-convex (the line exits and re-enters).
    clipped = full_line.intersection(geom_486)
    if clipped.is_empty:
        return

    if hasattr(clipped, "geoms"):
        segments = list(clipped.geoms)
    else:
        segments = [clipped]
    for seg in segments:
        if not isinstance(seg, LineString):
            continue
        xs, ys = seg.xy
        ax.plot(
            list(xs), list(ys),
            color="0.7",
            linestyle=":",
            linewidth=mpl.rcParams["axes.linewidth"],
            transform=DATA_TRANSFORM,
            zorder=2,
        )


def add_domain_boundary(ax) -> None:
    """Plot the model domain boundary as a dashed light-gray line."""
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


def add_ccamlr_labels(ax, ccamlr: gpd.GeoDataFrame) -> None:
    """Plot CCAMLR area code labels at polygon centroids, in semi-transparent
    framed boxes (so the underlying polygon and bathymetry remain visible).

    Subarea 48.6 gets two labels (48.6N, 48.6S) at the centroids of its
    two halves split at 60°S.
    """
    bbox = make_label_bbox(alpha=0.7)
    for _, row in ccamlr.iterrows():
        code = row["GAR_Long_L"]
        if code == "48.6":
            north_half, south_half = split_486(row.geometry)
            for half_geom, half_label in [(north_half, "48.6N"),
                                          (south_half, "48.6S")]:
                if half_geom.is_empty:
                    continue
                centroid = half_geom.centroid
                ax.text(
                    centroid.x, centroid.y, half_label,
                    transform=DATA_TRANSFORM,
                    ha="center", va="center",
                    bbox=bbox,
                    zorder=10,
                )
        else:
            centroid = row.geometry.centroid
            ax.text(
                centroid.x, centroid.y, code,
                transform=DATA_TRANSFORM,
                ha="center", va="center",
                bbox=bbox,
                zorder=10,
            )


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(aggregated_path: Path, out_path: Path) -> None:
    ds = xr.open_dataset(aggregated_path)
    hdept = ds["bathymetry_masked"].values
    nav_lon = ds["nav_lon"].values
    nav_lat = ds["nav_lat"].values

    # Build the two bathymetry zones from the masked field.
    shelf = np.where((hdept > 0) & (hdept < SHELF_MAX), hdept, np.nan)
    spawning = np.where(
        (hdept >= SPAWNING_MIN) & (hdept <= SPAWNING_MAX), hdept, np.nan
    )

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig = plt.figure()
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)
    ax = fig.add_subplot(1, 1, 1, projection=PROJECTION)
    setup_polar_axes(ax)

    # Bathymetry zones at zorder=1 (below boundaries and land).
    # Spawning zone in C0 (blue), continental shelf in C1 (orange).
    ax.contourf(
        nav_lon, nav_lat, spawning,
        levels=[SPAWNING_MIN, SPAWNING_MAX],
        colors=["C0"],
        alpha=0.5,
        transform=DATA_TRANSFORM,
        zorder=1,
    )
    ax.contourf(
        nav_lon, nav_lat, shelf,
        levels=[0, SHELF_MAX],
        colors=["C1"],
        alpha=0.5,
        transform=DATA_TRANSFORM,
        zorder=1,
    )

    # Boundaries.
    ccamlr = load_ccamlr_geometries()
    add_ccamlr_outlines(ax, ccamlr)
    add_486_split_line(ax, ccamlr)
    add_domain_boundary(ax)

    # CCAMLR labels (above land, above coastlines).
    add_ccamlr_labels(ax, ccamlr)

    # ------------------------------------------------------------------
    # Legends
    # ------------------------------------------------------------------
    # Legend 1 (upper right): zone / domain key.
    domain_patch = mpatches.Patch(
        facecolor="none", edgecolor="0.7", linestyle="--",
        label="Computational domain",
    )
    spawning_patch = mpatches.Patch(
        facecolor="C0", alpha=0.5, edgecolor="none",
        label="Spawning zone (1000–2000 m)",
    )
    shelf_patch = mpatches.Patch(
        facecolor="C1", alpha=0.5, edgecolor="none",
        label="Continental shelf (<1000 m)",
    )
    legend_zones = ax.legend(
        handles=[domain_patch, spawning_patch, shelf_patch],
        loc="upper right",
        framealpha=1.0,
    )

    # Legend 2 (lower left): CCAMLR subarea code -> geographic name.
    # Code-only labels with no handle, so the legend is a clean two-column
    # mapping table rather than a list of duplicate "none" patches.
    area_handles = [
        mpatches.Patch(facecolor="none", edgecolor="none",
                       label=f"{code}: {name}")
        for code, name in CCAMLR_AREA_NAMES.items()
    ]
    legend_areas = ax.legend(
        handles=area_handles,
        loc="lower left",
        handlelength=0,
        handletextpad=0,
        framealpha=1.0,
    )
    # Re-add the first legend (matplotlib drops it when a second is added).
    ax.add_artist(legend_zones)

    fig.savefig(out_path, dpi=500)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    here = Path(__file__).resolve().parent
    aggregated_path = here / "data" / "aggregated.nc"
    out_path = here / "domain_map.png"

    if not aggregated_path.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {aggregated_path}\n"
            f"Run aggregate.py first."
        )

    plot(aggregated_path, out_path)


if __name__ == "__main__":
    main()