"""
FS2 plot: mean circulation over the release depth band.

Reads data/aggregated.nc and produces currents.png:
  - Single polar-stereographic panel on the F4/F5 projection and extent.
  - Speed of the time-mean velocity vector as a log-scaled color field.
  - Streamlines of the same mean field over it.
  - CCAMLR outlines, 48.6 divider and domain boundary as in F1/F4/F5.

What is plotted, and what it does not mean:

  - The **speed of the mean vector**, |(u_bar, v_bar)|, not the mean of the
    speed. This measures persistent net transport, which is what displaces
    particles over 200 days. Regions that appear weak are regions of little
    *net* transport, not necessarily quiet regions: an energetic eddy field
    produces little net displacement and renders pale here.
  - **Velocity, not depth-integrated transport.** Transport scales with the
    available water column and would understate fast flow over the shelf,
    which is where the retention argument lives.
  - A thickness-weighted mean over **50-200 m, the release depth band**.
    Particles are advected in three dimensions over 200 days and sample well
    outside it; this is where every trajectory starts, not where larvae live.
  - Blank inside the 50 m isobath, where no model level falls in the band.
    White is "not represented", not "no flow".

The color scale is logarithmic because the field spans more than a decade:
the domain median is ~0.03 m/s against ~0.28 m/s at the 99th percentile and
0.72 m/s at the maximum. On a linear scale the ACC would be the only visible
feature and the entire shelf would render as uniform white.

Only the annual mean is plotted. The seasonal fields in the aggregation
differ from it by a few percent in every domain-wide statistic (medians
0.0288-0.0312 m/s), so four panels would show four nearly identical maps.
Pass --period to plot one of them anyway.

Inputs:  data/aggregated.nc  (from aggregate.py, run on the HPC)
Output:  currents.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import LogNorm
from shapely.geometry import LineString
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Style (see ../FIGURE_STYLE.md)
# ---------------------------------------------------------------------------

# One panel covering this extent is ~3.3 in tall at 6.5 in wide (the natural
# aspect used by the F4/F5 panels). Plus the colorbar (~0.2 in) and the
# bottom margin for its tick labels and label (~0.5 in).
mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "figure.figsize": (6.5, 4.0),
})


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Domain extent, matching the aggregation grid and F1/F4/F5.
LON_MIN, LON_MAX = -115.0, 40.0
LAT_MIN, LAT_MAX = -78.0, -40.0

# Cartopy projection (matches F1/F4/F5).
PROJECTION = ccrs.SouthPolarStereo(central_longitude=-37.5)
DATA_TRANSFORM = ccrs.PlateCarree()

# Latitude split for 48.6, consistent with F1/F4/F5.
SPLIT_LAT_486 = -60.0
N_DIVIDER_POINTS = 200

# Colour scale. Logarithmic; see the module docstring. VMIN is set an octave
# below the domain median so the shelf is resolved, VMAX just above the 99th
# percentile so only a fraction of a percent of cells saturate.
CMAP = "Blues"
VMIN, VMAX = 0.005, 0.5

# Streamlines. The mean field is subsampled before streamplot: at the native
# 1/12 degree the integration is slow and the result is illegible at domain
# scale.
STREAM_SUBSAMPLE = 6
STREAM_DENSITY = 2.5
STREAM_COLOR = "0.25"
STREAM_LINEWIDTH = 0.4
STREAM_ARROWSIZE = 0.6

PERIODS = ("annual", "DJF", "MAM", "JJA", "SON")


# ---------------------------------------------------------------------------
# Helpers shared with F4/F5 (kept identical so the maps read the same)
# ---------------------------------------------------------------------------

def resolve_legend_color(rc_key: str, axes_attr: str) -> str:
    """Resolve matplotlib's 'inherit' sentinel for legend colors."""
    val = mpl.rcParams[rc_key]
    if val == "inherit":
        return mpl.rcParams[axes_attr]
    return val


def setup_polar_axes(ax):
    """Configure a polar-stereographic axis with land + coastline.

    Zorder stack (bottom to top):
      data (1) < streamlines (1.5) < CCAMLR + domain boundaries (2) <
      land (3) < coastlines (4) < axis spines (5)
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
    repo_root = Path(__file__).resolve().parent.parent
    shp_path = repo_root / "ccamlr-data" / "CCAMLR_ASD_EPSG4326.shp"
    if not shp_path.exists():
        raise FileNotFoundError(f"CCAMLR shapefile not found: {shp_path}")
    ccamlr = gpd.read_file(shp_path)
    areas = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3"]
    return ccamlr[ccamlr["GAR_Long_L"].isin(areas)]


def add_486_split_line(ax, ccamlr: gpd.GeoDataFrame) -> None:
    """Dotted 60 degrees S divider inside 48.6, as in F1/F4/F5."""
    rows = ccamlr[ccamlr["GAR_Long_L"] == "48.6"]
    if rows.empty:
        return
    geom_486 = unary_union(rows.geometry.values)

    minx, _, maxx, _ = geom_486.bounds
    pad = 1.0
    lons = np.linspace(minx - pad, maxx + pad, N_DIVIDER_POINTS)
    coords = list(zip(lons, np.full_like(lons, SPLIT_LAT_486)))
    clipped = LineString(coords).intersection(geom_486)
    if clipped.is_empty:
        return

    segments = list(clipped.geoms) if hasattr(clipped, "geoms") else [clipped]
    for seg in segments:
        if not isinstance(seg, LineString):
            continue
        xs, ys = seg.xy
        ax.plot(
            list(xs), list(ys),
            color="0.7", linestyle=":",
            linewidth=mpl.rcParams["axes.linewidth"],
            transform=DATA_TRANSFORM, zorder=2,
        )


def add_overlays(ax, ccamlr: gpd.GeoDataFrame) -> None:
    """CCAMLR outlines (plain), 48.6 divider (dotted), domain (dashed)."""
    for _, row in ccamlr.iterrows():
        ax.add_geometries(
            [row.geometry], crs=DATA_TRANSFORM,
            facecolor="none", edgecolor="0.7",
            linewidth=mpl.rcParams["axes.linewidth"], zorder=2,
        )

    add_486_split_line(ax, ccamlr)

    n = 100
    lon_boundary = np.concatenate([
        np.linspace(LON_MIN, LON_MAX, n), np.full(n, LON_MAX),
        np.linspace(LON_MAX, LON_MIN, n), np.full(n, LON_MIN),
    ])
    lat_boundary = np.concatenate([
        np.full(n, LAT_MIN), np.linspace(LAT_MIN, LAT_MAX, n),
        np.full(n, LAT_MAX), np.linspace(LAT_MAX, LAT_MIN, n),
    ])
    ax.plot(
        lon_boundary, lat_boundary,
        color="0.7", linestyle="--",
        linewidth=mpl.rcParams["axes.linewidth"],
        transform=DATA_TRANSFORM, zorder=2,
    )


# ---------------------------------------------------------------------------
# Grid handling
# ---------------------------------------------------------------------------

def grid_1d(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """
    Reduce the 2-D nav_lon / nav_lat coordinate arrays to 1-D axes.

    GLORYS12 is a regular longitude-latitude grid, so every row of nav_lon is
    identical and every column of nav_lat is. Verified rather than assumed,
    because streamplot needs 1-D axes and would silently distort the field if
    the grid were curvilinear.
    """
    lon2d = ds["nav_lon"].values
    lat2d = ds["nav_lat"].values
    lon = lon2d[0, :]
    lat = lat2d[:, 0]

    dlon = float(np.nanmax(np.abs(lon2d - lon[None, :])))
    dlat = float(np.nanmax(np.abs(lat2d - lat[:, None])))
    tol = 1e-4
    if dlon > tol or dlat > tol:
        raise ValueError(
            "nav_lon/nav_lat are not separable to %.1e degrees "
            "(max deviation %.3e / %.3e): the grid is curvilinear and "
            "streamplot cannot take 1-D axes." % (tol, dlon, dlat)
        )
    return lon, lat


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot(aggregated_path: Path, out_path: Path, period: str) -> None:
    ds = xr.open_dataset(aggregated_path)

    u = ds[f"u_{period}"].values
    v = ds[f"v_{period}"].values
    n_days = ds[f"u_{period}"].attrs.get("n_days", "?")
    lon, lat = grid_1d(ds)

    speed = np.hypot(u, v)
    valid = np.isfinite(speed)

    print(f"Period {period}: n_days = {n_days}")
    print("  speed: median %.4f  p90 %.4f  p99 %.4f  max %.4f m/s" % (
        np.nanmedian(speed), np.nanpercentile(speed, 90),
        np.nanpercentile(speed, 99), np.nanmax(speed)))
    print("  colour scale: log, %.3f to %.3f m/s" % (VMIN, VMAX))
    print("  cells below vmin: %.1f%%   above vmax: %.2f%%" % (
        100 * np.mean(speed[valid] < VMIN), 100 * np.mean(speed[valid] > VMAX)))

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    fig = plt.figure()
    fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.12)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 0.05], hspace=0.05)

    ax = fig.add_subplot(gs[0, 0], projection=PROJECTION)
    setup_polar_axes(ax)

    mappable = ax.pcolormesh(
        lon, lat, np.where(valid, speed, np.nan),
        cmap=CMAP, norm=LogNorm(vmin=VMIN, vmax=VMAX),
        transform=DATA_TRANSFORM, shading="nearest", zorder=1,
    )

    # Streamlines of the same mean field. Subsampled, and NaN replaced by zero
    # so the integrator stops at the mask rather than failing on it.
    s = STREAM_SUBSAMPLE
    ax.streamplot(
        lon[::s], lat[::s],
        np.where(valid, u, 0.0)[::s, ::s],
        np.where(valid, v, 0.0)[::s, ::s],
        transform=DATA_TRANSFORM,
        density=STREAM_DENSITY,
        color=STREAM_COLOR,
        linewidth=STREAM_LINEWIDTH,
        arrowsize=STREAM_ARROWSIZE,
        zorder=1.5,
    )

    ccamlr = load_ccamlr_geometries()
    add_overlays(ax, ccamlr)

    # ------------------------------------------------------------------
    # Colorbar, inset to ~90% of the panel width so the end tick labels have
    # room to render (same device as F4/F5).
    # ------------------------------------------------------------------
    inner = gs[1, 0].subgridspec(
        nrows=1, ncols=3, width_ratios=[0.05, 0.9, 0.05], wspace=0.0)
    cax = fig.add_subplot(inner[0, 1])
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", extend="both")
    cb.set_label("Mean current speed (m s$^{-1}$)")

    fig.savefig(out_path, dpi=500)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    here = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--period", default="annual", choices=PERIODS,
                        help="which mean to plot (default: annual)")
    args = parser.parse_args()

    aggregated_path = here / "data" / "aggregated.nc"
    if not aggregated_path.exists():
        raise FileNotFoundError(
            f"Aggregated file not found: {aggregated_path}\n"
            f"Run aggregate.py on the HPC and copy data/aggregated.nc back."
        )

    suffix = "" if args.period == "annual" else f"_{args.period}"
    plot(aggregated_path, here / f"currents{suffix}.png", args.period)


if __name__ == "__main__":
    main()
