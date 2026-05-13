"""
F1 aggregation: mask GLORYS12 bathymetry within CCAMLR statistical areas.

The pixel-by-pixel point-in-polygon test takes ~1 minute on the full
GLORYS12 grid. Doing it once and writing the masked field to NetCDF
makes the plotting script fast for iterations.

Inputs:
  $KRICO_GLORYS12/glorys12_bathymetry.nc
  ccamlr-data/CCAMLR_ASD_EPSG4326.shp

Output:
  data/aggregated.nc
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Point
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CCAMLR_AREAS = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3"]


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(bathy_path: Path, ccamlr_shp_path: Path) -> xr.Dataset:
    """
    Load bathymetry and CCAMLR polygons, build a point-in-polygon mask
    for the project areas, and return the masked bathymetry.
    """
    print(f"Reading bathymetry: {bathy_path}")
    bathy_ds = xr.open_dataset(bathy_path)
    hdept = bathy_ds["hdept"].isel(t=0).values.astype(np.float32)
    nav_lon = bathy_ds["nav_lon"].values.astype(np.float32)
    nav_lat = bathy_ds["nav_lat"].values.astype(np.float32)
    bathy_ds.close()

    print(f"Reading CCAMLR shapefile: {ccamlr_shp_path}")
    ccamlr = gpd.read_file(ccamlr_shp_path)
    ccamlr = ccamlr[ccamlr["GAR_Long_L"].isin(CCAMLR_AREAS)]
    ccamlr_union = unary_union(ccamlr.geometry)

    print("Masking bathymetry pixels within CCAMLR areas (this takes a while)...")
    mask = np.zeros_like(hdept, dtype=bool)
    n_rows = nav_lon.shape[0]
    for i in range(n_rows):
        for j in range(nav_lon.shape[1]):
            if ccamlr_union.contains(Point(nav_lon[i, j], nav_lat[i, j])):
                mask[i, j] = True
        if (i + 1) % 100 == 0:
            print(f"  row {i + 1} / {n_rows}")
    print(f"  {int(mask.sum()):,} cells inside CCAMLR areas")

    masked_hdept = np.where(mask, hdept, np.nan).astype(np.float32)

    ds_out = xr.Dataset(
        data_vars={
            "bathymetry_masked": (
                ("y", "x"),
                masked_hdept,
                {
                    "long_name": "GLORYS12 bathymetry, NaN outside CCAMLR project areas",
                    "units": "m",
                    "positive": "down",
                },
            ),
            "nav_lon": (
                ("y", "x"),
                nav_lon,
                {"long_name": "longitude", "units": "degrees_east"},
            ),
            "nav_lat": (
                ("y", "x"),
                nav_lat,
                {"long_name": "latitude", "units": "degrees_north"},
            ),
        },
        attrs={
            "title": "F1 domain map aggregation — masked GLORYS12 bathymetry",
            "description": (
                "GLORYS12 bathymetry masked to KRICO project CCAMLR areas "
                "(48.1–48.6 and 88.3) for use as input to the F1 domain map plot."
            ),
            "ccamlr_areas": ", ".join(CCAMLR_AREAS),
        },
    )
    return ds_out


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # Bathymetry: external dependency, set via KRICO_GLORYS12 env var.
    # Only needed to re-run F1 aggregation from scratch; data/aggregated.nc
    # is committed to the repo for figure reproduction without HPC access.
    krico_glorys12 = os.environ.get("KRICO_GLORYS12")
    if not krico_glorys12:
        print(
            "ERROR: KRICO_GLORYS12 environment variable not set.\n"
            "Set it to the GLORYS12 preprocessing output directory, e.g.:\n"
            "  export KRICO_GLORYS12=/path/to/Pre/GLORYS12",
            file=sys.stderr,
        )
        sys.exit(1)

    bathy_path = Path(krico_glorys12) / "glorys12_bathymetry.nc"
    if not bathy_path.exists():
        print(f"ERROR: bathymetry file not found: {bathy_path}", file=sys.stderr)
        sys.exit(1)

    # CCAMLR shapefile: bundled in the repo at ccamlr-data/.
    repo_root = Path(__file__).resolve().parent.parent
    ccamlr_shp_path = repo_root / "ccamlr-data" / "CCAMLR_ASD_EPSG4326.shp"
    if not ccamlr_shp_path.exists():
        print(f"ERROR: CCAMLR shapefile not found: {ccamlr_shp_path}", file=sys.stderr)
        sys.exit(1)

    ds = aggregate(bathy_path, ccamlr_shp_path)

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "aggregated.nc"

    print(f"Writing {out_path}")
    encoding = {var: {"zlib": True, "complevel": 5} for var in ds.data_vars}
    ds.to_netcdf(out_path, encoding=encoding)

    print("Done.")


if __name__ == "__main__":
    main()
