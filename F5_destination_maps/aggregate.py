"""
F5 aggregation: per-outcome 2D fate-position density.

Reads per-particle recruitment outcomes from all cohort files in the
recruitment data directory. For each particle, takes its final_lon /
final_lat (position at the moment fate was determined) and bins onto a
regular 1° × 1° lon/lat grid, separately for each outcome category.

Output is normalized within each outcome (densities sum to 1 across cells
for each outcome): "of all particles that ended in this fate, this is the
fractional spatial fate-position distribution."

Inputs:  $KRICO_ROOT/Post/Production/recruitment/data/YYYY_MM_DD.nc
Output:  data/aggregated.nc
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Outcome integer flags from the recruitment pipeline (CF flag values).
OUTCOME_NAMES = (
    "success",                 # 0
    "censored",                # 1
    "killed_M1",               # 2
    "killed_M4",               # 3
    "killed_M5_no_FIV",        # 4
    "killed_M5_not_on_shelf",  # 5
    "killed_M6_no_advance",    # 6
    "exited_domain",           # 7
)
N_OUTCOMES = len(OUTCOME_NAMES)

# 1° × 1° lon/lat grid covering the KRICO domain.
# Domain bounds (from methodology / domain map): 115°W to 40°E, 78°S to 40°S.
LON_MIN, LON_MAX = -115.0, 40.0
LAT_MIN, LAT_MAX = -78.0, -40.0
DLON = DLAT = 1.0


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

COHORT_FILENAME_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


def find_cohort_files(data_dir: Path) -> list[Path]:
    out = []
    for path in sorted(data_dir.iterdir()):
        if COHORT_FILENAME_RE.match(path.name):
            out.append(path)
    return out


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(data_dir: Path) -> xr.Dataset:
    """
    Walk all cohort files and accumulate per-outcome counts on the lon/lat
    grid using final_lon / final_lat.

    Returns a dataset with dimensions (outcome, lat, lon) containing both
    raw counts and within-outcome normalized density.
    """
    # Cell edges (n+1) and centers (n).
    lon_edges = np.arange(LON_MIN, LON_MAX + DLON, DLON)
    lat_edges = np.arange(LAT_MIN, LAT_MAX + DLAT, DLAT)
    lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2.0
    lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2.0
    n_lon = len(lon_centers)
    n_lat = len(lat_centers)

    counts = np.zeros((N_OUTCOMES, n_lat, n_lon), dtype=np.int64)
    n_dropped_nan = 0

    files = find_cohort_files(data_dir)
    if not files:
        raise FileNotFoundError(f"No cohort files found in {data_dir}")

    print(f"Found {len(files)} cohort files in {data_dir}")
    print(f"Grid: lon {LON_MIN}..{LON_MAX} ({n_lon} cells), "
          f"lat {LAT_MIN}..{LAT_MAX} ({n_lat} cells)")

    for n, path in enumerate(files, start=1):
        with xr.open_dataset(path) as ds:
            outcome = ds["outcome"].values
            final_lon = ds["final_lon"].values
            final_lat = ds["final_lat"].values

        # Drop particles with NaN final positions. These can occur from
        # Parcels field-interpolation issues at the very last step even
        # for non-deleted particles; the recruitment pipeline preserves
        # the NaN rather than papering over it.
        valid = ~(np.isnan(final_lon) | np.isnan(final_lat))
        n_dropped_nan += int((~valid).sum())

        # Bin separately per outcome.
        for code in range(N_OUTCOMES):
            mask = valid & (outcome == code)
            if not mask.any():
                continue
            h, _, _ = np.histogram2d(
                final_lat[mask],
                final_lon[mask],
                bins=[lat_edges, lon_edges],
            )
            counts[code] += h.astype(np.int64)

        if n % 200 == 0:
            print(f"  processed {n} / {len(files)}")

    print("Done.")
    if n_dropped_nan > 0:
        print(f"  dropped {n_dropped_nan:,} particles with NaN final position")

    # Normalize within each outcome: density sums to 1 per outcome.
    totals = counts.sum(axis=(1, 2)).astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        density = np.where(
            totals[:, None, None] > 0,
            counts / totals[:, None, None],
            0.0,
        ).astype(np.float32)

    ds_out = xr.Dataset(
        data_vars={
            "counts": (
                ("outcome", "lat", "lon"),
                counts,
                {"long_name": "particle count per outcome per 1° × 1° fate-position cell"},
            ),
            "density": (
                ("outcome", "lat", "lon"),
                density,
                {
                    "long_name": "fraction of particles of given outcome with fate position in this cell",
                    "description": "Sums to 1 across (lat, lon) for each outcome.",
                },
            ),
            "total": (
                ("outcome",),
                totals.astype(np.int64),
                {"long_name": "total particles in each outcome category"},
            ),
        },
        coords={
            "outcome": ("outcome", np.array(OUTCOME_NAMES)),
            "lat": ("lat", lat_centers.astype(np.float32),
                    {"units": "degrees_north", "long_name": "fate latitude (cell center)"}),
            "lon": ("lon", lon_centers.astype(np.float32),
                    {"units": "degrees_east", "long_name": "fate longitude (cell center)"}),
        },
        attrs={
            "title": "F5 destination maps by fate — climatological aggregation",
            "description": (
                "Per-outcome density of fate positions (final_lon / final_lat) "
                "on a 1° × 1° lon/lat grid, aggregated over all 32 spawning years "
                "and all release dates."
            ),
            "grid_lon_min": LON_MIN,
            "grid_lon_max": LON_MAX,
            "grid_lat_min": LAT_MIN,
            "grid_lat_max": LAT_MAX,
            "grid_resolution_deg": DLON,
        },
    )
    return ds_out


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    krico_root = os.environ.get("KRICO_ROOT")
    if not krico_root:
        print("ERROR: KRICO_ROOT environment variable not set.", file=sys.stderr)
        sys.exit(1)

    data_dir = Path(krico_root) / "Post" / "Production" / "recruitment" / "data"
    if not data_dir.is_dir():
        print(f"ERROR: recruitment data dir not found: {data_dir}", file=sys.stderr)
        sys.exit(1)

    ds = aggregate(data_dir)

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "aggregated.nc"

    print(f"Writing {out_path}")
    encoding = {
        var: {"zlib": True, "complevel": 5}
        for var in ds.data_vars if ds[var].dtype.kind != "U"
    }
    ds.to_netcdf(out_path, encoding=encoding)

    print()
    print("Outcome totals:")
    for name, total in zip(OUTCOME_NAMES, ds["total"].values):
        print(f"  {name:28s} {int(total):14,d}")


if __name__ == "__main__":
    main()
