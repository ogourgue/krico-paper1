"""
F2 aggregation: 30-year (year, season_day) success rate for the phenology curve.

Reads per-particle recruitment outcomes from all cohort files in the
recruitment data directory and aggregates to a (year, season_day) grid:

  - season_day = 0 corresponds to Nov 15
  - season_day = 120 corresponds to Mar 15
  - Feb 29 cohorts in leap years are excluded (~0.2% data loss)
  - Cohorts that do not exist (Nov 31, Feb 30) are missing days, kept as NaN

For each (year, season_day) cell, the script records:
  - total particle count
  - success particle count

Success rate = success / total.

Inputs:  $KRICO_ROOT/Post/Production/recruitment/data/YYYY_MM_DD.nc (~3800 files)
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

# Spawning years span 1994..2025 (32 years).
# A spawning year Y corresponds to the season Nov 15 (Y-1) -> Mar 15 (Y).
SPAWNING_YEARS = np.arange(1994, 2026)
N_YEARS = len(SPAWNING_YEARS)
N_SEASON_DAYS = 121  # Nov 15 -> Mar 15 inclusive

# Outcome flag for "success" in recruitment NetCDFs (CF integer flag).
OUTCOME_SUCCESS = 0


# ---------------------------------------------------------------------------
# Date <-> season-day mapping
# ---------------------------------------------------------------------------

def season_day_index(year: int, month: int, day: int, spawning_year: int) -> int | None:
    """
    Convert a calendar date into a season-day index 0..120.

    A spawning year Y starts on Nov 15 (Y-1) and ends on Mar 15 (Y).

    Returns None for:
      - Feb 29 (leap years): excluded uniformly
      - Dates outside the Nov 15 -> Mar 15 window
      - Dates inconsistent with the given spawning_year

    Parameters
    ----------
    year, month, day : int
        Calendar date of the cohort release.
    spawning_year : int
        Spawning year Y; the season runs Nov 15 (Y-1) to Mar 15 (Y).
    """
    # Skip Feb 29 entirely.
    if month == 2 and day == 29:
        return None

    # November/December must belong to spawning_year - 1.
    if month in (11, 12):
        if year != spawning_year - 1:
            return None
    elif month in (1, 2, 3):
        if year != spawning_year:
            return None
    else:
        return None  # not in the season

    # Build a per-month offset table for non-leap-year season counting.
    # Days from Nov 15 to start of each month, treating Feb as 28 days.
    # season_day = 0 on Nov 15.
    if month == 11:
        if day < 15:
            return None
        idx = day - 15                              # 0..15
    elif month == 12:
        idx = (30 - 15) + day                       # 16..46  (Nov has 30 days)
    elif month == 1:
        idx = (30 - 15) + 31 + day                  # 47..77  (Dec 31 days)
    elif month == 2:
        if day > 28:
            return None                             # Feb 29 already filtered above
        idx = (30 - 15) + 31 + 31 + day             # 78..105 (Jan 31 days)
    elif month == 3:
        if day > 15:
            return None
        idx = (30 - 15) + 31 + 31 + 28 + day        # 106..120 (Feb 28 days)
    else:
        return None

    if idx < 0 or idx >= N_SEASON_DAYS:
        return None
    return idx


# ---------------------------------------------------------------------------
# File discovery and parsing
# ---------------------------------------------------------------------------

# Recruitment files are named YYYY_MM_DD.nc (matches release date).
COHORT_FILENAME_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


def find_cohort_files(data_dir: Path) -> list[tuple[int, int, int, Path]]:
    """
    Return [(year, month, day, path), ...] for all recruitment cohort files.
    """
    out = []
    for path in sorted(data_dir.iterdir()):
        m = COHORT_FILENAME_RE.match(path.name)
        if not m:
            continue
        y, mo, d = (int(g) for g in m.groups())
        out.append((y, mo, d, path))
    return out


def spawning_year_for_date(year: int, month: int) -> int:
    """
    Map a calendar (year, month) to the spawning year it belongs to.
    Nov-Dec belong to spawning year (year + 1); Jan-Mar belong to spawning year (year).
    """
    if month in (11, 12):
        return year + 1
    return year


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(data_dir: Path) -> xr.Dataset:
    """
    Walk all cohort files and accumulate counts on the (year, season_day) grid.
    """
    total = np.zeros((N_YEARS, N_SEASON_DAYS), dtype=np.int64)
    success = np.zeros((N_YEARS, N_SEASON_DAYS), dtype=np.int64)
    has_data = np.zeros((N_YEARS, N_SEASON_DAYS), dtype=bool)

    files = find_cohort_files(data_dir)
    if not files:
        raise FileNotFoundError(f"No cohort files found in {data_dir}")

    print(f"Found {len(files)} cohort files in {data_dir}")

    n_processed = 0
    n_skipped = 0

    for year, month, day, path in files:
        sy = spawning_year_for_date(year, month)
        sy_idx = sy - SPAWNING_YEARS[0]
        if sy_idx < 0 or sy_idx >= N_YEARS:
            n_skipped += 1
            continue

        sd_idx = season_day_index(year, month, day, sy)
        if sd_idx is None:
            n_skipped += 1
            continue

        # Read just the outcome variable; everything else is unused for F2.
        with xr.open_dataset(path) as ds:
            outcome = ds["outcome"].values

        n_total = outcome.size
        n_success = int((outcome == OUTCOME_SUCCESS).sum())

        total[sy_idx, sd_idx] = n_total
        success[sy_idx, sd_idx] = n_success
        has_data[sy_idx, sd_idx] = True
        n_processed += 1

        if n_processed % 200 == 0:
            print(f"  processed {n_processed} / {len(files)}")

    print(f"Done: {n_processed} processed, {n_skipped} skipped (e.g. Feb 29)")

    # Build output dataset.
    success_rate = np.where(has_data, success / np.maximum(total, 1), np.nan)

    ds = xr.Dataset(
        data_vars={
            "total": (("year", "season_day"), total),
            "success": (("year", "season_day"), success),
            "success_rate": (("year", "season_day"), success_rate.astype(np.float32)),
            "has_data": (("year", "season_day"), has_data),
        },
        coords={
            "year": ("year", SPAWNING_YEARS.astype(np.int32),
                     {"long_name": "spawning year",
                      "description": "Spawning year Y: season Nov 15 (Y-1) -> Mar 15 (Y)"}),
            "season_day": ("season_day", np.arange(N_SEASON_DAYS, dtype=np.int32),
                           {"long_name": "season day index",
                            "description": "0 = Nov 15, 120 = Mar 15; Feb 29 excluded"}),
        },
        attrs={
            "title": "F2 phenology curve aggregation",
            "description": "Per-(year, season_day) success counts for KRICO Paper 1 figure F2.",
        },
    )
    return ds


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
    encoding = {var: {"zlib": True, "complevel": 5} for var in ds.data_vars}
    ds.to_netcdf(out_path, encoding=encoding)

    # Quick summary to stdout.
    print()
    print("Summary:")
    print(f"  shape: {ds['total'].shape}")
    print(f"  total cells with data: {int(ds['has_data'].sum())} / {ds['has_data'].size}")
    print(f"  total particles:       {int(ds['total'].sum()):,}")
    print(f"  success particles:     {int(ds['success'].sum()):,}")

    sr = ds["success_rate"].values
    sr_valid = sr[ds["has_data"].values]
    print(f"  success rate: min {np.nanmin(sr_valid):.3f}, "
          f"mean {np.nanmean(sr_valid):.3f}, "
          f"max {np.nanmax(sr_valid):.3f}")


if __name__ == "__main__":
    main()
