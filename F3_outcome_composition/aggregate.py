"""
F3 aggregation: per-(year, season_day) counts of all 8 recruitment outcomes.

Reads per-particle recruitment outcomes from all cohort files in the
recruitment data directory and aggregates to a (year, season_day, outcome)
grid:

  - season_day = 0 corresponds to Nov 15
  - season_day = 120 corresponds to Mar 15
  - Feb 29 cohorts in leap years are excluded
  - Cohorts that do not exist (Nov 31, Feb 30) are kept as zero counts

For each (year, season_day) cell, the script records counts for each of the
8 outcome categories (success, censored, killed_M1, killed_M4, killed_M5_no_FIV,
killed_M5_not_on_shelf, killed_M6_no_advance, exited_domain) plus a total.

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

# Spawning years span 1994..2025 (32 years).
# A spawning year Y corresponds to the season Nov 15 (Y-1) -> Mar 15 (Y).
SPAWNING_YEARS = np.arange(1994, 2026)
N_YEARS = len(SPAWNING_YEARS)
N_SEASON_DAYS = 121  # Nov 15 -> Mar 15 inclusive

# Outcome integer flags from the recruitment pipeline (CF flag values).
# Order matches outcome.py in krico-post-production.
OUTCOME_NAMES = (
    "success",
    "censored",
    "killed_M1",
    "killed_M4",
    "killed_M5_no_FIV",
    "killed_M5_not_on_shelf",
    "killed_M6_no_advance",
    "exited_domain",
)
N_OUTCOMES = len(OUTCOME_NAMES)


# ---------------------------------------------------------------------------
# Date <-> season-day mapping (identical to F2; could be factored out later)
# ---------------------------------------------------------------------------

def season_day_index(year: int, month: int, day: int, spawning_year: int) -> int | None:
    """Convert a calendar date into a season-day index 0..120, or None if out of range."""
    if month == 2 and day == 29:
        return None

    if month in (11, 12):
        if year != spawning_year - 1:
            return None
    elif month in (1, 2, 3):
        if year != spawning_year:
            return None
    else:
        return None

    if month == 11:
        if day < 15:
            return None
        idx = day - 15
    elif month == 12:
        idx = (30 - 15) + day
    elif month == 1:
        idx = (30 - 15) + 31 + day
    elif month == 2:
        if day > 28:
            return None
        idx = (30 - 15) + 31 + 31 + day
    elif month == 3:
        if day > 15:
            return None
        idx = (30 - 15) + 31 + 31 + 28 + day
    else:
        return None

    if idx < 0 or idx >= N_SEASON_DAYS:
        return None
    return idx


def spawning_year_for_date(year: int, month: int) -> int:
    """Map a calendar (year, month) to the spawning year it belongs to."""
    if month in (11, 12):
        return year + 1
    return year


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

COHORT_FILENAME_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


def find_cohort_files(data_dir: Path) -> list[tuple[int, int, int, Path]]:
    out = []
    for path in sorted(data_dir.iterdir()):
        m = COHORT_FILENAME_RE.match(path.name)
        if not m:
            continue
        y, mo, d = (int(g) for g in m.groups())
        out.append((y, mo, d, path))
    return out


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def aggregate(data_dir: Path) -> xr.Dataset:
    """Walk all cohort files and accumulate per-outcome counts on the grid."""
    counts = np.zeros((N_YEARS, N_SEASON_DAYS, N_OUTCOMES), dtype=np.int64)
    total = np.zeros((N_YEARS, N_SEASON_DAYS), dtype=np.int64)
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

        with xr.open_dataset(path) as ds:
            outcome = ds["outcome"].values

        # Count occurrences of each outcome code (0..7).
        binc = np.bincount(outcome, minlength=N_OUTCOMES)
        counts[sy_idx, sd_idx, :] = binc
        total[sy_idx, sd_idx] = outcome.size
        has_data[sy_idx, sd_idx] = True
        n_processed += 1

        if n_processed % 200 == 0:
            print(f"  processed {n_processed} / {len(files)}")

    print(f"Done: {n_processed} processed, {n_skipped} skipped")

    ds_out = xr.Dataset(
        data_vars={
            "counts": (("year", "season_day", "outcome"), counts,
                       {"long_name": "number of particles in each outcome category"}),
            "total": (("year", "season_day"), total,
                      {"long_name": "total number of particles released"}),
            "has_data": (("year", "season_day"), has_data),
        },
        coords={
            "year": ("year", SPAWNING_YEARS.astype(np.int32),
                     {"long_name": "spawning year"}),
            "season_day": ("season_day", np.arange(N_SEASON_DAYS, dtype=np.int32),
                           {"long_name": "season day index",
                            "description": "0 = Nov 15, 120 = Mar 15; Feb 29 excluded"}),
            "outcome": ("outcome", np.array(OUTCOME_NAMES)),
        },
        attrs={
            "title": "F3 outcome composition aggregation",
            "description": "Per-(year, season_day) outcome counts for KRICO Paper 1 figure F3.",
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
    encoding = {var: {"zlib": True, "complevel": 5}
                for var in ds.data_vars if ds[var].dtype.kind != "U"}
    ds.to_netcdf(out_path, encoding=encoding)

    # Quick summary to stdout.
    print()
    print("Summary (overall fractions across the full 30-year dataset):")
    total_counts = ds["counts"].sum(dim=("year", "season_day")).values
    grand_total = int(total_counts.sum())
    for name, c in zip(OUTCOME_NAMES, total_counts):
        pct = 100.0 * int(c) / grand_total if grand_total else 0.0
        print(f"  {name:28s} {int(c):14,d}  ({pct:6.3f}%)")
    print(f"  {'TOTAL':28s} {grand_total:14,d}")


if __name__ == "__main__":
    main()
