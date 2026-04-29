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


def season_day_to_label(day: int) -> str:
    """Calendar label for a season day index (e.g., 0 -> 'Nov 15')."""
    if day <= 15:
        return f"Nov {15 + day}"
    if day <= 46:
        return f"Dec {day - 15}"
    if day <= 77:
        return f"Jan {day - 46}"
    if day <= 105:
        return f"Feb {day - 77}"
    return f"Mar {day - 105}"


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
# Summary statistics for the manuscript Results section
# ---------------------------------------------------------------------------

def print_summary_stats(ds: xr.Dataset) -> None:
    """
    Print extended summary statistics used in the manuscript Results
    section for F3.

    Computes:
      - Per-outcome climatological fraction at season start (Nov 15),
        climatological success peak (Jan 19, sd 65), and season end (Mar 14)
      - For each mortality category, the season day at which its
        climatological mean share is highest, and the value at that peak
      - The M1↔M6 crossover day (where M6 share first exceeds M1 share in
        the climatological mean)
      - At each anchor day, the total mortality fraction and the success
        fraction (with success fraction ≈ value from F2 at that day)

    Censored is folded into M6 (consistent with F3 plot.py).
    """
    import warnings

    counts = ds["counts"].values.copy()                  # (year, season_day, outcome)
    total = ds["total"].values                            # (year, season_day)
    has_data = ds["has_data"].values
    outcome_names = list(ds["outcome"].values.astype(str))

    # Fold censored into killed_M6_no_advance.
    censored_idx = outcome_names.index("censored")
    m6_idx = outcome_names.index("killed_M6_no_advance")
    counts[:, :, m6_idx] += counts[:, :, censored_idx]
    counts[:, :, censored_idx] = 0

    # Per-(year, season_day, outcome) fractions, in percent.
    with np.errstate(invalid="ignore", divide="ignore"):
        frac_pct = np.where(
            has_data[:, :, None],
            100.0 * counts / np.maximum(total[:, :, None], 1),
            np.nan,
        )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        # Climatological mean over years -> (season_day, outcome).
        mean_frac = np.nanmean(frac_pct, axis=0)

    # ---- Anchor days: Nov 15 (sd 0), success peak (sd 65, Jan 19), Mar 14 (sd 119)
    anchor_days = [(0, "Nov 15 (season start)"),
                   (65, "Jan 19 (climatological peak success)"),
                   (119, "Mar 14 (season end)")]

    print()
    print("=" * 72)
    print("F3 climatological breakdown by season day (mean over 32 years, in %)")
    print("=" * 72)
    print(f"{'Outcome':28s} {'Nov 15':>8s} {'Jan 19':>8s} {'Mar 14':>8s}")
    print("-" * 72)

    # Display order: top of stack to bottom of stack (success at top).
    display_order = [
        "success",
        "killed_M6_no_advance",
        "killed_M5_not_on_shelf",
        "killed_M5_no_FIV",
        "killed_M4",
        "killed_M1",
        "exited_domain",
    ]
    for name in display_order:
        oi = outcome_names.index(name)
        v0 = mean_frac[anchor_days[0][0], oi]
        v1 = mean_frac[anchor_days[1][0], oi]
        v2 = mean_frac[anchor_days[2][0], oi]
        print(f"{name:28s} {v0:8.2f} {v1:8.2f} {v2:8.2f}")

    # Total mortality (all killed_M*) plus exited_domain (modeling-domain
    # limitation, not a biological outcome -- separate row).
    mortality_indices = [outcome_names.index(n) for n in (
        "killed_M1", "killed_M4", "killed_M5_no_FIV",
        "killed_M5_not_on_shelf", "killed_M6_no_advance",
    )]
    total_mortality = mean_frac[:, mortality_indices].sum(axis=1)
    print("-" * 72)
    for label, vec in [
        ("Total mortality (M1+M4+M5a+M5b+M6)",
         [total_mortality[d] for d, _ in anchor_days]),
    ]:
        print(f"{label:28s} {vec[0]:8.2f} {vec[1]:8.2f} {vec[2]:8.2f}")

    # ---- Per-mortality-category climatological peak day & value
    print()
    print("=" * 72)
    print("Per-category climatological peak (over season days):")
    print("=" * 72)
    print(f"{'Outcome':28s} {'Peak day':>20s} {'Peak %':>8s}")
    print("-" * 72)
    for name in ("success",
                 "killed_M1", "killed_M4",
                 "killed_M5_no_FIV", "killed_M5_not_on_shelf",
                 "killed_M6_no_advance"):
        oi = outcome_names.index(name)
        col = mean_frac[:, oi]
        # Mask all-NaN columns (cohorts missing across all years).
        if np.all(np.isnan(col)):
            print(f"{name:28s} {'(no data)':>20s} {'-':>8s}")
            continue
        peak_sd = int(np.nanargmax(col))
        peak_val = float(np.nanmax(col))
        print(f"{name:28s} {f'sd {peak_sd} ({season_day_to_label(peak_sd)})':>20s} "
              f"{peak_val:8.2f}")

    # ---- M1 -> M6 crossover (climatological)
    print()
    print("=" * 72)
    print("M1 ↔ M6 crossover (climatological mean):")
    print("=" * 72)
    m1 = mean_frac[:, outcome_names.index("killed_M1")]
    m6 = mean_frac[:, outcome_names.index("killed_M6_no_advance")]
    diff = m6 - m1                           # positive when M6 > M1
    sign_changes = np.where(np.diff(np.signbit(diff).astype(int)))[0]
    if sign_changes.size > 0:
        # The first sign change (M6 starts exceeding M1).
        cross_sd = int(sign_changes[0]) + 1   # +1 because diff index k is between sd k and k+1
        print(f"  First season day where M6 share > M1 share: "
              f"sd {cross_sd} ({season_day_to_label(cross_sd)})")
        print(f"    M1 at sd {cross_sd-1}: {m1[cross_sd-1]:.2f}%   "
              f"M6 at sd {cross_sd-1}: {m6[cross_sd-1]:.2f}%")
        print(f"    M1 at sd {cross_sd}:   {m1[cross_sd]:.2f}%   "
              f"M6 at sd {cross_sd}:   {m6[cross_sd]:.2f}%")
    else:
        print("  No sign change detected — one category dominates throughout.")

    # ---- M1 reduction from season start to climatological peak
    print()
    print("=" * 72)
    print("M1 reduction from Nov 15 to Jan 19 (climatological mean):")
    print("=" * 72)
    m1_nov15 = m1[0]
    m1_jan19 = m1[65]
    print(f"  Nov 15: M1 = {m1_nov15:.2f}%")
    print(f"  Jan 19: M1 = {m1_jan19:.2f}%")
    print(f"  Reduction: {m1_nov15 - m1_jan19:+.2f} percentage points "
          f"({100 * (1 - m1_jan19 / max(m1_nov15, 1e-9)):.1f}% relative reduction)")

    # ---- M6 increase from climatological peak to season end
    print()
    print("=" * 72)
    print("M6 increase from Jan 19 to Mar 14 (climatological mean):")
    print("=" * 72)
    m6_jan19 = m6[65]
    m6_mar14 = m6[119]
    print(f"  Jan 19: M6 = {m6_jan19:.2f}%")
    print(f"  Mar 14: M6 = {m6_mar14:.2f}%")
    print(f"  Increase: {m6_mar14 - m6_jan19:+.2f} percentage points")


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

    # Extended statistics for the Results section.
    print_summary_stats(ds)


if __name__ == "__main__":
    main()