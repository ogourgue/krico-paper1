"""
F2 aggregation: 32-year (year, season_day) success rate for the phenology curve.

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


def season_day_to_label(day: int) -> str:
    """Calendar label for a season day index (e.g., 0 -> 'Nov 15', 47 -> 'Jan 1').

    Mirrors the helper in plot.py — kept local here so the aggregate script
    is self-contained for printing summary statistics.
    """
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
# Summary statistics for the manuscript Results section
# ---------------------------------------------------------------------------

def print_summary_stats(ds: xr.Dataset) -> None:
    """
    Print extended summary statistics used in the manuscript Results section.

    Computes climatological mean curve, peak, season endpoints, per-year
    peak distribution, per-year season-mean variability, and the width of
    the optimal release window. Suppresses the all-NaN warnings that come
    from cohorts missing across all years (e.g. Nov 31, Feb 30).
    """
    import warnings

    sr_pct = ds["success_rate"].values * 100.0   # (year, season_day), in %
    has_data = ds["has_data"].values
    sr_pct = np.where(has_data, sr_pct, np.nan)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        # Climatological curve (mean across years for each season_day).
        mean_per_day = np.nanmean(sr_pct, axis=0)
        p5_per_day = np.nanpercentile(sr_pct, 5, axis=0)
        p95_per_day = np.nanpercentile(sr_pct, 95, axis=0)
        std_per_day = np.nanstd(sr_pct, axis=0)

        # Climatological peak.
        peak_idx = int(np.nanargmax(mean_per_day))
        peak_mean = float(mean_per_day[peak_idx])
        peak_p5 = float(p5_per_day[peak_idx])
        peak_p95 = float(p95_per_day[peak_idx])
        peak_std = float(std_per_day[peak_idx])
        peak_label = season_day_to_label(peak_idx)

        print()
        print("=" * 60)
        print("Climatological mean curve (32-year mean by season day):")
        print("=" * 60)
        print(f"  Peak:  season_day {peak_idx} ({peak_label})  "
              f"mean {peak_mean:.2f}%  envelope [{peak_p5:.2f}%, {peak_p95:.2f}%]  "
              f"std {peak_std:.2f}%")

        # Season endpoints — handle missing cohorts (e.g. Nov 31 if mislabeled
        # as season_day 16) by using the first/last sd with non-NaN mean.
        valid_sd = np.where(~np.isnan(mean_per_day))[0]
        first_sd, last_sd = int(valid_sd[0]), int(valid_sd[-1])
        for label, sd in [("First (Nov 15)", first_sd),
                          ("Last  (Mar 14)", last_sd)]:
            print(f"  {label}: sd {sd:3d} ({season_day_to_label(sd)})  "
                  f"mean {mean_per_day[sd]:.2f}%  "
                  f"envelope [{p5_per_day[sd]:.2f}%, {p95_per_day[sd]:.2f}%]")

        # Width of the optimal release window.
        for frac in (0.90, 0.75, 0.50):
            threshold = frac * peak_mean
            above = mean_per_day >= threshold
            if above.any():
                first = int(np.argmax(above))
                # last True index, scanning from the right.
                last = len(above) - 1 - int(np.argmax(above[::-1]))
                width = last - first + 1
                print(f"  Window mean >= {frac*100:.0f}% of peak ({threshold:.2f}%): "
                      f"sd {first}-{last} ({season_day_to_label(first)} to "
                      f"{season_day_to_label(last)}), {width} days")

        # Per-year peak day & per-year peak value.
        # We mask all-NaN year rows out of nanargmax to avoid errors.
        years_with_any_data = np.any(~np.isnan(sr_pct), axis=1)
        peak_day_per_year = np.full(N_YEARS, -1, dtype=int)
        peak_val_per_year = np.full(N_YEARS, np.nan)
        for i in range(N_YEARS):
            if years_with_any_data[i]:
                peak_day_per_year[i] = int(np.nanargmax(sr_pct[i]))
                peak_val_per_year[i] = float(np.nanmax(sr_pct[i]))

        valid_years = years_with_any_data
        pd_valid = peak_day_per_year[valid_years]
        pv_valid = peak_val_per_year[valid_years]

        print()
        print("=" * 60)
        print(f"Per-year peak day distribution ({valid_years.sum()} years):")
        print("=" * 60)
        pd_min, pd_max = int(np.min(pd_valid)), int(np.max(pd_valid))
        pd_p25, pd_p50, pd_p75 = np.percentile(pd_valid, [25, 50, 75])
        print(f"  earliest peak: sd {pd_min} ({season_day_to_label(pd_min)})  "
              f"in year {SPAWNING_YEARS[np.where(valid_years)[0][np.argmin(pd_valid)]]}")
        print(f"  latest peak:   sd {pd_max} ({season_day_to_label(pd_max)})  "
              f"in year {SPAWNING_YEARS[np.where(valid_years)[0][np.argmax(pd_valid)]]}")
        print(f"  P25 peak:      sd {int(pd_p25)} ({season_day_to_label(int(pd_p25))})")
        print(f"  P50 peak:      sd {int(pd_p50)} ({season_day_to_label(int(pd_p50))})")
        print(f"  P75 peak:      sd {int(pd_p75)} ({season_day_to_label(int(pd_p75))})")
        print(f"  inter-quartile range: {int(pd_p75) - int(pd_p25)} days")

        print()
        print("=" * 60)
        print("Per-year peak value (best release date per year):")
        print("=" * 60)
        worst_year_idx = np.where(valid_years)[0][int(np.argmin(pv_valid))]
        best_year_idx = np.where(valid_years)[0][int(np.argmax(pv_valid))]
        print(f"  min:  {np.min(pv_valid):.2f}% in year {SPAWNING_YEARS[worst_year_idx]}")
        print(f"  max:  {np.max(pv_valid):.2f}% in year {SPAWNING_YEARS[best_year_idx]}")
        print(f"  mean: {np.mean(pv_valid):.2f}%")
        print(f"  std:  {np.std(pv_valid):.2f}%")

        # Per-year season-mean recruitment success.
        year_means = np.nanmean(sr_pct, axis=1)
        ym_valid = year_means[valid_years]
        worst_seasonal = np.where(valid_years)[0][int(np.argmin(ym_valid))]
        best_seasonal = np.where(valid_years)[0][int(np.argmax(ym_valid))]

        print()
        print("=" * 60)
        print("Per-year season-mean recruitment success:")
        print("=" * 60)
        print(f"  min:  {np.min(ym_valid):.2f}% in year {SPAWNING_YEARS[worst_seasonal]}")
        print(f"  max:  {np.max(ym_valid):.2f}% in year {SPAWNING_YEARS[best_seasonal]}")
        print(f"  mean: {np.mean(ym_valid):.2f}%")
        print(f"  std:  {np.std(ym_valid):.2f}%")
        print(f"  best/worst ratio: {np.max(ym_valid) / max(np.min(ym_valid), 1e-9):.2f}x")


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

    # Extended statistics for the Results section.
    print_summary_stats(ds)


if __name__ == "__main__":
    main()