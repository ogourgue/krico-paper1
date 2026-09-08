"""
S2 sweep: outcome composition under a common advance-detection cutoff.

Advance detection opens on April 1 for every cohort and closes at the end of
tracking, 200 days after release. The opening date is fixed but the closing
date is not, so the interval over which advance can be detected lengthens
through the release season: 63 days for a November 15 release (closing
June 3), 182 days for a March 14 release (closing September 30). M6 is
therefore evaluated over an unequal window across release dates, and its
reported seasonal decline is an upper bound on the underlying signal.

This sweep re-tallies the archived per-particle outcomes under a *common*
calendar cutoff. For each cohort and each cutoff date, any particle whose
advance event fell after the cutoff is demoted to M6, as it would have been
had tracking stopped there. The reported classification is recovered exactly
at the September 30 cutoff, which truncates no cohort.

What this does and does not establish
-------------------------------------
Truncating to a common cutoff removes the bookkeeping advantage that late
releases have over early ones, but it does so by discarding real advance
events that late-release larvae did experience. The resulting seasonal slope
is therefore biased flat, just as the reported slope is biased steep. The
sweep brackets the M6 decline from below; it does not measure it.

The unbiased test would extend early-season cohorts to a common calendar end
rather than truncating late-season ones, which requires re-running the
trajectories beyond 200 days and is out of scope here.

Only the outcomes evaluated at the advance event can be demoted: recruitment
success, M5a (under-developed at advance) and M5b (off-shelf at advance).
M1 is fixed at release, M4 at the starvation threshold, and M6, censored and
domain exit at the end of tracking, so none of them move.

Inputs:  $KRICO_POST/recruitment/data/YYYY_MM_DD.nc  (3848 cohort files)
Output:  data/m6_window_sweep.csv

Usage:
    export KRICO_POST=/path/to/krico-post-production
    python sweep.py
    python sweep.py --cutoffs 06-03,07-03,08-03,09-30
    python sweep.py --limit 40          # quick check on the first 40 cohorts
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
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

CODE = {name: i for i, name in enumerate(OUTCOME_NAMES)}

# Outcomes determined at the sea-ice advance event. Only these can be demoted.
ADVANCE_OUTCOMES = (
    CODE["success"],
    CODE["killed_M5_no_FIV"],
    CODE["killed_M5_not_on_shelf"],
)

# Common cutoffs to sweep, as (month, day) in the spawning year.
#   06-03  end of tracking for a Nov 15 release: the shortest window available
#   07-03  end of tracking for a Dec 15 release
#   08-03  end of tracking for a Jan 15 release
#   09-30  end of tracking for a Mar 14 release: truncates nothing, so this
#          row must reproduce the reported classification exactly
DEFAULT_CUTOFFS = ("06-03", "07-03", "08-03", "09-30")

COHORT_FILENAME_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


# ---------------------------------------------------------------------------
# Date handling
# ---------------------------------------------------------------------------

def spawning_year_for_date(year: int, month: int) -> int:
    """Map a calendar (year, month) to the spawning year it belongs to."""
    if month in (11, 12):
        return year + 1
    return year


def season_day_index(release: dt.date, spawning_year: int) -> int | None:
    """
    Season-day index 0..119 (Nov 15 -> Mar 14), or None if out of range.

    Feb 29 is retained here, unlike F2 and F3: this sweep reports per cohort
    and does not build a uniform (year, season_day) grid, so the leap-day
    cohorts carry a season_day of None rather than being dropped.
    """
    y, m, d = release.year, release.month, release.day
    if m in (11, 12):
        if y != spawning_year - 1:
            return None
    elif m in (1, 2, 3):
        if y != spawning_year:
            return None
    else:
        return None

    if m == 2 and d == 29:
        return None

    if m == 11:
        if d < 15:
            return None
        return d - 15
    if m == 12:
        return 15 + d
    if m == 1:
        return 15 + 31 + d
    if m == 2:
        return 15 + 31 + 31 + d
    if m == 3:
        if d > 14:
            return None
        return 15 + 31 + 31 + 28 + d
    return None


def parse_cutoffs(spec: str) -> list[tuple[int, int]]:
    """Parse 'MM-DD,MM-DD,...' into a list of (month, day) pairs."""
    out = []
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        try:
            month, day = (int(p) for p in item.split("-"))
            dt.date(2001, month, day)  # validate; 2001 is non-leap
        except ValueError:
            raise SystemExit(f"bad cutoff {item!r}: expected MM-DD, e.g. 06-03")
        if (month, day) < (4, 1):
            raise SystemExit(
                f"cutoff {item!r} precedes the April 1 detection gate, which "
                "would demote every advance event and is not a meaningful test"
            )
        out.append((month, day))
    if not out:
        raise SystemExit("no cutoffs given")
    return out


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------

def find_cohort_files(data_dir: Path) -> list[tuple[dt.date, Path]]:
    out = []
    for path in sorted(data_dir.iterdir()):
        m = COHORT_FILENAME_RE.match(path.name)
        if not m:
            continue
        y, mo, d = (int(g) for g in m.groups())
        out.append((dt.date(y, mo, d), path))
    return out


# ---------------------------------------------------------------------------
# Per-cohort re-tally
# ---------------------------------------------------------------------------

def retally(outcome: np.ndarray,
            travel_time: np.ndarray,
            release: dt.date,
            cutoff: dt.date) -> tuple[np.ndarray, int]:
    """
    Re-tally one cohort under a common advance-detection cutoff.

    Parameters
    ----------
    outcome : ndarray, int
        Per-particle outcome flag, 0..7.
    travel_time : ndarray, int
        Days from release to the moment the fate was determined. For the
        three advance-evaluated outcomes this is the advance day.
    release : datetime.date
        Release date of the cohort.
    cutoff : datetime.date
        Common calendar cutoff. Advance events after this date are demoted.

    Returns
    -------
    counts : ndarray of shape (8,), int64
        Outcome counts after demotion.
    n_demoted : int
        Number of particles moved into M6.
    """
    days_to_cutoff = (cutoff - release).days

    # Particles whose fate was set at the advance event, after the cutoff.
    at_advance = np.isin(outcome, ADVANCE_OUTCOMES)
    demote = at_advance & (travel_time > days_to_cutoff)

    adjusted = outcome.copy()
    adjusted[demote] = CODE["killed_M6_no_advance"]

    counts = np.bincount(adjusted, minlength=N_OUTCOMES).astype(np.int64)
    return counts, int(demote.sum())


# ---------------------------------------------------------------------------
# Sweep
# ---------------------------------------------------------------------------

def run_sweep(data_dir: Path, cutoffs: list[tuple[int, int]],
              out_path: Path, limit: int | None = None) -> None:
    files = find_cohort_files(data_dir)
    if not files:
        raise FileNotFoundError(f"No cohort files found in {data_dir}")
    if limit is not None:
        files = files[:limit]

    print(f"Found {len(files)} cohort files in {data_dir}")
    print(f"Cutoffs: {', '.join(f'{m:02d}-{d:02d}' for m, d in cutoffs)}")
    print()

    out_path.parent.mkdir(parents=True, exist_ok=True)

    header = (["release_date", "spawning_year", "season_day",
               "cutoff", "cutoff_date", "days_to_cutoff",
               "n_particles", "n_demoted"]
              + [f"n_{name}" for name in OUTCOME_NAMES])

    n_processed = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)

        for release, path in files:
            sy = spawning_year_for_date(release.year, release.month)
            sd = season_day_index(release, sy)

            with xr.open_dataset(path) as ds:
                outcome = ds["outcome"].values.astype(np.int64)
                travel_time = ds["travel_time"].values.astype(np.int64)

            if outcome.shape != travel_time.shape:
                raise ValueError(f"{path.name}: outcome and travel_time differ in shape")

            for month, day in cutoffs:
                cutoff_date = dt.date(sy, month, day)
                counts, n_demoted = retally(outcome, travel_time,
                                            release, cutoff_date)
                writer.writerow(
                    [release.isoformat(), sy,
                     "" if sd is None else sd,
                     f"{month:02d}-{day:02d}", cutoff_date.isoformat(),
                     (cutoff_date - release).days,
                     int(outcome.size), n_demoted]
                    + [int(c) for c in counts]
                )

            n_processed += 1
            if n_processed % 200 == 0:
                print(f"  processed {n_processed} / {len(files)}")

    print(f"Done: {n_processed} cohorts -> {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Re-tally KRICO outcomes under a common advance-detection cutoff."
    )
    parser.add_argument(
        "--post", type=str, default=os.environ.get("KRICO_POST"),
        help="krico-post-production root (default: $KRICO_POST). "
             "Cohort files are read from <post>/recruitment/data.",
    )
    parser.add_argument(
        "--cutoffs", type=str, default=",".join(DEFAULT_CUTOFFS),
        help=f"comma-separated MM-DD cutoffs (default: {','.join(DEFAULT_CUTOFFS)})",
    )
    parser.add_argument(
        "--out", type=str, default="data/m6_window_sweep.csv",
        help="output CSV path (default: data/m6_window_sweep.csv)",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="process only the first N cohorts (for a quick check)",
    )
    args = parser.parse_args()

    if not args.post:
        sys.exit("KRICO_POST is not set and --post was not given.")

    data_dir = Path(args.post) / "recruitment" / "data"
    if not data_dir.is_dir():
        sys.exit(f"Not a directory: {data_dir}")

    run_sweep(data_dir, parse_cutoffs(args.cutoffs),
              Path(args.out), args.limit)


if __name__ == "__main__":
    main()
