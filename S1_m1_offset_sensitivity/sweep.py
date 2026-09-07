"""
Sensitivity of the M1 classification to the descent-ascent offset.

M1 (ice at spawning) is evaluated at the spawning position and date. Spawning
precedes release by the descent-ascent interval, which Thorpe et al. (2019)
give as 23-26 days and which the pipeline applies as a constant. This script
sweeps that constant and reports how much the M1 classification depends on it.

An offset of 0 days samples the sea-ice field on the release date itself.
That case is included as a reference rather than as a candidate value: with
--verify, its classification is compared against the sea-ice concentration
sampled along each trajectory at day 0, which tests the field lookup end to
end.

Whether the two agree exactly depends on the date, for a reason external to
this code. GLORYS12 ice fields are stamped at 12:00 up to 2018-12-18 and at
11:52:30 from 2018-12-19 onward; the velocity, temperature and mixing fields
are stamped at 12:00 throughout. Before the step, particle times land on the
field timestamps, so Parcels' temporal interpolation returns the stored value
and the lookup reproduces it bit for bit. After it, every sample sits 7.5
minutes off the stamp, so Parcels returns a 450/86400 blend of two adjacent
daily fields while the lookup returns the containing day. Adjacent daily ice
fields differ by order 1e-2, so the blend differs by order 1e-4 for most
particles.

Measured over the full record: 1.8e-05 of particles classified differently,
mean absolute concentration difference 1.2e-05. By spawning year the rate is
about 1e-06 through 2018 and 5e-05 to 1.2e-04 from 2019, tracking the
timestamp step rather than anything physical.

Neither behaviour is an error. Taking the field for the day containing the
spawning date is the intended semantics for M1, and Parcels interpolates
correctly for the trajectory-sampled filters. The tolerance below is set
well above the observed rate: it exists to catch a genuine indexing failure,
which would produce percent-level disagreement, not 1e-04.

M1 is computed here by the same functions the pipeline uses, imported from
krico_recruitment, so this analysis cannot drift from the classification it
describes.

Usage
-----
    python sweep.py --runs   $KRICO_ROOT/Runs \
                    --glorys $KRICO_ROOT/Pre/GLORYS12 \
                    --offsets 0 23 24 25 26 \
                    --verify \
                    --out data/m1_offset_sweep.csv

Output is one row per cohort per offset.
"""

import argparse
import datetime
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from netCDF4 import Dataset

from krico_recruitment.filters import evaluate_M1
from krico_recruitment.sea_ice import spawning_sic

COHORT_RE = re.compile(r'^(\d{4})_(\d{2})_(\d{2})\.nc$')

# Fraction of particles whose M1 classification may differ between the
# offset-0 field lookup and the trajectory sample before the difference
# stops being attributable to the timestamp step documented above. The
# observed full-record rate is 1.8e-05, with a worst spawning year of
# 1.2e-04; a real indexing failure would be percent-level.
VERIFY_TOLERANCE = 1e-3


def find_cohorts(runs_root, first_year, last_year, stride):
    """List (release_date, path, spawning_year) for the cohort files."""
    cohorts = []
    for path in sorted(Path(runs_root).glob('KRICO_*/*.nc')):
        m = COHORT_RE.match(path.name)
        if not m:
            continue
        date = datetime.date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
        # A spawning year runs 15 Nov (Y-1) to 15 Mar (Y).
        spawning_year = date.year + 1 if date.month >= 11 else date.year
        if first_year is not None and spawning_year < first_year:
            continue
        if last_year is not None and spawning_year > last_year:
            continue
        cohorts.append((date, path, spawning_year))
    cohorts.sort()
    return cohorts[::stride] if stride > 1 else cohorts


def season_day(date):
    """Day index into the release season, 0 = 15 November."""
    year = date.year if date.month >= 11 else date.year - 1
    return (date - datetime.date(year, 11, 15)).days


def read_release_day(path, want_sic):
    """Release positions, and the trajectory sea-ice sample if requested."""
    with Dataset(path) as nc:
        lon = np.asarray(nc.variables['lon'][:, 0], dtype='float64')
        lat = np.asarray(nc.variables['lat'][:, 0], dtype='float64')
        sic = (np.asarray(nc.variables['sea_ice_area_fraction'][:, 0],
                          dtype='float64') if want_sic else None)
    return lon, lat, sic


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--runs', required=True)
    p.add_argument('--glorys', required=True)
    p.add_argument('--offsets', type=int, nargs='+',
                   default=[0, 23, 24, 25, 26])
    p.add_argument('--first-year', type=int, default=1994)
    p.add_argument('--last-year', type=int, default=None)
    p.add_argument('--stride', type=int, default=1)
    p.add_argument('--verify', action='store_true',
                   help='check the offset-0 lookup against the trajectory '
                        'sample at day 0')
    p.add_argument('--out', default='data/m1_offset_sweep.csv')
    args = p.parse_args()

    cohorts = find_cohorts(args.runs, args.first_year, args.last_year,
                           args.stride)
    if not cohorts:
        sys.exit('No cohort files found under --runs.')
    print(f'{len(cohorts)} cohorts, offsets {args.offsets}\n')

    verify_offset = 0 if (args.verify and 0 in args.offsets) else None
    if args.verify and verify_offset is None:
        print('--verify needs offset 0 in --offsets; verification skipped\n')

    rows = []
    for k, (date, path, spawning_year) in enumerate(cohorts, 1):
        lon, lat, sic_traj = read_release_day(path, verify_offset is not None)

        line = f'[{k}/{len(cohorts)}] {date}'
        for offset in args.offsets:
            try:
                sic = spawning_sic(lon, lat, date, args.glorys, offset)
            except FileNotFoundError:
                print(f'{line}  skip offset {offset}: sea-ice file missing')
                continue

            killed = evaluate_M1(sic)
            rows.append({
                'release_date': date.isoformat(),
                'spawning_year': spawning_year,
                'season_day': season_day(date),
                'offset_days': offset,
                'n_particles': int(killed.size),
                'n_M1': int(np.count_nonzero(killed)),
                'frac_M1': float(np.count_nonzero(killed) / killed.size),
                'n_missing': int(np.count_nonzero(np.isnan(sic))),
            })
            line += f'  {offset}d {rows[-1]["frac_M1"]:.2%}'

            # Regression test against the trajectory sample at day 0. Also
            # record the concentration differences, which separate
            # interpolation error at sharp gradients (a few coastal cells)
            # from a systematic offset in the lookup (a bias everywhere).
            if offset == verify_offset:
                killed_traj = evaluate_M1(sic_traj)
                mismatch = int(np.count_nonzero(killed != killed_traj))
                with np.errstate(invalid='ignore'):
                    diff = np.abs(sic - sic_traj)
                rows[-1]['n_verify_mismatch'] = mismatch
                rows[-1]['verify_max_diff'] = float(np.nanmax(diff))
                rows[-1]['verify_mean_diff'] = float(np.nanmean(diff))
                if mismatch:
                    line += f'  [verify {mismatch}]'

        print(line)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f'\nWritten: {args.out}')

    if 'n_verify_mismatch' in df.columns:
        ref = df[df.offset_days == verify_offset]
        total = int(ref['n_verify_mismatch'].sum())
        n = int(ref['n_particles'].sum())
        rate = total / n
        print('\nVerification against the trajectory sample at day 0:')
        print(f'  {total:,} of {n:,} particles ({rate:.4%}) '
              f'classified differently')
        print(f'  concentration difference: '
              f'mean {ref["verify_mean_diff"].mean():.2e}, '
              f'max {ref["verify_max_diff"].max():.2e}')
        if rate > VERIFY_TOLERANCE:
            print(f'  ABOVE TOLERANCE ({VERIFY_TOLERANCE:.0e}): the field '
                  f'lookup and the trajectory sampling disagree by more than '
                  f'the timestamp step accounts for. Investigate before using '
                  f'the sweep.')
        else:
            print(f'  Within tolerance ({VERIFY_TOLERANCE:.0e}). The residual '
                  f'follows the 2018-12-19 ice-field timestamp step: exact '
                  f'agreement before it, order 1e-04 after.')


if __name__ == '__main__':
    main()
