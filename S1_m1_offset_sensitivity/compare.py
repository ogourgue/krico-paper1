"""
Summarize the M1 offset sweep for the Supporting Information.

Reads the long-format CSV written by sweep.py and reports how much the M1
classification depends on the descent-ascent offset. The reference case,
offset 0, samples the sea-ice field on the release date and is reported
separately from the 23-26 day range that Thorpe et al. (2019) give.

Usage
-----
    python compare.py --sweep data/m1_offset_sweep.csv
"""

import argparse

import pandas as pd


def weighted(frame, column='n_M1'):
    """Particle-weighted fraction."""
    return frame[column].sum() / frame['n_particles'].sum()


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sweep', default='data/m1_offset_sweep.csv')
    args = p.parse_args()

    df = pd.read_csv(args.sweep)
    df['month'] = pd.to_datetime(df.release_date).dt.strftime('%m %b')

    offsets = sorted(int(o) for o in df.offset_days.unique())
    # Offset 0 is the reference; the descent-ascent range is everything else.
    range_offsets = [o for o in offsets if o > 0]

    print(f'offsets: {offsets}')
    print(f'cohorts: {df.release_date.nunique()}')
    print(f'spawning years: {df.spawning_year.min()}-{df.spawning_year.max()}')
    print(f'particles per offset: '
          f'{int(df[df.offset_days == offsets[0]].n_particles.sum()):,}\n')

    # --- overall -------------------------------------------------------
    print('M1, particle-weighted')
    for o in offsets:
        tag = '  (reference: release date)' if o == 0 else ''
        print(f'  {o:>2}d  {weighted(df[df.offset_days == o]):>7.2%}{tag}')
    lo = weighted(df[df.offset_days == min(range_offsets)])
    hi = weighted(df[df.offset_days == max(range_offsets)])
    print(f'\n  spread across {min(range_offsets)}-{max(range_offsets)}d: '
          f'{abs(hi - lo):.2%}')

    # --- by month ------------------------------------------------------
    print('\nM1 by month of release')
    by_month = df.pivot_table(index='month', columns='offset_days',
                              values=['n_M1', 'n_particles'], aggfunc='sum')
    table = pd.DataFrame({o: by_month[('n_M1', o)] / by_month[('n_particles', o)]
                          for o in offsets})
    table['spread'] = (table[range_offsets].max(axis=1)
                       - table[range_offsets].min(axis=1))
    print(table.map('{:.2%}'.format))

    # --- by spawning year ----------------------------------------------
    by_year = df.pivot_table(index='spawning_year', columns='offset_days',
                             values=['n_M1', 'n_particles'], aggfunc='sum')
    yearly = pd.DataFrame({o: by_year[('n_M1', o)] / by_year[('n_particles', o)]
                           for o in offsets})
    spread = (yearly[range_offsets].max(axis=1)
              - yearly[range_offsets].min(axis=1))
    print('\nM1 by spawning year, spread across the descent-ascent range')
    print(f'  mean {spread.mean():.2%}   '
          f'max {spread.max():.2%} in {spread.idxmax()}')
    above = list(spread.index[spread > 0.01])
    print(f'  years above 1%: {above if above else "none"}')

    # --- per cohort ----------------------------------------------------
    wide = df[df.offset_days.isin(range_offsets)].pivot(
        index='release_date', columns='offset_days', values='frac_M1')
    cohort_spread = wide.max(axis=1) - wide.min(axis=1)
    print('\nper-cohort spread across the descent-ascent range')
    print(f'  mean {cohort_spread.mean():.2%}   '
          f'median {cohort_spread.median():.2%}   '
          f'max {cohort_spread.max():.2%}')
    print(f'  cohorts above 5%: {(cohort_spread > 0.05).sum()} of '
          f'{len(cohort_spread)} ({(cohort_spread > 0.05).mean():.1%})')
    print('  largest:')
    for date, s in cohort_spread.nlargest(5).items():
        values = '  '.join(f'{o}d {wide.loc[date, o]:.1%}'
                           for o in range_offsets)
        print(f'    {date}  {s:.2%}   {values}')

    if 'n_verify_mismatch' in df.columns:
        ref = df[df.offset_days == 0]
        total = int(ref.n_verify_mismatch.sum())
        n = int(ref.n_particles.sum())
        print(f'\nOffset-0 verification against the trajectory sample: '
              f'{total:,} of {n:,} ({total / n:.3%}) differ.')


if __name__ == '__main__':
    main()
