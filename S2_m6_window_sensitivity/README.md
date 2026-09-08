# S2 — Sensitivity of M6 to the advance-detection window

Supporting Information analysis for Paper 1.

## What this answers

Sea-ice advance is detected from April 1 onward, within a tracking period
fixed at 200 days from release. The opening date is common to every cohort;
the closing date is not. The interval over which advance can be detected
therefore lengthens through the release season:

| Release | April 1 at day | Tracking ends | Detection window |
|---|---|---|---|
| Nov 15 | 137 | Jun 3 | 63 days |
| Dec 15 | 107 | Jul 3 | 93 days |
| Jan 15 | 76 | Aug 3 | 124 days |
| Feb 15 | 45 | Sep 3 | 155 days |
| Mar 14 | 18 | Sep 30 | 182 days |

M6 (no winter ice) is defined as no advance detected before the end of
tracking, so it is not measured over a comparable interval across release
dates. A particle drifting where ice advances in August is necessarily M6
for a November cohort and can be a success, M5a or M5b for a March cohort,
from identical physics. The reported seasonal decline in M6 is therefore an
upper bound on the underlying signal.

This analysis measures how much of that decline survives when every cohort
is held to the same detection window.

## Method

For each cohort and each common calendar cutoff, any particle whose advance
event fell after the cutoff is demoted to M6, as it would have been had
tracking stopped there. Only the three outcomes evaluated at the advance
event can move — recruitment success, M5a and M5b — because `travel_time`
in the recruitment output is the advance day for exactly those. M1 is fixed
at release, M4 at the starvation threshold, and M6, censored and domain exit
at the end of tracking.

The sweep runs entirely on the archived per-particle recruitment outputs. It
does not read the raw trajectories, does not use the GLORYS12 fields, and
does not modify the classification: the v2.0.0 dataset and every figure and
table in the paper are unchanged.

## What this establishes, and what it does not

Truncating to a common cutoff removes the bookkeeping advantage that late
releases hold over early ones, but it does so by discarding real advance
events that late-release larvae did experience. The resulting seasonal slope
is biased flat, just as the reported slope is biased steep.

**The sweep brackets the M6 decline from below. It does not measure it.**

The unbiased test would extend early-season cohorts to a common calendar end
rather than truncating late-season ones — tracking a November cohort to
September 30 as a March cohort already is — and requires re-running the
trajectories beyond 200 days. That is out of scope here.

Read the result accordingly. A slope that stays clearly positive across the
sweep indicates a real seasonal relaxation whose steepness is uncertain. A
slope collapsing toward zero at the shortest window would indicate that the
reported decline is largely an artefact of the tracking length.

## Reference case and verification

The September 30 cutoff is the end of tracking for a March 14 release, the
longest window any cohort has, so it truncates nothing. Its row must
reproduce the reported classification exactly: zero demotions, and M6 at
45.8% on November 15 falling to 33.6% on March 14. `compare.py` checks this
rather than assuming it, and reports `[ok]` or `[FAIL]`.

A second, weaker check comes free: the November 15 column is invariant
across cutoffs, because a November cohort's tracking already ends on June 3
and nothing can be removed from it.

## Second test: how much of M6 is beyond the reach of ice

The truncation sweep bounds the decline from below, but it is a blunt
instrument: it demotes real advance events wherever they occurred, including
in regions the ice never reaches. `ice_reachability.py` asks the
complementary question — what share of M6 sits north of the winter sea-ice
edge, and so could not have met ice however long tracking continued.

The edge is the 15% contour of the multi-year mean September concentration
field (NOAA/NSIDC G02202 v6), read from `F1_domain_map/data/aggregated.nc`:
the same field and threshold Figure 1 draws, so the two are consistent by
construction and nothing needs downloading. It is reduced to a northern
ice-edge latitude per 1-degree longitude bin, and a particle counts as
beyond reach if its fate latitude lies north of that line. Three edges are
reported — the two epochs and their northernmost envelope, the envelope
being the most conservative because it counts fewest particles as beyond
reach.

Two approximations, both stated in the script docstring: the fate position
stands in for the September position, which is reasonable because the ACC is
largely zonal so latitude is approximately conserved; and the edge is a
multi-year mean, about which individual years vary.

Particles beyond reach are M6 for a physical reason rather than a
bookkeeping one, and the truncation sweep cannot distinguish them. The
larger their share, the less of M6 is window-sensitive.

```bash
python ice_reachability.py --limit 40
python ice_reachability.py
```

Writes `data/m6_ice_reachability.csv`, one row per cohort, and prints the
share by release month.

## Result

Full 32-year record, 3 848 cohorts, 2 100 265 336 particles. Reproduce with
`compare.py` for the truncation sweep and `ice_reachability.py` for the
reach of winter ice.

**Truncation sweep.** The September 30 cutoff truncates no cohort and
reproduces the reported classification exactly, which `compare.py` checks
rather than assumes. The November 15 column is invariant across cutoffs, as
it must be: a November cohort's tracking already ends on June 3.

| Cutoff | Window | Demoted | M6 Nov 15 | M6 Mar 14 | Slope |
|---|---:|---:|---:|---:|---:|
| Jun 3 | 63 d | 10.03% | 45.77% | 49.76% | **-3.99** |
| Jul 3 | 93 d | 3.89% | 45.77% | 41.83% | 3.94 |
| Aug 3 | 124 d | 0.83% | 45.77% | 36.79% | 8.98 |
| Sep 30 | 182 d | 0.00% | 45.77% | 33.61% | **12.16** |

Slope is Nov 15 minus Mar 14, in percentage points. The seasonal decline in
M6 does not survive equalisation: it narrows monotonically as the common
window shortens, and reverses at 63 days.

Recruitment success behaves differently. Its peak release date moves only
from January 29 to January 21 across the whole bracket, well inside the
46-day interval over which success stays within 10% of its peak. The peak
value falls from 7.12% to 4.70%. **The January optimum is not an artefact of
the tracking window, though its amplitude is sensitive to it.**

**Reach of winter ice.** 88.9% of all M6 particles lie north of the
September sea-ice edge (envelope of the two epochs, the most conservative of
the three edges reported). The share is not uniform across the season, and
the decomposition is what matters:

| Release month | M6 (% released) | beyond reach | within reach |
|---|---:|---:|---:|
| November | 45.79 | 35.76 | 10.03 |
| December | 43.32 | 35.95 | 7.36 |
| January | 38.61 | 35.52 | 3.09 |
| February | 35.58 | 34.37 | 1.21 |
| March | 34.07 | 33.32 | 0.75 |

The beyond-reach component is nearly flat, 35.8% to 33.3%. The within-reach
component falls from 10.0% to 0.75% and carries 79% of the 11.7-point
seasonal decline in M6. Since a particle within reach is scored M6 only
because tracking stopped before the ice arrived, that component is largely a
window effect.

**Conclusion.** Report the magnitude of M6, not its seasonal trend. The
magnitude is robust and physical: most of the category sits beyond the
northern limit of winter sea ice, consistent with the absence of extensive
shelf-slope habitat in 48.3 and 48.6N. The trend is not recoverable from a
fixed 200-day tracking period.

Note that the two tests report the decline on different bases. The
truncation sweep uses release days (Nov 15 to Mar 14, 12.2 points), the
reachability test uses release months (November to March, 11.7 points). Both
are stated where used.

## Requirements

One environment variable:

- `KRICO_POST` — path to the krico-post-production root. Cohort files are
  read from `$KRICO_POST/recruitment/data`.

Unlike `S1_m1_offset_sensitivity/`, this analysis needs neither `KRICO_RUNS`
nor `KRICO_GLORYS12`, and does not import `krico_recruitment`. It reads only
the recruitment outputs, as the `F*/aggregate.py` scripts do, and runs
comfortably on a laptop.

## Running

```bash
export KRICO_POST="/path/to/krico-post-production"

cd S2_m6_window_sensitivity
python sweep.py --limit 40      # quick check on the first 40 cohorts
python sweep.py                 # full record, 3848 cohorts
python compare.py
```

`sweep.py` writes `data/m6_window_sweep.csv`, one row per (cohort, cutoff):
release date, spawning year, season day, cutoff, days from release to
cutoff, particle count, number demoted, and the eight outcome counts after
demotion. The CSV is committed; it is small and stable, like S1's.

To sweep different cutoffs:

```bash
python sweep.py --cutoffs 06-03,07-03,08-03,09-30
```

Cutoffs are given as `MM-DD` and applied within the spawning year. A cutoff
earlier than April 1 is rejected, since it would demote every advance event.
