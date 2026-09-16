# S3 — Distance travelled between release and fate

Supporting analysis for the Discussion.

## What this answers

The Discussion opens by noting that advection displaces particles a
substantial distance between release and sea-ice advance, and that it could
therefore have displaced or eliminated the January optimum. That claim was
stated qualitatively. This analysis measures it.

Two quantities, both read from the archived per-particle outputs:

- **Displacement** — great-circle distance from the release position to the
  fate position. Net transport: where the particle ended up relative to where
  it started.
- **Path length** — along-track distance actually covered, from the archived
  `trajectory_path_length` variable, which the recruitment pipeline
  accumulates along each trajectory.

Their ratio is a straightness index. A value of 1 would be a straight line;
lower values mean meandering, recirculation, or reversals along the way.

## Population

The headline figure pools **recruitment success, M5a and M5b** — the three
outcomes evaluated at the sea-ice advance event, which is the moment the
Discussion sentence refers to.

**M1 is excluded.** Its fate position is the release position by
construction, so its displacement is identically zero. Including it would
pull any pooled statistic toward zero without meaning, and M1 is roughly 28%
of all particles.

The remaining outcomes are reported separately. M6 is the interesting one:
those particles are the export cases, drifting in deep water for the full 200
days, and should travel furthest.

## Method

Displacement is computed with the haversine formula on the WGS84 mean radius
(6371 km). Distances are accumulated into 1 km histograms per outcome rather
than held in memory, so percentiles are exact to within one bin across all
2.1 billion particles. Means are accumulated exactly.

Reads only `$KRICO_POST/recruitment/data`. No trajectories, no GLORYS12
fields, no Parcels. Runs on a laptop.

## Requirements

- `KRICO_POST` — path to the krico-post-production root.

As with `S2_m6_window_sensitivity/`, and unlike `S1_m1_offset_sensitivity/`,
this needs neither `KRICO_RUNS` nor `KRICO_GLORYS12`.

## Running

```bash
export KRICO_POST="/path/to/krico-post-production"

cd S3_displacement
python displacement.py --limit 40      # quick check on the first 40 cohorts
python displacement.py                 # full record, 3848 cohorts
```

Writes `data/displacement.csv`: one row per outcome plus the pooled
at-advance population, with mean and 5th, 25th, 50th, 75th and 95th
percentiles of both distances, and the straightness index.

## Result

Full 32-year record, 3 848 cohorts, 2 100 265 336 particles.

| Population | n | Displacement (km) | | Path (km) | Straightness |
|---|---:|---:|---|---:|---:|
| | | median | IQR | median | |
| **at advance** (success + M5a + M5b) | 523 773 649 | **224** | 112–414 | 562 | 0.43 |
| Recruitment success | 110 417 165 | 154 | 76–272 | 440 | 0.39 |
| M1 (ice at spawning) | 595 485 067 | 0 | 0–0 | 0 | — |
| M4 (calyptopis starvation) | 134 504 220 | 48 | 22–108 | 72 | 0.64 |
| M5a (under-developed) | 131 306 367 | 138 | 72–250 | 296 | 0.53 |
| M5b (off-shelf) | 282 050 117 | 318 | 178–546 | 834 | 0.42 |
| M6 (no winter ice) | 803 211 965 | 1 262 | 778–1 652 | 2 676 | 0.46 |
| Censored | 27 221 228 | 588 | 344–960 | 1 522 | 0.42 |
| Domain exit | 16 069 207 | 1 780 | 1 472–1 904 | 3 232 | 0.55 |

Particles evaluated at sea-ice advance are displaced a median of 224 km from
their release position, with an interquartile range of 112–414 km. They cover
562 km along track to get there, giving a straightness index of 0.43: less
than half the distance travelled becomes net displacement.

Successes travel least of the three, a median of 154 km, consistent with
recruitment resolving close to its origin at subarea resolution (Table S2).
M5b, the off-shelf failures, travel twice as far at 318 km — the same origins,
a different destination. M6 travels furthest of the classified outcomes at
1 262 km, eight times the success median, which is the export regime as a
number. Domain exit is larger still at 1 780 km, as expected for particles
that reach the boundary.

M1 displacement is identically zero, confirming that its fate position is its
release position by construction. Its straightness is undefined and reported
as such.

M4 is the shortest non-zero displacement at 48 km, and the straightest at
0.64. Both follow from its fate moment: the starvation threshold is crossed
during the calyptopis stages, within roughly the first month, so these
particles have had far less time to travel or to meander than any other
outcome.

## What this does and does not establish

The displacement is measured to the moment each particle's fate is
determined, which differs by outcome: the advance event for success, M5a and
M5b, the starvation threshold for M4, and the end of tracking for M6.
Distances are therefore not directly comparable between outcomes without also
considering `travel_time`, which the archived dataset carries.

The straightness index is a property of the flow sampled over these
trajectories, not a measure of how far larvae could disperse in principle. A
low value reflects the eddying, recirculating character of the Antarctic
Circumpolar Current over the release depth band (Figure S2), not any
behaviour of the particles, which are passive.
