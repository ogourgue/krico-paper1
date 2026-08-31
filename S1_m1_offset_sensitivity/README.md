# S1 — Sensitivity of M1 to the descent-ascent offset

Supporting Information analysis for Paper 1.

## What this answers

M1 (ice at spawning) is a constraint acting on the spawning adult. The
23-26 day descent-ascent cycle from spawning to calyptopis I is not
simulated, so particles enter the model as calyptopis I and their release
date postdates spawning. The recruitment pipeline therefore evaluates M1 at
the release position on the spawning date, sampling the GLORYS12 sea-ice
field directly rather than reading the trajectory at day 0.

That requires a value for the offset. Thorpe et al. (2019) give 23-26 days
and the pipeline applies the midpoint, 24 days, as a constant
(`krico_recruitment.sea_ice.SPAWNING_OFFSET_DAYS`). This analysis measures
how much the M1 classification depends on that choice.

## Result

Across the full 32-year record (3 848 cohorts), varying the offset over
23-26 days changes the domain-wide M1 fraction by less than one percentage
point, against a much larger difference between evaluating M1 at spawning
and at release. Spread is largest for December releases, whose spawning
dates fall in mid-to-late November on the steepest part of the seasonal
retreat, and no single spawning year is materially more sensitive than the
record as a whole.

The constant is therefore adequate, and a temperature-dependent,
per-particle descent-ascent duration is not needed.

Run `compare.py` for the current numbers.

## Reference case and verification

The sweep includes offset 0, which samples the sea-ice field on the release
date itself. It is not a candidate value: it is a regression test. With
`--verify`, the offset-0 classification is compared against the sea-ice
concentration sampled along each trajectory at day 0.

Whether the two agree exactly depends on the date, for a reason external to
this analysis. The GLORYS12 ice fields are stamped at 12:00 up to
2018-12-18 and at 11:52:30 from 2018-12-19 onward. The velocity,
temperature and mixing fields are stamped at 12:00 throughout, so only the
ice product is affected.

Before the step, particle times coincide with the field timestamps: Parcels'
temporal interpolation returns the stored value, and the nearest-neighbour
lookup reproduces it **bit for bit**. After the step, every sample sits 7.5
minutes off the stamp, so Parcels returns a 450/86400 blend of two adjacent
daily fields while the lookup returns the day containing the sample. Adjacent
daily ice fields differ by order 10⁻², so the blend differs by order 10⁻⁴ for
most particles, with a consistent sign.

Measured over the full 32-year record:

| Quantity | Value |
|---|---|
| Particles classified differently by M1 | 1.8 × 10⁻⁵ |
| Mean absolute concentration difference | 1.2 × 10⁻⁵ |
| Rate, spawning years 1994-2018 | ~1 × 10⁻⁶ |
| Rate, spawning years 2019-2025 | 5 × 10⁻⁵ to 1.2 × 10⁻⁴ |

The rate follows the timestamp step, not anything physical: it is flat on
both sides and jumps between consecutive spawning years. Cohorts released
before the step but tracked past it show exact agreement at day 0 and order
10⁻⁴ disagreement at day 100, which is what identified the mechanism.

Neither behaviour is an error. Taking the field for the day containing the
spawning date is the intended semantics for M1, and Parcels interpolates
correctly for the filters that read sea ice along the trajectory (M4, M5,
M6). One particle in roughly 57 000 is classified differently, in a category
holding 28% of the dataset, so no reported quantity is affected.

`sweep.py` flags the verification above `VERIFY_TOLERANCE` (10⁻³), set well
above the observed rate: it exists to catch a genuine indexing failure, which
would produce percent-level disagreement.

## Files

| File | Description |
|---|---|
| `sweep.py` | Computes M1 at each offset for every cohort |
| `sweep.sh` | SLURM driver for the full-record sweep |
| `compare.py` | Summarizes the sweep into the SI numbers |
| `data/m1_offset_sweep.csv` | Sweep output, one row per cohort per offset |

M1 is computed by importing `evaluate_M1` and `spawning_sic` from
`krico_recruitment`, so this analysis cannot drift from the classification
it describes.

## Reproducing

Unlike the figure folders, this analysis reads the raw trajectories and the
GLORYS12 sea-ice fields, and imports the recruitment package. Three
environment variables are required:

- `KRICO_POST` — path to a [krico-post-production](https://github.com/ogourgue/krico-post-production) **clone** (not merely a directory of recruitment outputs), used to import `krico_recruitment`.
- `KRICO_RUNS` — path to the raw trajectory simulations.
- `KRICO_GLORYS12` — path to the GLORYS12 preprocessing output directory containing `glorys12_ice_YYYY_MM.nc`.

```bash
cd S1_m1_offset_sensitivity
sbatch sweep.sh                              # ~3 848 cohorts, 5 offsets
python compare.py --sweep data/m1_offset_sweep.csv
```

`sweep.sh` sets `PYTHONPATH` from `KRICO_POST`; its `#SBATCH` directives are
specific to ECMWF's Atos HPCF. To run without SLURM, or on a subset:

```bash
export PYTHONPATH="$KRICO_POST/recruitment:$PYTHONPATH"
python sweep.py --runs "$KRICO_RUNS" --glorys "$KRICO_GLORYS12" \
                --offsets 0 23 24 25 26 --verify --stride 30 \
                --out data/subset.csv
```

The sweep reads only day 0 of each trajectory file, so it does not require
the full trajectory archive to be staged in memory, but it does require the
release positions and the day-0 sea-ice sample.

Spawning year 1994 requires the October 1993 sea-ice field. The Mercator
OPeNDAP server moved behind an SSO portal whose sessions are not usable from
the ECMWF login nodes, so that month was retrieved by hand through the
browser data-access form and decoded from the DAP2 binary response. Note
that the form's `start:stride:stop` indexing is 1-based and inclusive,
whereas the download scripts slice with 0-based exclusive `isel`; the
resulting file was verified cell-by-cell against an existing month before
use.
