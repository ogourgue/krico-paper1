# MS1 — Outcome movie

Supporting Information movie: one release-day cohort followed from release to the end of tracking, with every particle turning green or red on the day its outcome is decided. It shows in motion what Figures 3 to 5 show in aggregate: the drift during summer, the sea-ice advance from April, and where larvae end up ready or not ready to overwinter.

## Scope

- **Cohort:** a single release day, chosen from F2's aggregation (see below). By default 29 January 2006, spawned on 5 January 2006, all ~546,000 particles.
- **Particles:**
  - gray while drifting, at their daily position;
  - green at their fate position from the day recruitment success is decided, and they stay;
  - red at their fate position from the day any mortality outcome (M1, M4, M5a, M5b, M6) is decided, fading out over `FADE_DAYS` (5) days so the map stays readable;
  - domain exits drift in gray until deleted, then disappear. They are a modelling limitation, not an outcome.
- **Sea ice:** daily GLORYS12 concentration from the 15% advance threshold upward, block-averaged to 1/4°. These are the fields the simulation and the classification use, not the observed climatologies of F1.
- **Isobath:** 2000 m, the bathymetric limit of recruitment habitat (the M5b criterion), from F1's masked bathymetry.
- **Base map:** F1's land, coastline, CCAMLR outlines, 48.6 split and domain boundary, imported from `../F1_domain_map/plot.py` rather than copied. F1's bathymetry zones, sea-ice climatologies and legends are left out, and region names replace the subarea codes.
- **Format:** 1080 × 1080 px frames, one per tracking day (201) plus `FADE_DAYS` hold frames at the end.

Outcome timing comes from the archived fate day (`travel_time`) of the recruitment classification, so the movie cannot disagree with it:

- M1 is decided at spawning, before release: those particles are red on the first frame and never drift.
- M4 is decided when the starvation threshold is crossed, during the calyptopis stages.
- Success, M5a and M5b are decided at each particle's own sea-ice advance, from 1 April onward, so they appear as the ice edge moves north.
- M6 (with censored folded in, as in the paper) is decided at the end of tracking: those particles all turn red on the last day, and the hold frames let them fade while success stays.

## Cohort choice

`aggregate.py` picks the cohort from `../F2_phenology_curve/data/aggregated.nc` instead of a hand-picked date:

1. the climatological peak release date, the maximum of the 32-year mean curve (29 January, as in Figure 2);
2. among spawning years up to 2016, the year whose own peak release date is closest to it;
3. ties broken by the peak value closest to the median of per-year peaks.

The cohort is that year's own peak release date. Years after 2016 are excluded so the movie cannot be read as a statement about the post-2016 low-ice regime. With the current data this gives spawning year 2006, whose own peak falls exactly on 29 January at 7.98% (median per-year peak 7.71%, season mean 5.67%). `--date` overrides the choice.

## Files

- `aggregate.py` — reads the cohort's raw trajectory, its recruitment file and the daily GLORYS12 sea ice, and writes `data/aggregated.nc`: daily positions up to each particle's fate day (int16-packed, fill afterwards), outcome, fate day and fate position per particle, and the coarsened daily sea ice. Runs on the HPC.
- `job.sh` — SLURM driver for `aggregate.py`.
- `plot.py` — reads `data/aggregated.nc` and writes `frames/frame_NNN.png`. Runs locally.
- `make_movie.sh` — assembles the frames into `outcome_movie.mp4` with ffmpeg, holding the first and last frames.

`data/aggregated.nc` and `frames/` are gitignored: at full resolution the aggregation is a few hundred MB, and the frames are regenerated from it.

## Inputs (read by `aggregate.py`)

- `$KRICO_RUNS/KRICO_*/YYYY_MM_DD.nc` — raw trajectory of the cohort (not publicly archived; see krico-post-production)
- `$KRICO_POST/recruitment/data/YYYY_MM_DD.nc` — recruitment outcomes, v2.0.0
- `$KRICO_GLORYS12/glorys12_ice_YYYY_MM.nc` — daily sea-ice concentration, for the months covered by tracking
- `../F2_phenology_curve/data/aggregated.nc` — cohort choice (committed)

The trajectory and recruitment files are matched row by row, as written by `process_cohort.py`. `aggregate.py` checks this by comparing each particle's trajectory position on its fate day with the archived fate position, and stops if they disagree.

## Run

On the HPC:

```bash
cd MS1_outcome_movie
sbatch job.sh                          # or: sbatch job.sh --date 2006-01-29
```

Then copy `data/aggregated.nc` to a local clone, and:

```bash
python plot.py --frames 0:1 --stride 10    # quick look at one frame
python plot.py                             # all frames, a few seconds each
./make_movie.sh                            # FPS=12 by default
```

Frame ranges are independent (`--frames 0:100`, `--frames 100:`), so a long render can be split across processes.
