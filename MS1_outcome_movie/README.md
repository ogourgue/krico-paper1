# MS1 — Outcome movie

Supporting Information movie: one release-day cohort followed from release to the end of tracking, with particles turning green or red on the day their outcome is decided during tracking. It shows in motion what Figures 3 to 5 show in aggregate: the drift during summer, the sea-ice advance from April, and where larvae end up ready or not ready to overwinter.

## Scope

- **Cohort:** a single release day, chosen from F2's aggregation (see below). By default 29 January 2006, spawned on 5 January 2006, all ~546,000 particles.
- **Particles:** only particles that drift are shown.
  - translucent gray while drifting, at their daily position;
  - green at their fate position from the day recruitment success is decided, and they stay;
  - red at their fate position from the day M4, M5a or M5b is decided, fading out over `FADE_DAYS` (5) days so the map stays readable;
  - M1 particles are not shown: they are never spawned;
  - M6 particles (with censored folded in, as in the paper) stay gray to the end;
  - domain exits drift in gray until deleted, then disappear. They are a modelling limitation, not an outcome.
- **Sea ice:** daily GLORYS12 concentration above the 15% advance threshold, block-averaged to 1/4°, as a single colour, matching the threshold (`SIC_BINARY = True`; `False` shades by concentration instead). These are the fields the simulation and the classification use, not the observed climatologies of F1.
- **Shelf-slope zone:** bathymetry shallower than 2000 m, the bathymetric condition for recruitment success (its complement is M5b), as a translucent light gray drawn over the sea ice. It comes from F1's bathymetry, which is masked to the CCAMLR subareas, so the zone is drawn inside them only. The strip south of the GLORYS12 grid (about 77°S), where the bathymetry ends, is filled with the same colour across the domain's longitudes, since it is shelf in front of the Weddell ice shelves.
- **Base map:** F1's land, coastline, CCAMLR outlines and domain boundary, imported from `../F1_domain_map/plot.py` rather than copied. F1's subarea labels, 48.6 split, bathymetry zones, sea-ice climatologies and legends are left out; the movie's legend names the outlines as CCAMLR subareas.
- **Format:** 1080 × 1080 px frames, one per tracking day, from release (frame 0) to the end of tracking (frame 200).

Outcome timing comes from the archived fate day (`travel_time`) of the recruitment classification, so the movie cannot disagree with it:

- M4 is decided when the starvation threshold is crossed, during the calyptopis stages.
- Success, M5a and M5b are decided at each particle's own sea-ice advance, from 1 April onward, so they appear as the ice edge moves north.
- M6 is only decided at the end of tracking, which is where the movie ends: those particles are still gray on the last frame, and that gray cloud is the larvae winter sea ice never reached.
- M1 is decided at spawning, before release, so those particles never enter the movie.

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
- `make_movie.sh` — assembles the frames into `outcome_movie.mp4` with ffmpeg, at 10 frames per second (about 20 s), with no frame held at the start or the end.
- `job_series.sh` — SLURM array driver for the per-year series (below): runs the three steps above on the HPC for every spawning year.

`data/aggregated.nc`, `frames/` and `series/` are gitignored: at full resolution the aggregation is a few hundred MB, and the frames are regenerated from it.

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
./make_movie.sh                            # FPS=10 by default
```

Frame ranges are independent (`--frames 0:100`, `--frames 100:`), so a long render can be split across processes.

## Per-year series (not part of the paper)

A variant made at the co-authors' request for a stakeholder meeting: one movie per spawning year, 1994 to 2025, all for the cohort released on the climatological peak day (29 January, the same step 1 as the cohort choice above, read from F2 by `aggregate.py --year`). The design is identical to `outcome_movie.mp4`; only the cohort changes, so the years can be compared side by side. Only the scripts are committed; the movies are kept locally and passed on directly. The whole chain runs on the HPC, one array task per year, with the frames split across the task's cores:

```bash
cd MS1_outcome_movie
bash job_series.sh --preflight        # on a login node: environment, packages, Natural Earth files, ffmpeg
mkdir -p logs
sbatch job_series.sh                  # 32 tasks; sbatch --array=16 job_series.sh for 2010 only
```

Each task writes `series/YYYY/aggregated.nc`, `series/YYYY/frames/` and `series/outcome_movie_YYYY.mp4`; logs go to `logs/ms1_series_<jobid>_<task>.out`. A re-run reuses an existing `aggregated.nc` and re-renders the frames (`REDO=1 sbatch job_series.sh` to re-aggregate as well). The frames need cartopy and geopandas from the `python3` module, the Natural Earth 50m land and coastline files, which cartopy downloads once (the tasks serialise on a lock), and an ffmpeg binary, found on the PATH, in an `ffmpeg` module, or bundled with `imageio-ffmpeg`.
