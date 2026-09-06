# FS2 — Mean circulation

Supporting Information figure (Figure S2): the 32-year mean circulation over the
depth band in which particles are released, as speed of the mean velocity vector
with streamlines of the same field.

Cited from the Results as the physical basis for the eastward transport that
moves off-shelf failures and no-winter-ice particles away from their release
subareas.

## Scope

- **Period:** November 1993 – October 2025, 384 monthly files, 11 688 daily
  fields. The four seasonal means (DJF, MAM, JJA, SON) are computed in the same
  pass and stored, but only the annual mean is plotted — see "Seasonality" below.
- **Depth band:** 50–200 m, thickness-weighted. Nine GLORYS12 levels (centres
  47.4–186.1 m), with the weights of the top and bottom cells clipped to the
  band, summing to 150.0 m.
- **Domain and projection:** 115°W to 40°E, 78°S to 40°S; South Polar
  Stereographic, central_longitude = −37.5°, identical to F1 / F4 / F5.
- Single panel, 6.5 × 4 in.

## What is plotted, and what it does not mean

Four points that the caption also has to carry:

- **The speed of the mean vector**, |(u̅, v̅)|, not the mean of the speed. This
  measures persistent net transport, which is what displaces a particle over 200
  days. The mean of |(u, v)| would instead measure how energetic the flow is.
  Consequence: pale regions are regions of weak *net* transport, not necessarily
  of weak flow — an energetic eddy field produces little net displacement and
  renders pale here.
- **Velocity, not depth-integrated transport (U·H).** Transport is the right
  quantity for volume fluxes, but it scales with the available water column and
  so understates fast flow over the shelf, which is exactly where the retention
  argument lives. Particles are advected by velocity.
- **50–200 m is the release band, not where larvae live.** Particles are advected
  in three dimensions with vertical mixing over 200 days and sample well outside
  it. This is the band every trajectory starts in.
- **Blank where the seabed is above 50 m**, because no model level falls in the
  band there: 23.2% of the domain is land or too shallow, 1.7% is partial (the
  seabed lies within the band, so fewer levels are averaged), 75.1% is the full
  150 m. White means "not represented", not "no flow". The per-cell averaging
  thickness is written to the aggregation as `depth_averaged`.

## Design

- **Logarithmic colour scale**, `LogNorm(0.005, 0.5)` on `Blues`. The field spans
  more than a decade: domain median 0.030 m/s, 90th percentile 0.120, 99th 0.283,
  maximum 0.720. On a linear scale the ACC would be the only visible feature and
  the entire shelf would render as uniform white. `vmin` sits below the median so
  the shelf resolves; `vmax` sits just above the 99th percentile so almost
  nothing saturates. The colorbar is drawn with `extend="both"` so the clipping
  at each end is visible rather than hidden, and `plot.py` prints the fraction of
  cells outside both bounds.
- **Streamlines** of the same mean field, subsampled 6× before integration (at
  the native 1/12° the result is illegible at domain scale and slow to compute),
  `density=2.5`, dark gray at linewidth 0.4. Masked cells are set to zero
  velocity so the integrator stops at the mask rather than failing on it.
- **Boundary overlays**, linestyle convention consistent with F1 / F4 / F5:
  - CCAMLR subarea outlines: plain light gray.
  - 48.6 split at 60°S: dotted light gray — densified to follow the curved
    parallel and clipped to the 48.6 polygon.
  - Model computational domain: dashed light gray — external boundary.
- Coastlines and gray land fill via cartopy NaturalEarth (50 m resolution); axis
  spines pushed above them so the panel border is never covered.
- `grid_1d()` verifies that `nav_lon` / `nav_lat` are separable to 1e-4 degrees
  before reducing them to the 1-D axes `streamplot` requires. A curvilinear grid
  would otherwise be silently distorted rather than rejected.

## Seasonality

The seasonal means differ from the annual mean by a few percent in every
domain-wide statistic (medians 0.0288–0.0312 m/s, 90th percentiles 0.120–0.123),
so four panels would show four nearly identical maps. The annual mean is
therefore the one plotted. The seasonal fields remain in the aggregation and can
be plotted with `--period`; note the domain-wide statistics do not exclude
regional seasonality, only rule it out as a first-order effect at this scale.

## Reproducibility: this figure is the repo's one exception

`data/aggregated.nc` is **not committed** (~58 MB on disk, 94 MB in memory), so
unlike every other figure in this repo, `plot.py` alone is not enough. To
reproduce the figure you must first re-run `aggregate.py` on a machine with the
GLORYS12 velocity fields, then copy `data/aggregated.nc` back.

`aggregate.py` reads roughly 1.5 TB over 384 month-pairs and takes ~4.5 h on one
ECMWF node, so it is a batch job, not an interactive one; `job.sh` is the SLURM
script used. It checkpoints every 12 month-pairs to `data/checkpoint.npz`
(gitignored) and skips completed months on restart, so a timeout costs at most 12
month-pairs. `--no-resume` forces a clean run; `--limit N` processes the first N
month-pairs for testing.

## Files

- `aggregate.py` — reads the reformatted GLORYS12 velocity files, takes the
  thickness-weighted depth mean over 50–200 m of each daily field, accumulates
  sums and counts per period, writes `data/aggregated.nc` (gitignored). Runs on
  the HPC.
- `job.sh` — SLURM submission script for `aggregate.py`.
- `plot.py` — reads `data/aggregated.nc`, writes `currents.png`. Runs locally.

## Inputs (read by `aggregate.py`)

- `$KRICO_ROOT/Pre/GLORYS12/glorys12_u_YYYY_MM.nc`
- `$KRICO_ROOT/Pre/GLORYS12/glorys12_v_YYYY_MM.nc`

Both components must be present for a month to be processed; a `u` file without
its `v` is skipped with a warning. The reformatted GLORYS12 fields are produced
by the preprocessing step in
[krico-templates](https://github.com/ogourgue/krico-templates).

## Run

On the HPC, once:

```bash
python3 aggregate.py --limit 2      # check the level table and paths
sbatch job.sh                        # full run, ~4.5 h
```

Then copy `data/aggregated.nc` back and, locally:

```bash
python plot.py
python plot.py --period DJF          # optional, any of DJF/MAM/JJA/SON
```
