# F4 — Source maps by fate

Multi-panel grid of release-position concentration maps, one panel per outcome. For each fate (success, killed_M1, killed_M4, killed_M5_no_FIV, killed_M5_not_on_shelf, exited_domain), shows where in the domain the particles ending in that fate were *released from*, normalized within the outcome.

## Scope

- Aggregates over all 32 spawning years (1994–2025) and all release dates within each season — climatological.
- 6 panels in a 3 × 2 grid; success top-left.
- Per-cell quantity: fraction of particles of a given outcome that were released in that cell. Sums to 1 within each panel.
- Grid: 1° × 1° lon/lat.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

## Design

- **Colormaps (echoing F3):**
  - `Greens` for success.
  - `Reds` for the four killed_M* panels (shared colormap, comparable within the killed group).
  - `Greys` for exited_domain.
- **Three colorbars** at the bottom of the figure, one per color group.
- **Panel labels:** letter (a–f) at top-left, fate name at top-center, both in framed boxes matching the legend style from earlier figures.
- Coastlines + gray land fill via cartopy NaturalEarth.
- `killed_M6_no_advance` and `censored` are dropped (each well below 0.05% over the full dataset).

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, computes per-outcome 2D density on a 1° × 1° grid using `release_lon` / `release_lat`, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, writes `source_maps.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.
