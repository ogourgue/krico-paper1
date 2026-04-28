# F4 — Source maps by fate

Multi-panel grid of release-position concentration maps, one panel per outcome. For each fate (success and the five killed_M* categories), shows where in the domain the particles ending in that fate were *released from*, normalized within the outcome.

## Scope

- Aggregates over all 32 spawning years (1994–2025) and all release dates within each season — climatological.
- 6 panels in a 3 × 2 grid; success top-left.
- Per-cell quantity: fraction of particles of a given outcome that were released in that cell. Sums to 1 within each panel.
- Grid: 1° × 1° lon/lat.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

## Design

- **Colormaps (echoing F3):**
  - `Greens` for success.
  - `Reds` for the five killed_M* panels (shared colormap, comparable within the killed group).
- **Single shared vmax** across both colormaps, computed as the 99th percentile of all nonzero density cells pooled together.
- **Two colorbars** at the bottom of the figure, one per color group, centered.
- **Panel labels:** combined letter + fate name in a framed box, upper-right.
- Coastlines + gray land fill via cartopy NaturalEarth.
- CCAMLR area outlines (light gray) and model domain boundary (dashed light gray) overlaid.
- Censored particles are folded into killed_M6 at plot time (consistent with F3 and the methodology decision); they share the same end-of-tracking, no-advance condition. The killed_M6 vs censored breakdown lives in the SI.
- `exited_domain` is not shown — it is a modeling-domain limitation (~1% of particles dataset-wide). A spatial map is in the SI.

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, computes per-outcome 2D density on a 1° × 1° grid using `release_lon` / `release_lat`, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, folds censored into killed_M6, writes `source_maps.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.