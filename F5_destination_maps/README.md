# F5 — Destination maps by fate

Multi-panel grid of fate-position concentration maps, one panel per outcome. For each fate (success and the five killed_M* categories), shows where in the domain the particles were *at the moment fate was determined*, normalized within the outcome.

## Fate position by outcome

The "destination" interpretation depends on what defines the moment of fate (`fate_day` in the recruitment pipeline):

- **Success / Killed M5 (no FIV) / Killed M5 (off-shelf):** position when the winter sea-ice advance event was detected. For successes, this is the recruitment moment — particle had reached FIV and was on shelf when ice arrived. For both M5 categories, particle had ice arrive but failed one of the recruitment criteria.
- **Killed M4:** position when the 10-consecutive-days SIC > 40% threshold was crossed during the calyptope window.
- **Killed M1:** position at release (fate is determined at day 0). Panel (b) is therefore identical to F4 panel (b) by construction.
- **Killed M6 (no advance):** position at the last valid trajectory day (typically day 200 for particles surviving the full tracking window). After folding censored particles into killed_M6, this panel reflects the spatial distribution of all particles that survived the calyptope window but did not encounter winter sea-ice advance within the tracking window.

## Scope

- Aggregates over all 32 spawning years (1994–2025) and all release dates within each season — climatological.
- 6 panels in a 3 × 2 grid; success top-left.
- Per-cell quantity: fraction of particles of a given outcome whose fate position fell in that cell. Sums to 1 within each panel.
- Grid: 1° × 1° lon/lat.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

## Design

Identical to F4 (see `../F4_source_maps/README.md`):

- **Colormaps:** `Greens` for success, `Reds` for the five killed_M* panels (shared scale across all six panels).
- **Single shared vmax** computed as the 99th percentile of all nonzero density cells pooled together.
- **Two colorbars** at the bottom of the figure, one per color group, centered under their column.
- **Panel labels:** combined letter + fate name in a framed box, upper-right.
- Coastlines + gray land fill via cartopy NaturalEarth.
- CCAMLR area outlines (light gray) and model domain boundary (dashed light gray) overlaid.
- Censored particles are folded into killed_M6 at plot time (consistent with F3 and F4); they share the same end-of-tracking, no-advance condition. The killed_M6 vs censored breakdown lives in the SI.
- `exited_domain` is not shown — it is a modeling-domain limitation (~1% of particles dataset-wide). A spatial map is in the SI.

F5's color scale is independent of F4's. Source positions are constrained to the narrow bathymetric release zone (1000–2000 m), while destinations spread across the entire domain — the inherent dynamic range differs, so each figure is normalized to its own data.

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, computes per-outcome 2D density on a 1° × 1° grid using `final_lon` / `final_lat`, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, folds censored into killed_M6, writes `destination_maps.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.