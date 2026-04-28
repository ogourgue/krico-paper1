# F5 — Destination maps by fate

Multi-panel grid of fate-position concentration maps, one panel per outcome. For each fate (success, killed_M1, killed_M4, killed_M5_no_FIV, killed_M5_not_on_shelf, exited_domain), shows where in the domain the particles were *at the moment fate was determined*, normalized within the outcome.

## Fate position by outcome

The "destination" interpretation depends on what defines the moment of fate (`fate_day` in the recruitment pipeline):

- **Success / Killed M5 (no FIV) / Killed M5 (off-shelf):** position when the winter sea-ice advance event was detected. For successes, this is the recruitment moment — particle had reached FIV and was on shelf when ice arrived. For both M5 categories, particle had ice arrive but failed one of the recruitment criteria.
- **Killed M4:** position when the 10-consecutive-days SIC > 40% threshold was crossed during the calyptope window.
- **Killed M1:** position at release (fate is determined at day 0). Panel (b) is therefore identical to F4 panel (b) by construction.
- **Killed M6 (no advance):** position at end of tracking (day 200). Negligible in practice (~0.01% of dataset, dropped from the figure).
- **Exited domain:** position at deletion (last valid trajectory day).

## Scope

- Aggregates over all 32 spawning years (1994–2025) and all release dates within each season — climatological.
- 6 panels in a 3 × 2 grid; success top-left.
- Per-cell quantity: fraction of particles of a given outcome whose fate position fell in that cell. Sums to 1 within each panel.
- Grid: 1° × 1° lon/lat.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

## Design

Identical to F4 (see `../F4_source_maps/README.md`):

- **Colormaps:** `Greens` for success, `Reds` for the four killed_M* panels (shared scale), `Greys` for exited_domain.
- **Single shared vmax** across all colormaps, computed as the 99th percentile of all nonzero density cells pooled together.
- **Three colorbars** at the bottom of the figure, one per color group.
- **Panel labels:** combined letter + fate name in a framed box, upper-right.
- Coastlines + gray land fill via cartopy NaturalEarth.
- CCAMLR area outlines (light gray) and model domain boundary (dashed light gray) overlaid.
- `killed_M6_no_advance` and `censored` are dropped (each well below 0.05% over the full dataset).

F5's color scale is independent of F4's. Source positions are constrained to the narrow bathymetric release zone (1000–2000 m), while destinations spread across the entire domain — the inherent dynamic range differs, so each figure is normalized to its own data.

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, computes per-outcome 2D density on a 1° × 1° grid using `final_lon` / `final_lat`, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, writes `destination_maps.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.
