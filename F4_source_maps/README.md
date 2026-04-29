# F4 — Source maps by fate

Multi-panel grid of release-position concentration maps, one panel per outcome. For each fate (recruitment success and the five mortality categories M1, M4, M5a, M5b, M6), shows where in the domain the particles ending in that fate were *released from*, normalized within the outcome.

## Scope

- Aggregates over all 32 spawning years (1994–2025) and all release dates within each season — climatological.
- 6 panels in a 3 × 2 grid; success top-left.
- Per-cell quantity: fraction of particles of a given outcome that were released in that cell. Sums to 1 within each panel.
- Grid: 1° × 1° lon/lat.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

## Design

- **Colormaps (echoing F3):**
  - `Greens` for recruitment success.
  - `Reds` for the five mortality panels (shared colormap, comparable within the mortality group).
- **Single shared vmax** across both colormaps, computed as the 99th percentile of all nonzero density cells pooled together.
- **Two colorbars** at the bottom of the figure, one per color group, centered. Labels: "Recruitment success density" (under col 0) and "Mortality density" (under col 1).
- **Panel labels:** combined letter + fate name in a framed box, upper-right. Manuscript labeling convention from `../METHODOLOGY.md` §2.7: plain-English name first, M-code in parentheses (e.g., "Calyptope starvation (M4)"). M5a and M5b are abbreviated to "Under-developed (M5a)" and "Off-shelf (M5b)" for compactness.
- Coastlines + gray land fill via cartopy NaturalEarth.
- CCAMLR subarea outlines (light gray) and model domain boundary (dashed light gray) overlaid.
- Censored particles are folded into M6 at plot time (consistent with F3 and the methodology decision); they share the same end-of-tracking, no-advance condition. The `killed_M6` vs censored breakdown lives in the SI.
- `exited_domain` is not shown — it is a modeling-domain limitation (~1% of particles dataset-wide). A spatial map is in the SI. (F3 retains it as a noise-floor band in the stacked-area for accounting completeness.)

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, computes per-outcome 2D density on a 1° × 1° grid using `release_lon` / `release_lat`, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, folds censored into `killed_M6`, writes `source_maps.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.