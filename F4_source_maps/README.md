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
- **Panel labels:** combined letter + fate name in a framed box, upper-right. Manuscript labeling convention: plain-English name first, M-code in parentheses (e.g., "Calyptopis starvation (M4)"). M5a and M5b are abbreviated to "Under-developed (M5a)" and "Off-shelf (M5b)" for compactness.
- Coastlines + gray land fill via cartopy NaturalEarth.
- **Boundary overlays** (linestyle convention consistent with F1):
  - CCAMLR subarea outlines: plain light gray.
  - 48.6 split at 60°S: dotted light gray — internal subdivision used in the per-CCAMLR-subarea breakdown (`ccamlr_summary.py`); densified to follow the curved parallel and clipped to the 48.6 polygon.
  - Model computational domain: dashed light gray — external boundary.
- Censored particles are folded into M6 at plot time (consistent with F3 and the methodology decision); they share the same end-of-tracking, no-advance condition. The `killed_M6` vs censored breakdown lives in the SI.
- `exited_domain` is not shown — it is a modeling-domain limitation (~1% of particles dataset-wide). A spatial map is in the SI. (F3 retains it as a noise-floor band in the stacked-area for accounting completeness.)

## Per-CCAMLR-subarea breakdown

The companion utility `ccamlr_summary.py` reads the gridded `data/aggregated.nc` and computes a per-fate × per-CCAMLR-subarea breakdown. Subarea 48.6 is split at 60°S into Northern Bouvet (48.6N) and Southern Bouvet (48.6S). The script writes:

- A formatted table to stdout (per-fate columns × per-subarea rows, percentages summing to 100% within each fate).
- `ccamlr_summary.csv` — wide-format CSV with the same data, `#`-prefixed comment header. **Committed to the repo as the source for SI Table S1.**

The "outside" column captures particles whose release position fell outside any CCAMLR project subarea. For F4 this is essentially zero by construction (releases are constrained to bathymetric spawning zones inside CCAMLR); the column is kept for symmetry with F5 where it carries real signal (M6 destinations exported by the ACC).

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, computes per-outcome 2D density on a 1° × 1° grid using `release_lon` / `release_lat`, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, folds censored into `killed_M6`, writes `source_maps.png`.
- `ccamlr_summary.py` — reads `data/aggregated.nc`, computes per-fate × per-subarea breakdown, writes `ccamlr_summary.csv` (committed) and prints the same data to stdout.

## Run

```bash
python aggregate.py
python plot.py
python ccamlr_summary.py
```

Running `aggregate.py` requires `$KRICO_POST` to be set; see the repo README.