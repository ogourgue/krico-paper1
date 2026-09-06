# FS1 — Per-year phenology

Supporting Information figure (Figure S1), companion to F2: the 32 individual
per-year recruitment success curves, and the distribution of the release date at
which each year peaks.

## Scope

- Same underlying quantity as F2 — recruitment success rate as a function of
  release date — resolved by spawning year instead of averaged over years.
- Panel (a): all 32 per-year curves (1994–2025), the 32-year climatological mean
  overlaid, each year's own peak marked, and two vertical lines for the
  mean-curve peak and the median per-year peak.
- Panel (b): histogram of the 32 per-year peak release dates, with median and
  inter-quartile range.
- X-axis shared between panels, so the peak scatter in (b) aligns with the
  curves in (a). Release date, Nov 15 → Mar 14.
- Feb 29 cohorts excluded, as in F2 and F3.

## Why this figure exists

F2's shaded envelope describes the spread of success *rates* at a given release
date. It says nothing about the spread of *peak dates*, and the two are easily
conflated: because `max()` is nonlinear, the maximum of the mean curve is not
the mean of the per-year maxima. This figure separates them.

The climatological curve peaks on January 29 at 7.12%; the per-year peaks have a
median of January 24 and a mean value of 7.61%. The mean curve therefore peaks
five days later and half a point lower than a typical individual year, because
averaging curves whose maxima fall on different dates broadens and flattens the
composite.

The per-year peak-date distribution is **unimodal**: one KDE mode at bandwidths
from 0.8× to 1.5× Scott's factor, stable under rebinning at 5, 7, 10 and 14 days
and under bin-origin shifts, with mean minus median of 0.5 days and skew −0.08.
The apparent dip in the plotted 7-day bins is a single count at the edge of the
distribution, and the central run is flat within Poisson noise.

Keeping this in Supporting Information rather than as an inset on F2 holds the
boundary between this paper and the interannual-variability paper: the main text
stays climatological.

## Design

- Two stacked panels of equal height, 6.5 × 8 in total, sized to fit a US Letter
  page with 1 in margins (6.5 × 9 in live area) with room for the caption.
- **Colour roles**, consistent across both panels: black for the raw per-year
  curves, `C0` for climatological-mean quantities, `C1` for per-year quantities.
  An earlier one-colour-per-panel scheme was rejected because it could not
  express that the median peak in (a) and the histogram in (b) are the same
  quantity.
  - (a) per-year curves black, linewidth 1.0, alpha 0.1; climatological mean
    `C0` at the default linewidth; per-year peaks as open markers (white face,
    `C1` edge, opaque); mean-curve peak `C0` dashed; median peak `C1` dashed.
  - (b) histogram `C1`, 7-day bins, `rwidth=0.8`, alpha 0.5; IQR span `C1` at
    alpha 0.2; median peak `C1` dashed.
- Line width 1.0 for the per-year curves is the only departure from the rcParams
  defaults: at the default 1.5 the 32 curves fill the panel, and alpha alone does
  not separate them from the mean.
- Legends opaque (`framealpha=1.0`), as elsewhere in the manuscript. Panel (a)
  legend order is set explicitly by handle rather than by draw order.
- Y-grid at the matplotlib default grey (`#b0b0b0`), alpha 0.5, matching F2.
- X-axis ticks identical to F2 (Dec 1, Jan 1, Feb 1, Mar 1), so the two figures
  read the same way.
- **Panel (a) z-order**, bottom to top: per-year curves (2), median peak line
  (3), mean-curve peak line (4), climatological mean (5), per-year peak markers
  (6). Numbering starts at 2 because the default `axes.axisbelow="line"` draws
  the grid at z-order 1.5; starting at 1 would put gridlines over the curves.

## Files

- `plot.py` — reads the **F2** aggregation and writes `phenology_per_year.png`.

There is no `aggregate.py` and no `data/` here. This figure needs exactly the
quantities F2 already computes, so `plot.py` reads
`../F2_phenology_curve/data/aggregated.nc` directly rather than re-walking ~3800
cohort NetCDFs. This is a deliberate departure from the self-contained-folder
convention used elsewhere in the repo: duplicating the aggregation would let the
two figures disagree about the same numbers. Pass `--aggregated` to point at a
different file.

## Run

```bash
python plot.py
```

Requires `../F2_phenology_curve/data/aggregated.nc`, which is committed. No
environment variables and no external data are needed.

`plot.py` also prints both peak estimators, the per-year peaks sorted by date,
and a text histogram of the peak-date distribution — the numbers quoted above
and in the figure caption.
