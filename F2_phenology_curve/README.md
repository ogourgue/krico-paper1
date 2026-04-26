# F2 — Phenology curve

Headline figure for KRICO Paper 1: 30-year climatological mean recruitment success rate as a function of release date, with 5th–95th percentile ribbon across years.

## Scope

- Aggregates over all source regions (no spatial decomposition).
- Aggregates over all 32 spawning years (1994–2025).
- Y-axis: success rate (%).
- X-axis: release season-day, 0 (Nov 15) → 120 (Mar 15).
- Annotations: climatological peak day, optimum window width (days within 80% of peak).

Feb 29 cohorts are excluded from the climatology to keep the season-day grid uniformly 121 days across leap and non-leap years.

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs from `$KRICO_ROOT/Post/Production/recruitment/data/`, aggregates to `(year, season_day)` grid, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, writes `phenology_curve.png`.

## Run

```bash
python aggregate.py    # ~3800 cohort NetCDFs → small aggregated.nc
python plot.py         # aggregated.nc → phenology_curve.png
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.
