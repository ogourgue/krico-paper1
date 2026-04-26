# F3 — Outcome composition

Stacked-area showing what happens to each particle as a function of release date: 32-year climatological fractional breakdown into the recruitment outcome states.

## Scope

- Aggregates over all source regions (no spatial decomposition).
- Aggregates over all 32 spawning years (1994–2025).
- X-axis: release date, Nov 15 → Mar 14.
- Y-axis: cumulative fraction (%).
- Layers, bottom to top:
  - `exited_domain` — particle left the model domain (a modeling-domain limitation, not a biological outcome; placed at the bottom as a "noise floor")
  - `killed_M1` — at release (sea ice too dense to spawn)
  - `killed_M4` — during calyptope window (starvation under sea ice)
  - `killed_M5_no_FIV` — at sea-ice advance event (development too slow)
  - `killed_M5_not_on_shelf` — at sea-ice advance event (off-shelf at advance)
  - `success` — at sea-ice advance event (top of stack)

`killed_M6_no_advance` and `censored` are negligible (each well below 0.05% over the full 30 years) and dropped from the stack at plot time. The `aggregate.py` script always computes all 8 counts; `plot.py` reports which categories were dropped to stdout.

Feb 29 cohorts are excluded (consistent with F2).

## Design notes

- Colors: `success` is matplotlib `C2` (green); the four `killed_M*` layers use a perceptually uniform gradient sampled from `YlOrRd`, lightest = earliest filter (M1), darkest = latest filter (M5_not_on_shelf); `exited_domain` is a light gray (`"0.7"`).
- No y-axis grid: grid lines compete with filled bands and would only show through the lightest layers, providing inconsistent reference. Y-values are read from band thicknesses, not axis position.
- Legend inside, lower-left corner, on top of the `exited_domain` band where contrast is good and there is no biological signal to obscure.

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, aggregates per-(year, season_day) outcome counts, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, writes `outcome_composition.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.