# F3 — Outcome composition

Stacked-area showing what happens to each particle as a function of release date: 32-year climatological fractional breakdown into the recruitment outcome states.

## Scope

- Aggregates over all source regions (no spatial decomposition).
- Aggregates over all 32 spawning years (1994–2025).
- X-axis: release date, Nov 15 → Mar 14.
- Y-axis: cumulative fraction (%).
- Layers, bottom to top (with manuscript labels and internal flag names):
  - Domain exit — `exited_domain` — particle left the model domain (a modeling-domain limitation, not a biological outcome; placed at the bottom as a "noise floor" so the y-axis sums to 100%)
  - Ice at spawning (M1) — `killed_M1` — at release (sea ice too dense to spawn)
  - Calyptope starvation (M4) — `killed_M4` — during calyptope window (starvation under sea ice)
  - Under-developed (M5a) — `killed_M5_no_FIV` — at sea-ice advance event (development too slow)
  - Off-shelf (M5b) — `killed_M5_not_on_shelf` — at sea-ice advance event (off-shelf at advance)
  - No winter ice (M6) — `killed_M6_no_advance` — alive at end of tracking, no winter sea-ice advance detected
  - Recruitment success — `success` — at sea-ice advance event (top of stack)

Censored particles (alive at end of tracking, no advance detected, but SIC rising at cutoff) are folded into `killed_M6_no_advance` for this figure. They share the same end-of-tracking, no-advance condition; the "censored" label only marks the subset whose classification is provisional. The `killed_M6` vs censored breakdown lives in the SI.

Feb 29 cohorts are excluded (consistent with F2).

## Design notes

- Colors: `success` is matplotlib `C2` (green); the five `killed_M*` layers use a perceptually uniform gradient sampled from `YlOrRd`, lightest = earliest filter (M1), darkest = latest filter (M6_no_advance); `exited_domain` is a light gray (`"0.7"`).
- Legend labels follow the manuscript convention: plain-English name first, M-code in parentheses (e.g., "Calyptope starvation (M4)"). M5a and M5b are abbreviated to "Under-developed (M5a)" and "Off-shelf (M5b)" for legend compactness; full "at sea-ice advance" definitions are in the figure caption, with the M-code connecting them.
- No y-axis grid: grid lines compete with filled bands and would only show through the lightest layers, providing inconsistent reference. Y-values are read from band thicknesses, not axis position.
- Legend inside, lower-left corner. The legend covers the early-season M1↔M4 boundary in a region where M1 dominates uniformly and the boundary is roughly horizontal — least informational cost of available placements. The upper-right would obscure the late-season M6 dominance and the success peak, both of which carry more of the paper's story.

## Files

- `aggregate.py` — reads recruitment cohort NetCDFs, aggregates per-(year, season_day) outcome counts, writes `data/aggregated.nc`.
- `plot.py` — reads `data/aggregated.nc`, writes `outcome_composition.png`.

## Run

```bash
python aggregate.py
python plot.py
```

Both scripts assume `$KRICO_ROOT` is set and that the recruitment pipeline output exists at the expected location.