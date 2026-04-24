# KRICO Paper 1

Figures and analysis for KRICO Paper 1: *Spawning phenology and upstream connectivity jointly determine larval recruitment success in Antarctic krill*.

Each figure lives in its own self-contained subfolder with aggregation and plotting scripts.

## Dependencies

Reads recruitment outcome data produced by [krico-post-production](https://github.com/ogourgue/krico-post-production), specifically the `recruitment/` pipeline. Expected location: `$KRICO_ROOT/Post/Production/recruitment/data/`.

## Layout

```
F1_domain_map/
F2_phenology_curve/
F3_outcome_composition/
F4_thermal_ice_mechanism/
F5_spatial_recruitment/
```

Each figure folder contains its own `README.md`, `aggregate.py`, `plot.py`, and `data/` subdirectory.
