# KRICO Paper 1

Figures and analysis for KRICO Paper 1: *Spawning phenology and upstream connectivity jointly determine larval recruitment success in Antarctic krill*.

Each figure lives in its own self-contained subfolder with aggregation and plotting scripts.

## Dependencies

Reads recruitment outcome data produced by [krico-post-production](https://github.com/ogourgue/krico-post-production), specifically the `recruitment/` pipeline. Expected location: `$KRICO_ROOT/Post/Production/recruitment/data/`.

## Layout

```
F1_domain_map/                  # methods: domain, bathymetry zones, CCAMLR areas
F2_phenology_curve/             # headline: 32-year mean recruitment success vs release date
F3_outcome_composition/         # 32-year mean fractional outcome breakdown vs release date
F4_source_maps/                 # release-position density per outcome (where particles came from)
F5_destination_maps/            # fate-position density per outcome (where particles ended up)
```

Each figure folder contains its own `README.md`, `aggregate.py`, `plot.py`, and `data/` subdirectory.