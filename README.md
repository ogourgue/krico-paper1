# KRICO: Paper 1

Figures and analysis for KRICO Paper 1: *Spawning phenology and upstream connectivity jointly determine larval recruitment success in Antarctic krill*. Publication-ready visualizations of 32-year hindcast (1994–2025) larval dispersal outcomes and spatial connectivity across CCAMLR Areas 48 and 88.

Author: Olivier Gourgue (RBINS)

Related repositories:

* __[krico-templates](https://github.com/ogourgue/krico-templates)__ — Simulation templates (Parcels + GLORYS12v1)
* __[krico-post-production](https://github.com/ogourgue/krico-post-production)__ — Post-processing and recruitment classification pipeline

## Dependencies

Reads recruitment outcome data produced by [krico-post-production](https://github.com/ogourgue/krico-post-production), specifically the `recruitment/` pipeline. Expected location: `$KRICO_ROOT/Post/Production/recruitment/data/`.

## Layout

```
F1_domain_map/                  # methods: domain, bathymetry zones, CCAMLR subareas
F2_phenology_curve/             # headline: 32-year mean recruitment success vs release date
F3_outcome_composition/         # 32-year mean fractional outcome breakdown vs release date
F4_source_maps/                 # release-position density per outcome (where particles came from)
F5_destination_maps/            # fate-position density per outcome (where particles ended up)
```

Each figure folder contains its own `README.md`, `aggregate.py`, `plot.py`, and `data/` subdirectory.