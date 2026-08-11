# KRICO: Paper 1

Figures and analysis for KRICO Paper 1: *Sea-ice timing governs the January peak and regional asymmetries in Antarctic krill larval recruitment*. Publication-ready visualizations of 32-year hindcast (1994–2025) larval dispersal outcomes and spatial connectivity across CCAMLR Areas 48 and 88.

Author: Olivier Gourgue (RBINS)

Related repositories:

* __[krico-templates](https://github.com/ogourgue/krico-templates)__ — Simulation templates (Parcels + GLORYS12v1)
* __[krico-post-production](https://github.com/ogourgue/krico-post-production)__ — Post-processing and recruitment classification pipeline

## Reproducing the figures

The aggregated NetCDF files (`F*/data/aggregated.nc`) and the rendered figures (`F*/*.png`) are committed to the repo, so reproducing each figure requires only:

```bash
cd F1_domain_map
python plot.py
```

Repeat in `F2_phenology_curve`, `F3_outcome_composition`, `F4_source_maps`, and `F5_destination_maps`. The CCAMLR Statistical Areas shapefile is bundled at `ccamlr-data/`. No environment variables or external data are needed for this workflow.

## Re-running the aggregation from scratch

To regenerate the `aggregated.nc` files from upstream sources, two environment variables are required:

- `KRICO_POST` — path to the krico-post-production root directory (e.g., `/path/to/krico-post-production`). Used by F2, F3, F4, F5 `aggregate.py`; recruitment outcome data is read from `$KRICO_POST/recruitment/data`.
- `KRICO_GLORYS12` — path to the GLORYS12 preprocessing output directory containing `glorys12_bathymetry.nc` (e.g., `/path/to/Pre/GLORYS12`). Used by F1 `aggregate.py` only.

The recruitment outcome data is archived on Zenodo (DOI: [10.5281/zenodo.20101159](https://doi.org/10.5281/zenodo.20101159)); see [krico-post-production](https://github.com/ogourgue/krico-post-production) for the download workflow. The GLORYS12 bathymetry is produced by the preprocessing step in [krico-templates](https://github.com/ogourgue/krico-templates).

F1 additionally needs observed sea-ice concentration, which is not committed (~22 MB). Run `F1_domain_map/download_sic.sh` first; it fetches the February and September monthly files of the NOAA/NSIDC Climate Data Record (G02202 v6, DOI: [10.7265/b18j-z797](https://doi.org/10.7265/b18j-z797)) into `F1_domain_map/data/sic/` (gitignored) over anonymous HTTPS, and records DOI, source URL, access date and per-file checksums in `F1_domain_map/sic_provenance.txt` (committed).

After each aggregation run, re-run the corresponding `plot.py` to regenerate the figure.

## Layout

```
F1_domain_map/                  # methods: domain, bathymetry zones, CCAMLR subareas, sea-ice edges
F2_phenology_curve/             # headline: 32-year mean recruitment success vs release date
F3_outcome_composition/         # 32-year mean fractional outcome breakdown vs release date
F4_source_maps/                 # release-position density per outcome (where particles came from)
F5_destination_maps/            # fate-position density per outcome (where particles ended up)
ccamlr-data/                    # CCAMLR Statistical Areas shapefile (bundled; see ccamlr-data/README.md)
```

Each figure folder contains its own `README.md`, `aggregate.py`, `plot.py`, `data/aggregated.nc`, and the rendered PNG(s). F1 additionally contains `download_sic.sh` and `sic_provenance.txt`; F4 and F5 additionally contain `ccamlr_summary.py` and `data/ccamlr_summary.csv`.