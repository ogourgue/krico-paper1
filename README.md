# KRICO: Paper 1

Figures and analysis for KRICO Paper 1: *Timing of sea-ice retreat and advance governs the January peak and regional asymmetries in Antarctic krill larval recruitment*. Publication-ready visualizations of 32-year hindcast (1994–2025) larval dispersal outcomes and spatial connectivity across CCAMLR Areas 48 and 88.

Author: Olivier Gourgue (RBINS)

Related repositories:

* __[krico-templates](https://github.com/ogourgue/krico-templates)__ — Simulation templates (Parcels + GLORYS12v1)
* __[krico-post-production](https://github.com/ogourgue/krico-post-production)__ — Post-processing and recruitment classification pipeline

## Reproducing the figures

The aggregated NetCDF files (`*/data/aggregated.nc`) and the rendered figures (`*/*.png`) are committed to the repo, so reproducing each figure requires only:

```bash
cd F1_domain_map
python plot.py
```

Repeat in `F2_phenology_curve`, `F3_outcome_composition`, `F4_source_maps`, `F5_destination_maps`, and `FS1_phenology_per_year`. The CCAMLR Statistical Areas shapefile is bundled at `ccamlr-data/`. No environment variables or external data are needed for this workflow.

Two folders depart from that pattern:

- **`FS1_phenology_per_year`** has no `aggregate.py` and no `data/`. It needs exactly the quantities F2 already computes, so its `plot.py` reads `../F2_phenology_curve/data/aggregated.nc` directly. Duplicating the aggregation would let the two figures disagree about the same numbers.
- **`FS2_currents`** is the one figure that cannot be reproduced from this repo alone. Its `data/aggregated.nc` is a mean field on the full 1/12° grid (~58 MB on disk), too large to track, so it is gitignored. Regenerating it means re-running `FS2_currents/aggregate.py` on a machine holding the GLORYS12 velocity fields — ~1.5 TB of reads over 384 monthly files, ~4.5 h as a batch job — and copying the result back. See `FS2_currents/README.md`.

## Re-running the aggregation from scratch

To regenerate the `aggregated.nc` files from upstream sources, three environment variables are required:

- `KRICO_POST` — path to the krico-post-production root directory (e.g., `/path/to/krico-post-production`). Used by F2, F3, F4, F5 `aggregate.py`; recruitment outcome data is read from `$KRICO_POST/recruitment/data`.
- `KRICO_GLORYS12` — path to the GLORYS12 preprocessing output directory containing `glorys12_bathymetry.nc` (e.g., `/path/to/Pre/GLORYS12`). Used by F1 `aggregate.py` only.
- `KRICO_ROOT` — path to the KRICO project root, from which `Pre/GLORYS12/` is derived. Used by FS2 `aggregate.py` only, to read the reformatted velocity fields (`glorys12_u_YYYY_MM.nc`, `glorys12_v_YYYY_MM.nc`).

The recruitment outcome data is archived on Zenodo (DOI: [10.5281/zenodo.20101159](https://doi.org/10.5281/zenodo.20101159)); see [krico-post-production](https://github.com/ogourgue/krico-post-production) for the download workflow. The GLORYS12 bathymetry and the reformatted velocity fields are produced by the preprocessing step in [krico-templates](https://github.com/ogourgue/krico-templates).

F1 additionally needs observed sea-ice concentration, which is not committed (~22 MB). Run `F1_domain_map/download_sic.sh` first; it fetches the February and September monthly files of the NOAA/NSIDC Climate Data Record (G02202 v6, DOI: [10.7265/b18j-z797](https://doi.org/10.7265/b18j-z797)) into `F1_domain_map/data/sic/` (gitignored) over anonymous HTTPS, and records DOI, source URL, access date and per-file checksums in `F1_domain_map/sic_provenance.txt` (committed).

FS2 is a batch job rather than an interactive one: submit `FS2_currents/job.sh` rather than running its `aggregate.py` on a login node. It checkpoints every 12 month-pairs and resumes, so an interrupted run does not start over.

After each aggregation run, re-run the corresponding `plot.py` to regenerate the figure.

## Layout

```
F1_domain_map/                  # methods: domain, bathymetry zones, CCAMLR subareas, sea-ice edges
F2_phenology_curve/             # headline: 32-year mean recruitment success vs release date
F3_outcome_composition/         # 32-year mean fractional outcome breakdown vs release date
F4_source_maps/                 # release-position density per outcome (where particles came from)
F5_destination_maps/            # fate-position density per outcome (where particles ended up)
FS1_phenology_per_year/         # SI: per-year success curves and the distribution of peak dates
FS2_currents/                   # SI: 32-year mean circulation over the 50-200 m release depth band
ccamlr-data/                    # CCAMLR Statistical Areas shapefile (bundled; see ccamlr-data/README.md)
```

Each figure folder contains its own `README.md`, `aggregate.py`, `plot.py`, `data/aggregated.nc`, and the rendered PNG(s), with three exceptions: FS1 has no `aggregate.py` or `data/` (it reads F2's aggregation), and FS2's `data/aggregated.nc` is gitignored rather than committed. F1 additionally contains `download_sic.sh` and `sic_provenance.txt`; F4 and F5 additionally contain `ccamlr_summary.py` and `data/ccamlr_summary.csv`; FS2 additionally contains `job.sh`.
