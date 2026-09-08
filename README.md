# KRICO: Paper 1

Figures and analysis for KRICO Paper 1: *Timing of sea-ice retreat and advance governs the January peak and regional asymmetries in Antarctic krill larval recruitment*. Publication-ready visualizations of 32-year hindcast (1994–2025) larval dispersal outcomes and spatial connectivity across CCAMLR Areas 48 and 88.

Author: Olivier Gourgue (RBINS)

Related repositories:

* __[krico-templates](https://github.com/ogourgue/krico-templates)__ — Simulation templates (Parcels + GLORYS12v1)
* __[krico-post-production](https://github.com/ogourgue/krico-post-production)__ — Post-processing and recruitment classification pipeline

## Folder naming

Three kinds of folder, distinguished by prefix:

- `F*` — main-text figures (F1–F5).
- `FS*` — Supporting Information figures (FS1, FS2). These behave like the `F*` folders: each produces a rendered figure.
- `S*` — Supporting Information analyses (S1, S2). These produce no figure. Their requirements differ: S1 reads the raw trajectories and imports the recruitment package rather than only its outputs, so it is heavier than any figure folder; S2 reads only the recruitment outputs, like the `F*` aggregations.

Each folder maps to one or two items in the manuscript:

| Folder | Manuscript item |
|---|---|
| `F1_domain_map` | Figure 1 |
| `F2_phenology_curve` | Figure 2 |
| `F3_outcome_composition` | Figure 3 |
| `F4_source_maps` | Figure 4, Table S1 |
| `F5_destination_maps` | Figure 5, Table S2 |
| `FS1_phenology_per_year` | Figure S1 |
| `FS2_currents` | Figure S2 |
| `S1_m1_offset_sensitivity` | Text S1 |
| `S2_m6_window_sensitivity` | Text S2 |

AGU numbers the Text, Figure and Table series separately in the Supporting Information, so a shared number does not mean a shared item: `S1` is Text S1 and `FS1` is Figure S1, and the two are unrelated. Tables S1 and S2 are produced by `ccamlr_summary.py` in the F4 and F5 folders rather than by a folder of their own.

## Environment

Most of this repository runs locally on Python 3.14.7. Two workloads run on the ECMWF HPC instead, on Python 3.13.13: `FS2_currents/aggregate.py`, submitted via its `job.sh`, and `S1_m1_offset_sensitivity/`, submitted via its `sweep.sh`. Neither environment requires Parcels — every folder here reads simulation output rather than producing it.

The split follows the environment variables described below: folders needing `KRICO_GLORYS12`, `KRICO_ROOT` or `KRICO_RUNS` run on the HPC, everything else runs locally. In particular `S2_m6_window_sensitivity/` runs on a laptop despite reading the same recruitment outputs as S1.

Dependencies are `numpy`, `xarray`, `pandas`, `netCDF4`, `matplotlib`, `cartopy`, `geopandas` and `shapely`. S1 additionally imports the `krico_recruitment` package from [krico-post-production](https://github.com/ogourgue/krico-post-production), which is why its `KRICO_POST` must point at a clone of that repository rather than at a directory of outputs.

## Reproducing the figures

The aggregated NetCDF files (`F*/data/aggregated.nc`) and the rendered figures (`F*/*.png`, `FS*/*.png`) are committed to the repo, so reproducing each figure requires only:

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

The recruitment outcome data is archived on Zenodo (DOI: [10.5281/zenodo.22548763](https://doi.org/10.5281/zenodo.22548763)); see [krico-post-production](https://github.com/ogourgue/krico-post-production) for the download workflow. That is version 2.0.0 of the dataset, the M1-at-spawning classification these figures are built on; version 1.0.0 evaluated M1 at release and will not reproduce them. The GLORYS12 bathymetry and the reformatted velocity fields are produced by the preprocessing step in [krico-templates](https://github.com/ogourgue/krico-templates).

F1 additionally needs observed sea-ice concentration, which is not committed (~22 MB). Run `F1_domain_map/download_sic.sh` first; it fetches the February and September monthly files of the NOAA/NSIDC Climate Data Record (G02202 v6, DOI: [10.7265/b18j-z797](https://doi.org/10.7265/b18j-z797)) into `F1_domain_map/data/sic/` (gitignored) over anonymous HTTPS, and records DOI, source URL, access date and per-file checksums in `F1_domain_map/sic_provenance.txt` (committed).

FS2 is a batch job rather than an interactive one: submit `FS2_currents/job.sh` rather than running its `aggregate.py` on a login node. It checkpoints every 12 month-pairs and resumes, so an interrupted run does not start over.

After each aggregation run, re-run the corresponding `plot.py` to regenerate the figure.

## Supporting Information analyses

`S*` folders hold analyses that support the Supporting Information rather than a figure.

- `S1_m1_offset_sensitivity/` — how much the M1 classification depends on the descent-ascent offset between spawning and calyptopis I.
- `S2_m6_window_sensitivity/` — how much the seasonal decline in M6 depends on the advance-detection window, which lengthens through the release season.

The two differ in what they need. In addition to `KRICO_POST` and `KRICO_GLORYS12` above, S1 requires:

- `KRICO_RUNS` — path to the raw trajectory simulations (e.g., `/path/to/KRICO/Runs`). The raw trajectories are not publicly archived; see [krico-post-production](https://github.com/ogourgue/krico-post-production).

Note that S1 uses `KRICO_POST` to import `krico_recruitment`, so it must point at a clone of the pipeline repository, not merely at a directory of recruitment outputs. It also reads the GLORYS12 monthly sea-ice files (`glorys12_ice_YYYY_MM.nc`) from `KRICO_GLORYS12`, alongside the bathymetry that F1 uses.

S2 needs only `KRICO_POST`, and reads recruitment outputs from `$KRICO_POST/recruitment/data` exactly as the `F*` aggregations do. Its second test reads the committed `F1_domain_map/data/aggregated.nc` for the sea-ice edge, so nothing is downloaded. It runs on a laptop.

See each folder's own README for what the analysis establishes and how to run it.

## Layout

```
F1_domain_map/                  # methods: domain, bathymetry zones, CCAMLR subareas, sea-ice edges
F2_phenology_curve/             # headline: 32-year mean recruitment success vs release date
F3_outcome_composition/         # 32-year mean fractional outcome breakdown vs release date
F4_source_maps/                 # release-position density per outcome (where particles came from)
F5_destination_maps/            # fate-position density per outcome (where particles ended up)
FS1_phenology_per_year/         # SI figure: per-year success curves and the distribution of peak dates
FS2_currents/                   # SI figure: 32-year mean circulation over the 50-200 m release depth band
S1_m1_offset_sensitivity/       # SI analysis: sensitivity of M1 to the descent-ascent offset
S2_m6_window_sensitivity/       # SI analysis: sensitivity of M6 to the advance-detection window
ccamlr-data/                    # CCAMLR Statistical Areas shapefile (bundled; see ccamlr-data/README.md)
```

Each figure folder contains its own `README.md`, `aggregate.py`, `plot.py`, `data/aggregated.nc`, and the rendered PNG(s), with two exceptions: FS1 has no `aggregate.py` or `data/` (it reads F2's aggregation), and FS2's `data/aggregated.nc` is gitignored rather than committed. F1 additionally contains `download_sic.sh` and `sic_provenance.txt`; F4 and F5 additionally contain `ccamlr_summary.py` and `ccamlr_summary.csv`; FS2 additionally contains `job.sh`. The `S*` folders are not figure folders: S1 contains `sweep.py`, `sweep.sh`, `compare.py` and its output CSV; S2 contains `sweep.py`, `compare.py`, `ice_reachability.py` and two output CSVs.
