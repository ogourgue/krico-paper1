# F1 — Domain map

Method figure showing the computational domain, the bathymetry-defined spawning and continental-shelf zones, the CCAMLR Statistical Subareas covered by the simulation, and mean sea-ice edges for two months and two periods.

## Scope

- **Computational domain:** 115°W to 40°E, 78°S to 40°S.
- **Spawning zone:** bathymetry 1000–2000 m, masked to CCAMLR project subareas.
- **Continental shelf:** bathymetry < 1000 m, masked to CCAMLR project subareas.
- **CCAMLR Statistical Subareas shown:** 48.1, 48.2, 48.3, 48.4, 48.5, 48.6, 88.3. Subarea 48.6 spans a wide latitude range from the Antarctic shelf to the open Southern Ocean and is reported as two halves (48.6N, 48.6S) split at 60°S in the F4 / F5 spatial analyses; F1 visualizes this split.
- **Sea-ice edges:** 15% contour of the multi-year mean monthly sea-ice concentration, for February and September, each for 1979–2015 and 2016–2025 (the recent period ending in 2025 to match the simulated spawning years).
- Projection: South Polar Stereographic, central_longitude = −37.5°.

The "continental shelf" naming (rather than "shelf recruitment zone") avoids conflation with the < 2000 m recruitment criterion used in M5b. The recruitment habitat is the combined shelf + slope (< 2000 m); F1 shows the shelf component alone for geographic context.

The sea-ice edges are illustrative context for the post-2016 low-extent regime referred to in the introduction; no quantity derived from them enters the results. They are observations, from NOAA/NSIDC G02202 Version 6 (Meier et al., 2026; [doi:10.7265/b18j-z797](https://doi.org/10.7265/b18j-z797)) — chosen as a single continuous 1978-present record under one DOI, with the SSMIS → AMSR2 transition of January 2025 handled inside the product. The simulation's own sea-ice criteria are applied to GLORYS12 concentration, not to this product, and the caption states so. February is the ice minimum, September the maximum — and also the northernmost the ice reaches all year, so it determines whether winter ice arrives at all in 48.3 and 48.6N, the M6 outcome the Discussion turns on. Each mean concentration field is computed on the native NSIDC grid and the 15% contour taken from that mean (not as the median of annual contours); the threshold matches per-particle sea-ice advance detection.

## Design

- **Bathymetry zones:**
  - Spawning zone: matplotlib `C0` (blue) at 50% opacity.
  - Continental shelf: matplotlib `C1` (orange) at 50% opacity.
- **Boundaries (consistent with F4 / F5):**
  - CCAMLR subarea outlines: plain light-gray (`"0.7"`).
  - Model computational domain: dashed light-gray (`"0.7"`) — external boundary.
  - 48.6 split at 60°S: dotted light-gray (`"0.7"`) — internal subdivision. Densified to follow the curved parallel on the stereographic projection, then clipped to the 48.6 polygon.
- **Sea-ice edges:** colour for the period (`C2` green 1979–2015, `C3` red 2016–2025), linestyle for the month (solid September, dashed February), at `axes.linewidth` and alpha 0.7, above the bathymetry fills but below the boundaries. Masked to the domain before contouring, so each edge terminates at the domain boundary rather than trailing off into the East Antarctic sector.
- **CCAMLR subarea labels:** centered on polygon centroids in framed semi-transparent boxes (style matching panel labels in F4 / F5), 48.6 getting one per half. `CCAMLR_LABEL_ANCHORS` overrides placement where a centroid would obscure data; 48.1 is anchored on its northern boundary.
- **Three legends,** all left-aligned and opaque, one per free corner:
  - Upper right: zone / domain key (3 entries — domain, spawning zone, continental shelf).
  - Lower right, titled *Sea-ice edge:* — the four contours, ordered north to south as on the map.
  - Lower left: CCAMLR subarea code → geographic name lookup (8 entries — 48.1 Antarctic Peninsula, 48.2 South Orkney Islands, 48.3 South Georgia, 48.4 South Sandwich Islands, 48.5 Weddell Sea, 48.6N Northern Bouvet, 48.6S Southern Bouvet, 88.3 Amundsen Sea).
- Coastlines and gray land fill via cartopy NaturalEarth (50 m resolution).
- Axis spines pushed above land/coastlines so the panel border is never covered.

## Linestyle convention

In F1 (and consistent with the rest of the manuscript figures), gray lines distinguish boundary types:

- **Dashed gray:** external boundary (computational domain).
- **Dotted gray:** internal subdivision (60°S split within 48.6).
- **Solid gray:** CCAMLR subarea outlines.

Gray is reserved for geography, colour for data, solid black for the coastline.

## Files

- `download_sic.sh` — fetches the 94 February and September G02202 v6 files into `data/sic/` (gitignored) and writes `sic_provenance.txt` (DOI, URL, access date, per-file sha256; committed). Only needed to re-run aggregation.
- `aggregate.py` — reads the GLORYS12 bathymetry NetCDF and the CCAMLR shapefile, masks bathymetry pixels to the project subareas via point-in-polygon, computes the four sea-ice climatologies, writes `data/aggregated.nc`. Slow step (~1 minute).
- `plot.py` — reads `data/aggregated.nc`, re-loads the CCAMLR shapefile for polygon outlines and centroids, splits 48.6 at 60°S for label placement, writes `domain_map.png`. Fast.

## Inputs (read by `aggregate.py`)

- `$KRICO_GLORYS12/glorys12_bathymetry.nc`
- `ccamlr-data/CCAMLR_ASD_EPSG4326.shp` (bundled in repo)
- `data/sic/*.nc` (downloaded by `download_sic.sh`; gitignored)

## Run

```bash
./download_sic.sh
python aggregate.py
python plot.py
```