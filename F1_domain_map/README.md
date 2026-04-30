# F1 — Domain map

Method figure showing the computational domain, the bathymetry-defined spawning and continental-shelf zones, and the CCAMLR Statistical Subareas covered by the simulation.

## Scope

- **Computational domain:** 115°W to 40°E, 78°S to 40°S.
- **Spawning zone:** bathymetry 1000–2000 m, masked to CCAMLR project subareas.
- **Continental shelf:** bathymetry < 1000 m, masked to CCAMLR project subareas.
- **CCAMLR Statistical Subareas shown:** 48.1, 48.2, 48.3, 48.4, 48.5, 48.6, 88.3. Subarea 48.6 spans a wide latitude range from the Antarctic shelf to the open Southern Ocean and is reported as two halves (48.6N, 48.6S) split at 60°S in the F4 / F5 spatial analyses; F1 visualizes this split. See `../METHODOLOGY.md` §0 for the rationale.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

The "continental shelf" naming (rather than "shelf recruitment zone") avoids conflation with the < 2000 m recruitment criterion used in M5b. The recruitment habitat is the combined shelf + slope (< 2000 m); F1 shows the shelf component alone for geographic context. See `../METHODOLOGY.md` §4 for the bathymetric vocabulary.

## Design

- **Bathymetry zones:**
  - Spawning zone: matplotlib `C0` (blue) at 50% opacity.
  - Continental shelf: matplotlib `C1` (orange) at 50% opacity.
- **Boundaries (consistent with F4 / F5):**
  - CCAMLR subarea outlines: plain light-gray (`"0.7"`).
  - Model computational domain: dashed light-gray (`"0.7"`) — external boundary.
  - 48.6 split at 60°S: dotted light-gray (`"0.7"`) — internal subdivision. Densified to follow the curved parallel on the stereographic projection, then clipped to the 48.6 polygon.
- **CCAMLR subarea labels:** centered on polygon centroids in framed semi-transparent boxes (style matching panel labels in F4 / F5; the underlying polygon and bathymetry remain visible through the label frame). Subarea 48.6 gets two labels (48.6N, 48.6S) at the centroids of its halves north and south of 60°S.
- **Two legends:**
  - Upper right, opaque: zone / domain key (3 entries — domain, spawning zone, continental shelf).
  - Lower left, opaque: CCAMLR subarea code → geographic name lookup (8 entries — 48.1 Antarctic Peninsula, 48.2 South Orkney Islands, 48.3 South Georgia, 48.4 South Sandwich Islands, 48.5 Weddell Sea, 48.6N Northern Bouvet, 48.6S Southern Bouvet, 88.3 Amundsen Sea).
- Coastlines and gray land fill via cartopy NaturalEarth (50 m resolution).
- Axis spines pushed above land/coastlines so the panel border is never covered.

## Linestyle convention

In F1 (and consistent with the rest of the manuscript figures), gray lines distinguish boundary types:

- **Dashed gray:** external boundary (computational domain).
- **Dotted gray:** internal subdivision (60°S split within 48.6).
- **Solid gray:** CCAMLR subarea outlines.

## Files

- `aggregate.py` — reads the GLORYS12 bathymetry NetCDF and the CCAMLR shapefile, masks bathymetry pixels to the project subareas via point-in-polygon, writes `data/aggregated.nc`. Slow step (~1 minute).
- `plot.py` — reads `data/aggregated.nc`, re-loads the CCAMLR shapefile for polygon outlines and centroids, splits 48.6 at 60°S for label placement, writes `domain_map.png`. Fast.

## Inputs (read by `aggregate.py`)

- `$KRICO_ROOT/Pre/glorys12/glorys12_bathymetry.nc`
- `$KRICO_ROOT/Pre/ccamlr-data/geographical_data/asd/CCAMLR_ASD_EPSG4326.shp`

## Run

```bash
python aggregate.py
python plot.py
```