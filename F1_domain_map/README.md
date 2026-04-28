# F1 — Domain map

Method figure showing the computational domain, the bathymetry-defined release and recruitment zones, and the CCAMLR statistical areas covered by the simulation.

## Scope

- **Computational domain:** 115°W to 40°E, 78°S to 40°S.
- **Spawning zone:** bathymetry 1000–2000 m, masked to CCAMLR project areas.
- **Shelf recruitment zone:** bathymetry < 1000 m, masked to CCAMLR project areas.
- **CCAMLR Subareas shown:** 48.1, 48.2, 48.3, 48.4, 48.5, 48.6, 88.3.
- Projection: South Polar Stereographic, central_longitude = −37.5°.

## Design

- **Bathymetry zones:**
  - Spawning zone: matplotlib `C1` (orange) at 50% opacity.
  - Shelf recruitment zone: matplotlib `C2` (green) at 50% opacity.
- **Boundaries (consistent with F4 / F5):**
  - CCAMLR area outlines: plain light-gray (`"0.7"`).
  - Model computational domain: dashed light-gray (`"0.7"`).
- **CCAMLR area labels:** centered on polygon centroids in framed boxes (style matching panel labels in F4 / F5).
- **Legend:** upper-right, opaque, three entries (domain, spawning zone, shelf zone).
- Coastlines and gray land fill via cartopy NaturalEarth (50 m resolution).
- Axis spines pushed above land/coastlines so the panel border is never covered.

## Files

- `aggregate.py` — reads the GLORYS12 bathymetry NetCDF and the CCAMLR shapefile, masks bathymetry pixels to the project areas via point-in-polygon, writes `data/aggregated.nc`. Slow step (~1 minute).
- `plot.py` — reads `data/aggregated.nc`, re-loads the CCAMLR shapefile for polygon outlines and centroids, writes `domain_map.png`. Fast.

## Inputs (read by `aggregate.py`)

- `$KRICO_ROOT/Pre/glorys12/glorys12_bathymetry.nc`
- `$KRICO_ROOT/Pre/ccamlr-data/geographical_data/asd/CCAMLR_ASD_EPSG4326.shp`

## Run

```bash
python aggregate.py
python plot.py
```
