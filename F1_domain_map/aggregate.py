"""
F1 aggregation: mask GLORYS12 bathymetry within CCAMLR statistical areas, and
compute sea-ice concentration climatologies for two months and two periods.

The pixel-by-pixel point-in-polygon test takes ~1 minute on the full
GLORYS12 grid. Doing it once and writing the masked field to NetCDF
makes the plotting script fast for iterations.

The sea-ice climatologies are stored on the native NSIDC polar-stereographic
grid with their own 2-D lon/lat coordinates, so that plot.py can contour them
with the same PlateCarree transform used for the bathymetry zones. The 15%
contour that defines the ice edge is applied at plot time, not here.

Inputs:
  $KRICO_GLORYS12/glorys12_bathymetry.nc
  ccamlr-data/CCAMLR_ASD_EPSG4326.shp   (bundled in repo)
  data/sic/*.nc                          (downloaded by download_sic.sh; gitignored)

Output:
  data/aggregated.nc
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Point
from shapely.ops import unary_union


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CCAMLR_AREAS = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3"]

# Sea-ice climatology: February (southern-hemisphere minimum) and September
# (maximum), each split into the pre- and post-2016 regimes. The recent period
# ends in 2025 to match the span of simulated spawning years (1994-2025).
#
# September is not merely a generic winter month: it is the northernmost the
# ice reaches all year, and so determines whether winter ice arrives at all in
# 48.3 and 48.6N -- the M6 outcome the Discussion turns on.
SIC_MONTHS = {"February": 2, "September": 9}
SIC_PERIODS = {
    "1979-2015": (1979, 2015),
    "2016-2025": (2016, 2025),
}

# Candidate names for the concentration variable across G02202 versions.
SIC_VARIABLE_CANDIDATES = [
    "cdr_seaice_conc_monthly",
    "cdr_seaice_conc",
    "seaice_conc_monthly_cdr",
]

REPO_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

def display_path(path: Path) -> str:
    """
    Render a path for logging: repo-relative when it lies inside the repo,
    otherwise absolute.

    In-repo paths are short and identical on every machine, which is what makes
    a log comparable between runs. External inputs stay absolute, because a
    relative path would misrepresent where they actually live; the GLORYS12
    bathymetry is shown through its environment variable so the log also records
    how it was located.
    """
    path = Path(path).resolve()
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        pass
    krico_glorys12 = os.environ.get("KRICO_GLORYS12")
    if krico_glorys12:
        try:
            relative = path.relative_to(Path(krico_glorys12).resolve())
            return f"$KRICO_GLORYS12/{relative}"
        except ValueError:
            pass
    return str(path)


# ---------------------------------------------------------------------------
# Bathymetry aggregation
# ---------------------------------------------------------------------------

def aggregate_bathymetry(bathy_path: Path, ccamlr_shp_path: Path) -> xr.Dataset:
    """
    Load bathymetry and CCAMLR polygons, build a point-in-polygon mask
    for the project areas, and return the masked bathymetry.
    """
    print(f"Reading bathymetry: {display_path(bathy_path)}")
    bathy_ds = xr.open_dataset(bathy_path)
    hdept = bathy_ds["hdept"].isel(t=0).values.astype(np.float32)
    nav_lon = bathy_ds["nav_lon"].values.astype(np.float32)
    nav_lat = bathy_ds["nav_lat"].values.astype(np.float32)
    bathy_ds.close()

    print(f"Reading CCAMLR shapefile: {display_path(ccamlr_shp_path)}")
    ccamlr = gpd.read_file(ccamlr_shp_path)
    ccamlr = ccamlr[ccamlr["GAR_Long_L"].isin(CCAMLR_AREAS)]
    ccamlr_union = unary_union(ccamlr.geometry)

    print("Masking bathymetry pixels within CCAMLR areas (this takes a while)...")
    mask = np.zeros_like(hdept, dtype=bool)
    n_rows = nav_lon.shape[0]
    for i in range(n_rows):
        for j in range(nav_lon.shape[1]):
            if ccamlr_union.contains(Point(nav_lon[i, j], nav_lat[i, j])):
                mask[i, j] = True
        if (i + 1) % 100 == 0:
            print(f"  row {i + 1} / {n_rows}")
    print(f"  {int(mask.sum()):,} cells inside CCAMLR areas")

    masked_hdept = np.where(mask, hdept, np.nan).astype(np.float32)

    return xr.Dataset(
        data_vars={
            "bathymetry_masked": (
                ("y", "x"),
                masked_hdept,
                {
                    "long_name": "GLORYS12 bathymetry, NaN outside CCAMLR project areas",
                    "units": "m",
                    "positive": "down",
                },
            ),
            "nav_lon": (
                ("y", "x"),
                nav_lon,
                {"long_name": "longitude", "units": "degrees_east"},
            ),
            "nav_lat": (
                ("y", "x"),
                nav_lat,
                {"long_name": "latitude", "units": "degrees_north"},
            ),
        },
    )


# ---------------------------------------------------------------------------
# Sea-ice aggregation
# ---------------------------------------------------------------------------

def find_sic_variable(ds: xr.Dataset) -> str:
    """Locate the sea-ice concentration variable in a G02202 file."""
    for name in SIC_VARIABLE_CANDIDATES:
        if name in ds.variables:
            return name
    raise KeyError(
        "No recognised sea-ice concentration variable. "
        f"Variables present: {sorted(ds.data_vars)}"
    )


def find_sic_coordinates(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """
    Return 2-D longitude and latitude arrays for the G02202 grid.

    Recent versions ship `longitude` / `latitude` as 2-D coordinate variables.
    Older layouts keep them in the ancillary file only; in that case the
    projected `xgrid` / `ygrid` axes are converted with pyproj.
    """
    for lon_name, lat_name in (("longitude", "latitude"), ("lon", "lat")):
        if lon_name in ds.variables and lat_name in ds.variables:
            lon = np.asarray(ds[lon_name].values, dtype=np.float32)
            lat = np.asarray(ds[lat_name].values, dtype=np.float32)
            if lon.ndim == 2:
                return lon, lat

    print("  2-D lon/lat not found in file; deriving from projected axes")
    from pyproj import CRS, Transformer

    x_name = "xgrid" if "xgrid" in ds.variables else "x"
    y_name = "ygrid" if "ygrid" in ds.variables else "y"
    x2d, y2d = np.meshgrid(ds[x_name].values, ds[y_name].values)

    crs = None
    for name in ("crs", "projection", "spatial_ref"):
        if name in ds.variables:
            attrs = ds[name].attrs
            for key in ("proj4text", "spatial_ref", "crs_wkt"):
                if key in attrs:
                    crs = CRS.from_user_input(attrs[key])
                    break
        if crs is not None:
            break
    if crs is None:
        # NSIDC Sea Ice Polar Stereographic South.
        crs = CRS.from_user_input("EPSG:3412")
        print("  projection metadata absent; assuming EPSG:3412")

    transformer = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    lon, lat = transformer.transform(x2d, y2d)
    return lon.astype(np.float32), lat.astype(np.float32)


def to_fraction(da: xr.DataArray) -> xr.DataArray:
    """
    Mask flag values (land, missing, pole hole) and return concentration as a
    fraction in [0, 1].

    Two traps here, both silent if got wrong:

    `valid_range` is stored in the raw integer space of the file, while xarray
    has already applied `scale_factor` by the time this runs. G02202 packs
    concentration as unsigned bytes with a 0.01 scale, so an unscaled
    `valid_range` of [0, 100] would pass the 2.51-2.55 flag values straight
    through. It is scaled here to match the decoded data.

    Flags must also be removed before testing the magnitude of the data, and a
    *missing* `units` attribute must not be read as "fraction": a percent file
    without units would then skip conversion, and the whole field would fall
    below the 15% contour.
    """
    values = da.astype("float32")

    valid_range = da.attrs.get("valid_range")
    if valid_range is not None and len(valid_range) == 2:
        low, high = (float(v) for v in valid_range)
        scale = da.encoding.get("scale_factor")
        offset = da.encoding.get("add_offset", 0.0)
        if scale is not None:
            low = low * float(scale) + float(offset)
            high = high * float(scale) + float(offset)
            print(f"  valid_range scaled by {float(scale):g} to match decoded data")
        print(f"  masking to valid_range [{low:g}, {high:g}]")
        values = values.where((values >= low) & (values <= high))
    else:
        print("  WARNING: no valid_range attribute; flags may survive masking",
              file=sys.stderr)

    units = str(da.attrs.get("units", "")).strip().lower()
    if units in {"percent", "%"}:
        is_percent = True
    elif units in {"1", "fraction"}:
        is_percent = False
    else:
        is_percent = None

    vmax = float(values.max(skipna=True).compute())
    if is_percent is None:
        is_percent = vmax > 1.5
        label = "absent" if not units else f"unrecognised ({units!r})"
        print(f"  units attribute {label}; inferring "
              f"{'percent' if is_percent else 'fraction'} from maximum {vmax:.3f}")

    if is_percent:
        print(f"  converting from percent (maximum {vmax:.1f})")
        values = values / 100.0

    return values.where((values >= 0.0) & (values <= 1.0))


def aggregate_sea_ice(sic_dir: Path) -> xr.Dataset:
    """
    Compute the mean sea-ice concentration field for each month and period.

    Returns a dataset with dimensions (month, period, sic_y, sic_x) plus the
    2-D lon/lat coordinates of the native NSIDC grid.
    """
    files = sorted(sic_dir.glob("*.nc"))
    if not files:
        raise FileNotFoundError(
            f"No NetCDF files found in {display_path(sic_dir)}. "
            "Run ./download_sic.sh first."
        )
    print(f"Reading {len(files)} sea-ice files from {display_path(sic_dir)}")

    # download_sic.sh fetches only the wanted months, but the month filter below
    # is kept so the script stays correct if the full record is downloaded.

    # data_vars="minimal" concatenates only variables carrying the time
    # dimension; grid metadata is taken from the first file rather than
    # stacked 47 times. Set explicitly because the xarray default is changing.
    ds = xr.open_mfdataset(
        files,
        combine="by_coords",
        mask_and_scale=True,
        data_vars="minimal",
        coords="minimal",
        compat="override",
    )
    sic_name = find_sic_variable(ds)
    print(f"  concentration variable: {sic_name}")

    sic_lon, sic_lat = find_sic_coordinates(ds)

    concentration = to_fraction(ds[sic_name])

    months = list(SIC_MONTHS)
    periods = list(SIC_PERIODS)
    fields = np.empty((len(months), len(periods)), dtype=object)
    counts = np.zeros((len(months), len(periods)), dtype=np.int16)

    for i, month_name in enumerate(months):
        month_number = SIC_MONTHS[month_name]
        monthly = concentration.sel(
            time=concentration["time"].dt.month == month_number
        )
        for j, period in enumerate(periods):
            year_min, year_max = SIC_PERIODS[period]
            subset = monthly.sel(
                time=slice(f"{year_min}-01-01", f"{year_max}-12-31")
            )
            n_fields = int(subset.sizes["time"])
            expected = year_max - year_min + 1
            print(f"  {month_name} {period}: averaging {n_fields} fields")
            if n_fields != expected:
                print(
                    f"    WARNING: expected {expected} fields; "
                    f"check for gaps in the record",
                    file=sys.stderr,
                )
            fields[i, j] = (
                subset.mean("time", skipna=True).compute().values.astype(np.float32)
            )
            counts[i, j] = n_fields

    stacked = np.stack([np.stack(list(row)) for row in fields])

    ds.close()

    return xr.Dataset(
        data_vars={
            "sic_climatology": (
                ("month", "period", "sic_y", "sic_x"),
                stacked,
                {
                    "long_name": (
                        "Multi-year mean monthly sea-ice concentration, "
                        "NOAA/NSIDC G02202 v6"
                    ),
                    "units": "1",
                    "source_doi": "10.7265/b18j-z797",
                },
            ),
            "sic_n_fields": (
                ("month", "period"),
                counts,
                {"long_name": "number of monthly fields averaged"},
            ),
            "sic_lon": (
                ("sic_y", "sic_x"),
                sic_lon,
                {"long_name": "longitude", "units": "degrees_east"},
            ),
            "sic_lat": (
                ("sic_y", "sic_x"),
                sic_lat,
                {"long_name": "latitude", "units": "degrees_north"},
            ),
        },
        coords={
            "month": np.array(months, dtype=object),
            "period": np.array(periods, dtype=object),
        },
    )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    # Bathymetry: external dependency, set via KRICO_GLORYS12 env var.
    # Only needed to re-run F1 aggregation from scratch; data/aggregated.nc
    # is committed to the repo for figure reproduction without HPC access.
    krico_glorys12 = os.environ.get("KRICO_GLORYS12")
    if not krico_glorys12:
        print(
            "ERROR: KRICO_GLORYS12 environment variable not set.\n"
            "Set it to the GLORYS12 preprocessing output directory, e.g.:\n"
            "  export KRICO_GLORYS12=/path/to/Pre/GLORYS12",
            file=sys.stderr,
        )
        sys.exit(1)

    bathy_path = Path(krico_glorys12) / "glorys12_bathymetry.nc"
    if not bathy_path.exists():
        print(f"ERROR: bathymetry file not found: {display_path(bathy_path)}", file=sys.stderr)
        sys.exit(1)

    # CCAMLR shapefile: bundled in the repo at ccamlr-data/.
    ccamlr_shp_path = REPO_ROOT / "ccamlr-data" / "CCAMLR_ASD_EPSG4326.shp"
    if not ccamlr_shp_path.exists():
        print(f"ERROR: CCAMLR shapefile not found: {display_path(ccamlr_shp_path)}", file=sys.stderr)
        sys.exit(1)

    # Sea-ice concentration: downloaded into data/sic/ by download_sic.sh and
    # gitignored. Not needed to reproduce the figure -- data/aggregated.nc
    # carries the climatologies -- only to re-run this aggregation.
    sic_dir = Path(__file__).resolve().parent / "data" / "sic"
    if not sic_dir.is_dir():
        print(
            f"ERROR: sea-ice directory not found: {display_path(sic_dir)}\n"
            "Run ./download_sic.sh first.",
            file=sys.stderr,
        )
        sys.exit(1)

    ds_bathy = aggregate_bathymetry(bathy_path, ccamlr_shp_path)
    ds_sic = aggregate_sea_ice(sic_dir)

    ds = xr.merge([ds_bathy, ds_sic])
    ds.attrs = {
        "title": "F1 domain map aggregation - masked bathymetry and sea-ice climatology",
        "description": (
            "GLORYS12 bathymetry masked to KRICO project CCAMLR areas "
            "(48.1-48.6 and 88.3), and multi-year mean sea-ice concentration "
            "for two months and two periods, for use as input to the F1 "
            "domain map plot."
        ),
        "ccamlr_areas": ", ".join(CCAMLR_AREAS),
        "sic_source": (
            "Meier, W. N., Fetterer, F., Windnagel, A. K., Stewart, J. S. & "
            "Stafford, T. (2026). NOAA/NSIDC Climate Data Record of Passive "
            "Microwave Sea Ice Concentration (G02202, Version 6). NSIDC. "
            "https://doi.org/10.7265/b18j-z797"
        ),
    }

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "aggregated.nc"

    print(f"Writing {display_path(out_path)}")
    encoding = {
        var: {"zlib": True, "complevel": 5}
        for var in ds.data_vars
        if np.issubdtype(ds[var].dtype, np.number)
    }
    ds.to_netcdf(out_path, encoding=encoding)

    print("Done.")


if __name__ == "__main__":
    main()