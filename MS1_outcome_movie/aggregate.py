"""
MS1 aggregation: one cohort, followed day by day, for Movie S1.

The movie follows a single release-day cohort from release to the end of
tracking. Every particle drifts in gray until its outcome is decided, then
turns green (recruitment success) or red (any mortality outcome) at its fate
position and stops moving. This script extracts what the frames need, so that
plot.py can render them on a laptop:

  - daily positions of every particle, up to and including its fate day
    (positions after the fate day are written as fill values: the particle
    no longer moves, and the fill values compress to almost nothing),
  - per-particle outcome, fate day and fate position, read from the archived
    recruitment classification rather than recomputed,
  - daily GLORYS12 sea-ice concentration over the domain for the tracking
    window, block-averaged to a coarser grid.

Cohort choice
-------------
By default the cohort is chosen from F2's aggregation, so the choice is
reproducible rather than hand-picked:

  1. the climatological peak release date (the max of the 32-year mean curve,
     January 29, as annotated in Figure 2);
  2. among spawning years up to LAST_YEAR_BEFORE_REGIME, the year whose own
     peak release date is closest to it;
  3. ties broken by the peak value closest to the median of per-year peaks.

The cohort is that year's own peak release date. Spawning years after 2016
are excluded so the movie cannot be read as a statement about the post-2016
low-ice regime. --date overrides the choice.

Particles are matched between the trajectory file and the recruitment file by
position in the file: process_cohort.py writes one row per trajectory, in
trajectory order. This is checked rather than assumed, by comparing the
trajectory position on each particle's fate day with the archived fate
position.

Inputs
------
  $KRICO_RUNS/KRICO_*/YYYY_MM_DD.nc                raw trajectories
  $KRICO_POST/recruitment/data/YYYY_MM_DD.nc       recruitment outcomes (v2.0.0)
  $KRICO_GLORYS12/glorys12_ice_YYYY_MM.nc          daily sea-ice concentration
  ../F2_phenology_curve/data/aggregated.nc         cohort choice

Output
------
  data/aggregated.nc     (gitignored; a few hundred MB at full resolution)

Usage
-----
    python aggregate.py                      # cohort chosen from F2
    python aggregate.py --date 2006-01-29    # explicit cohort
    python aggregate.py --stride 4           # every 4th particle, for tests

Submit via job.sh on the HPC rather than running on a login node.
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
F2_AGGREGATED = HERE.parent / "F2_phenology_curve" / "data" / "aggregated.nc"
OUT_PATH = HERE / "data" / "aggregated.nc"

# Domain extent (matches F1 / F4 / F5).
LON_MIN, LON_MAX = -115.0, 40.0
LAT_MIN, LAT_MAX = -78.0, -40.0

# Last spawning year eligible for the default cohort choice. The low-extent
# regime begins in late 2016, after the 2016 spawning season (Nov 2015 to
# Mar 2016) has been tracked.
LAST_YEAR_BEFORE_REGIME = 2016

# Descent-ascent offset between spawning and release (calyptopis I). Matches
# krico_recruitment.sea_ice.SPAWNING_OFFSET_DAYS; only used for the subtitle.
SPAWNING_OFFSET_DAYS = 24

# Sea ice is block-averaged by this factor on each axis (1/12 deg -> 1/4 deg).
# The frames are 1080 px wide for a ~155 deg wide domain, so a 1/4 deg cell is
# still about one pixel; the full grid would only multiply the file size.
SIC_COARSEN = 3

# Positions are packed as int16 with this scale (degrees). 0.005 deg is
# ~0.5 km in latitude, two orders of magnitude below a pixel, and the int16
# range then spans +-163 deg, which covers the domain.
POSITION_SCALE = np.float32(0.005)
POSITION_FILL = np.int16(-32768)

# Sea-ice concentration is stored as integer percent; 255 marks land.
SIC_FILL = np.uint8(255)

GLORYS_ICE_PATTERN = "glorys12_ice_{year:04d}_{month:02d}.nc"


# ---------------------------------------------------------------------------
# Cohort choice
# ---------------------------------------------------------------------------

def season_day_to_date(spawning_year: int, season_day: int) -> pd.Timestamp:
    """
    Calendar date of an F2 season-day index.

    F2 counts days from Nov 15 of the year before the spawning year, with
    February always 28 days long (Feb 29 cohorts are excluded). In leap years
    every index after Feb 28 is therefore one calendar day short.
    """
    start = pd.Timestamp(year=spawning_year - 1, month=11, day=15)
    date = start + pd.Timedelta(days=int(season_day))
    feb28 = pd.Timestamp(year=spawning_year, month=2, day=28)
    if pd.Timestamp(year=spawning_year, month=1, day=1).is_leap_year and date > feb28:
        date += pd.Timedelta(days=1)
    return date


def choose_cohort() -> tuple[pd.Timestamp, str]:
    """Pick the cohort from F2's aggregation; return (date, rationale)."""
    ds = xr.open_dataset(F2_AGGREGATED)
    rate = ds["success_rate"].where(ds["has_data"]).values * 100.0
    years = ds["year"].values

    # The last season-day slot (Mar 15) is never populated, hence the
    # all-NaN column warning silenced here.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        clim_peak = int(np.nanargmax(np.nanmean(rate, axis=0)))

    peak_day = np.full(years.size, -1)
    peak_value = np.full(years.size, np.nan)
    for i in range(years.size):
        if not np.all(np.isnan(rate[i])):
            peak_day[i] = int(np.nanargmax(rate[i]))
            peak_value[i] = float(np.nanmax(rate[i]))

    eligible = (peak_day >= 0) & (years <= LAST_YEAR_BEFORE_REGIME)
    median_peak = float(np.nanmedian(peak_value[peak_day >= 0]))

    # Lexicographic: distance to the climatological peak first, then peak
    # value distance to the median per-year peak.
    idx = np.flatnonzero(eligible)
    order = np.lexsort((np.abs(peak_value[idx] - median_peak),
                        np.abs(peak_day[idx] - clim_peak)))
    best = idx[order[0]]

    year = int(years[best])
    date = season_day_to_date(year, peak_day[best])
    clim_date = season_day_to_date(year, clim_peak)
    rationale = (
        f"spawning year {year}: own peak on {date:%b %d} at "
        f"{peak_value[best]:.2f}%; climatological peak {clim_date:%b %d}; "
        f"median per-year peak {median_peak:.2f}%; "
        f"years after {LAST_YEAR_BEFORE_REGIME} excluded"
    )
    return date, rationale


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def require_dir(var: str) -> Path:
    value = os.environ.get(var)
    if not value:
        sys.exit(f"ERROR: {var} environment variable not set. See the repo README.")
    path = Path(value)
    if not path.is_dir():
        sys.exit(f"ERROR: {var}={value} is not a directory.")
    return path


def find_trajectory(runs: Path, date: pd.Timestamp) -> Path:
    """Locate the cohort's trajectory file among the KRICO_NNNN folders."""
    name = f"{date:%Y_%m_%d}.nc"
    hits = sorted(runs.glob(f"KRICO_*/{name}"))
    if len(hits) != 1:
        sys.exit(f"ERROR: expected one {name} under {runs}/KRICO_*, found {len(hits)}.")
    return hits[0]


def read_trajectory(path: Path, stride: int):
    """Daily lon/lat as (particle, day) float32 arrays, NaN where invalid."""
    with Dataset(path) as nc:
        lon = np.ma.filled(nc["lon"][::stride, :].astype(np.float32), np.nan)
        lat = np.ma.filled(nc["lat"][::stride, :].astype(np.float32), np.nan)
    # Decode the release time as process_cohort.py does, through xarray.
    with xr.open_dataset(path, decode_timedelta=False) as ds:
        t0 = pd.Timestamp(ds["time"].isel(trajectory=0, obs=0).values)
    return lon, lat, t0.normalize()


def read_outcomes(path: Path, stride: int) -> xr.Dataset:
    ds = xr.open_dataset(path)
    return ds.isel(particle=slice(None, None, stride)).load()


def regular_axes(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """1-D lon/lat axes from the 2-D NEMO coordinates, checked to be regular."""
    lon = np.asarray(ds["nav_lon"].values)
    lat = np.asarray(ds["nav_lat"].values)
    if lon.ndim == 2:
        if not np.allclose(lon, lon[0:1, :], equal_nan=True):
            raise ValueError("nav_lon varies along y: grid is not regular.")
        lon = lon[0, :]
    if lat.ndim == 2:
        if not np.allclose(lat, lat[:, 0:1], equal_nan=True):
            raise ValueError("nav_lat varies along x: grid is not regular.")
        lat = lat[:, 0]
    return lon, lat


def block_mean(field: np.ndarray, k: int) -> np.ndarray:
    """k x k block average ignoring NaN; all-NaN blocks stay NaN (land)."""
    ny, nx = (field.shape[0] // k) * k, (field.shape[1] // k) * k
    blocks = field[:ny, :nx].reshape(ny // k, k, nx // k, k)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(blocks, axis=(1, 3))


def read_sea_ice(glorys: Path, dates: pd.DatetimeIndex):
    """
    Daily sea-ice concentration over the domain for each date, coarsened.

    Daily means are stamped at midday (12:00 up to 2018-12-18, 11:52:30
    afterwards), so midday is requested and the returned day is checked, as
    in krico_recruitment.sea_ice.spawning_sic.
    """
    out = None
    lon_c = lat_c = None
    for month_start in pd.date_range(dates[0].replace(day=1), dates[-1], freq="MS"):
        path = glorys / GLORYS_ICE_PATTERN.format(year=month_start.year,
                                                  month=month_start.month)
        with xr.open_dataset(path) as ds:
            lon, lat = regular_axes(ds)
            ii = np.flatnonzero((lon >= LON_MIN) & (lon <= LON_MAX))
            jj = np.flatnonzero((lat >= LAT_MIN) & (lat <= LAT_MAX))
            if lon_c is None:
                k = SIC_COARSEN
                lon_c = lon[ii][: (ii.size // k) * k].reshape(-1, k).mean(axis=1)
                lat_c = lat[jj][: (jj.size // k) * k].reshape(-1, k).mean(axis=1)
                out = np.full((dates.size, lat_c.size, lon_c.size), SIC_FILL, np.uint8)

            in_month = np.flatnonzero((dates.year == month_start.year)
                                      & (dates.month == month_start.month))
            for t in in_month:
                target = np.datetime64(dates[t].strftime("%Y-%m-%dT12:00:00"))
                da = ds["ileadfra"].sel(time_counter=target, method="nearest")
                got = pd.Timestamp(da["time_counter"].values).normalize()
                if got != dates[t]:
                    raise ValueError(f"asked sea ice for {dates[t].date()}, "
                                     f"got {got.date()} from {path.name}")
                field = np.squeeze(np.asarray(da.values, dtype=np.float64))
                field = field[np.ix_(jj, ii)]
                coarse = block_mean(field, SIC_COARSEN)
                pct = np.where(np.isnan(coarse), SIC_FILL,
                               np.clip(np.rint(coarse * 100.0), 0, 100))
                out[t] = pct.astype(np.uint8)
        print(f"  sea ice: {path.name}")
    return out, lon_c.astype(np.float32), lat_c.astype(np.float32)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--date", help="release date YYYY-MM-DD (default: chosen from F2)")
    parser.add_argument("--stride", type=int, default=1,
                        help="keep every Nth particle (default 1, all)")
    parser.add_argument("--out", type=Path, default=OUT_PATH)
    args = parser.parse_args()

    runs = require_dir("KRICO_RUNS")
    post = require_dir("KRICO_POST")
    glorys = require_dir("KRICO_GLORYS12")

    if args.date:
        date = pd.Timestamp(args.date).normalize()
        rationale = "release date given on the command line"
    else:
        date, rationale = choose_cohort()
    print(f"Cohort: {date.date()} ({rationale})")

    traj_path = find_trajectory(runs, date)
    rec_path = post / "recruitment" / "data" / f"{date:%Y_%m_%d}.nc"
    print(f"Trajectories: {traj_path}")
    print(f"Outcomes:     {rec_path}")

    lon, lat, t0 = read_trajectory(traj_path, args.stride)
    if t0 != date:
        sys.exit(f"ERROR: trajectory starts on {t0.date()}, expected {date.date()}.")
    rec = read_outcomes(rec_path, args.stride)
    n_particles, n_obs = lon.shape
    if rec.sizes["particle"] != n_particles:
        sys.exit("ERROR: trajectory and recruitment files differ in particle count.")
    print(f"  {n_particles} particles, {n_obs} daily positions")

    outcome = rec["outcome"].values.astype(np.int8)
    fate_day = rec["travel_time"].values.astype(np.int16)
    final_lon = rec["final_lon"].values.astype(np.float32)
    final_lat = rec["final_lat"].values.astype(np.float32)

    # Row alignment check: the trajectory position on the fate day must be the
    # archived fate position. Exact up to float32 round-trip.
    rows = np.arange(n_particles)
    traj_fate_lon = lon[rows, np.clip(fate_day, 0, n_obs - 1)]
    both = np.isfinite(traj_fate_lon) & np.isfinite(final_lon)
    mismatch = np.mean(np.abs(traj_fate_lon[both] - final_lon[both]) > 1e-3)
    print(f"  fate-position mismatch: {mismatch:.2e} of particles")
    if mismatch > 1e-4:
        sys.exit("ERROR: trajectory and recruitment rows do not line up.")

    # Positions after the fate day are not needed: the particle stops there.
    after = np.arange(n_obs)[None, :] > fate_day[:, None]
    lon[after] = np.nan
    lat[after] = np.nan

    dates = pd.date_range(date, periods=n_obs, freq="D")
    sic, sic_lon, sic_lat = read_sea_ice(glorys, dates)

    flag_meanings = rec["outcome"].attrs["flag_meanings"]
    flag_values = rec["outcome"].attrs["flag_values"]

    # Packing to int16 is done by xarray from this encoding (NaN -> fill).
    position_encoding = {
        "dtype": "int16", "scale_factor": POSITION_SCALE,
        "_FillValue": POSITION_FILL, "zlib": True, "complevel": 4,
        "chunksizes": (1, n_particles),
    }
    ds = xr.Dataset(
        data_vars={
            # Stored (day, particle) and chunked by day: plot.py reads one day
            # per frame.
            "lon": (("day", "particle"), lon.T,
                    {"units": "degrees_east", "long_name": "daily longitude, "
                     "fill after the fate day"}),
            "lat": (("day", "particle"), lat.T,
                    {"units": "degrees_north", "long_name": "daily latitude, "
                     "fill after the fate day"}),
            "outcome": ("particle", outcome,
                        {"long_name": "recruitment outcome",
                         "flag_values": flag_values,
                         "flag_meanings": flag_meanings}),
            "fate_day": ("particle", fate_day,
                         {"long_name": "day index at which the outcome is decided",
                          "units": "days since release"}),
            "final_lon": ("particle", final_lon, {"units": "degrees_east"}),
            "final_lat": ("particle", final_lat, {"units": "degrees_north"}),
            "sic": (("day", "sic_lat", "sic_lon"), sic,
                    {"long_name": "GLORYS12 sea-ice concentration, "
                     f"{SIC_COARSEN}x{SIC_COARSEN} block mean",
                     "units": "percent", "_FillValue": SIC_FILL}),
        },
        coords={
            "day": ("day", np.arange(n_obs, dtype=np.int16)),
            "sic_lon": ("sic_lon", sic_lon, {"units": "degrees_east"}),
            "sic_lat": ("sic_lat", sic_lat, {"units": "degrees_north"}),
        },
        attrs={
            "title": "MS1 movie aggregation: one cohort followed day by day",
            "release_date": str(date.date()),
            "spawning_offset_days": SPAWNING_OFFSET_DAYS,
            "cohort_choice": rationale,
            "particle_stride": args.stride,
            "trajectory_file": str(traj_path),
            "recruitment_file": str(rec_path),
        },
    )
    encoding = {
        "lon": position_encoding,
        "lat": position_encoding,
        "sic": {"zlib": True, "complevel": 4,
                "chunksizes": (1, sic_lat.size, sic_lon.size)},
        "outcome": {"zlib": True}, "fate_day": {"zlib": True},
        "final_lon": {"zlib": True}, "final_lat": {"zlib": True},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(args.out, encoding=encoding)
    size_mb = args.out.stat().st_size / 1e6
    print(f"Wrote {args.out} ({size_mb:.0f} MB)")

    print("Outcome fractions in this cohort:")
    for code, label in zip(np.atleast_1d(flag_values), flag_meanings.split()):
        print(f"  {label:24s} {100 * np.mean(outcome == code):6.2f}%")


if __name__ == "__main__":
    main()
