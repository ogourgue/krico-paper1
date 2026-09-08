"""
S2, second part: how much of M6 lies beyond the reach of winter sea ice.

The truncation sweep (sweep.py) brackets the seasonal decline in M6 from
below by holding every cohort to a common detection window. That bound is
deliberately pessimistic: it demotes real advance events that late-release
larvae did experience, wherever they occurred, and so removes successes from
regions the ice never reaches as readily as from regions it does.

This test asks the complementary question. A particle whose fate position
lies north of the winter sea-ice edge could not have met ice however long
tracking continued. Such particles are M6 for a physical reason rather than
a bookkeeping one, and the truncation sweep cannot distinguish them. The
larger their share, the less of M6 is window-sensitive, and the closer the
truth sits to the reported classification rather than to the truncated bound.

Method
------
The sea-ice edge is the 15% contour of the multi-year mean September
concentration field (NOAA/NSIDC G02202 v6), the same field and threshold
Figure 1 uses, read from F1_domain_map/data/aggregated.nc. September is the
climatological maximum, so this is the furthest north ice reaches in a
typical year.

The contour is reduced to a northern ice-edge latitude per 1-degree
longitude bin: the northernmost cell at that longitude reaching 15%. A
particle is counted as beyond the reach of ice if its fate latitude lies
north of that line. Three edges are reported: the 1979-2015 mean, the
2016-2025 mean, and their northernmost envelope, which is the most
conservative of the three because it counts fewest particles as unreachable.

Two approximations, both stated rather than hidden
--------------------------------------------------
1. The fate position is where the particle sat at the end of tracking, not
   where it would have been in September. Particles continue to drift. The
   Antarctic Circumpolar Current is largely zonal over this domain, so
   latitude is approximately conserved and the fate latitude is a reasonable
   proxy, but a particle drifting south after the cutoff would be misscored
   as unreachable.
2. The edge is a multi-year mean. Individual years vary about it, so a
   particle just north of the mean edge may meet ice in some years. Using
   the northernmost envelope of the two epochs mitigates this in the
   conservative direction.

Neither approximation is load-bearing for the qualitative result, which is
about the share of M6 sitting far north rather than marginally north.

Inputs:  F1_domain_map/data/aggregated.nc  (committed; no download needed)
         $KRICO_POST/recruitment/data/YYYY_MM_DD.nc
Output:  data/m6_ice_reachability.csv

Usage:
    export KRICO_POST=/path/to/krico-post-production
    python ice_reachability.py
    python ice_reachability.py --limit 40
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Same threshold Figure 1 contours (F1_domain_map/plot.py, SIC_THRESHOLD).
SIC_THRESHOLD = 0.15

# Month of the climatological sea-ice maximum, as stored in the F1 aggregation.
SIC_MONTH = "September"

# Outcome flags folded together as M6 (consistent with F3, F4 and F5).
M6_CODES = (1, 6)  # censored, killed_M6_no_advance

# Longitude bin width for the ice-edge lookup, in degrees.
LON_BIN_DEG = 1.0

F1_AGGREGATE = Path("../F1_domain_map/data/aggregated.nc")

COHORT_FILENAME_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")


# ---------------------------------------------------------------------------
# Sea-ice edge
# ---------------------------------------------------------------------------

def build_ice_edges(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """
    Northern sea-ice edge latitude per longitude bin, from the F1 aggregation.

    Returns
    -------
    bin_edges : ndarray
        Longitude bin boundaries, length n_bins + 1.
    edges : dict of str -> ndarray
        Edge latitude per bin for each period, plus "envelope", the
        northernmost of the periods. NaN where no cell at that longitude
        reaches the threshold.
    """
    with xr.open_dataset(path) as ds:
        lon = ds["sic_lon"].values.ravel()
        lat = ds["sic_lat"].values.ravel()
        periods = [str(p) for p in ds["period"].values]
        fields = {p: ds["sic_climatology"].sel(month=SIC_MONTH, period=p).values.ravel()
                  for p in periods}

    bin_edges = np.arange(-180.0, 180.0 + LON_BIN_DEG, LON_BIN_DEG)
    n_bins = len(bin_edges) - 1

    edges: dict[str, np.ndarray] = {}
    for period, sic in fields.items():
        iced = np.isfinite(sic) & (sic >= SIC_THRESHOLD)
        # Initialise to -inf, not NaN: np.maximum propagates NaN.
        e = np.full(n_bins, -np.inf)
        idx = np.clip(np.digitize(lon[iced], bin_edges) - 1, 0, n_bins - 1)
        np.maximum.at(e, idx, lat[iced])
        e[np.isinf(e)] = np.nan
        edges[period] = e

    stacked = np.vstack([edges[p] for p in periods])
    with np.errstate(invalid="ignore"):
        edges["envelope"] = np.nanmax(stacked, axis=0)

    return bin_edges, edges


def north_of_edge(lon: np.ndarray, lat: np.ndarray,
                  bin_edges: np.ndarray, edge: np.ndarray) -> np.ndarray:
    """
    True where a position lies north of the ice edge at its longitude.

    Bins with no ice at any latitude (NaN edge) count as north, since no ice
    was ever present there.
    """
    idx = np.clip(np.digitize(lon, bin_edges) - 1, 0, len(edge) - 1)
    edge_lat = edge[idx]
    return np.isnan(edge_lat) | (lat > edge_lat)


# ---------------------------------------------------------------------------
# Cohort walk
# ---------------------------------------------------------------------------

def spawning_year_for_date(year: int, month: int) -> int:
    return year + 1 if month in (11, 12) else year


def find_cohort_files(data_dir: Path) -> list[tuple[dt.date, Path]]:
    out = []
    for path in sorted(data_dir.iterdir()):
        m = COHORT_FILENAME_RE.match(path.name)
        if m:
            y, mo, d = (int(g) for g in m.groups())
            out.append((dt.date(y, mo, d), path))
    return out


def run(data_dir: Path, f1_path: Path, out_path: Path,
        limit: int | None = None) -> None:
    bin_edges, edges = build_ice_edges(f1_path)
    edge_names = [k for k in edges if k != "envelope"] + ["envelope"]

    print(f"Sea-ice edge from {f1_path} ({SIC_MONTH}, {SIC_THRESHOLD:.0%} contour)")
    for name in edge_names:
        e = edges[name]
        sel = (bin_edges[:-1] >= -115) & (bin_edges[:-1] <= 40)
        print(f"  {name:12s} over domain: {np.nanmin(e[sel]):.1f} to "
              f"{np.nanmax(e[sel]):.1f} degrees latitude")
    print()

    files = find_cohort_files(data_dir)
    if not files:
        raise FileNotFoundError(f"No cohort files found in {data_dir}")
    if limit is not None:
        files = files[:limit]
    print(f"Found {len(files)} cohort files in {data_dir}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    header = (["release_date", "spawning_year", "release_month",
               "n_particles", "n_M6"]
              + [f"n_M6_north_{name}" for name in edge_names])

    n_processed = 0
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)

        for release, path in files:
            with xr.open_dataset(path) as ds:
                outcome = ds["outcome"].values
                flon = ds["final_lon"].values
                flat = ds["final_lat"].values

            is_m6 = np.isin(outcome, M6_CODES)
            n_m6 = int(is_m6.sum())
            lon_m6 = flon[is_m6].astype("float64")
            lat_m6 = flat[is_m6].astype("float64")

            counts = []
            for name in edge_names:
                if n_m6 == 0:
                    counts.append(0)
                    continue
                counts.append(int(north_of_edge(lon_m6, lat_m6,
                                                bin_edges, edges[name]).sum()))

            writer.writerow(
                [release.isoformat(),
                 spawning_year_for_date(release.year, release.month),
                 release.month, int(outcome.size), n_m6] + counts
            )

            n_processed += 1
            if n_processed % 200 == 0:
                print(f"  processed {n_processed} / {len(files)}")

    print(f"Done: {n_processed} cohorts -> {out_path}")
    summarise(out_path, edge_names)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

MONTH_NAME = {11: "November", 12: "December", 1: "January",
              2: "February", 3: "March"}


def summarise(csv_path: Path, edge_names: list[str]) -> None:
    df = pd.read_csv(csv_path)

    print()
    print("=" * 78)
    print("Share of M6 lying north of the September sea-ice edge")
    print("=" * 78)
    print("A particle north of the edge could not have met ice however long")
    print("tracking continued: its M6 classification is physical, not an")
    print("artefact of the 200-day window.")
    print()

    cols = [f"n_M6_north_{n}" for n in edge_names]
    print(f"{'release month':>15s} {'M6 particles':>14s} "
          + " ".join(f"{n:>14s}" for n in edge_names))
    print("-" * 78)

    for month in (11, 12, 1, 2, 3):
        sub = df[df["release_month"] == month]
        if sub.empty:
            continue
        n_m6 = sub["n_M6"].sum()
        pcts = [100.0 * sub[c].sum() / n_m6 if n_m6 else float("nan") for c in cols]
        print(f"{MONTH_NAME[month]:>15s} {n_m6:14,d} "
              + " ".join(f"{p:13.1f}%" for p in pcts))

    n_m6 = df["n_M6"].sum()
    pcts = [100.0 * df[c].sum() / n_m6 for c in cols]
    print("-" * 78)
    print(f"{'all releases':>15s} {n_m6:14,d} "
          + " ".join(f"{p:13.1f}%" for p in pcts))

    print()
    print("The 'envelope' column is the most conservative: it uses the")
    print("northernmost ice edge of the two epochs, so it counts the fewest")
    print("particles as beyond reach.")
    print()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Share of M6 beyond the reach of winter sea ice."
    )
    parser.add_argument("--post", type=str, default=os.environ.get("KRICO_POST"))
    parser.add_argument("--f1", type=str, default=str(F1_AGGREGATE),
                        help="path to F1_domain_map/data/aggregated.nc")
    parser.add_argument("--out", type=str, default="data/m6_ice_reachability.csv")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    if not args.post:
        sys.exit("KRICO_POST is not set and --post was not given.")

    data_dir = Path(args.post) / "recruitment" / "data"
    if not data_dir.is_dir():
        sys.exit(f"Not a directory: {data_dir}")

    f1_path = Path(args.f1)
    if not f1_path.is_file():
        sys.exit(f"Not found: {f1_path}. Run from inside S2_m6_window_sensitivity/.")

    run(data_dir, f1_path, Path(args.out), args.limit)


if __name__ == "__main__":
    main()
