"""
FS1 aggregation: time-mean GLORYS12 currents over the release depth band.

Reduces the reformatted GLORYS12 velocity files to depth-averaged, time-mean
horizontal velocity fields over the 50-200 m band, for the annual mean and
for each of the four seasons. Runs on the HPC where Pre/GLORYS12/ lives; the
output is small enough to copy back and plot locally.

What is averaged, and why:

  - The mean of the velocity *vector*, not of the speed. |(u_bar, v_bar)|
    measures persistent net transport, which is what displaces particles
    over 200 days; the mean of |(u, v)| would instead measure how energetic
    the flow is and would light up eddy fields that produce little net
    displacement. Regions that look weak here are regions of little net
    transport, not necessarily quiet regions.

  - Velocity, not depth-integrated transport. Transport scales with the
    available water column, so it understates fast flow over the shelf,
    which is precisely where the retention argument lives. Particles are
    advected by velocity.

  - 50-200 m is the release depth band. Particles are advected in three
    dimensions and sample well outside it over 200 days, so this is the
    band every trajectory starts in, not the depth larvae occupy. Say so in
    the caption.

  - Depth averaging is weighted by layer thickness, since the GLORYS
    vertical grid is stretched (about 11 m per level at the top of the band,
    35 m at the bottom). Thicknesses come from the midpoints between level
    centres, clipped to the band; the files carry no depth bounds and no
    e3t, so NEMO partial bottom cells are not represented. Where the seabed
    is shallower than 200 m only the levels present are averaged, so the
    effective averaging depth varies in space; it is written to the output
    as `depth_averaged` and should be stated in the caption.

  - Time averaging accumulates daily fields and divides by the count, so
    months are weighted by their length automatically.

Inputs:  $KRICO_ROOT/Pre/GLORYS12/glorys12_u_YYYY_MM.nc
         $KRICO_ROOT/Pre/GLORYS12/glorys12_v_YYYY_MM.nc
Output:  data/aggregated.nc

Usage:
    python aggregate.py                 # all files, resuming if interrupted
    python aggregate.py --no-resume     # ignore any checkpoint and restart
    python aggregate.py --limit 3       # first 3 month-pairs, for testing
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path

import numpy as np
import xarray as xr


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Release depth band, in metres.
DEPTH_MIN = 50.0
DEPTH_MAX = 200.0

# Periods accumulated in a single pass over the files. "annual" takes every
# day; the others select on calendar month.
PERIODS = {
    "annual": None,
    "DJF": (12, 1, 2),
    "MAM": (3, 4, 5),
    "JJA": (6, 7, 8),
    "SON": (9, 10, 11),
}

# Checkpoint cadence, in month-pairs.
CHECKPOINT_EVERY = 12

FILENAME_RE = re.compile(r"^glorys12_u_(\d{4})_(\d{2})\.nc$")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def find_month_pairs(forcing_dir: Path) -> list[tuple[int, int, Path, Path]]:
    """All (year, month, u_path, v_path) with both components present."""
    pairs = []
    for path in sorted(forcing_dir.iterdir()):
        m = FILENAME_RE.match(path.name)
        if not m:
            continue
        year, month = int(m.group(1)), int(m.group(2))
        v_path = path.with_name(f"glorys12_v_{year:04d}_{month:02d}.nc")
        if not v_path.exists():
            print(f"  WARNING: no v file for {path.name}, skipping",
                  file=sys.stderr)
            continue
        pairs.append((year, month, path, v_path))
    return pairs


def velocity_var(ds: xr.Dataset) -> str:
    """Name of the single 4-D velocity variable in a reformatted file."""
    candidates = [name for name, da in ds.data_vars.items() if da.ndim == 4]
    if len(candidates) != 1:
        raise ValueError(
            f"expected exactly one 4-D variable, found {candidates}")
    return candidates[0]


def layer_weights(depth: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Thickness weights for the levels overlapping [DEPTH_MIN, DEPTH_MAX].

    The files carry level centres only, so cell edges are taken as the
    midpoints between centres, with the top edge at the surface and the
    bottom edge extrapolated symmetrically. Each cell is then clipped to the
    band, so the partially overlapping cells at the top and bottom
    contribute only their overlapping thickness.

    Returns the level indices and their weights, in metres.
    """
    edges = np.empty(depth.size + 1)
    edges[0] = 0.0
    edges[1:-1] = 0.5 * (depth[:-1] + depth[1:])
    edges[-1] = depth[-1] + (depth[-1] - edges[-2])

    top = np.maximum(edges[:-1], DEPTH_MIN)
    bot = np.minimum(edges[1:], DEPTH_MAX)
    w = np.maximum(bot - top, 0.0)

    idx = np.flatnonzero(w > 0.0)
    return idx, w[idx]


def depth_average(block: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """
    Thickness-weighted depth mean of a (level, y, x) block, ignoring levels
    that are missing because the seabed is shallower than the band.

    Returns (y, x); cells with no valid level are NaN.
    """
    valid = np.isfinite(block)
    w = weights[:, None, None] * valid
    wsum = w.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(wsum > 0,
                       np.nansum(np.where(valid, block, 0.0) * w, axis=0)
                       / np.maximum(wsum, 1e-12),
                       np.nan)
    return out


def effective_thickness(block: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Total weighted thickness present at each cell, in metres."""
    valid = np.isfinite(block)
    return (weights[:, None, None] * valid).sum(axis=0)


# ---------------------------------------------------------------------------
# Accumulation
# ---------------------------------------------------------------------------

class Accumulator:
    """Running sums of the depth-averaged daily fields, per period."""

    def __init__(self, shape: tuple[int, int]):
        self.shape = shape
        self.u = {k: np.zeros(shape) for k in PERIODS}
        self.v = {k: np.zeros(shape) for k in PERIODS}
        self.n = {k: 0 for k in PERIODS}
        self.done: list[str] = []
        self.depth_averaged: np.ndarray | None = None

    def add_day(self, month: int, u2d: np.ndarray, v2d: np.ndarray) -> None:
        for name, months in PERIODS.items():
            if months is None or month in months:
                # NaN over land; treat as zero in the sum and rely on the
                # static mask, which is identical for every timestep.
                self.u[name] += np.nan_to_num(u2d)
                self.v[name] += np.nan_to_num(v2d)
                self.n[name] += 1

    def to_npz(self, path: Path) -> None:
        np.savez_compressed(
            path,
            **{f"u_{k}": v for k, v in self.u.items()},
            **{f"v_{k}": v for k, v in self.v.items()},
            **{f"n_{k}": np.array(v) for k, v in self.n.items()},
            done=np.array(self.done),
            depth_averaged=(self.depth_averaged
                            if self.depth_averaged is not None
                            else np.zeros(self.shape)),
        )

    @classmethod
    def from_npz(cls, path: Path) -> "Accumulator":
        z = np.load(path, allow_pickle=False)
        shape = z["u_annual"].shape
        acc = cls(shape)
        for k in PERIODS:
            acc.u[k] = z[f"u_{k}"]
            acc.v[k] = z[f"v_{k}"]
            acc.n[k] = int(z[f"n_{k}"])
        acc.done = [str(s) for s in z["done"]]
        acc.depth_averaged = z["depth_averaged"]
        return acc


# ---------------------------------------------------------------------------
# Main pass
# ---------------------------------------------------------------------------

def aggregate(forcing_dir: Path, out_dir: Path, resume: bool,
              limit: int | None) -> xr.Dataset:
    pairs = find_month_pairs(forcing_dir)
    if not pairs:
        raise FileNotFoundError(f"No glorys12_u_*.nc files in {forcing_dir}")
    if limit is not None:
        pairs = pairs[:limit]
    print(f"Found {len(pairs)} month-pairs in {forcing_dir}")

    checkpoint = out_dir / "checkpoint.npz"
    acc: Accumulator | None = None
    if resume and checkpoint.exists():
        acc = Accumulator.from_npz(checkpoint)
        print(f"Resuming from checkpoint: {len(acc.done)} month-pairs done")

    lon = lat = None
    levels = weights = None

    for i, (year, month, u_path, v_path) in enumerate(pairs, start=1):
        tag = f"{year:04d}_{month:02d}"
        if acc is not None and tag in acc.done:
            continue

        with xr.open_dataset(u_path) as du, xr.open_dataset(v_path) as dv:
            if levels is None:
                depth = du["deptht"].values
                levels, weights = layer_weights(depth)
                print(f"Depth band {DEPTH_MIN:.0f}-{DEPTH_MAX:.0f} m: "
                      f"{levels.size} levels")
                for k, w in zip(levels, weights):
                    print(f"  level {k:2d}  centre {depth[k]:7.2f} m  "
                          f"weight {w:6.2f} m")
                lon = du["nav_lon"].values
                lat = du["nav_lat"].values

            uname, vname = velocity_var(du), velocity_var(dv)
            u_block = du[uname].isel(deptht=levels)
            v_block = dv[vname].isel(deptht=levels)

            if acc is None:
                acc = Accumulator(u_block.shape[-2:])

            n_days = u_block.sizes["time_counter"]
            for t in range(n_days):
                ub = u_block.isel(time_counter=t).values
                vb = v_block.isel(time_counter=t).values
                if acc.depth_averaged is None or not acc.depth_averaged.any():
                    acc.depth_averaged = effective_thickness(ub, weights)
                acc.add_day(month,
                            depth_average(ub, weights),
                            depth_average(vb, weights))

        acc.done.append(tag)
        print(f"  [{i:3d}/{len(pairs)}] {tag}  {n_days} days  "
              f"(annual n = {acc.n['annual']})")

        if i % CHECKPOINT_EVERY == 0:
            acc.to_npz(checkpoint)
            print(f"    checkpoint written ({len(acc.done)} done)")

    assert acc is not None

    # ------------------------------------------------------------------
    # Means
    # ------------------------------------------------------------------
    data_vars = {}
    for name in PERIODS:
        n = max(acc.n[name], 1)
        # Restore the land mask, which nan_to_num removed from the sums.
        mask = np.isfinite(acc.depth_averaged) & (acc.depth_averaged > 0)
        u_mean = np.where(mask, acc.u[name] / n, np.nan).astype(np.float32)
        v_mean = np.where(mask, acc.v[name] / n, np.nan).astype(np.float32)
        data_vars[f"u_{name}"] = (
            ("y", "x"), u_mean,
            {"long_name": f"{name} mean zonal velocity, "
                          f"{DEPTH_MIN:.0f}-{DEPTH_MAX:.0f} m",
             "units": "m s-1", "n_days": acc.n[name]})
        data_vars[f"v_{name}"] = (
            ("y", "x"), v_mean,
            {"long_name": f"{name} mean meridional velocity, "
                          f"{DEPTH_MIN:.0f}-{DEPTH_MAX:.0f} m",
             "units": "m s-1", "n_days": acc.n[name]})

    data_vars["depth_averaged"] = (
        ("y", "x"), acc.depth_averaged.astype(np.float32),
        {"long_name": "total layer thickness averaged over, within the band",
         "units": "m"})

    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            "nav_lon": (("y", "x"), lon.astype(np.float32),
                        {"units": "degrees_east"}),
            "nav_lat": (("y", "x"), lat.astype(np.float32),
                        {"units": "degrees_north"}),
        },
        attrs={
            "title": "FS1 mean GLORYS12 circulation over the release depth band",
            "description": (
                "Thickness-weighted depth mean over "
                f"{DEPTH_MIN:.0f}-{DEPTH_MAX:.0f} m of the daily GLORYS12v1 "
                "velocity fields, then time-averaged. Mean of the velocity "
                "vector, not of the speed. Where the seabed is shallower "
                "than the band, only the levels present are averaged; see "
                "depth_averaged."),
            "depth_min_m": DEPTH_MIN,
            "depth_max_m": DEPTH_MAX,
            "n_month_pairs": len(acc.done),
            "first_month": acc.done[0] if acc.done else "",
            "last_month": acc.done[-1] if acc.done else "",
        },
    )

    print()
    print("Days accumulated per period:")
    for name in PERIODS:
        print(f"  {name:7s} {acc.n[name]:6d}")

    return ds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-resume", action="store_true",
                        help="ignore any existing checkpoint")
    parser.add_argument("--limit", type=int, default=None,
                        help="process only the first N month-pairs")
    args = parser.parse_args()

    krico_root = os.environ.get("KRICO_ROOT")
    if not krico_root:
        print("ERROR: KRICO_ROOT environment variable not set.",
              file=sys.stderr)
        sys.exit(1)

    forcing_dir = Path(krico_root) / "Pre" / "GLORYS12"
    if not forcing_dir.is_dir():
        print(f"ERROR: forcing dir not found: {forcing_dir}", file=sys.stderr)
        sys.exit(1)

    out_dir = Path(__file__).resolve().parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    ds = aggregate(forcing_dir, out_dir, resume=not args.no_resume,
                   limit=args.limit)

    out_path = out_dir / "aggregated.nc"
    print(f"Writing {out_path}")
    encoding = {v: {"zlib": True, "complevel": 5} for v in ds.data_vars}
    ds.to_netcdf(out_path, encoding=encoding)
    print("Done. The checkpoint in data/ can be deleted once this is copied "
          "back.")


if __name__ == "__main__":
    main()
