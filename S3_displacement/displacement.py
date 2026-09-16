"""
S3: how far particles travel between release and the moment their fate is set.

Supports the Discussion claim that advection displaces particles a substantial
distance between release and sea-ice advance. Two quantities are reported:

  displacement  great-circle distance from the release position to the fate
                position, i.e. net transport
  path length   along-track distance actually covered, read from the archived
                `trajectory_path_length` variable

Their ratio is a straightness index: 1 would be a straight line, lower values
mean meandering or recirculation.

Headline population is the three outcomes evaluated at sea-ice advance —
recruitment success, M5a and M5b — because that is what the Discussion
sentence is about. M1 is excluded: its fate position is the release position
by construction, so its displacement is identically zero and would drag any
pooled statistic down without meaning. The remaining outcomes are reported
separately for context, M6 in particular, since those are the export cases.

Reads only the archived recruitment outputs. No trajectories, no GLORYS12.

Usage:
    export KRICO_POST=/path/to/krico-post-production
    python displacement.py
    python displacement.py --limit 40
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


OUTCOME_NAMES = (
    "success", "censored", "killed_M1", "killed_M4",
    "killed_M5_no_FIV", "killed_M5_not_on_shelf",
    "killed_M6_no_advance", "exited_domain",
)
CODE = {n: i for i, n in enumerate(OUTCOME_NAMES)}

# Evaluated at the sea-ice advance event: the population the Discussion means.
AT_ADVANCE = (CODE["success"], CODE["killed_M5_no_FIV"], CODE["killed_M5_not_on_shelf"])

EARTH_RADIUS_KM = 6371.0
COHORT_RE = re.compile(r"^(\d{4})_(\d{2})_(\d{2})\.nc$")

# Percentiles reported. The median is the headline; the tails say how much
# spread sits behind it.
PCTL = (5, 25, 50, 75, 95)


def great_circle_km(lon1, lat1, lon2, lat2):
    """Great-circle distance in km. Haversine, vectorised."""
    lon1, lat1, lon2, lat2 = (np.radians(np.asarray(v, dtype="float64"))
                              for v in (lon1, lat1, lon2, lat2))
    dlon, dlat = lon2 - lon1, lat2 - lat1
    a = np.sin(dlat / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2) ** 2
    return 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def run(data_dir: Path, out_path: Path, limit: int | None = None) -> None:
    files = sorted(p for p in data_dir.iterdir() if COHORT_RE.match(p.name))
    if not files:
        raise FileNotFoundError(f"No cohort files found in {data_dir}")
    if limit is not None:
        files = files[:limit]
    print(f"Found {len(files)} cohort files in {data_dir}")

    # Accumulate per-outcome histograms rather than holding 2.1e9 values in
    # memory. 1 km bins to 20000 km is ample and keeps percentiles exact to
    # within one bin.
    edges = np.arange(0, 20001, 1.0)
    hist_disp = {c: np.zeros(len(edges) - 1, dtype=np.int64) for c in range(8)}
    hist_path = {c: np.zeros(len(edges) - 1, dtype=np.int64) for c in range(8)}
    counts = np.zeros(8, dtype=np.int64)
    sum_disp = np.zeros(8)
    sum_path = np.zeros(8)

    for n, path in enumerate(files, 1):
        with xr.open_dataset(path) as ds:
            outcome = ds["outcome"].values
            disp = great_circle_km(ds["release_lon"].values, ds["release_lat"].values,
                                   ds["final_lon"].values, ds["final_lat"].values)
            plen = ds["trajectory_path_length"].values.astype("float64") / 1000.0

        ok = np.isfinite(disp) & np.isfinite(plen)
        for c in range(8):
            m = ok & (outcome == c)
            if not m.any():
                continue
            counts[c] += m.sum()
            sum_disp[c] += disp[m].sum()
            sum_path[c] += plen[m].sum()
            hist_disp[c] += np.histogram(disp[m], bins=edges)[0]
            hist_path[c] += np.histogram(plen[m], bins=edges)[0]

        if n % 200 == 0:
            print(f"  processed {n} / {len(files)}")

    centres = 0.5 * (edges[:-1] + edges[1:])

    def pctls(h):
        tot = h.sum()
        if tot == 0:
            return {p: float("nan") for p in PCTL}
        cum = np.cumsum(h) / tot
        return {p: float(centres[np.searchsorted(cum, p / 100.0)]) for p in PCTL}

    rows = []
    for c in range(8):
        if counts[c] == 0:
            continue
        d, p = pctls(hist_disp[c]), pctls(hist_path[c])
        rows.append({
            "population": OUTCOME_NAMES[c], "n": int(counts[c]),
            "disp_mean_km": sum_disp[c] / counts[c],
            **{f"disp_p{q}_km": d[q] for q in PCTL},
            "path_mean_km": sum_path[c] / counts[c],
            **{f"path_p{q}_km": p[q] for q in PCTL},
            # M1 never moves, so its mean path length is zero and the ratio
            # is undefined rather than infinite.
            "straightness": ((sum_disp[c] / counts[c]) / (sum_path[c] / counts[c])
                             if sum_path[c] > 0 else float("nan")),
        })

    # Pooled population evaluated at advance: the Discussion's subject.
    hd = sum(hist_disp[c] for c in AT_ADVANCE)
    hp = sum(hist_path[c] for c in AT_ADVANCE)
    n_adv = int(counts[list(AT_ADVANCE)].sum())
    d, p = pctls(hd), pctls(hp)
    md = sum_disp[list(AT_ADVANCE)].sum() / n_adv
    mp = sum_path[list(AT_ADVANCE)].sum() / n_adv
    rows.insert(0, {
        "population": "at_advance (success + M5a + M5b)", "n": n_adv,
        "disp_mean_km": md, **{f"disp_p{q}_km": d[q] for q in PCTL},
        "path_mean_km": mp, **{f"path_p{q}_km": p[q] for q in PCTL},
        "straightness": md / mp,
    })

    df = pd.DataFrame(rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False, float_format="%.1f")
    print(f"\nWrote {out_path}\n")

    print("=" * 96)
    print("Distance from release to the moment fate is determined")
    print("=" * 96)
    print(f"{'population':34s} {'n':>14s} {'disp med':>9s} {'disp IQR':>15s} "
          f"{'path med':>9s} {'straight':>9s}")
    print("-" * 96)
    for r in rows:
        print(f"{r['population']:34s} {r['n']:14,d} {r['disp_p50_km']:8.0f}  "
              f"{r['disp_p25_km']:6.0f}-{r['disp_p75_km']:<7.0f} "
              f"{r['path_p50_km']:8.0f}  {r['straightness']:8.2f}")
    print("-" * 96)
    print("All distances in km. 'straight' is mean displacement / mean path length.")
    print("M1 displacement is zero by construction and is excluded from the")
    print("at_advance population.\n")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Displacement and path length from release to fate.")
    ap.add_argument("--post", type=str, default=os.environ.get("KRICO_POST"))
    ap.add_argument("--out", type=str, default="data/displacement.csv")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    if not args.post:
        sys.exit("KRICO_POST is not set and --post was not given.")
    data_dir = Path(args.post) / "recruitment" / "data"
    if not data_dir.is_dir():
        sys.exit(f"Not a directory: {data_dir}")
    run(data_dir, Path(args.out), args.limit)


if __name__ == "__main__":
    main()
