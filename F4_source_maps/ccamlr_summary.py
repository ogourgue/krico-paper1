"""
F4 CCAMLR summary: per-fate, per-CCAMLR-subarea breakdown.

Reads data/aggregated.nc (the 1° x 1° per-fate density grid produced by
aggregate.py) and reports, for each outcome, the fraction of that
outcome's particles released in each CCAMLR subarea. Subarea 48.6 is
split at 60°S into northern and southern halves.

Produces two outputs:
  - stdout: a formatted table for human reading.
  - ccamlr_summary.csv (sibling of this script): wide-format CSV for SI
    Table S1, with rows = outcomes and columns = subareas. The CSV is
    committed to the repo (it is a small, stable text artifact, unlike
    the gitignored aggregated.nc grid).

Usage:
  python ccamlr_summary.py

Assumes:
  - $KRICO_ROOT is set
  - data/aggregated.nc exists in the same folder
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr
from shapely.geometry import Point


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# CCAMLR subareas of interest, with full geographic names.
# 48.6 is split at 60°S into northern and southern halves.
SUBAREA_CODES = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6", "88.3"]
SUBAREA_NAMES = {
    "48.1": "Antarctic Peninsula",
    "48.2": "South Orkney Islands",
    "48.3": "South Georgia",
    "48.4": "South Sandwich Islands",
    "48.5": "Weddell Sea",
    "48.6N": "Bouvet Island (north of 60°S)",
    "48.6S": "Bouvet Island (south of 60°S)",
    "88.3": "Amundsen Sea",
}
# Order for the per-fate report (same as F1 legend; 48.6 split into N/S).
REPORT_ORDER = ["48.1", "48.2", "48.3", "48.4", "48.5", "48.6N", "48.6S", "88.3"]

# Latitude split for 48.6.
SPLIT_LAT_486 = -60.0

# Outcome display order (top of F3 stack to bottom).
DISPLAY_ORDER = [
    "success",
    "killed_M6_no_advance",
    "killed_M5_not_on_shelf",
    "killed_M5_no_FIV",
    "killed_M4",
    "killed_M1",
    "exited_domain",
]
DISPLAY_LABELS = {
    "success":                "Recruitment success",
    "killed_M6_no_advance":   "No winter ice (M6)",
    "killed_M5_not_on_shelf": "Off-shelf (M5b)",
    "killed_M5_no_FIV":       "Under-developed (M5a)",
    "killed_M4":              "Calyptope starvation (M4)",
    "killed_M1":              "Ice at spawning (M1)",
    "exited_domain":          "Domain exit",
}

# Fates that should be folded together before reporting (consistent with F4 plot).
FOLD_INTO = {
    "censored": "killed_M6_no_advance",
}


# ---------------------------------------------------------------------------
# CCAMLR mapping
# ---------------------------------------------------------------------------

def load_ccamlr_geometries() -> gpd.GeoDataFrame:
    """Load CCAMLR statistical area polygons from the bundled shapefile."""
    repo_root = Path(__file__).resolve().parent.parent
    shp_path = repo_root / "ccamlr-data" / "CCAMLR_ASD_EPSG4326.shp"
    if not shp_path.exists():
        raise FileNotFoundError(f"CCAMLR shapefile not found: {shp_path}")
    ccamlr = gpd.read_file(shp_path)
    return ccamlr[ccamlr["GAR_Long_L"].isin(SUBAREA_CODES)]


def build_cell_subarea_map(lon_centers: np.ndarray,
                           lat_centers: np.ndarray) -> np.ndarray:
    """
    Return a (n_lat, n_lon) array of subarea labels for each grid cell.

    Cells outside any project subarea are labeled "outside". Cells within
    48.6 are split into "48.6N" / "48.6S" based on cell-center latitude
    (split at SPLIT_LAT_486 = -60°).
    """
    ccamlr = load_ccamlr_geometries()

    n_lon = len(lon_centers)
    n_lat = len(lat_centers)
    labels = np.full((n_lat, n_lon), "outside", dtype=object)

    # For each subarea, build a single MultiPolygon and test cell centers.
    # At 1° resolution this is fast (~16k cells × 7 subareas).
    from shapely.ops import unary_union
    for code in SUBAREA_CODES:
        rows = ccamlr[ccamlr["GAR_Long_L"] == code]
        if rows.empty:
            continue
        geom = unary_union(rows.geometry.values)
        for i, lat in enumerate(lat_centers):
            for j, lon in enumerate(lon_centers):
                if labels[i, j] != "outside":
                    continue
                if geom.contains(Point(lon, lat)):
                    if code == "48.6":
                        labels[i, j] = "48.6N" if lat >= SPLIT_LAT_486 else "48.6S"
                    else:
                        labels[i, j] = code

    return labels


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def fold_outcomes(counts: np.ndarray, outcome_names: list[str]) -> np.ndarray:
    """Fold outcomes (e.g. censored -> killed_M6_no_advance) at the count level."""
    out = counts.copy()
    for src, dst in FOLD_INTO.items():
        if src in outcome_names and dst in outcome_names:
            si = outcome_names.index(src)
            di = outcome_names.index(dst)
            out[di] = out[di] + out[si]
            out[si] = 0
    return out


def per_fate_per_subarea(
    counts: np.ndarray,
    cell_labels: np.ndarray,
    outcome_names: list[str],
) -> dict[str, dict[str, int]]:
    """
    Build a nested dict {outcome: {subarea: count}} from the per-cell counts.
    Each cell is assigned to exactly one subarea (or "outside").
    """
    n_outcomes, n_lat, n_lon = counts.shape
    result: dict[str, dict[str, int]] = {
        oname: {label: 0 for label in REPORT_ORDER + ["outside"]}
        for oname in outcome_names
    }

    # Iterate per outcome (small N).
    for oi, oname in enumerate(outcome_names):
        layer = counts[oi]
        for label in REPORT_ORDER + ["outside"]:
            mask = (cell_labels == label)
            result[oname][label] = int(layer[mask].sum())

    return result


def to_percent_table(per_fate: dict[str, dict[str, int]]) -> dict[str, dict[str, float]]:
    """Convert raw counts to within-outcome percentages."""
    pct: dict[str, dict[str, float]] = {}
    for oname, row in per_fate.items():
        total = sum(row.values())
        pct[oname] = (
            {k: 100.0 * v / total for k, v in row.items()}
            if total > 0
            else {k: 0.0 for k in row}
        )
    return pct


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def print_summary(pct_table: dict[str, dict[str, float]],
                  outcome_names: list[str],
                  figure_label: str) -> None:
    """Print the per-fate, per-subarea breakdown to stdout."""
    print()
    print("=" * 84)
    print(f"{figure_label}: per-fate breakdown by CCAMLR subarea (within-outcome %)")
    print("=" * 84)

    cols = REPORT_ORDER + ["outside"]
    header = f"{'Outcome':28s} " + " ".join(f"{c:>7s}" for c in cols)
    print(header)
    print("-" * len(header))

    for oname in DISPLAY_ORDER:
        if oname not in outcome_names or oname not in pct_table:
            continue
        row = pct_table[oname]
        cells = " ".join(f"{row[c]:7.2f}" for c in cols)
        print(f"{DISPLAY_LABELS.get(oname, oname):28s} {cells}")

    print()
    print("Subarea names:")
    for code in cols:
        if code == "outside":
            print(f"  {code:7s} : outside CCAMLR project subareas")
        else:
            print(f"  {code:7s} : {SUBAREA_NAMES.get(code, '?')}")


def write_csv(pct_table: dict[str, dict[str, float]],
              outcome_names: list[str],
              csv_path: Path,
              figure_label: str) -> None:
    """
    Write the per-fate × per-subarea percentage table as a wide-format CSV.

    Rows are outcomes (in DISPLAY_ORDER), columns are subareas + "outside".
    Includes a header comment line identifying the source figure and the
    convention (within-outcome percentages, rows sum to 100).
    """
    cols = REPORT_ORDER + ["outside"]
    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        # Comment header, prefixed with '#'. Most CSV readers will treat
        # the first non-comment row as the header. Pandas can ingest this
        # via pd.read_csv(..., comment='#').
        writer.writerow([
            f"# {figure_label}: per-fate breakdown by CCAMLR subarea "
            f"(within-outcome %, rows sum to 100)"
        ])
        writer.writerow(["Outcome"] + cols)
        for oname in DISPLAY_ORDER:
            if oname not in outcome_names or oname not in pct_table:
                continue
            row = pct_table[oname]
            writer.writerow(
                [DISPLAY_LABELS.get(oname, oname)]
                + [f"{row[c]:.2f}" for c in cols]
            )
    print(f"Wrote {csv_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    here = Path(__file__).resolve().parent
    aggregated_path = here / "data" / "aggregated.nc"
    if not aggregated_path.exists():
        print(f"ERROR: aggregated file not found: {aggregated_path}", file=sys.stderr)
        print("Run aggregate.py first.", file=sys.stderr)
        sys.exit(1)

    ds = xr.open_dataset(aggregated_path)
    counts = ds["counts"].values                       # (outcome, lat, lon)
    lon_centers = ds["lon"].values
    lat_centers = ds["lat"].values
    outcome_names = list(ds["outcome"].values.astype(str))

    print("Building cell -> subarea mapping...")
    cell_labels = build_cell_subarea_map(lon_centers, lat_centers)
    print(f"  {(cell_labels != 'outside').sum()} cells inside project subareas, "
          f"{(cell_labels == 'outside').sum()} outside.")

    counts = fold_outcomes(counts, outcome_names)
    per_fate = per_fate_per_subarea(counts, cell_labels, outcome_names)
    pct_table = to_percent_table(per_fate)

    figure_label = "F4 (release positions)"
    print_summary(pct_table, outcome_names, figure_label=figure_label)
    write_csv(pct_table, outcome_names, here / "ccamlr_summary.csv",
              figure_label=figure_label)


if __name__ == "__main__":
    main()