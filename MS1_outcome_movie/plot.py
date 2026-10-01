"""
MS1 plot: render the frames of Movie S1.

One frame per tracking day, plus FADE_DAYS hold frames at the end. Each frame
shows, on the F1 base map:

  - daily GLORYS12 sea-ice concentration above the 15% advance threshold,
  - the 2000 m isobath, the bathymetric limit of recruitment habitat (M5b),
  - every particle still drifting, in gray,
  - particles whose outcome is decided, frozen at their fate position:
    green for recruitment success, which stays; red for every mortality
    outcome, which fades out over FADE_DAYS days so the map stays readable.

Outcome timing follows the classification exactly, through the archived fate
day. M1 particles are decided at spawning, before release, so they are red on
the first frame and fade from there; they never drift. M6 and censored
particles are decided at the end of tracking, so they turn red together on the
last day, and the hold frames let them fade while success stays on screen.
Domain exits are a modelling limitation rather than an outcome: those
particles drift in gray until deleted, then disappear.

The base map reuses the F1 helpers (land, coastline, CCAMLR outlines, 48.6
split, domain boundary) by importing ../F1_domain_map/plot.py, so the two
cannot drift apart. The bathymetry zones, sea-ice climatology and legends of
F1 are left out; region names replace the subarea codes.

Reads data/aggregated.nc produced by aggregate.py, and the bathymetry in
../F1_domain_map/data/aggregated.nc for the isobath.

Usage
-----
    python plot.py                    # all frames into frames/
    python plot.py --frames 0:10      # a range, python slice semantics
    python plot.py --stride 10        # every 10th particle, for quick tests

Frame ranges are independent, so long renders can be split across processes.
Assemble the movie with make_movie.sh.
"""

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent

# F1 helpers, imported by path: the module is named plot.py like this one.
_spec = importlib.util.spec_from_file_location(
    "f1_plot", REPO_ROOT / "F1_domain_map" / "plot.py")
f1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(f1)

PROJECTION = f1.PROJECTION
DATA_TRANSFORM = f1.DATA_TRANSFORM


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Square 1080 x 1080 px: shown at full width in the LinkedIn feed on both
# desktop and mobile. The map is about 2:1, so it takes the middle band and the
# text sits above and below it.
FIG_SIZE_IN = 5.4
DPI = 200
MAP_RECT = [0.0, 0.25, 1.0, 0.54]      # left, bottom, width, height (figure)
TITLE_Y, SUBTITLE_Y, DATE_Y = 0.945, 0.9, 0.83
LEGEND_TOP_Y = 0.255
EXPLANATION_Y = 0.105

# Days over which a red particle fades from opaque to invisible.
FADE_DAYS = 5

# Particles.
MARKER_SIZE = 0.35                     # points^2; ~546,000 particles at 1080 px
DRIFTING_COLOR = mpl.colors.to_rgba("0.3", alpha=0.45)
SUCCESS_COLOR = mpl.colors.to_rgba("C2", alpha=1.0)
LOST_COLOR = mpl.colors.to_rgba("C3", alpha=1.0)

# Sea ice: transparent below the 15% advance threshold, then from a pale to a
# light blue with concentration. Kept light so particles read on top of it and
# distinct from the gray land.
SIC_THRESHOLD = 0.15
SIC_CMAP = LinearSegmentedColormap.from_list(
    "ice", [(0.0, "#e6f2fb"), (1.0, "#9cc8e8")])
SIC_CMAP.set_under((0, 0, 0, 0))
SIC_CMAP.set_bad((0, 0, 0, 0))

# Recruitment habitat limit (M5b criterion).
ISOBATH_DEPTH = 2000.0
ISOBATH_COLOR = "0.45"
ISOBATH_WIDTH = 0.4

# Region names in place of the CCAMLR codes. Placement reuses F1's anchors.
REGION_NAMES = {
    "48.1": "Antarctic\nPeninsula",
    "48.2": "South\nOrkney Is.",
    "48.3": "South\nGeorgia",
    "48.4": "South\nSandwich Is.",
    "48.5": "Weddell\nSea",
    "48.6N": "Bouvet\n(north)",
    "48.6S": "Bouvet\n(south)",
    "88.3": "Amundsen\nSea",
}
REGION_FONTSIZE = 5.5
REGION_BOX_ALPHA = 0.6

# Text.
FONT_SIZE = 8
TITLE = "Following virtual krill larvae to their first winter"
TITLE_SIZE = 12
SUBTITLE_SIZE = 8.5
DATE_SIZE = 15
LEGEND_SIZE = 8
EXPLANATION_SIZE = 7.5
# Most viewers will watch without sound and without reading the paper, so the
# success criterion is spelled out once, in plain words.
EXPLANATION = (
    "A larva is ready to overwinter if, when winter sea ice reaches it,\n"
    "it is developed enough and above the shelf or slope (shallower than 2000 m).\n"
    "On the first day, red marks larvae that could not be spawned because of sea ice."
)
LEGEND_LABELS = {
    "drifting": "Drifting larvae",
    "success": "Ready to overwinter",
    "lost": "Did not make it",
    "ice": "Sea ice",
    "isobath": "2000 m depth",
}

# Zorder stack, slotted around F1's: ice and isobath under the CCAMLR
# outlines (2), particles above them, all under land (3).
ZORDER_ICE = 1.2
ZORDER_ISOBATH = 1.6
ZORDER_DRIFTING = 2.4
ZORDER_LOST = 2.5
ZORDER_SUCCESS = 2.6

AGGREGATED = HERE / "data" / "aggregated.nc"
F1_AGGREGATED = REPO_ROOT / "F1_domain_map" / "data" / "aggregated.nc"
FRAMES_DIR = HERE / "frames"


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

class Cohort:
    """Lazy access to the movie aggregation, one day at a time."""

    def __init__(self, path: Path, stride: int):
        self.ds = xr.open_dataset(path)
        self.stride = stride
        sel = slice(None, None, stride)

        meanings = self.ds["outcome"].attrs["flag_meanings"].split()
        values = np.atleast_1d(self.ds["outcome"].attrs["flag_values"])
        code = dict(zip(meanings, values))
        outcome = self.ds["outcome"].values[sel]

        self.success = outcome == code["success"]
        self.exited = outcome == code["exited_domain"]
        self.lost = ~self.success & ~self.exited
        self.fate_day = self.ds["fate_day"].values[sel].astype(int)

        # Fate positions are fixed: project them once.
        xy = PROJECTION.transform_points(
            DATA_TRANSFORM,
            self.ds["final_lon"].values[sel].astype(np.float64),
            self.ds["final_lat"].values[sel].astype(np.float64))
        self.fate_xy = xy[:, :2]

        self.n_obs = self.ds.sizes["day"]
        self.release_date = pd.Timestamp(self.ds.attrs["release_date"])
        self.spawning_date = self.release_date - pd.Timedelta(
            days=int(self.ds.attrs["spawning_offset_days"]))
        self.sic_lon = self.ds["sic_lon"].values
        self.sic_lat = self.ds["sic_lat"].values

    def positions(self, day: int, mask: np.ndarray) -> np.ndarray:
        """Projected positions on `day` of the particles selected by `mask`."""
        sel = slice(None, None, self.stride)
        lon = self.ds["lon"].isel(day=day).values[sel][mask]
        lat = self.ds["lat"].isel(day=day).values[sel][mask]
        ok = np.isfinite(lon) & np.isfinite(lat)
        xy = PROJECTION.transform_points(DATA_TRANSFORM,
                                         lon[ok].astype(np.float64),
                                         lat[ok].astype(np.float64))
        return xy[:, :2]

    def sea_ice(self, day: int) -> np.ndarray:
        """Concentration as a fraction, NaN over land."""
        return self.ds["sic"].isel(day=day).values / 100.0


def frame_state(cohort: Cohort, frame: int):
    """
    Which particles are drawn how on a given frame.

    Frames past the last tracking day are hold frames: nothing moves, red
    particles keep fading, success stays.
    """
    age = frame - cohort.fate_day                     # days since the fate
    drifting = np.where(cohort.exited, age <= 0, age < 0)
    success = cohort.success & (age >= 0)
    lost = cohort.lost & (age >= 0) & (age < FADE_DAYS)
    lost_alpha = 1.0 - age[lost] / FADE_DAYS
    return drifting, success, lost, lost_alpha


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def add_region_labels(ax, ccamlr) -> None:
    """Region names at F1's label positions (centroids, 48.1 on its north edge)."""
    bbox = f1.make_label_bbox(alpha=REGION_BOX_ALPHA)

    def place(geom, code):
        if geom.is_empty or code not in REGION_NAMES:
            return
        anchor = f1.CCAMLR_LABEL_ANCHORS.get(code, "centroid")
        lat = geom.bounds[3] if anchor == "north" else geom.centroid.y
        ax.text(geom.centroid.x, lat, REGION_NAMES[code],
                transform=DATA_TRANSFORM, ha="center", va="center",
                fontsize=REGION_FONTSIZE, linespacing=1.0,
                bbox=bbox, zorder=10)

    for _, row in ccamlr.iterrows():
        code = row["GAR_Long_L"]
        if code == "48.6":
            north, south = f1.split_486(row.geometry)
            place(north, "48.6N")
            place(south, "48.6S")
        else:
            place(row.geometry, code)


def add_isobath(ax) -> None:
    ds = xr.open_dataset(F1_AGGREGATED)
    ax.contour(ds["nav_lon"].values, ds["nav_lat"].values,
               ds["bathymetry_masked"].values,
               levels=[ISOBATH_DEPTH], colors=[ISOBATH_COLOR],
               linewidths=ISOBATH_WIDTH, transform=DATA_TRANSFORM,
               zorder=ZORDER_ISOBATH)


def add_legend(fig) -> None:
    def dot(color, label):
        return Line2D([], [], linestyle="none", marker="o", markersize=5,
                      markerfacecolor=color, markeredgecolor="none",
                      label=label)

    handles = [
        dot(mpl.colors.to_rgba(DRIFTING_COLOR, 1.0), LEGEND_LABELS["drifting"]),
        dot(SUCCESS_COLOR, LEGEND_LABELS["success"]),
        dot(LOST_COLOR, LEGEND_LABELS["lost"]),
        Patch(facecolor=SIC_CMAP(0.5), edgecolor="none", label=LEGEND_LABELS["ice"]),
        Line2D([], [], color=ISOBATH_COLOR, linewidth=1.0,
               label=LEGEND_LABELS["isobath"]),
    ]
    fig.legend(handles=handles, loc="upper center",
               bbox_to_anchor=(0.5, LEGEND_TOP_Y),
               ncol=3, frameon=False, fontsize=LEGEND_SIZE,
               handletextpad=0.4, columnspacing=1.6)


def build_figure(cohort: Cohort):
    """Static layers once; returns the figure and the artists updated per frame."""
    mpl.rcParams.update({"font.size": FONT_SIZE})
    fig = plt.figure(figsize=(FIG_SIZE_IN, FIG_SIZE_IN), dpi=DPI)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes(MAP_RECT, projection=PROJECTION)
    f1.setup_polar_axes(ax)

    ccamlr = f1.load_ccamlr_geometries()
    f1.add_ccamlr_outlines(ax, ccamlr)
    f1.add_486_split_line(ax, ccamlr)
    f1.add_domain_boundary(ax)
    add_isobath(ax)
    add_region_labels(ax, ccamlr)

    ice = ax.pcolormesh(cohort.sic_lon, cohort.sic_lat,
                        np.ma.masked_invalid(cohort.sea_ice(0)),
                        cmap=SIC_CMAP, vmin=SIC_THRESHOLD, vmax=1.0,
                        shading="nearest", transform=DATA_TRANSFORM,
                        zorder=ZORDER_ICE, rasterized=True)

    def layer(color, zorder):
        return ax.scatter(np.empty(0), np.empty(0), s=MARKER_SIZE,
                          color=color, linewidths=0, transform=PROJECTION,
                          zorder=zorder, rasterized=True)

    drifting = layer(DRIFTING_COLOR, ZORDER_DRIFTING)
    lost = layer(LOST_COLOR, ZORDER_LOST)
    success = layer(SUCCESS_COLOR, ZORDER_SUCCESS)

    fig.text(0.5, TITLE_Y, TITLE, ha="center", va="center",
             fontsize=TITLE_SIZE, fontweight="bold")
    d = cohort.spawning_date
    fig.text(0.5, SUBTITLE_Y,
             f"Spawned around {d.day} {d:%B %Y}, followed for "
             f"{cohort.n_obs - 1} days with the ocean currents",
             ha="center", va="center", fontsize=SUBTITLE_SIZE, color="0.3")
    date_text = fig.text(0.5, DATE_Y, "", ha="center", va="center",
                         fontsize=DATE_SIZE, fontweight="bold")
    add_legend(fig)
    fig.text(0.5, EXPLANATION_Y, EXPLANATION, ha="center", va="center",
             fontsize=EXPLANATION_SIZE, color="0.3", linespacing=1.5)

    return fig, {"ice": ice, "drifting": drifting, "lost": lost,
                 "success": success, "date": date_text}


def update(cohort: Cohort, artists: dict, frame: int) -> None:
    day = min(frame, cohort.n_obs - 1)
    date = cohort.release_date + pd.Timedelta(days=day)
    artists["date"].set_text(f"{date.day} {date:%B %Y}")

    sic = cohort.sea_ice(day)
    artists["ice"].set_array(np.ma.masked_invalid(sic))

    drifting, success, lost, lost_alpha = frame_state(cohort, frame)
    artists["drifting"].set_offsets(
        cohort.positions(day, drifting) if frame < cohort.n_obs
        else np.empty((0, 2)))
    artists["success"].set_offsets(cohort.fate_xy[success])
    artists["lost"].set_offsets(cohort.fate_xy[lost])
    colors = np.tile(LOST_COLOR, (lost_alpha.size, 1))
    colors[:, 3] = lost_alpha
    artists["lost"].set_facecolor(colors)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_frames(spec: str | None, n_frames: int) -> range:
    if not spec:
        return range(n_frames)
    start, _, stop = spec.partition(":")
    return range(n_frames)[slice(int(start) if start else None,
                                 int(stop) if stop else None)]


def main():
    parser = argparse.ArgumentParser(description="Render the Movie S1 frames.")
    parser.add_argument("--frames", help="frame range, e.g. 0:50 (default: all)")
    parser.add_argument("--stride", type=int, default=1,
                        help="draw every Nth particle (default 1, all)")
    parser.add_argument("--data", type=Path, default=AGGREGATED)
    parser.add_argument("--out", type=Path, default=FRAMES_DIR)
    args = parser.parse_args()

    if not args.data.exists():
        raise FileNotFoundError(f"{args.data} not found. Run aggregate.py "
                                f"on the HPC first and copy it here.")

    cohort = Cohort(args.data, args.stride)
    n_frames = cohort.n_obs + FADE_DAYS
    frames = parse_frames(args.frames, n_frames)

    fig, artists = build_figure(cohort)
    args.out.mkdir(parents=True, exist_ok=True)
    for frame in frames:
        update(cohort, artists, frame)
        path = args.out / f"frame_{frame:03d}.png"
        fig.savefig(path, dpi=DPI)
        print(f"Wrote {f1.display_path(path)}")
    plt.close(fig)


if __name__ == "__main__":
    main()
