#!/bin/bash
#SBATCH --job-name=krico_ms1_series
#SBATCH --qos=nf
#SBATCH --array=0-31
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --output=logs/ms1_series_%A_%a.out

# ============================================================================
# MS1 per-year movie series — SLURM array driver
# ============================================================================
#
# One task per spawning year (1994 to 2025, array index 0 to 31). Each task
# runs the whole MS1 chain for the cohort released on the climatological peak
# day of its year (January 29, read from F2 by aggregate.py --year):
#
#   aggregate.py  ->  series/YYYY/aggregated.nc
#   plot.py       ->  series/YYYY/frames/frame_NNN.png   (split over the cores)
#   make_movie.sh ->  series/outcome_movie_YYYY.mp4
#
# The movies are identical in design to outcome_movie.mp4; only the cohort
# changes. The frames are rendered on the HPC rather than locally, which
# needs cartopy and geopandas from the python3 module, the Natural Earth
# 50m land and coastline files (downloaded once by cartopy, under a lock so
# the tasks do not race), and an ffmpeg binary (PATH, an ffmpeg module, or
# the one bundled with imageio-ffmpeg).
#
# Check all of that from a login node before queueing 32 tasks:
#   bash job_series.sh --preflight
#
# The SBATCH directives above are specific to ECMWF's Atos HPCF.
#
# Submit with:
#   cd MS1_outcome_movie
#   mkdir -p logs
#   sbatch job_series.sh                  # all 32 years
#   sbatch --array=16 job_series.sh       # one year: 1994 + 16 = 2010
#   sbatch --array=0-3,31 job_series.sh   # a subset
#
# A task that is re-run reuses an existing series/YYYY/aggregated.nc and
# re-renders the frames; REDO=1 sbatch job_series.sh re-aggregates too.
# ============================================================================

set -euo pipefail

FIRST_YEAR=1994        # spawning years as in F2's aggregation
LAST_YEAR=2025
SERIES_DIR="series"    # relative to this folder; gitignored

cd "$(dirname "${BASH_SOURCE[0]}")"

module load python3

export PYTHONUNBUFFERED=1
export MPLBACKEND=Agg          # no display on the compute nodes
export OMP_NUM_THREADS=1       # several plot.py processes per task, see below

# ----------------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------------

die() { echo "ERROR: $*" >&2; exit 1; }

check_env() {
    for var in KRICO_POST KRICO_RUNS KRICO_GLORYS12; do
        [[ -n "${!var:-}" ]] || die "${var} environment variable not set. See the repo README."
        [[ -d "${!var}" ]] || die "${var}=${!var} is not a directory."
    done
}

check_python() {
    python3 -c "import cartopy, geopandas, shapely, matplotlib, xarray, netCDF4, pandas" \
        || die "the python3 module lacks a package plot.py needs (see above)."
}

# Print the path of an ffmpeg binary, trying the PATH, an ffmpeg module, then
# the binary that imageio-ffmpeg bundles (pip install --user imageio-ffmpeg
# if none of these works).
find_ffmpeg() {
    if command -v ffmpeg >/dev/null 2>&1; then
        command -v ffmpeg; return 0
    fi
    module load ffmpeg >/dev/null 2>&1 || true
    if command -v ffmpeg >/dev/null 2>&1; then
        command -v ffmpeg; return 0
    fi
    python3 -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())" 2>/dev/null \
        && return 0
    return 1
}

# F1's base map draws Natural Earth 50m land and coastlines, which cartopy
# downloads on first use into its data directory. Done here once, under a
# lock, so that 32 tasks starting together do not all download (and clobber)
# the same files. A no-op once the files are present.
warm_cartopy() {
    mkdir -p "${SERIES_DIR}"
    flock "${SERIES_DIR}/.cartopy.lock" python3 - <<'EOF'
import cartopy.io.shapereader as shp
for name in ("land", "coastline"):
    shp.natural_earth(resolution="50m", category="physical", name=name)
EOF
}

preflight() {
    check_env
    check_python
    FFMPEG="$(find_ffmpeg)" || die "no ffmpeg found (PATH, module, or imageio-ffmpeg)."
    warm_cartopy
    echo "Preflight OK: ffmpeg = ${FFMPEG}"
}

if [[ "${1:-}" == "--preflight" ]]; then
    preflight
    exit 0
fi

[[ -n "${SLURM_ARRAY_TASK_ID:-}" ]] \
    || die "not an array task: submit with sbatch, or run 'bash job_series.sh --preflight'."

YEAR=$(( FIRST_YEAR + SLURM_ARRAY_TASK_ID ))
(( YEAR <= LAST_YEAR )) || die "array index ${SLURM_ARRAY_TASK_ID} is past ${LAST_YEAR}."

preflight

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------

WORK="${SERIES_DIR}/${YEAR}"
AGG="${WORK}/aggregated.nc"
FRAMES="${WORK}/frames"
MOVIE="${SERIES_DIR}/outcome_movie_${YEAR}.mp4"
mkdir -p "${WORK}" logs

echo "=== Spawning year ${YEAR} (task ${SLURM_ARRAY_TASK_ID} of job ${SLURM_ARRAY_JOB_ID:-?}) ==="

# ----------------------------------------------------------------------------
# 1. Aggregate
# ----------------------------------------------------------------------------

if [[ -f "${AGG}" && "${REDO:-0}" != "1" ]]; then
    echo "--- Aggregation: reusing ${AGG} (REDO=1 to recompute)"
else
    echo "--- Aggregation"
    # Written to a temporary name and renamed, so a task killed mid-write
    # cannot leave a truncated file that a re-run would then reuse.
    python3 aggregate.py --year "${YEAR}" --out "${AGG}.part"
    mv "${AGG}.part" "${AGG}"
fi

# ----------------------------------------------------------------------------
# 2. Frames, one plot.py process per core on disjoint frame ranges. The ranges
#    are independent by design (see plot.py), so this is a plain split.
# ----------------------------------------------------------------------------

echo "--- Frames"
N_FRAMES=$(python3 -c "import xarray as xr, sys; \
    print(xr.open_dataset(sys.argv[1], decode_times=False).sizes['day'])" "${AGG}")
N_PROC="${SLURM_CPUS_PER_TASK:-1}"
rm -rf "${FRAMES}"
mkdir -p "${FRAMES}"

pids=()
for (( i = 0; i < N_PROC; i++ )); do
    start=$(( i * N_FRAMES / N_PROC ))
    stop=$(( (i + 1) * N_FRAMES / N_PROC ))
    (( stop > start )) || continue
    python3 plot.py --data "${AGG}" --out "${FRAMES}" --frames "${start}:${stop}" \
        > "${WORK}/plot_${i}.log" 2>&1 &
    pids+=("$!")
done

failed=0
for pid in "${pids[@]}"; do
    wait "${pid}" || failed=1
done
if (( failed )); then
    cat "${WORK}"/plot_*.log >&2
    die "frame rendering failed for ${YEAR}; see ${WORK}/plot_*.log."
fi

n_png=$(find "${FRAMES}" -name 'frame_*.png' | wc -l)
(( n_png == N_FRAMES )) || die "expected ${N_FRAMES} frames, found ${n_png} in ${FRAMES}."
echo "    ${n_png} frames in ${FRAMES}"

# ----------------------------------------------------------------------------
# 3. Movie
# ----------------------------------------------------------------------------

echo "--- Movie"
FFMPEG="${FFMPEG}" FRAMES="$(realpath "${FRAMES}")" OUT="$(realpath -m "${MOVIE}")" \
    ./make_movie.sh

echo "=== Done: ${MOVIE} ==="
