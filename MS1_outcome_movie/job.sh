#!/bin/bash
#SBATCH --job-name=krico_ms1
#SBATCH --qos=nf
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/ms1_%j.out

# ============================================================================
# MS1 movie aggregation — SLURM driver
# ============================================================================
#
# Extracts one cohort (trajectories, outcomes, daily sea ice) into
# data/aggregated.nc. Pass aggregate.py options after the script name, e.g.
#   sbatch job.sh --date 2006-01-29
#
# The SBATCH directives above are specific to ECMWF's Atos HPCF.
#
# Submit with:
#   cd MS1_outcome_movie
#   sbatch job.sh
# ============================================================================

set -euo pipefail

module load python3

for var in KRICO_POST KRICO_RUNS KRICO_GLORYS12; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: ${var} environment variable not set. See the repo README." >&2
        exit 1
    fi
    if [[ ! -d "${!var}" ]]; then
        echo "ERROR: ${var}=${!var} is not a directory." >&2
        exit 1
    fi
done

mkdir -p logs data
export PYTHONUNBUFFERED=1

python3 aggregate.py "$@"
