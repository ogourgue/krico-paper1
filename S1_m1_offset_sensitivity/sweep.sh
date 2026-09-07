#!/bin/bash
#SBATCH --job-name=m1_sweep
#SBATCH --qos=nf
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/m1_sweep_%j.out

# ============================================================================
# M1 descent-ascent offset sweep — SLURM driver
# ============================================================================
#
# Computes M1 at each offset for every cohort in the 32-year hindcast, and
# verifies the offset-0 field lookup against the sea-ice concentration
# sampled along each trajectory at day 0.
#
# All five offsets are evaluated in a single pass: the cost is dominated by
# reading day 0 of ~3 848 trajectory files, while each additional offset is
# one lookup in an already-open monthly ice file. Splitting the offsets
# across array tasks would multiply the trajectory reads for no gain.
#
# The SBATCH directives above are specific to ECMWF's Atos HPCF.
#
# Submit with:
#   cd S1_m1_offset_sensitivity
#   sbatch sweep.sh
# ============================================================================

set -euo pipefail

module load python3

# Check that each variable is set and points somewhere real: an unset
# variable inside a path expansion yields a plausible-looking but wrong
# directory, which would otherwise surface as an import error after the job
# has queued.
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

if [[ ! -d "${KRICO_POST}/recruitment/krico_recruitment" ]]; then
    echo "ERROR: KRICO_POST=${KRICO_POST} does not look like a" >&2
    echo "krico-post-production clone (no recruitment/krico_recruitment)." >&2
    exit 1
fi

mkdir -p logs data

# M1 is computed with the same functions the pipeline uses, imported from
# krico_recruitment, so this analysis cannot drift from the classification it
# describes. KRICO_POST must therefore point at a krico-post-production
# clone, not merely at a directory of recruitment outputs.
export PYTHONPATH="${KRICO_POST}/recruitment:${PYTHONPATH:-}"
export PYTHONUNBUFFERED=1

python3 sweep.py \
    --runs    "${KRICO_RUNS}" \
    --glorys  "${KRICO_GLORYS12}" \
    --offsets 0 23 24 25 26 \
    --verify \
    --out     data/m1_offset_sweep.csv
