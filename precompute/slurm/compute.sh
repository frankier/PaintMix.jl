#!/bin/bash
#SBATCH --job-name=paintmix
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
set -euo pipefail
: "${PAINTMIX_WORK:?Set PAINTMIX_WORK to your isolated disk-backed workspace}"
: "${PAINTMIX_REPO:=$PAINTMIX_WORK/repo}"
export TMPDIR="$PAINTMIX_WORK/tmp"
export JULIA_DEPOT_PATH="$PAINTMIX_WORK/depot"
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export JULIA_NUM_PRECOMPILE_TASKS="${SLURM_CPUS_PER_TASK:-1}"
mkdir -p "$TMPDIR"
cd "$PAINTMIX_REPO"
exec "${PAINTMIX_JULIA:-$PAINTMIX_WORK/julia-1.13.3/bin/julia}" "$@"
