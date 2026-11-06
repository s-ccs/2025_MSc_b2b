#!/bin/bash
#SBATCH --job-name=roamm_rerp_rms_inner5
#SBATCH --account=st
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=36:00:00
#SBATCH --output=slurm_results/logs/roamm/rerp_rms_inner5/roamm_rerp_%A_%a.log
#SBATCH --error=slurm_results/logs/roamm/rerp_rms_inner5/roamm_rerp_%A_%a.log

set -e

# Go to the project root
cd "$SLURM_SUBMIT_DIR"

# Avoid CPU oversubscription
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

# julia 1.12.3
JULIA="$HOME/.juliaup/bin/julia"

echo "Job ID: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"

# Run the existing pipeline
"$JULIA" --threads=1 --project=. scripts/different_pipelines_on_roamm/02_rerp_roamm_slurm.jl

echo "Finished: $(date)"