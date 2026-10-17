#!/bin/bash
#SBATCH --job-name=roamm_two_step
#SBATCH --account=st
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=roamm_%A_%a.out

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
"$JULIA" --threads=1 --project=. scripts/roamm/run_roamm_all_subjects.jl

echo "Finished: $(date)"