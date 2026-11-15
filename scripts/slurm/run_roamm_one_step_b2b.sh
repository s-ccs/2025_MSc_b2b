#!/bin/bash
#SBATCH --job-name=roamm_one_step_b2b
#SBATCH --account=st
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=72:00:00
#SBATCH --output=slurm_results/logs/roamm/one_step_b2b/%A_%a.out
#SBATCH --error=slurm_results/logs/roamm/one_step_b2b/%A_%a.err

set -e

# Go to the project root
cd "$SLURM_SUBMIT_DIR"

# Avoid CPU oversubscription
export OPENBLAS_NUM_THREADS=8
export OMP_NUM_THREADS=1

# julia 1.12.3
JULIA="$HOME/.juliaup/bin/julia"

echo "Job ID: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"

# Run the existing pipeline
"$JULIA" --threads=1 --project=. scripts/different_pipelines_on_roamm/04_one_step_b2b_roamm_slurm.jl

echo "Finished: $(date)" 

# sbatch --array=1 --cpus-per-task=8 scripts/slurm/run_roamm_one_step_b2b.sh