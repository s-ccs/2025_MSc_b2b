#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using MScB2B
using CSV, DataFrames
using Serialization
using Dates
using LinearAlgebra

# Avoid CPU oversubscription
LinearAlgebra.BLAS.set_num_threads(1)

# Configuration
eeg_dir = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/eeg"
events_path = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/preprocessing 05_roamm_b2b_events.csv"
#eeg_dir = "/scratch/data/ROAMM/outputs/eeg"
#events_path = "/scratch/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"

PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
models_dir = joinpath(PROJECT_ROOT, "slurm_results", "models", "roamm", "two_step_b2b_evaluation_G_H")
logs_dir = joinpath(PROJECT_ROOT, "slurm_results", "logs", "roamm", "two_step_b2b_evaluation_G_H")
mkpath(models_dir)
mkpath(logs_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
CROSS_VAL_REPS = 5
SCALING = :center
OUTER_NFOLDS = 5
SEED = 12
SFREQ = 256.0


all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
@assert 1 <= task_id <= length(subjects)
subjects = [subjects[task_id]]


for subject in subjects
    result_path = joinpath(models_dir, "$(subject)_$(SCALING)_two_step_b2b_evaluation_G_H.jls")
    isfile(result_path) && (println("SKIPPED: $subject"); continue)


    t0 = time()
    println("START: $subject"); flush(stdout)

    eeg, evts = MScB2B.concatenate_runs(subject, all_events; eeg_dir = eeg_dir)
    evts[!, :event] = fill("fixation", nrow(evts))
    n_samples = size(eeg, 2)
    n_fixations = nrow(evts)
    println("DATA: $subject | Fixations: $n_fixations | Samples: $n_samples")


    # Decode each predictor
    result = MScB2B.predictive_evaluate_two_step_b2b_roamm(
        eeg,
        evts;
        predictors = PREDICTORS,
        sfreq = SFREQ,
        nfolds = OUTER_NFOLDS,
        cross_val_reps = CROSS_VAL_REPS,
        scaling = SCALING,
        seed = SEED
    )
    
    temp_path = result_path * ".tmp"

    serialize(temp_path, (
        subject = subject,
        predictors = PREDICTORS,
        scaling = SCALING,
        outer_nfolds = OUTER_NFOLDS,
        cross_val_reps = CROSS_VAL_REPS,
        seed = SEED, 
        evaluation = result
    ))
    

    mv(temp_path, result_path; force = true)
    elapsed = round(time() - t0 / 60; digits = 2)
    println("SUCCESS: $subject | Time: $(elapsed) min")
end 
println("Finished: $(now())")