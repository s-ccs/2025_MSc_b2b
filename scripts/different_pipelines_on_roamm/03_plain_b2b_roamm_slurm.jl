#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using MScB2B
using CSV, DataFrames
using Serialization
using Dates
using LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(1)

# Configuration
eeg_dir = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/eeg"
events_path = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"
#eeg_dir = "/scratch/data/ROAMM/outputs/eeg"
#events_path = "/scratch/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"

PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
models_dir = joinpath(PROJECT_ROOT, "slurm_results", "models", "roamm", "plain_b2b")
mkpath(models_dir)


PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
CROSS_VAL_REPS = 5
SCALING = :center


all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
@assert 1 <= task_id <= length(subjects)
subjects = [subjects[task_id]]


for subject in subjects
    result_path = joinpath(models_dir, "$(subject)_$(SCALING)_roamm_plain_b2b.jls")
    isfile(result_path) && (println("SKIPPED: $subject"); continue)


    t0 = time()
    println("START: $subject"); flush(stdout)

    
    # 1. Concatenate continuous EEG runs
    eeg, evts = MScB2B.concatenate_runs(subject, all_events; eeg_dir = eeg_dir)
    evts[!, :event] = fill("fixation", nrow(evts))

    # prepare epochs
    dat_epoched, evts_epoched, times = MScB2B.prepare_roamm_epochs(eeg, evts)
    println("EEG size: $(size(eeg)) | Epochs: $(nrow(evts_epoched))")

    result = MScB2B.fit_plain_b2b_roamm(
        dat_epoched,
        evts_epoched,
        times;
        cross_val_reps = CROSS_VAL_REPS,
        scaling = SCALING
    )

    temp_path = result_path * ".tmp"

    serialize(temp_path, (
        subject = subject,
        predictors = PREDICTORS,
        scaling = SCALING,
        cross_val_reps = CROSS_VAL_REPS,
        b2b_fit_plain = result
    ))
    

    mv(temp_path, result_path; force = true)
    elapsed = round((time() - t0) / 60; digits = 1)
    println("SUCCESS: $subject | Time: $(elapsed) min")
end 
println("Finished: $(now())")