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
#eeg_dir = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/eeg"
# events_path = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/preprocessing 05_roamm_b2b_events.csv"
eeg_dir = "/scratch/data/ROAMM/outputs/eeg"
events_path = "/scratch/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"

PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
models_dir = joinpath(PROJECT_ROOT, "slurm_results", "models", "roamm", "two_step_b2b_evaluation")
logs_dir = joinpath(PROJECT_ROOT, "slurm_results", "logs", "roamm", "two_step_b2b_evaluation")
mkpath(models_dir)
mkpath(logs_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
CROSS_VAL_REPS = 5
SCALING = :center
OUTRTER_NFOLDS = 5
SEED = 12
SFREQ = 256.0


all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))

log_path = joinpath(logs_dir, "two_step_b2b_evaluation.log")
log_lock = ReentrantLock()
open(log_path, "a") do log
    logmsg(msg) = lock(() -> (println(log, msg); flush(log)), log_lock)
    logmsg("Started: $(now()) | Subjects: $(length(subjects)) | Threads: $(Threads.nthreads())")
    logmsg("Outer folds: $OUTER_NFOLDS | B2B reps: $CROSS_VAL_REPS | Seed: $SEED")

    Threads.@threads for i in eachindex(subjects)
        subject = subjects[i]
        result_path = joinpath(models_dir, "$(subject)_$(SCALING)_roamm_two_step_b2b_evaluation.jls")
        if isfile(result_path)
            logmsg("SKIPPED: $subject")
            continue
        end

        t0 = time()
        logmsg("START: $subject | Thread $(Threads.threadid())")

        try
            # 1. Concatenate continuous EEG runs
            eeg, evts = MScB2B.concatenate_runs(subject, all_events; eeg_dir = eeg_dir)
            evts[!, :event] = fill("fixation", nrow(evts))
            n_samples = size(eeg, 2)
            n_fixations = nrow(evts)
            logmsg("DATA: $subject | Epochs: $(nrow(evts_epoched)) | EEG size: $(size(dat_epoched))")


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
            
            # 5. Save results
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
            elapsed = round(time() - t0; digits = 2)
            logmsg("SUCCESS: $subject | Epochs: $(size(dat_epoched, 3)) | Time: $(elapsed)s")

        catch err 
            logmsg("FAILED: $subject")
            logmsg(sprint(showerror, err, catch_backtrace()))
        end
    end 

    logmsg("Finished: $(now())")
end