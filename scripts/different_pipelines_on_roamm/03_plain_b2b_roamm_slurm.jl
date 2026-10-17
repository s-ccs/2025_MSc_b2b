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
logs_dir = joinpath(PROJECT_ROOT, "slurm_results", "logs", "roamm", "plain_b2b")

mkpath(models_dir)
mkpath(logs_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
CROSS_VAL_REPS = 5


all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
@assert 1 <= task_id <= length(subjects)
subjects = [subjects[task_id]]

log_path = joinpath(logs_dir, "roamm_$(subjects[1])_plain_b2b.log")
log_lock = ReentrantLock()
open(log_path, "a") do log
    logmsg(msg) = lock(() -> (println(log, msg); flush(log)), log_lock)
    logmsg("Started: $(now()) | Subjects: $(length(subjects)) | Cross_val_reps: $CROSS_VAL_REPS")
    #logmsg("Inner folds: $INNER_NFOLDS | Outer folds: $OUTER_NFOLDS | Resolution: $RESOLUTION | Seed: $SEED")

    for i in eachindex(subjects)
        subject = subjects[i]
        result_path = joinpath(models_dir, "$(subject)_plain_b2b.jls")
        if isfile(result_path)
            logmsg("SKIPPED: $subject")
            continue
        end

        t0 = time()
        temp_path = result_path * ".tmp"

        try
            # 1. Concatenate continuous EEG runs
            eeg, evts = MScB2B.concatenate_runs(subject, all_events; eeg_dir = eeg_dir)
            evts[!, :event] = fill("fixation", nrow(evts))
 
            dat_epoched, evts_epoched, times = MScB2B.prepare_roamm_epochs(eeg, evts)
            result = MScB2B.fit_plain_b2b_roamm(
                dat_epoched,
                evts_epoched,
                times;
                cross_val_reps = CROSS_VAL_REPS
            )

            serialize(temp_path, (
                subject = subject,
                predictors = PREDICTORS,
                cross_val_reps = CROSS_VAL_REPS,
                b2b_fit_plain = result
            ))
            
        
            mv(temp_path, result_path; force = true)
            elapsed = round((time() - t0) / 60; digits = 1)
            logmsg("SUCCESS: $subject | Fixations: $(nrow(evts))| Time: $(elapsed) min")

        catch err 
            logmsg("FAILED: $subject")
            logmsg(sprint(showerror, err, catch_backtrace()))
            isfile(temp_path) && rm(temp_path)
            rethrow()
        end
    end 

    logmsg("Finished: $(now())")
end