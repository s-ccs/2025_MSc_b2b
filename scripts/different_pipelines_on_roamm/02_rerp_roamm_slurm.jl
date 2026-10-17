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
models_dir = joinpath(PROJECT_ROOT, "slurm_results", "models", "roamm", "rerp")
logs_dir = joinpath(PROJECT_ROOT, "slurm_results", "logs", "roamm", "rerp")

mkpath(models_dir)
mkpath(logs_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
INNER_NFOLDS = 3
OUTER_NFOLDS = 5
RESOLUTION = 10
SEED = 12 

all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
@assert 1 <= task_id <= length(subjects)
subjects = [subjects[task_id]]

log_path = joinpath(logs_dir, "roamm_$(subjects[1])_rerp_decoding.log")
log_lock = ReentrantLock()
open(log_path, "a") do log
    logmsg(msg) = lock(() -> (println(log, msg); flush(log)), log_lock)
    logmsg("Started: $(now()) | Subjects: $(length(subjects))")
    #logmsg("Started: $(now()) | Subjects: $(length(subjects)) | Threads: $(Threads.nthreads())")
    #logmsg("Inner folds: $INNER_NFOLDS | Outer folds: $OUTER_NFOLDS | Resolution: $RESOLUTION | Seed: $SEED")

    for i in eachindex(subjects)
        subject = subjects[i]
        result_path = joinpath(models_dir, "$(subject)_roamm_rerp_decoding.jls")
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
 
            # 3. Tuned ridge model
            model = MScB2B.make_ridge_tuned_model(inner_nfolds = INNER_NFOLDS, resolution = RESOLUTION)

            # Decode each predictor
            result = MScB2B.fit_rerp_decoding_roamm(
                eeg,
                evts;
                model = model,
                targets = PREDICTORS,
                nfolds = OUTER_NFOLDS
            )
            # 5. Save results
            temp_path = result_path * ".tmp"

            serialize(temp_path, (
                subject = subject,
                predictors = PREDICTORS,
                inner_nfolds = INNER_NFOLDS,
                outer_nfolds = OUTER_NFOLDS,
                scores = result.scores
            ))
            
        
            mv(temp_path, result_path; force = true)
            elapsed = round((time() - t0) / 60; digits = 1)
            logmsg("SUCCESS: $subject | Fixations: $(nrow(evts))| Time: $(elapsed) min")

        catch err 
            logmsg("FAILED: $subject")
            logmsg(sprint(showerror, err, catch_backtrace()))
            rethrow()
        end
    end 

    logmsg("Finished: $(now())")
end