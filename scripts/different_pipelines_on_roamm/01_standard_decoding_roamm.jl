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
#eeg_dir = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/eeg"
#events_path = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"
eeg_dir = "/scratch/data/ROAMM/outputs/eeg"
events_path = "/scratch/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"

PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
models_dir = joinpath(PROJECT_ROOT, "slurm_results", "models", "roamm", "standard_foldmean")
logs_dir = joinpath(PROJECT_ROOT, "slurm_results", "logs", "roamm", "standard_foldmean")

mkpath(models_dir)
mkpath(logs_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
INNER_NFOLDS = 3
OUTER_NFOLDS = 5
RESOLUTION = 10
SEED = 12 
SFREQ = 256.0

all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))

log_path = joinpath(logs_dir, "standard_decoding_foldmean.log")
log_lock = ReentrantLock()
open(log_path, "a") do log
    logmsg(msg) = lock(() -> (println(log, msg); flush(log)), log_lock)
    logmsg("Started: $(now()) | Subjects: $(length(subjects)) | Threads: $(Threads.nthreads())")
    logmsg("Inner folds: $INNER_NFOLDS | Outer folds: $OUTER_NFOLDS | Resolution: $RESOLUTION | Seed: $SEED")

    Threads.@threads for i in eachindex(subjects)
        subject = subjects[i]
        result_path = joinpath(models_dir, "$(subject)_standard_ridge_foldmean.jls")
        if isfile(result_path)
            logmsg("SKIPPED: $subject")
            continue
        end

        t0 = time()
        logmsg("START: $subject | Thread $(Threads.threadid())")

        try
            # 1. Concatenate continuous EEG runs
            eeg, evts = MScB2B.concatenate_runs(subject, all_events; eeg_dir = eeg_dir)

            # 2. Extract fixation-locked epochs
            dat_epoched, events_epoched, times = MScB2B.prepare_roamm_epochs(eeg, evts; sfreq = SFREQ)
            @assert size(dat_epoched, 3) == nrow(events_epoched)
            @assert size(dat_epoched, 2) == length(times)
 
            # 3. Tuned ridge model
            model = MScB2B.make_ridge_tuned_model(inner_nfolds = INNER_NFOLDS, resolution = RESOLUTION)

            # Decode each predictor
            result = MScB2B.fit_standard_decoding_roamm(
                dat_epoched,
                events_epoched,
                times;
                model = model,
                targets = PREDICTORS,
                nfolds = OUTER_NFOLDS,
                seed = SEED
            )
            @assert :r in propertynames(result.scores)
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
            elapsed = round(time() - t0; digits = 2)
            logmsg("SUCCESS: $subject | Epochs: $(size(dat_epoched, 3)) | Time: $(elapsed)s")

        catch err 
            logmsg("FAILED: $subject")
            logmsg(sprint(showerror, err, catch_backtrace()))
        end
    end 

    logmsg("Finished: $(now())")
end