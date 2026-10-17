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
models_dir = joinpath(PROJECT_ROOT, "results", "models", "roamm")
logs_dir = joinpath(PROJECT_ROOT, "results", "logs", "roamm")
mkpath(models_dir)
mkpath(logs_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
SCALING = :center
CROSS_VAL_REPS = 5
SFREQ = 256.0

all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))
task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
@assert 1 <= task_id <= length(subjects)
subjects = [subjects[task_id]]

log_path = joinpath(
    logs_dir,
    "roamm_$(subjects[1])_$(SCALING)_cv$(CROSS_VAL_REPS).log"
)
log_lock = ReentrantLock()

open(log_path, "a") do log

    function logmsg(msg)
        lock(log_lock) do
            println(log, msg)
            flush(log)
        end
    end

    logmsg("\n========================================")
    logmsg("ROAMM Two-step B2B")
    logmsg("Started: $(now())")
    logmsg("Subjects: $(length(subjects))")
    logmsg("Predictors: $PREDICTORS")
    logmsg("Scaling: $SCALING")
    logmsg("Cross-validation repetitions: $CROSS_VAL_REPS")
    logmsg("Sampling rate: $SFREQ Hz")
    logmsg("Julia threads: $(Threads.nthreads())")
    logmsg("BLAS threads: $(LinearAlgebra.BLAS.get_num_threads())")
    logmsg("========================================")

    for i in eachindex(subjects)

        subject = subjects[i]

        result_path = joinpath(
            models_dir,
            "$(subject)_$(SCALING)_cv$(CROSS_VAL_REPS).jls"
        )

        if isfile(result_path)
            logmsg("SKIPPED: $subject")
            continue
        end

        t0 = time()
        logmsg("START: $subject | Thread $(Threads.threadid())")

        try
            eeg, evts = MScB2B.concatenate_runs(
                subject,
                all_events;
                eeg_dir = eeg_dir
            )

            n_samples = size(eeg, 2)
            n_fixations = nrow(evts)

            result = MScB2B.fit_two_step_b2b_roamm(
                eeg,
                evts;
                predictors = PREDICTORS,
                scaling = SCALING,
                sfreq = SFREQ,
                cross_val_reps = CROSS_VAL_REPS
            )

            temp_path = result_path * ".tmp"

            serialize(
                temp_path,
                (
                    subject = subject,
                    predictors = result.predictors,
                    scaling = SCALING,
                    cross_val_reps = CROSS_VAL_REPS,
                    sfreq = SFREQ,
                    b2b_fit_corrected = result.b2b_fit_corrected
                )
            )

            mv(temp_path, result_path; force = true)

            elapsed = round(time() - t0; digits = 2)

            logmsg(
                "SUCCESS: $subject | " *
                "Samples: $n_samples | " *
                "Fixations: $n_fixations | " *
                "Time: $(elapsed)s"
            )

        catch err
            logmsg("FAILED: $subject")
            logmsg(sprint(showerror, err, catch_backtrace()))
            rethrow()
        end
    end

    logmsg("========================================")
    logmsg("Batch finished: $(now())")
end