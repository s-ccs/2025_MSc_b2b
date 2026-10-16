using MScB2B
using CSV, DataFrames
using Serialization
using Dates
using Random

# ------------------------------------------
# Configuration
# ------------------------------------------

subject = "sub-10183"

eeg_dir = "/scratch/data/ROAMM/outputs/eeg"
events_path = "/scratch/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"
PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
models_dir = joinpath(PROJECT_ROOT, "results", "models", "roamm")
logs_dir = joinpath(PROJECT_ROOT, "results", "logs", "roamm")
mkpath(models_dir)
mkpath(logs_dir)

predictors = [
    :word_length,
    :frequency_zipf,
    :word_surprisal
]

scaling = :center
cross_val_reps = 5
sfreq = 256.0
seed = 42

# ------------------------------------------
# Logging
# ------------------------------------------

timestamp = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

log_path = joinpath(
    logs_dir,
     "$(subject)_$(scaling)_cv$(cross_val_reps)_memory_test.log"
)

result_path = joinpath(
    models_dir,
    "$(subject)_$(scaling)_cv$(cross_val_reps).jls"
)
@assert !isfile(result_path) "Result already exists: $result_path"
open(log_path, "w") do log
    # Print to both terminal and log file
    function logmsg(msg)
        println(msg)
        println(log, msg)
        flush(log)
    end

    t0 = time()

    logmsg("========================================")
    logmsg("ROAMM Two-step B2B")
    logmsg("========================================")
    logmsg("Subject: $subject")
    logmsg("Date: $(now())")
    logmsg("Scaling: $scaling")
    logmsg("Predictors: $predictors")
    logmsg("Sampling rate: $sfreq Hz")
    logmsg("Cross-validation repetitions: $cross_val_reps")
    logmsg("Random seed: $seed")
    logmsg("")

    try
        Random.seed!(seed)

        # Load events
        all_events = CSV.read(events_path, DataFrame)

        # Concatenate five runs
        logmsg("Loading and concatenating 5 runs...")

        eeg, evts = MScB2B.concatenate_runs(
            subject,
            all_events;
            eeg_dir = eeg_dir
        )

        logmsg("EEG shape: $(size(eeg))")
        logmsg("Fixations: $(nrow(evts))")

        # Fit pooled rERP + B2B
        logmsg("Starting pooled rERP + B2B...")

        result = MScB2B.fit_two_step_b2b_roamm(
            eeg,
            evts;
            predictors = predictors,
            scaling = scaling,
            sfreq = sfreq,
            cross_val_reps = cross_val_reps
        )
        
        # Save
        serialize(
            result_path,
            (
                subject = subject,
                scaling = scaling,
                predictors = result.predictors,
                model_predictors = result.model_predictors,
                corrected_trials_shape = result.corrected_trials_shape,
                cross_val_reps = cross_val_reps,
                sfreq = sfreq,
                seed = seed,
                rerp_model = result.rerp_model,
                b2b_fit_corrected = result.b2b_fit_corrected
            )
        )
        #saved = deserialize(result_path)

        #Unfold.coeftable(saved.rerp_model)
        #Unfold.coeftable(saved.b2b_fit_corrected)   
        
        logmsg("Saved result: $result_path")
        logmsg("Status: SUCCESS")

    catch err
        logmsg("Status: FAILED")
        logmsg("Error: $(sprint(showerror, err, catch_backtrace()))")
        rethrow()
    finally
        elapsed = round(time() - t0; digits = 2)
        logmsg("Elapsed time: $elapsed seconds")
    end
end


