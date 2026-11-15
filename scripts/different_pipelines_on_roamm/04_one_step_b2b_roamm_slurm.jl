
#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using MScB2B
using CSV, DataFrames
using Serialization
using Dates
using LinearAlgebra

LinearAlgebra.BLAS.set_num_threads(parse(Int, ENV["SLURM_CPUS_PER_TASK"]))

# Configuration
eeg_dir = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/eeg"
events_path = "/home/st/st_us-051900/st_st181570/data/ROAMM/outputs/preprocessing/05_roamm_b2b_events.csv"

PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

models_dir = joinpath(
    PROJECT_ROOT,
    "slurm_results", "models", "roamm", "one_step_b2b"
)

mkpath(models_dir)

PREDICTORS = [:word_length, :frequency_zipf, :word_surprisal]
CROSS_VAL_REPS = 5
SCALING = :center
SFREQ = 256.0

# Load events and select subject
all_events = CSV.read(events_path, DataFrame)
subjects = sort(unique(all_events.subject_id))

task_id = parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
@assert 1 <= task_id <= length(subjects)

subject = subjects[task_id]

result_path = joinpath(
    models_dir,
    "$(subject)_$(SCALING)_one_step_b2b.jls"
)

if isfile(result_path)
    println("SKIPPED: $subject")
else
    t0 = time()
    temp_path = result_path * ".tmp"

    println("START: $subject")
    flush(stdout)

    try
        # 1. Concatenate continuous EEG runs
        eeg, evts = MScB2B.concatenate_runs(
            subject,
            all_events;
            eeg_dir = eeg_dir
        )

        println(
            "EEG size: $(size(eeg)) | Fixations: $(nrow(evts))"
        )
        flush(stdout)

        # 2. Fit one-step B2B
        result = MScB2B.fit_one_step_b2b_roamm(
            eeg,
            evts;
            predictors = PREDICTORS,
            sfreq = SFREQ,
            cross_val_reps = CROSS_VAL_REPS,
            scaling = SCALING
        )

        # 3. Save result
        serialize(temp_path, (
            subject = subject,
            predictors = PREDICTORS,
            scaling = SCALING,
            sfreq = SFREQ,
            cross_val_reps = CROSS_VAL_REPS,
            b2b_fit_one_step = result
        ))

        mv(temp_path, result_path; force = true)

        elapsed = round((time() - t0) / 60; digits = 2)

        println(
            "SUCCESS: $subject | Time: $(elapsed) min"
        )

    catch err
        println("FAILED: $subject")
        showerror(stdout, err, catch_backtrace())
        println()

        isfile(temp_path) && rm(temp_path)
        rethrow()
    end
end

println("Finished: $(now())")
