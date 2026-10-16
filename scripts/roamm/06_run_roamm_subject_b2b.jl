#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MScB2B
using CSV
using DataFrames
using Unfold
using Random

Random.seed!(42)


const DATA_DIR =
    "/scratch/data/ROAMM/outputs/eeg"

const PREPROC_DIR =
    "/scratch/data/ROAMM/outputs/preprocessing"

const RESULTS_DIR =
    joinpath(
        @__DIR__,
        "..",
        "results",
        "roamm"
    )

const PREDICTORS = [:word_length]

const CROSS_VAL_REPS = 20

mkpath(RESULTS_DIR)


# ==============================
# Load event table
# ==============================

all_events = CSV.read(
    joinpath(
        PREPROC_DIR,
        "roamm_b2b_events.csv"
    ),
    DataFrame
)


# ==============================
# Subjects
# ==============================

# Run one subject for timing test
subjects = ["sub-10121"]

# Later, for all subjects:
# subjects = sort(unique(all_events.subject_id))


println(
    "Number of subjects: ",
    length(subjects)
)


# ==============================
# Batch loop
# ==============================

batch_start = time()


for subject in subjects

    subject_start = time()

    println()
    println("==============================")
    println("Processing ", subject)
    println("==============================")


    output_file = joinpath(
        RESULTS_DIR,
        "$(subject)_word_length_b2b.csv"
    )


    # ==============================
    # Allows restarting the batch
    # ==============================

    if isfile(output_file)
        println("Already finished, skipping.")
        continue
    end


    try

        # ==============================
        # Fit two-step B2B
        # ==============================

        fit_time = @elapsed begin

            fit =
                MScB2B.fit_two_step_b2b_roamm_subject(
                    subject,
                    all_events;
                    data_dir = DATA_DIR,
                    predictors = PREDICTORS,
                    sfreq = 256.0,
                    cross_val_reps = CROSS_VAL_REPS,
                )

        end


        println()
        println("------------------------------")
        println("Timing")
        println("------------------------------")

        println(
            "B2B fitting time: ",
            round(fit_time, digits=2),
            " s"
        )

        println(
            "B2B fitting time: ",
            round(fit_time / 60, digits=2),
            " min"
        )


        # ==============================
        # Convert result to DataFrame
        # ==============================

        result = DataFrame(
            Unfold.coeftable(
                fit.b2b_fit_corrected
            )
        )


        # Remove unused columns containing `nothing`
        select!(
            result,
            Not([:group, :stderror])
        )


        # ==============================
        # Store signed values
        # ==============================

        result[!, :estimate_signed] =
            copy(result.estimate)


        # ==============================
        # B2B magnitude
        # ==============================

        result.estimate .=
            abs.(result.estimate)


        # ==============================
        # Remove intercept
        # ==============================

        result = result[
            result.coefname .!= "(Intercept)",
            :
        ]


        # ==============================
        # Add metadata
        # ==============================

        result[!, :subject] =
            fill(
                subject,
                nrow(result)
            )


        result[!, :n_trials] =
            fill(
                fit.n_trials,
                nrow(result)
            )


        # ==============================
        # Save
        # ==============================

        CSV.write(
            output_file,
            result
        )


        println()
        println(
            "Saved: ",
            output_file
        )


        # ==============================
        # Total subject time
        # ==============================

        subject_time =
            time() - subject_start


        println()
        println("------------------------------")
        println("Subject summary")
        println("------------------------------")

        println(
            "Subject: ",
            subject
        )

        println(
            "Trials: ",
            fit.n_trials
        )

        println(
            "B2B fitting time: ",
            round(fit_time, digits=2),
            " s (",
            round(fit_time / 60, digits=2),
            " min)"
        )

        println(
            "Total subject time: ",
            round(subject_time, digits=2),
            " s (",
            round(subject_time / 60, digits=2),
            " min)"
        )


    catch err

        # Remove incomplete result file
        if isfile(output_file)
            rm(
                output_file;
                force = true
            )
        end


        println()
        println(
            "FAILED: ",
            subject
        )


        showerror(
            stdout,
            err,
            catch_backtrace()
        )

        println()

    end


    GC.gc()

end


# ==============================
# Batch timing
# ==============================

batch_time =
    time() - batch_start


println()
println("==============================")
println("Finished batch.")
println("==============================")

println(
    "Total batch time: ",
    round(batch_time, digits=2),
    " s"
)

println(
    "Total batch time: ",
    round(batch_time / 60, digits=2),
    " min"
)