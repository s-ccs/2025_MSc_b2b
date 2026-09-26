#!/usr/bin/env julia

using MScB2B
using CSV
using DataFrames
using NPZ
using Statistics
using StatsModels: @formula
using Unfold
using UnfoldDecode
using Serialization


# ============================================================
# Config
# ============================================================

const DATA_DIR = "/scratch/data/ROAMM/outputs"
const EEG_DIR = joinpath(DATA_DIR, "eeg")

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(PROJECT_ROOT, "results", "roamm_two_step_lsq")
const CORRECTED_TRIALS_DIR = joinpath(RESULTS_DIR, "corrected_trials")

mkpath(RESULTS_DIR)
mkpath(CORRECTED_TRIALS_DIR)

const SFREQ = 256.0
const CROSS_VAL_REPS = 10


# ============================================================
# Step 1: continuous EEG -> corrected trials
# ============================================================

function correct_roamm_run_local(
    dat,
    evts;
    sfreq = SFREQ,
)
    evts = copy(evts)
    sort!(evts, :latency)

    evts[!, :event] = fill("stimulus", nrow(evts))

    design_rerp = [
        "stimulus" => (
            @formula(0 ~ 1 + surprisal_z),
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = sfreq,
            ),
        )
    ]

    uf_rerp = Unfold.fit(
        UnfoldLinearModelContinuousTime,
        design_rerp,
        evts,
        dat;
        eventcolumn = :event,
    )

    X_corrected = UnfoldDecode.singletrials(
        dat,
        uf_rerp,
        evts,
        "stimulus",
        :event,
    )

    times = Unfold.times(uf_rerp)[1]

    @assert size(X_corrected, 3) == nrow(evts)

    return (
        corrected_trials = X_corrected,
        events = evts,
        times = times,
    )
end


# ============================================================
# Step 2: corrected trials -> B2B
# ============================================================

function fit_subject_b2b_local(
    X_corrected,
    evts;
    times,
    cross_val_reps = CROSS_VAL_REPS,
)
    design_b2b = [
        "stimulus" => (
            @formula(0 ~ 1 + surprisal_z),
            times,
        )
    ]

    b2b_solver = (X, y) -> begin

        X_clean, y_clean =
            Unfold.drop_missing_epochs(X, y)

        X_clean = Float64.(X_clean)

        UnfoldDecode.solver_b2b(
            X_clean,
            y_clean;
            cross_val_reps = cross_val_reps,
            solver_G = UnfoldDecode.model_lsq,
            solver_H = UnfoldDecode.model_lsq,
            show_progress = false,
        )
    end

    return Unfold.fit(
        UnfoldModel,
        design_b2b,
        evts,
        X_corrected;
        solver = b2b_solver,
    )
end


function add_zscores!(df)

    μ = mean(df.surprisal_bits)
    σ = std(df.surprisal_bits)

    @assert σ > 0

    df[!, :surprisal_z] =
        (df.surprisal_bits .- μ) ./ σ

    return df
end


# ============================================================
# Load events
# ============================================================

evts_all = CSV.read(
    joinpath(DATA_DIR, "roamm_b2b_events.csv"),
    DataFrame,
)

# left eye only
filter!(:eye => ==("L"), evts_all)

dropmissing!(evts_all, :surprisal_bits)


# ============================================================
# Find EEG files
# ============================================================

eeg_pattern = r"^(sub-\d+)_run(\d+)_eeg\.npy$"

subject_runs =
    Dict{String, Vector{Tuple{Int,String}}}()

for filename in readdir(EEG_DIR)

    m = match(eeg_pattern, filename)

    isnothing(m) && continue

    subject = m.captures[1]
    run = parse(Int, m.captures[2])

    push!(
        get!(
            subject_runs,
            subject,
            Tuple{Int,String}[],
        ),
        (
            run,
            joinpath(EEG_DIR, filename),
        ),
    )
end

for runs in values(subject_runs)
    sort!(runs, by = first)
end


println("EEG subjects: ", length(subject_runs))
println("EEG runs: ", sum(length, values(subject_runs)))


# ============================================================
# TEST FIRST:
# only sub-10014
# ============================================================

# subjects = ["sub-10014"]

# Later, for all subjects, replace the line above with:
subjects = sort(collect(keys(subject_runs)))


all_scores = DataFrame[]


# ============================================================
# Main analysis
# ============================================================

for subject in subjects

    println()
    println("============================")
    println("Subject: ", subject)
    println("============================")

    evts_subject = filter(
        :subject_id => ==(subject),
        evts_all,
    )

    isempty(evts_subject) && continue

    # one z-score across all runs of this subject
    add_zscores!(evts_subject)

    X_runs = []
    evt_runs = DataFrame[]

    times = nothing

    for (run, eeg_file) in subject_runs[subject]

        evts_run = filter(
            :run_num => ==(run),
            evts_subject,
        )

        isempty(evts_run) && continue

        println(
            "run ", run,
            " | events = ", nrow(evts_run)
        )

        dat = npzread(eeg_file)

        println("EEG size: ", size(dat))

        corrected = correct_roamm_run_local(
            dat,
            evts_run,
        )

        println(
            "Corrected size: ",
            size(corrected.corrected_trials)
        )

        # save expensive overlap-corrected trials
        serialize(
            joinpath(
                CORRECTED_TRIALS_DIR,
                "$(subject)_run$(run)_corrected.jls",
            ),
            corrected,
        )

        push!(
            X_runs,
            corrected.corrected_trials,
        )

        push!(
            evt_runs,
            corrected.events,
        )

        if isnothing(times)
            times = corrected.times
        else
            @assert times ≈ corrected.times
        end
    end


    isempty(X_runs) && continue


    # combine corrected trials from all runs
    X_subject = cat(
        X_runs...;
        dims = 3,
    )

    evts_b2b = vcat(
        evt_runs...
    )

    @assert size(X_subject, 3) == nrow(evts_b2b)

    println(
        "Subject trials: ",
        size(X_subject, 3)
    )


    # subject-level B2B
    fit = fit_subject_b2b_local(
        X_subject,
        evts_b2b;
        times = times,
    )


    # extract B2B scores
    scores = DataFrame(
        Unfold.coeftable(fit)
    )

    scores[!, :estimate_signed] =
        copy(scores.estimate)

    scores.estimate .=
        abs.(scores.estimate)

    filter!(
        :coefname => !=("(Intercept)"),
        scores,
    )

    scores[!, :subject] =
        fill(subject, nrow(scores))

    push!(
        all_scores,
        scores,
    )


    # save subject result
    serialize(
        joinpath(
            RESULTS_DIR,
            "$(subject)_two_step_lsq.jls",
        ),
        (
            model = fit,
            scores = scores,
        ),
    )

    println("Finished subject: ", subject)
end


# ============================================================
# Save combined results
# ============================================================

@assert !isempty(all_scores) "No subject results were generated."

scores_all = vcat(all_scores...)

CSV.write(
    joinpath(
        RESULTS_DIR,
        "roamm_two_step_lsq_all_subjects.csv",
    ),
    select(
        scores_all,
        :subject,
        :coefname,
        :time,
        :estimate,
        :estimate_signed,
    ),
)

serialize(
    joinpath(
        RESULTS_DIR,
        "roamm_two_step_lsq_all_subjects.jls",
    ),
    scores_all,
)

println()
println("Finished.")