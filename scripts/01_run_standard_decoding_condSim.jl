#!/usr/bin/env julia

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MScB2B
using Serialization
using Dates


# ============================================================
# Paths
# ============================================================

const PROJECT_ROOT =
    normpath(
        joinpath(@__DIR__, "..")
    )

const RESULTS_DIR =
    joinpath(
        PROJECT_ROOT,
        "results",
        "final",
        "simulation_condition_continuous"
    )

mkpath(RESULTS_DIR)


# ============================================================
# Configuration
# ============================================================

n_trials = 1500
sfreq = 100.0

outer_nfolds = 5
inner_nfolds = 3
ridge_resolution = 10

seed = 12


cfg_cond =
    MScB2B.ConditionContinuousConfig(
        n_trials = n_trials,
        sfreq = sfreq
    )


ridge_model =
    MScB2B.make_ridge_tuned_model(
        inner_nfolds = inner_nfolds,
        resolution = ridge_resolution
    )


# ============================================================
# Log header
# ============================================================

println("==============================")
println("Standard ridge decoding")
println("==============================")

println("Active project: ", Base.active_project())
println("Threads: ", Threads.nthreads())
println("Started: ", now())

println()

println("Configuration:")
println("  n_trials           = ", cfg_cond.n_trials)
println("  sfreq              = ", cfg_cond.sfreq)
println("  noiselevel         = ", cfg_cond.noiselevel)
println("  channel_noise_sd   = ", cfg_cond.channel_noise_sd)
println("  rho                = ", cfg_cond.rho)
println("  outer nfolds       = ", outer_nfolds)
println("  ridge nfolds       = ", inner_nfolds)
println("  ridge resolution   = ", ridge_resolution)
println("  seed               = ", seed)

println()


# ============================================================
# Run
# ============================================================

t_start = time()


println("Simulating cases...")

sim_cond =
    MScB2B.simulate_cond_cont_cases(
        cfg_cond;
        seed = seed
    )

println("Done simulating cases.")
println()


# ------------------------------------------------------------
# Condition decoding
# ------------------------------------------------------------

println("Running standard decoding...")
println("Starting condition decoding...")

standard_condition =
    MScB2B.run_standard_decoding(
        cfg_cond,
        sim_cond;
        model = ridge_model,
        target = :condition_num,
        nfolds = outer_nfolds,
        seed = seed
    )

println("Finished condition decoding.")


# ------------------------------------------------------------
# Continuous decoding
# ------------------------------------------------------------

println("Starting continuous decoding...")

standard_continuous =
    MScB2B.run_standard_decoding(
        cfg_cond,
        sim_cond;
        model = ridge_model,
        target = :continuous,
        nfolds = outer_nfolds,
        seed = seed
    )

println("Finished continuous decoding.")
println()


# ============================================================
# Save
# ============================================================

outfile =
    joinpath(
        RESULTS_DIR,
        "standard_decoding_results.jls"
    )


serialize(
    outfile,
    (
        cfg_cond = cfg_cond,
        sim_cond = sim_cond,

        standard_condition = standard_condition,
        standard_continuous = standard_continuous,

        outer_nfolds = outer_nfolds,
        inner_nfolds = inner_nfolds,
        ridge_resolution = ridge_resolution,

        seed = seed
    )
)


# ============================================================
# Finish
# ============================================================

elapsed_min =
    (time() - t_start) / 60


println("Saved to:")
println(outfile)

println()

println(
    "Total time: ",
    round(elapsed_min; digits = 2),
    " min"
)

println("Finished: ", now())