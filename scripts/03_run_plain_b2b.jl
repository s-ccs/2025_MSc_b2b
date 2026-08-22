#!/usr/bin/env julia

using MScB2B
using Serialization
using Dates

const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(PROJECT_ROOT, "results")

mkpath(RESULTS_DIR)


# ----------------------------------------
# Load configuration
# ----------------------------------------
cfg = SimulationConfig(
    n_trials = 1000,
    sfreq = 100.0,
    n_channels = 20,
    
    noiselevel = 0.3,
    channel_noise_sd = 0.3,

    β0_n170 = 5.0,
    β_condition = 3.0,

    β0_p300 = 5.0,
    β_continuous = 1.0,

    rho = 0.8,

    overlap_interval_ms = 250.0,
    onset_condition_bias = -0.6,
    shift_onset = true
)

seed = 12
cross_val_reps = 3


println("===============================")
println("03 - plain B2B")
println("===============================")
println("Active project: ", Base.active_project())
println("Threads: ", Threads.nthreads())
println("Started: ", Dates.now())
println()

println("Configuration:")
println("  n_trials           = ", cfg.n_trials)
println("  noiselevel         = ", cfg.noiselevel)
println("  channel_noise_sd   = ", cfg.channel_noise_sd)
println("  rho                = ", cfg.rho)
println("  cross val reps       = ", cross_val_reps)
println("  seed               = ", seed)
println()

flush(stdout)

# -----------------------------------------
# Simulation    
# -----------------------------------------
println("Simulating cases...")
flush(stdout)

cases = MScB2B.simulate_cases(cfg; seed = seed)

println("Done simulating cases.")
flush(stdout)

# -----------------------------------------
# Run plain B2B
# -----------------------------------------
println()
println("Running plain B2B...")
flush(stdout)

elapsed = @elapsed begin
    results = MScB2B.run_plain_b2b(
        cases,
        cross_val_reps = cross_val_reps,
    )

    # -----------------------------------------
    # Save results
    # -----------------------------------------
    output_file = joinpath(RESULTS_DIR, "03_plain_b2b_results.jls")

    serialize(
        output_file,
        (
            cfg = cfg,
            seed = seed,
            cross_val_reps = cross_val_reps,
            score_tables = results.score_tables,
            fitted_models = results.fitted_models,
        )
    )

    println() 
    println("Saved to:")
    println(output_file)
end 

println()
println("Total time: ", round(elapsed / 60; digits = 2), " min")
println("Finished: ", Dates.now())
flush(stdout)