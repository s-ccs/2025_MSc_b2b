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

# outer CV: evaluate decoding performance
nfolds = 2 

# inner CV: tune hyperparameters
ridge_tuning_nfolds = 2
ridge_tuning_resolution = 4

println("===============================")
println("02 - rerp decoding")
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
println("  outer nfolds       = ", nfolds)
println("  ridge nfolds       = ", ridge_tuning_nfolds)
println("  ridge resolution   = ", ridge_tuning_resolution)
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
# Ridge model
# -----------------------------------------
model = MScB2B.make_ridge_tuned_model(
    inner_nfolds = ridge_tuning_nfolds,
    resolution = ridge_tuning_resolution
)


# -----------------------------------------
# Run rerp decoding
# -----------------------------------------
println()
println("Running rerp decoding...")
flush(stdout)

elapsed = @elapsed begin 
    println("Starting condition decoding...")
    flush(stdout)

    condition_results = MScB2B.run_rerp_decoding(
        cfg,
        cases;
        model = model,
        target = :condition_num,
        nfolds = nfolds
    )

    println("Finished condition decoding.")
    flush(stdout)

    println("Starting continuous decoding...")
    flush(stdout)

    continuous_results = MScB2B.run_rerp_decoding(
        cfg,
        cases;
        model = model,
        target = :continuous,
        nfolds = nfolds
    )

    println("Finished continuous decoding.")
    flush(stdout)

    # -----------------------------------------
    # Save results
    # -----------------------------------------
    output_file = joinpath(RESULTS_DIR, "02_rerp_decoidng_results.jls")

    serialize(
        output_file,
        (
            cfg = cfg,
            seed = seed,
            nfolds = nfolds,
            ridge_tuning_nfolds = ridge_tuning_nfolds,
            ridge_tuning_resolution = ridge_tuning_resolution,
            condition_results = condition_results,
            continuous_results = continuous_results,
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
