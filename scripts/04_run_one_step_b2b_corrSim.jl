#!/usr/bin/env julia

using MScB2B
using Serialization

println("==============================")
println("Correlated continuous one-step B2B")
println("==============================")

cfg_corr = CorrelatedContinuousConfig(
    n_trials = 1500,
    sfreq = 100.0
)

println("Config:")
println(cfg_corr)

println("\nSimulating data...")

sim_corr = simulate_corr_cont_cases(
    cfg_corr;
    seed = 12
)

println("Simulation finished.")

println("\nRunning one-step B2B...")

one_step_corr = run_one_step_b2b(
    cfg_corr,
    sim_corr;
    cross_val_reps = 3
)

println("One-step B2B finished.")

result_dir = joinpath(
    "results",
    "final",
    "simulation_correlated_continuous"
)

mkpath(result_dir)

outfile = joinpath(
    result_dir,
    "one_step_b2b_1500trials_100hz.jls"
)

serialize(
    outfile,
    (
        cfg_corr = cfg_corr,
        sim_corr = sim_corr,
        one_step_corr = one_step_corr
    )
)

println("\nSaved to:")
println(outfile)

println("==============================")
println("Done")
println("==============================")