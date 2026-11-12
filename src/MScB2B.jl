module MScB2B

using Random
using LinearAlgebra
using Statistics
using DataFrames
using StatsModels
using MLJ
using MLJLinearModels: RidgeRegressor
using Unfold
using UnfoldSim
using UnfoldDecode
using UnfoldMakie
using CairoMakie
using NPZ
using Serialization

import Distributions

# Simulations
include("simulations/ConditionContinuousSim.jl")
include("simulations/CorrelatedContinuousSim.jl")

# Simulation plotting
include("plotting/simulation/plot_decoding.jl")
include("plotting/simulation/plot_b2b.jl")
include("plotting/simulation/plot_correlated_decoding.jl")
include("plotting/simulation/plot_correlated_b2b.jl")

# ROAMM plotting and results
include("plotting/roamm/load_results.jl")
include("plotting/roamm/plot_roamm_b2b.jl")
include("plotting/roamm/plot_roamm_comparisons.jl")
include("plotting/roamm/plot_roamm_decoding.jl")
include("plotting/roamm/plot_roamm_evaluation.jl")

# Data preparation
include("pipelines/00_prepare_roamm_data.jl")

# Five analysis pipelines (Simulation + ROAMM)
include("pipelines/01_standard_decoding.jl")
include("pipelines/02_rerp_decoding.jl")
include("pipelines/03_plain_b2b.jl")
include("pipelines/04_one_step_b2b.jl")
include("pipelines/05_two_step_b2b.jl")

# B2B evaluation
include("pipelines/B2B_R_evaluation.jl")

# Exports
export ConditionContinuousConfig, CorrelatedContinuousConfig
export simulate_cond_cont_cases, simulate_corr_cont_cases

export run_standard_decoding, run_rerp_decoding
export run_plain_b2b, run_one_step_b2b, run_two_step_b2b

export plot_standard_decoding_grid, plot_one_step_b2b_grid
export plot_correlated_results, plot_correlated_decoding

export concatenate_runs, prepare_roamm_epochs

export fit_standard_decoding_roamm
export fit_rerp_decoding_roamm
export fit_plain_b2b_roamm
export fit_one_step_b2b_roamm
export fit_two_step_b2b_roamm

export predictive_evaluate_plain_b2b_roamm
export predictive_evaluate_two_step_b2b_roamm
export calculate_b2b_R

end