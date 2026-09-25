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

import Distributions

# Simulations
include("simulations/ConditionContinuousSim.jl")
include("simulations/CorrelatedContinuousSim.jl")

# Plotting
include("plotting/plot_decoding.jl")
include("plotting/plot_b2b.jl")
include("plotting/plot_correlated_decoding.jl")
include("plotting/plot_correlated_b2b.jl")


# Simulation pipelines
include("pipelines/01_standard_decoding.jl")
include("pipelines/02_rerp_decoding.jl")
include("pipelines/03_plain_b2b.jl")
include("pipelines/04_one_step_b2b.jl")
include("pipelines/05_two_step_b2b.jl")

# ROAMM 
include("roamm/two_step_b2b_roamm.jl")



export ConditionContinuousConfig
export CorrelatedContinuousConfig
export simulate_cond_cont_cases
export simulate_corr_cont_cases

export run_standard_decoding
export run_rerp_decoding
export run_plain_b2b
export run_one_step_b2b
export run_two_step_b2b

export plot_standard_decoding_grid
export plot_one_step_b2b_grid
export plot_correlated_results
export plot_correlated_decoding
end