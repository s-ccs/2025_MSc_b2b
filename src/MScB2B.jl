module MScB2B

using Random
using Statistics
using DataFrames
using StatsModels

using MLJ
using MLJLinearModels: RidgeRegressor
using Unfold
using UnfoldSim
using UnfoldDecode
using CairoMakie


include("simulation.jl")
include("Pipelines.jl")
include("plotting.jl")

export SimulationConfig
export simulate_cases

export run_standard_decoding
#export run_rerp_decoding
export run_plain_b2b
export run_debug_pipelines

export plot_standard_decoding_grid
end