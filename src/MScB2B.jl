module MScB2B

using Random
using Statistics
using DataFrames
using StatsModels
using UnfoldSim

include("simulation.jl")

export SimulationConfig
export simulate_cases
export simulate_cases_continuous

end