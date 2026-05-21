# src/B2BSim.jl
# Main entry point for the B2B simulation helper code.

module B2BSim

using Random
using DataFrames
using Statistics
using StatsModels
using UnfoldSim
using UnfoldDecode

include("EventDesign.jl")
include("OnsetModels.jl")
include("Components.jl")

end # module B2BSim
