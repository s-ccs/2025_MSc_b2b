# scripts/run_single_case.jl
# Run one minimal simulation case for debugging the new src structure.
#
# Usage from the b2bsim_template directory:
#     julia --project=. scripts/run_single_case.jl
#
# Usage from the repository root:
#     julia --project=./b2bsim_template \
#         ./b2bsim_template/scripts/run_single_case.jl

using Random
using DataFrames
using Statistics
using UnfoldSim
using UnfoldDecode

include(joinpath(@__DIR__, "..", "src", "B2BSim.jl"))
using .B2BSim

# -----------------------------
# 1. Basic settings
# -----------------------------

sfreq = 20
seed = 42
n_trials = 400

collinearity_strength = 0.8
overlap_level = :high
noiselevel = 0.1

rng = MersenneTwister(seed)

println("\n=== Single-case B2B simulation sanity check ===")
println("seed = $seed")
println("sfreq = $sfreq")
println("n_trials = $n_trials")
println("collinearity_strength = $collinearity_strength")
println("overlap_level = $overlap_level")
println("noiselevel = $noiselevel")

# -----------------------------
# 2. Event design
# -----------------------------

base_design = make_base_design(;
    condition_name = :condition,
    low_level = "car",
    high_level = "face",
)

design = CollinearCovariateDesign(;
    design = base_design,
    n_trials = n_trials,
    collinearity_strength = collinearity_strength,
    condition_name = :condition,
    continuous_name = :continuous,
    low_level = "car",
    high_level = "face",
)

# Generate events once for sanity checks.
preview_events = generate_events(MersenneTwister(seed), design)

println("\n--- Event preview ---")
show(stdout, first(preview_events, 10); allcols = true)
println()

println("\n--- Condition-continuous sanity check ---")
actual_cor = condition_continuous_cor(preview_events)
println("actual cor(condition == face, continuous) = ", round(actual_cor, digits = 3))
show(stdout, summarize_condition_continuous(preview_events); allcols = true)
println()

# -----------------------------
# 3. Components and onset model
# -----------------------------

components = make_components(;
    sfreq = sfreq,
    use_positive_kernels = true,
)

onset_model = make_lognormal_onset_from_level(
    overlap_level;
    sfreq = sfreq,
    σ = 0.35,
    predictor_dependent = false,
)

noise = PinkNoise(; noiselevel = noiselevel)

# -----------------------------
# 4. Simulate continuous EEG
# -----------------------------

println("\n--- Simulating continuous EEG ---")

dat, evts = simulate(
    rng,
    design,
    components,
    onset_model,
    noise,
)

println("continuous data size = ", size(dat))
println("number of events = ", nrow(evts))
println("first event rows:")
show(stdout, first(evts, 10); allcols = true)
println()

# -----------------------------
# 5. Epoch simulated data
# -----------------------------

println("\n--- Epoching simulated data ---")

dat_e, times = UnfoldDecode.Unfold.epoch(
    dat,
    evts,
    [-0.1, 0.5],
    sfreq,
)

evts_e, dat_e = UnfoldDecode.Unfold.drop_missing_epochs(
    evts,
    dat_e,
)

println("epoched data size = ", size(dat_e))
println("number of complete epochs = ", nrow(evts_e))
println("time range = ", (first(times), last(times)))

println("\n--- Epoched event sanity check ---")
println("actual cor(condition == face, continuous) after epoch drop = ",
    round(condition_continuous_cor(evts_e), digits = 3))
show(stdout, summarize_condition_continuous(evts_e); allcols = true)
println()

# -----------------------------
# 6. Simple raw-epoch averages
# -----------------------------

condition_string = string.(evts_e.condition)
car_ix = condition_string .== "car"
face_ix = condition_string .== "face"

car_avg = vec(mean(dat_e[1, :, car_ix], dims = 2))
face_avg = vec(mean(dat_e[1, :, face_ix], dims = 2))

q25 = quantile(evts_e.continuous, 0.25)
q75 = quantile(evts_e.continuous, 0.75)
low_ix = evts_e.continuous .<= q25
high_ix = evts_e.continuous .>= q75

low_avg = vec(mean(dat_e[1, :, low_ix], dims = 2))
high_avg = vec(mean(dat_e[1, :, high_ix], dims = 2))

println("\n--- Raw average sanity values ---")
println("mean(face_avg - car_avg) = ", round(mean(face_avg .- car_avg), digits = 4))
println("mean(high_continuous_avg - low_continuous_avg) = ", round(mean(high_avg .- low_avg), digits = 4))

println("\nDone. This script only checks that the design, onset model, components, simulation, and epoching run together.")
