# src/OnsetModels.jl
# Onset and temporal-overlap utilities.
# This file is included inside the B2BSim module.

export DEFAULT_OVERLAP_TARGETS_MS,
       target_mean_ms_from_overlap,
       make_lognormal_onset,
       make_lognormal_onset_from_level

const DEFAULT_OVERLAP_TARGETS_MS = Dict(
    :low => 800.0,
    :medium => 400.0,
    :high => 200.0,
)

"""
    target_mean_ms_from_overlap(overlap_level; targets = DEFAULT_OVERLAP_TARGETS_MS)

Map an overlap label to a target mean inter-event interval in milliseconds.

Default interpretation:
- `:low`    -> 800 ms, weak temporal overlap
- `:medium` -> 400 ms, moderate temporal overlap
- `:high`   -> 200 ms, strong temporal overlap
"""
function target_mean_ms_from_overlap(
    overlap_level::Symbol;
    targets::Dict{Symbol, <:Real} = DEFAULT_OVERLAP_TARGETS_MS,
)
    haskey(targets, overlap_level) ||
        error("Unknown overlap_level: $overlap_level. Available levels: $(collect(keys(targets)))")

    return Float64(targets[overlap_level])
end

"""
    make_lognormal_onset(; sfreq = 20,
                           target_mean_ms = 250.0,
                           σ = 0.35,
                           predictor_dependent = false)

Create a LogNormal onset model for UnfoldSim.

By default, onset timing is independent of predictors:

    μ_formula = @formula(0 ~ 1)

This is the recommended clean setting for the first overlap simulation.
"""
function make_lognormal_onset(;
    sfreq::Real = 20,
    target_mean_ms::Real = 250.0,
    σ::Real = 0.35,
    predictor_dependent::Bool = false,
    condition_effect_on_onset::Real = 0.10,
    continuous_effect_on_onset::Real = 0.05,
)
    target_mean_samples = target_mean_ms / 1000 * sfreq

    target_mean_samples > 0 ||
        error("target_mean_ms and sfreq must imply a positive target_mean_samples")

    # LogNormal mean = exp(μ + σ^2 / 2), so this μ gives the desired mean.
    μ0 = log(target_mean_samples) - (σ^2) / 2

    if predictor_dependent
        return LogNormalOnsetFormula(
            μ_formula = @formula(0 ~ 1 + condition + continuous),
            μ_β = [μ0, condition_effect_on_onset, continuous_effect_on_onset],
            σ_β = [σ],
            offset_β = [0.0],
            truncate_upper = nothing,
        )
    else
        return LogNormalOnsetFormula(
            μ_formula = @formula(0 ~ 1),
            μ_β = [μ0],
            σ_β = [σ],
            offset_β = [0.0],
            truncate_upper = nothing,
        )
    end
end

"""
    make_lognormal_onset_from_level(overlap_level; sfreq = 20, σ = 0.35)

Convenience wrapper that creates a LogNormal onset model from `:low`, `:medium`, or `:high`.
"""
function make_lognormal_onset_from_level(
    overlap_level::Symbol;
    sfreq::Real = 20,
    σ::Real = 0.35,
    predictor_dependent::Bool = false,
)
    target_mean_ms = target_mean_ms_from_overlap(overlap_level)

    return make_lognormal_onset(;
        sfreq = sfreq,
        target_mean_ms = target_mean_ms,
        σ = σ,
        predictor_dependent = predictor_dependent,
    )
end
