# src/EventDesign.jl
# Event design utilities for B2B simulation.
# This file is included inside the B2BSim module.

export make_base_design,
       CollinearCovariateDesign,
       condition_continuous_cor,
       summarize_condition_continuous

"""
    make_base_design(; condition_name = :condition,
                       low_level = "car",
                       high_level = "face",
                       event_order_function = shuffle)

Create a simple categorical design with one condition variable.

Default output has one column called `:condition` with two levels: `"car"` and `"face"`.
"""
function make_base_design(;
    condition_name::Symbol = :condition,
    low_level = "car",
    high_level = "face",
    event_order_function = shuffle,
)
    return SingleSubjectDesign(;
        conditions = Dict(condition_name => [low_level, high_level]),
        event_order_function = event_order_function,
    )
end

"""
    CollinearCovariateDesign

A custom UnfoldSim design that first generates a categorical design and then adds
one continuous covariate whose distribution depends on the categorical condition.

The parameter `collinearity_strength` controls how strongly the continuous
covariate differs between the two condition levels.

Important:
- `collinearity_strength` is not the exact Pearson correlation.
- It controls the mean separation between the continuous distributions.
"""
UnfoldSim.@with_kw struct CollinearCovariateDesign <: UnfoldSim.AbstractDesign
    design
    n_trials::Int

    collinearity_strength::Float64 = 0.0

    condition_name::Symbol = :condition
    continuous_name::Symbol = :continuous

    low_level = "car"
    high_level = "face"

    center::Float64 = 2.5
    noise_sd::Float64 = 1.0
    max_separation::Float64 = 3.0

    lower::Float64 = 0.0
    upper::Float64 = 5.0
end

Base.size(d::CollinearCovariateDesign) = (d.n_trials,)
UnfoldSim.size(d::CollinearCovariateDesign) = d.n_trials

"""
    UnfoldSim.generate_events(rng, d::CollinearCovariateDesign)

Generate events for a `CollinearCovariateDesign`.

The categorical condition is generated from `d.design`.
Then a continuous covariate is generated conditionally on the condition level:

- high level, by default `"face"`, receives a higher mean
- low level, by default `"car"`, receives a lower mean
"""
function UnfoldSim.generate_events(
    rng::Random.AbstractRNG,
    d::CollinearCovariateDesign,
)
    0.0 <= d.collinearity_strength <= 1.0 ||
        error("collinearity_strength must be between 0 and 1")

    base_n = size(d.design)[1]

    d.n_trials % base_n == 0 ||
        error("n_trials must be divisible by size(d.design) = $base_n")

    n_rep = div(d.n_trials, base_n)

    categorical_events = UnfoldSim.generate_events(
        deepcopy(rng),
        RepeatDesign(d.design, n_rep),
    )

    is_high =
        string.(categorical_events[!, d.condition_name]) .== string(d.high_level)

    separation = d.max_separation * d.collinearity_strength

    μ = ifelse.(
        is_high,
        d.center + separation / 2,
        d.center - separation / 2,
    )

    continuous = μ .+ d.noise_sd .* randn(rng, nrow(categorical_events))
    continuous = clamp.(continuous, d.lower, d.upper)

    categorical_events[!, d.continuous_name] = continuous

    return categorical_events
end

"""
    condition_continuous_cor(events; condition_name = :condition,
                              continuous_name = :continuous,
                              high_level = "face")

Compute the actual Pearson correlation between the categorical condition and the
continuous covariate. The condition is dummy-coded as 1 for `high_level` and 0 otherwise.
"""
function condition_continuous_cor(
    events;
    condition_name::Symbol = :condition,
    continuous_name::Symbol = :continuous,
    high_level = "face",
)
    condition_numeric =
        Float64.(string.(events[!, condition_name]) .== string(high_level))

    return cor(condition_numeric, events[!, continuous_name])
end

"""
    summarize_condition_continuous(events; condition_name = :condition,
                                    continuous_name = :continuous)

Return a small sanity-check table showing the mean, standard deviation, and
number of trials for the continuous covariate by condition.
"""
function summarize_condition_continuous(
    events;
    condition_name::Symbol = :condition,
    continuous_name::Symbol = :continuous,
)
    return combine(
        groupby(events, condition_name),
        continuous_name => mean,
        continuous_name => std,
        nrow => :n,
    )
end
