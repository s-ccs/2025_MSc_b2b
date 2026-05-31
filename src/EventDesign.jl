UnfoldSim.@with_kw struct CollinearCovariateDesign <: UnfoldSim.AbstractDesign
    Design 
    n_trails::Int

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

Base.size(d::CollinearCovariateDesign) = (d.n_trails,)
UnfoldSim.size(d::CollinearCovariateDesign) = d.n_trails

function UnfoldSim.generate_events(
    rng::Random.AbstractRNG,
    d::CollinearCovariateDesign, 
)
    0.0 <= d.collinearity_strength <= 1.0 ||
        error("collinearity_strength must be within [0, 1].")

    base_n = size(d.design)[1]
    d.n_trails % base_n == 0 ||
        error("n_trials must be divisible by size(design) = $base_n")

    n_rep = div(d.n_trials, base_n)
    events = UnfoldSim.generate_events(deepcopy(rng), RepeatDesign(d.design, n_rep))

    is_high = string.(events[!, d.condition_name]) .== string(d.high_level) 
    separation = d.max_separation * d.collinearity_strength # The higher the collinearity strength, the more separation between the two conditions

    μ = ifelse.(
        is_high,
        d.center + separation / 2,
        d.center - separation / 2,
    )

    continuous = μ .+ d.noise_sd .* randn(rng, nrow(events))
    continuous = clamp.(continuous, d.lower, d.upper) 
    events[!, d.continuous_name] = continuous

    return events
end 

function collinearity_strength(level::Symbol)
    level === :none && return 0.0
    level === :low && return 0.4
    level === :high && return 0.8
    error("Unknown collinearity level: $level")
end


function make_base_design(cfg::SimConfig)
    validate(cfg) # cfg should have condition_levels field, which is a vector of condition levels (e.g., ["car", "face"])
    return SingleSubjectDesign(;
        conditions = Dict(:condition => collect(cfg.condition_levels)),
        event_order_function = shuffle,
    )
end


function make_design(cfg::SimConfig)
    return CollinearCovariateDesign(;
        design = make_base_design(cfg),
        n_trails = cfg.n_trials,
        collinearity_strength = collinearity_strength(cfg.collinearity_level),
        low_level = cfg.condition_levels[1],
        high_level = cfg.condition_levels[2],
        center = cfg.continuous_center,
        noise_sd = cfg.continuous_noise_sd,
        max_separation = cfg.continuous_max_separation,
        lower = cfg.continuous_bounds[1],
        upper = cfg.continuous_bounds[2],
    )
end


function make_events(cfg::SimConfig; rng = MersenneTwister(cfg.seed))
    return generate_events(rng, make_design(cfg))
end


function condition_continuous_cor(events; condition_name::Symbol = :condition)
    levels = unique(string.(events[!, condition_name]))
    length(levels) == 2 || error("Expected two condition levels, got $(levels).")
    high_level = levels[end]
    condition_numeric = Float64.(string.(events[!, condition_name])) .== high_level
    return cor(condition_numeric, events.continuous)
end

function summarize_events(events; condition_name::Symbol = :condition)
    summary = combine(
        groupby(events, condition_name),
        :continuous => mean => :continuous_mean,
        :continuous => std => :continuous_sd,
        nrow => :n,
    )
    summary[!, :condition_continuous_cor] .= condition_continuous_cor(events)
    return summary
end