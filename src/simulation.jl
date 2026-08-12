Base.@kwdef struct SimulationConfig
    n_trials::Int = 400
    sfreq::Float64 = 100.0

    n_channels::Int = 20
    noiselevel::Float64 = 0.1
    channel_noise_sd::Float64 = 0.1

    β0_n170::Float64 = 5.0
	β_condition::Float64 = 3.0
	
    β0_p300::Float64 = 5.0
	β_continuous::Float64 = 1.0
    
    rho::Float64 = 0.8

    overlap_interval_ms::Float64 = 250.0
    onset_condition_bias::Float64 = -0.6
    shift_onset::Bool = true
end

struct SimulationDesign <: UnfoldSim.AbstractDesign
    n_trials::Int
    confounded::Bool
    rho::Float64
end

UnfoldSim.size(design::SimulationDesign) = (design.n_trials,)
Base.length(design::SimulationDesign) = design.n_trials

function UnfoldSim.generate_events(
    rng::AbstractRNG,
    design::SimulationDesign,
)
    @assert iseven(design.n_trials)

    n_half = design.n_trials ÷ 2

    condition = vcat(fill("car", n_half), fill("face", n_half))

    shuffle!(rng, condition)

    condition_num = ifelse.(condition .== "face", 1.0, -1.0)

    if design.confounded
        continuous = 
            design.rho * condition_num .+ 
            sqrt(1 - design.rho^2) .* 
            randn(rng, design.n_trials)
    else
        continuous = randn(rng, design.n_trials)
    end

    continuous = (continuous .- mean(continuous)) ./ std(continuous)

    return DataFrame(
        condition = condition,
        condition_num = condition_num,
        continuous = continuous,
    )
end

function make_components(cfg::SimulationConfig)

    n170 = UnfoldSim.LinearModelComponent(
        basis = UnfoldSim.n170(; sfreq = cfg.sfreq,),
        formula = @formula(0 ~ 1 + condition),
        β = [cfg.β0_n170, cfg.β_condition],
        contrasts = Dict(),
    )

    p300 = UnfoldSim.LinearModelComponent(
        basis = UnfoldSim.p300(; sfreq = cfg.sfreq,),
        formula = @formula(0 ~ 1 + continuous),
        β = [cfg.β0_p300, cfg.β_continuous],
        contrasts = Dict(),
    )

    return [n170, p300]
end 

function make_onset(
    cfg::SimulationConfig;
    overlap::Bool = false,
    continuous::Bool = false,
)
    # No-overlap case
    if !overlap

        if !continuous
            return UnfoldSim.NoOnset()
        end
        # Continuous data connot use NoOnset().
        # Give events enough spacing that the components do not overlap.
        max_component_length = maximum([
            length(UnfoldSim.n170(; sfreq = cfg.sfreq)),
            length(UnfoldSim.p300(; sfreq = cfg.sfreq)),
        ])

        return UnfoldSim.UniformOnset(
            width = 1,
            offset = max_component_length + 1, 
        )
    end

    # Overlap case
    σ = 0.35

    target_mean_samples = cfg.overlap_interval_ms / 1000 * cfg.sfreq
    μ0 = log(target_mean_samples) - σ^2 / 2

    onset = UnfoldSim.LogNormalOnsetFormula(
        μ_formula = @formula(0 ~ 1 + condition),
        μ_β = [μ0,    # car/reference-level log-onset mean
            overlap ? cfg.onset_condition_bias : 0.0, ], # face minus car difference
        σ_β = [σ],
        offset_β = [0.0],
        truncate_upper = nothing,
    )
    
    if cfg.shift_onset
        return UnfoldSim.ShiftOnsetByOne(onset)
    else
        return onset
    end
end


function simulate_condition_events(
    cfg::SimulationConfig;
    confounded::Bool = false,
    overlap::Bool = false,
    seed::Int = 12,
)
    rng = MersenneTwister(seed)
    design = SimulationDesign(cfg.n_trials, confounded, cfg.rho,)
    components = make_components(cfg)
    onset = make_onset(cfg; overlap = overlap,)
    noise = UnfoldSim.PinkNoise(noiselevel = cfg.noiselevel,)

    dat, evts = UnfoldSim.simulate(rng, design, components, onset, noise; return_epoched = true,)
    evts[!, :event] = fill("stimulus", nrow(evts),)

    dat_3d = permutedims(repeat(dat, 1, 1, cfg.n_channels), (3, 1, 2),)
    dat_3d .+= cfg.channel_noise_sd .* randn(rng, size(dat_3d))

    return dat_3d, evts
end

function simulate_continuous_events(
    cfg::SimulationConfig;
    confounded::Bool = false,
    overlap::Bool = false,
    seed::Int = 12,
)
    rng = MersenneTwister(seed)
    design = SimulationDesign(cfg.n_trials, confounded, cfg.rho,)
    components = make_components(cfg)
    onset = make_onset(cfg; overlap = overlap, continuous = true,)
    noise = UnfoldSim.PinkNoise(noiselevel = cfg.noiselevel,)

    dat_1ch, evts = UnfoldSim.simulate(rng, design, components, onset, noise; return_epoched = false,)

    evts[!, :event] = fill("stimulus", nrow(evts),)

    dat_vector = vec(dat_1ch)

	# channels x continuous samples
	dat_cont = repeat(reshape(dat_vector, 1, :), cfg.n_channels, 1,)

	dat_cont .+= cfg.channel_noise_sd .* randn(rng, size(dat_cont))

	return dat_cont, evts 
end


function simulate_cases(
    cfg::SimulationConfig;
    continuous::Bool = false,
    seed::Int = 12,
)
    simulator = continuous ? simulate_continuous_events : simulate_condition_events

    return (
        clean = simulator(cfg; confounded = false, overlap = false, seed = seed,),
        overlap = simulator(cfg; confounded = false, overlap = true, seed = seed,),
        confounded = simulator(cfg; confounded = true, overlap = false, seed = seed,),
        both = simulator(cfg; confounded = true, overlap = true, seed = seed,),
    )
end
