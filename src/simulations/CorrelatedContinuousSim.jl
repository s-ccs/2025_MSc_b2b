# 1. Configuration
Base.@kwdef struct CorrelatedContinuousConfig
    n_trials::Int = 1500
    sfreq::Float64 = 100.0
    rho::Float64 = 0.8

    n_channels::Int = 20
    noiselevel::Float64 = 0.3
    channel_noise_sd::Float64 = 0.3

    """
    β0_p100::Float64 = 5.0
    β_p100::Float64 = 2.0

    β0_n170::Float64 = 5.0
    β_n170::Float64 = 3.0

    β0_p300::Float64 = 5.0
    β_p300::Float64 = 1.5
    """

    component_width::Float64 = 0.15

    peak1::Float64 = 0.15
    peak2::Float64 = 0.45
    peak3::Float64 = 0.75

    β0::Float64 = 5.0
    β::Float64 = 2.0
    

    overlap_interval_ms::Float64 = 250.0
    onset_predictor_bias::Float64 = -0.6
    shift_onset::Bool = true
end

# 2. Design
Base.@kwdef struct CorrelatedContinuousDesign <: UnfoldSim.AbstractDesign
    n_trials::Int 
    rho::Float64 
end

UnfoldSim.size(design::CorrelatedContinuousDesign) = (design.n_trials,)
Base.length(design::CorrelatedContinuousDesign) = design.n_trials

function UnfoldSim.generate_events(
    rng::AbstractRNG,
    design::CorrelatedContinuousDesign
)
    ρ = design.rho
    # positive definite covariance matrix
    Σ = [
        1.0  ρ    ρ
        ρ    1.0  ρ
        ρ    ρ    1.0
    ]

    @assert isposdef(Symmetric(Σ))

    X = rand(
        rng,
        Distributions.MvNormal(zeros(3), Σ), # multivariate normal distribution
        design.n_trials
    )

    return DataFrame(
        continuous1 = X[1, :],
        continuous2 = X[2, :],
        continuous3 = X[3, :]
    )
end

function make_corr_cont_components(cfg::CorrelatedContinuousConfig)
    """
    p100 = UnfoldSim.LinearModelComponent(
        basis = UnfoldSim.p100(; sfreq = cfg.sfreq),
        formula = @formula(0 ~ 1 + continuous1),
        β = [cfg.β0_p100, cfg.β_p100],
        contrasts = Dict()
    )

    n170 = UnfoldSim.LinearModelComponent(
        basis = -UnfoldSim.n170(; sfreq = cfg.sfreq),
        formula = @formula(0 ~ 1 + continuous2),
        β = [cfg.β0_n170, cfg.β_n170],
        contrasts = Dict()
    )

    p300 = UnfoldSim.LinearModelComponent(
        # basis = UnfoldSim.p300(; sfreq = cfg.sfreq),
        basis = UnfoldSim.hanning(0.30, 0.35, cfg.sfreq),
        formula = @formula(0 ~ 1 + continuous3),
        β = [cfg.β0_p300, cfg.β_p300],
        contrasts = Dict()
    )

    """

    component1 = UnfoldSim.LinearModelComponent(
        basis = UnfoldSim.hanning(cfg.component_width, cfg.peak1, cfg.sfreq),
        formula = @formula(0 ~ 1 + continuous1),
        β = [cfg.β0, cfg.β],
    )

    component2 = UnfoldSim.LinearModelComponent(
        basis = UnfoldSim.hanning(cfg.component_width, cfg.peak2, cfg.sfreq),
        formula = @formula(0 ~ 1 + continuous2),
        β = [cfg.β0, cfg.β],
    )

    component3 = UnfoldSim.LinearModelComponent(
        basis = UnfoldSim.hanning(cfg.component_width, cfg.peak3, cfg.sfreq),
        formula = @formula(0 ~ 1 + continuous3),
        β = [cfg.β0, cfg.β],
    )
    return [component1, component2, component3]
    
    # return [p100, n170, p300]
end

function make_corr_cont_onset(
    cfg::CorrelatedContinuousConfig;
    overlap::Bool = false
)
    if !overlap
        max_component_length = maximum([
            length(UnfoldSim.p100(; sfreq = cfg.sfreq)),
            length(UnfoldSim.n170(; sfreq = cfg.sfreq)),
            length(UnfoldSim.p300(; sfreq = cfg.sfreq))
        ])

        return UnfoldSim.UniformOnset(
            width = 1,
            offset = max_component_length + 1
        )
    end

    # biased overlap
    σ = 0.35
    target_mean_samples = cfg.overlap_interval_ms / 1000 * cfg.sfreq

    μ0 = log(target_mean_samples) - σ^2 / 2
    onset = UnfoldSim.LogNormalOnsetFormula(
        μ_formula = @formula(0 ~ 1 + continuous1),
        μ_β = [
            μ0,
            cfg.onset_predictor_bias
        ],

        σ_β = [σ],
        offset_β = [0.0],
        truncate_upper = nothing
    )
    return cfg.shift_onset ? 
        UnfoldSim.ShiftOnsetByOne(onset) : onset
end

function simulator_corr_cont(
    cfg::CorrelatedContinuousConfig;
    overlap::Bool = false,
    seed::Int = 12,
    epoch_window = (-0.1, 1.0)
)
    rng = MersenneTwister(seed)
    design = CorrelatedContinuousDesign(cfg.n_trials, cfg.rho)
    components = make_corr_cont_components(cfg)
    onset = make_corr_cont_onset(cfg; overlap = overlap)
    noise = UnfoldSim.PinkNoise(noiselevel = cfg.noiselevel)

    # simulate continuous data (1 channel) and events
    dat_1ch, evts = UnfoldSim.simulate(rng, design, components, onset, noise; return_epoched = false)
    evts[!, :event] = fill("stimulus",nrow(evts))

    # make multichannel
    dat_cont = repeat(reshape(vec(dat_1ch), 1, :), cfg.n_channels, 1)
    dat_cont .+= cfg.channel_noise_sd .* randn(rng, size(dat_cont))

    # cut epochs from dat_cont
    dat_epoched, times = Unfold.epoch(dat_cont, evts, epoch_window, cfg.sfreq)
    evts_epoched, dat_epoched = Unfold.drop_missing_epochs(evts, dat_epoched)

    return (
        continuous = dat_cont,
        events_continuous = evts,
        epoched = dat_epoched,
        events_epoched = evts_epoched,
        times = times
    )
end

function simulate_corr_cont_cases(
    cfg::CorrelatedContinuousConfig;
    seed::Int = 12
)
    return (
        no_overlap = simulator_corr_cont(cfg; overlap = false, seed = seed),
        overlap = simulator_corr_cont(cfg; overlap = true, seed = seed)
    )
end