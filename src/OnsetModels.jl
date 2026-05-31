function overlap_parameters(level::Symbol)
    level === :none && return (; target_mean_ms = 800.0, sigma = 0.20, offset_ms = 0.0)
    level === :high && return (; target_mean_ms = 250.0, sigma = 0.35, offset_ms = 0.0)
    error("Unknown overlap level: $level")
end

function make_onset_model(cfg::SimConfig)
    params = overlap_parameters(cfg.overlap_level)
    target_mean_samples = params.target_mean_ms / 1000.0 * cfg.sfreq
    offset_samples = params.offset_ms / 1000.0 * cfg.sfreq

    target_mean_samples > offset_samples || 
        error("Onset offset must be smaller than target mean interval.")

    μ0 = log(target_mean_samples - offset_samples) - (params.sigma^2) / 2

    return LogNormalOnsetFormula(
        μ_formula = @formula(0 ~ 1),
        μ_μ_β = [μ0],
        σ_β = [params.sigma],
        offset_β = [offset_samples],
        truncate_upper = nothing,
    )
end


function simulate_isis(cfg::SimConfig, rng = MersenneTwister(cfg.seed))
    design = make_design(cfg)
    onset_model = make_onset_model(cfg)
    events = UnfoldSim.generate_events(deepcopy(rng), design)
    onsets = UnfoldSim.simulate_interonset_distances(deepcopy(rng), onset_model, design)
    return (; events, onsets)
end


function summarize_isis(onsets)
    return DataFrame(;
        mean_isi = [mean(onsets)],
        std_isi = [std(onsets)],
        min_isi = [minimum(onsets)],
        max_isi = [maximum(onsets)],
    )
end