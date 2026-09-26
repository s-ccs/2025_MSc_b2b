### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 37ed115a-9582-11f1-bb1d-1dd80c25901c
begin

	using Pkg
	Pkg.activate(mktempdir())

	#Pkg.develop(path = "/home/xu/dev/UnfoldDecode_debug",)
	
	Pkg.add(url="https://github.com/unfoldtoolbox/UnfoldDecode.jl",rev="predict_type")
	Pkg.add(["UnfoldSim","UnfoldMakie","CairoMakie","Unfold","MLJ","MultivariateStats","MLJMultivariateStatsInterface","PlutoUI", "DataFrames", "StatsModels", "MLJLinearModels", 
			])
end


# ╔═╡ 76fc2129-60c2-4416-8e13-9381ad8ee05a
begin
	using Random
	using Statistics
	using DataFrames
	using StatsModels
	using PlutoUI
	using UnfoldDecode
	using LinearAlgebra
	using CairoMakie
	using UnfoldMakie
	using UnfoldSim
	using UnfoldMakie
	using CairoMakie
	using Unfold
end

# ╔═╡ 466fd4a3-b29d-443b-83b2-09acbc2ef590
# =====================
# Design struct
# =====================

begin
	UnfoldSim.@with_kw struct Decoding_Design <: UnfoldSim.AbstractDesign
		n_trials::Int = 400
		confounded::Bool = false
		rho::Float64 = 0.8
	end

	UnfoldSim.size(design::Decoding_Design) = (design.n_trials,)
	Base.length(design::Decoding_Design) = design.n_trials

	function UnfoldSim.generate_events(
		rng::UnfoldSim.AbstractRNG, 
		design::Decoding_Design,
	)
		@assert iseven(design.n_trials)
	
		n_half = div(design.n_trials, 2)
	
		condition = vcat(
			fill("car", n_half),
			fill("face", n_half),
		)
	
		condition = shuffle(rng, condition)
	
		# local numeric version, only used inside this function
		# not returned to event table
		condition_num = ifelse.(condition .== "face", 1.0, -1.0)
	
		if design.confounded
			continuous =
				design.rho .* condition_num .+
				sqrt(1 - design.rho^2) .* randn(rng, design.n_trials)
		else
			continuous = randn(rng, design.n_trials)
		end
	
		# z-score continuous
		continuous = (continuous .- mean(continuous)) ./ std(continuous)
		return DataFrame(
			condition = condition,
			condition_num = condition_num,
			continuous = continuous,
		)
	end
end

# ╔═╡ 5033a34c-f279-4080-ac72-b0b184c9acf5
# =============================
# Simulation functions
# ============================
begin
	# ------------------
	# onset function
	# ------------------
	
	function make_onset(cfg; overlap = false)
		if !overlap
			return UnfoldSim.NoOnset()
		end
	
		σ = 0.35
		target_mean_ms = 250.0
		target_mean_samples = cfg.onset_interval_ms / 1000 * cfg.sfreq
		μ0 = log(target_mean_samples) - σ^2 / 2
	
		return UnfoldSim.LogNormalOnsetFormula(
			μ_formula = @formula(0 ~ 1 + condition),
			μ_β = [μ0,    # car/reference-level log-onset mean
				   cfg.onset_condition_bias], # face minus car difference
			σ_β = [σ],
			offset_β = [0.0],
			truncate_upper = nothing,
		) |> ShiftOnsetByOne
	end

	function make_onset_continuous(
		cfg::Decoding_Config;
		overlap::Bool = false,
	)
		if !overlap
			max_component_length = maximum([
				length(UnfoldSim.n170(; sfreq = cfg.sfreq)),
				length(UnfoldSim.p300(; sfreq = cfg.sfreq)),
			])

			return UnfoldSim.UniformOnset(;
				width = 1,
				offset = max_component_length + 1,
			)
		end

		# with overlap
		σ = 0.35

		target_mean_samples = cfg.onset_interval_ms / 1000 * cfg.sfreq

		μ0 = log(target_mean_samples) - σ^2 / 2

		return UnfoldSim.LogNormalOnsetFormula(;
			μ_formula = @formula(0 ~ 1 + condition),

			μ_β = [
				μ0,
				cfg.onset_condition_bias,
			],

			σ_β = [σ],
			offset_β = [0.0],
			truncate_upper= nothing,	
		)
	end
	# --------------------
	# make components
	# -------------------
	
	function make_components(cfg::Decoding_Config)
		n1 = UnfoldSim.LinearModelComponent(;
			basis = UnfoldSim.n170(; sfreq = cfg.sfreq),
			formula = @formula(0 ~ 1 + condition),
			β = [cfg.β0_n170, cfg.β_condition],
			contrasts = Dict(),
		)
	
		p3 = UnfoldSim.LinearModelComponent(;
			basis = UnfoldSim.p300(; sfreq = cfg.sfreq),
			formula = @formula(0 ~ 1 + continuous),
			β = [cfg.β0_p300, cfg.β_continuous],
			contrasts = Dict(),
		)
	
		return [n1, p3]
	end
	
	# -------------------------------------
	# simulate cases for standard decoding
	# -------------------------------------
	
	function simulate_case(
		cfg::Decoding_Config;
		confounded::Bool = false,
		overlap::Bool = false,
		seed::Int = 12,
	)
		rng = MersenneTwister(seed)
	
		design = Decoding_Design(;
			n_trials = cfg.n_trials,
			confounded = confounded,
			rho = cfg.rho,
		)
	
		components = make_components(cfg)
	
		onset = make_onset(cfg; overlap = overlap)
	
		noise = UnfoldSim.PinkNoise(; noiselevel = cfg.noiselevel)
	
		dat, evts = UnfoldSim.simulate(
			rng,
			design,
			components,
			onset,
			noise;
			return_epoched = true,
		)
	
		dat_3d = permutedims(
			repeat(dat, 1, 1, cfg.n_channels),
			(3, 1, 2)
		)
		dat_3d .+= cfg.channel_noise_sd .* randn(rng, size(dat_3d))
	
		return dat_3d, evts
		# return channels x timepoints x trials
	end

	# ------------------------------
	# Continuous simulation for rERP
	# ------------------------------

	function simulate_case_continuous(
		cfg::Decoding_Config;
		confounded::Bool = false,
		overlap::Bool = false,
		seed::Int = 12, 
	)
		rng = MersenneTwister(seed)

		design = Decoding_Design(;
			n_trials = cfg.n_trials,
			confounded = confounded,
			rho = cfg.rho,
		)

		components = make_components(cfg)
		onset = make_onset_continuous(cfg; overlap = overlap)
		noise = UnfoldSim.PinkNoise(
			noiselevel = cfg.noiselevel,
		)

		dat_1ch, evts = UnfoldSim.simulate(
			rng,
			design,
			components,
			onset,
			noise;
			return_epoched = false,
		)

		evts[!, :event] = fill("stimulus", nrow(evts), )

		dat_vector = vec(dat_1ch)

		# channels x continuous samples
		dat_cont = repeat(
			reshape(dat_vector, 1, :),
			cfg.n_channels,
			1,
		)

		dat_cont .+= 
			cfg.channel_noise_sd .* randn(
				rng,
				size(dat_cont),
			)

		return dat_cont, evts 
		# return channels x continuous samples
	end
end

# ╔═╡ 463f3a3b-0111-4dc7-8f61-58295bc11303
cross_val_reps = 3

# ╔═╡ 2e08bc1f-3910-4fbd-a5b5-a582ef186f36
b2b_solver = (x, y) -> UnfoldDecode.solver_b2b(x, y; cross_val_reps=cross_val_reps);

# ╔═╡ 732598fb-cb29-41cb-adfd-4c7a982b0b63
function fit_b2b_only_cases(
    cases,
    design,
)
    result_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for (case_name, (dat, evts)) in pairs(cases)

        uf_b2b = Unfold.fit(
            UnfoldDecode.UnfoldModel,
            design,
            evts,
            dat;
            solver = b2b_solver,
        )

        result = coeftable(uf_b2b)

        result = result[
            result.coefname .!= "(Intercept)",
            :
        ]

        result[!, :case] .= String(case_name)

        push!(result_tables, result)
        fitted_models[case_name] = uf_b2b
    end

    return vcat(result_tables...), fitted_models
end

# ╔═╡ 111cb69b-997f-4f54-aad4-fa78e9b00be3
begin
	md"""
	### Simulation parameters

	**Mean inter-event interval:**
	$(@bind onset_interval_ms_slider PlutoUI.Slider(
		150.0:25.0:800.0;
		default = 250.0,
		show_value = true,
	))

	**Onset condition bias:**
	$(@bind biased_overlap_slider PlutoUI.Slider(
		-1.0:0.05:1.0;
		default = -0.6,
		show_value = true,
	))
	
	**Confound strength ρ:**  
	$(@bind rho_slider PlutoUI.Slider(
		0.0:0.05:0.95;
		default = 0.8,
		show_value = true,
	))

	**EEG noise level:**  
	$(@bind noiselevel_slider PlutoUI.Slider(
		0.0:0.05:1.0;
		default = 0.1,
		show_value = true,
	))

	**Channel noise:**  
	$(@bind channel_noise_slider PlutoUI.Slider(
		0.0:0.02:0.5;
		default = 0.1,
		show_value = true,
	))

	**Common ERP intercept β₀:**  
	$(@bind β0_slider PlutoUI.Slider(
		0.0:0.5:10.0;
		default = 5.0,
		show_value = true,
	))

	**Condition effect βcondition:**  
	$(@bind β_condition_slider PlutoUI.Slider(
		0.0:0.5:6.0;
		default = 3.0,
		show_value = true,
	))

	**Continuous effect βcontinuous:**  
	$(@bind β_continuous_slider PlutoUI.Slider(
		0.0:0.1:3.0;
		default = 1.0,
		show_value = true,
	))
	"""
end

# ╔═╡ 11b64a1a-9527-4470-b3c1-8d5da3e6924c
# =====================
# Config
# =====================

begin
	Base.@kwdef struct Decoding_Config
		n_trials::Int = 1000
		sfreq::Float64 = 100.0
	
		n_channels::Int = 20
		channel_noise_sd::Float64 = 0.1
		noiselevel::Float64 = 0.1
	
		β0_n170::Float64 = 5.0
		β_condition::Float64 = 3.0
		β0_p300::Float64 = 5.0
		β_continuous::Float64 = 1.0 
			
		rho::Float64 = 0.8

		onset_interval_ms::Float64 = 250.0
		onset_condition_bias::Float64 = -0.6
	end


	cfg = Decoding_Config(
		rho = rho_slider,
		onset_interval_ms = onset_interval_ms_slider,
		onset_condition_bias = biased_overlap_slider,
		noiselevel = noiselevel_slider,
		channel_noise_sd = channel_noise_slider,
		β0_n170 = β0_slider,
		β0_p300 = β0_slider,
		β_condition = β_condition_slider,
		β_continuous = β_continuous_slider,
		
	)
end

# ╔═╡ e21ea44a-20ad-41e5-8217-2dc21bacf622
# =================================
# Four cases for normal epoched data decoding
# =================================

begin
	standard_cases = (
		clean = simulate_case(
			cfg;
			overlap = false,
			confounded = false,
			seed = 12,
		),
		overlap = simulate_case(
			cfg;
			overlap = true,
			confounded = false,
			seed = 12,
		),
		confound = simulate_case(
			cfg;
			overlap = false,
			confounded = true,
			seed = 12,
		),
		both = simulate_case(
			cfg;
			overlap = true,
			confounded = true,
			seed = 12,
		),
	)

end 

# ╔═╡ c5ed5df3-4f34-46eb-927c-4ae04dbe5d3c
# ==========================
# rerp design
# =======================

des_rerp = des_rerp = [
    "stimulus" => (
        @formula(0 ~ 1 + condition + continuous),
        Unfold.firbasis(
            τ = [-0.1, 1.0],
            sfreq = cfg.sfreq,
            name = "stimulus",
        ),
    ),
]

# ╔═╡ 7fa27818-b03f-4001-a3ba-77ca71140d51
begin
    fo_b2b = @formula(0 ~ 1 + condition + continuous)

    n_timepoints_b2b = size(standard_cases.clean[1], 2)
    times_b2b = (0:n_timepoints_b2b-1) ./ cfg.sfreq

    des_b2b_plain = [
        Any => (
            fo_b2b,
            times_b2b,
        )
    ]
end

# ╔═╡ ea173ba6-1f55-4b11-b88d-ce21474d1bbd
b2b_only_scores, b2b_only_models =
    fit_b2b_only_cases(
        standard_cases,
        des_b2b_plain,
    )

# ╔═╡ 1d8a4dc5-f7dc-4f21-9927-d6748c0f4c0f
begin
	b2b_condition_scores =
	    b2b_only_scores[
	        string.(b2b_only_scores.coefname) .== "condition",
	        :
	    ]
	
	b2b_continuous_scores =
	    b2b_only_scores[
	        string.(b2b_only_scores.coefname) .== "continuous",
	        :
	    ]
end

# ╔═╡ f8e12b0c-61c8-4a0b-bd9b-0e627da42381
begin
    @show unique(b2b_only_scores.coefname)
    @show unique(b2b_only_scores.case)
    @show propertynames(b2b_only_scores)

    first(b2b_only_scores, 20)
end


# ╔═╡ f54e9a49-fb88-4414-94fb-ba1b4496a07f
function plot_b2b_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    score_col::Symbol = :estimate,
    use_abs::Bool = true,
    guide_peak::Real = 0.45,
)

    columns = propertynames(scores)

    :case in columns ||
        error("`scores` must contain a :case column.")

    :coefname in columns ||
        error("`scores` must contain a :coefname column.")

    score_col in columns ||
        error(
            "`scores` has no column $score_col. " *
            "Available columns: $columns"
        )

    time_col =
        if :time in columns
            :time
        elseif :timepoint in columns
            :timepoint
        else
            error("`scores` must contain :time or :timepoint.")
        end

    # ------------------------------------------
    # choose the coefficient rows for this target
    # ------------------------------------------

    coefnames = string.(scores.coefname)

    target_mask =
        if target == :continuous
            occursin.("continuous", coefnames)

        elseif target in (:condition, :condition_num)
            occursin.("condition", coefnames)

        else
            error("`target` must be :continuous or :condition")
        end

    target_scores = scores[target_mask, :]

    isempty(target_scores) &&
        error(
            "No rows found for target=$target. " *
            "Available coefnames: $(unique(coefnames))"
        )

    # ------------------------------------------
    # helper: extract one case curve
    # ------------------------------------------

    function get_curve(case_name::String)
        mask = string.(target_scores.case) .== case_name

        any(mask) ||
            error(
                "No case named \"$case_name\". " *
                "Available cases: $(unique(string.(target_scores.case)))"
            )

        time =
            if time_col == :time
                Float64.(target_scores[mask, :time])
            else
                (Float64.(target_scores[mask, :timepoint]) .- 1) ./ cfg.sfreq
            end

        values = Float64.(target_scores[mask, score_col])

        if use_abs
            values = abs.(values)
        end

        order = sortperm(time)
        return time[order], values[order]
    end

    # ------------------------------------------
    # timing guides
    # ------------------------------------------

    function scale_guide(effect)
        effect_abs = abs.(Float64.(effect))
        peak = maximum(effect_abs)

        return peak == 0 ?
               zeros(length(effect_abs)) :
               guide_peak .* effect_abs ./ peak
    end

    n170_effect =
        cfg.β_condition .*
        UnfoldSim.n170(; sfreq = cfg.sfreq)

    p300_effect =
        cfg.β_continuous .*
        UnfoldSim.p300(; sfreq = cfg.sfreq)

    n170_guide = scale_guide(n170_effect)
    p300_guide = scale_guide(p300_effect)

    n170_time = (0:length(n170_guide)-1) ./ cfg.sfreq
    p300_time = (0:length(p300_guide)-1) ./ cfg.sfreq

    # ------------------------------------------
    # labels and colors
    # ------------------------------------------

    condition_target = target in (:condition, :condition_num)

    figure_title =
        condition_target ?
        "Condition B2B" :
        "Continuous B2B"

    ylabel_text =
        if use_abs
            "B2B estimate magnitude"
        else
            "B2B estimate"
        end

    case_order = ("clean", "overlap", "confound", "both")

    case_titles = Dict(
        "clean" => "Clean",
        "overlap" => "Overlap",
        "confound" => "Confound",
        "both" => "Both",
    )

    colors = Dict(
        "clean" => :dodgerblue,
        "overlap" => :darkorange,
        "confound" => :seagreen,
        "both" => :deeppink,
    )

    # ------------------------------------------
    # figure
    # ------------------------------------------

    fig = Figure(
        size = (1000, 720),
        backgroundcolor = :white,
    )

    axes = Axis[]

    for (index, case_name) in enumerate(case_order)

        row = index <= 2 ? 1 : 2
        col = isodd(index) ? 1 : 2

        show_x = row == 2
        show_y = col == 1

        ax = Axis(
            fig[row, col];
            title = case_titles[case_name],
            xlabel = show_x ? "Time [s]" : "",
            ylabel = show_y ? ylabel_text : "",
            xticklabelsvisible = show_x,
            xticksvisible = show_x,
            yticklabelsvisible = show_y,
            yticksvisible = show_y,
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
        )

        time, values = get_curve(case_name)

        lines!(
            ax,
            time,
            values;
            color = colors[case_name],
            linewidth = 2.5,
        )

        # --------------------------------------
        # decide which guides to show
        # --------------------------------------

        confounded_case = case_name in ("confound", "both")

        show_n170 =
            if condition_target
                true
            else
                confounded_case
            end

        show_p300 =
            if condition_target
                confounded_case
            else
                true
            end

        if show_n170
            lines!(
                ax,
                n170_time,
                n170_guide;
                color = :black,
                linestyle = :dash,
                linewidth = 2.5,
            )
        end

        if show_p300
            lines!(
                ax,
                p300_time,
                p300_guide;
                color = :gray40,
                linestyle = :dot,
                linewidth = 2.5,
            )
        end

        hlines!(
            ax,
            [0.0];
            color = (:gray, 0.35),
            linewidth = 1,
        )

        push!(axes, ax)
    end

    linkxaxes!(axes...)
    linkyaxes!(axes...)

    Label(
        fig[0, 1:2],
        figure_title;
        fontsize = 24,
        font = :bold,
    )

    legend_elements = [
        LineElement(color = colors["clean"], linewidth = 2.5),
        LineElement(color = colors["overlap"], linewidth = 2.5),
        LineElement(color = colors["confound"], linewidth = 2.5),
        LineElement(color = colors["both"], linewidth = 2.5),
        LineElement(color = :black, linestyle = :dash, linewidth = 2.5),
        LineElement(color = :gray40, linestyle = :dot, linewidth = 2.5),
    ]

    legend_labels = [
        "Clean",
        "Overlap",
        "Confound",
        "Both",
        "N170 timing guide",
        "P300 timing guide",
    ]

    Legend(
        fig[3, 1:2],
        legend_elements,
        legend_labels;
        orientation = :horizontal,
        framevisible = false,
        nbanks = 2,
    )

    colgap!(fig.layout, 24)
    rowgap!(fig.layout, 16)

    return fig
end

# ╔═╡ c037c6af-9ec5-416f-be6a-53a0ad6660a9
fig_b2b_condition = plot_b2b_grid(
    b2b_only_scores,
    cfg;
    target = :condition,
	use_abs = false,
)

# ╔═╡ 86d0ac54-2e52-4c0f-bff9-06d3d8a0291c
fig_b2b_continuous = plot_b2b_grid(
    b2b_only_scores,
    cfg;
    target = :continuous,
	use_abs = false,
)

# ╔═╡ c766cc45-769b-4c25-b3a3-fd29a11761c3
begin
    clean_b2b = b2b_only_scores[
        b2b_only_scores.case .== "clean",
        :
    ]

    @show propertynames(clean_b2b)
    @show unique(clean_b2b.coefname)

    clean_b2b
end

# ╔═╡ 981f8e30-8e8a-490f-8008-8fd4f03be7b5
begin
    for coef_name in unique(string.(clean_b2b.coefname))

        tmp = clean_b2b[
            string.(clean_b2b.coefname) .== coef_name,
            :
        ]

        i = argmax(abs.(tmp.estimate))

        println(
            coef_name,
            " | peak time = ",
            round(tmp.time[i]; digits = 3),
            " | estimate = ",
            round(tmp.estimate[i]; digits = 3),
        )
    end
end

# ╔═╡ 38f5d3a8-2616-46e8-a297-d31453020f02
# ╠═╡ disabled = true
#=╠═╡
begin
    dat_clean, evts_clean = standard_cases.clean

    @show size(dat_clean)
    @show nrow(evts_clean)
    @show cor(evts_clean.condition_num, evts_clean.continuous)
end
  ╠═╡ =#

# ╔═╡ 48e6370b-8685-4641-ba68-0f3c792048ea
#=╠═╡
begin
    dat_clean, evts_clean = standard_cases.clean

    # average channels
    # time × trials
    Y_time_trial =
        dropdims(
            mean(dat_clean; dims = 1),
            dims = 1,
        )

    # trials × time
    Y = permutedims(Y_time_trial)

    # simple design:
    # intercept + condition(-1/+1) + continuous
    X = hcat(
        ones(nrow(evts_clean)),
        evts_clean.condition_num,
        evts_clean.continuous,
    )

    # predictor × time
    B = X \ Y

    times_debug = (0:size(Y, 2)-1) ./ cfg.sfreq

    fig = Figure()
    ax = Axis(
        fig[1, 1],
        xlabel = "Time [s]",
        ylabel = "OLS coefficient",
        title = "Sanity check BEFORE B2B",
    )

    lines!(
        ax,
        times_debug,
        B[2, :],
        label = "condition",
    )

    lines!(
        ax,
        times_debug,
        B[3, :],
        label = "continuous",
    )

    axislegend(ax)

    fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═37ed115a-9582-11f1-bb1d-1dd80c25901c
# ╠═76fc2129-60c2-4416-8e13-9381ad8ee05a
# ╠═11b64a1a-9527-4470-b3c1-8d5da3e6924c
# ╠═5033a34c-f279-4080-ac72-b0b184c9acf5
# ╠═466fd4a3-b29d-443b-83b2-09acbc2ef590
# ╠═e21ea44a-20ad-41e5-8217-2dc21bacf622
# ╠═c5ed5df3-4f34-46eb-927c-4ae04dbe5d3c
# ╠═463f3a3b-0111-4dc7-8f61-58295bc11303
# ╠═2e08bc1f-3910-4fbd-a5b5-a582ef186f36
# ╠═7fa27818-b03f-4001-a3ba-77ca71140d51
# ╠═732598fb-cb29-41cb-adfd-4c7a982b0b63
# ╠═ea173ba6-1f55-4b11-b88d-ce21474d1bbd
# ╠═1d8a4dc5-f7dc-4f21-9927-d6748c0f4c0f
# ╟─111cb69b-997f-4f54-aad4-fa78e9b00be3
# ╠═c037c6af-9ec5-416f-be6a-53a0ad6660a9
# ╠═86d0ac54-2e52-4c0f-bff9-06d3d8a0291c
# ╠═f8e12b0c-61c8-4a0b-bd9b-0e627da42381
# ╠═38f5d3a8-2616-46e8-a297-d31453020f02
# ╠═48e6370b-8685-4641-ba68-0f3c792048ea
# ╟─f54e9a49-fb88-4414-94fb-ba1b4496a07f
# ╠═c766cc45-769b-4c25-b3a3-fd29a11761c3
# ╠═981f8e30-8e8a-490f-8008-8fd4f03be7b5
