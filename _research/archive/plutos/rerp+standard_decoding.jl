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

# ╔═╡ f43ecec0-8a81-11f1-a925-3d063fa9b39c
begin

	using Pkg
	Pkg.activate(mktempdir())

	#Pkg.develop(path = "/home/xu/dev/UnfoldDecode_debug",)
	
	Pkg.add(url="https://github.com/unfoldtoolbox/UnfoldDecode.jl",rev="predict_type")
	Pkg.add(["UnfoldSim","UnfoldMakie","CairoMakie","Unfold","MLJ","MultivariateStats","MLJMultivariateStatsInterface","PlutoUI", "DataFrames", "StatsModels", "MLJLinearModels", 
			])
end


# ╔═╡ c8a40da4-af0c-470a-9227-2073e5bfc622
using MLJ, MultivariateStats, MLJMultivariateStatsInterface


# ╔═╡ 07238baa-e2bb-4c65-9f65-4506170dbd83
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

# ╔═╡ c9daa551-e9be-4ff2-96c0-86d857189650
RidgeRegressor = @load RidgeRegressor pkg=MLJLinearModels 

# ╔═╡ e0c3a358-844c-4ff4-931d-b83a512fc8cf
ridgeModel = RidgeRegressor()

# ╔═╡ 0a43fabc-9209-40da-b568-51f0ff423954
MLJ.machine(model::RidgeRegressor, X::AbstractMatrix{Float64},y::SubArray{Float64}; kwargs...,) = MLJ.machine(model,MLJ.table(X),y;kwargs...,)

# ╔═╡ ad6634be-abea-42b0-bd59-c24ffd40106f
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

# ╔═╡ b937fc2e-42bb-4ad4-83d2-bb29c00ccce0
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

# ╔═╡ bd65f918-f274-4d2c-a51d-dad55da44133
# ============================
# Standard decoding
# ============================

function decode_standard_ridge(
    dat,
    evts,
    model;
    target::Symbol = :continuous,
    nfolds::Int = 2,
    seed::Int = 12,
    sfreq::Float64 = 100.0,
)
    # dat: channels × timepoints × trials
    @assert ndims(dat) == 3
    @assert size(dat, 3) == nrow(evts)

    y = Float64.(evts[!, target])

    n_trials = size(dat, 3)
    n_timepoints = size(dat, 2)

    folds = UnfoldDecode.MLJBase.train_test_pairs(
        UnfoldDecode.MLJBase.CV(
            nfolds = nfolds,
            shuffle = true,
            rng = seed,
        ),
        1:n_trials,
    )

    # rows = trials, columns = timepoints
    yhat = fill(NaN, n_trials, n_timepoints)

    for (train_indices, test_indices) in folds

        # Reuse the function already implemented in UnfoldDecode
        machines = UnfoldDecode.fit_timepoints(
            model,
            @view(dat[:, :, train_indices]),
            @view(y[train_indices]),
        )

        predictions = UnfoldDecode.predict_timepoints(
            machines,
            @view(dat[:, :, test_indices]),
        )

        # predictions[t] contains predictions for one timepoint
        for t in eachindex(predictions)
            yhat[test_indices, t] .= Float64.(predictions[t])
        end
    end

    scores = DataFrame(
        timepoint = collect(1:n_timepoints),
        time = (collect(1:n_timepoints) .- 1) ./ sfreq,

        r = [
            cor(yhat[:, t], y)
            for t in 1:n_timepoints
        ],

        r_squared = [
            RSquared()(yhat[:, t], y)
            for t in 1:n_timepoints
        ],
    )

    return scores, yhat
end

# ╔═╡ ac6cb027-550b-4a93-a2d4-2ec3ad3325ac
function decode_standard_cases(
	cases,
	model;
	target::Symbol = :continuous,
	nfolds::Int = 2,
	seed::Int = 2,
	sfreq::Real = 100.0,
)

	score_tables = DataFrame[]
	yhats = Dict{Symbol, Matrix{Float64}}()

	for(case_name, (dat, evts)) in pairs(cases)
		score_case, yhat_case = decode_standard_ridge(
			dat,
			evts,
			model;
			target = target,
			nfolds = nfolds,
			seed = seed,
			sfreq = sfreq,
		)

		score_case.case = fill(String(case_name), nrow(score_case))

		push!(score_tables, score_case)
		yhats[case_name] = yhat_case
	end

	return vcat(score_tables...), yhats
end

# ╔═╡ b0d97bbf-6e1b-496d-8ac7-0dc9688cc309
begin

	# ========================================
	# Standard_decoding_plot
	# ========================================

	function _sd_get_curve(
		scores::DataFrame,
		case_name::AbstractString;
		sfreq::Real,
	)
		required_columns = Set([:timepoint, :r, :case])
		available_columns = Set(propertynames(scores))

		issubset(required_columns, available_columns) ||
			error(
				"`scores` must contain :timepoint, :r and :case. " *
				"Available columns: $(propertynames(scores))"
			)

		case_mask =
			string.(scores.case) .== case_name


		# rERP results have a real :time column
		if :time in propertynames(scores)
			sub = scores[
				case_mask,
				[:time, :r],
			]

			isempty(sub) && 
				error("No rows found for case = \"$case_name\".")

			sort!(sub, :time)

			return Float64.(sub.time), Float16.(sub.r)

		# Standard decoding uses :timepoint
		elseif :timepoint in propertynames(scores)
			
			sub = scores[
				case_mask,
				[:timepoint, :r],
			]
	
			isempty(sub) &&
				error(
					"No rows found for case = \"$case_name\". " *
					"Available cases: $(unique(string.(scores.case)))"
				)
	
			sort!(sub, :timepoint)

			time = 
				(Float64.(sub.timepoint) .- 1.0) ./ sfreq

			return time, Float64.(sub.r)

		else 
			error(
				"`scores` must contain either :time or :timepoint."
			)
		end 
		
		# Convert sample indices to seconds.
		# timepoint 1 becomes time 0 s.
		time =
			(Float64.(sub.timepoint) .- 1.0) ./ sfreq

		values =
			Float64.(sub.r)

		return time, values
	end


	# ============================================================
	# Create the corresponding ground-truth effect
	# ============================================================

	function _sd_true_effect(
		cfg,
		target::Symbol,
		n_timepoints::Integer;
		target_peak::Real = 0.45,
		use_magnitude::Bool = true,
	)
		raw_effect, label =
			if target == :condition
				(
					cfg.β_condition .*
					UnfoldSim.n170(; sfreq = cfg.sfreq),

					use_magnitude ?
					"True N170 condition-effect magnitude (scaled)" :
					"True N170 condition effect (scaled)",
				)

			elseif target == :continuous
				(
					cfg.β_continuous .*
					UnfoldSim.p300(; sfreq = cfg.sfreq),

					use_magnitude ?
					"True P300 continuous-effect magnitude (scaled)" :
					"True P300 continuous effect (scaled)",
				)

			else
				error(
					"Unknown target: $target. " *
					"Use :condition or :continuous."
				)
			end

		raw_effect =
			Float64.(raw_effect)

		# For decoding correlation, the temporal magnitude is generally
		# more directly comparable than the signed ERP polarity.
		effect_for_plot =
			use_magnitude ?
			abs.(raw_effect) :
			raw_effect

		# Protect against different waveform/data lengths.
		n =
			min(
				Int(n_timepoints),
				length(effect_for_plot),
			)

		effect_for_plot =
			effect_for_plot[1:n]

		time =
			(0:n-1) ./ cfg.sfreq

		maximum_absolute =
			maximum(abs.(effect_for_plot))

		scaled_effect =
			maximum_absolute == 0 ?
			zeros(Float64, n) :
			Float64(target_peak) .*
			effect_for_plot ./ maximum_absolute

		return time, scaled_effect, label
	end


	# ============================================================
	# Make one consistently styled axis
	# ============================================================

	function _sd_make_axis(
		position;
		title::AbstractString,
		ylabel::AbstractString,
		show_x::Bool,
		show_y::Bool,
		x_limits,
		y_limits,
		x_ticks,
	)
		ax = Axis(
			position;

			title = title,
			titlesize = 16,
			titlefont = :bold,
			titlegap = 8,

			xlabel = show_x ? "Time [s]" : "",
			ylabel = show_y ? ylabel : "",

			xlabelsize = 14,
			ylabelsize = 14,

			xticklabelsize = 12,
			yticklabelsize = 12,

			xticks = x_ticks,

			xgridvisible = false,
			ygridvisible = false,

			topspinevisible = false,
			rightspinevisible = false,

			xticklabelsvisible = show_x,
			xticksvisible = show_x,
			xlabelvisible = show_x,

			yticklabelsvisible = show_y,
			yticksvisible = show_y,
			ylabelvisible = show_y,

			backgroundcolor = :white,
		)

		xlims!(ax, x_limits...)
		ylims!(ax, y_limits...)

		# Zero-correlation reference line
		hlines!(
			ax,
			[0.0];
			color = (:gray45, 0.55),
			linestyle = :dash,
			linewidth = 1.0,
		)

		return ax
	end


	# ============================================================
	# Main plotting function
	# ============================================================

	function plot_standard_decoding_grid(
		scores::DataFrame,
		cfg;
		target::Symbol,
		true_effect_peak::Real = 0.45,
		use_effect_magnitude::Bool = true,
	)
		# --------------------------------------------------------
		# Target-specific labels
		# --------------------------------------------------------

		figure_title,
		ylabel_text,
		bottom_left_title =
			if target == :condition
				(
					"Condition decoding",
					"corr(predicted condition, true condition)",
					"Confound and true condition effect",
				)

			elseif target == :continuous
				(
					"Continuous decoding",
					"corr(predicted continuous, true continuous)",
					"Confound and true continuous effect",
				)

			else
				error(
					"Unknown target: $target. " *
					"Use :condition or :continuous."
				)
			end


		# --------------------------------------------------------
		# Fixed visual identity
		#
		# Only `both` is dashed.
		# --------------------------------------------------------

		case_colors = Dict(
			"clean"    => :dodgerblue,
			"overlap"  => :darkorange,
			"confound" => :seagreen,
			"both"     => :deeppink,
		)

		case_linestyles = Dict(
			"clean"    => :solid,
			"overlap"  => :solid,
			"confound" => :solid,
			"both"     => :dash,
		)

		case_linewidths = Dict(
			"clean"    => 2.6,
			"overlap"  => 2.6,
			"confound" => 2.6,
			"both"     => 2.8,
		)

		true_effect_color = :black
		true_effect_linestyle = :dashdot
		true_effect_linewidth = 2.2


		# --------------------------------------------------------
		# Extract four decoding curves
		# --------------------------------------------------------

		t_clean, y_clean =
			_sd_get_curve(
				scores,
				"clean";
				sfreq = cfg.sfreq,
			)

		t_overlap, y_overlap =
			_sd_get_curve(
				scores,
				"overlap";
				sfreq = cfg.sfreq,
			)

		t_confound, y_confound =
			_sd_get_curve(
				scores,
				"confound";
				sfreq = cfg.sfreq,
			)

		t_both, y_both =
			_sd_get_curve(
				scores,
				"both";
				sfreq = cfg.sfreq,
			)


		# --------------------------------------------------------
		# Corresponding ground-truth effect
		# --------------------------------------------------------

		n_timepoints =
			Int(maximum(scores.timepoint))

		t_true,
		y_true,
		true_effect_label =
			_sd_true_effect(
				cfg,
				target,
				n_timepoints;
				target_peak = true_effect_peak,
				use_magnitude = use_effect_magnitude,
			)


		# --------------------------------------------------------
		# Shared x limits
		# --------------------------------------------------------

		x_min = 0.0

		x_max =
			maximum(
				[
					maximum(t_clean),
					maximum(t_overlap),
					maximum(t_confound),
					maximum(t_both),
				]
			)

		x_limits =
			(x_min, x_max)

		x_ticks =
			collect(
				range(
					x_min,
					x_max;
					length = 5,
				)
			)


		# --------------------------------------------------------
		# Shared y limits
		# --------------------------------------------------------

		all_decoding_values =
			Float64.(scores.r)

		all_decoding_values =
			all_decoding_values[
				isfinite.(all_decoding_values)
			]

		isempty(all_decoding_values) &&
			error("No finite decoding scores were found.")

		all_y_values =
			vcat(
				all_decoding_values,
				y_true,
			)

		y_min =
			minimum(all_y_values)

		y_max =
			maximum(all_y_values)

		y_span =
			y_max - y_min

		y_padding =
			y_span == 0 ?
			0.1 :
			0.06 * y_span

		y_limits = (
			y_min - y_padding,
			y_max + y_padding,
		)


		# --------------------------------------------------------
		# Figure and four axes
		# --------------------------------------------------------

		fig = Figure(
			size = (1080, 800),
			backgroundcolor = :white,
		)

		ax_clean = _sd_make_axis(
			fig[1, 1];

			title = "Clean",
			ylabel = ylabel_text,

			show_x = false,
			show_y = true,

			x_limits = x_limits,
			y_limits = y_limits,
			x_ticks = x_ticks,
		)

		ax_overlap = _sd_make_axis(
			fig[1, 2];

			title = "Overlap",
			ylabel = ylabel_text,

			show_x = false,
			show_y = false,

			x_limits = x_limits,
			y_limits = y_limits,
			x_ticks = x_ticks,
		)

		ax_confound_true = _sd_make_axis(
			fig[2, 1];

			title = bottom_left_title,
			ylabel = ylabel_text,

			show_x = true,
			show_y = true,

			x_limits = x_limits,
			y_limits = y_limits,
			x_ticks = x_ticks,
		)

		ax_biased_cases = _sd_make_axis(
			fig[2, 2];

			title = "Overlap, Confound, and Both",
			ylabel = ylabel_text,

			show_x = true,
			show_y = false,

			x_limits = x_limits,
			y_limits = y_limits,
			x_ticks = x_ticks,
		)


		# --------------------------------------------------------
		# Top-left: clean only
		# --------------------------------------------------------

		lines!(
			ax_clean,
			t_clean,
			y_clean;

			color = case_colors["clean"],
			linestyle = case_linestyles["clean"],
			linewidth = case_linewidths["clean"],
		)


		# --------------------------------------------------------
		# Top-right: overlap only
		# --------------------------------------------------------

		lines!(
			ax_overlap,
			t_overlap,
			y_overlap;

			color = case_colors["overlap"],
			linestyle = case_linestyles["overlap"],
			linewidth = case_linewidths["overlap"],
		)


		# --------------------------------------------------------
		# Bottom-left:
		# confound decoding + corresponding true effect
		#
		# No `both` curve here.
		# --------------------------------------------------------

		lines!(
			ax_confound_true,
			t_confound,
			y_confound;

			color = case_colors["confound"],
			linestyle = case_linestyles["confound"],
			linewidth = case_linewidths["confound"],
		)

		lines!(
			ax_confound_true,
			t_true,
			y_true;

			color = true_effect_color,
			linestyle = true_effect_linestyle,
			linewidth = true_effect_linewidth,
		)


		# --------------------------------------------------------
		# Bottom-right:
		# overlap + confound + both
		# --------------------------------------------------------

		lines!(
			ax_biased_cases,
			t_overlap,
			y_overlap;

			color = case_colors["overlap"],
			linestyle = case_linestyles["overlap"],
			linewidth = case_linewidths["overlap"],
		)

		lines!(
			ax_biased_cases,
			t_confound,
			y_confound;

			color = case_colors["confound"],
			linestyle = case_linestyles["confound"],
			linewidth = case_linewidths["confound"],
		)

		lines!(
			ax_biased_cases,
			t_both,
			y_both;

			color = case_colors["both"],
			linestyle = case_linestyles["both"],
			linewidth = case_linewidths["both"],
		)


		# --------------------------------------------------------
		# Shared axes
		# --------------------------------------------------------

		linkxaxes!(
			ax_clean,
			ax_overlap,
			ax_confound_true,
			ax_biased_cases,
		)

		linkyaxes!(
			ax_clean,
			ax_overlap,
			ax_confound_true,
			ax_biased_cases,
		)


		# --------------------------------------------------------
		# Main title
		# --------------------------------------------------------

		Label(
			fig[0, 1:2],
			figure_title;

			fontsize = 25,
			font = :bold,
			padding = (0, 0, 8, 8),
		)


		# --------------------------------------------------------
		# One shared legend
		# --------------------------------------------------------

		legend_elements = LineElement[
			LineElement(
				color = case_colors["clean"],
				linestyle = case_linestyles["clean"],
				linewidth = case_linewidths["clean"],
			),

			LineElement(
				color = case_colors["overlap"],
				linestyle = case_linestyles["overlap"],
				linewidth = case_linewidths["overlap"],
			),

			LineElement(
				color = case_colors["confound"],
				linestyle = case_linestyles["confound"],
				linewidth = case_linewidths["confound"],
			),

			LineElement(
				color = case_colors["both"],
				linestyle = case_linestyles["both"],
				linewidth = case_linewidths["both"],
			),

			LineElement(
				color = true_effect_color,
				linestyle = true_effect_linestyle,
				linewidth = true_effect_linewidth,
			),
		]

		legend_labels = [
			"Clean",
			"Overlap",
			"Confound",
			"Both",
			true_effect_label,
		]

		Legend(
			fig[3, 1:2],
			legend_elements,
			legend_labels;

			orientation = :horizontal,
			framevisible = false,

			labelsize = 13,
			patchsize = (30, 12),

			tellheight = true,
		)


		# --------------------------------------------------------
		# Layout spacing
		# --------------------------------------------------------

		colgap!(fig.layout, 28)
		rowgap!(fig.layout, 18)

		rowsize!(
			fig.layout,
			0,
			Auto(0.10),
		)

		rowsize!(
			fig.layout,
			3,
			Auto(0.10),
		)

		return fig
	end
end

# ╔═╡ c42e65fe-d053-4c0c-b13a-5a222bbd31b8


# ╔═╡ fb06e408-5790-4e80-a0b0-1c21eab7e057


# ╔═╡ 7a6db845-47aa-40ca-9d5d-2cbf802cf544


# ╔═╡ 190a86aa-0bf3-4f41-a95d-343db2e3a718


# ╔═╡ e088c8cc-1bdd-45bc-bac0-eeecf2477cac


# ╔═╡ 5a5c1d01-74c9-422e-9511-762536fd212f


# ╔═╡ 9048dc35-3bec-4e19-87c1-ca163f38f0a5


# ╔═╡ 474bb3b3-64d8-4f94-8f8c-e6d9f3723246
"""rERP decoding"""

# ╔═╡ 17caab34-da12-453b-87d1-679196ac5c9e
begin
	ri= range(
		ridgeModel,
		:lambda;
		lower = 1e-6,
		upper = 1e2,
		scale = :log,
	)
	ridgeTunedModel = TunedModel(
		model=ridgeModel,
		resampling=CV(nfolds=3),
		range = ri,
		tuning=Grid(resolution=4), # test 4 λ values on a logarithmic scale
		measure = RSquared(),
	)
end

# ╔═╡ 4593f063-419f-4898-a202-ba6ddf320656
function fit_rerp_cases(
    cases,
    design,
	model;
	target::Symbol = :continuous,
	nfolds::Int = 3,
)
    scores_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for (case_name, (dat_cont, evts)) in pairs(cases)

        uf_rerp = Unfold.fit(
            UnfoldDecodingModel,
            design,
            evts,
            dat_cont,
			model,
			"stimulus" => target;
			nfolds = nfolds,
			predict_type = Continuous,
            eventcolumn = :event,
			multithreading = false,
        )

		score_case = coeftable(uf_rerp; measure = Statistics.cor, averaged = true,)

		score_case[!, :case] = fill(String(case_name), nrow(score_case))

        push!(scores_tables, score_case)
        fitted_models[case_name] = uf_rerp
    end

    return vcat(scores_tables...), fitted_models
end

# ╔═╡ 660f5fae-737f-401d-9929-a8d1f761de2b


# ╔═╡ 848c8f62-09e2-49b2-bf19-c737725e9127


# ╔═╡ 73e2179d-5369-4745-b7b4-7bc4c6408130


# ╔═╡ 003fd384-9c6e-4f98-88a2-77132225ec3c
"""b2b only"""

# ╔═╡ 30cb6203-d9a4-40a6-b802-f3f00fa3d76a
cross_val_reps = 3

# ╔═╡ 62a8848c-6ae8-41a8-923a-05265d9b1bbd
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

# ╔═╡ 49778986-c23e-4905-9991-aef0e1e76322
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

# ╔═╡ 30324158-c36d-488f-af8b-225d6ce9d7d2
# =================================
# Four cases for standard decoding
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


# ╔═╡ 161394db-5dfb-4fcb-bf62-6326f5ca06b6
# ╠═╡ disabled = true
#=╠═╡
standard_condition_scores, standard_condition_yhats = 
	decode_standard_cases(
		standard_cases,
		ridgeTunedModel;
		target = :condition_num,
		nfolds = 2,
		seed = 12,
		sfreq = cfg.sfreq,
	)
  ╠═╡ =#

# ╔═╡ 5716eddb-27fb-43de-bb5c-e69f99c3a511
# ╠═╡ disabled = true
#=╠═╡
standard_continuous_scores, standard_continuous_yhats = 
	decode_standard_cases(
		standard_cases,
		ridgeTunedModel;
		target = :continuous,
		nfolds = 2,
		seed = 12,
		sfreq = cfg.sfreq,
	)
  ╠═╡ =#

# ╔═╡ 9b09f723-7140-4dc6-adf0-359cbe3b0732
# ╠═╡ disabled = true
#=╠═╡
begin
	fig_condition = plot_standard_decoding_grid(
	standard_condition_scores,
	cfg;
	target = :condition,
)

fig_condition
end 
  ╠═╡ =#

# ╔═╡ c531111b-ed31-4389-a868-03fe28f96e50
# ╠═╡ disabled = true
#=╠═╡
begin
	fig_continuous = plot_standard_decoding_grid(
	standard_continuous_scores,
	cfg;
	target = :continuous,
)

fig_continuous
end 
  ╠═╡ =#

# ╔═╡ a7c4d66b-6e63-4a68-b344-10f19f8ab823
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

# ╔═╡ d24b5c10-9243-4787-99ca-9240a9f215ab
# ================================
# rERP decoding continuous cases
# ================================

begin
	rerp_cases = (
		clean = simulate_case_continuous(
			cfg;
			overlap = false,
			confounded = false,
			seed = 12,
		),

		overlap = simulate_case_continuous(
			cfg;
			overlap = true,
			confounded = false,
			seed = 12,
		),

		confound = simulate_case_continuous(
			cfg;
			overlap = false,
			confounded = true,
			seed = 12,
		),

		both = simulate_case_continuous(
			cfg;
			overlap = true,
			confounded = true,
			seed = 12,
		),
	)
end

# ╔═╡ 9fdd09b1-a4c9-4f2b-a5e5-0ecbce5a5e54
# ╠═╡ disabled = true
#=╠═╡
rerp_continuous_scores, rerp_continuous_models =
    fit_rerp_cases(
        rerp_cases,
        des_rerp,
        ridgeTunedModel;
        target = :continuous,
    )
  ╠═╡ =#

# ╔═╡ 0cdb3a71-0706-4af8-b892-579dbab4526f
# ╠═╡ disabled = true
#=╠═╡
combine(
    groupby(rerp_continuous_scores, :case),
    :estimate => minimum => :minimum_R2,
    :estimate => maximum => :maximum_R2,
    nrow => :nrows,
)
  ╠═╡ =#

# ╔═╡ 63f6376c-cb17-4112-b6b8-c9ac82771178
# ╠═╡ disabled = true
#=╠═╡
rerp_condition_scores, rerp_condition_models =
    fit_rerp_cases(
        rerp_cases,
        des_rerp,
        ridgeTunedModel;
        target = :condition_num,
    )
  ╠═╡ =#

# ╔═╡ 7e8b486f-f14a-4c38-b005-901d2c581928
# ╠═╡ disabled = true
#=╠═╡
fig_rerp_condition = plot_decoding_grid(
    rerp_condition_scores,
    cfg;
    target = :condition_num,
    score_col = :estimate,
)
  ╠═╡ =#

# ╔═╡ f2ad34ad-3408-4b7e-87a6-1be4871084f9
# ╠═╡ disabled = true
#=╠═╡
begin
	fig_rerp_continuous = plot_decoding_grid(
	    rerp_continuous_scores,
	    cfg;
	    target = :continuous,
	    score_col = :estimate,
	)
	save("fig_rerp_continuous.svg", fig_rerp_continuous)
end
  ╠═╡ =#

# ╔═╡ e6941d15-90d8-4d28-98f5-be90b3d91efd
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

# ╔═╡ 774da26a-ea96-4780-b8c8-6fbb96644c98
begin
	using Serialization
	
	mkpath("results")
	
	serialize(
	    "results/one_step_input.jls",
	    (
	        rerp_cases = rerp_cases,
	        sfreq = cfg.sfreq,
	        tau = [-0.1, 1.0],
	        cross_val_reps = cross_val_reps,
	    ),
	)
end

# ╔═╡ 111b0a7c-50b9-435b-a82d-f69cc032ba66
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

# ╔═╡ 62606fef-1c93-475e-a643-3477528222d6


# ╔═╡ 117e5fed-526d-4880-9a0f-ef3fafb81d20


# ╔═╡ 766a7af2-0899-41fb-962d-8d48e23a7713
"""single trial overlap corrected + b2b"""

# ╔═╡ 4c6ba9b2-c8d8-4d35-ba7e-5d07cf059729
# ╠═╡ disabled = true
#=╠═╡
function plot_b2b_grid_pretty_window(
    scores,
    cfg;
    target::Symbol = :condition,
    title = nothing,
    show_aux_window::Bool = true,
    use_abs::Bool = true,
)

    # --------------------------------------------------
    # 1) choose coefficient rows
    # --------------------------------------------------
    coefnames = string.(scores.coefname)

    target_mask =
        if target == :continuous
            occursin.("continuous", coefnames)
        elseif target in (:condition, :condition_num)
            occursin.("condition", coefnames)
        else
            error("target must be :condition or :continuous")
        end

    df = scores[target_mask, :]

    # make sure case column is string
    case_strings = string.(df.case)

    case_order = ["clean", "overlap", "confound", "both"]

    panel_titles = Dict(
        "clean" => "Clean",
        "overlap" => "Overlap",
        "confound" => "Confound",
        "both" => "Both",
    )

    case_colors = Dict(
        "clean" => :dodgerblue,
        "overlap" => :darkorange,
        "confound" => :seagreen,
        "both" => :deeppink,
    )

    # --------------------------------------------------
    # 2) timing windows
    #    adjust if needed
    # --------------------------------------------------
    n170_window = (0.10, 0.22)
    p300_window = (0.20, 0.42)

    if target in (:condition, :condition_num)
        fig_title = isnothing(title) ?
            "Condition predictor recovery (Two-step rERP + B2B)" :
            title

        ylims_use = (0.0, 0.9)

        main_window = n170_window
        aux_window  = p300_window

        main_label = "N170 timing window"
        aux_label  = "P300 timing window"

    else
        fig_title = isnothing(title) ?
            "Continuous predictor recovery (Two-step rERP + B2B)" :
            title

        ylims_use = (0.0, 0.18)

        main_window = p300_window
        aux_window  = n170_window

        main_label = "P300 timing window"
        aux_label  = "N170 timing window"
    end

    # --------------------------------------------------
    # 3) figure
    # --------------------------------------------------
    fig = Figure(size = (1150, 820))
    Label(fig[0, 1:2], fig_title, fontsize = 28, font = :bold)

    legend_handles = Any[]
    legend_labels = String[]

    # --------------------------------------------------
    # 4) panels
    # --------------------------------------------------
    for (i, case_name) in enumerate(case_order)
        row = i <= 2 ? 1 : 2
        col = i <= 2 ? i : i - 2

        ax = Axis(
            fig[row, col],
            title = panel_titles[case_name],
            xlabel = row == 2 ? "Time [s]" : "",
            ylabel = col == 1 ? "B2B recoverability estimate" : "",
            topspinevisible = false,
            rightspinevisible = false,
        )

        subdf = df[case_strings .== case_name, :]

        # if empty, still show panel and continue
        if nrow(subdf) == 0
            @warn "No rows found for case = $case_name and target = $target"
            xlims!(ax, -0.1, 1.0)
            ylims!(ax, ylims_use...)
            continue
        end

        # sort by time
        time = Float64.(subdf.time)
        values = Float64.(subdf.estimate)

        if use_abs
            values = abs.(values)
        end

        ord = sortperm(time)
        time = time[ord]
        values = values[ord]

        y0, y1 = ylims_use

        # ----------------------------------------------
        # main timing window (slightly darker)
        # ----------------------------------------------
        poly!(
            ax,
            Point2f[
                (main_window[1], y0),
                (main_window[2], y0),
                (main_window[2], y1),
                (main_window[1], y1),
            ],
            color = (:black, 0.08),
            strokecolor = :transparent,
        )

        # ----------------------------------------------
        # auxiliary timing window (lighter)
        # ----------------------------------------------
        if show_aux_window
            poly!(
                ax,
                Point2f[
                    (aux_window[1], y0),
                    (aux_window[2], y0),
                    (aux_window[2], y1),
                    (aux_window[1], y1),
                ],
                color = (:gray, 0.05),
                strokecolor = :transparent,
            )
        end

        # boundary lines
        vlines!(
            ax,
            [main_window[1], main_window[2]],
            color = (:black, 0.35),
            linestyle = :dash,
            linewidth = 1.5,
        )

        if show_aux_window
            vlines!(
                ax,
                [aux_window[1], aux_window[2]],
                color = (:gray, 0.35),
                linestyle = :dot,
                linewidth = 1.2,
            )
        end

        # baseline
        hlines!(ax, [0.0], color = (:gray, 0.35), linewidth = 1)

        # main B2B curve
        lines!(
            ax,
            time,
            values,
            color = case_colors[case_name],
            linewidth = 3,
        )

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, ylims_use...)

        # collect legend items once
        if i == 1
            push!(legend_handles, LineElement(color = :dodgerblue, linewidth = 3))
            push!(legend_labels, "Clean")

            push!(legend_handles, LineElement(color = :darkorange, linewidth = 3))
            push!(legend_labels, "Overlap")

            push!(legend_handles, LineElement(color = :seagreen, linewidth = 3))
            push!(legend_labels, "Confound")

            push!(legend_handles, LineElement(color = :deeppink, linewidth = 3))
            push!(legend_labels, "Both")

            push!(legend_handles, PolyElement(color = (:black, 0.08)))
            push!(legend_labels, main_label)

            if show_aux_window
                push!(legend_handles, PolyElement(color = (:gray, 0.05)))
                push!(legend_labels, aux_label)
            end
        end
    end

    # --------------------------------------------------
    # 5) legend
    # --------------------------------------------------
    Legend(
        fig[3, 1:2],
        legend_handles,
        legend_labels,
        orientation = :horizontal,
        framevisible = false,
    )

    rowsize!(fig.layout, 3, Auto(0.12))

    return fig
end
  ╠═╡ =#

# ╔═╡ 9f8610a5-621e-406a-84f4-b3631e62af85
function UnfoldDecode.solver_b2b(
    X,
    data::AbstractArray{T,2};
    cross_val_reps = 1,
    multithreading = true,
    show_progress = true,
    solver_G = UnfoldDecode.model_ridge,
    solver_H = UnfoldDecode.model_ridge,
) where {T<:Number}

    E = zeros(
        size(X, 2),
        size(X, 2),
    )

    prog = UnfoldDecode.Progress(
        cross_val_reps;
        dt = 0.1,
        enabled = show_progress,
    )

    UnfoldDecode.@maybe_threads multithreading for m = 1:cross_val_reps

        k_ix = collect(
            UnfoldDecode.Kfold(
                size(data, 2),
                2,
            )
        )

        Y1 = data[:, k_ix[1]]
        Y2 = data[:, k_ix[2]]

        X1 = X[k_ix[1], :]
        X2 = X[k_ix[2], :]

        # backward model
        Y2G = solver_G(
            Y1',
            X1,
            Y2',
        )

        # forward model
        H = solver_H(
            X2,
            Y2G,
        )

        E[:, :] =
            E[:, :] +
            UnfoldDecode.Diagonal(
                H[UnfoldDecode.diagind(H)]
            )

        UnfoldDecode.ProgressMeter.next!(
            prog;
            showvalues = [
                (:cross_val_rep, m)
            ],
        )
    end

    E[:, :] ./= cross_val_reps

    W =
        (X * E[:, :])' /
        data[:, :]

    # B2B S = diag(E)
    beta = UnfoldDecode.diag(E)

    # Unfold expects coef dimension
    beta = reshape(beta, 1, :)

    modelinfo = Dict(
        "W" => W,
        "E" => E,
        "cross_val_reps" => cross_val_reps,
    )

    return UnfoldDecode.Unfold.LinearModelFit{
        eltype(beta),
        2,
    }(
        beta,
        modelinfo,
    )
end

# ╔═╡ 5dd04bb5-edf0-4cad-8460-8fb4afd6256b
b2b_solver = (x, y) -> UnfoldDecode.solver_b2b(x, y; cross_val_reps=cross_val_reps);

# ╔═╡ b2acf9d6-3715-4a67-9f60-c50b92996d47
#=╠═╡
b2b_only_scores, b2b_only_models =
    fit_b2b_only_cases(
        standard_cases,
        des_b2b_plain,
    )
  ╠═╡ =#

# ╔═╡ fa1c36ea-14c9-40ff-bff3-34f6de835b1d
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ d52118da-353c-419c-8f02-939725d21570
#=╠═╡
fig_b2b_continuous = plot_b2b_grid(
    b2b_only_scores,
    cfg;
    target = :continuous,
)
  ╠═╡ =#

# ╔═╡ 810aeb3f-ab03-485a-9faf-c4a76da3ad20
#=╠═╡
fig_b2b_condition = plot_b2b_grid(
    b2b_only_scores,
    cfg;
    target = :condition,
)
  ╠═╡ =#

# ╔═╡ eb646ead-ca2b-4b60-94f4-1d341bba4ae8
function fit_two_step_b2b_cases(
    cases,
    des_rerp;
    eventname = "stimulus",
    eventcolumn = :event,
)
    result_tables = DataFrame[]

    for (case_name, (dat_cont, evts)) in pairs(cases)

		# ==========================================
		# Step 1: rERP
		# ==========================================
		
		uf_rerp = Unfold.fit(
		    UnfoldLinearModelContinuousTime,
		    des_rerp,
		    evts,
		    dat_cont;
		    eventcolumn = eventcolumn,
		)
		
		
		# ==========================================
		# Step 2: corrected single trials
		# ==========================================
		
		dat_corrected = UnfoldDecode.singletrials(
		    dat_cont,
		    uf_rerp,
		    evts,
		    eventname,
		    eventcolumn,
		)
		
		
		# ==========================================
		# Step 2.5: remove incomplete boundary epochs
		# ==========================================
		
		has_missing = dropdims(
		    any(
		        ismissing.(dat_corrected),
		        dims = (1, 2),
		    ),
		    dims = (1, 2),
		)
		
		good_trials = .!has_missing
		
		dat_corrected =
		    Float64.(
		        dat_corrected[:, :, good_trials]
		    )
		
		evts_corrected =
		    evts[good_trials, :]
		
		@assert size(dat_corrected, 3) == nrow(evts_corrected)
		@assert count(ismissing, dat_corrected) == 0
		
		
		# ==========================================
		# Step 3: epoched B2B design
		# ==========================================
				times_corrected =
		    Unfold.times(uf_rerp)[1]
		
		@assert length(times_corrected) ==
		        size(dat_corrected, 2)
		
		des_b2b_corrected = [
		    Any => (
		        @formula(0 ~ 1 + condition + continuous),
		        times_corrected,
		    )
		]
		
		
		# ==========================================
		# Step 4: B2B
		# ==========================================
		
		uf_b2b = Unfold.fit(
    UnfoldDecode.UnfoldModel,
    des_b2b_corrected,
    evts_corrected,
    dat_corrected;
    solver = b2b_solver,
)

        # ------------------------------------------
        # 6. collect coefficients
        # ------------------------------------------

        result = coeftable(uf_b2b)

        result = result[
            result.coefname .!= "(Intercept)",
            :
        ]

        result[!, :case] =
            fill(String(case_name), nrow(result))

        push!(result_tables, result)
    end

    return vcat(result_tables...)
end

# ╔═╡ 107be910-a2e4-4c86-9c39-7a6bfee91a14
# ╠═╡ disabled = true
#=╠═╡
two_step_scores =
    fit_two_step_b2b_cases(
        rerp_cases,
        des_rerp,
    )
  ╠═╡ =#

# ╔═╡ ef5d1d63-3110-43b1-b9a4-649f4786bea7
# ╠═╡ disabled = true
#=╠═╡
fig_two_step_continuous = plot_b2b_grid(
    two_step_scores,
    cfg;
    target = :continuous,
)
  ╠═╡ =#

# ╔═╡ e9d86b5d-c4ff-49bb-a551-3cb7e963f576
# ╠═╡ disabled = true
#=╠═╡
fig_two_step_condition = plot_b2b_grid(
    two_step_scores,
    cfg;
    target = :condition,
)
  ╠═╡ =#

# ╔═╡ 1d39e8ad-5509-46e2-90a7-c72540337fd5
function fit_one_step_b2b_cases(
    cases,
    design,
)
    result_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for (case_name, (dat_cont, evts)) in pairs(cases)

        uf_b2b = Unfold.fit(
            UnfoldDecode.UnfoldModel,
            design,
            evts,
            dat_cont;
            solver = b2b_solver,
        )

        result = coeftable(uf_b2b)

        result = result[
            result.coefname .!= "(Intercept)",
            :
        ]

        result[!, :case] =
            fill(String(case_name), nrow(result))

        push!(result_tables, result)
        fitted_models[case_name] = uf_b2b
    end

    return vcat(result_tables...), fitted_models
end

# ╔═╡ 61f2f8a2-4c14-4914-b5d5-b8d3d23e73b7
pwd()

# ╔═╡ 4531fbfb-24b0-46c0-a16f-46b203b24e0e
# ╠═╡ disabled = true
#=╠═╡
one_step_scores,
one_step_models =
    fit_one_step_b2b_cases(
        rerp_cases,
        des_rerp,
    )
  ╠═╡ =#

# ╔═╡ 5e8adcf9-b8de-42fa-9aca-ac802a10880d
# ╠═╡ disabled = true
#=╠═╡
fig_one_step_continuous = plot_b2b_grid(
    one_step_scores,
    cfg;
    target = :continuous,
)
  ╠═╡ =#

# ╔═╡ 1eaf9e64-7bbc-41cb-975b-cb7c5d8239d9
#=╠═╡
fig_one_step_condition = plot_b2b_grid(
    one_step_scores,
    cfg;
    target = :condition,
)
  ╠═╡ =#

# ╔═╡ 70130678-72e7-437e-8580-5245ee9d57bb
cfg_debug = Decoding_Config(
    n_trials = 200,
    sfreq = 50.0,
    n_channels = 5,
    noiselevel = cfg.noiselevel,
    channel_noise_sd = cfg.channel_noise_sd,
    β0_n170 = cfg.β0_n170,
    β_condition = cfg.β_condition,
    β0_p300 = cfg.β0_p300,
    β_continuous = cfg.β_continuous,
    rho = cfg.rho,
    onset_interval_ms = cfg.onset_interval_ms,
    onset_condition_bias = cfg.onset_condition_bias,
)

# ╔═╡ d0d1859d-304e-46f6-8d20-d7f685418082


# ╔═╡ 7cddeb26-ea36-45db-aa69-364f50148c1b
# ╠═╡ disabled = true
#=╠═╡
function plot_decoding_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    score_col::Symbol,
    true_effect_peak::Real = 0.45,
)
    columns = propertynames(scores)

    :case in columns ||
        error("`scores` must contain a :case column.")

    score_col in columns ||
        error(
            "`scores` has no column $score_col. " *
            "Available columns: $columns"
        )

    (:time in columns || :timepoint in columns) ||
        error("`scores` must contain either :time or :timepoint.")

    # -------------------------------------------------
    # Extract one case
    # -------------------------------------------------

    function get_curve(case_name)
        mask = string.(scores.case) .== case_name

        any(mask) ||
            error(
                "No case named \"$case_name\". " *
                "Available cases: $(unique(string.(scores.case)))"
            )

        time =
            if :time in columns
                Float64.(scores[mask, :time])
            else
                (
                    Float64.(scores[mask, :timepoint]) .- 1
                ) ./ cfg.sfreq
            end

        values = Float64.(scores[mask, score_col])

        order = sortperm(time)

        return time[order], values[order]
    end

    curves = Dict(
        name => get_curve(name)
        for name in ("clean", "overlap", "confound", "both")
    )

    # -------------------------------------------------
    # Ground-truth temporal shape
    # -------------------------------------------------

    raw_effect, true_label =
        if target in (:condition, :condition_num)
            (
                cfg.β_condition .*
                UnfoldSim.n170(; sfreq = cfg.sfreq),

                "True N170 condition effect (scaled)",
            )

        elseif target == :continuous
            (
                cfg.β_continuous .*
                UnfoldSim.p300(; sfreq = cfg.sfreq),

                "True P300 continuous effect (scaled)",
            )

        else
            error(
                "target must be :condition_num or :continuous."
            )
        end

    raw_effect = abs.(Float64.(raw_effect))

    true_values =
        maximum(raw_effect) == 0 ?
        zeros(length(raw_effect)) :
        true_effect_peak .* raw_effect ./ maximum(raw_effect)

    true_time =
        (0:length(true_values)-1) ./ cfg.sfreq

    # -------------------------------------------------
    # Labels
    # -------------------------------------------------

    condition_target =
        target in (:condition, :condition_num)

    figure_title =
        condition_target ?
        "Condition decoding" :
        "Continuous decoding"

    ylabel_text =
        if score_col == :r
            "Cross-validated correlation r"
        elseif score_col in (:r_squared, :estimate)
            "Cross-validated R²"
        else
            string(score_col)
        end

    confound_title =
        condition_target ?
        "Confound and true condition effect" :
        "Confound and true continuous effect"

    colors = Dict(
        "clean" => :dodgerblue,
        "overlap" => :darkorange,
        "confound" => :seagreen,
        "both" => :deeppink,
    )

    # -------------------------------------------------
    # Figure
    # -------------------------------------------------

    fig = Figure(
        size = (1000, 720),
        backgroundcolor = :white,
    )

    function make_axis(
        position,
        title;
        show_x::Bool,
        show_y::Bool,
    )
        ax = Axis(
            position;
            title = title,
            xlabel = show_x ? "Time [s]" : "",
            ylabel = show_y ? ylabel_text : "",
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
            xticklabelsvisible = show_x,
            xticksvisible = show_x,
            yticklabelsvisible = show_y,
            yticksvisible = show_y,
        )

        hlines!(
            ax,
            [0.0];
            color = (:gray, 0.5),
            linestyle = :dash,
        )

        return ax
    end

    ax_clean = make_axis(
        fig[1, 1],
        "Clean";
        show_x = false,
        show_y = true,
    )

    ax_overlap = make_axis(
        fig[1, 2],
        "Overlap";
        show_x = false,
        show_y = false,
    )

    ax_confound = make_axis(
        fig[2, 1],
        confound_title;
        show_x = true,
        show_y = true,
    )

    ax_both = make_axis(
        fig[2, 2],
        "Overlap, Confound, and Both";
        show_x = true,
        show_y = false,
    )

    # Clean
    lines!(
        ax_clean,
        curves["clean"]...;
        color = colors["clean"],
        linewidth = 2.5,
    )

    # Overlap
    lines!(
        ax_overlap,
        curves["overlap"]...;
        color = colors["overlap"],
        linewidth = 2.5,
    )

    # Confound + truth
    lines!(
        ax_confound,
        curves["confound"]...;
        color = colors["confound"],
        linewidth = 2.5,
    )

    lines!(
        ax_confound,
        true_time,
        true_values;
        color = :black,
        linestyle = :dashdot,
        linewidth = 2,
    )

    # Biased cases
    lines!(
        ax_both,
        curves["overlap"]...;
        color = colors["overlap"],
        linewidth = 2.5,
    )

    lines!(
        ax_both,
        curves["confound"]...;
        color = colors["confound"],
        linewidth = 2.5,
    )

    lines!(
        ax_both,
        curves["both"]...;
        color = colors["both"],
        linestyle = :dash,
        linewidth = 2.5,
    )

    linkxaxes!(
        ax_clean,
        ax_overlap,
        ax_confound,
        ax_both,
    )

    linkyaxes!(
        ax_clean,
        ax_overlap,
        ax_confound,
        ax_both,
    )

    Label(
        fig[0, 1:2],
        figure_title;
        fontsize = 24,
        font = :bold,
    )

    legend_elements = [
        LineElement(
            color = colors["clean"],
            linewidth = 2.5,
        ),
        LineElement(
            color = colors["overlap"],
            linewidth = 2.5,
        ),
        LineElement(
            color = colors["confound"],
            linewidth = 2.5,
        ),
        LineElement(
            color = colors["both"],
            linestyle = :dash,
            linewidth = 2.5,
        ),
        LineElement(
            color = :black,
            linestyle = :dashdot,
            linewidth = 2,
        ),
    ]

    Legend(
        fig[3, 1:2],
        legend_elements,
        [
            "Clean",
            "Overlap",
            "Confound",
            "Both",
            true_label,
        ];
        orientation = :horizontal,
        framevisible = false,
    )

    colgap!(fig.layout, 24)
    rowgap!(fig.layout, 16)

    return fig
end
  ╠═╡ =#

# ╔═╡ ee0f9f8f-ecdb-475d-b601-0d600bda4a87
# ╠═╡ disabled = true
#=╠═╡
function fit_b2b_only_cases(
    cases,
	design;
)
    scores_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for (case_name, (dat, evts)) in pairs(cases)

        uf_b2b_only = Unfold.fit(
            UnfoldDecode.UnfoldModel,
            design,
            evts,
            dat;
			solver = b2b_solver,
        )

		b2b_only_results = coeftable(uf_b2b_only)
		b2b_only_results.estimate = abs.(b2b_only_results.estimate)
		b2b_only_results = b2b_only_results[b2b_only_results.coefname .!="(Intercept)", :]
		score_case = b2b_only_results

		score_case[!, :case] = fill(String(case_name), nrow(score_case))

        push!(scores_tables, score_case)
        fitted_models[case_name] = uf_b2b_only
    end

    return vcat(scores_tables...), fitted_models
end
  ╠═╡ =#

# ╔═╡ ae25b8b8-8c27-47d6-ae91-0bc5425d80e0
#=╠═╡
function plot_decoding_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    score_col::Symbol = :estimate,
    guide_peak::Real = 0.45,
)
    columns = propertynames(scores)

    :case in columns ||
        error("`scores` must contain a :case column.")

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

    # ==================================================
    # Extract one decoding curve
    # ==================================================

    function get_curve(case_name::String)
        mask = string.(scores.case) .== case_name

        any(mask) ||
            error(
                "No case named \"$case_name\". " *
                "Available cases: $(unique(string.(scores.case)))"
            )

        time =
            if time_col == :time
                Float64.(scores[mask, :time])
            else
                (
                    Float64.(scores[mask, :timepoint]) .- 1
                ) ./ cfg.sfreq
            end

        values =
            Float64.(scores[mask, score_col])

        order = sortperm(time)

        return time[order], values[order]
    end

    # ==================================================
    # Ground-truth timing guides
    #
    # These curves only indicate expected timing.
    # They are not theoretical R² values.
    # ==================================================

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

    n170_guide =
        scale_guide(n170_effect)

    p300_guide =
        scale_guide(p300_effect)

    n170_time =
        (0:length(n170_guide)-1) ./ cfg.sfreq

    p300_time =
        (0:length(p300_guide)-1) ./ cfg.sfreq

    condition_target =
        target in (:condition, :condition_num)

    if condition_target
        figure_title = "Condition decoding"
    elseif target == :continuous
        figure_title = "Continuous decoding"
    else
        error(
            "`target` must be :condition, " *
            ":condition_num, or :continuous."
        )
    end

    # ==================================================
    # Panel settings
    # ==================================================

    case_order = (
        "clean",
        "overlap",
        "confound",
        "both",
    )

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

    ylabel_text =
        if score_col == :r
            "Cross-validated correlation r"
        else
            string(score_col)
        end

    # ==================================================
    # Figure
    # ==================================================

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

        # ----------------------------------------------
        # Decoding score
        # ----------------------------------------------

        time, values =
            get_curve(case_name)

        lines!(
            ax,
            time,
            values;
            color = colors[case_name],
            linewidth = 2.5,
        )

        # ----------------------------------------------
        # Decide which timing guides belong in this panel
        # ----------------------------------------------

        confounded_case =
            case_name in ("confound", "both")

        show_n170 =
            if condition_target
                # Direct condition effect in every case
                true
            else
                # N170 only becomes a proxy for continuous
                # when condition and continuous are correlated
                confounded_case
            end

        show_p300 =
            if condition_target
                # P300 only becomes a proxy for condition
                # when condition and continuous are correlated
                confounded_case
            else
                # Direct continuous effect in every case
                true
            end

        # ----------------------------------------------
        # N170 style is always:
        # black + dashed
        # ----------------------------------------------

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

        # ----------------------------------------------
        # P300 style is always:
        # gray + dotted
        # ----------------------------------------------

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

        # Zero reference line
        hlines!(
            ax,
            [0.0];
            color = (:gray, 0.35),
            linestyle = :solid,
            linewidth = 1,
        )

        push!(axes, ax)
    end

    # Same x- and y-axis ranges across panels
    linkxaxes!(axes...)
    linkyaxes!(axes...)

    Label(
        fig[0, 1:2],
        figure_title;
        fontsize = 24,
        font = :bold,
    )

    # ==================================================
    # Legend
    # ==================================================

    legend_elements = [
        LineElement(
            color = colors["clean"],
            linewidth = 2.5,
        ),
        LineElement(
            color = colors["overlap"],
            linewidth = 2.5,
        ),
        LineElement(
            color = colors["confound"],
            linewidth = 2.5,
        ),
        LineElement(
            color = colors["both"],
            linewidth = 2.5,
        ),
        LineElement(
            color = :black,
            linestyle = :dash,
            linewidth = 2.5,
        ),
        LineElement(
            color = :gray40,
            linestyle = :dot,
            linewidth = 2.5,
        ),
    ]

    Legend(
        fig[3, 1:2],
        legend_elements,
        [
            "Clean",
            "Overlap",
            "Confound",
            "Both",
            "N170 timing guide",
            "P300 timing guide",
        ];
        orientation = :horizontal,
        framevisible = false,
        nbanks = 2,
    )

    colgap!(fig.layout, 24)
    rowgap!(fig.layout, 16)

    return fig
end
  ╠═╡ =#

# ╔═╡ eebded5d-134d-43fe-9107-9ce960396de7
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═f43ecec0-8a81-11f1-a925-3d063fa9b39c
# ╠═c8a40da4-af0c-470a-9227-2073e5bfc622
# ╠═c9daa551-e9be-4ff2-96c0-86d857189650
# ╠═e0c3a358-844c-4ff4-931d-b83a512fc8cf
# ╠═0a43fabc-9209-40da-b568-51f0ff423954
# ╠═07238baa-e2bb-4c65-9f65-4506170dbd83
# ╟─49778986-c23e-4905-9991-aef0e1e76322
# ╟─b937fc2e-42bb-4ad4-83d2-bb29c00ccce0
# ╟─bd65f918-f274-4d2c-a51d-dad55da44133
# ╠═ad6634be-abea-42b0-bd59-c24ffd40106f
# ╠═30324158-c36d-488f-af8b-225d6ce9d7d2
# ╟─ac6cb027-550b-4a93-a2d4-2ec3ad3325ac
# ╟─b0d97bbf-6e1b-496d-8ac7-0dc9688cc309
# ╠═161394db-5dfb-4fcb-bf62-6326f5ca06b6
# ╠═5716eddb-27fb-43de-bb5c-e69f99c3a511
# ╠═9b09f723-7140-4dc6-adf0-359cbe3b0732
# ╠═c531111b-ed31-4389-a868-03fe28f96e50
# ╠═c42e65fe-d053-4c0c-b13a-5a222bbd31b8
# ╠═fb06e408-5790-4e80-a0b0-1c21eab7e057
# ╠═7a6db845-47aa-40ca-9d5d-2cbf802cf544
# ╠═190a86aa-0bf3-4f41-a95d-343db2e3a718
# ╠═e088c8cc-1bdd-45bc-bac0-eeecf2477cac
# ╠═5a5c1d01-74c9-422e-9511-762536fd212f
# ╠═9048dc35-3bec-4e19-87c1-ca163f38f0a5
# ╠═474bb3b3-64d8-4f94-8f8c-e6d9f3723246
# ╠═a7c4d66b-6e63-4a68-b344-10f19f8ab823
# ╠═d24b5c10-9243-4787-99ca-9240a9f215ab
# ╟─17caab34-da12-453b-87d1-679196ac5c9e
# ╠═4593f063-419f-4898-a202-ba6ddf320656
# ╠═ae25b8b8-8c27-47d6-ae91-0bc5425d80e0
# ╠═9fdd09b1-a4c9-4f2b-a5e5-0ecbce5a5e54
# ╠═63f6376c-cb17-4112-b6b8-c9ac82771178
# ╠═7e8b486f-f14a-4c38-b005-901d2c581928
# ╠═f2ad34ad-3408-4b7e-87a6-1be4871084f9
# ╠═0cdb3a71-0706-4af8-b892-579dbab4526f
# ╠═7cddeb26-ea36-45db-aa69-364f50148c1b
# ╠═660f5fae-737f-401d-9929-a8d1f761de2b
# ╠═848c8f62-09e2-49b2-bf19-c737725e9127
# ╠═73e2179d-5369-4745-b7b4-7bc4c6408130
# ╠═003fd384-9c6e-4f98-88a2-77132225ec3c
# ╠═30cb6203-d9a4-40a6-b802-f3f00fa3d76a
# ╠═5dd04bb5-edf0-4cad-8460-8fb4afd6256b
# ╠═e6941d15-90d8-4d28-98f5-be90b3d91efd
# ╠═ee0f9f8f-ecdb-475d-b601-0d600bda4a87
# ╠═eebded5d-134d-43fe-9107-9ce960396de7
# ╠═b2acf9d6-3715-4a67-9f60-c50b92996d47
# ╠═fa1c36ea-14c9-40ff-bff3-34f6de835b1d
# ╠═62a8848c-6ae8-41a8-923a-05265d9b1bbd
# ╠═111b0a7c-50b9-435b-a82d-f69cc032ba66
# ╠═d52118da-353c-419c-8f02-939725d21570
# ╠═810aeb3f-ab03-485a-9faf-c4a76da3ad20
# ╠═62606fef-1c93-475e-a643-3477528222d6
# ╠═117e5fed-526d-4880-9a0f-ef3fafb81d20
# ╠═766a7af2-0899-41fb-962d-8d48e23a7713
# ╠═eb646ead-ca2b-4b60-94f4-1d341bba4ae8
# ╠═107be910-a2e4-4c86-9c39-7a6bfee91a14
# ╠═ef5d1d63-3110-43b1-b9a4-649f4786bea7
# ╠═e9d86b5d-c4ff-49bb-a551-3cb7e963f576
# ╟─4c6ba9b2-c8d8-4d35-ba7e-5d07cf059729
# ╠═1d39e8ad-5509-46e2-90a7-c72540337fd5
# ╠═9f8610a5-621e-406a-84f4-b3631e62af85
# ╠═774da26a-ea96-4780-b8c8-6fbb96644c98
# ╠═61f2f8a2-4c14-4914-b5d5-b8d3d23e73b7
# ╠═4531fbfb-24b0-46c0-a16f-46b203b24e0e
# ╠═5e8adcf9-b8de-42fa-9aca-ac802a10880d
# ╠═1eaf9e64-7bbc-41cb-975b-cb7c5d8239d9
# ╠═70130678-72e7-437e-8580-5245ee9d57bb
# ╠═d0d1859d-304e-46f6-8d20-d7f685418082
