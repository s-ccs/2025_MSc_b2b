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

# ╔═╡ ae25b8b8-8c27-47d6-ae91-0bc5425d80e0
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
standard_condition_scores, standard_condition_yhats = 
	decode_standard_cases(
		standard_cases,
		ridgeTunedModel;
		target = :condition_num,
		nfolds = 2,
		seed = 12,
		sfreq = cfg.sfreq,
	)

# ╔═╡ 5716eddb-27fb-43de-bb5c-e69f99c3a511
standard_continuous_scores, standard_continuous_yhats = 
	decode_standard_cases(
		standard_cases,
		ridgeTunedModel;
		target = :continuous,
		nfolds = 2,
		seed = 12,
		sfreq = cfg.sfreq,
	)

# ╔═╡ 9b09f723-7140-4dc6-adf0-359cbe3b0732
begin
	fig_condition = plot_standard_decoding_grid(
	standard_condition_scores,
	cfg;
	target = :condition,
)

fig_condition
end 

# ╔═╡ c531111b-ed31-4389-a868-03fe28f96e50
begin
	fig_continuous = plot_standard_decoding_grid(
	standard_continuous_scores,
	cfg;
	target = :continuous,
)

fig_continuous
end 

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
rerp_continuous_scores, rerp_continuous_models =
    fit_rerp_cases(
        rerp_cases,
        des_rerp,
        ridgeTunedModel;
        target = :continuous,
        nfolds = 3,
    )

# ╔═╡ 63f6376c-cb17-4112-b6b8-c9ac82771178
rerp_condition_scores, rerp_condition_models =
    fit_rerp_cases(
        rerp_cases,
        des_rerp,
        ridgeTunedModel;
        target = :condition_num,
        nfolds = 3,
    )

# ╔═╡ 7e8b486f-f14a-4c38-b005-901d2c581928
fig_rerp_condition = plot_decoding_grid(
    rerp_condition_scores,
    cfg;
    target = :condition_num,
    score_col = :estimate,
)

# ╔═╡ f2ad34ad-3408-4b7e-87a6-1be4871084f9
fig_rerp_continuous = plot_decoding_grid(
    rerp_continuous_scores,
    cfg;
    target = :continuous,
    score_col = :estimate,
)

# ╔═╡ 0cdb3a71-0706-4af8-b892-579dbab4526f
combine(
    groupby(rerp_continuous_scores, :case),
    :estimate => minimum => :minimum_R2,
    :estimate => maximum => :maximum_R2,
    nrow => :nrows,
)

# ╔═╡ 660f5fae-737f-401d-9929-a8d1f761de2b


# ╔═╡ bbf7c325-9db8-4598-8456-2e050124a94a
# ╠═╡ disabled = true
#=╠═╡
begin

    fig_rerp_continuous = plot_decoding_grid(
        rerp_continuous_scores,
        cfg;
        target = :continuous,
        score_col = :estimate,
    )

    fig_rerp_continuous
end
  ╠═╡ =#

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

# ╔═╡ Cell order:
# ╠═f43ecec0-8a81-11f1-a925-3d063fa9b39c
# ╠═c8a40da4-af0c-470a-9227-2073e5bfc622
# ╠═c9daa551-e9be-4ff2-96c0-86d857189650
# ╠═e0c3a358-844c-4ff4-931d-b83a512fc8cf
# ╠═0a43fabc-9209-40da-b568-51f0ff423954
# ╠═07238baa-e2bb-4c65-9f65-4506170dbd83
# ╠═49778986-c23e-4905-9991-aef0e1e76322
# ╠═b937fc2e-42bb-4ad4-83d2-bb29c00ccce0
# ╠═bd65f918-f274-4d2c-a51d-dad55da44133
# ╠═ad6634be-abea-42b0-bd59-c24ffd40106f
# ╠═30324158-c36d-488f-af8b-225d6ce9d7d2
# ╠═ac6cb027-550b-4a93-a2d4-2ec3ad3325ac
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
# ╠═17caab34-da12-453b-87d1-679196ac5c9e
# ╠═4593f063-419f-4898-a202-ba6ddf320656
# ╟─ae25b8b8-8c27-47d6-ae91-0bc5425d80e0
# ╠═9fdd09b1-a4c9-4f2b-a5e5-0ecbce5a5e54
# ╠═63f6376c-cb17-4112-b6b8-c9ac82771178
# ╟─62a8848c-6ae8-41a8-923a-05265d9b1bbd
# ╠═7e8b486f-f14a-4c38-b005-901d2c581928
# ╠═f2ad34ad-3408-4b7e-87a6-1be4871084f9
# ╠═0cdb3a71-0706-4af8-b892-579dbab4526f
# ╠═660f5fae-737f-401d-9929-a8d1f761de2b
# ╟─bbf7c325-9db8-4598-8456-2e050124a94a
# ╟─7cddeb26-ea36-45db-aa69-364f50148c1b
