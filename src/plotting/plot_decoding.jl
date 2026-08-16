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