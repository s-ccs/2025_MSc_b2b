# ========================================
# Standard decoding plotting
# ========================================


# ============================================================
# Extract one decoding curve
# ============================================================

function _sd_get_curve(
    scores::AbstractDataFrame,
    case_name::AbstractString;
    sfreq::Real
)
    required_columns = Set([:r, :case])
    available_columns = Set(propertynames(scores))

    issubset(required_columns, available_columns) ||
        error(
            "`scores` must contain :r and :case. " *
            "Available columns: $(propertynames(scores))"
        )

    case_mask =
        string.(scores.case) .== case_name

    any(case_mask) ||
        error(
            "No rows found for case = \"$case_name\". " *
            "Available cases: $(unique(string.(scores.case)))"
        )

    # --------------------------------------------------------
    # Prefer real time if available
    # --------------------------------------------------------

    if :time in propertynames(scores)

        sub = DataFrame(
            scores[
                case_mask,
                [:time, :r]
            ]
        )

        sort!(sub, :time)

        return (
            Float64.(sub.time),
            Float64.(sub.r)
        )

    # --------------------------------------------------------
    # Otherwise convert timepoint to seconds
    # --------------------------------------------------------

    elseif :timepoint in propertynames(scores)

        sub = DataFrame(
            scores[
                case_mask,
                [:timepoint, :r]
            ]
        )

        sort!(sub, :timepoint)

        time =
            (Float64.(sub.timepoint) .- 1.0) ./ sfreq

        return (
            time,
            Float64.(sub.r)
        )

    else
        error(
            "`scores` must contain either :time or :timepoint."
        )
    end
end


# ============================================================
# Trim curve to plotting window
# ============================================================

function _sd_trim_curve(
    time::AbstractVector,
    values::AbstractVector,
    x_window
)
    x_min, x_max = x_window

    mask =
        (time .>= x_min) .&
        (time .<= x_max)

    return (
        Float64.(time[mask]),
        Float64.(values[mask])
    )
end


# ============================================================
# Create ground-truth target effect
# ============================================================

function _sd_true_effect(
    cfg,
    target::Symbol;
    target_peak::Real = 0.45,
    use_magnitude::Bool = true
)
    raw_effect, label =
        if target == :condition

            (
                cfg.β_condition .*
                    UnfoldSim.n170(; sfreq = cfg.sfreq),

                "Ground-truth N170 effect waveform (scaled)"
            )

        elseif target == :continuous

            (
                cfg.β_continuous .*
                    UnfoldSim.p300(; sfreq = cfg.sfreq),

                "Ground-truth P300 effect waveform (scaled)"
            )

        else
            error(
                "Unknown target: $target. " *
                "Use :condition or :continuous."
            )
        end

    raw_effect =
        Float64.(raw_effect)

    # Decoding correlation has no ERP polarity,
    # so magnitude is easier to compare visually.
    effect_for_plot =
        use_magnitude ?
            abs.(raw_effect) :
            raw_effect

    maximum_absolute =
        maximum(abs.(effect_for_plot))

    scaled_effect =
        maximum_absolute == 0 ?
            zeros(Float64, length(effect_for_plot)) :
            Float64(target_peak) .*
            effect_for_plot ./ maximum_absolute

    # ERP waveform starts at stimulus onset = 0 s
    time =
        (0:length(scaled_effect)-1) ./ cfg.sfreq

    return (
        Float64.(time),
        Float64.(scaled_effect),
        label
    )
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
    x_ticks
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

        backgroundcolor = :white
    )

    xlims!(ax, x_limits...)
    ylims!(ax, y_limits...)

    # Zero correlation reference line
    hlines!(
        ax,
        [0.0];
        color = (:gray45, 0.55),
        linestyle = :dash,
        linewidth = 1.0
    )

    return ax
end


# ============================================================
# Main plotting function
# ============================================================

function plot_standard_decoding_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    x_window = (-0.1, 0.43),
    true_effect_peak::Real = 0.45,
    use_effect_magnitude::Bool = true,
    show_true_effect::Bool = true,
    true_effect_cases = ("overlap", "confound")
)

    # --------------------------------------------------------
    # Target-specific labels
    # --------------------------------------------------------

    figure_title, ylabel_text =
        if target == :condition

            (
                "Condition decoding",
                "corr(predicted condition, true condition)"
            )

        elseif target == :continuous

            (
                "Continuous decoding",
                "corr(predicted continuous, true continuous)"
            )

        else
            error(
                "Unknown target: $target. " *
                "Use :condition or :continuous."
            )
        end


    # --------------------------------------------------------
    # Visual identity
    # --------------------------------------------------------

    case_colors = Dict(
        "clean"    => :dodgerblue,
        "overlap"  => :darkorange,
        "confound" => :seagreen,
        "both"     => :deeppink
    )

    case_linestyles = Dict(
        "clean"    => :solid,
        "overlap"  => :solid,
        "confound" => :solid,
        "both"     => :dash
    )

    case_linewidths = Dict(
        "clean"    => 2.6,
        "overlap"  => 2.6,
        "confound" => 2.6,
        "both"     => 2.8
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
            sfreq = cfg.sfreq
        )

    t_overlap, y_overlap =
        _sd_get_curve(
            scores,
            "overlap";
            sfreq = cfg.sfreq
        )

    t_confound, y_confound =
        _sd_get_curve(
            scores,
            "confound";
            sfreq = cfg.sfreq
        )

    t_both, y_both =
        _sd_get_curve(
            scores,
            "both";
            sfreq = cfg.sfreq
        )


    # --------------------------------------------------------
    # Trim decoding curves to plotting window
    # --------------------------------------------------------

    t_clean, y_clean =
        _sd_trim_curve(
            t_clean,
            y_clean,
            x_window
        )

    t_overlap, y_overlap =
        _sd_trim_curve(
            t_overlap,
            y_overlap,
            x_window
        )

    t_confound, y_confound =
        _sd_trim_curve(
            t_confound,
            y_confound,
            x_window
        )

    t_both, y_both =
        _sd_trim_curve(
            t_both,
            y_both,
            x_window
        )


    # --------------------------------------------------------
    # Ground-truth effect
    # --------------------------------------------------------

    t_true,
    y_true,
    true_effect_label =
        _sd_true_effect(
            cfg,
            target;
            target_peak = true_effect_peak,
            use_magnitude = use_effect_magnitude
        )

    t_true, y_true =
        _sd_trim_curve(
            t_true,
            y_true,
            x_window
        )


    # --------------------------------------------------------
    # Shared x limits
    # --------------------------------------------------------

    x_min, x_max =
        Float64.(x_window)

    x_limits =
        (x_min, x_max)

        x_ticks = (
            [-0.1, 0.0, 0.1, 0.2, 0.3, 0.4],
            ["-0.1", "0", "0.1", "0.2", "0.3", "0.4"]
        )


    # --------------------------------------------------------
    # Shared y limits
    # Only use visible data
    # --------------------------------------------------------

    all_y_values =
        vcat(
            y_clean,
            y_overlap,
            y_confound,
            y_both
        )

    if show_true_effect
        all_y_values =
            vcat(
                all_y_values,
                y_true
            )
    end

    all_y_values =
        all_y_values[
            isfinite.(all_y_values)
        ]

    isempty(all_y_values) &&
        error(
            "No finite decoding scores were found " *
            "inside x_window = $x_window."
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
        y_max + y_padding
    )


    # --------------------------------------------------------
    # Figure and four axes
    # --------------------------------------------------------

    fig = Figure(
        size = (1080, 800),
        backgroundcolor = :white
    )


    ax_clean = _sd_make_axis(
        fig[1, 1];

        title = "Clean",
        ylabel = ylabel_text,

        show_x = false,
        show_y = true,

        x_limits = x_limits,
        y_limits = y_limits,
        x_ticks = x_ticks
    )


    ax_overlap = _sd_make_axis(
        fig[1, 2];

        title = "Overlap",
        ylabel = ylabel_text,

        show_x = false,
        show_y = false,

        x_limits = x_limits,
        y_limits = y_limits,
        x_ticks = x_ticks
    )


    ax_confound = _sd_make_axis(
        fig[2, 1];

        title = "Confound",
        ylabel = ylabel_text,

        show_x = true,
        show_y = true,

        x_limits = x_limits,
        y_limits = y_limits,
        x_ticks = x_ticks
    )


    ax_biased_cases = _sd_make_axis(
        fig[2, 2];

        title = "Overlap, Confound, and Both",
        ylabel = ylabel_text,

        show_x = true,
        show_y = false,

        x_limits = x_limits,
        y_limits = y_limits,
        x_ticks = x_ticks
    )


    # --------------------------------------------------------
    # Top-left: Clean
    # --------------------------------------------------------

    lines!(
        ax_clean,
        t_clean,
        y_clean;

        color = case_colors["clean"],
        linestyle = case_linestyles["clean"],
        linewidth = case_linewidths["clean"]
    )

    if show_true_effect &&
       ("clean" in true_effect_cases)

        lines!(
            ax_clean,
            t_true,
            y_true;

            color = true_effect_color,
            linestyle = true_effect_linestyle,
            linewidth = true_effect_linewidth
        )
    end


    # --------------------------------------------------------
    # Top-right: Overlap
    # --------------------------------------------------------

    lines!(
        ax_overlap,
        t_overlap,
        y_overlap;

        color = case_colors["overlap"],
        linestyle = case_linestyles["overlap"],
        linewidth = case_linewidths["overlap"]
    )

    if show_true_effect &&
       ("overlap" in true_effect_cases)

        lines!(
            ax_overlap,
            t_true,
            y_true;

            color = true_effect_color,
            linestyle = true_effect_linestyle,
            linewidth = true_effect_linewidth
        )
    end


    # --------------------------------------------------------
    # Bottom-left: Confound
    # --------------------------------------------------------

    lines!(
        ax_confound,
        t_confound,
        y_confound;

        color = case_colors["confound"],
        linestyle = case_linestyles["confound"],
        linewidth = case_linewidths["confound"]
    )

    if show_true_effect &&
       ("confound" in true_effect_cases)

        lines!(
            ax_confound,
            t_true,
            y_true;

            color = true_effect_color,
            linestyle = true_effect_linestyle,
            linewidth = true_effect_linewidth
        )
    end


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
        linewidth = case_linewidths["overlap"]
    )

    lines!(
        ax_biased_cases,
        t_confound,
        y_confound;

        color = case_colors["confound"],
        linestyle = case_linestyles["confound"],
        linewidth = case_linewidths["confound"]
    )

    lines!(
        ax_biased_cases,
        t_both,
        y_both;

        color = case_colors["both"],
        linestyle = case_linestyles["both"],
        linewidth = case_linewidths["both"]
    )

    if show_true_effect &&
       ("both" in true_effect_cases)

        lines!(
            ax_biased_cases,
            t_true,
            y_true;

            color = true_effect_color,
            linestyle = true_effect_linestyle,
            linewidth = true_effect_linewidth
        )
    end


    # --------------------------------------------------------
    # Shared axes
    # --------------------------------------------------------

    linkxaxes!(
        ax_clean,
        ax_overlap,
        ax_confound,
        ax_biased_cases
    )

    linkyaxes!(
        ax_clean,
        ax_overlap,
        ax_confound,
        ax_biased_cases
    )


    # --------------------------------------------------------
    # Main title
    # --------------------------------------------------------

    Label(
        fig[0, 1:2],
        figure_title;

        fontsize = 25,
        font = :bold,
        padding = (0, 0, 8, 8)
    )


    # --------------------------------------------------------
    # Shared legend
    # --------------------------------------------------------

    legend_elements =
        LineElement[
            LineElement(
                color = case_colors["clean"],
                linestyle = case_linestyles["clean"],
                linewidth = case_linewidths["clean"]
            ),

            LineElement(
                color = case_colors["overlap"],
                linestyle = case_linestyles["overlap"],
                linewidth = case_linewidths["overlap"]
            ),

            LineElement(
                color = case_colors["confound"],
                linestyle = case_linestyles["confound"],
                linewidth = case_linewidths["confound"]
            ),

            LineElement(
                color = case_colors["both"],
                linestyle = case_linestyles["both"],
                linewidth = case_linewidths["both"]
            )
        ]


    legend_labels = [
        "Clean",
        "Overlap",
        "Confound",
        "Both"
    ]


    if show_true_effect

        push!(
            legend_elements,

            LineElement(
                color = true_effect_color,
                linestyle = true_effect_linestyle,
                linewidth = true_effect_linewidth
            )
        )

        push!(
            legend_labels,
            true_effect_label
        )
    end


    Legend(
        fig[3, 1:2],
        legend_elements,
        legend_labels;

        orientation = :horizontal,
        framevisible = false,

        labelsize = 13,
        patchsize = (30, 12),

        tellheight = true
    )


    # --------------------------------------------------------
    # Layout spacing
    # --------------------------------------------------------

    colgap!(
        fig.layout,
        28
    )

    rowgap!(
        fig.layout,
        18
    )

    rowsize!(
        fig.layout,
        0,
        Auto(0.10)
    )

    rowsize!(
        fig.layout,
        3,
        Auto(0.10)
    )


    return fig
end