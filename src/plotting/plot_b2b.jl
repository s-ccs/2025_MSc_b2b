# ========================================
# B2B plotting
# ========================================


# ============================================================
# 1. Select condition or continuous coefficients
# ============================================================

function _b2b_target_mask(
    scores::AbstractDataFrame,
    target::Symbol
)
    :coefname in propertynames(scores) ||
        error(
            "`scores` must contain :coefname. " *
            "Available columns: $(propertynames(scores))"
        )

    coefnames = string.(scores.coefname)

    mask =
        if target == :continuous

            occursin.(
                "continuous",
                coefnames
            )

        elseif target in (:condition, :condition_num)

            occursin.(
                "condition",
                coefnames
            )

        else

            error(
                "`target` must be " *
                ":condition, :condition_num, or :continuous."
            )
        end

    any(mask) ||
        error(
            "No rows found for target = $target. " *
            "Available coefnames: $(unique(coefnames))"
        )

    return mask
end


# ============================================================
# 2. Extract one case curve
# ============================================================

function _b2b_get_curve(
    scores::AbstractDataFrame,
    case_name::AbstractString;
    score_col::Symbol,
    sfreq::Real,
    use_abs::Bool
)
    :case in propertynames(scores) ||
        error("`scores` must contain :case.")

    score_col in propertynames(scores) ||
        error(
            "`scores` must contain $score_col. " *
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
    # Prefer real time
    # --------------------------------------------------------

    if :time in propertynames(scores)

        sub =
            DataFrame(
                scores[
                    case_mask,
                    [:time, score_col]
                ]
            )

        sort!(
            sub,
            :time
        )

        time =
            Float64.(sub.time)


    # --------------------------------------------------------
    # Otherwise derive time from timepoint
    # --------------------------------------------------------

    elseif :timepoint in propertynames(scores)

        sub =
            DataFrame(
                scores[
                    case_mask,
                    [:timepoint, score_col]
                ]
            )

        sort!(
            sub,
            :timepoint
        )

        time =
            (
                Float64.(sub.timepoint) .- 1.0
            ) ./ sfreq

    else

        error(
            "`scores` must contain either :time or :timepoint."
        )

    end


    values =
        Float64.(sub[!, score_col])

    if use_abs
        values =
            abs.(values)
    end


    return (
        time = time,
        values = values
    )
end


# ============================================================
# 3. Trim curve to plotting window
# ============================================================

function _b2b_trim_curve(
    curve,
    x_window
)
    x_min, x_max =
        x_window

    mask =
        (curve.time .>= x_min) .&
        (curve.time .<= x_max)

    return (
        time =
            Float64.(curve.time[mask]),

        values =
            Float64.(curve.values[mask])
    )
end


# ============================================================
# 4. Get ground-truth effect timing window
#
# threshold = 0.20 means:
# keep samples where the true component is >= 20% of its peak.
# ============================================================

function _b2b_effect_window(
    cfg,
    target::Symbol;
    threshold::Real = 0.20
)
    waveform, label =
        if target in (:condition, :condition_num)

            (
                cfg.β_condition .*
                UnfoldSim.n170(
                    ;
                    sfreq = cfg.sfreq
                ),

                "Ground-truth N170 timing"
            )

        elseif target == :continuous

            (
                cfg.β_continuous .*
                UnfoldSim.p300(
                    ;
                    sfreq = cfg.sfreq
                ),

                "Ground-truth P300 timing"
            )

        else

            error(
                "`target` must be " *
                ":condition, :condition_num, or :continuous."
            )

        end


    waveform =
        abs.(Float64.(waveform))

    peak =
        maximum(waveform)

    peak > 0 ||
        error(
            "Ground-truth target effect has zero amplitude."
        )


    time =
        (0:length(waveform)-1) ./ cfg.sfreq


    mask =
        waveform .>=
        threshold * peak

    any(mask) ||
        error(
            "No waveform samples passed threshold = $threshold."
        )


    return (
        start =
            minimum(time[mask]),

        stop =
            maximum(time[mask]),

        label =
            label
    )
end


# ============================================================
# 5. Nice x-axis ticks
# ============================================================

function _b2b_xticks(
    x_window
)
    x_min, x_max =
        Float64.(x_window)

    first_tick =
        ceil(x_min * 10) / 10

    last_tick =
        floor(x_max * 10) / 10

    values =
        collect(
            first_tick:0.1:last_tick
        )

    labels =
        [
            isapprox(
                value,
                0.0;
                atol = 1e-8
            ) ?
                "0" :
                string(
                    round(
                        value;
                        digits = 1
                    )
                )

            for value in values
        ]

    return (
        values,
        labels
    )
end


# ============================================================
# 6. Main B2B plotting function
# ============================================================

function plot_b2b_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    score_col::Symbol = :estimate,
    use_abs::Bool = true,
    x_window = (-0.1, 0.43),
    timing_threshold::Real = 0.20,
    show_timing_window::Bool = true,
    timing_window_cases = ("overlap", "confound")
)

    # --------------------------------------------------------
    # Select target
    # --------------------------------------------------------

    target_scores =
        DataFrame(
            scores[
                _b2b_target_mask(
                    scores,
                    target
                ),
                :
            ]
        )


    # --------------------------------------------------------
    # Extract four cases
    # --------------------------------------------------------

    curves =
        Dict{String, Any}()

    for case_name in (
        "clean",
        "overlap",
        "confound",
        "both"
    )

        curve =
            _b2b_get_curve(
                target_scores,
                case_name;

                score_col = score_col,
                sfreq = cfg.sfreq,
                use_abs = use_abs
            )

        curves[case_name] =
            _b2b_trim_curve(
                curve,
                x_window
            )
    end


    # --------------------------------------------------------
    # Ground-truth target timing
    # --------------------------------------------------------

    timing =
        _b2b_effect_window(
            cfg,
            target;
            threshold = timing_threshold
        )


    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    condition_target =
        target in (
            :condition,
            :condition_num
        )

    figure_title =
        condition_target ?
            "Condition B2B" :
            "Continuous B2B"

    ylabel_text =
        use_abs ?
            "B2B estimate magnitude" :
            "B2B estimate"


    # --------------------------------------------------------
    # Colors / styles
    # --------------------------------------------------------

    case_colors =
        Dict(
            "clean" =>
                :dodgerblue,

            "overlap" =>
                :darkorange,

            "confound" =>
                :seagreen,

            "both" =>
                :deeppink
        )

    case_linestyles =
        Dict(
            "clean" =>
                :solid,

            "overlap" =>
                :solid,

            "confound" =>
                :solid,

            "both" =>
                :dash
        )

    case_linewidths =
        Dict(
            "clean" =>
                2.6,

            "overlap" =>
                2.6,

            "confound" =>
                2.6,

            "both" =>
                2.8
        )


    # --------------------------------------------------------
    # X axis
    # --------------------------------------------------------

    x_limits =
        (
            Float64(x_window[1]),
            Float64(x_window[2])
        )

    x_ticks =
        _b2b_xticks(
            x_window
        )


    # --------------------------------------------------------
    # Y limits
    #
    # IMPORTANT:
    # only B2B estimates determine the y-axis.
    # The shaded timing window does NOT affect scaling.
    # --------------------------------------------------------

    all_y_values =
        vcat(
            curves["clean"].values,
            curves["overlap"].values,
            curves["confound"].values,
            curves["both"].values
        )

    all_y_values =
        all_y_values[
            isfinite.(all_y_values)
        ]

    isempty(all_y_values) &&
        error(
            "No finite B2B values found."
        )


    y_min =
        use_abs ?
            0.0 :
            minimum(all_y_values)

    y_max =
        maximum(all_y_values)


    if y_max == y_min

        y_limits =
            use_abs ?
                (0.0, y_max + 0.1) :
                (y_min - 0.1, y_max + 0.1)

    else

        y_padding =
            0.06 * (y_max - y_min)

        y_limits =
            use_abs ?
                (
                    0.0,
                    y_max + y_padding
                ) :
                (
                    y_min - y_padding,
                    y_max + y_padding
                )

    end


    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig =
        Figure(
            size = (1080, 800),
            backgroundcolor = :white
        )


    # --------------------------------------------------------
    # Axis helper
    # --------------------------------------------------------

    function make_axis(
        position;
        title,
        show_x,
        show_y
    )

        ax =
            Axis(
                position;

                title =
                    title,

                titlesize =
                    16,

                titlefont =
                    :bold,

                titlegap =
                    8,

                xlabel =
                    show_x ?
                        "Time [s]" :
                        "",

                ylabel =
                    show_y ?
                        ylabel_text :
                        "",

                xlabelsize =
                    14,

                ylabelsize =
                    14,

                xticklabelsize =
                    12,

                yticklabelsize =
                    12,

                xticks =
                    x_ticks,

                xticklabelsvisible =
                    show_x,

                xticksvisible =
                    show_x,

                xlabelvisible =
                    show_x,

                yticklabelsvisible =
                    show_y,

                yticksvisible =
                    show_y,

                ylabelvisible =
                    show_y,

                xgridvisible =
                    false,

                ygridvisible =
                    false,

                topspinevisible =
                    false,

                rightspinevisible =
                    false,

                backgroundcolor =
                    :white
            )


        xlims!(
            ax,
            x_limits...
        )

        ylims!(
            ax,
            y_limits...
        )


        hlines!(
            ax,
            [0.0];

            color =
                (:gray45, 0.55),

            linestyle =
                :dash,

            linewidth =
                1.0
        )


        return ax
    end


    # --------------------------------------------------------
    # Four axes
    # --------------------------------------------------------

    ax_clean =
        make_axis(
            fig[1, 1];

            title =
                "Clean",

            show_x =
                false,

            show_y =
                true
        )


    ax_overlap =
        make_axis(
            fig[1, 2];

            title =
                "Overlap",

            show_x =
                false,

            show_y =
                false
        )


    ax_confound =
        make_axis(
            fig[2, 1];

            title =
                "Confound",

            show_x =
                true,

            show_y =
                true
        )


    ax_biased_cases =
        make_axis(
            fig[2, 2];

            title =
                "Overlap, Confound, and Both",

            show_x =
                true,

            show_y =
                false
        )


    # --------------------------------------------------------
    # Add timing window
    # --------------------------------------------------------

    function add_timing_window!(ax)

        y_top = y_limits[2]
        y_bar = y_top * 0.96
    
        lines!(
            ax,
            [timing.start, timing.stop],
            [y_bar, y_bar];
    
            color = (:gray30, 0.55),
            linewidth = 6
        )
    end


    # --------------------------------------------------------
    # Clean
    # --------------------------------------------------------

    if show_timing_window &&
       ("clean" in timing_window_cases)

        add_timing_window!(
            ax_clean
        )
    end


    lines!(
        ax_clean,
        curves["clean"].time,
        curves["clean"].values;

        color =
            case_colors["clean"],

        linestyle =
            case_linestyles["clean"],

        linewidth =
            case_linewidths["clean"]
    )


    # --------------------------------------------------------
    # Overlap
    # --------------------------------------------------------

    if show_timing_window &&
       ("overlap" in timing_window_cases)

        add_timing_window!(
            ax_overlap
        )
    end


    lines!(
        ax_overlap,
        curves["overlap"].time,
        curves["overlap"].values;

        color =
            case_colors["overlap"],

        linestyle =
            case_linestyles["overlap"],

        linewidth =
            case_linewidths["overlap"]
    )


    # --------------------------------------------------------
    # Confound
    # --------------------------------------------------------

    if show_timing_window &&
       ("confound" in timing_window_cases)

        add_timing_window!(
            ax_confound
        )
    end


    lines!(
        ax_confound,
        curves["confound"].time,
        curves["confound"].values;

        color =
            case_colors["confound"],

        linestyle =
            case_linestyles["confound"],

        linewidth =
            case_linewidths["confound"]
    )


    # --------------------------------------------------------
    # Bottom-right:
    # overlap + confound + both
    # --------------------------------------------------------

    if show_timing_window &&
       ("both" in timing_window_cases)

        add_timing_window!(
            ax_biased_cases
        )
    end


    lines!(
        ax_biased_cases,
        curves["overlap"].time,
        curves["overlap"].values;

        color =
            case_colors["overlap"],

        linestyle =
            case_linestyles["overlap"],

        linewidth =
            case_linewidths["overlap"]
    )


    lines!(
        ax_biased_cases,
        curves["confound"].time,
        curves["confound"].values;

        color =
            case_colors["confound"],

        linestyle =
            case_linestyles["confound"],

        linewidth =
            case_linewidths["confound"]
    )


    lines!(
        ax_biased_cases,
        curves["both"].time,
        curves["both"].values;

        color =
            case_colors["both"],

        linestyle =
            case_linestyles["both"],

        linewidth =
            case_linewidths["both"]
    )


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

        fontsize =
            25,

        font =
            :bold,

        padding =
            (0, 0, 8, 8)
    )


    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    legend_elements =
        Any[
            LineElement(
                color =
                    case_colors["clean"],

                linestyle =
                    case_linestyles["clean"],

                linewidth =
                    case_linewidths["clean"]
            ),

            LineElement(
                color =
                    case_colors["overlap"],

                linestyle =
                    case_linestyles["overlap"],

                linewidth =
                    case_linewidths["overlap"]
            ),

            LineElement(
                color =
                    case_colors["confound"],

                linestyle =
                    case_linestyles["confound"],

                linewidth =
                    case_linewidths["confound"]
            ),

            LineElement(
                color =
                    case_colors["both"],

                linestyle =
                    case_linestyles["both"],

                linewidth =
                    case_linewidths["both"]
            )
        ]


    legend_labels =
        [
            "Clean",
            "Overlap",
            "Confound",
            "Both"
        ]


    if show_timing_window

        push!(
            legend_elements,

            LineElement(
                color =
                    (:gray30, 0.55),
                linewidth = 6
            )
        )

        push!(
            legend_labels,
            timing.label
        )

    end


    Legend(
        fig[3, 1:2],
        legend_elements,
        legend_labels;

        orientation =
            :horizontal,

        framevisible =
            false,

        labelsize =
            13,

        patchsize =
            (30, 12),

        tellheight =
            true
    )


    # --------------------------------------------------------
    # Layout
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