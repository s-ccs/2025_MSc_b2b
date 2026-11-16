# ========================================
# Standard decoding plotting
# ========================================


# ============================================================
# Prepare decoding scores
# ============================================================

function _sd_prepare_scores(
    scores::AbstractDataFrame,
    cfg,
    x_window,
)
    df = DataFrame(scores)

    if !(:time in propertynames(df))
        :timepoint in propertynames(df) ||
            error("scores must contain :time or :timepoint")

        df.time =
            (Float64.(df.timepoint) .- 1.0) ./ cfg.sfreq
    end

    x_min, x_max = x_window

    filter!(
        :time => t -> x_min <= t <= x_max,
        df,
    )

    sort!(df, [:case, :time])

    return df
end


# ============================================================
# Ground-truth timing window
#
# Uses the actual non-zero support of the simulated
# UnfoldSim component. No arbitrary amplitude threshold.
# ============================================================

function _sd_groundtruth_window(
    cfg,
    target::Symbol,
)
    basis, β, label =
        if target == :condition

            (
                UnfoldSim.n170(; sfreq = cfg.sfreq),
                cfg.β_condition,
                "Ground-truth N170 timing",
            )

        elseif target == :continuous

            (
                UnfoldSim.p300(; sfreq = cfg.sfreq),
                cfg.β_continuous,
                "Ground-truth P300 timing",
            )

        else
            error(
                "Unknown target: $target. " *
                "Use :condition or :continuous."
            )
        end

    # If the simulated effect is actually zero,
    # there is no ground-truth effect window to display.
    iszero(β) && return nothing

    active =
        findall(
            x -> !iszero(x),
            basis,
        )

    isempty(active) && return nothing

    return (
        start =
            (first(active) - 1) / cfg.sfreq,

        stop =
            (last(active) - 1) / cfg.sfreq,

        label =
            label,
    )
end


# ============================================================
# Convert ground-truth timing to the format expected by
# UnfoldMakie significance / vspan plotting
# ============================================================

function _sd_window_dataframe(window)

    isnothing(window) &&
        return nothing

    return DataFrame(
        from = [window.start],
        to = [window.stop],
    )
end


# ============================================================
# Plot one panel with UnfoldMakie
# ============================================================

function _sd_plot_panel!(
    position,
    df;
    title,
    color,
    x_limits,
    y_limits,
    x_ticks,
    show_x = true,
    show_y = true,
    significance = nothing,
)

    UnfoldMakie.plot_erp!(
        position,
        df;

        mapping = (;
            x = :time,
            y = :r,
        ),

        visual = (;
            color = color,
            linewidth = 2.6,
        ),

        significance =
            significance,

        # Yes, this keyword is currently spelled
        # "sigifnicance_visual" in UnfoldMakie.
        sigifnicance_visual =
            :vspan,

        significance_vspan = (;
            alpha = 0.10,
        ),

        axis = (;
            title = title,
            titlesize = 16,
            titlefont = :bold,

            xlabel =
                show_x ?
                "Time [s]" :
                "",

            ylabel =
                show_y ?
                "Decoding r" :
                "",

            xticks =
                x_ticks,

            limits = (
                x_limits[1],
                x_limits[2],
                y_limits[1],
                y_limits[2],
            ),

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
        ),

        layout = (;
            show_legend = false,
        ),
    )
end


# ============================================================
# Plot overlap + confound + both together
# ============================================================

function _sd_plot_biased_panel!(
    position,
    df;
    x_limits,
    y_limits,
    x_ticks,
    significance = nothing,
)

    # Prefixes make the categorical plotting order explicit.
    case_names = Dict(
        "overlap"  => "1 Overlap",
        "confound" => "2 Confound",
        "both"     => "3 Both",
    )

    plot_df = copy(df)

    plot_df.case_plot =
        [
            case_names[string(case)]
            for case in plot_df.case
        ]

    UnfoldMakie.plot_erp!(
        position,
        plot_df;

        mapping = (;
            x = :time,
            y = :r,
            color = :case_plot,
            group = :case_plot,
        ),

        visual = (;
            color = [
                :darkorange,
                :seagreen,
                :deeppink,
            ],
            linewidth = 2.6,
        ),

        significance =
            significance,

        sigifnicance_visual =
            :vspan,

        significance_vspan = (;
            alpha = 0.10,
        ),

        axis = (;
            title =
                "Overlap, Confound, and Both",

            titlesize =
                16,

            titlefont =
                :bold,

            xlabel =
                "Time [s]",

            ylabel =
                "",

            xticks =
                x_ticks,

            limits = (
                x_limits[1],
                x_limits[2],
                y_limits[1],
                y_limits[2],
            ),

            yticklabelsvisible =
                false,

            yticksvisible =
                false,

            ylabelvisible =
                false,

            xgridvisible =
                false,

            ygridvisible =
                false,

            topspinevisible =
                false,

            rightspinevisible =
                false,
        ),

        layout = (;
            show_legend = false,
        ),
    )
end


# ============================================================
# Main plotting function
# ============================================================

function plot_standard_decoding_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    x_window = (-0.1, 0.43),
    show_true_window::Bool = true,
    true_window_cases = (
        "clean",
        "overlap",
        "confound",
        "both",
    ),
)

    # --------------------------------------------------------
    # Labels
    # --------------------------------------------------------

    figure_title =
        if target == :condition
            "Condition decoding"

        elseif target == :continuous
            "Continuous decoding"

        else
            error(
                "Unknown target: $target. " *
                "Use :condition or :continuous."
            )
        end


    # --------------------------------------------------------
    # Prepare tidy plotting data
    # --------------------------------------------------------

    df =
        _sd_prepare_scores(
            scores,
            cfg,
            x_window,
        )


    # --------------------------------------------------------
    # Shared x axis
    # --------------------------------------------------------

    x_min, x_max =
        Float64.(x_window)

    x_limits =
        (x_min, x_max)

    first_tick =
        ceil(x_min * 10) / 10

    last_tick =
        floor(x_max * 10) / 10

    x_ticks =
        first_tick:0.1:last_tick


    # --------------------------------------------------------
    # Shared y axis
    # --------------------------------------------------------

    y =
        Float64.(
            df.r[
                isfinite.(df.r)
            ]
        )

    isempty(y) &&
        error(
            "No finite decoding values inside x_window."
        )

    y_min =
        min(
            minimum(y),
            0.0,
        )

    y_max =
        max(
            maximum(y),
            0.0,
        )

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
    # Ground-truth timing window
    # --------------------------------------------------------

    true_window =
        show_true_window ?
        _sd_groundtruth_window(
            cfg,
            target,
        ) :
        nothing

    timing_df =
        _sd_window_dataframe(
            true_window
        )

    function timing_for(case)

        if show_true_window &&
           !isnothing(timing_df) &&
           case in true_window_cases

            return timing_df

        else
            return nothing
        end
    end


    # --------------------------------------------------------
    # Extract plotting subsets
    # --------------------------------------------------------

    clean =
        filter(
            :case => ==("clean"),
            df,
        )

    overlap =
        filter(
            :case => ==("overlap"),
            df,
        )

    confound =
        filter(
            :case => ==("confound"),
            df,
        )

    both =
        filter(
            :case => ==("both"),
            df,
        )

    biased =
        vcat(
            overlap,
            confound,
            both,
        )


    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig =
        Figure(
            size = (1080, 800),
            backgroundcolor = :white,
        )


    # --------------------------------------------------------
    # Four panels
    # --------------------------------------------------------

    _sd_plot_panel!(
        fig[1, 1],
        clean;

        title =
            "Clean",

        color =
            :dodgerblue,

        x_limits =
            x_limits,

        y_limits =
            y_limits,

        x_ticks =
            x_ticks,

        show_x =
            false,

        show_y =
            true,

        significance =
            timing_for("clean"),
    )


    _sd_plot_panel!(
        fig[1, 2],
        overlap;

        title =
            "Overlap",

        color =
            :darkorange,

        x_limits =
            x_limits,

        y_limits =
            y_limits,

        x_ticks =
            x_ticks,

        show_x =
            false,

        show_y =
            false,

        significance =
            timing_for("overlap"),
    )


    _sd_plot_panel!(
        fig[2, 1],
        confound;

        title =
            "Confound",

        color =
            :seagreen,

        x_limits =
            x_limits,

        y_limits =
            y_limits,

        x_ticks =
            x_ticks,

        show_x =
            true,

        show_y =
            true,

        significance =
            timing_for("confound"),
    )


    _sd_plot_biased_panel!(
        fig[2, 2],
        biased;

        x_limits =
            x_limits,

        y_limits =
            y_limits,

        x_ticks =
            x_ticks,

        significance =
            timing_for("both"),
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
            (0, 0, 8, 8),
    )


    # --------------------------------------------------------
    # Shared legend
    # --------------------------------------------------------

    legend_elements =
        Any[
            LineElement(
                color = :dodgerblue,
                linewidth = 2.6,
            ),

            LineElement(
                color = :darkorange,
                linewidth = 2.6,
            ),

            LineElement(
                color = :seagreen,
                linewidth = 2.6,
            ),

            LineElement(
                color = :deeppink,
                linewidth = 2.6,
            ),
        ]

    legend_labels =
        [
            "Clean",
            "Overlap",
            "Confound",
            "Both",
        ]


    if !isnothing(true_window)

        push!(
            legend_elements,

            PolyElement(
                color =
                    (:gray30, 0.10),
            )
        )

        push!(
            legend_labels,
            true_window.label,
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
            true,
    )


    # --------------------------------------------------------
    # Layout
    # --------------------------------------------------------

    colgap!(
        fig.layout,
        28,
    )

    rowgap!(
        fig.layout,
        18,
    )

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