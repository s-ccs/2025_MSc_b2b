# ========================================
# Correlated continuous decoding plotting
# ========================================

function _correlated_truth_window(
    cfg::CorrelatedContinuousConfig,
    target::AbstractString
)
    peak =
        if target == "continuous1"
            cfg.peak1
        elseif target == "continuous2"
            cfg.peak2
        elseif target == "continuous3"
            cfg.peak3
        else
            error("Unknown target: $target")
        end

    half_width = cfg.component_width / 2

    return (
        start = peak - half_width,
        stop = peak + half_width
    )
end

function plot_correlated_decoding(
    scores::AbstractDataFrame,
    cfg::CorrelatedContinuousConfig;
    title::String = "Decoding on correlated predictors",
    x_window = (-0.1, 1.0),
    #truth_threshold::Real = 0.20,
    show_truth::Bool = true
)

    # --------------------------------------------------------
    # Checks
    # --------------------------------------------------------

    for col in (:time, :r, :case, :target)
        col in propertynames(scores) ||
            error(
                "`scores` must contain $col. " *
                "Available columns: $(propertynames(scores))"
            )
    end


    # --------------------------------------------------------
    # Target order
    # --------------------------------------------------------

    targets = [
        "continuous1",
        "continuous2",
        "continuous3"
    ]
    """
    target_titles = Dict(
        "continuous1" => "Continuous 1 (P100)",
        "continuous2" => "Continuous 2 (N170)",
        "continuous3" => "Continuous 3 (P300)"
    )
    """

    target_titles = Dict(
        "continuous1" => "Continuous 1 ($(round(Int, cfg.peak1 * 1000)) ms)",
        "continuous2" => "Continuous 2 ($(round(Int, cfg.peak2 * 1000)) ms)",
        "continuous3" => "Continuous 3 ($(round(Int, cfg.peak3 * 1000)) ms)"
    )

    available_targets =
        unique(
            string.(scores.target)
        )

    targets =
        [
            t for t in targets
            if t in available_targets
        ]

    isempty(targets) &&
        error(
            "No continuous1/2/3 targets found."
        )


    # --------------------------------------------------------
    # Cases
    # --------------------------------------------------------

    cases = [
        "no_overlap",
        "overlap"
    ]


    # --------------------------------------------------------
    # Shared y-axis
    # --------------------------------------------------------

    r_values =
        Float64.(scores.r)

    r_values =
        r_values[
            isfinite.(r_values)
        ]

    isempty(r_values) &&
        error(
            "No finite correlation values found."
        )

    y_min =
        min(
            0.0,
            minimum(r_values)
        )

    y_max =
        maximum(r_values)

    y_range =
        max(
            y_max - y_min,
            0.1
        )

    y_upper =
        y_max +
        0.12 * y_range

    y_bar =
        y_max +
        0.07 * y_range


    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    fig =
        Figure(
            size = (1500, 450)
        )

    axes = Axis[]


    # --------------------------------------------------------
    # Panels
    # --------------------------------------------------------

    for (i, target_name) in enumerate(targets)

        ax =
            Axis(
                fig[1, i],

                title =
                    target_titles[target_name],

                xlabel = "Time [s]",

                ylabel =
                    i == 1 ?
                    "Correlation (r)" :
                    "",

                topspinevisible = false,
                rightspinevisible = false,

                xgridvisible = false,
                ygridvisible = false,

                xminorgridvisible = false,
                yminorgridvisible = false,
            )


        target_mask =
            string.(scores.target) .==
            target_name


        for case_name in cases

            mask =
                target_mask .&
                (
                    string.(scores.case) .==
                    case_name
                )

            any(mask) ||
                continue

            sub =
                DataFrame(
                    scores[
                        mask,
                        :
                    ]
                )

            sort!(
                sub,
                :time
            )

            lines!(
                ax,
                sub.time,
                sub.r,
                label =
                    replace(
                        case_name,
                        "_" => " "
                    )
            )
        end


        # ----------------------------------------------------
        # Ground-truth timing bar
        # ----------------------------------------------------

        if show_truth

            truth =
                _correlated_truth_window(
                    cfg,
                    target_name
                    #threshold =
                        #truth_threshold
                )

            lines!(
                ax,
                [
                    truth.start,
                    truth.stop
                ],
                [
                    y_bar,
                    y_bar
                ],
                color =
                    (:gray30, 0.45),
                linewidth = 7,
                label =
                    "ground-truth timing"
            )
        end


        # ----------------------------------------------------
        # Limits
        # ----------------------------------------------------

        xlims!(
            ax,
            x_window...
        )

        ylims!(
            ax,
            y_min,
            y_upper
        )

        push!(
            axes,
            ax
        )
    end


    # --------------------------------------------------------
    # Shared axes
    # --------------------------------------------------------

    linkxaxes!(
        axes...
    )

    linkyaxes!(
        axes...
    )


    # --------------------------------------------------------
    # Legend
    # --------------------------------------------------------

    Legend(
        fig[2, 1:length(targets)],
        axes[1];
        orientation = :horizontal,
        framevisible = false
    )


    # --------------------------------------------------------
    # Main title
    # --------------------------------------------------------

    Label(
        fig[0, 1:length(targets)],
        title,
        fontsize = 24,
        font = :bold,
        padding = (0, 0, 10, 0)
    )


    return fig
end