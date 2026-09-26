
function plot_correlated_b2b(
    scores,
    cfg::CorrelatedContinuousConfig;
    y::Symbol = :estimate,
    panel::Symbol = :coefname,
    ylabel::String = "B2B estimate magnitude",
    x_window = (-0.1, 1.0),
    title = "Two-step B2B on correlated predictors"
)

    panels = [
        "continuous1",
        "continuous2",
        "continuous3"
    ]

    panel_titles = Dict(
        "continuous1" => "Continuous 1 ($(round(Int, cfg.peak1 * 1000)) ms)",
        "continuous2" => "Continuous 2 ($(round(Int, cfg.peak2 * 1000)) ms)",
        "continuous3" => "Continuous 3 ($(round(Int, cfg.peak3 * 1000)) ms)"
    )

    cases = [
        "no_overlap",
        "overlap"
    ]

    model_values =
        Float64.(scores[!, y])

    model_values =
        model_values[
            isfinite.(model_values)
        ]

    y_min =
        min(
            0.0,
            minimum(model_values)
        )

    y_max =
        maximum(model_values)

    y_range =
        max(
            y_max - y_min,
            eps()
        )

    y_upper =
        y_max +
        0.12 * y_range

    y_bar =
        y_max +
        0.07 * y_range

    fig =
        Figure(
            size = (1500, 450)
        )

    axes = Axis[]

    for (i, panel_name) in enumerate(panels)

        ax =
            Axis(
                fig[1, i],
                title = panel_titles[panel_name],
                xlabel = "Time [s]",
                ylabel =
                    i == 1 ?
                    ylabel :
                    "",
                topspinevisible = false,
                rightspinevisible = false,
                xgridvisible = false,
                ygridvisible = false,
                xminorgridvisible = false,
                yminorgridvisible = false,
            )

        panel_mask =
            string.(scores[!, panel]) .==
            panel_name

        for case_name in cases

            mask =
                panel_mask .&
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
                sub[!, y],
                label =
                    replace(
                        case_name,
                        "_" => " "
                    )
            )
        end

        truth =
            _correlated_truth_window(
                cfg,
                panel_name
            )

        lines!(
            ax,
            [truth.start, truth.stop],
            [y_bar, y_bar],
            color = (:gray30, 0.45),
            linewidth = 7,
            label = "ground-truth timing"
        )

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

    linkyaxes!(
        axes...
    )

    Legend(
        fig[2, 1:3],
        axes[1];
        orientation = :horizontal,
        framevisible = false
    )

    Label(
        fig[0, 1:3],
        title,
        fontsize = 24,
        font = :bold,
        padding = (0, 0, 10, 0)
    )

    return fig
end