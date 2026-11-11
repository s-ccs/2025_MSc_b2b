function plot_roamm_b2b(
    scores_all,
    scores_mean;
    title = "ROAMM — B2B"
)

    targets = [
        "word_length",
        "frequency_zipf",
        "word_surprisal"
    ]


    titles = [
        "Word length",
        "Zipf frequency",
        "Word surprisal"
    ]

    n_subjects = length(unique(scores_all.subject))

    fig = Figure(size = (1500, 460))
    colgap!(fig.layout, 30)
    for (i, target) in enumerate(targets)

        sub = scores_all[scores_all.coefname .== target, :]
        sm = sort(
            scores_mean[scores_mean.coefname .== target, :],
            :time
        )

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Absolute B2B coefficient (log scale)" : "",
            xticks = -0.1:0.1:1.0,

            # Shared log scale with scientific notation
            yscale = log10,
            yticks = (
                [1e-8, 1e-6, 1e-4, 1e-2],
                ["10⁻⁸", "10⁻⁶", "10⁻⁴", "10⁻²"]
            ),
            yticklabelsvisible = i == 1,
            yticksvisible = i == 1,

            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        # Individual subjects
        for subject in unique(sub.subject)
            ss = sort(
                sub[sub.subject .== subject, :],
                :time
            )

            lines!(
                ax,
                ss.time,
                ss.estimate,
                color = ("#85B8D8", 0.18),
                linewidth = 1
            )
        end

        # 10% winsorized mean
        lines!(
            ax,
            sm.time,
            sm.winsor20,
            color = "#246B9A",
            linewidth = 3
        )
        
        # Median
        lines!(
            ax,
            sm.time,
            sm.estimate_median,
            #color = :darkorange,
            color = "#D66A58",
            linestyle = :dash,
            linewidth = 2
        )
        # Fixation onset
        vlines!(
            ax,
            [0.0],
            color = (:black, 0.55),
            linestyle = :dash,
            linewidth = 1
        )

        # Identical axis limits for all three predictors
        xlims!(ax, -0.1, 1.0)
        ylims!(ax, 1e-8, 1.0)
    end

    Label(
        fig[0, 1:3],
        "$title (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

    return fig
end

function plot_roamm_b2b_heatmap(scores_all; title = "ROAMM — B2B — subject-level heatmaps")

    targets = [
        "word_length_c",
        "frequency_zipf_c",
        "word_surprisal_c"
    ]

    titles = [
        "Word length",
        "Zipf frequency",
        "Word surprisal"
    ]

    fig = Figure(size = (1500, 650))
    subjects = sort(unique(scores_all.subject))

        # ---- global color scale across all 3 predictors ----
        vals = Float64[]
        for target in targets
            sub = scores_all[scores_all.coefname .== target, :]
            append!(vals, sub.estimate[sub.estimate .> 0])
        end
    
        cmin = minimum(vals)
        cmax = maximum(vals)
    
        hms = Any[]
    for (i, target) in enumerate(targets)

        sub = scores_all[scores_all.coefname .== target, :]
        times = sort(unique(sub.time))

        # Matrix: time × subjects
        mat = Matrix{Float64}(undef, length(times), length(subjects))

        for (j, subject) in enumerate(subjects)
            ss = sort(sub[sub.subject .== subject, :], :time)

            @assert ss.time == times "Time mismatch for $subject, $target"

            mat[:, j] = ss.estimate
        end

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Subject index" : "",
            xticks = -0.1:0.1:1.0,
            yticks = 1:5:length(subjects),
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )
        mat_log = map(x -> x > 0 ? x : NaN, mat)
        hm = heatmap!(
            ax,
            times,
            1:length(subjects),
            mat_log,
            colormap = :viridis,
            colorscale=log10
        )

        vlines!(
            ax, [0.0],
            color = :white,
            linestyle = :dash,
            linewidth = 1.5
        )

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, 0.5, length(subjects) + 0.5)

        Colorbar(
            fig[2, i],
            hm,
            vertical = false,
            label = "Absolute B2B coefficient"
        )
    end

    Label(
        fig[0, 1:3],
        title,
        fontsize = 22,
        font = :bold
    )

    return fig
end

function plot_roamm_b2b_mean(
    scores_summary;
    title = "ROAMM — B2B",
    ylim = nothing
)
    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    # Handle centered predictor names in two-step B2B
    coefnames = unique(String.(scores_summary.coefname))
    if "word_length_c" in coefnames
        targets = targets .* "_c"
    end

    fig = Figure(size = (1500, 420))
    colgap!(fig.layout, 35)

    # Shared y-limits across the 3 predictors
    if ylim === nothing
        ymax = maximum(scores_summary.estimate_mean)
        ylim = (0.0, 1.1 * ymax)
    end

    for (i, target) in enumerate(targets)

        sm = sort(
            scores_summary[scores_summary.coefname .== target, :],
            :time
        )

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Absolute B2B coefficient" : "",
            xticks = -0.1:0.2:1.0,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        # Group mean only
        lines!(
            ax,
            sm.time,
            sm.estimate_mean,
            color = "#0072B2",
            linewidth = 2.5
        )

        # Fixation onset
        vlines!(
            ax, [0.0],
            color = (:black, 0.4),
            linestyle = :dash,
            linewidth = 1
        )

        # Zero baseline
        hlines!(
            ax, [0.0],
            color = (:gray, 0.4),
            linestyle = :dash,
            linewidth = 1
        )

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, ylim...)
    end

    Label(
        fig[0, 1:3],
        title,
        fontsize = 22,
        font = :bold
    )

    return fig
end