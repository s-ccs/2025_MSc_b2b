function plot_roamm_decoding(scores_all, scores_summary;
    title = "ROAMM — rERP decoding"
)
    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    n_subjects = length(unique(scores_all.subject))
    fig = Figure(size = (1500, 460))

    for (i, target) in enumerate(targets)

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Correlation (r)" : "",
            xticks = -0.1:0.1:1.0,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        sub = scores_all[scores_all.target .== target, :]

        # Individual subjects
        for subject in unique(sub.subject)
            ss = sort(sub[sub.subject .== subject, :], :time)

            lines!(
                ax,
                ss.time,
                ss.r,
                color = ("#85B8D8", 0.18),
                linewidth = 1
            )
        end

        # Group summary
        sm = sort(
            scores_summary[scores_summary.target .== target, :],
            :time
        )

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
            sm.r_median,
            #color = :darkorange,
            color = "#D66A58",
            linestyle = :dash,
            linewidth = 2
        )
        # Fixation onset
        vlines!(
            ax, [0.0],
            color = (:black, 0.55),
            linestyle = :dash,
            linewidth = 1
        )

        # Chance-level reference
        hlines!(
            ax, [0.0],
            color = (:gray, 0.5),
            linestyle = :dash,
            linewidth = 1
        )

        xlims!(ax, -0.1, 1.0)
    end

    Label(
        fig[0, 1:3],
        "$title (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

    return fig
end

function plot_roamm_decoding_log(
    scores_all,
    scores_summary;
    title = "ROAMM — Decoding"
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

    # Average repeated scores within each subject/time
    subject_scores = combine(
        groupby(scores_all, [:subject, :target, :time]),
        :r => mean => :r
    )

    for (i, target) in enumerate(targets)

        sub = subject_scores[
            subject_scores.target .== target, :
        ]

        sm = sort(
            scores_summary[
                scores_summary.target .== target, :
            ],
            :time
        )

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Absolute decoding correlation |r| (log scale)" : "",
            xticks = -0.1:0.1:1.0,

            yscale = log10,

            yticks = (
                [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0],
                ["10⁻⁵", "10⁻⁴", "10⁻³", "10⁻²", "10⁻¹", "10⁰"]
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

            values = abs.(ss.r)
            values[values .<= 0] .= NaN

            lines!(
                ax,
                ss.time,
                values,
                color = ("#B6D8EF", 0.15),
                linewidth = 1
            )
        end

        # 20% winsorized mean
        winsor_values = abs.(sm.winsor10)
        winsor_values[winsor_values .<= 0] .= NaN

        lines!(
            ax,
            sm.time,
            winsor_values,
            color = "#246B9A",
            linewidth = 3
        )

        # Median
        median_values = abs.(sm.r_median)
        median_values[median_values .<= 0] .= NaN

        lines!(
            ax,
            sm.time,
            median_values,
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

        # Shared axis limits
        xlims!(ax, -0.1, 1.0)
        ylims!(ax, 1e-5, 1.0)
    end

    Label(
        fig[0, 1:3],
        "$title (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

    return fig
end

function plot_roamm_decoding_mean(
    scores_summary;
    title = "ROAMM — Decoding",
    ylim = nothing
)
    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    fig = Figure(size = (1500, 420))
    colgap!(fig.layout, 35)

    # Shared y-axis limits across all predictors
    if ylim === nothing
        rmax = maximum(abs.(scores_summary.r_mean))
        ylim = (-1.1 * rmax, 1.1 * rmax)
    end

    for (i, target) in enumerate(targets)

        sm = sort(
            scores_summary[scores_summary.target .== target, :],
            :time
        )

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Decoding correlation (r)" : "",
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
            sm.r_mean,
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

        # Zero correlation
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

function plot_roamm_decoding_heatmap(
    scores_all;
    title = "ROAMM — Decoding — subject-level heatmaps"
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

    # Average repeated observations within each subject and time
    scores = combine(
        groupby(scores_all, [:subject, :target, :time]),
        :r => mean => :r
    )

    fig = Figure(size = (1500, 650))
    subjects = sort(unique(scores.subject))

    for (i, target) in enumerate(targets)

        sub = scores[scores.target .== target, :]
        times = sort(unique(sub.time))

        # Matrix: time × subjects
        mat = Matrix{Float64}(
            undef,
            length(times),
            length(subjects)
        )

        for (j, subject) in enumerate(subjects)

            ss = sort(
                sub[sub.subject .== subject, :],
                :time
            )

            @assert ss.time == times "Time mismatch for $subject, $target"

            mat[:, j] = ss.r
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

        # Symmetric color range centered at zero
        rmax = maximum(abs.(mat))

        hm = heatmap!(
            ax,
            times,
            1:length(subjects),
            mat,
            colormap = :balance,
            colorrange = (-rmax, rmax)
        )

        # Fixation onset
        vlines!(
            ax, [0.0],
            color = :black,
            linestyle = :dash,
            linewidth = 1.5
        )

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, 0.5, length(subjects) + 0.5)

        Colorbar(
            fig[2, i],
            hm,
            vertical = false,
            label = "Decoding correlation (r)"
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