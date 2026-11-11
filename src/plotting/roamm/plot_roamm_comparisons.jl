function plot_standard_vs_rerp(
    scores_standard,
    scores_rerp
)
    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    # Only compare subjects available in BOTH pipelines
    subjects = intersect(
        unique(scores_standard.subject),
        unique(scores_rerp.subject)
    )

    @info "Matched subjects" n = length(subjects)
    isempty(subjects) && error("No matching subjects found.")

    function summarize_matched(df)
        matched = filter(row -> row.subject in subjects, df)

        # Average any repeated scores within subject/time
        subject_mean = combine(
            groupby(matched, [:subject, :target, :time]),
            :r => mean => :r
        )

        # Mean across subjects
        combine(
            groupby(subject_mean, [:target, :time]),
            :r => mean => :winsor20
        )
    end

    standard_mean = summarize_matched(scores_standard)
    rerp_mean = summarize_matched(scores_rerp)

    fig = Figure(size = (1500, 460))
    colgap!(fig.layout, 35)

    for (i, target) in enumerate(targets)
        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Decoding correlation (r)" : "",
            xticks = -0.1:0.2:1.0,
            yticklabelsvisible = i == 1,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        std = sort(
            standard_mean[standard_mean.target .== target, :],
            :time
        )

        rer = sort(
            rerp_mean[rerp_mean.target .== target, :],
            :time
        )

        lines!(
            ax, std.time, std.winsor20;
            color = "#0072B2",
            linewidth = 3,
            label = "Standard decoding"
        )

        lines!(
            ax, rer.time, rer.winsor20;
            color = "#D55E00",
            linewidth = 3,
            label = "rERP decoding"
        )

        hlines!(ax, [0.0], color = (:gray, 0.4), linestyle = :dash)
        vlines!(ax, [0.0], color = (:black, 0.5), linestyle = :dash)

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, -0.05, 0.15)

        if i == 1
            axislegend(ax, position = :rt, framevisible = false)
        end
    end

    Label(
        fig[0, 1:3],
        "ROAMM — Standard vs rERP Decoding (N = $(length(subjects)))",
        fontsize = 22,
        font = :bold
    )

    return fig
end

function plot_plain_vs_two_step_b2b(
    scores_plain,
    scores_two_step
)
    targets = [
        "word_length_c",
        "frequency_zipf_c",
        "word_surprisal_c"
    ]

    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    # Only compare subjects present in both methods
    subjects = intersect(
        unique(scores_plain.subject),
        unique(scores_two_step.subject)
    )

    @info "Matched subjects" n = length(subjects)
    isempty(subjects) && error("No matching subjects found.")

    function summarize_matched(df)
        matched = filter(
            row -> row.subject in subjects &&
                   string(row.coefname) in targets,
            df
        )

        # One score per subject, predictor, time point
        subject_mean = combine(
            groupby(matched, [:subject, :coefname, :time]),
            :estimate => mean => :estimate
        )

        # Mean across subjects
        combine(
            groupby(subject_mean, [:coefname, :time]),
            :estimate => mean => :estimate_mean
        )
    end

    plain_mean = summarize_matched(scores_plain)
    two_step_mean = summarize_matched(scores_two_step)

    fig = Figure(size = (1500, 460))
    colgap!(fig.layout, 35)

    for (i, target) in enumerate(targets)
        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Absolute B2B coefficient" : "",
            xticks = -0.1:0.2:1.0,
            yticklabelsvisible = i == 1,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        plain = sort(
            plain_mean[string.(plain_mean.coefname) .== target, :],
            :time
        )

        two = sort(
            two_step_mean[string.(two_step_mean.coefname) .== target, :],
            :time
        )

        lines!(
            ax, plain.time, plain.estimate_mean;
            color = "#0072B2",
            linewidth = 3,
            label = "Plain B2B"
        )

        lines!(
            ax, two.time, two.estimate_mean;
            color = "#D55E00",
            linewidth = 3,
            label = "Two-step B2B"
        )

        vlines!(ax, [0.0], color = (:black, 0.5), linestyle = :dash)

        xlims!(ax, -0.1, 1.0)

        if i == 1
            axislegend(ax, position = :rt, framevisible = false)
        end
    end

    Label(
        fig[0, 1:3],
        "ROAMM — Plain vs Two-step B2B (N = $(length(subjects)))",
        fontsize = 22,
        font = :bold
    )

    return fig
end