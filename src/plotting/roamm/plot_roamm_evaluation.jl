function plot_group_b2b_delta(summary; title="ROAMM — B2B ΔR", n_subjects=nothing)

    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    fig = Figure(size=(1500, 460))
    colgap!(fig.layout, 35)

    for (i, target) in enumerate(targets)

        sm = sort(summary[summary.target .== target, :], :time)

        isempty(sm) && error("Missing target: $target")

        ax = Axis(
            fig[1, i],
            title=titles[i],
            xlabel="Time relative to fixation onset (s)",
            ylabel=i == 1 ? "ΔR" : "",
            xticks=-0.1:0.2:1.0,
            topspinevisible=false,
            rightspinevisible=false,
            xgridvisible=false,
            ygridvisible=false
        )

        lines!(
            ax,
            sm.time,
            sm.delta_mean;
            color="#0072B2",
            linewidth=3,
            label="20% winsorized mean"
        )

        lines!(
            ax,
            sm.time,
            sm.delta_median;
            color="#D55E00",
            linestyle=:dash,
            linewidth=2,
            label="Median"
        )

        vlines!(ax, [0.0], color=(:black, 0.5), linestyle=:dash)
        hlines!(ax, [0.0], color=(:gray, 0.5), linestyle=:dash)

        xlims!(ax, -0.1, 1.0)

        if i == 1
            axislegend(ax, position=:rt, framevisible=false)
        end
    end

    label = n_subjects === nothing ? title : "$title (N = $n_subjects)"

    Label(fig[0, 1:3], label, fontsize=22, font=:bold)

    return fig
end


function plot_plain_vs_two_step_delta(scores_plain, scores_two_step)

    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    # Match subjects
    subjects = intersect(
        string.(unique(scores_plain.subject)),
        string.(unique(scores_two_step.subject))
    )

    isempty(subjects) && error("No matching subjects found")

    # Retain shared subjects
    plain_matched = filter(
        row -> string(row.subject) in subjects,
        scores_plain
    )

    two_matched = filter(
        row -> string(row.subject) in subjects,
        scores_two_step
    )

    # Match exact subject, target and time observations
    plain_subject = combine(
        groupby(plain_matched, [:subject, :target, :time]),
        :delta_R => mean => :plain_delta
    )

    two_subject = combine(
        groupby(two_matched, [:subject, :target, :time]),
        :delta_R => mean => :two_delta
    )

    paired = innerjoin(
        plain_subject,
        two_subject;
        on=[:subject, :target, :time]
    )

    isempty(paired) && error("No matching subject/target/time observations")

    # Group summary for each pipeline, using matched observations
    plain_mean = combine(
        groupby(paired, [:target, :time]),
		:plain_delta =>  mean => :delta_mean
        #:plain_delta => (x -> winsorized_mean(x; proportion=0.2)) => :delta_mean
    )

    two_mean = combine(
        groupby(paired, [:target, :time]),
        :two_delta => mean => :delta_mean
		#:two_delta => (x -> winsorized_mean(x; proportion=0.2)) => :delta_mean
    )

    n_subjects = length(unique(paired.subject))

    fig = Figure(size=(1500, 460))
    colgap!(fig.layout, 35)
	axes = Axis[]

    for (i, target) in enumerate(targets)

        plain = sort(plain_mean[plain_mean.target .== target, :], :time)
        two = sort(two_mean[two_mean.target .== target, :], :time)

        isempty(plain) && error("Missing Plain B2B target: $target")
        isempty(two) && error("Missing Two-step B2B target: $target")

        ax = Axis(
            fig[1, i],
            title=titles[i],
            xlabel="Time relative to fixation onset (s)",
            ylabel=i == 1 ? " 20% winsorized mean ΔR" : "",
            xticks=-0.1:0.2:1.0,
            topspinevisible=false,
            rightspinevisible=false,
            xgridvisible=false,
            ygridvisible=false
        )
		push!(axes, ax)
		
        lines!(
            ax,
            plain.time,
            plain.delta_mean;
            color="#0072B2",
            linewidth=3,
            label="Plain B2B"
        )

        lines!(
            ax,
            two.time,
            two.delta_mean;
            color="#D55E00",
            linewidth=3,
            label="Two-step B2B"
        )

        vlines!(ax, [0.0], color=(:black, 0.5), linestyle=:dash)
        hlines!(ax, [0.0], color=(:gray, 0.5), linestyle=:dash)

        xlims!(ax, -0.1, 1.0)

        if i == 1
            axislegend(ax, position=:rt, framevisible=false)
			
        end
    end

	linkyaxes!(axes...)
    Label(
        fig[0, 1:3],
        "ROAMM — Plain vs Two-step B2B ΔR (N = $n_subjects)",
        fontsize=22,
        font=:bold
    )

    return fig
end

function plot_four_pipeline_R(
    scores_standard,
    scores_rerp,
    scores_plain,
    scores_two_step
)
    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    datasets = [
        scores_standard,
        scores_rerp,
        scores_plain,
        scores_two_step
    ]

    # Match subjects across all four pipelines
    subjects = reduce(intersect, [
        Set(string.(unique(df.subject))) for df in datasets
    ])

    isempty(subjects) && error("No matching subjects")

    n_subjects = length(subjects)

    # Group mean across subjects, no smoothing
    function summarize_matched(df)
        matched = filter(
            row -> string(row.subject) in subjects,
            df
        )

        # Make subject IDs consistent
        matched.subject = string.(matched.subject)

        subject_mean = combine(
            groupby(matched, [:subject, :target, :time]),
            :r => mean => :r
        )

        return combine(
            groupby(subject_mean, [:target, :time]),
            #:r => (x -> winsorized_mean(x; proportion=0.2)) => :winsor20
            :r => (x -> begin
                v = filter(isfinite, x)
                isempty(v) ? NaN : mean(v)
            end) => :mean_r
        )
    end

    summaries = [summarize_matched(df) for df in datasets]

    labels = [
        "Standard decoding (r)",
        "rERP decoding (r)",
        "Plain B2B (R)",
        "Two-step B2B (R)"
    ]

    colors = [
        "#0072B2",
        "#D55E00",
        "#009E73",
        "#CC79A7"
    ]

    fig = Figure(size=(1500, 480))
    colgap!(fig.layout, 35)

    axes = Axis[]
    legend_lines = []

    for (i, target) in enumerate(targets)

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "ΔR" : "",
            xticks = -0.1:0.2:1.0,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        push!(axes, ax)

        for j in 1:4
            sm = sort(
                summaries[j][summaries[j].target .== target, :],
                :time
            )

            line = lines!(
                ax,
                sm.time,
                sm.mean_r;
                color = colors[j],
                linewidth = 2.5
            )

            if i == 1
                push!(legend_lines, line)
            end
        end

        vlines!(
            ax, [0.0];
            color = (:black, 0.5),
            linestyle = :dash
        )

        hlines!(
            ax, [0.0];
            color = (:gray, 0.4),
            linestyle = :dash
        )

        xlims!(ax, -0.1, 1.0)
    end

    # Same y-axis scale across all predictors
    linkyaxes!(axes...)

    # Shared legend
    Legend(
        fig[2, 1:3],
        legend_lines,
        labels;
        orientation = :horizontal,
        framevisible = false
    )

    Label(
        fig[0, 1:3],
        "ROAMM — Four-pipeline comparison (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

    return fig
end