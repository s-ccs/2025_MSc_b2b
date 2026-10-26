### A Pluto.jl notebook ###
# v1.0.4

using Markdown
using InteractiveUtils

# ╔═╡ e9ce96e4-ca36-11f1-84a8-df0069b540ab
begin
    import Pkg
	Pkg.activate(joinpath(@__DIR__, "..", ".."))
end

# ╔═╡ 4001b213-854f-4305-8238-7a1742945c9a
begin
	using Serialization, Unfold, DataFrames
	using CairoMakie, Statistics, CSV
end

# ╔═╡ 4f8ae38b-5aff-4ab5-b546-9b1e88c1dc85
import MScB2B

# ╔═╡ d63e9e0e-1a2a-4832-be56-1a5de5085976
md"""
# Step1: load all subject results"""

# ╔═╡ 74612582-bead-443b-95b0-d123fbb2bb62
begin
	base_dir = "/home/st/st_us-051900/st_st181570/2025_MSc_b2b/slurm_results/models/roamm"
	#base_dir = "/store/users/xu/2025_MSc_b2b/slurm_results/models/roamm"
	
	readdir(base_dir)
end

# ╔═╡ 28c97b20-219b-4913-b147-e4e422a0bca0
md"""
# Helper functions
"""

# ╔═╡ 75e697c5-db1c-4d53-a95d-9f5615ffd99a
function load_b2b_results(results_dir, model_field::Symbol)
	files = filter(f -> endswith(f, ".jls"), readdir(results_dir; join = true))
	isempty(files) && error("No .jls files found in $results_dir")

	tables = DataFrame[]
	for file in files
		result = deserialize(file)
		df = DataFrame(Unfold.coeftable(getproperty(result, model_field)))
		df.subject = fill(String(result.subject), nrow(df))
		df.estimate_signed = copy(df.estimate)
		df.estimate = abs.(df.estimate_signed)
		filter!(:coefname => !=("(Intercept)"), df)
		push!(tables, df)
	end

	return vcat(tables...)
end

# ╔═╡ 01fcdcce-3ece-45b8-9e25-7101fa06c6bc
function load_roamm_decoding_scores(results_dir)
	files = filter(f -> endswith(f, ".jls"), readdir(results_dir; join = true))
	isempty(files) && error("No .jls files found in $results_dir")
	tables = DataFrame[]
	for file in files
		result = deserialize(file)
		df = DataFrame(result.scores)
		df[!, :subject] = fill(String(result.subject), nrow(df))
		push!(tables, df)
	end

	return vcat(tables...)
end

# ╔═╡ 57ec5561-d008-407f-99b7-846ef9b7c265
function winsorized_mean(x; proportion=0.1)
	v = sort(collect(x))
	k = floor(Int, proportion * length(v))
	k == 0 && return mean(v)
	v[1:k] .= v[k+1]
	v[end-k+1:end] .= v[end-k]
	return mean(v)
end

# ╔═╡ 1f536d81-ad42-4fa7-b5e3-9cc5d430cabd
function summarize_b2b(scores)
	summary = combine(
		groupby(scores, [:coefname, :time]),
		:estimate => mean => :estimate_mean,
		:estimate => median => :estimate_median,
		:estimate => std => :estimate_std,
		:estimate => length => :n,
		:estimate => (x -> winsorized_mean(x; proportion=0.1)) => :winsor10,
		:estimate => (x -> winsorized_mean(x; proportion=0.2)) => :winsor20	
	)
	summary.estimate_sem = summary.estimate_std ./ sqrt.(summary.n)
	return summary
end

# ╔═╡ 58460a3e-9c05-46b4-a402-f7f65ed79b15
function summarize_decoding(scores)
	# One value per participant, predictor and time point
	subject_scores = combine(groupby(scores, [:subject, :target, :time]), :r => mean => :r)
	summary = combine(
		groupby(subject_scores, [:target, :time]),
		:r => mean => :r_mean,
		:r => median => :r_median,
		:r => std => :r_std,
		:r => length => :n,
		:r => (x -> winsorized_mean(x; proportion=0.1)) => :winsor10,
		:r => (x -> winsorized_mean(x; proportion=0.2)) => :winsor20	
	)
	summary.r_sem = summary.r_std ./ sqrt.(summary.n)
	return summary
end

# ╔═╡ 57b707fc-6089-485a-acfb-139ee57f52ce
md"""
# Plotting"""

# ╔═╡ 30d5aa98-a963-44ab-bded-24058687f5df
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

# ╔═╡ fadd9408-edeb-4b56-bf33-95238fc4adfb
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

# ╔═╡ 9b318652-6a11-417a-a795-f2c119562e57
function plot_roamm_b2b(
    scores_all,
    scores_mean;
    title = "ROAMM — B2B"
)

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

# ╔═╡ 2624f753-ee4e-4c04-867f-7e1028045d88
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

# ╔═╡ ac49403c-97b7-4628-a4bd-1e57ae8b1863
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

# ╔═╡ c8088f5b-18b7-448e-aa08-c8cfb0b34066
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

# ╔═╡ e9bf43b8-fb22-4979-ae96-c40717e43f9e
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

# ╔═╡ 4685dc86-21ec-4a03-ac5f-41f7ba9d46b7
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

# ╔═╡ 4ea40ee8-fdd7-47c8-8dc5-24d21c6f86d9
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

# ╔═╡ 3503d980-be17-42ef-88cd-eefc4ff429c6
md"""
# experimental
"""

# ╔═╡ b12bdc36-92da-420a-a0ed-025f24ae3b3f
# ╠═╡ disabled = true
#=╠═╡
function plot_roamm_b2b_heatmap(
    scores_all;
    title = "ROAMM — B2B — subject-level heatmaps",
    colorranges = nothing
)
    targets = [
        "word_length",
        "frequency_zipf",
        "word_surprisal"
    ]

    # Handle centered predictors in two-step B2B
    if "word_length_c" in scores_all.coefname
        targets = targets .* "_c"
    end

    titles = [
        "Word length",
        "Zipf frequency",
        "Word surprisal"
    ]

    fig = Figure(size = (1500, 650))
    colgap!(fig.layout, 30)

    subjects = sort(unique(scores_all.subject))

    for (i, target) in enumerate(targets)

        sub = scores_all[scores_all.coefname .== target, :]
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

        # Linear color scale
        crange = if colorranges === nothing
            (0.0, maximum(mat))
        else
            (0.0, colorranges[i])
        end

        hm = heatmap!(
            ax,
            times,
            1:length(subjects),
            mat,
            colormap = :viridis,
            colorrange = crange
        )

        # Fixation onset
        vlines!(
            ax,
            [0.0],
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
  ╠═╡ =#

# ╔═╡ 3e044ac4-b157-4a0a-97d7-8de8049e6c5f
# ╠═╡ disabled = true
#=╠═╡
function plot_roamm_decoding(
    scores_all,
    scores_summary;
    title = "ROAMM — decoding",
    ylim = (-0.2, 0.2)
)
    targets = ["word_length", "frequency_zipf", "word_surprisal"]
    titles = ["Word length", "Zipf frequency", "Word surprisal"]

    n_subjects = length(unique(scores_all.subject))

    # One curve per subject, target, and time point
    subject_scores = combine(
        groupby(scores_all, [:subject, :target, :time]),
        :r => mean => :r
    )

    fig = Figure(size = (1500, 460))
    colgap!(fig.layout, 40)

    for (i, target) in enumerate(targets)

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Decoding correlation (r)" : "",
            xticks = -0.1:0.1:1.0,
            yticks = -0.2:0.1:0.2,
            yticklabelsvisible = i == 1,
            yticksvisible = i == 1,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        sub = subject_scores[subject_scores.target .== target, :]

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

        # Group mean
        lines!(
            ax,
            sm.time,
            sm.r_mean,
            color = "#155A91",
            linewidth = 3
        )

        # 10% winsorized mean
        lines!(
            ax,
            sm.time,
            sm.winsor10,
            color = "#D55E00",
            linewidth = 2,
            linestyle = :dash
        )

        # Fixation onset
        vlines!(
            ax, [0.0],
            color = (:black, 0.55),
            linestyle = :dash,
            linewidth = 1
        )

        # Zero correlation
        hlines!(
            ax, [0.0],
            color = (:gray, 0.5),
            linestyle = :dash,
            linewidth = 1
        )

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, ylim...)
    end

    Label(
        fig[0, 1:3],
        "$title (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

    return fig
end
  ╠═╡ =#

# ╔═╡ 1f9045d8-d9c7-41ef-aa56-034804f6ce23
# ╠═╡ disabled = true
#=╠═╡
function plot_roamm_b2b(
    scores_all,
    scores_mean;
    title = "ROAMM — B2B",
    ylims_by_target = nothing
)
    targets = [
        "word_length",
        "frequency_zipf",
        "word_surprisal"
    ]

    # Handle centered predictors in two-step B2B
    if "word_length_c" in scores_all.coefname
        targets = targets .* "_c"
    end

    titles = [
        "Word length",
        "Zipf frequency",
        "Word surprisal"
    ]

    n_subjects = length(unique(scores_all.subject))

    fig = Figure(size = (1500, 460))
    colgap!(fig.layout, 40)

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
            ylabel = i == 1 ? "Absolute B2B coefficient" : "",
            xticks = -0.1:0.1:1.0,

            # Linear scale
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

        # Group mean
        lines!(
            ax,
            sm.time,
            sm.estimate_mean,
            color = "#155A91",
            linewidth = 3
        )

        # Group median
        lines!(
            ax,
            sm.time,
            sm.estimate_median,
            color = :darkorange,
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

        xlims!(ax, -0.1, 1.0)

        # Optional fixed limits for each predictor
        if ylims_by_target !== nothing
            ylims!(ax, 0, ylims_by_target[i])
        else
            ylims!(ax, low = 0)
        end
    end

    Label(
        fig[0, 1:3],
        "$title (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

    return fig
end
  ╠═╡ =#

# ╔═╡ 37593f70-0cda-4901-bbb2-a69ebc384f3f
# ╠═╡ disabled = true
#=╠═╡
function plot_roamm_decoding_heatmap(
    scores_all;
    title = "ROAMM — decoding — subject-level heatmaps"
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

    fig = Figure(size = (1500, 650))
    colgap!(fig.layout, 30)

    subjects = sort(unique(scores_all.subject))

    # ---- global symmetric color scale across all 3 predictors ----
    all_r = Float64.(scores_all.r)
    all_r = all_r[isfinite.(all_r)]
    rmax = maximum(abs.(all_r))
    crange = (-rmax, rmax)

    hms = Any[]

    for (i, target) in enumerate(targets)

        sub = scores_all[string.(scores_all.target) .== target, :]
        times = sort(unique(sub.time))

        # Matrix: time × subjects
        mat = Matrix{Float64}(undef, length(times), length(subjects))

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

        hm = heatmap!(
            ax,
            times,
            1:length(subjects),
            mat,
            colormap = :balance,
            colorrange = crange
        )

        push!(hms, hm)

        vlines!(
            ax, [0.0],
            color = :white,
            linestyle = :dash,
            linewidth = 1.5
        )

        hlines!(
            ax, [0.0],
            color = (:white, 0.0)  # do nothing, just keep structure optional
        )

        xlims!(ax, -0.1, 1.0)
        ylims!(ax, 0.5, length(subjects) + 0.5)
    end

    # one shared colorbar
    Colorbar(
        fig[2, 1:3],
        hms[1],
        vertical = false,
        label = "Decoding correlation (r)"
    )

    Label(
        fig[0, 1:3],
        title,
        fontsize = 22,
        font = :bold
    )

    return fig
end
  ╠═╡ =#

# ╔═╡ 3eacb1ef-27af-40bf-a738-b4aad27fa9b9
md"""
# Plot Standard ridge decoding
"""

# ╔═╡ 326c7e77-55e6-48a4-b834-31b4a635ec12
begin
	scores_standard_decoding = load_roamm_decoding_scores(joinpath(base_dir, "standard_centering_foldmean"))
	summary_standard_decoding = summarize_decoding(scores_standard_decoding)
	fig_standard = plot_roamm_decoding(
		scores_standard_decoding,
		summary_standard_decoding;
		title = "ROAMM - Standard decoding"
	)
	fig_standard
end

# ╔═╡ f68010fb-0e8f-442f-96d4-e1b37903be28
propertynames(scores_standard_decoding)

# ╔═╡ c82a398e-8581-455c-b97d-7884f7db39aa
# ╠═╡ disabled = true
#=╠═╡
begin
	fig_standard_log = plot_roamm_decoding_log(
		scores_standard_decoding,
		summary_standard_decoding;
		#title = "ROAMM - Standard decoding"
	)
	fig_standard_log
end
  ╠═╡ =#

# ╔═╡ 54a711cb-5284-4a0a-8314-3e6cf9dc2bc7
# ╠═╡ disabled = true
#=╠═╡
begin
	fig_standard_mean = plot_roamm_decoding_mean(
		summary_standard_decoding;
		title = "ROAMM - Standard decoding"
	)
	fig_standard_mean
end
  ╠═╡ =#

# ╔═╡ d6e439f9-b27d-4bcd-ae99-e06b061f1264
begin
    fig_standard_decoding_heatmap = plot_roamm_decoding_heatmap(
        scores_standard_decoding;
        title = "ROAMM — standard decoding — subject-level heatmaps"
    )

    fig_standard_decoding_heatmap
end

# ╔═╡ 9236fe57-780a-4028-bef6-a90399ef492b
md"""
# Plot rerp"""

# ╔═╡ 8f43b5ca-8c27-4180-9fbd-975fd881c762
begin
    scores_rerp = load_roamm_decoding_scores(
        joinpath("/home/st/st_us-051900/st_st181570/2025_MSc_b2b/slurm_results/models/roamm/rerp")
    )

    summary_rerp = summarize_decoding(scores_rerp)

    fig_rerp = plot_roamm_decoding(
        scores_rerp,
        summary_rerp;
        title = "ROAMM — rERP decoding"
    )

    fig_rerp
end

# ╔═╡ d8c8e4f3-833c-44cf-9c48-814d42fd87a6
begin
	names(scores_rerp)
	first(scores_rerp, 5)
end

# ╔═╡ 6c9341b3-076c-41a7-99c3-c154aa13e8f5
# ╠═╡ disabled = true
#=╠═╡
begin
    fig_comparison = plot_standard_vs_rerp(
        scores_standard_decoding,
        scores_rerp
    )
    
    fig_comparison
end
  ╠═╡ =#

# ╔═╡ 62b100f0-f9d2-4a2a-964c-86179e9546e6
names(scores_rerp)

# ╔═╡ 0e68529d-5fc6-4a16-b151-8558fe9d52bb
length(unique(scores_rerp.subject))

# ╔═╡ 6c5b118d-5f09-40ed-9805-d49085dfa0f1
begin
    fig_comparison = plot_standard_vs_rerp(
        scores_standard_decoding,
        scores_rerp
    )
    
    fig_comparison
end

# ╔═╡ ffcacd1d-38b4-4f34-8f9b-dd90606d7f85
begin
    @show extrema(scores_standard_decoding.time)
    @show extrema(scores_rerp.time)

    @show length(unique(scores_standard_decoding.time))
    @show length(unique(scores_rerp.time))
end

# ╔═╡ 07b34fdd-1e7d-4307-aea3-0098d2f75ce3
begin
    fig_rERP_decoding_heatmap = plot_roamm_decoding_heatmap(
        scores_rerp;
        title = "ROAMM — rERP decoding — subject-level heatmaps"
    )

    fig_rERP_decoding_heatmap
end

# ╔═╡ c3cec440-5592-4f29-8dcc-c3ce58b3a1d0
md"""
# Standard decoding check
"""


# ╔═╡ d3503f74-b9d5-4ce3-b553-fe647d3e3921


# ╔═╡ 130e150c-a908-4d70-9b1d-356a6009416b
md"""
# Plot plain_b2b
"""

# ╔═╡ 7d93eff2-0bec-4423-8db1-0439527dd2fa
begin
	scores_plain_b2b = load_b2b_results(joinpath(base_dir, "plain_b2b"), :b2b_fit_plain)
	summary_plain_b2b = summarize_b2b(scores_plain_b2b)
	fig_plain_b2b = plot_roamm_b2b(scores_plain_b2b, summary_plain_b2b; title = "ROAMM - Plain B2B")
	fig_plain_b2b
end

# ╔═╡ 0d646518-4b97-42c6-ac1e-b7c8f9a36df1
begin
    fig_plain_heatmap = plot_roamm_b2b_heatmap(
        scores_plain_b2b;
        title = "ROAMM — Plain B2B — subject-level heatmaps"
    )

    fig_plain_heatmap
end

# ╔═╡ d05990b6-7035-478a-889a-73dd8be464d9
md"""
# Plot two-step b2b
"""

# ╔═╡ 65bd6d0c-7db3-4087-ae5c-33433b67d613
begin
	scores_two_step = load_b2b_results(joinpath(base_dir, "two_step_roamm"), :b2b_fit_corrected)
	summary_two_step = summarize_b2b(scores_two_step)
	fig_two_step = plot_roamm_b2b(scores_two_step, summary_two_step; title = "ROAMM - Two-step B2B")
	
	#save("roamm_two_step_final.svg", fig_two_step)
	fig_two_step
	
end

# ╔═╡ f619cfe4-6758-4535-9056-7debf7739312
begin 
    fig = plot_plain_vs_two_step_b2b(
    scores_plain_b2b,
    scores_two_step
)

fig
#save("roamm_plain_vs_two_step_b2b.svg", fig)
end

# ╔═╡ a7da75d3-85af-4c50-a580-a28cdd9b9d6a
begin
	fig_two_step_mean = plot_roamm_b2b_mean(summary_two_step; title = "ROAMM - two_step B2B")
	fig_two_step_mean
end

# ╔═╡ 7b1ddef0-eacf-4b04-9bbe-c08151b7c1c6
begin
    fig_two_step_heatmap = plot_roamm_b2b_heatmap(
        scores_two_step;
        title = "ROAMM — Two-step B2B — subject-level heatmaps"
    )

    fig_two_step_heatmap
end

# ╔═╡ c17bfc57-2d51-493d-9da8-a72b04e5f3bb
# ╠═╡ disabled = true
#=╠═╡
begin
	save("01_ROAMM_standard_decoding.svg", fig_standard)
	save("01_ROAMM_standard_decoding_heatmap.svg", fig_standard_decoding_heatmap)
	save("02_ROAMM_rerp_decoding.svg", fig_rerp)
	save("standard_and_rerp_comparison.svg", fig_comparison)
	save("02_ROAMM_rerp_decoding_heatmap.svg", fig_rERP_decoding_heatmap)
	save("03_ROAMM_plain_b2b.svg", fig_plain_b2b)
	save("03_ROAMM_plain_b2b_heatmap.svg", fig_plain_heatmap)
	save("05_ROAMM_two_step_b2b.svg", fig_two_step)
	save("05_ROAMM_two_step_b2b_heatmap.svg", fig_two_step_heatmap)
	
end
  ╠═╡ =#

# ╔═╡ 7770b6ef-7c84-4d81-a4e9-f35a1b6b0d11
md"""
# Plot one-step b2b
"""

# ╔═╡ 8cfdac79-78da-40c9-9e5d-4ae43bb0cb18
md"""
# Are there negative B2B coefficients?
"""

# ╔═╡ f8aa3617-9110-4d71-bbeb-ee58101f73bf
begin
    targets = ["word_length", "frequency_zipf", "word_surprisal"]

    for target in targets

        standard = sort(
            summary_standard_decoding[
                summary_standard_decoding.target .== target, :
            ],
            :time
        )

        # Try centered or uncentered predictor name
        b2b_target = if target * "_c" in summary_plain_b2b.coefname
            target * "_c"
        else
            target
        end

        b2b = sort(
            summary_plain_b2b[
                summary_plain_b2b.coefname .== b2b_target, :
            ],
            :time
        )

        println("$target: standard=$(nrow(standard)), B2B=$(nrow(b2b))")

        if isempty(standard) || isempty(b2b)
            continue
        end

        @assert length(standard.time) == length(b2b.time)
        @assert all(isapprox.(standard.time, b2b.time; atol=1e-6))

        r = cor(standard.r_mean, b2b.estimate_mean)

        println("Temporal correlation = $(round(r, digits=3))")
    end
end

# ╔═╡ 28f461e7-4e8c-472a-accc-1afb6cf8d95a
begin
	early_b2b = filter(
    row -> 0.1 <= row.time <= 0.35,
    scores_two_step
)

subject_effects = combine(
    groupby(early_b2b, [:subject, :coefname]),
    :estimate => mean => :early_mean,
    :estimate => maximum => :early_peak
)

sort(subject_effects, [:coefname, :early_mean])
end

# ╔═╡ Cell order:
# ╠═e9ce96e4-ca36-11f1-84a8-df0069b540ab
# ╠═4f8ae38b-5aff-4ab5-b546-9b1e88c1dc85
# ╠═4001b213-854f-4305-8238-7a1742945c9a
# ╟─d63e9e0e-1a2a-4832-be56-1a5de5085976
# ╟─74612582-bead-443b-95b0-d123fbb2bb62
# ╟─f68010fb-0e8f-442f-96d4-e1b37903be28
# ╟─28c97b20-219b-4913-b147-e4e422a0bca0
# ╟─75e697c5-db1c-4d53-a95d-9f5615ffd99a
# ╟─01fcdcce-3ece-45b8-9e25-7101fa06c6bc
# ╟─57ec5561-d008-407f-99b7-846ef9b7c265
# ╟─1f536d81-ad42-4fa7-b5e3-9cc5d430cabd
# ╟─58460a3e-9c05-46b4-a402-f7f65ed79b15
# ╟─57b707fc-6089-485a-acfb-139ee57f52ce
# ╟─30d5aa98-a963-44ab-bded-24058687f5df
# ╟─fadd9408-edeb-4b56-bf33-95238fc4adfb
# ╟─9b318652-6a11-417a-a795-f2c119562e57
# ╟─2624f753-ee4e-4c04-867f-7e1028045d88
# ╟─ac49403c-97b7-4628-a4bd-1e57ae8b1863
# ╟─c8088f5b-18b7-448e-aa08-c8cfb0b34066
# ╟─e9bf43b8-fb22-4979-ae96-c40717e43f9e
# ╟─4685dc86-21ec-4a03-ac5f-41f7ba9d46b7
# ╟─4ea40ee8-fdd7-47c8-8dc5-24d21c6f86d9
# ╠═f619cfe4-6758-4535-9056-7debf7739312
# ╟─3503d980-be17-42ef-88cd-eefc4ff429c6
# ╟─b12bdc36-92da-420a-a0ed-025f24ae3b3f
# ╟─3e044ac4-b157-4a0a-97d7-8de8049e6c5f
# ╟─1f9045d8-d9c7-41ef-aa56-034804f6ce23
# ╟─37593f70-0cda-4901-bbb2-a69ebc384f3f
# ╟─3eacb1ef-27af-40bf-a738-b4aad27fa9b9
# ╠═326c7e77-55e6-48a4-b834-31b4a635ec12
# ╠═c82a398e-8581-455c-b97d-7884f7db39aa
# ╠═54a711cb-5284-4a0a-8314-3e6cf9dc2bc7
# ╠═d6e439f9-b27d-4bcd-ae99-e06b061f1264
# ╠═d8c8e4f3-833c-44cf-9c48-814d42fd87a6
# ╟─9236fe57-780a-4028-bef6-a90399ef492b
# ╠═8f43b5ca-8c27-4180-9fbd-975fd881c762
# ╠═6c9341b3-076c-41a7-99c3-c154aa13e8f5
# ╠═62b100f0-f9d2-4a2a-964c-86179e9546e6
# ╠═0e68529d-5fc6-4a16-b151-8558fe9d52bb
# ╠═6c5b118d-5f09-40ed-9805-d49085dfa0f1
# ╠═ffcacd1d-38b4-4f34-8f9b-dd90606d7f85
# ╠═07b34fdd-1e7d-4307-aea3-0098d2f75ce3
# ╟─c3cec440-5592-4f29-8dcc-c3ce58b3a1d0
# ╠═d3503f74-b9d5-4ce3-b553-fe647d3e3921
# ╟─130e150c-a908-4d70-9b1d-356a6009416b
# ╠═7d93eff2-0bec-4423-8db1-0439527dd2fa
# ╠═0d646518-4b97-42c6-ac1e-b7c8f9a36df1
# ╟─d05990b6-7035-478a-889a-73dd8be464d9
# ╠═65bd6d0c-7db3-4087-ae5c-33433b67d613
# ╠═a7da75d3-85af-4c50-a580-a28cdd9b9d6a
# ╠═7b1ddef0-eacf-4b04-9bbe-c08151b7c1c6
# ╠═c17bfc57-2d51-493d-9da8-a72b04e5f3bb
# ╟─7770b6ef-7c84-4d81-a4e9-f35a1b6b0d11
# ╠═8cfdac79-78da-40c9-9e5d-4ae43bb0cb18
# ╠═f8aa3617-9110-4d71-bbeb-ee58101f73bf
# ╠═28f461e7-4e8c-472a-accc-1afb6cf8d95a
