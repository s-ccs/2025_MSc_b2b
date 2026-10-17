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
using Serialization, Unfold, DataFrames

# ╔═╡ 8bdced5b-26e5-435c-b3dc-6710fb27c5be
using CairoMakie, Statistics

# ╔═╡ d63e9e0e-1a2a-4832-be56-1a5de5085976
md"""
# Step1: load all subject results"""

# ╔═╡ 486482b5-2120-41bf-b749-21f314180fc4
begin 
	function load_roamm_two_step_b2b_results(results_dir)
		files = filter(f -> endswith(f, ".jls"), readdir(results_dir; join = true))
		isempty(files) && error("No .jls files found in $results_dir")

		tables = DataFrame[]
		for file in files
			result = deserialize(file)
			coef_df = DataFrame(Unfold.coeftable(result.b2b_fit_corrected))
			coef_df[!, :subject] = fill(String(result.subject), nrow(coef_df))
			push!(tables, coef_df)
		end

		return vcat(tables...)
	end
	
end

# ╔═╡ e1890edc-0177-4a96-ab87-b0e392a07b22
scores_all = load_roamm_two_step_b2b_results("/home/st/st_us-051900/st_st181570/2025_MSc_b2b/slurm_results/models/roamm/two_step_roamm")

# ╔═╡ 0fcebb9e-ede3-4814-963e-103e32a2fb5c
md"""
# Step2: average across subjects
"""

# ╔═╡ 8d5901a1-72b7-4fba-83c1-9b938426ca54
function summarize_roamm_scores(scores_all)
	g = groupby(scores_all, [:coefname, :time])

	summary = combine(
		g,
		:estimate => mean => :estimate_mean,
		:estimate => std => :estimate_std,
		:estimate => median => :estimate_median,
		:estimate => length => :n
	)

	summary[!, :estimate_sem] = summary.estimate_std ./ sqrt.(summary.n)
	return summary
end

# ╔═╡ 725d8080-5f43-412a-8384-eee1d0dfb2ce
scores_mean = summarize_roamm_scores(scores_all)

# ╔═╡ af88e308-78e8-4fc3-9fdd-34a70b4d33db
# ╠═╡ disabled = true
#=╠═╡


function plot_roamm_b2b(scores_mean)

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

    fig = Figure(size = (1500, 450))

    for (i, target) in enumerate(targets)

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time [s]",
            ylabel = i == 1 ? "B2B coefficient" : ""
        )

        sub = sort(
            scores_mean[scores_mean.coefname .== target, :],
            :time
        )

        band!(
            ax,
            sub.time,
            sub.estimate_mean .- sub.estimate_sem,
            sub.estimate_mean .+ sub.estimate_sem,
            color = (:steelblue, 0.2)
        )

        lines!(
            ax,
            sub.time,
            sub.estimate_mean,
            color = :steelblue,
            linewidth = 2
        )

        hlines!(ax, [0.0], color = :gray, linestyle = :dash)

        xlims!(ax, -0.1, 1.0)
    end

    Label(
        fig[0, 1:3],
        "ROAMM — Two-step B2B",
        fontsize = 24,
        font = :bold
    )

    return fig
end
  ╠═╡ =#

# ╔═╡ ddd6f2d0-7630-4919-bc7e-d9bd25146c45
#=╠═╡
begin
	fig = plot_roamm_b2b(scores_mean)
fig
end
  ╠═╡ =#

# ╔═╡ 57b707fc-6089-485a-acfb-139ee57f52ce
md"""
# Plotting"""

# ╔═╡ e0b2923d-ee37-406e-92cd-5ce65c46b039
function plot_roamm_b2b_individual(scores_all, scores_mean)

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

    for (i, target) in enumerate(targets)

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "B2B coefficient" : "",
            xticks = -0.1:0.1:1.0,
            topspinevisible = false,
            rightspinevisible = false,
            xgridvisible = false,
            ygridvisible = false
        )

        sub = scores_all[scores_all.coefname .== target, :]

        # Individual subjects
        for subject in unique(sub.subject)
            ss = sub[sub.subject .== subject, :]
            sort!(ss, :time)

            lines!(
                ax,
                ss.time,
                ss.estimate,
                color = ("#85B8D8", 0.18),
                linewidth = 1
            )
        end

        # Group mean
        sm = scores_mean[scores_mean.coefname .== target, :]
        sort!(sm, :time)

        lines!(
            ax,
            sm.time,
            sm.estimate_mean,
            color = "#155A91",
            linewidth = 3,
            label = "Mean"
        )
        # median
        lines!(ax, sm.time, sm.estimate_median, color = :darkorange, linewidth = 2, linestyle = :dash, label = "Median")
        
        # Fixation onset
        vlines!(
            ax, [0.0],
            color = (:black, 0.55),
            linestyle = :dash,
            linewidth = 1
        )

        # Zero baseline
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
        "ROAMM — Two-step B2B (N = $n_subjects)",
        fontsize = 22,
        font = :bold
    )

       

    return fig
end

# ╔═╡ 9d08c2d9-cfe2-4be2-afc3-75d88b5ed4b8
begin
	fig_individual = plot_roamm_b2b_individual(scores_all, scores_mean)
fig_individual
end

# ╔═╡ d7e216c3-d3f9-430c-8304-fc82c66b9653
md"""
# Heatmap
"""

# ╔═╡ 1772e0ce-78ea-416b-b1ab-a0a38eb06696

function plot_roamm_b2b_heatmap(scores_all)

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

    for (i, target) in enumerate(targets)

        sub = scores_all[scores_all.coefname .== target, :]
        times = sort(unique(sub.time))

        # Matrix: time × subjects
        mat = Matrix{Float64}(undef, length(times), length(subjects))

        for (j, subject) in enumerate(subjects)
            ss = sort(sub[sub.subject .== subject, :], :time)

            # Check that subjects have the same time points
            @assert ss.time == times "Time mismatch for $subject, $target"

            mat[:, j] = ss.estimate
        end

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time (s)",
            ylabel = i == 1 ? "Subject index" : "",
            xticks = -0.1:0.1:1.0,
            yticks = 1:5:length(subjects),
            topspinevisible = false,
            rightspinevisible = false
        )

        hm = heatmap!(
            ax,
            times,
            1:length(subjects),
            mat,
            colormap = :viridis
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
            label = "B2B coefficient"
        )
    end

    Label(
        fig[0, 1:3],
        "ROAMM — Two-step B2B — subject-level heatmaps",
        fontsize = 22,
        font = :bold
    )

    return fig
end


# ╔═╡ b320e128-2a0a-4db4-83f3-794b7c2cf4c7
begin
	fig_heatmap = plot_roamm_b2b_heatmap(scores_all)
fig_heatmap
end

# ╔═╡ 8cfdac79-78da-40c9-9e5d-4ae43bb0cb18
md"""
# Are there negative B2B coefficients?
"""

# ╔═╡ ba0d9fa9-036c-4240-9c99-c2d7cae4a245

combine(
    groupby(scores_all, :coefname),
    :estimate => minimum => :min,
    :estimate => maximum => :max,
    :estimate => mean => :mean,
    :estimate => median => :median,
    :estimate => (x -> mean(x .< 0)) => :fraction_negative,
    :estimate => (x -> mean(x .== 0)) => :fraction_zero
)


# ╔═╡ 944b1cad-843d-46a7-8d9a-956eaa557a65
# ╠═╡ disabled = true
#=╠═╡
function plot_roamm_standard_decoding(
	scores::AbstractDataFrame;
	title::String = "ROAMM standard decoding",
	x_window = (-0.1, 1.0),
	metric::Symbol = :r_mean,
	show_band::Bool = true
)
	for col in (:time, :target, metric)
		col in propertynames(scores) || error("Missing column: $col")
	end

	targets = ["word_length", "frequency_zipf", "word_surprisal"]
	available_targets = unique(string.(scores.target))
	targets = [t for t in targets if t in available_targets]

	isempty(targets) && error("No expected targets found.")

	fig = Figure(size = (1500, 450))
	axes = Axis[]

	# shared y range
	yvals = Float64.(scores[!, metric])
	yvals = yvals[isfinite.(yvals)]

	y_min = min(0.0, minimum(yvals))
	y_max = maximum(yvals)
	y_range = max(y_max - y_min, 0.1)
	y_upper = y_max + 0.1 * y_range

	target_titles = Dict(
		"word_length" => "Word length",
		"frequency_zipf" => "Frequency (Zipf)",
		"word_surprisal" => "Word surprisal"
	)

	for (i, target_name) in enumerate(targets)
		ax = Axis(
			fig[1, i],
			title = target_titles[target_name],
			xlabel = "Time [s]",
			ylabel = i == 1 ? "Correlation (r)" : "",
			topspinevisible = false,
			rightspinevisible = false,
			xgridvisible = false,
			ygridvisible = false
		)
		sub = DataFrame(scores[string.(scores.target) .== target_name, :])
		sort!(sub, :time)

		lines!(ax, sub.time, sub[!, metric])

		if show_band && (:r_sem in propertynames(sub))
			band!(
				ax,
				sub.time,
				sub[!, metric] .- sub.r_sem,
				sub[!, metric] .+ sub.r_sem,
				color = (:steelblue, 0.2)
			)
			lines!(ax, sub.time, sub[!, metric])
		end
		xlims!(ax, x_window...)
        ylims!(ax, y_min, y_upper)
		push!(axes, ax)
	end
	linkxaxes!(axes...)
	linkyaxes!(axes...)
	Label(
        fig[0, 1:length(targets)],
        title,
        fontsize = 24,
        font = :bold,
        padding = (0, 0, 10, 0)
    )

    return fig
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═e9ce96e4-ca36-11f1-84a8-df0069b540ab
# ╠═4001b213-854f-4305-8238-7a1742945c9a
# ╠═d63e9e0e-1a2a-4832-be56-1a5de5085976
# ╠═486482b5-2120-41bf-b749-21f314180fc4
# ╠═e1890edc-0177-4a96-ab87-b0e392a07b22
# ╠═0fcebb9e-ede3-4814-963e-103e32a2fb5c
# ╠═8bdced5b-26e5-435c-b3dc-6710fb27c5be
# ╠═8d5901a1-72b7-4fba-83c1-9b938426ca54
# ╠═725d8080-5f43-412a-8384-eee1d0dfb2ce
# ╠═af88e308-78e8-4fc3-9fdd-34a70b4d33db
# ╠═ddd6f2d0-7630-4919-bc7e-d9bd25146c45
# ╠═57b707fc-6089-485a-acfb-139ee57f52ce
# ╠═e0b2923d-ee37-406e-92cd-5ce65c46b039
# ╠═9d08c2d9-cfe2-4be2-afc3-75d88b5ed4b8
# ╠═d7e216c3-d3f9-430c-8304-fc82c66b9653
# ╠═1772e0ce-78ea-416b-b1ab-a0a38eb06696
# ╠═b320e128-2a0a-4db4-83f3-794b7c2cf4c7
# ╟─8cfdac79-78da-40c9-9e5d-4ae43bb0cb18
# ╠═ba0d9fa9-036c-4240-9c99-c2d7cae4a245
# ╠═944b1cad-843d-46a7-8d9a-956eaa557a65
