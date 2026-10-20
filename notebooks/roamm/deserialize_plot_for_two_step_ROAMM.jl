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
	using CairoMakie, Statistics
end

# ╔═╡ d63e9e0e-1a2a-4832-be56-1a5de5085976
md"""
# Step1: load all subject results"""

# ╔═╡ 74612582-bead-443b-95b0-d123fbb2bb62
begin
	base_dir = "/home/st/st_us-051900/st_st181570/2025_MSc_b2b/slurm_results/models/roamm"
	
	readdir(base_dir)
end

# ╔═╡ f68010fb-0e8f-442f-96d4-e1b37903be28
# ╠═╡ disabled = true
#=╠═╡
propertynames(scores_rerp)
  ╠═╡ =#

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
		df.coefname = replace.(string.(df.coefname), r"_c$" => "")
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

# ╔═╡ 7dc15294-de07-48d0-89eb-6b618a7663fe
function plot_roamm_decoding(
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

# ╔═╡ e0b2923d-ee37-406e-92cd-5ce65c46b039
function plot_roamm_b2b(scores_all, scores_mean; title = "ROAMM - B2B")

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

    for (i, target) in enumerate(targets)

        ax = Axis(
            fig[1, i],
            title = titles[i],
            xlabel = "Time relative to fixation onset (s)",
            ylabel = i == 1 ? "Absolute B2B coefficient (log_scale)" : "",
            yscale = log10,
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
            sm.winsor10,
            color = "#155A91",
            linewidth = 3,
            label = "Winsorized mean (10%)"
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
        #hlines!(
         #   ax, [0.0],
          #  color = (:gray, 0.5),
          #  linestyle = :dash,
           # linewidth = 1
        #)

        #xlims!(ax, -0.1, 1.0)
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
    subjects = sort(unique(scores_all.subject))

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

# ╔═╡ 3eacb1ef-27af-40bf-a738-b4aad27fa9b9
md"""
# Plot Standard ridge decoding
"""

# ╔═╡ 9236fe57-780a-4028-bef6-a90399ef492b
md"""
# Plot rerp"""

# ╔═╡ 30d5aa98-a963-44ab-bded-24058687f5df
function plot_roamm_decoding_individual(scores_all, scores_summary;
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

# ╔═╡ e8c2b476-b90c-4939-8360-9e871aac7118
unique(scores_plain_b2b.coefname)

# ╔═╡ 44cdb525-2581-43f5-9137-127fb2d8d2ef
combine(
    groupby(scores_plain_b2b, :coefname),
    :estimate => minimum => :min,
    :estimate => maximum => :max,
    :estimate => (x -> count(isfinite, x)) => :n_finite
)

# ╔═╡ 0d646518-4b97-42c6-ac1e-b7c8f9a36df1
begin
    fig_plain_heatmap = plot_roamm_b2b_heatmap(
        scores_plain_b2b;
        title = "ROAMM — Plain B2B — subject-level heatmaps"
    )

    fig_plain_heatmap
end

# ╔═╡ 6b99c440-4177-45ea-8c12-3c1faa0db6fd
begin
	save("final_roamm_plain_b2b.svg", fig_plain_b2b)
	save("final_roamm_plain_b2b_heatmap.svg", fig_plain_heatmap)
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

# ╔═╡ 6ff590fa-5e6b-4cd3-bcb0-063c9ee4cc95
unique(scores_two_step.coefname)

# ╔═╡ 7b1ddef0-eacf-4b04-9bbe-c08151b7c1c6
begin
    fig_two_step_heatmap = plot_roamm_b2b_heatmap(
        scores_two_step;
        title = "ROAMM — Two-step B2B — subject-level heatmaps"
    )

    fig_two_step_heatmap
end

# ╔═╡ c520c71d-9d96-4883-93ae-a691d964970b
begin
	save("final_roamm_two_step_b2b.svg", fig_two_step)
	save("final_roamm_two_step_b2b_heatmap.svg", fig_two_step_heatmap)
end

# ╔═╡ 7770b6ef-7c84-4d81-a4e9-f35a1b6b0d11
md"""
# Plot one-step b2b
"""

# ╔═╡ 8cfdac79-78da-40c9-9e5d-4ae43bb0cb18
md"""
# Are there negative B2B coefficients?
"""

# ╔═╡ 4ce4326c-69cf-4b51-82a3-8b2e627c78a2
# ╠═╡ disabled = true
#=╠═╡
begin
	scores_rerp_subject = combine(
    groupby(scores_rerp, [:subject, :target, :time]),
    :r => mean => :r
)

fig_rerp = plot_roamm_decoding_individual(
    scores_rerp_subject,
    summary_rerp;
    title = "ROAMM — rERP decoding"
)
end
  ╠═╡ =#

# ╔═╡ 97b5db89-64e3-4d78-be29-89d57f4ebc28
# ╠═╡ disabled = true
#=╠═╡
begin
    scores_rerp = load_roamm_decoding_scores(joinpath(base_dir, "rerp"))
    summary_rerp = summarize_decoding(scores_rerp)
    
    fig_rerp = plot_roamm_decoding_individual(
        summary_rerp;
        title = "ROAMM — rERP decoding (group mean)"
    )

    
    #save("roamm_rerp_decoding.svg", fig_rerp)

    fig_rerp
end
  ╠═╡ =#

# ╔═╡ 8f43b5ca-8c27-4180-9fbd-975fd881c762
begin
    scores_rerp = load_roamm_decoding_scores(
        joinpath(base_dir, "rerp")
    )

    summary_rerp = summarize_decoding(scores_rerp)

    fig_rerp = plot_roamm_decoding_individual(
        scores_rerp,
        summary_rerp;
        title = "ROAMM — rERP decoding"
    )

    # save("roamm_rerp_decoding.svg", fig_rerp)

    fig_rerp
end

# ╔═╡ Cell order:
# ╠═e9ce96e4-ca36-11f1-84a8-df0069b540ab
# ╠═4001b213-854f-4305-8238-7a1742945c9a
# ╟─d63e9e0e-1a2a-4832-be56-1a5de5085976
# ╠═74612582-bead-443b-95b0-d123fbb2bb62
# ╠═f68010fb-0e8f-442f-96d4-e1b37903be28
# ╠═28c97b20-219b-4913-b147-e4e422a0bca0
# ╠═75e697c5-db1c-4d53-a95d-9f5615ffd99a
# ╠═01fcdcce-3ece-45b8-9e25-7101fa06c6bc
# ╠═57ec5561-d008-407f-99b7-846ef9b7c265
# ╠═1f536d81-ad42-4fa7-b5e3-9cc5d430cabd
# ╠═58460a3e-9c05-46b4-a402-f7f65ed79b15
# ╟─57b707fc-6089-485a-acfb-139ee57f52ce
# ╟─7dc15294-de07-48d0-89eb-6b618a7663fe
# ╠═e0b2923d-ee37-406e-92cd-5ce65c46b039
# ╠═2624f753-ee4e-4c04-867f-7e1028045d88
# ╟─3eacb1ef-27af-40bf-a738-b4aad27fa9b9
# ╟─9236fe57-780a-4028-bef6-a90399ef492b
# ╠═30d5aa98-a963-44ab-bded-24058687f5df
# ╠═8f43b5ca-8c27-4180-9fbd-975fd881c762
# ╠═4ce4326c-69cf-4b51-82a3-8b2e627c78a2
# ╠═97b5db89-64e3-4d78-be29-89d57f4ebc28
# ╟─130e150c-a908-4d70-9b1d-356a6009416b
# ╠═e8c2b476-b90c-4939-8360-9e871aac7118
# ╠═6ff590fa-5e6b-4cd3-bcb0-063c9ee4cc95
# ╠═44cdb525-2581-43f5-9137-127fb2d8d2ef
# ╠═7d93eff2-0bec-4423-8db1-0439527dd2fa
# ╠═0d646518-4b97-42c6-ac1e-b7c8f9a36df1
# ╠═6b99c440-4177-45ea-8c12-3c1faa0db6fd
# ╟─d05990b6-7035-478a-889a-73dd8be464d9
# ╠═65bd6d0c-7db3-4087-ae5c-33433b67d613
# ╠═7b1ddef0-eacf-4b04-9bbe-c08151b7c1c6
# ╠═c520c71d-9d96-4883-93ae-a691d964970b
# ╟─7770b6ef-7c84-4d81-a4e9-f35a1b6b0d11
# ╠═8cfdac79-78da-40c9-9e5d-4ae43bb0cb18
