### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ e9ce96e4-ca36-11f1-84a8-df0069b540ab
begin
    import Pkg
	Pkg.activate(joinpath(@__DIR__, "..", ".."))
end

# ╔═╡ 4001b213-854f-4305-8238-7a1742945c9a
using Serialization, DataFrames

# ╔═╡ 8bdced5b-26e5-435c-b3dc-6710fb27c5be
using CairoMakie, Statistics

# ╔═╡ f77b3a99-50f6-40bd-ba12-d56dcd107242
begin
	standard_decoding = deserialize("/store/users/xu/2025_MSc_b2b/results/models/roamm/standard/sub-10014_standard_ridge.jls")
	keys(standard_decoding)
end

# ╔═╡ f4b1719a-d6ba-4172-8d77-45b3ad560f71
first(standard_decoding.scores, 5)

# ╔═╡ c8f859bc-ab5b-4a03-b4b1-60206a4b22e1
 names(standard_decoding.scores)

# ╔═╡ dcf2f4e4-04d1-4b09-8582-1859843c072d
unique(standard_decoding.scores.target)

# ╔═╡ d63e9e0e-1a2a-4832-be56-1a5de5085976
md"""
# Step1: load all subject results"""

# ╔═╡ 486482b5-2120-41bf-b749-21f314180fc4
begin 
	function load_roamm_standard_scores(results_dir)
		files = filter(f -> endswith(f, ".jls"), readdir(results_dir; join = true))
		isempty(files) && error("No .jls files found in $results_dir")

		tables = DataFrame[]
		for file in files
			result = deserialize(file)
			scores = DataFrame(result.scores)
			scores[!, :subject] = fill(String(result.subject), nrow(scores))
			push!(tables, scores)
		end

		return vcat(tables...)
	end
end

# ╔═╡ e1890edc-0177-4a96-ab87-b0e392a07b22
scores_all = load_roamm_standard_scores("/store/users/xu/2025_MSc_b2b/results/models/roamm/standard")

# ╔═╡ 0fcebb9e-ede3-4814-963e-103e32a2fb5c
md"""
# Step2: average across subjects
"""

# ╔═╡ 8d5901a1-72b7-4fba-83c1-9b938426ca54
function summarize_roamm_scores(scores_all)
	g = groupby(scores_all, [:target, :time])

	combine(
		g,
		:r => mean => :r_mean,
		:r => std => :r_std,
		nrow => :n
	)
end

# ╔═╡ 725d8080-5f43-412a-8384-eee1d0dfb2ce
scores_mean = summarize_roamm_scores(scores_all)

# ╔═╡ 57b707fc-6089-485a-acfb-139ee57f52ce
md"""
# Plotting"""

# ╔═╡ 944b1cad-843d-46a7-8d9a-956eaa557a65
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

# ╔═╡ c507a9a1-e725-4085-8a30-9c58c422d6da
begin
	
	fig = plot_roamm_standard_decoding(scores_mean)
	#save("results/plots/roamm_standard_decoding_mean.svg", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═e9ce96e4-ca36-11f1-84a8-df0069b540ab
# ╠═4001b213-854f-4305-8238-7a1742945c9a
# ╠═f77b3a99-50f6-40bd-ba12-d56dcd107242
# ╠═f4b1719a-d6ba-4172-8d77-45b3ad560f71
# ╠═c8f859bc-ab5b-4a03-b4b1-60206a4b22e1
# ╠═dcf2f4e4-04d1-4b09-8582-1859843c072d
# ╠═d63e9e0e-1a2a-4832-be56-1a5de5085976
# ╠═486482b5-2120-41bf-b749-21f314180fc4
# ╠═e1890edc-0177-4a96-ab87-b0e392a07b22
# ╠═0fcebb9e-ede3-4814-963e-103e32a2fb5c
# ╠═8bdced5b-26e5-435c-b3dc-6710fb27c5be
# ╠═8d5901a1-72b7-4fba-83c1-9b938426ca54
# ╠═725d8080-5f43-412a-8384-eee1d0dfb2ce
# ╠═57b707fc-6089-485a-acfb-139ee57f52ce
# ╠═944b1cad-843d-46a7-8d9a-956eaa557a65
# ╠═c507a9a1-e725-4085-8a30-9c58c422d6da
