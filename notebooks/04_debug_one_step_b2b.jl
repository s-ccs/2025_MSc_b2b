### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 8925770f-fa3b-440b-8ff1-9a29ee70d93e
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ bdf08ee2-4ec8-4fa2-87fb-9bbc7ce69029
using PlutoLinks: @revise, @ingredients

# ╔═╡ 1464865c-a796-49dc-a125-c333be11117d
@revise using MScB2B

# ╔═╡ 140dbed4-c56d-4dad-b162-1308d743f2b3
begin
    using Serialization
    using DataFrames
    using CairoMakie
end

# ╔═╡ 42b1c426-bd86-48ab-912e-17f038efbde1
begin
	project_root = normpath(joinpath(@__DIR__, ".."))
	
	results_dir = joinpath(project_root, "results")
end

# ╔═╡ aa265660-52c6-4509-ba9f-8acc0791b234
one_step_b2b = deserialize(
	    joinpath(results_dir, "one_step_b2b_hanning_1500trials_results.jls"))

# ╔═╡ fde9f5c6-d65d-49af-b1ce-f8fb3a2d1c5d
scores_one_step_b2b = one_step_b2b.score_tables

# ╔═╡ 0198155c-b34f-4d63-b6a6-41e5b4465717
fig_one_step_b2b_condition = MScB2B.plot_b2b_grid(
	scores_one_step_b2b,
	one_step_b2b.cfg;
	target = :condition,
	x_window= (-0.1, 0.6)
)

# ╔═╡ 5aa42d95-ca20-45ac-85b8-1a2c519e9b5b
fig_one_step_b2b_continuous = MScB2B.plot_b2b_grid(
		scores_one_step_b2b,
		one_step_b2b.cfg;
		target = :continuous,
		x_window = (-0.1, 0.6)
	)

# ╔═╡ Cell order:
# ╠═8925770f-fa3b-440b-8ff1-9a29ee70d93e
# ╠═bdf08ee2-4ec8-4fa2-87fb-9bbc7ce69029
# ╠═1464865c-a796-49dc-a125-c333be11117d
# ╠═140dbed4-c56d-4dad-b162-1308d743f2b3
# ╠═42b1c426-bd86-48ab-912e-17f038efbde1
# ╠═aa265660-52c6-4509-ba9f-8acc0791b234
# ╠═fde9f5c6-d65d-49af-b1ce-f8fb3a2d1c5d
# ╠═0198155c-b34f-4d63-b6a6-41e5b4465717
# ╠═5aa42d95-ca20-45ac-85b8-1a2c519e9b5b
