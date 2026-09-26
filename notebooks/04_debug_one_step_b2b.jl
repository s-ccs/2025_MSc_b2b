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
	    joinpath(results_dir, "0_rho_one_step_b2b_hanning_1500trials_results.jls"))

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

# ╔═╡ 80a7fc5b-5e66-44eb-a5f7-d79c63bf0bd7
begin
	save("0_rho_one_step_b2b_condition_hanning_1500_trials.svg",
		 fig_one_step_b2b_condition)

	save("0_rho_one_step_b2b_continuous_hanning_1500_trials.svg",
		fig_one_step_b2b_continuous)
end

# ╔═╡ d359f16d-b719-4ea9-863e-7b511512b373
scores_one_step_b2b

# ╔═╡ 04490951-027d-4e66-b585-f908a17b4d88
AoG = MScB2B.UnfoldMakie.AlgebraOfGraphics

# ╔═╡ f2927a02-ebe5-490d-ad7c-e0a7200c89e5
begin
	h = AoG.data(subset(scores_one_step_b2b,:coefname=> x-> x.=="continuous"))*
	AoG.mapping(:time,:estimate,color="case",layout="case"=>AoG.sorter(["clean","overlap","confound","both"]))*AoG.visual(AoG.Lines)|>AoG.draw
end

# ╔═╡ e6070312-cbb4-4689-9e30-7c01c6df6299
h.figure

# ╔═╡ 2c21af61-407d-4fdb-b328-4733fc4f51ab
MScB2B.plot_erp(subset(scores_one_step_b2b,:coefname=> x-> x.=="continuous"),mapping=(;color=:case,layout=:case=>AoG.sorter(["clean","overlap","confound","both"])))

# ╔═╡ 8b8c3271-fe3a-40cb-815a-36f40356b518
scores_one_step_b2b

# ╔═╡ Cell order:
# ╠═8925770f-fa3b-440b-8ff1-9a29ee70d93e
# ╠═bdf08ee2-4ec8-4fa2-87fb-9bbc7ce69029
# ╠═1464865c-a796-49dc-a125-c333be11117d
# ╠═140dbed4-c56d-4dad-b162-1308d743f2b3
# ╠═42b1c426-bd86-48ab-912e-17f038efbde1
# ╠═aa265660-52c6-4509-ba9f-8acc0791b234
# ╠═fde9f5c6-d65d-49af-b1ce-f8fb3a2d1c5d
# ╠═0198155c-b34f-4d63-b6a6-41e5b4465717
# ╠═80a7fc5b-5e66-44eb-a5f7-d79c63bf0bd7
# ╠═5aa42d95-ca20-45ac-85b8-1a2c519e9b5b
# ╠═d359f16d-b719-4ea9-863e-7b511512b373
# ╠═04490951-027d-4e66-b585-f908a17b4d88
# ╠═f2927a02-ebe5-490d-ad7c-e0a7200c89e5
# ╠═e6070312-cbb4-4689-9e30-7c01c6df6299
# ╠═2c21af61-407d-4fdb-b328-4733fc4f51ab
# ╠═8b8c3271-fe3a-40cb-815a-36f40356b518
