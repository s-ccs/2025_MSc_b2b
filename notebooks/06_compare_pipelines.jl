### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ f390cf4c-9e73-11f1-b1f4-1dfed743fa21
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 607a836b-a82f-4e35-9980-51fab46d8c03
using PlutoLinks: @revise, @ingredients

# ╔═╡ f529b8c2-b7d8-450a-b29b-b53310702a3b
@revise using MScB2B

# ╔═╡ 48f54015-00dd-4a42-aef4-ea36cf3baa46
begin
    using Serialization
    using DataFrames
    using CairoMakie
end

# ╔═╡ a8b46923-0f34-45ed-a965-6a1ca27c498b
begin
	project_root = normpath(joinpath(@__DIR__, ".."))
	
	results_dir = joinpath(project_root, "results")
end

# ╔═╡ 0d11f181-fe2f-4796-8981-f9fa9ef1d761
begin
	standard = deserialize(
	    joinpath(results_dir, "01_standard_decoding_results.jls")
	)
	
	rerp = deserialize(
	    joinpath(results_dir, "02_rerp_decoding_results.jls")
	)
	
	plain_b2b = deserialize(
	    joinpath(results_dir, "03_plain_b2b_results.jls")
	)
	
	one_step_b2b = deserialize(
	    joinpath(results_dir, "04_one_step_b2b_results.jls")
	)
	
	two_step_b2b = deserialize(
	    joinpath(results_dir, "05_two_step_b2b_results.jls")
	)
end

# ╔═╡ b3f15889-7799-4ab1-b4f5-f9f337a6e876
# ==========================================
# 01 - Standard decoding
# ==========================================
begin
	scores_standard_condition = standard.condition_results.scores
	scores_standard_continuous = standard.continuous_results.scores
end

# ╔═╡ 1a83acf7-fa70-4124-b54f-3ab55e9ea419
fig_standard_condition = MScB2B.plot_standard_decoding_grid(
	scores_standard_condition,
	standard.cfg;
	target = :condition
)

# ╔═╡ 4eafcc8f-5ee6-44ee-b3c8-f9a3301b899c
fig_standard_continuous = MScB2B.plot_standard_decoding_grid(
	scores_standard_continuous,
	standard.cfg;
	target = :continuous
)
#save("standard_decoding_continuous_tuned_lambda.svg", fig_standard_continuous)

# ╔═╡ fa98fc46-6ec8-44f5-a9d5-a1b3a8bdfd41
# =========================================
# 02 - rerp decoding
# ========================================
begin
	scores_rerp_condition = rerp.condition_results.scores
	scores_rerp_continuous = rerp.continuous_results.scores
end

# ╔═╡ a09e150c-63d4-41e7-aae7-42d6b95ecad9
fig_rerp_condition = MScB2B.plot_standard_decoding_grid(
		scores_rerp_condition,
		standard.cfg;
		target = :condition
	)

# ╔═╡ 7ecf5bc1-6198-4e2e-aab6-c457cd485101
fig_rerp_continuous = MScB2B.plot_standard_decoding_grid(
	scores_rerp_continuous,
	standard.cfg;
	target = :continuous
)

# ╔═╡ a1b8bacc-67d1-400f-a796-4f90b3a3db07
# =======================================
# 03 - plain B2B
# ========================================
scores_plain_b2b = plain_b2b.score_tables

# ╔═╡ 4754f6e5-9092-48f3-a455-15f76b5f2b9d
fig_plain_b2b_condition = MScB2B.plot_b2b_grid(
	scores_plain_b2b,
	plain_b2b.cfg;
	target = :condition
)

# ╔═╡ 9e5bed10-48c6-472c-96a5-ef2be8747cd4
fig_plain_b2b_continuous = MScB2B.plot_b2b_grid(
		scores_plain_b2b,
		plain_b2b.cfg;
		target = :continuous
	)

# ╔═╡ 494ad80a-7133-42f6-be53-a7a43243740a
# ===================================
# 04 - one step B2B
# =====================================
scores_one_step_b2b = one_step_b2b.score_tables

# ╔═╡ 123ed81f-a95f-4861-ab48-e3a6fb736b06
fig_one_step_b2b_condition = MScB2B.plot_b2b_grid(
	scores_one_step_b2b,
	one_step_b2b.cfg;
	target = :condition
)

# ╔═╡ 95c902cf-19c7-4748-bc9f-9e549c11c98e
fig_one_step_b2b_continuous = MScB2B.plot_b2b_grid(
		scores_one_step_b2b,
		one_step_b2b.cfg;
		target = :continuous
	)

# ╔═╡ 03a985fd-74ee-4d3c-ba7f-38252e50c0f7
# ======================================
# 05 -  two step B2B
# ======================================
scores_two_step_b2b = two_step_b2b.score_tables

# ╔═╡ 02cb3245-4a38-47f0-847f-e91dec3e61c0
fig_two_step_b2b_condition = MScB2B.plot_b2b_grid(
	scores_two_step_b2b,
	two_step_b2b.cfg;
	target = :condition
)

# ╔═╡ e0b923ca-5989-44a4-a2e7-7d8fbe246d58
fig_two_step_b2b_continuous = MScB2B.plot_b2b_grid(
		scores_two_step_b2b,
		two_step_b2b.cfg;
		target = :continuous
	)

# ╔═╡ 25ab2016-3821-4783-bbed-5d1deb1af7d5
begin
    plot_dir = joinpath(@__DIR__, "..", "plots")

    # 01 - Standard decoding
    save(
        joinpath(plot_dir, "01_standard_condition.svg"),
        fig_standard_condition
    )

    save(
        joinpath(plot_dir, "01_standard_continuous.svg"),
        fig_standard_continuous
    )

    # 02 - rERP decoding
    save(
        joinpath(plot_dir, "02_rerp_condition.svg"),
        fig_rerp_condition
    )

    save(
        joinpath(plot_dir, "02_rerp_continuous.svg"),
        fig_rerp_continuous
    )

    # 03 - Plain B2B
    save(
        joinpath(plot_dir, "03_plain_b2b_condition.svg"),
        fig_plain_b2b_condition
    )

    save(
        joinpath(plot_dir, "03_plain_b2b_continuous.svg"),
        fig_plain_b2b_continuous
    )

    # 04 - One-step B2B
    save(
        joinpath(plot_dir, "04_one_step_b2b_condition.svg"),
        fig_one_step_b2b_condition
    )

    save(
        joinpath(plot_dir, "04_one_step_b2b_continuous.svg"),
        fig_one_step_b2b_continuous
    )

    # 05 - Two-step B2B
    save(
        joinpath(plot_dir, "05_two_step_b2b_condition.svg"),
        fig_two_step_b2b_condition
    )

    save(
        joinpath(plot_dir, "05_two_step_b2b_continuous.svg"),
        fig_two_step_b2b_continuous
    )

    "Plots saved to: $plot_dir"
end

# ╔═╡ Cell order:
# ╠═f390cf4c-9e73-11f1-b1f4-1dfed743fa21
# ╠═607a836b-a82f-4e35-9980-51fab46d8c03
# ╠═f529b8c2-b7d8-450a-b29b-b53310702a3b
# ╠═48f54015-00dd-4a42-aef4-ea36cf3baa46
# ╠═a8b46923-0f34-45ed-a965-6a1ca27c498b
# ╠═0d11f181-fe2f-4796-8981-f9fa9ef1d761
# ╠═b3f15889-7799-4ab1-b4f5-f9f337a6e876
# ╠═1a83acf7-fa70-4124-b54f-3ab55e9ea419
# ╠═4eafcc8f-5ee6-44ee-b3c8-f9a3301b899c
# ╠═fa98fc46-6ec8-44f5-a9d5-a1b3a8bdfd41
# ╠═a09e150c-63d4-41e7-aae7-42d6b95ecad9
# ╟─7ecf5bc1-6198-4e2e-aab6-c457cd485101
# ╠═a1b8bacc-67d1-400f-a796-4f90b3a3db07
# ╠═4754f6e5-9092-48f3-a455-15f76b5f2b9d
# ╠═9e5bed10-48c6-472c-96a5-ef2be8747cd4
# ╠═494ad80a-7133-42f6-be53-a7a43243740a
# ╠═123ed81f-a95f-4861-ab48-e3a6fb736b06
# ╠═95c902cf-19c7-4748-bc9f-9e549c11c98e
# ╠═03a985fd-74ee-4d3c-ba7f-38252e50c0f7
# ╠═02cb3245-4a38-47f0-847f-e91dec3e61c0
# ╠═e0b923ca-5989-44a4-a2e7-7d8fbe246d58
# ╠═25ab2016-3821-4783-bbed-5d1deb1af7d5
