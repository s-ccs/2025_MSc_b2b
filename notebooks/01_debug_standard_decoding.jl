### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 68307fde-996d-11f1-9b14-6d01c4c95307
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ cdffc3fd-e905-43f0-862e-fbc3e00b879b
using PlutoLinks: @revise, @ingredients

# ╔═╡ 75b38e4e-1946-4868-90da-acea59c546f5
@revise using MScB2B

# ╔═╡ 75e898d1-34ef-4cd0-8ce1-54460ace049a
begin
	using CairoMakie
	using DataFrames
	using Statistics
end

# ╔═╡ c9a67dce-c511-440a-8e55-94582196c4ff
Controls = @ingredients(joinpath(@__DIR__, "simulation_controls.jl"))

# ╔═╡ b6475db8-fcdd-4993-8b73-6de01629a7c9
@bind sim Controls.simulation_controls()

# ╔═╡ fab7bf06-f3ea-4fc9-bb14-605ffe3e34a7
common_cfg = (
    n_trials = sim.n_trials,
    sfreq = 100.0,
    n_channels = 20,
    noiselevel = sim.noiselevel,
    channel_noise_sd = sim.channel_noise_sd,
    rho = sim.rho,
    overlap_interval_ms = sim.overlap_interval_ms,
    shift_onset = sim.shift_onset,
)

# ╔═╡ cdf46f17-1510-4218-b180-7d88e10678e6
ridge_model = MScB2B.make_ridge_tuned_model(
	inner_nfolds = 2,
	resolution = 4
)

# ╔═╡ c1dbc92d-9f76-46b9-bc66-3b20d9f4a06c
md"""
# Condition + Continuous Simulation
"""

# ╔═╡ 42c0f09f-66e8-46e3-9c6d-4fc539044c30
cfg_cond = MScB2B.ConditionContinuousConfig(
    ;
    common_cfg...,

    n_trials = 1500,
    sfreq = 100.0,

    β0_n170 = sim.β0_n170,
    β_condition = sim.β_condition,

    β0_p300 = sim.β0_p300,
    β_continuous = sim.β_continuous,

    onset_condition_bias = sim.onset_condition_bias,
)

# ╔═╡ 4e43ae69-2fd9-4e83-9e32-3889d8daad6a
sim_cond = MScB2B.simulate_cond_cont_cases(cfg_cond; seed = 12)

# ╔═╡ e89f6ee9-f7cd-4922-9b1d-c9fafaaad1dc
standard_condition =
    MScB2B.run_standard_decoding(
		cfg_cond,
		sim_cond;
		model = ridge_model,
        target = :condition_num,
        nfolds = 2,
        seed = 12,
    )

# ╔═╡ e2ce9fc1-bfca-4faa-8e8f-91de1c58163d
standard_continuous =
    MScB2B.run_standard_decoding(
		cfg_cond,
		sim_cond;
		model = ridge_model,
        target = :continuous,
        nfolds = 2,
        seed = 12,
    )

# ╔═╡ 98e0a496-7f9f-4a23-8944-f56782378f8a
standard_condition.scores

# ╔═╡ 7e8d007b-4e32-475b-82c3-b665988254c6
keys(standard_condition.yhats)

# ╔═╡ 48a50dc7-3aa0-4bf9-a5f1-10622688016c
MScB2B.plot_standard_decoding_grid(
    standard_condition.scores,
    cfg_cond;
	target = :condition
)

# ╔═╡ 1e576cb3-f9eb-4311-b3c2-37e4a126897b
MScB2B.plot_standard_decoding_grid(
    standard_continuous.scores,
    cfg_cond;
	target = :continuous
)

# ╔═╡ a1171a55-019d-45c1-9ae8-b70b45faaf27
md"""
# Correlated Continuous Simulation
"""

# ╔═╡ 3b577bdc-068e-4ff9-bca6-7b9904df5b38
cfg_corr = MScB2B.CorrelatedContinuousConfig(
	;
    common_cfg...,

    β0_p100 = 5.0,
    β_p100 = 2.0,

    β0_n170 = 5.0,
    β_n170 = 3.0,

    β0_p300 = 5.0,
    β_p300 = 1.5,

    onset_predictor_bias = -0.6,
)

# ╔═╡ 9a49f770-4e1d-431c-91bc-63f0c017672f
sim_corr = MScB2B.simulate_corr_cont_cases(cfg_corr; seed = 12)

# ╔═╡ d87c3fe8-cb0d-4dc6-96f2-8a860391d77b
standard_corr =
    MScB2B.run_standard_decoding(
		cfg_corr,
		sim_corr;
		model = ridge_model,
        nfolds = 2,
        seed = 12,
    )

# ╔═╡ a3f109c1-7b09-477b-a310-a576b17ff360
MScB2B.plot_standard_decoding_grid(
    standard_corr.scores,
    cfg_corr;
	target = :continuous1)

# ╔═╡ Cell order:
# ╠═68307fde-996d-11f1-9b14-6d01c4c95307
# ╠═cdffc3fd-e905-43f0-862e-fbc3e00b879b
# ╠═75b38e4e-1946-4868-90da-acea59c546f5
# ╠═75e898d1-34ef-4cd0-8ce1-54460ace049a
# ╠═c9a67dce-c511-440a-8e55-94582196c4ff
# ╠═fab7bf06-f3ea-4fc9-bb14-605ffe3e34a7
# ╠═b6475db8-fcdd-4993-8b73-6de01629a7c9
# ╠═cdf46f17-1510-4218-b180-7d88e10678e6
# ╟─c1dbc92d-9f76-46b9-bc66-3b20d9f4a06c
# ╠═42c0f09f-66e8-46e3-9c6d-4fc539044c30
# ╠═4e43ae69-2fd9-4e83-9e32-3889d8daad6a
# ╠═e89f6ee9-f7cd-4922-9b1d-c9fafaaad1dc
# ╠═e2ce9fc1-bfca-4faa-8e8f-91de1c58163d
# ╠═98e0a496-7f9f-4a23-8944-f56782378f8a
# ╠═7e8d007b-4e32-475b-82c3-b665988254c6
# ╠═48a50dc7-3aa0-4bf9-a5f1-10622688016c
# ╠═1e576cb3-f9eb-4311-b3c2-37e4a126897b
# ╟─a1171a55-019d-45c1-9ae8-b70b45faaf27
# ╠═3b577bdc-068e-4ff9-bca6-7b9904df5b38
# ╠═9a49f770-4e1d-431c-91bc-63f0c017672f
# ╠═d87c3fe8-cb0d-4dc6-96f2-8a860391d77b
# ╠═a3f109c1-7b09-477b-a310-a576b17ff360
