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

# ╔═╡ ebf4dc54-fcd1-4c8c-bd6d-d12a1dd8dfdc
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 930ebe90-9996-11f1-9773-990c6272f37b
using PlutoLinks: @revise, @ingredients

# ╔═╡ 5600884f-4488-4313-bf13-fd782d5ce705
@revise using MScB2B

# ╔═╡ 952709b5-78a6-4fb6-91b1-a352cbce67d0
begin
	using CairoMakie
	using DataFrames
	using Statistics
end

# ╔═╡ afb7b18a-ea7d-4ca7-ae76-5736fc8052fd
Controls = @ingredients(joinpath(@__DIR__, "simulation_controls.jl"))

# ╔═╡ e6a94e55-5591-44ce-be7d-9e26c1ceaa20
@bind sim Controls.simulation_controls()

# ╔═╡ 25844241-84ec-4997-a577-cfceade26837
cfg = MScB2B.SimulationConfig(
    n_trials = sim.n_trials,

    noiselevel = sim.noiselevel,
    channel_noise_sd = sim.channel_noise_sd,

    β0_n170 = sim.β0_n170,
    β0_p300 = sim.β0_p300,

    β_condition = sim.β_condition,
    β_continuous = sim.β_continuous,

    rho = sim.rho,

    overlap_interval_ms = sim.overlap_interval_ms,
    onset_condition_bias = sim.onset_condition_bias,
    shift_onset = sim.shift_onset,
)

# ╔═╡ d1c7c42c-7a67-4eb3-8588-7cc3397aa60a
cases = MScB2B.simulate_cases(cfg)

# ╔═╡ f1ec1326-55c0-463a-b396-01b9f375e9aa
begin
	clean = cases.clean
	
	(
	    data_size = size(clean.continuous),
	    n_events = nrow(clean.events_continuous),
	    first_latency = first(clean.events_continuous.latency),
	    last_latency = last(clean.events_continuous.latency),
	)
end

# ╔═╡ 1945d498-3c32-4306-ada6-0ac3d397144e
rerp_clean =
    MScB2B.fit_rerp_case(
        cfg,
        cases.clean,
        MScB2B.make_ridge_tuned_model();
        target = :continuous,
        nfolds = 3,
    )

# ╔═╡ c5b6523f-601b-4501-8e83-ae4e4dbb9eda
rerp_clean_scores =
    MScB2B.Unfold.coeftable(
        rerp_clean;
        measure = Statistics.cor,
        averaged = true,
    )

# ╔═╡ 31078d82-0855-4667-b802-ffc3922e871e
rerp_condition =
    MScB2B.run_rerp_decoding(
        cfg,
        cases;
        target = :condition_num,
        nfolds = 3,
    )

# ╔═╡ ebfc4ca9-88a5-4171-814f-f8fffbe29f2e
rerp_continuous =
    MScB2B.run_rerp_decoding(
        cfg,
        cases;
        target = :continuous,
        nfolds = 3,
    )

# ╔═╡ 602208c3-e1f9-4254-ac31-340476470056
MScB2B.plot_standard_decoding_grid(
    rerp_condition.scores,
    cfg;
    target = :condition,
)

# ╔═╡ 3484cd86-3245-4fa6-88fd-7774996ccd82
MScB2B.plot_standard_decoding_grid(
    rerp_continuous.scores,
    cfg;
    target = :continuous,
)

# ╔═╡ Cell order:
# ╠═ebf4dc54-fcd1-4c8c-bd6d-d12a1dd8dfdc
# ╠═930ebe90-9996-11f1-9773-990c6272f37b
# ╠═5600884f-4488-4313-bf13-fd782d5ce705
# ╠═952709b5-78a6-4fb6-91b1-a352cbce67d0
# ╠═afb7b18a-ea7d-4ca7-ae76-5736fc8052fd
# ╠═e6a94e55-5591-44ce-be7d-9e26c1ceaa20
# ╠═25844241-84ec-4997-a577-cfceade26837
# ╠═d1c7c42c-7a67-4eb3-8588-7cc3397aa60a
# ╠═f1ec1326-55c0-463a-b396-01b9f375e9aa
# ╠═1945d498-3c32-4306-ada6-0ac3d397144e
# ╠═c5b6523f-601b-4501-8e83-ae4e4dbb9eda
# ╠═31078d82-0855-4667-b802-ffc3922e871e
# ╠═ebfc4ca9-88a5-4171-814f-f8fffbe29f2e
# ╠═602208c3-e1f9-4254-ac31-340476470056
# ╠═3484cd86-3245-4fa6-88fd-7774996ccd82
