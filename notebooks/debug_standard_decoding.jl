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

# ╔═╡ bfb4a3f6-4983-43cd-bc0b-7c4bb8f22669
cfg = MScB2B.SimulationConfig(
    n_trials = sim.n_trials,

    noiselevel = sim.noiselevel,
    channel_noise_sd = sim.channel_noise_sd,

    β0_n170 = sim.β0,
    β0_p300 = sim.β0,

    β_condition = sim.β_condition,
    β_continuous = sim.β_continuous,

    rho = sim.rho,

    overlap_interval_ms = sim.overlap_interval_ms,
    onset_condition_bias = sim.onset_condition_bias,
    shift_onset = sim.shift_onset,
)

# ╔═╡ c166cfe9-28b0-47b1-9e1c-4f8fefcfb08e
cases = MScB2B.simulate_cases(cfg)

# ╔═╡ e89f6ee9-f7cd-4922-9b1d-c9fafaaad1dc
standard_condition =
    MScB2B.run_standard_decoding(
        cfg,
        cases;
        target = :condition_num,
        nfolds = 2,
        seed = 12,
    )

# ╔═╡ e2ce9fc1-bfca-4faa-8e8f-91de1c58163d
standard_continuous =
    MScB2B.run_standard_decoding(
        cfg,
        cases;
        target = :continuous,
        nfolds = 2,
        seed = 12,
    )

# ╔═╡ 48a50dc7-3aa0-4bf9-a5f1-10622688016c
MScB2B.plot_standard_decoding_grid(
    standard_condition.scores,
    cfg;
    target = :condition,
)

# ╔═╡ 1e576cb3-f9eb-4311-b3c2-37e4a126897b
MScB2B.plot_standard_decoding_grid(
    standard_continuous.scores,
    cfg;
    target = :continuous,
)

# ╔═╡ Cell order:
# ╠═68307fde-996d-11f1-9b14-6d01c4c95307
# ╠═cdffc3fd-e905-43f0-862e-fbc3e00b879b
# ╠═75b38e4e-1946-4868-90da-acea59c546f5
# ╠═75e898d1-34ef-4cd0-8ce1-54460ace049a
# ╠═c9a67dce-c511-440a-8e55-94582196c4ff
# ╠═bfb4a3f6-4983-43cd-bc0b-7c4bb8f22669
# ╠═c166cfe9-28b0-47b1-9e1c-4f8fefcfb08e
# ╠═e89f6ee9-f7cd-4922-9b1d-c9fafaaad1dc
# ╠═e2ce9fc1-bfca-4faa-8e8f-91de1c58163d
# ╠═48a50dc7-3aa0-4bf9-a5f1-10622688016c
# ╠═1e576cb3-f9eb-4311-b3c2-37e4a126897b
# ╠═b6475db8-fcdd-4993-8b73-6de01629a7c9
