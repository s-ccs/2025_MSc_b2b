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

# ╔═╡ 8b8cc14c-9d99-11f1-93c4-09291a590ac4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 49496a9e-03ac-42fd-a7c0-76ee3aac3c0c
using PlutoLinks: @revise, @ingredients

# ╔═╡ 58ed75f6-d941-4083-8e81-f1c4de3f1c03
@revise using MScB2B

# ╔═╡ ffa9faea-e1e1-4ad0-956c-d16b77dd1bf2
begin
	using CairoMakie
	using DataFrames
	using Statistics
    using Unfold
    using UnfoldDecode
    using StatsModels: @formula
end

# ╔═╡ 577347be-c791-4972-8e2c-099137d80077
Controls = @ingredients(joinpath(@__DIR__, "simulation_controls.jl"))

# ╔═╡ af6d5abb-efbf-475e-9efd-765f3c234ddc
@bind sim Controls.simulation_controls()

# ╔═╡ 266b38ea-3fee-45dd-bf3a-5504e4177227
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

# ╔═╡ 4b92f32a-8b10-4a3e-9577-5ad859ee251b
cases = MScB2B.simulate_cases(cfg)

# ╔═╡ 10ce8434-6e0a-4459-a25b-189abd874c08
begin
	case_data = cases.clean
	
	dat_cont = case_data.continuous
	evts = case_data.events_continuous
end

# ╔═╡ 6e663d4d-a433-4b14-ae94-3dafbdeef40e
(
	size(dat_cont),
	nrow(evts),
	names(evts),
	unique(evts.event)
)

# ╔═╡ f954ca62-610b-4c4b-8966-35dcdfe37471
begin
	design_rerp = [
	    "stimulus" => (
	        @formula(0 ~ 1 + condition + continuous),
	        Unfold.firbasis(
	            τ = [-0.1, 1.0],
	            sfreq = cfg.sfreq,
	            name = "stimulus",
	        ),
	    ),
	]
	
	uf_rerp = Unfold.fit(
	    UnfoldLinearModelContinuousTime,
	    design_rerp,
	    evts,
	    dat_cont;
	    eventcolumn = :event,
	)
end

# ╔═╡ a03528f8-af5b-4f1e-88d0-8d30cadc76f9
(
	typeof(uf_rerp),
	Unfold.times(uf_rerp)[1],
	Unfold.basisname(uf_rerp)
)

# ╔═╡ 134faf2e-4c93-4531-9edc-19ee0bd19352
times = Unfold.times(uf_rerp)[1]

# ╔═╡ b13d9cfb-e897-4760-af30-5577f485a6a3
X_corrected = UnfoldDecode.singletrials(
    dat_cont,
    uf_rerp,
    evts,
    "stimulus",
    :event,
)

# ╔═╡ 9df66431-c596-4c79-9aa6-dc52f4ba600a
size(X_corrected)

# ╔═╡ b6a65070-3cd6-4aca-bafc-516c84e28342
@assert size(X_corrected, 1) == cfg.n_channels

# ╔═╡ a8cd339e-0329-4baf-ab39-170020186472
@assert size(X_corrected, 2) == length(times)

# ╔═╡ 3529f7cf-df20-4547-b103-4dfc3fff01b1
@assert size(X_corrected, 3) == nrow(evts)

# ╔═╡ 3054b4b7-52e0-4577-8dc5-a9ca949965ee
eltype(X_corrected)

# ╔═╡ 07f071f5-8d0d-4a11-a2b2-6bd27a85720e
extrema(skipmissing(vec(X_corrected)))

# ╔═╡ 7d18bdff-f362-4f77-802a-c6c6cc731b05
design_b2b = [
    "stimulus" => (
        @formula(0 ~ 1 + condition + continuous),
        times,
    ),
]

# ╔═╡ 4c0e2532-9b80-4ae6-a74e-a890e621d1e6
b2b_solver = (X, y) ->
    UnfoldDecode.solver_b2b(
        X,
        y;
        cross_val_reps = 3,
    )

# ╔═╡ 8be4de65-74ff-43b9-a017-21e2e6aea411
X_test = Unfold.modelmatrix(
    Unfold.fit(
        UnfoldModel,
        design_b2b,
        evts,
        X_corrected;
        solver = (X, y) -> begin
            @show size(X)
            @show size(y)
            @show eltype(y)

            X_clean, y_clean =
                Unfold.drop_missing_epochs(X, y)

            @show size(X_clean)
            @show size(y_clean)
            @show eltype(y_clean)

            UnfoldDecode.solver_b2b(
                X_clean,
                y_clean;
                cross_val_reps = 1,
            )
        end,
        eventcolumn = :event,
    )
)

# ╔═╡ 8b01d59b-8db8-4e46-bfe3-30e76dfa2821
uf_b2b = Unfold.fit(
    UnfoldModel,
    design_b2b,
    evts,
    X_corrected;
    solver = b2b_solver,
    eventcolumn = :event,
)

# ╔═╡ Cell order:
# ╠═8b8cc14c-9d99-11f1-93c4-09291a590ac4
# ╠═49496a9e-03ac-42fd-a7c0-76ee3aac3c0c
# ╠═58ed75f6-d941-4083-8e81-f1c4de3f1c03
# ╠═ffa9faea-e1e1-4ad0-956c-d16b77dd1bf2
# ╠═577347be-c791-4972-8e2c-099137d80077
# ╠═af6d5abb-efbf-475e-9efd-765f3c234ddc
# ╠═266b38ea-3fee-45dd-bf3a-5504e4177227
# ╠═4b92f32a-8b10-4a3e-9577-5ad859ee251b
# ╠═10ce8434-6e0a-4459-a25b-189abd874c08
# ╠═6e663d4d-a433-4b14-ae94-3dafbdeef40e
# ╠═a03528f8-af5b-4f1e-88d0-8d30cadc76f9
# ╠═f954ca62-610b-4c4b-8966-35dcdfe37471
# ╠═134faf2e-4c93-4531-9edc-19ee0bd19352
# ╠═b13d9cfb-e897-4760-af30-5577f485a6a3
# ╠═9df66431-c596-4c79-9aa6-dc52f4ba600a
# ╠═b6a65070-3cd6-4aca-bafc-516c84e28342
# ╠═a8cd339e-0329-4baf-ab39-170020186472
# ╠═3529f7cf-df20-4547-b103-4dfc3fff01b1
# ╠═3054b4b7-52e0-4577-8dc5-a9ca949965ee
# ╠═07f071f5-8d0d-4a11-a2b2-6bd27a85720e
# ╠═7d18bdff-f362-4f77-802a-c6c6cc731b05
# ╠═4c0e2532-9b80-4ae6-a74e-a890e621d1e6
# ╠═8be4de65-74ff-43b9-a017-21e2e6aea411
# ╠═8b01d59b-8db8-4e46-bfe3-30e76dfa2821
