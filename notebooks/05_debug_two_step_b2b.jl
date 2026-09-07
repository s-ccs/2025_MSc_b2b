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
	using Serialization
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

# ╔═╡ 4d25fb12-c5b3-48fa-a637-5ac4f7a9fe7e
two_step = MScB2B.run_two_step_b2b(
    cfg,
    cases;
    cross_val_reps = 3,
)

# ╔═╡ 7c2e277a-36a1-4606-acf5-71cb32c489e4
fig_two_step_b2b_condition_hanning = MScB2B.plot_b2b_grid(
    two_step.score_tables,
    cfg;
    target = :condition,
	x_window = (-0.1, 0.6),
)

# ╔═╡ 3d5e3d62-c554-4256-9672-0042c4204a68
fig_two_step_b2b_continuous_hanning = MScB2B.plot_b2b_grid(
    two_step.score_tables,
    cfg;
    target = :continuous,
	x_window = (-0.1, 0.6),
)

# ╔═╡ 2472b378-d948-40b0-98e1-08f54dd2ee40
# ╠═╡ disabled = true
#=╠═╡
begin
save( "05_two_step_b2b_condition_hanning_1500trials_long_windows.svg",
        fig_two_step_b2b_condition_hanning
    )
save( "05_two_step_b2b_continuous_hanning_1500trials_long_windows.svg",
        fig_two_step_b2b_continuous_hanning
    )
end
  ╠═╡ =#

# ╔═╡ 75df5261-b6b4-482a-aeb5-a9c4c5cb99f4


# ╔═╡ ddbfb942-9591-40e0-bd98-75d995fcbd02


# ╔═╡ 5de2d518-5bd3-404d-81bd-9a8f0b852a01


# ╔═╡ f7690397-cecd-4184-8da4-962d50be7d16


# ╔═╡ a0cd698c-0a73-4918-91d5-920a3fb45177


# ╔═╡ 8a6416e9-2c7d-49fc-9d31-fbc09003adcd


# ╔═╡ 9fe10962-9032-4a3a-a91c-3dc63491d7ba


# ╔═╡ 6f7dc6f2-f4ab-4a15-9da2-94e189b2aebb


# ╔═╡ 43fe4334-da86-4f11-85c0-d50e5998fbac


# ╔═╡ 98f309fc-1d8c-4b18-a65c-1412e541b81a


# ╔═╡ 5864388e-3fb0-4723-86d2-df6662efb0d3


# ╔═╡ f27fc2cf-411a-4544-ba94-15153f885da2


# ╔═╡ f9eb4eec-a141-4574-bdc4-a39e9d9f36e6


# ╔═╡ 7edc508a-5147-4760-855b-5ceaf5228212


# ╔═╡ b76b52da-3df5-43b0-89ee-51836a021c47


# ╔═╡ 434a0d07-58aa-44ea-a030-09d947b08bc8


# ╔═╡ 4bb1837a-179b-4cd1-8aa4-0bdb98d47704


# ╔═╡ d020d01d-c910-4833-a633-320279a01253


# ╔═╡ b8757652-0ce5-425c-a47e-40e36f37cc4f


# ╔═╡ b0a9a532-8856-4ca4-bd7f-18314ff2f87d


# ╔═╡ 3b55d105-a21a-4604-9224-128bf3867ef8


# ╔═╡ b86bb227-b008-41ad-ba74-faf8c18c9104


# ╔═╡ 10ce8434-6e0a-4459-a25b-189abd874c08
# ╠═╡ disabled = true
#=╠═╡
begin
	case_data = cases.clean
	
	dat_cont = case_data.continuous
	evts = case_data.events_continuous
end
  ╠═╡ =#

# ╔═╡ 6e663d4d-a433-4b14-ae94-3dafbdeef40e
# ╠═╡ disabled = true
#=╠═╡
(
	size(dat_cont),
	nrow(evts),
	names(evts),
	unique(evts.event)
)
  ╠═╡ =#

# ╔═╡ a03528f8-af5b-4f1e-88d0-8d30cadc76f9
# ╠═╡ disabled = true
#=╠═╡
(
	typeof(uf_rerp),
	Unfold.times(uf_rerp)[1],
	Unfold.basisname(uf_rerp)
)
  ╠═╡ =#

# ╔═╡ f954ca62-610b-4c4b-8966-35dcdfe37471
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 134faf2e-4c93-4531-9edc-19ee0bd19352
# ╠═╡ disabled = true
#=╠═╡
times = Unfold.times(uf_rerp)[1]
  ╠═╡ =#

# ╔═╡ b13d9cfb-e897-4760-af30-5577f485a6a3
# ╠═╡ disabled = true
#=╠═╡
X_corrected = UnfoldDecode.singletrials(
    dat_cont,
    uf_rerp,
    evts,
    "stimulus",
    :event,
)
  ╠═╡ =#

# ╔═╡ 9df66431-c596-4c79-9aa6-dc52f4ba600a
# ╠═╡ disabled = true
#=╠═╡
size(X_corrected)
  ╠═╡ =#

# ╔═╡ b6a65070-3cd6-4aca-bafc-516c84e28342
# ╠═╡ disabled = true
#=╠═╡
@assert size(X_corrected, 1) == cfg.n_channels
  ╠═╡ =#

# ╔═╡ a8cd339e-0329-4baf-ab39-170020186472
# ╠═╡ disabled = true
#=╠═╡
@assert size(X_corrected, 2) == length(times)
  ╠═╡ =#

# ╔═╡ 3529f7cf-df20-4547-b103-4dfc3fff01b1
# ╠═╡ disabled = true
#=╠═╡
@assert size(X_corrected, 3) == nrow(evts)
  ╠═╡ =#

# ╔═╡ 3054b4b7-52e0-4577-8dc5-a9ca949965ee
# ╠═╡ disabled = true
#=╠═╡
eltype(X_corrected)
  ╠═╡ =#

# ╔═╡ 07f071f5-8d0d-4a11-a2b2-6bd27a85720e
# ╠═╡ disabled = true
#=╠═╡
extrema(skipmissing(vec(X_corrected)))
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═8b8cc14c-9d99-11f1-93c4-09291a590ac4
# ╠═49496a9e-03ac-42fd-a7c0-76ee3aac3c0c
# ╠═58ed75f6-d941-4083-8e81-f1c4de3f1c03
# ╠═ffa9faea-e1e1-4ad0-956c-d16b77dd1bf2
# ╠═577347be-c791-4972-8e2c-099137d80077
# ╠═266b38ea-3fee-45dd-bf3a-5504e4177227
# ╠═4b92f32a-8b10-4a3e-9577-5ad859ee251b
# ╠═4d25fb12-c5b3-48fa-a637-5ac4f7a9fe7e
# ╠═af6d5abb-efbf-475e-9efd-765f3c234ddc
# ╠═7c2e277a-36a1-4606-acf5-71cb32c489e4
# ╠═3d5e3d62-c554-4256-9672-0042c4204a68
# ╠═2472b378-d948-40b0-98e1-08f54dd2ee40
# ╠═75df5261-b6b4-482a-aeb5-a9c4c5cb99f4
# ╠═ddbfb942-9591-40e0-bd98-75d995fcbd02
# ╠═5de2d518-5bd3-404d-81bd-9a8f0b852a01
# ╠═f7690397-cecd-4184-8da4-962d50be7d16
# ╠═a0cd698c-0a73-4918-91d5-920a3fb45177
# ╠═8a6416e9-2c7d-49fc-9d31-fbc09003adcd
# ╠═9fe10962-9032-4a3a-a91c-3dc63491d7ba
# ╠═6f7dc6f2-f4ab-4a15-9da2-94e189b2aebb
# ╠═43fe4334-da86-4f11-85c0-d50e5998fbac
# ╠═98f309fc-1d8c-4b18-a65c-1412e541b81a
# ╠═5864388e-3fb0-4723-86d2-df6662efb0d3
# ╠═f27fc2cf-411a-4544-ba94-15153f885da2
# ╠═f9eb4eec-a141-4574-bdc4-a39e9d9f36e6
# ╠═7edc508a-5147-4760-855b-5ceaf5228212
# ╠═b76b52da-3df5-43b0-89ee-51836a021c47
# ╠═434a0d07-58aa-44ea-a030-09d947b08bc8
# ╠═4bb1837a-179b-4cd1-8aa4-0bdb98d47704
# ╠═d020d01d-c910-4833-a633-320279a01253
# ╠═b8757652-0ce5-425c-a47e-40e36f37cc4f
# ╠═b0a9a532-8856-4ca4-bd7f-18314ff2f87d
# ╠═3b55d105-a21a-4604-9224-128bf3867ef8
# ╠═b86bb227-b008-41ad-ba74-faf8c18c9104
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
