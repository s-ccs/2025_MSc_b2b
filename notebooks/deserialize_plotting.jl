### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ bb713dfd-bdb9-46c8-b1eb-b0c5872d7bd6
begin
    import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 33175d2a-4fb8-4471-a5a2-61ec2411ba71
using PlutoLinks: @revise, @ingredients

# ╔═╡ d9bdc570-3f31-40df-baf8-676a19ba5960
@revise using MScB2B

# ╔═╡ e19c8aa5-5bb2-40e0-935d-cf90182ba8cd
begin
	using Serialization
	using CairoMakie
end

# ╔═╡ 35d251b4-caf8-4f68-ad99-f09d8d6172ed
result_one_step = deserialize(
        joinpath(
            @__DIR__,
            "..",
            "results",
            "final",
            "simulation_correlated_continuous",
            "one_step_b2b_1500trials_100hz.jls"
        )
    )

# ╔═╡ 317b2b1f-11dc-49f7-8b85-ff3829ac4cd8
propertynames(result_one_step)

# ╔═╡ 22e19442-4b57-4400-a047-4188e20bfd04
begin
    cfg_corr = result_one_step.cfg_corr
    sim_corr = result_one_step.sim_corr
    one_step_corr = result_one_step.one_step_corr
end

# ╔═╡ 7134dc44-d793-44f3-aee1-dec128ead46b
begin
	fig = MScB2B.plot_correlated_b2b(
	    one_step_corr.score_tables,
	    cfg_corr;
	    title = "One-step B2B on correlated continuous predictors"
	)
	
	fig
end

# ╔═╡ 065f24d7-a953-4024-9f14-88692dc74e64
# ╠═╡ disabled = true
#=╠═╡
cfg_cond = ConditionContinuousConfig()
  ╠═╡ =#

# ╔═╡ bb6e801d-7023-4a4d-a0a8-1932fc8b8597
save("04_one_step_b2b_correlatedSim.svg", fig)

# ╔═╡ d9b9b9b8-a493-45b3-902d-35acd919d57e
#=╠═╡
begin
using UnfoldSim

sfreq = cfg_cond.sfreq

n170_truth = UnfoldSim.n170(; sfreq = sfreq)
p300_truth = UnfoldSim.p300(; sfreq = sfreq)

t_n170 = (0:length(n170_truth)-1) ./ sfreq
t_p300 = (0:length(p300_truth)-1) ./ sfreq

fig = Figure(size = (900, 400))

ax = Axis(
    fig[1, 1],
    xlabel = "Time [s]",
    ylabel = "Ground-truth effect",
	xgridvisible = false,
    ygridvisible = false,
	topspinevisible = false,
	rightspinevisible = false,
)

lines!(
    ax,
	
    t_n170,
    n170_truth,
    label = "Condition → N170"
)

lines!(
    ax,
    t_p300,
    p300_truth,
    label = "Continuous → P300"
)

axislegend(ax; framevisible = false)

fig
end
  ╠═╡ =#

# ╔═╡ bd8e70cc-4861-42ae-8d62-37f7ed27d8e8
#=╠═╡
save("CondSim_ground_truth.svg", fig)
  ╠═╡ =#

# ╔═╡ a43f38b2-61d6-4a21-b3e9-e1a5abaa5324
# ╠═╡ disabled = true
#=╠═╡
cfg_corr = CorrelatedContinuousConfig()
  ╠═╡ =#

# ╔═╡ 1110114e-2954-40ea-82ca-2e43303b43dd
# ╠═╡ disabled = true
#=╠═╡
begin
	using UnfoldSim
    sfreq = cfg_corr.sfreq

    truth1 = UnfoldSim.hanning(
        cfg_corr.component_width,
        cfg_corr.peak1,
        sfreq
    )

    truth2 = UnfoldSim.hanning(
        cfg_corr.component_width,
        cfg_corr.peak2,
        sfreq
    )

    truth3 = UnfoldSim.hanning(
        cfg_corr.component_width,
        cfg_corr.peak3,
        sfreq
    )

    t1 = (0:length(truth1)-1) ./ sfreq
    t2 = (0:length(truth2)-1) ./ sfreq
    t3 = (0:length(truth3)-1) ./ sfreq

    fig = Figure(size = (1000, 400))

ax = Axis(
    fig[1, 1],
    xlabel = "Time [s]",
    ylabel = "Ground-truth effect",
    xgridvisible = false,
    ygridvisible = false,
    topspinevisible = false,
    rightspinevisible = false,
)

l1 = lines!(ax, t1, truth1, label = "Continuous 1 → 150 ms")
l2 = lines!(ax, t2, truth2, label = "Continuous 2 → 450 ms")
l3 = lines!(ax, t3, truth3, label = "Continuous 3 → 750 ms")

Legend(
    fig[1, 2],
    ax,
    framevisible = false
)

fig
    
end
  ╠═╡ =#

# ╔═╡ ceef775a-ff52-4341-9348-e623f5e43ce9
# ╠═╡ disabled = true
#=╠═╡
save("CorrSim_ground_truth.svg", fig)
  ╠═╡ =#

# ╔═╡ be2998b0-d558-4f25-bed3-095e152a6d2c
# ╠═╡ disabled = true
#=╠═╡
begin
	const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
	
	result_dir = joinpath(
	    PROJECT_ROOT,
	    "results",
	    "final",
	    "simulation_condition_continuous"
	)
end
  ╠═╡ =#

# ╔═╡ 61a1a804-8a8c-4918-a778-bfab3dc46ddd
# ╠═╡ disabled = true
#=╠═╡
standard = deserialize(
    joinpath(
        result_dir,
        "standard_decoding_results.jls"
    )
)
  ╠═╡ =#

# ╔═╡ 9be70209-4859-4893-9054-af572442bf60
# ╠═╡ disabled = true
#=╠═╡
fig_standard_condition_midterm_1500trials =
    MScB2B.plot_standard_decoding_grid(
        standard.standard_condition.scores,
        standard.cfg_cond;
        target = :condition
    )
  ╠═╡ =#

# ╔═╡ b2a15db1-aa13-4f0e-a648-e30c2f1a98e6
# ╠═╡ disabled = true
#=╠═╡
fig_standard_continuous_midterm_1500trials =
    MScB2B.plot_standard_decoding_grid(
        standard.standard_continuous.scores,
        standard.cfg_cond;
        target = :continuous
    )
  ╠═╡ =#

# ╔═╡ 340b8ba9-b8d2-4f53-b928-07744075bfc2
# ╠═╡ disabled = true
#=╠═╡
begin
	save("fig_correlatedSim_standard_condition_midterm_1500trials.svg",fig_standard_condition_midterm_1500trials)

	save("fig_correlatedSim_standard_continuous_midterm_1500trials.svg",fig_standard_continuous_midterm_1500trials)
end
  ╠═╡ =#

# ╔═╡ 8df82456-15ab-470b-97b3-d7acd53ada91
# ╠═╡ disabled = true
#=╠═╡
rerp = deserialize(
    joinpath(
        result_dir,
        "rerp_decoding_results.jls"
    )
)
  ╠═╡ =#

# ╔═╡ b465428f-9dca-42e4-9533-22a1fa9b53c4
# ╠═╡ disabled = true
#=╠═╡
(
    standard = (
        n_trials = standard.cfg_cond.n_trials,
        sfreq = standard.cfg_cond.sfreq
    ),

    rerp = (
        n_trials = rerp.cfg.n_trials,
        sfreq = rerp.cfg.sfreq
    )
)
  ╠═╡ =#

# ╔═╡ fb95b5e3-a4a3-46da-9442-f79c336bcb47
# ╠═╡ disabled = true
#=╠═╡
fig_rerp_condition_midterm_1500trials =
    MScB2B.plot_standard_decoding_grid(
        rerp.condition_results.scores,
        rerp.cfg;
        target = :condition
    )
  ╠═╡ =#

# ╔═╡ c06f4920-0963-461b-9655-ca84091c1d9d
# ╠═╡ disabled = true
#=╠═╡
fig_rerp_continuous_midterm_1500trials =
    MScB2B.plot_standard_decoding_grid(
        rerp.continuous_results.scores,
        rerp.cfg;
        target = :continuous
    )
  ╠═╡ =#

# ╔═╡ 47d570cf-fd20-4ee0-a3b2-067d452b9cb0
# ╠═╡ disabled = true
#=╠═╡
begin
	save("fig_correlatedSim_rerp_condition_midterm_1500trials.svg",fig_rerp_condition_midterm_1500trials)

	save("fig_correlatedSim_rerp_continuous_midterm_1500trials.svg",fig_rerp_continuous_midterm_1500trials)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═bb713dfd-bdb9-46c8-b1eb-b0c5872d7bd6
# ╠═33175d2a-4fb8-4471-a5a2-61ec2411ba71
# ╠═d9bdc570-3f31-40df-baf8-676a19ba5960
# ╠═e19c8aa5-5bb2-40e0-935d-cf90182ba8cd
# ╠═35d251b4-caf8-4f68-ad99-f09d8d6172ed
# ╠═317b2b1f-11dc-49f7-8b85-ff3829ac4cd8
# ╠═22e19442-4b57-4400-a047-4188e20bfd04
# ╠═7134dc44-d793-44f3-aee1-dec128ead46b
# ╠═065f24d7-a953-4024-9f14-88692dc74e64
# ╠═bb6e801d-7023-4a4d-a0a8-1932fc8b8597
# ╠═d9b9b9b8-a493-45b3-902d-35acd919d57e
# ╠═bd8e70cc-4861-42ae-8d62-37f7ed27d8e8
# ╠═a43f38b2-61d6-4a21-b3e9-e1a5abaa5324
# ╠═1110114e-2954-40ea-82ca-2e43303b43dd
# ╠═ceef775a-ff52-4341-9348-e623f5e43ce9
# ╠═be2998b0-d558-4f25-bed3-095e152a6d2c
# ╠═61a1a804-8a8c-4918-a778-bfab3dc46ddd
# ╠═9be70209-4859-4893-9054-af572442bf60
# ╠═b2a15db1-aa13-4f0e-a648-e30c2f1a98e6
# ╠═340b8ba9-b8d2-4f53-b928-07744075bfc2
# ╠═8df82456-15ab-470b-97b3-d7acd53ada91
# ╠═b465428f-9dca-42e4-9533-22a1fa9b53c4
# ╠═fb95b5e3-a4a3-46da-9442-f79c336bcb47
# ╠═c06f4920-0963-461b-9655-ca84091c1d9d
# ╠═47d570cf-fd20-4ee0-a3b2-067d452b9cb0
