### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 02ca7078-b808-11f1-bc85-313c1753f1c7
begin
    import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ d413a7cd-ac08-41b0-9269-c3419ffa889f
using PlutoLinks: @revise, @ingredients

# ╔═╡ 7e973892-038f-4b07-ab9e-3bb7a69ec90c
@revise using MScB2B

# ╔═╡ b28fa550-cbe6-42fa-ab30-a722ae038c01
begin
	using CairoMakie
	using DataFrames
	using Statistics
	using Serialization
end

# ╔═╡ df497ab9-b21e-400c-b7a0-ca6254ead1ed
using UnfoldSim


# ╔═╡ 1698ae93-5479-4132-b604-07e61299335e
begin
    using LinearAlgebra
    import Distributions
    using Random
    

    Σ = [
        1.0  0.8  0.3
        0.8  1.0  0.5
        0.3  0.5  1.0
    ]

    println("symmetric: ", issymmetric(Σ))
    println("positive definite: ", isposdef(Symmetric(Σ)))
    println("eigenvalues: ", eigvals(Symmetric(Σ)))

    dist = Distributions.MvNormal(zeros(3), Σ)

    rng = MersenneTwister(12)
    X = rand(rng, dist, 10_000)'   # observations × predictors

    println("sample correlations:")
    display(cor(X))
end

# ╔═╡ 2aad1cc2-f3b3-4c40-8ead-f6ccb415f1c0
begin
	# Output directory for this run
	const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
	const PLOT_DIR = joinpath(PROJECT_ROOT, "plots", "2026-09-27-Corr_Sim-rerun-with-correct-no-overlap")
end

# ╔═╡ 7bc141e7-6f73-45e7-83dc-a5e0495b8876
# ╠═╡ disabled = true
#=╠═╡
begin
	results_corr = deserialize("results/final/simulation_correlated_continuous/correlated_continuous_results.jls")

	cfg_corr = results_corr.cfg_corr
	sim_corr = results_corr.sim_corr
	standard_corr = results_corr.standard_corr
	rerp_corr = results_corr.rerp_corr
	plain_b2b_corr = results_corr.plain_b2b_corr
	two_step_corr = results_corr.two_step_corr
end
  ╠═╡ =#

# ╔═╡ 1f7b0c0c-a07d-41e9-bee0-5bb44500baf5
Controls = @ingredients(joinpath(@__DIR__, "simulation_controls.jl"))

# ╔═╡ fbaef126-7bed-46cd-98ab-9fc66b798c81
cfg_corr = CorrelatedContinuousConfig(
    n_trials = 1500,
    sfreq = 100.0,
	
  	rho12 = 0.8,
    rho13 = 0.8,
    rho23 = 0.8,

    noiselevel = 0.3,
    channel_noise_sd = 0.3,

	
	component_width = 0.15,
    peak1 = 0.15,
    peak2 = 0.45,
    peak3 = 0.75,

    β0 = 5.0,
    β = 2.0,
	
    #β0_p100 = 5.0,
    #β_p100 = 2.0,

    #β0_n170 = 5.0,
    #β_n170 = 3.0,

    #β0_p300 = 5.0,
    #β_p300 = 1.5,

    overlap_interval_ms = 250.0,
    onset_predictor_bias = -0.6,
    shift_onset = true,
)

# ╔═╡ 9b22a44d-60fa-4101-bdf3-01f3fa211e54
begin
	old_max = maximum([
	    length(UnfoldSim.p100(; sfreq = cfg_corr.sfreq)),
	    length(UnfoldSim.n170(; sfreq = cfg_corr.sfreq)),
	    length(UnfoldSim.p300(; sfreq = cfg_corr.sfreq))
	])
	
	new_max = maximum([
	    length(UnfoldSim.hanning(
	        cfg_corr.component_width,
	        cfg_corr.peak1,
	        cfg_corr.sfreq
	    )),
	    length(UnfoldSim.hanning(
	        cfg_corr.component_width,
	        cfg_corr.peak2,
	        cfg_corr.sfreq
	    )),
	    length(UnfoldSim.hanning(
	        cfg_corr.component_width,
	        cfg_corr.peak3,
	        cfg_corr.sfreq
	    ))
	])
	


	(old_max = old_max, new_max = new_max, safe = old_max >= new_max)
end

# ╔═╡ 0f02aff9-3dfc-421d-8cc6-e32e623bc0c2
sim_corr = simulate_corr_cont_cases(
    cfg_corr;
    seed = 12
)

# ╔═╡ cc11676d-b434-476f-ae24-6c35c4fa0333
# ╠═╡ disabled = true
#=╠═╡
begin
    X = Matrix(
        sim_corr.overlap.events_continuous[
            :,
            [:continuous1, :continuous2, :continuous3]
        ]
    )

    cor(X)
end
  ╠═╡ =#

# ╔═╡ ebf10fc0-f163-4fab-9565-2561c058c4fd
ridge_model = MScB2B.make_ridge_tuned_model(
    inner_nfolds = 2,
    resolution = 4
)

# ╔═╡ fc969ae8-fc8f-43d0-af80-959c481c30e8
begin
	ev = sim_corr.overlap.events_continuous
	iei = diff(ev.latency) ./ cfg_corr.sfreq
	x = ev.continuous1[1:end-1]
	
	cor(x, iei)
end

# ╔═╡ 235b519d-cf95-46ad-b6d1-c2a399a9bab6
dis = Distributions.MvNormal(zeros(3), Σ)

# ╔═╡ 9d97eb5d-b208-4e3f-8a82-d63caac9c1d7


# ╔═╡ 6688cf47-a90f-476d-9daf-1a9ba44e5f65
md"""
# Pipeline 01: Standard decoding
"""

# ╔═╡ 6abb536f-797d-4059-9503-0b71e78df713
standard_corr = 
	run_standard_decoding(
		cfg_corr,
		sim_corr;
		nfolds = 3
	)

# ╔═╡ eb71a409-b18c-4cf1-bbef-69f15205adc7
# ╠═╡ disabled = true
#=╠═╡
MScB2B.plot_erp(
	standard_corr.scores,
	mapping = (
		x = :time,
		y = :r,
		color = :case,
		layout = :target
	)
)
  ╠═╡ =#

# ╔═╡ eb968238-6113-441a-87f0-e67e058fc9dc
standard_decoding_correlatedSim = plot_correlated_decoding(
    standard_corr.scores,
    cfg_corr;
    title = "Standard decoding on correlated continuous predictors"
)

# ╔═╡ eaba59e5-76f4-41df-9727-f6d1e900d52d
# ╠═╡ disabled = true
#=╠═╡
begin
	outdir = joinpath(
		"results",
		"experiments",
		"simulation_correlated_continuous",
		"rho_sensitivity"
	)
	
	makpath(outdir)
	
	serialize(joinpath(outdir, "rho_0.8_0.5_0.3.jl"), rho_results)
end
  ╠═╡ =#

# ╔═╡ bcb8dd5a-e8bc-4b66-aafb-e2db9bfee107
# ╠═╡ disabled = true
#=╠═╡
begin
	ev = sim_corr.overlap.events_continuous
	
	iei = diff(ev.latency) ./ cfg_corr.sfreq
	x = ev.continuous1[1:end-1]
	
	cor(x, iei)
end
  ╠═╡ =#

# ╔═╡ 1d9076b0-332d-443a-bfd8-075b73a0d76d
# ╠═╡ disabled = true
#=╠═╡
begin
	cut = median(x)
	
	mean(iei[x .<= cut]),
	mean(iei[x .> cut])
end
  ╠═╡ =#

# ╔═╡ aa9eab56-9228-4fb6-958b-86393d84e206
md"""
# Pipeline 02: rERP decoding
"""

# ╔═╡ cd48cfd2-637a-437f-afc3-c0fc143af3b4
rerp_corr = 
	run_rerp_decoding(
		cfg_corr,
		sim_corr;
		nfolds = 3
	)

# ╔═╡ 85f152de-865b-4e8a-b6e8-8f06a88bc5bf
rerp_decoding_correlatedSim = plot_correlated_decoding(
    rerp_corr.scores,
    cfg_corr;
    title = "rERP decoding on correlated continuous predictors"
)

# ╔═╡ 07e71057-6107-48cc-a78e-c8746579271c
# ╠═╡ disabled = true
#=╠═╡
MScB2B.plot_erp(
	rerp_corr.scores,
	mapping = (
		x = :time,
		y = :r,
		color = :case,
		layout = :target
	)
)
  ╠═╡ =#

# ╔═╡ d15ff62c-b202-405a-8c22-0c88d5dd4071
md"""
# Pipeline 03: Plain B2B
"""

# ╔═╡ 1f277c0b-257e-4c68-81fc-9506fda5c43c
plain_b2b_corr = 
	run_plain_b2b(
		cfg_corr,
		sim_corr;
		cross_val_reps = 3
	)

# ╔═╡ 5156ed01-c13b-4c45-b985-e87f7ba30266
# ╠═╡ disabled = true
#=╠═╡
MScB2B.plot_erp(
	plain_b2b_corr.score_tables,
	mapping = (
		x = :time,
		y = :estimate,
		color = :case,
		layout = :coefname
	)
)
  ╠═╡ =#

# ╔═╡ b5b217df-e188-4430-8fc3-fe9cdd807bd4
plain_b2b_correlatedSim = MScB2B.plot_correlated_b2b(
    plain_b2b_corr.score_tables,
    cfg_corr;
    title = "Plain B2B on correlated continuous predictors"
)

# ╔═╡ 16e47145-992c-48e5-b9ef-e92505770cd1
md"""
# Pipeline 05: two step B2B
"""

# ╔═╡ 6ed849df-4c53-4162-9806-15a203652976
two_step_corr = 
	run_two_step_b2b(
		cfg_corr,
		sim_corr;
		cross_val_reps = 3
	)

# ╔═╡ 73d3520f-85b7-4b24-9dfd-25e6cbd63fea
begin
	result_dir = joinpath(
		"results",
		"simulation_correlated_continuous"
	)
	mkpath(result_dir)
	
	serialize(
		joinpath(result_dir, "rerun_with_correct_no_overlap_corrSim_results.jls"),
		(
			cfg_corr = cfg_corr,
			sim_corr = sim_corr,
			standard_corr = standard_corr,
			rerp_corr = rerp_corr,
			plain_b2b_corr = plain_b2b_corr,
			two_step_b2b_corr = two_step_corr
		)
	)
end

# ╔═╡ c2df51c4-0c44-4aa1-805a-2cb32c7bbfb6
# ╠═╡ disabled = true
#=╠═╡
two_step_b2b_correlatedSim = MScB2B.plot_erp(
	two_step_corr.score_tables,
	mapping = (
		x = :time,
		y = :estimate,
		color = :case,
		layout = :coefname
	)
)
  ╠═╡ =#

# ╔═╡ 7a338232-cd7f-478c-901d-163b5ebf7c81
two_step_b2b_correlatedSim = MScB2B.plot_correlated_b2b(
    two_step_corr.score_tables,
    cfg_corr;
    title = "Two-step B2B on correlated continuous predictors"
)

# ╔═╡ 515db147-a448-4ccf-8e93-24ba5513cac5
# ╠═╡ disabled = true
#=╠═╡
begin
using AlgebraOfGraphics


fig_rerp_corr =
    data(two_step_corr.score_tables) *
    mapping(
        :time,
        :estimate,
        color = :case,
        layout = :coefname
    ) *
    visual(Lines) |>
    draw
end
  ╠═╡ =#

# ╔═╡ 675a83ea-4160-46f2-91b3-2d44c6d1e336
begin
	save("03_plain_b2b_correlatedSim_correct_no_overlap.svg", plain_b2b_correlatedSim)
	save("05_two_step_b2b_correlatedSim_correct_no_overlap.svg", two_step_b2b_correlatedSim)
	save("01_standard_decoding_correlatedSim_correct_no_overlap.svg", standard_decoding_correlatedSim)
	save("02_rerp_decoding_correlatedSim_correct_no_overlap.svg", rerp_decoding_correlatedSim)

end

# ╔═╡ Cell order:
# ╠═02ca7078-b808-11f1-bc85-313c1753f1c7
# ╠═d413a7cd-ac08-41b0-9269-c3419ffa889f
# ╠═7e973892-038f-4b07-ab9e-3bb7a69ec90c
# ╠═2aad1cc2-f3b3-4c40-8ead-f6ccb415f1c0
# ╠═7bc141e7-6f73-45e7-83dc-a5e0495b8876
# ╠═b28fa550-cbe6-42fa-ab30-a722ae038c01
# ╠═73d3520f-85b7-4b24-9dfd-25e6cbd63fea
# ╠═df497ab9-b21e-400c-b7a0-ca6254ead1ed
# ╠═9b22a44d-60fa-4101-bdf3-01f3fa211e54
# ╠═1f7b0c0c-a07d-41e9-bee0-5bb44500baf5
# ╠═fbaef126-7bed-46cd-98ab-9fc66b798c81
# ╠═0f02aff9-3dfc-421d-8cc6-e32e623bc0c2
# ╠═cc11676d-b434-476f-ae24-6c35c4fa0333
# ╠═ebf10fc0-f163-4fab-9565-2561c058c4fd
# ╠═fc969ae8-fc8f-43d0-af80-959c481c30e8
# ╠═1698ae93-5479-4132-b604-07e61299335e
# ╠═235b519d-cf95-46ad-b6d1-c2a399a9bab6
# ╠═9d97eb5d-b208-4e3f-8a82-d63caac9c1d7
# ╟─6688cf47-a90f-476d-9daf-1a9ba44e5f65
# ╠═6abb536f-797d-4059-9503-0b71e78df713
# ╠═eb71a409-b18c-4cf1-bbef-69f15205adc7
# ╠═eb968238-6113-441a-87f0-e67e058fc9dc
# ╠═eaba59e5-76f4-41df-9727-f6d1e900d52d
# ╠═bcb8dd5a-e8bc-4b66-aafb-e2db9bfee107
# ╠═1d9076b0-332d-443a-bfd8-075b73a0d76d
# ╟─aa9eab56-9228-4fb6-958b-86393d84e206
# ╠═cd48cfd2-637a-437f-afc3-c0fc143af3b4
# ╠═85f152de-865b-4e8a-b6e8-8f06a88bc5bf
# ╠═07e71057-6107-48cc-a78e-c8746579271c
# ╟─d15ff62c-b202-405a-8c22-0c88d5dd4071
# ╠═1f277c0b-257e-4c68-81fc-9506fda5c43c
# ╠═5156ed01-c13b-4c45-b985-e87f7ba30266
# ╠═b5b217df-e188-4430-8fc3-fe9cdd807bd4
# ╟─16e47145-992c-48e5-b9ef-e92505770cd1
# ╠═6ed849df-4c53-4162-9806-15a203652976
# ╠═c2df51c4-0c44-4aa1-805a-2cb32c7bbfb6
# ╠═7a338232-cd7f-478c-901d-163b5ebf7c81
# ╠═515db147-a448-4ccf-8e93-24ba5513cac5
# ╠═675a83ea-4160-46f2-91b3-2d44c6d1e336
