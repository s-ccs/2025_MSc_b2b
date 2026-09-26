### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ f52ff27a-b5d3-11f1-9927-553946ff73b0
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ c1cd1ec3-6ce3-4b28-aeac-5e467926680e
using PlutoLinks: @revise, @ingredients

# ╔═╡ fcfce135-63a3-4329-94ce-7fda78c05c88
@revise using MScB2B

# ╔═╡ f810c2bd-c0f1-479c-8fe3-1f978f1de913
begin
    using CSV
    using DataFrames
    using Statistics
    using NPZ
    using CairoMakie
    using Unfold
	using UnfoldMakie
end

# ╔═╡ d2e282fb-f449-4714-8505-ab404ae6b239
begin
	DATA_DIR = "/scratch/data/ROAMM/outputs/eeg"
	PREPROC_DIR = "/scratch/data/ROAMM/outputs/preprocessing"
end

# ╔═╡ f832dbcf-29a4-4940-950b-208f80d8d3e6
EEG_CHANNELS = [
    "Fp1", "AF7", "AF3", "F1", "F3", "F5", "F7", "FT7",
    "FC5", "FC3", "FC1", "C1", "C3", "C5", "T7", "TP7",
    "CP5", "CP3", "CP1", "P1", "P3", "P5", "P7", "P9",
    "PO7", "PO3", "O1", "Iz", "Oz", "POz", "Pz", "CPz",
    "Fpz", "Fp2", "AF8", "AF4", "Afz", "Fz", "F2", "F4",
    "F6", "F8", "FT8", "FC6", "FC4", "FC2", "FCz", "Cz",
    "C2", "C4", "C6", "T8", "TP8", "CP6", "CP4", "CP2",
    "P2", "P4", "P6", "P8", "P10", "PO8", "PO4", "O2",
]

# ╔═╡ 1adc39fc-2d26-427b-9b87-39605b692b03
begin
	SUBJECT = "sub-10177"
	RUN = 3
	SFREQ = 256
end

# ╔═╡ 05303ca7-8c71-4474-bd74-96c9266f5ae5
begin
	eeg_file = joinpath(DATA_DIR, "$(SUBJECT)_run$(RUN)_eeg.npy")
	events_file = joinpath(PREPROC_DIR, "roamm_b2b_events.csv")
	all_events = CSV.read(events_file, DataFrame)
	dat = NPZ.npzread(eeg_file)
	evts = all_events[
		(all_events.subject_id .== SUBJECT) .&
		(all_events.run_num .== RUN), 
		:
	]
end

# ╔═╡ c7148930-cb7c-4d4a-89c9-23569e2c2c75
size(dat)

# ╔═╡ 7d8d5dfe-2184-4d32-8871-c47d8276c41f
size(evts)

# ╔═╡ 59abb5b9-a2dc-4e3d-b43e-91db46d0e964
first(evts, 5)

# ╔═╡ fce46637-e9fb-4c0d-b1a6-f2e1d6f9bfdd
extrema(evts.latency)

# ╔═╡ 711ef9be-66dd-4a9b-85fe-58ab69396fbc
size(dat, 2)

# ╔═╡ ec688884-3abd-4fb7-af69-13b4740555af
maximum(evts.latency) <= size(dat, 2)

# ╔═╡ ed70adf3-f076-4304-9ef9-1da20968dd9b
begin
	dat_uv = dat .* 1e6
	extrema(dat_uv)
end

# ╔═╡ 3a2ff9a7-952b-4e19-b578-0b3c181e6714
md"""
# 1. Raw EEG quick check"""

# ╔═╡ 852348ff-dabc-428f-b107-4901a01b4b70
begin
	oz_idx = findfirst(==("Oz"), EEG_CHANNELS)
	let
	    ix = 5000:6280 # 5 seconds at 256Hz
	
	    event_ix = findall(
	        (evts.latency .>= first(ix)) .&
	        (evts.latency .<= last(ix))
	    )
	
	    f = Figure()
	    ax = Axis(
	        f[1,1],
	        xlabel = "Sample",
	        ylabel = "Amplitude (μV)",
	        title = "$(SUBJECT), run $(RUN), channel Oz"
	    )
	
	    lines!(ax, ix, dat_uv[oz_idx, ix])
	    vlines!(ax, evts.latency[event_ix])
	
	    f
		save("raw_eeg_quick_check.svg", f)
	end
	
end

# ╔═╡ cff33dc9-01d5-49ab-831c-433d42e9a6db
let 
	ix = 4000:10000
	lines(ix, dat_uv[29, ix])
	vlines!(evts.latency[30:45])
	current_figure()
end

# ╔═╡ 2e23406e-3326-4e6e-b408-bb1812c3924a
md"""
# 2. Fixation timing sanity check
"""

# ╔═╡ ba344a6c-87b6-4294-9986-3034f391b492
hist(evts.duration_ms, bins = 100)

# ╔═╡ 3b55cfa6-c906-41c9-8021-a933c07d3d51
begin
	
	let
		f = hist(
		    evts.duration_ms,
		    bins = 100,
		    axis = (
		        xlabel = "Fixation duration [ms]",
		        ylabel = "Count",
		        title = "$(SUBJECT), run $(RUN)"
		    )
		)
		
		save("fixation_duration_hist.svg", f)
	end
end

# ╔═╡ 20d89361-5f71-45bd-b2dd-4072d8d99b87
md"""
# 3. All-channel FRP butterfly
"""

# ╔═╡ c8b9305a-5d5e-44d1-97fd-a26cc6dd95af
begin
	# epoch
	dat_e, times_e = Unfold.epoch(dat_uv, evts, [-0.2, 0.8], SFREQ)

	# raw FRP
	frp = dropdims(mean(dat_e, dims=3), dims=3)
end

# ╔═╡ a0bf8d24-1c6b-438a-b8c8-1898d87ab90a
begin
	series(frp, solid_color=:black)
	f_series = current_figure()
	save("frp_butterfly.svg", f_series)
	f_series
end

# ╔═╡ cbb0b283-0cf3-46ba-bd3d-152c1c911c82
begin
	heatmap(frp')
	f_heatmap = current_figure()
	save("frp_heatmap.svg", f_heatmap)
	f_heatmap
end

# ╔═╡ 70f3f436-62a0-4cc3-88ae-541aa36c2997
md"""
# 4. Posterior FRP"""

# ╔═╡ aaebce5e-7848-4856-8622-d90e089556c7
frp_df = DataFrame(
	time = repeat(collect(times_e), inner=length(EEG_CHANNELS)),
	channel = repeat(EEG_CHANNELS, length(times_e)),
	estimate = vec(frp)
)

# ╔═╡ 72cdb441-3662-4ca4-8874-0443e2c72722
begin
	posterior_channels = ["O1", "Oz", "O2", "POz"]
	posterior_df = frp_df[in.(frp_df.channel, Ref(posterior_channels)),:]

	posterior_idx = [findfirst(==(ch), EEG_CHANNELS) for ch in posterior_channels]
end

# ╔═╡ f1b1c379-85bd-4ccd-b099-c5522957148a
begin
	Posterior_FRP= plot_erp(
	    posterior_df;
	    mapping = (
	        color = :channel,
	    ),
	    axis = (
	        xlabel = "Time from fixation onset [s]",
	        ylabel = "Amplitude [μV]",
	        title = "$(SUBJECT), run $(RUN): raw FRP",
	    )
	)

	save("Posterior_FRP.svg", Posterior_FRP)
end

# ╔═╡ 6c014db7-4522-467f-8612-a56d68c167a3
md"""
# 5. Raw FRP vs deconvolved FRP
## intercept-only
"""

# ╔═╡ a759d4db-fb7f-460f-9717-4b7bd62512eb
m_epoched = fit(
    UnfoldModel,
    @formula(0 ~ 1),
    evts,
    dat_e,
    times_e
)

# ╔═╡ 5a03f494-a183-421a-aa14-74950094a61d
epoched_frp = dropdims(
    coef(m_epoched),
    dims = 3
)

# ╔═╡ a7a69670-e174-4055-848f-76577ea73a09
series(frp, solid_color = :black)

# ╔═╡ 45d30899-d7b3-4078-9618-42a2800f1b9d
series(
    epoched_frp,
    solid_color = :black
)

# ╔═╡ 54155cf9-bd4d-4cde-9269-5d9f60c8adb0
m_deconv = fit(
	UnfoldModel,
	@formula(0~1),
	evts,
	dat_uv,
	firbasis(
		sfreq = 256,
		τ = (-0.2, 0.8)
	)
)

# ╔═╡ 312d022e-5d24-4ae7-bee9-a291ecc8288b
deconv_frp = coef(m_deconv)

# ╔═╡ e2b6e0a4-5858-40b6-a600-313256ecaea6
series(deconv_frp, solid_color=:black)

# ╔═╡ 0b596b30-1d57-4092-b3a9-abaad2c59d63
begin
	raw_post = vec(mean(frp[posterior_idx, :], dims=1))
	epoched_post = vec(mean(epoched_frp[posterior_idx, :], dims=1))
	deconv_post = vec(mean(deconv_frp[posterior_idx, :], dims=1))
end

# ╔═╡ 606051e1-9f3c-47f9-893b-0af1fba24c6e
heatmap(deconv_frp')

# ╔═╡ b30a76ac-f64d-45b5-845d-18ebb95df9e6
let
    lines(
        collect(times_e),
        raw_post,
        label = "Raw mean FRP",
		linestyle = :dash
    )

	lines!(
        collect(times_e),
        epoched_post,
        label = "Epoched",
		linestyle = :dot
    )
	
    lines!(
        collect(times_e),
        deconv_post,
        label = "Deconvolved FRP"
    )

    axislegend()
    current_figure()
end

# ╔═╡ 3d798e42-0160-43d1-be89-c3d7a762d793


# ╔═╡ 6e2f6e5c-0e1e-40f5-ac47-33454b94d69b
md"""
# 6. Two-step B2B
## One subject, one run
"""

# ╔═╡ 1807f9ce-3e52-4274-bda4-240325715092
# ╠═╡ disabled = true
#=╠═╡
fit_length = MScB2B.fit_two_step_b2b_roamm(
	dat_uv,
	evts;
	predictors = [:word_length],
	cross_val_reps = 20
)
  ╠═╡ =#

# ╔═╡ dc6f7f07-34c1-4bae-b81f-29530416d1bc
# ╠═╡ disabled = true
#=╠═╡
begin
	result_length = DataFrame(
    Unfold.coeftable(fit_length.b2b_fit_corrected)
)

	result_length = result_length[
    result_length.coefname .!= "(Intercept)",
    :
]

	result_length[!, :estimate_signed] = copy(result_length.estimate)
	result_length.estimate .= abs.(result_length.estimate)
											 
end
  ╠═╡ =#

# ╔═╡ c1950600-59f5-476f-9877-8395a17d7e3c
#=╠═╡
one_subject_one_run_roamm = plot_erp(result_length)
  ╠═╡ =#

# ╔═╡ 32ea42a5-7f7c-4328-a630-5b4ffc2b5fc2
#=╠═╡
save("one_subject_one_run_roamm.svg", one_subject_one_run_roamm)
  ╠═╡ =#

# ╔═╡ 5c7984af-fb9b-42d0-88a9-33241af6c123
# ╠═╡ disabled = true
#=╠═╡
predictors = [
		:word_length,
		:frequency_zipf,
		:word_surprisal
	]
  ╠═╡ =#

# ╔═╡ ac6f6215-f719-4b8a-89a5-3dba3756525a
# ╠═╡ disabled = true
#=╠═╡
fit_all = MScB2B.fit_two_step_b2b_roamm(
	dat_uv,
	evts;
	predictors
)
  ╠═╡ =#

# ╔═╡ 42180da3-93eb-4d45-94e4-006d1bef3326
# ╠═╡ disabled = true
#=╠═╡
begin
	result_all = DataFrame(
    Unfold.coeftable(fit_all.b2b_fit_corrected)
)

#	result_all = result_all[
#   result_all.coefname .!= "(Intercept)",
#    :
#]

	result_all[!, :estimate_signed] = copy(result_all.estimate)
	result_all.estimate .= abs.(result_all.estimate)
											 
end
  ╠═╡ =#

# ╔═╡ 81f5bd9f-51d6-43d7-aae5-cb923e69d713
#=╠═╡
plot_erp(result_all, mapping=(color = :coefname,))
  ╠═╡ =#

# ╔═╡ d15064ee-9057-4638-843e-f92b72b83ccb
md"""
## Subject-level two-step b2b
"""

# ╔═╡ 857398fb-6ad9-4f6e-b670-af8f8bfe1baf
# ╠═╡ disabled = true
#=╠═╡
result_subject = CSV.read(
    "/store/users/xu/2025_MSc_b2b/results/roamm_b2b/sub-10177_word_length_b2b.csv",
    DataFrame
)
  ╠═╡ =#

# ╔═╡ 5e10757b-a882-4c6c-b22b-0650197d984b
# ╠═╡ disabled = true
#=╠═╡
one_subject_roamm = plot_erp(
    result_subject;
    axis = (
        xlabel = "Time [s]",
        ylabel = "B2B estimate",
        title = "sub-10177: word length",
		ytickformat = "{:.1e}"
    )
)
  ╠═╡ =#

# ╔═╡ e3ec59c4-0a3e-4e2a-b6ac-75541d03ffde
# ╠═╡ disabled = true
#=╠═╡
save("one_subject_roamm.svg", one_subject_roamm)
  ╠═╡ =#

# ╔═╡ 6611800c-cb5f-4906-8d22-a0da635e0bb6
#=╠═╡
extrema(result_length.estimate)
  ╠═╡ =#

# ╔═╡ 9df6a51e-0910-4677-b3a0-f3153d5a71f3
#=╠═╡
extrema(result_subject.estimate)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═f52ff27a-b5d3-11f1-9927-553946ff73b0
# ╠═c1cd1ec3-6ce3-4b28-aeac-5e467926680e
# ╠═fcfce135-63a3-4329-94ce-7fda78c05c88
# ╠═f810c2bd-c0f1-479c-8fe3-1f978f1de913
# ╠═d2e282fb-f449-4714-8505-ab404ae6b239
# ╠═f832dbcf-29a4-4940-950b-208f80d8d3e6
# ╠═1adc39fc-2d26-427b-9b87-39605b692b03
# ╠═05303ca7-8c71-4474-bd74-96c9266f5ae5
# ╠═c7148930-cb7c-4d4a-89c9-23569e2c2c75
# ╠═7d8d5dfe-2184-4d32-8871-c47d8276c41f
# ╠═59abb5b9-a2dc-4e3d-b43e-91db46d0e964
# ╠═fce46637-e9fb-4c0d-b1a6-f2e1d6f9bfdd
# ╠═711ef9be-66dd-4a9b-85fe-58ab69396fbc
# ╠═ec688884-3abd-4fb7-af69-13b4740555af
# ╠═ed70adf3-f076-4304-9ef9-1da20968dd9b
# ╟─3a2ff9a7-952b-4e19-b578-0b3c181e6714
# ╠═852348ff-dabc-428f-b107-4901a01b4b70
# ╠═cff33dc9-01d5-49ab-831c-433d42e9a6db
# ╟─2e23406e-3326-4e6e-b408-bb1812c3924a
# ╠═ba344a6c-87b6-4294-9986-3034f391b492
# ╠═3b55cfa6-c906-41c9-8021-a933c07d3d51
# ╟─20d89361-5f71-45bd-b2dd-4072d8d99b87
# ╠═c8b9305a-5d5e-44d1-97fd-a26cc6dd95af
# ╠═a0bf8d24-1c6b-438a-b8c8-1898d87ab90a
# ╠═cbb0b283-0cf3-46ba-bd3d-152c1c911c82
# ╟─70f3f436-62a0-4cc3-88ae-541aa36c2997
# ╠═aaebce5e-7848-4856-8622-d90e089556c7
# ╠═72cdb441-3662-4ca4-8874-0443e2c72722
# ╠═f1b1c379-85bd-4ccd-b099-c5522957148a
# ╟─6c014db7-4522-467f-8612-a56d68c167a3
# ╠═a759d4db-fb7f-460f-9717-4b7bd62512eb
# ╠═5a03f494-a183-421a-aa14-74950094a61d
# ╠═a7a69670-e174-4055-848f-76577ea73a09
# ╠═45d30899-d7b3-4078-9618-42a2800f1b9d
# ╠═e2b6e0a4-5858-40b6-a600-313256ecaea6
# ╠═54155cf9-bd4d-4cde-9269-5d9f60c8adb0
# ╠═312d022e-5d24-4ae7-bee9-a291ecc8288b
# ╠═0b596b30-1d57-4092-b3a9-abaad2c59d63
# ╠═606051e1-9f3c-47f9-893b-0af1fba24c6e
# ╠═b30a76ac-f64d-45b5-845d-18ebb95df9e6
# ╠═3d798e42-0160-43d1-be89-c3d7a762d793
# ╟─6e2f6e5c-0e1e-40f5-ac47-33454b94d69b
# ╠═1807f9ce-3e52-4274-bda4-240325715092
# ╠═dc6f7f07-34c1-4bae-b81f-29530416d1bc
# ╠═c1950600-59f5-476f-9877-8395a17d7e3c
# ╠═32ea42a5-7f7c-4328-a630-5b4ffc2b5fc2
# ╠═5c7984af-fb9b-42d0-88a9-33241af6c123
# ╠═ac6f6215-f719-4b8a-89a5-3dba3756525a
# ╠═42180da3-93eb-4d45-94e4-006d1bef3326
# ╠═81f5bd9f-51d6-43d7-aae5-cb923e69d713
# ╟─d15064ee-9057-4638-843e-f92b72b83ccb
# ╠═857398fb-6ad9-4f6e-b670-af8f8bfe1baf
# ╠═5e10757b-a882-4c6c-b22b-0650197d984b
# ╠═e3ec59c4-0a3e-4e2a-b6ac-75541d03ffde
# ╠═6611800c-cb5f-4906-8d22-a0da635e0bb6
# ╠═9df6a51e-0910-4677-b3a0-f3153d5a71f3
