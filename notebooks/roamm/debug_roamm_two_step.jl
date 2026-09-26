### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 3ee335c0-abc6-11f1-90c5-05876e25c543
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 7d609d2d-6e9b-4350-8cdb-52b3a4ded147
using PlutoLinks: @revise, @ingredients

# ╔═╡ 6886a972-e132-46e4-a2f9-530ce4dcdbeb
@revise using MScB2B

# ╔═╡ c64c8d29-d53c-4846-adc2-64313139e761
begin
	using CairoMakie
	using Statistics
	using NPZ
	using CSV
	using DataFrames

	using Unfold
	using UnfoldDecode
	using StatsModels: @formula
end

# ╔═╡ 3d7235bf-62ec-4cd0-baa4-32c1220f5983
begin
	DATA_DIR = "/scratch/data/ROAMM/outputs"
	eeg_file = joinpath(DATA_DIR, "sub-10014_run1_eeg.npy")
	events_files = joinpath(DATA_DIR, "sub-10014_run1_fixation_events.csv")

	dat = npzread(eeg_file)
	evts = CSV.read(events_files, DataFrame)
	evts = filter(:eye => ==("L"), evts)
	evts = filter(row -> !ismissing(row.word_key) && row.word_key != "**", evts)

	println(size(dat))
	println(size(evts))
	println("Valid word fixations: ", nrow(evts))
first(evts[:, [:word, :word_key, :duration_ms, :latency]], 10)
end	
	

# ╔═╡ f44e7a04-2b38-43db-9025-6a96343e1774
begin
	lexical_from_csv = CSV.read("/scratch/data/ROAMM/outputs/roamm_words_with_surprisal.csv", DataFrame)

	lexical = select(lexical_from_csv, :word_key, :word_length, :frequency_zipf, :surprisal_bits)

	println("Lexical rows: ", nrow(lexical))
	println("Unique word keys: ", length(unique(lexical.word_key)))

end

# ╔═╡ 7442f9a1-063d-4c95-9533-7b8d4050ab50
begin
	evts_b2b = leftjoin(
    evts,
    lexical,
    on = :word_key
)

println("Before join: ", nrow(evts))
println("After join:  ", nrow(evts_b2b))
end

# ╔═╡ 4e7aff75-362a-495c-89bb-0366036d1985
begin
	sort!(evts_b2b, :latency)
	
	first(
	select(
		evts_b2b,
		:word,
		:latency,
		:word_length,
		:frequency_zipf,
		:surprisal_bits,
	),
	100,
)
end

# ╔═╡ a185f004-16db-4fa2-9565-8b53a80e301d
for col in [:word_length, :frequency_zipf, :surprisal_bits]
    println(col, ": ", count(ismissing, evts_b2b[!, col]))
end

# ╔═╡ 7cf428e0-26f5-43d8-b7a5-896d324dd02c
begin
	b2b_events_file = joinpath(
    DATA_DIR,
    "sub-10014_run1_b2b_events.csv"
)

CSV.write(b2b_events_file, evts_b2b)
end

# ╔═╡ 887b69cc-4e18-4881-a28c-f5d9aa5a0776
begin
	dat_uv = dat .* 1e6
	extrema(dat_uv)
end

# ╔═╡ 4c9a8478-9f60-4504-b986-ec069e7df680
size(dat)

# ╔═╡ 01477af5-721c-4640-9a9a-c7c6cbd660d6
series(dat_uv[:,1:100],solid_color=:black)

# ╔═╡ e7e334f7-f5da-4476-9544-e8d2fa9d354d
let
	ix = 5000:10000
lines(ix,dat_uv[7,ix])
	vlines!(evts.latency[30:45])
current_figure()
end

# ╔═╡ 5820e6e4-f88a-419c-9506-e60b8b417f88
hist(evts.duration_ms,bins=100)

# ╔═╡ e96f8b0c-6ce0-4e06-8ce4-00d1531a2dfd
evts_b2b.frequency_zipf_c = evts_b2b.frequency_zipf .- mean(evts_b2b.frequency_zipf) 

# ╔═╡ 4c8863d8-4c8d-468c-964a-8735b9c3a3d0
evts_b2b.word_length_c = evts_b2b.word_length .- mean(evts_b2b.word_length) 

# ╔═╡ ccc13664-e3c1-4069-b663-d1003af3d36d
m = fit(UnfoldModel,@formula(0~1+word_length_c+frequency_zipf_c),evts_b2b,dat_uv,firbasis(sfreq=256,τ=(-0.2,0.8)))

# ╔═╡ 82be34e0-77bf-4c6d-ba4b-7b9cd6931485
begin
	dat_e,times_e = Unfold.epoch(dat_uv,evts_b2b,[-0.2,0.8],256)
	m_e = fit(UnfoldModel,@formula(0~1+word_length_c+frequency_zipf_c),evts_b2b,dat_e,times_e)
end

# ╔═╡ 46a50f0f-a02e-4f23-83f4-a0a645d51b08
heatmap(reshape(coef(m_e)[:,:,:],64,:)')

# ╔═╡ da678d98-4642-4723-8823-72f785426ab1
heatmap(coef(m)')

# ╔═╡ 7a7b404b-b821-4864-8afb-96409d2f213c
series(reshape(coef(m_e),64,:)[:,1:250],solid_color=:black)

# ╔═╡ 2974c740-746a-4c9b-a57d-b534b06428de
series(coef(m)[:,1:250],solid_color=:black)

# ╔═╡ 9c0c2fff-de08-4045-8f65-f2d4bbdf122f


# ╔═╡ 19f51df0-7536-4427-bb5e-93bd77cb7ef6
# ╠═╡ disabled = true
#=╠═╡
two_step_roamm = MScB2B.fit_two_step_b2b_roamm(
	dat_uv,
	evts_b2b,
	sfreq = 256.0,
	cross_val_reps = 1
)
  ╠═╡ =#

# ╔═╡ 444278c7-0000-4a0e-bb47-6fc83a0dd4a3
#=╠═╡
series(coef(two_step_roamm[4])[1,:,:]')
  ╠═╡ =#

# ╔═╡ f690a840-7160-44dc-91fe-75a4e57b8e8a
begin
scores_roamm_uv = DataFrame(
    Unfold.coeftable(
        two_step_roamm_uv.b2b_fit_corrected
    )
)

scores_roamm_uv[!, :estimate_signed] =
    copy(scores_roamm_uv.estimate)

scores_roamm_uv.estimate .= abs.(scores_roamm_uv.estimate)

scores_roamm_uv = subset(
    scores_roamm_uv,
    :coefname => x -> x .== "surprisal_z"
)
end

# ╔═╡ e8261bd2-cdd5-46cb-8d3e-856ac9bfcd4e
#=╠═╡
(
	corrected_trials = size(two_step_roamm.corrected_trials),
	n_events = nrow(two_step_roamm.events),
	rerp_model = typeof(two_step_roamm.rerp_model),
	b2b_model = typeof(two_step_roamm.b2b_fit_corrected)
)
  ╠═╡ =#

# ╔═╡ 7098ade1-258a-4ab8-a668-45ac3e7ead1b
#=╠═╡
result_raw = DataFrame(
	Unfold.coeftable(
		two_step_roamm.b2b_fit_corrected
	)
)
  ╠═╡ =#

# ╔═╡ 93b90f33-0d6f-4dff-bfe9-e2e894a606bc
#=╠═╡
begin
    result = copy(result_raw)

    result[!, :estimate_signed] =
        copy(result.estimate)

    result.estimate .= abs.(result.estimate)

    filter!(
        :coefname => !=("(Intercept)"),
        result
    )

    result
end
  ╠═╡ =#

# ╔═╡ f04ff098-800c-4cc4-b798-5fd6061595a5
#=╠═╡
begin
    scores_roamm = DataFrame(
        Unfold.coeftable(
            two_step_roamm.b2b_fit_corrected
        )
    )

    # keep signed estimate for debugging
    scores_roamm[!, :estimate_signed] = copy(scores_roamm.estimate)

    # B2B: use magnitude
    scores_roamm.estimate .= abs.(scores_roamm.estimate)

    # remove intercept
    scores_roamm = subset(
        scores_roamm,
        :coefname => x -> x .!= "(Intercept)"
    )
end
  ╠═╡ =#

# ╔═╡ e4903dc3-8749-46da-a925-3901b1e2e288
#=╠═╡
# uV+Ridge

MScB2B.plot_erp(
    subset(
        scores_roamm,
        :coefname => x -> x .== "surprisal_z"
    )
)
  ╠═╡ =#

# ╔═╡ 535f187b-50ce-4626-8a80-46bfcf981f9b


# ╔═╡ 9c4e91e7-ffac-4045-a7f6-3ecd6ecd18a7
#=╠═╡
size(two_step_roamm.corrected_trials)
  ╠═╡ =#

# ╔═╡ cb245384-d60f-469a-8bec-b95ab6581172
#=╠═╡
extrema(skipmissing(vec(two_step_roamm.corrected_trials)))
  ╠═╡ =#

# ╔═╡ 6f0b62fe-f34b-4059-9d73-e668f368276a
#=╠═╡
extrema(scores_roamm.estimate_signed)
  ╠═╡ =#

# ╔═╡ 76f53762-2cf7-47ee-a42d-b1658cca6abe
#=╠═╡
sort(scores_roamm, :estimate, rev=true)[1:10, [:time, :estimate_signed, :estimate]]
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═3ee335c0-abc6-11f1-90c5-05876e25c543
# ╠═7d609d2d-6e9b-4350-8cdb-52b3a4ded147
# ╠═6886a972-e132-46e4-a2f9-530ce4dcdbeb
# ╠═c64c8d29-d53c-4846-adc2-64313139e761
# ╠═3d7235bf-62ec-4cd0-baa4-32c1220f5983
# ╠═f44e7a04-2b38-43db-9025-6a96343e1774
# ╠═7442f9a1-063d-4c95-9533-7b8d4050ab50
# ╠═4e7aff75-362a-495c-89bb-0366036d1985
# ╠═a185f004-16db-4fa2-9565-8b53a80e301d
# ╠═7cf428e0-26f5-43d8-b7a5-896d324dd02c
# ╠═887b69cc-4e18-4881-a28c-f5d9aa5a0776
# ╠═4c9a8478-9f60-4504-b986-ec069e7df680
# ╠═01477af5-721c-4640-9a9a-c7c6cbd660d6
# ╠═e7e334f7-f5da-4476-9544-e8d2fa9d354d
# ╠═5820e6e4-f88a-419c-9506-e60b8b417f88
# ╠═e96f8b0c-6ce0-4e06-8ce4-00d1531a2dfd
# ╠═4c8863d8-4c8d-468c-964a-8735b9c3a3d0
# ╠═ccc13664-e3c1-4069-b663-d1003af3d36d
# ╠═82be34e0-77bf-4c6d-ba4b-7b9cd6931485
# ╠═46a50f0f-a02e-4f23-83f4-a0a645d51b08
# ╠═da678d98-4642-4723-8823-72f785426ab1
# ╠═7a7b404b-b821-4864-8afb-96409d2f213c
# ╠═2974c740-746a-4c9b-a57d-b534b06428de
# ╠═9c0c2fff-de08-4045-8f65-f2d4bbdf122f
# ╠═19f51df0-7536-4427-bb5e-93bd77cb7ef6
# ╠═444278c7-0000-4a0e-bb47-6fc83a0dd4a3
# ╠═f690a840-7160-44dc-91fe-75a4e57b8e8a
# ╠═e8261bd2-cdd5-46cb-8d3e-856ac9bfcd4e
# ╠═7098ade1-258a-4ab8-a668-45ac3e7ead1b
# ╠═93b90f33-0d6f-4dff-bfe9-e2e894a606bc
# ╠═f04ff098-800c-4cc4-b798-5fd6061595a5
# ╠═e4903dc3-8749-46da-a925-3901b1e2e288
# ╠═535f187b-50ce-4626-8a80-46bfcf981f9b
# ╠═9c4e91e7-ffac-4045-a7f6-3ecd6ecd18a7
# ╠═cb245384-d60f-469a-8bec-b95ab6581172
# ╠═6f0b62fe-f34b-4059-9d73-e668f368276a
# ╠═76f53762-2cf7-47ee-a42d-b1658cca6abe
