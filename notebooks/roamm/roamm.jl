### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 218de97b-c7cb-4667-b993-94031552dd63
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 71ac26a6-c9b0-4933-b9e7-b0a4fe37ad48
begin
	using CSV
	using DataFrames
	using NPZ
end

# ╔═╡ de9502bb-2a16-4517-966d-92d249981f6c
using Statistics

# ╔═╡ b9c0f854-9ebf-46ac-a4e4-6113d71c0161
const DATA_DIR = "/scratch/data/ROAMM/outputs/"

# ╔═╡ 8fd55001-633c-4abe-8fd4-29d6c64df17e
begin
	eeg_file = joinpath(DATA_DIR, "sub-10014_run1_eeg.npy")
	events_files = joinpath(DATA_DIR, "roamm_fixation_events.csv")
	dat = npzread(eeg_file)
	evts_all = CSV.read(events_files, DataFrame)

	println("EEG size: ", size(dat))
	println("All events: ", size(evts_all))
end

# ╔═╡ bd79afd9-56e9-48a7-a549-69b45b23207c
begin
	evts = filter(
		row -> row.subject_id == "sub-10014" && row.run_num == 1,
		evts_all
	)

	println("Matched events: ", nrow(evts))
	println(unique(evts_all.subject_id))
	println(unique(evts_all.run_num))
end

# ╔═╡ c502dafb-d7ce-477f-b23b-f13f06d04ff4
begin
	println("EEG samples: ", size(dat, 2))
println("Latency min: ", minimum(evts.latency))
println("Latency max: ", maximum(evts.latency))

@assert minimum(evts.latency) >= 1
@assert maximum(evts.latency) <= size(dat, 2)
end

# ╔═╡ 8248122d-57ee-4d86-8397-b60b7cd60a13
begin
println("EEG size: ", size(dat))
println("Matched events: ", nrow(evts))

println(unique(evts.subject_id))
println(unique(evts.run_num))
println(names(evts))
end

# ╔═╡ b4b7df86-2061-4b90-8cd2-f2eecf566058
first(
    select(
        evts,
        :onset_s,
        :eeg_sample,
        :latency,
        :eeg_time_s,
        :sync_error_ms,
    ),
    10,
)

# ╔═╡ 276c00ea-9489-4249-b860-e12bf1fddaa1
all(evts.latency .== evts.eeg_sample .+ 1)

# ╔═╡ dd0bfb1d-48a6-4c76-b83f-b5bb10537ad2
begin
	sfreq = 256.0

eeg_time_check = (evts.latency .- 1) ./ sfreq

maximum(abs.(eeg_time_check .- evts.eeg_time_s))

errors = abs.(evts.sync_error_ms)

println("Mean sync error: ", mean(errors), " ms")
println("Max sync error:  ", maximum(errors), " ms")
end


# ╔═╡ 03a4464b-70f0-4549-9fe3-014aca58217d
begin
	filtered_events_file = joinpath(
    DATA_DIR,
    "sub-10014_run1_fixation_events.csv"
)

CSV.write(filtered_events_file, evts)

println("Saved to: ", filtered_events_file)
end

# ╔═╡ 8018a40c-9156-42c1-8809-6a68285cff21
begin
	evts_check = CSV.read(filtered_events_file, DataFrame)
	
	println(size(evts_check))
	println(unique(evts_check.subject_id))
	println(unique(evts_check.run_num))
end

# ╔═╡ Cell order:
# ╠═218de97b-c7cb-4667-b993-94031552dd63
# ╠═71ac26a6-c9b0-4933-b9e7-b0a4fe37ad48
# ╠═b9c0f854-9ebf-46ac-a4e4-6113d71c0161
# ╠═8fd55001-633c-4abe-8fd4-29d6c64df17e
# ╠═bd79afd9-56e9-48a7-a549-69b45b23207c
# ╠═c502dafb-d7ce-477f-b23b-f13f06d04ff4
# ╠═8248122d-57ee-4d86-8397-b60b7cd60a13
# ╠═b4b7df86-2061-4b90-8cd2-f2eecf566058
# ╠═276c00ea-9489-4249-b860-e12bf1fddaa1
# ╠═de9502bb-2a16-4517-966d-92d249981f6c
# ╠═dd0bfb1d-48a6-4c76-b83f-b5bb10537ad2
# ╠═03a4464b-70f0-4549-9fe3-014aca58217d
# ╠═8018a40c-9156-42c1-8809-6a68285cff21
