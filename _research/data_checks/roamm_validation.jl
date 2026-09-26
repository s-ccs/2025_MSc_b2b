### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ bef71bca-a561-11f1-b7a9-3f7785b36ca5
begin 
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ aaa8a0f4-ec53-4353-950e-9f1b0703ce41
begin
	using PythonCall
	pd = pyimport("pandas")
end

# ╔═╡ bbd37711-e086-45d4-adb8-0e933c8f094d
using Statistics

# ╔═╡ 199165ad-d364-4303-a7dc-6a4154064b4c
begin
	using Unfold
	using StatsModels
	
	basis = firbasis(τ = (-0.2, 0.5), sfreq = 256)
	
	design = Dict(
	    "fixation" => (@formula(0 ~ 1), basis)
	)
end

# ╔═╡ 544061b2-1f01-4a14-a8f3-925aeb30bf34
begin
	path = "/scratch/data/ROAMM/derivatives/synced/sub-10014/sub-10014_task-ReMind_run-01_mldata.pkl"

	x = pd.read_pickle(path)
end

# ╔═╡ 4974bced-a007-4847-8991-f82571d3589e
x.shape

# ╔═╡ d955f259-bd5a-4a16-b513-c56bbae79b66
x.columns

# ╔═╡ ef7d587a-31d2-4fb4-96a8-f8d49dacae56
begin
	fix_cols = pybuiltins.list([
	    "time",
	    "is_fix",
	    "fix_L_tStart",
	    "fix_L_tEnd",
	    "fix_L_duration",
	    "fix_L_xAvg",
	    "fix_L_yAvg",
	    "fix_L_fixed_word",
	    "fix_L_fixed_word_key"
	])
	
	fix_samples = x.filter(items=fix_cols)
end

# ╔═╡ e1d661c3-bd90-4172-91e9-5c768b01fb4e
fix_samples.dropna(
    subset=pybuiltins.list(["fix_L_fixed_word"])
).head(20)

# ╔═╡ 5f3a04dc-a598-4989-80ff-581dd2ad3526
md"""
# Create one-row-per-fixation table
"""

# ╔═╡ 68c3714e-65fb-4b29-8e2a-3adf038e91e6
fixevents = fix_samples.dropna(
    subset=pybuiltins.list(["fix_L_tStart"])
).drop_duplicates(
    subset=pybuiltins.list(["fix_L_tStart"])
)

# ╔═╡ f48dc289-65d8-400b-a38d-a678a3e9f4db
fixevents.head(20)

# ╔═╡ 978505fd-6ab1-4ec0-8d1d-87f08283b684
fixevents.shape

# ╔═╡ 966d8389-fa0d-44c9-b631-acf61b19c52c


# ╔═╡ bb24cf7f-7440-4ae5-91c9-5cc74183e913


# ╔═╡ 4c6b31d9-c94e-41b8-a9e3-cd379a3c1bfb
begin
	eeg_time = pyconvert(Vector{Float64}, x["time"].to_numpy())

	dt = diff(eeg_time)

(
	first_time = first(eeg_time),
	last_time = last(eeg_time),
	median_dt = median(dt),
	sfreq = 1 / median(dt),
	n_samples = length(eeg_time)
)
end

# ╔═╡ f437c396-c726-463b-a2cb-fc1cbc78ac54
begin
	fix_onsets = pyconvert(
    Vector{Float64},
    fixevents["fix_L_tStart"].to_numpy()
)

latencies = [
    argmin(abs.(eeg_time .- onset))
    for onset in fix_onsets
]
end

# ╔═╡ 2c215824-e567-4979-b005-07868dfbfb3f
for i in 1:10
    println(
        "fixation onset = ", fix_onsets[i],
        " s   → EEG sample = ", latencies[i],
        "   → EEG time = ", eeg_time[latencies[i]]
    )
end

# ╔═╡ efc179d2-11bf-4b60-96a2-254bd09bff60
begin
    fix_duration = pyconvert(
        Vector{Float64},
        fixevents["fix_L_duration"].to_numpy()
    )

    fix_x = pyconvert(
        Vector{Float64},
        fixevents["fix_L_xAvg"].to_numpy()
    )

    fix_y = pyconvert(
        Vector{Float64},
        fixevents["fix_L_yAvg"].to_numpy()
    )
end

# ╔═╡ 75038155-a5d6-4ad2-ba4b-d90468bad2dc
begin
	using DataFrames
	
	events = DataFrame(
	    latency = latencies,
	    onset = fix_onsets,
	    duration = fix_duration,
	    x = fix_x,
	    y = fix_y,
	)
end

# ╔═╡ fa0058ca-4a1a-4f24-8e8a-d7465fa9c542
latencies_fast = round.(Int, fix_onsets .* 256) .+ 1

# ╔═╡ 68dff82a-045b-421d-9e5e-471cf5294651
all(latencies_fast .== latencies)

# ╔═╡ 2b834e6a-e7bb-45e9-95cd-e38b7bbcef71
cols = pyconvert(Vector{String}, x.columns.tolist())

# ╔═╡ 3273cc39-5b6c-468f-b723-5ba0b2694bd6
filter(cols) do c
    occursin("fix", lowercase(c))
end

# ╔═╡ 5818d59e-a8a7-4a6d-abcc-c999f37c1f21
begin
    L_start = pyconvert(Vector{Float64},
        x["fix_L_tStart"].dropna().drop_duplicates().to_numpy())

    R_start = pyconvert(Vector{Float64},
        x["fix_R_tStart"].dropna().drop_duplicates().to_numpy())

    (
        n_left_fixations = length(L_start),
        n_right_fixations = length(R_start),
    )
end

# ╔═╡ 6d8f6da5-24d3-4242-b1a3-3a4ee4b236d5
for i in 1:min(10, length(L_start), length(R_start))
    println(
        "L = ", L_start[i],
        "   R = ", R_start[i],
        "   difference = ", abs(L_start[i] - R_start[i]) * 1000, " ms"
    )
end

# ╔═╡ 7475bb6a-3bb3-4349-a1d0-03f36e7ca575
md"""
# For initial validation, use left-eye fixations only.
# Final analysis should define a consistent binocular/reference-eye rule.
"""

# ╔═╡ 824704d4-4a11-4e73-a96d-05598198844a
# EEG matrix
begin
    time_idx = findfirst(==("time"), cols)
    eeg_cols = cols[1:time_idx-1]

    length(eeg_cols), eeg_cols
end

# ╔═╡ d982c367-9596-4913-875f-db5c6214b164
begin
    eeg_py = x.filter(items=pybuiltins.list(eeg_cols)).to_numpy()

    # pandas: samples × channels
    eeg_samples_channels = pyconvert(Matrix{Float64}, eeg_py)

    # Unfold/EEG analysis usually更方便用 channels × samples
    eeg = permutedims(eeg_samples_channels)

    size(eeg)
end

# ╔═╡ c7e7f428-3fb6-485b-a415-1831e8fd99be
# Fixation-locked average
begin
	oz = findfirst(==("Oz"), eeg_cols)
	tmin = -0.2
	tmax = 0.5
	
	pre  = round(Int, abs(tmin) * 256)
	post = round(Int, tmax * 256)

	valid_latencies = filter(events.latency) do lat
    lat - pre >= 1 && lat + post <= size(eeg, 2)
	end
	epochs = [
    begin
        ep = copy(eeg[oz, lat-pre:lat+post])
        ep .-= mean(ep[1:pre])
        ep
    end
    for lat in valid_latencies
]

erp = mean(reduce(hcat, epochs), dims=2)[:, 1]

	times = collect(-pre:post) ./ 256
end

# ╔═╡ fb465b04-d3fa-43c5-91cf-9bfb6b2077ab
begin
	using CairoMakie

fig = Figure()
ax = Axis(
    fig[1, 1],
    xlabel = "Time from fixation onset (s)",
    ylabel = "Amplitude",
    title = "Fixation-locked average at Oz"
)

lines!(ax, times, erp)
vlines!(ax, [0])

fig
end

# Very good sanity check. There is a clear posterior fixation-locked response around 100 ms.

# ╔═╡ f4acd753-c8ac-4eb7-8f7d-9b0899b11ea9
(
    n_events = nrow(events),
    first_latency = minimum(events.latency),
    last_latency = maximum(events.latency),
    eeg_samples = size(eeg, 2),
    all_inside = all(1 .<= events.latency .<= size(eeg, 2)),
)

# ╔═╡ e8833aee-f15c-4c01-bf38-5a08364ae8bb
begin
	word_fixevents = fixevents.dropna(
    subset = pybuiltins.list(["fix_L_fixed_word"])
)

word_fixevents.shape
end

# ╔═╡ 1216151b-8fbf-4d07-bd36-f4f992ac5c33
begin
    word_onsets = pyconvert(
        Vector{Float64},
        word_fixevents["fix_L_tStart"].to_numpy()
    )

    word_durations = pyconvert(
        Vector{Float64},
        word_fixevents["fix_L_duration"].to_numpy()
    )

    word_x = pyconvert(
        Vector{Float64},
        word_fixevents["fix_L_xAvg"].to_numpy()
    )

    word_y = pyconvert(
        Vector{Float64},
        word_fixevents["fix_L_yAvg"].to_numpy()
    )

    words = pyconvert(
        Vector{String},
        word_fixevents["fix_L_fixed_word"].tolist()
    )

    word_latencies = round.(Int, word_onsets .* 256) .+ 1
end

# ╔═╡ 20c12800-b2ff-40a0-9323-65ecfd02a8bc
word_events = DataFrame(
    latency = word_latencies,
    onset = word_onsets,
    duration = word_durations,
    word = words,
    x = word_x,
    y = word_y,
)



# ╔═╡ 364b1427-7e39-49c8-8e60-b26bb9294039
begin
	valid_word_latencies = filter(word_events.latency) do lat
    lat - pre >= 1 && lat + post <= size(eeg, 2)
end

word_epochs = [
    begin
        ep = copy(eeg[oz, lat-pre:lat+post])
        ep .-= mean(ep[1:pre])
        ep
    end
    for lat in valid_word_latencies
]

word_erp = mean(reduce(hcat, word_epochs), dims=2)[:, 1]
	end

# ╔═╡ 568a0cb7-ec56-48d6-a220-d26c80fef6d5
begin
	fig_ = Figure()

ax_ = Axis(
    fig_[1, 1],
    xlabel = "Time from fixation onset (s)",
    ylabel = "Amplitude (V)",
    title = "Word-fixation locked average at Oz"
)

lines!(ax_, times, word_erp)
vlines!(ax_, [0])

fig_
	
end

# ╔═╡ 88e6f331-d011-437e-83c7-f6be1e08d8e6
word_events.event = fill("fixation", nrow(word_events))

# ╔═╡ Cell order:
# ╠═bef71bca-a561-11f1-b7a9-3f7785b36ca5
# ╠═aaa8a0f4-ec53-4353-950e-9f1b0703ce41
# ╠═544061b2-1f01-4a14-a8f3-925aeb30bf34
# ╠═4974bced-a007-4847-8991-f82571d3589e
# ╠═d955f259-bd5a-4a16-b513-c56bbae79b66
# ╠═ef7d587a-31d2-4fb4-96a8-f8d49dacae56
# ╠═e1d661c3-bd90-4172-91e9-5c768b01fb4e
# ╟─5f3a04dc-a598-4989-80ff-581dd2ad3526
# ╠═68c3714e-65fb-4b29-8e2a-3adf038e91e6
# ╠═f48dc289-65d8-400b-a38d-a678a3e9f4db
# ╠═978505fd-6ab1-4ec0-8d1d-87f08283b684
# ╠═966d8389-fa0d-44c9-b631-acf61b19c52c
# ╠═bb24cf7f-7440-4ae5-91c9-5cc74183e913
# ╠═bbd37711-e086-45d4-adb8-0e933c8f094d
# ╠═4c6b31d9-c94e-41b8-a9e3-cd379a3c1bfb
# ╠═f437c396-c726-463b-a2cb-fc1cbc78ac54
# ╠═2c215824-e567-4979-b005-07868dfbfb3f
# ╠═efc179d2-11bf-4b60-96a2-254bd09bff60
# ╠═75038155-a5d6-4ad2-ba4b-d90468bad2dc
# ╠═fa0058ca-4a1a-4f24-8e8a-d7465fa9c542
# ╠═68dff82a-045b-421d-9e5e-471cf5294651
# ╠═2b834e6a-e7bb-45e9-95cd-e38b7bbcef71
# ╠═3273cc39-5b6c-468f-b723-5ba0b2694bd6
# ╠═5818d59e-a8a7-4a6d-abcc-c999f37c1f21
# ╠═6d8f6da5-24d3-4242-b1a3-3a4ee4b236d5
# ╟─7475bb6a-3bb3-4349-a1d0-03f36e7ca575
# ╠═824704d4-4a11-4e73-a96d-05598198844a
# ╠═d982c367-9596-4913-875f-db5c6214b164
# ╠═f4acd753-c8ac-4eb7-8f7d-9b0899b11ea9
# ╠═c7e7f428-3fb6-485b-a415-1831e8fd99be
# ╠═fb465b04-d3fa-43c5-91cf-9bfb6b2077ab
# ╠═e8833aee-f15c-4c01-bf38-5a08364ae8bb
# ╠═1216151b-8fbf-4d07-bd36-f4f992ac5c33
# ╠═20c12800-b2ff-40a0-9323-65ecfd02a8bc
# ╠═364b1427-7e39-49c8-8e60-b26bb9294039
# ╠═568a0cb7-ec56-48d6-a220-d26c80fef6d5
# ╠═88e6f331-d011-437e-83c7-f6be1e08d8e6
# ╠═199165ad-d364-4303-a7dc-6a4154064b4c
