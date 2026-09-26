### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 032a0c56-b0dd-11f1-8b66-7d0d34715e27
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 7d21e4fc-0b3a-4725-b91a-ea41ff7b55ef
using PlutoLinks: @revise, @ingredients

# ╔═╡ cd005e17-5b37-45fb-9236-312a6c3a8ad7
@revise using MScB2B

# ╔═╡ 5f625afd-c55a-4ae4-939d-a8daa2ab355c
begin
	using Serialization
	using DataFrames
	using Statistics
	using NPZ
end

# ╔═╡ b1a2cf70-5d39-4c2a-b04f-4fd3957edcbf
pwd()

# ╔═╡ 8b54dd6e-4c34-4201-aade-5b4f53d83a8a
project_root = normpath(joinpath(@__DIR__, ".."))

# ╔═╡ c690b5d9-d770-45d9-8077-d937f1404134
begin 
	result = deserialize(
		joinpath(
			project_root,
			"results",
			"roamm_two_step_lsq",
			"sub-10014_two_step_lsq.jls"
		)
	)
	scores = result.scores
end

# ╔═╡ 30969bec-ead1-4d3a-9e8d-786a4f30101c
scores_surprisal = subset(
    scores,
    :coefname => x -> x .== "surprisal_z"
)

# ╔═╡ 3beb63f4-28f6-4070-80e9-c435804e8908
MScB2B.plot_erp(scores_surprisal)

# ╔═╡ 03a93966-6ccf-4e1a-81ea-62104b3c4875
dat = npzread("/scratch/data/ROAMM/outputs/eeg/sub-10014_run1_eeg.npy")

# ╔═╡ 86d4ab9a-367a-49cd-9634-1f94426d782f
size(dat)

# ╔═╡ 8374d32a-268b-475c-a4fc-5fa70e3db801
dat[1, 1:100]


# ╔═╡ a8afa9db-afe3-4f32-a1c2-e2c4d4c127be
channel_names = [
    "Fp1","AF7","AF3","F1","F3","F5","F7","FT7",
    "FC5","FC3","FC1","C1","C3","C5","T7","TP7",
    "CP5","CP3","CP1","P1","P3","P5","P7","P9",
    "PO7","PO3","O1","Iz","Oz","POz","Pz","CPz",
    "Fpz","Fp2","AF8","AF4","Afz","Fz","F2","F4",
    "F6","F8","FT8","FC6","FC4","FC2","FCz","Cz",
    "C2","C4","C6","T8","TP8","CP6","CP4","CP2",
    "P2","P4","P6","P8","P10","PO8","PO4","O2"
]

# ╔═╡ 104a5d38-394a-4e85-9ffe-cbfa45e557a5
begin
	channel_stats = DataFrame(
    channel = channel_names,

    mean_μV = [
        mean(dat[ch, :])
        for ch in 1:size(dat, 1)
    ],

    std_μV = [
        std(dat[ch, :])
        for ch in 1:size(dat, 1)
    ],

    min_μV = [
        minimum(dat[ch, :])
        for ch in 1:size(dat, 1)
    ],

    max_μV = [
        maximum(dat[ch, :])
        for ch in 1:size(dat, 1)
    ],

    p01_μV = [
        quantile(dat[ch, :], 0.01)
        for ch in 1:size(dat, 1)
    ],

    p99_μV = [
        quantile(dat[ch, :], 0.99)
        for ch in 1:size(dat, 1)
    ],
)

channel_stats
end

# ╔═╡ Cell order:
# ╠═032a0c56-b0dd-11f1-8b66-7d0d34715e27
# ╠═7d21e4fc-0b3a-4725-b91a-ea41ff7b55ef
# ╠═cd005e17-5b37-45fb-9236-312a6c3a8ad7
# ╠═5f625afd-c55a-4ae4-939d-a8daa2ab355c
# ╠═b1a2cf70-5d39-4c2a-b04f-4fd3957edcbf
# ╠═8b54dd6e-4c34-4201-aade-5b4f53d83a8a
# ╠═c690b5d9-d770-45d9-8077-d937f1404134
# ╠═30969bec-ead1-4d3a-9e8d-786a4f30101c
# ╠═3beb63f4-28f6-4070-80e9-c435804e8908
# ╠═03a93966-6ccf-4e1a-81ea-62104b3c4875
# ╠═86d4ab9a-367a-49cd-9634-1f94426d782f
# ╠═8374d32a-268b-475c-a4fc-5fa70e3db801
# ╠═a8afa9db-afe3-4f32-a1c2-e2c4d4c127be
# ╠═104a5d38-394a-4e85-9ffe-cbfa45e557a5
