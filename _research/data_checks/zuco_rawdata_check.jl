### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 5c0e85d2-a31f-11f1-aca2-4d7ba3bde0bd
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
end

# ╔═╡ 2254a607-5c0b-4098-b526-6d110a5ba284
using MAT

# ╔═╡ 63990c20-f487-4a23-a770-12ec5247c79e
path = "/scratch/data/ZuCo2/raw/YRP_NR1_EEG.mat"

# ╔═╡ adcbe945-2b03-4f4c-9ba7-7d49ddd8197f
f = matopen(path)

# ╔═╡ c9034718-207a-4034-a62a-03bad82ad421
keys(f)

# ╔═╡ d7de0bb5-fd93-41f2-8683-16192c0cc00e
EEG = read(f, "EEG")

# ╔═╡ a0c2220e-a2d9-4b8c-bf93-76995f780d79
close(f)

# ╔═╡ 7ef24a3d-44f8-4ccc-aa33-6ccde98d2404
typeof(EEG)

# ╔═╡ 837938d2-850a-4923-9b84-85213c82addc
keys(EEG)

# ╔═╡ 3adf3b51-b0b9-4069-b565-fa4f7772318b
typeof(EEG["data"])

# ╔═╡ 11d55963-9ee2-44ae-9bf7-a4a82606bc5d
size(EEG["data"])  # continous EEG

# ╔═╡ da844bb6-4ba1-420c-93da-c83d36dd2b31
(
	EEG["srate"],  # 500 hz
	EEG["nbchan"], 
	EEG["pnts"],   # 250062 samples
	EEG["trials"]
)

# ╔═╡ eed56ed2-ff8a-4af0-ac40-a74ef82000ad
ev = EEG["event"]

# ╔═╡ f5154cfe-cfb3-4ebc-9196-332b7c9cf7fe
typeof(ev)

# ╔═╡ 46d6e631-0a55-4c7a-81a7-f0d169387ce1
vec(ev["type"])[1:10]

# ╔═╡ dde5699d-0905-46c8-9515-86d10cc673dc
vec(ev["latency"])[1:10]

# ╔═╡ 27bd0663-5f2b-4b07-adaa-df57f680a75d


# ╔═╡ 1e4e5534-0c1d-4a44-9655-0a95ae5ce9d6


# ╔═╡ 35c9474a-a676-4861-a693-4ab5a7d3da62


# ╔═╡ 07cd5150-d83b-494e-9f4e-b2c53e3a8aff


# ╔═╡ 6ef8c338-6ff5-45a5-a416-be5c668b9ecd


# ╔═╡ b4d72282-5fd7-42c8-b92c-19ef1e768b51


# ╔═╡ 82a8b74b-2222-4a6d-bab1-410f7319e7e8
path_ET = "/scratch/data/ZuCo2/raw/YRP_NR1_ET.mat"

# ╔═╡ 1de3f0ec-925a-4a8e-abfd-641c2809b740
f_ET = matopen(path_ET)

# ╔═╡ 434664c0-c36f-4e35-810c-466cafe58def
keys(f_ET)

# ╔═╡ f2a4ef2e-f82b-43d7-a825-964d4a312b81
ET = matread(path_ET)

# ╔═╡ a2bb5f43-c402-4301-8885-e1fea58412fa
typeof(ET)

# ╔═╡ 446d43b8-77bf-410a-b3e0-0cf250d329d4
keys(ET)

# ╔═╡ 82b9e895-db48-4503-b66d-78ec178b0420
ET["event"]

# ╔═╡ 0831c71b-4979-43d3-a77a-966f4586c3b6
typeof(ET["event"])

# ╔═╡ 72aa35ff-af19-4ee0-b857-d582f4df70e8
size(ET["event"])

# ╔═╡ 18e0ee0f-882f-4aed-8e86-4906942037d8
ET["event"][1:10, :]

# ╔═╡ c56d43bb-3f9b-4e71-8789-44e2ea1a9282


# ╔═╡ 1c60bcff-5a60-41c4-bd12-3580a8e87680
begin
	et_time = ET["event"][:, 1]
	et_type = Int.(round.(ET["event"][:, 2]))

	eeg_type = strip.(string.(vec(ev["type"])))
	eeg_latency = Float64.(vec(ev["latency"]))
end

# ╔═╡ e08bc378-f23e-4ac4-b680-aa7133ab8a28
[
    (
        i = i,
        et_time = et_time[i],
        et_type = et_type[i],
        eeg_type = eeg_type[i],
        eeg_sample = eeg_latency[i]
    )
    for i in 1:10
]

# ╔═╡ Cell order:
# ╠═5c0e85d2-a31f-11f1-aca2-4d7ba3bde0bd
# ╠═2254a607-5c0b-4098-b526-6d110a5ba284
# ╠═63990c20-f487-4a23-a770-12ec5247c79e
# ╠═adcbe945-2b03-4f4c-9ba7-7d49ddd8197f
# ╠═c9034718-207a-4034-a62a-03bad82ad421
# ╠═d7de0bb5-fd93-41f2-8683-16192c0cc00e
# ╠═a0c2220e-a2d9-4b8c-bf93-76995f780d79
# ╠═7ef24a3d-44f8-4ccc-aa33-6ccde98d2404
# ╠═837938d2-850a-4923-9b84-85213c82addc
# ╠═3adf3b51-b0b9-4069-b565-fa4f7772318b
# ╠═11d55963-9ee2-44ae-9bf7-a4a82606bc5d
# ╠═da844bb6-4ba1-420c-93da-c83d36dd2b31
# ╠═eed56ed2-ff8a-4af0-ac40-a74ef82000ad
# ╠═f5154cfe-cfb3-4ebc-9196-332b7c9cf7fe
# ╠═46d6e631-0a55-4c7a-81a7-f0d169387ce1
# ╠═dde5699d-0905-46c8-9515-86d10cc673dc
# ╠═27bd0663-5f2b-4b07-adaa-df57f680a75d
# ╠═1e4e5534-0c1d-4a44-9655-0a95ae5ce9d6
# ╠═35c9474a-a676-4861-a693-4ab5a7d3da62
# ╠═07cd5150-d83b-494e-9f4e-b2c53e3a8aff
# ╠═6ef8c338-6ff5-45a5-a416-be5c668b9ecd
# ╠═b4d72282-5fd7-42c8-b92c-19ef1e768b51
# ╠═82a8b74b-2222-4a6d-bab1-410f7319e7e8
# ╠═1de3f0ec-925a-4a8e-abfd-641c2809b740
# ╠═434664c0-c36f-4e35-810c-466cafe58def
# ╠═f2a4ef2e-f82b-43d7-a825-964d4a312b81
# ╠═a2bb5f43-c402-4301-8885-e1fea58412fa
# ╠═446d43b8-77bf-410a-b3e0-0cf250d329d4
# ╠═82b9e895-db48-4503-b66d-78ec178b0420
# ╠═0831c71b-4979-43d3-a77a-966f4586c3b6
# ╠═72aa35ff-af19-4ee0-b857-d582f4df70e8
# ╠═18e0ee0f-882f-4aed-8e86-4906942037d8
# ╠═c56d43bb-3f9b-4e71-8789-44e2ea1a9282
# ╠═1c60bcff-5a60-41c4-bd12-3580a8e87680
# ╠═e08bc378-f23e-4ac4-b680-aa7133ab8a28
