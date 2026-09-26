### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ dadef766-4d7a-11f1-b6f6-dd96007b2989
function make_onset_model(; sfreq, overlap_strength = 1.0, predictor_overlap = true)

    # overlap_strength = 0 -> low overlap, long SOA
    # overlap_strength = 1 -> high overlap, short SOA
    target_mean_ms =
        (1 - overlap_strength) * 900.0 +
        overlap_strength * 250.0

    target_mean_samples = target_mean_ms / 1000 * sfreq

    σ = 0.35
    μ0 = log(target_mean_samples) - (σ^2) / 2

    if predictor_overlap
        μ_β = [
            μ0,
            0.10 * overlap_strength,
            0.05 * overlap_strength,
        ]
    else
        μ_β = [
            μ0,
            0.0,
            0.0,
        ]
    end

    return LogNormalOnsetFormula(
        μ_formula = @formula(0 ~ 1 + condition + continuous),
        μ_β = μ_β,
        σ_β = [σ],
        offset_β = [0.0],
        truncate_upper = nothing,
    )
end

# ╔═╡ ae251a7c-f0f0-47db-93d1-a8c162be2e23


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.3"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╠═dadef766-4d7a-11f1-b6f6-dd96007b2989
# ╠═ae251a7c-f0f0-47db-93d1-a8c162be2e23
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
