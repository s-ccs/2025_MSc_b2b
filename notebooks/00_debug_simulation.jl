### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ d8e43b36-9d02-4d65-88cd-c05a1089c9d7
begin
    import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 7d73fa7c-9591-11f1-ad7b-3fbe41a4c4d8
using PlutoLinks: @revise

# ╔═╡ e1f400d6-d913-4fc5-a547-f1447b7787b3
begin
	using Statistics
	using PlutoUI
end

# ╔═╡ 4147db38-b343-4af9-9c6e-4bfb4cda5c13
@revise using MScB2B

# ╔═╡ fea975cf-9492-42bc-8b3a-53d633763080
cfg = MScB2B.SimulationConfig()

# ╔═╡ 5fac3bbb-5a64-43cf-bf7e-ea9c05edce44
begin
    (
        active_project = Base.active_project(),
        package_path = pathof(MScB2B),
        config_type = typeof(cfg),
    )
end

# ╔═╡ 6efcdb91-39e7-47a4-84d5-bdf240a9e316
epoched_cases =
    MScB2B.simulate_cases(
        cfg;
        continuous = false,
    )

# ╔═╡ d34cfeb3-3f9a-4aef-91b7-915170d9325e
continuous_cases =
    MScB2B.simulate_cases(
        cfg;
        continuous = false,
    )

# ╔═╡ 445fe497-8d67-4733-9dfd-693e53ba8d76
begin
	size(epoched_cases.clean[1])
	size(continuous_cases.clean[1])
	
	first(epoched_cases.clean[2], 5)
	first(continuous_cases.clean[2], 5)
end

# ╔═╡ 519261a4-908a-464a-b286-b75732faab74
begin


clean_evts = epoched_cases.clean[2]
confound_evts = epoched_cases.confounded[2]

cor(clean_evts.condition_num, clean_evts.continuous)
cor(confound_evts.condition_num, confound_evts.continuous)
end

# ╔═╡ Cell order:
# ╠═5fac3bbb-5a64-43cf-bf7e-ea9c05edce44
# ╠═d8e43b36-9d02-4d65-88cd-c05a1089c9d7
# ╠═7d73fa7c-9591-11f1-ad7b-3fbe41a4c4d8
# ╠═e1f400d6-d913-4fc5-a547-f1447b7787b3
# ╠═4147db38-b343-4af9-9c6e-4bfb4cda5c13
# ╠═fea975cf-9492-42bc-8b3a-53d633763080
# ╠═6efcdb91-39e7-47a4-84d5-bdf240a9e316
# ╠═d34cfeb3-3f9a-4aef-91b7-915170d9325e
# ╠═445fe497-8d67-4733-9dfd-693e53ba8d76
# ╠═519261a4-908a-464a-b286-b75732faab74
