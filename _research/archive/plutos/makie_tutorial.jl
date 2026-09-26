### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ cc4d04c2-09f9-457f-8e0c-e122733a43f4
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 6c6a4c0c-b4de-4a33-aee0-b6845c820727
# ╠═╡ disabled = true
#=╠═╡
begin
	using Unfold
	using UnfoldMakie
	using DataFrames
	using CairoMakie
	using DataFramesMeta
	using UnfoldSim
end
  ╠═╡ =#

# ╔═╡ 6dc0ff3c-a6f1-11f1-834b-7b0799fdc2a6
# ╠═╡ disabled = true
#=╠═╡
begin
	using UnfoldMakie, CairoMakie, UnfoldSim,Unfold
results = Unfold.coeftable(UnfoldMakie.example_data("UnfoldLinearModel"))
f = Figure(size=(500, 350))
plot_erp!(f,
    results,
    mapping = (; row = :coefname, color = :coefname => "Conditions"),
    axis = (; xlabel = "Time [s]"),
    stderror = true,
)
f
end
  ╠═╡ =#

# ╔═╡ ef5262c2-31eb-42d4-88d1-63d06e06fcfd
#=╠═╡
plot_erp(results)
  ╠═╡ =#

# ╔═╡ e02753ae-61a4-4ce6-855d-97d2916d93ce
#=╠═╡
plot_erp(
	results,
	mapping = (; color = :coefname => "Conditions"),
	axis = (; xlabel = "Time[s]")
)
  ╠═╡ =#

# ╔═╡ 5f124d36-c1d9-4114-8f35-819f522b291c
#=╠═╡
begin
	m = UnfoldMakie.example_data("UnfoldLinearModel")
results = coeftable(m)
significancevalues = DataFrame(
    from = [0.01, 0.25],
    to = [0.2, 0.29],
    coefname = ["(Intercept)", "condition: face"], # if coefname not specified, line should be black
)
plot_erp(
    results;
    :significance => significancevalues,
    mapping = (; color = :coefname => "Conditions"),
    axis = (; xlabel = "Time [s]"),
)
end
  ╠═╡ =#

# ╔═╡ 07f19fa8-0d36-40ec-94a6-64741dc13431
#=╠═╡
begin
	erp_matrix, evts = UnfoldSim.predef_eeg(; noiselevel = 12, return_epoched = true)
	erp_matrix = reshape(erp_matrix, (1, size(erp_matrix)...))
	f = @formula 0 ~ 1 + condition + continuous
	se_solver = (x, y) -> Unfold.solver_default(x, y, stderror = true);
	
	m = fit(
	    UnfoldModel,
	    Dict(Any => (f, range(0, step = 1 / 100, length = size(erp_matrix, 2)))),
	    evts,
	    erp_matrix,
	    solver = se_solver,
	);
	results = coeftable(m)
	res_effects = effects(Dict(:continuous => -5:0.5:5), m);
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═cc4d04c2-09f9-457f-8e0c-e122733a43f4
# ╠═6dc0ff3c-a6f1-11f1-834b-7b0799fdc2a6
# ╠═6c6a4c0c-b4de-4a33-aee0-b6845c820727
# ╠═07f19fa8-0d36-40ec-94a6-64741dc13431
# ╠═ef5262c2-31eb-42d4-88d1-63d06e06fcfd
# ╠═e02753ae-61a4-4ce6-855d-97d2916d93ce
# ╠═5f124d36-c1d9-4114-8f35-819f522b291c
