### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 8342e65e-a91b-4124-9f4f-9ebcb501fbe2
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 50b45037-f786-445d-a739-150f5bf58d53
using PlutoLinks: @revise

# ╔═╡ 9fae2e81-09cf-479e-9b5a-c58f8f415cd1
using PlutoLinks: @ingredients

# ╔═╡ ec67fe10-3d5d-4900-a0b2-11c783f15912
begin
	using DataFrames
	using StatsModels
	using Unfold
	using UnfoldDecode
	using CairoMakie
	using UnfoldMakie
	using UnfoldSim
	using PlutoUI
end

# ╔═╡ f54f1b3a-44ac-4224-94d5-f4014a87d22a
import PlutoLinks

# ╔═╡ c270273c-96fb-464a-94ee-4484a879a478
@revise using MScB2B

# ╔═╡ 67cfe216-3b20-4bee-a063-7b607f7626f9
Controls = @ingredients "simulation_controls.jl"

# ╔═╡ 2bd9a33e-cdb4-4110-905d-d7eead0e64e2
typeof(Controls)

# ╔═╡ 454cce81-787c-42e5-bafb-d0cb791c8ced
Controls_test =
    PlutoLinks.ingredients(
        joinpath(@__DIR__, "simulation_controls.jl")
    )

# ╔═╡ 2ab62963-8107-4e3c-b995-80486df369c5
@bind sim Controls.simulation_controls()

# ╔═╡ 11750070-aa2a-4467-bc97-45f832515637
begin
	cross_val_reps = 5
	
	b2b_solver = (x, y) ->
	    UnfoldDecode.solver_b2b(
	        x,
	        y;
	        cross_val_reps = cross_val_reps,
	    )
end

# ╔═╡ 0779b840-29c1-42b5-837f-d22b30aa70e6
cfg_debug = MScB2B.SimulationConfig(
    n_trials = 2000,
    β_continuous = 3.0,
)

# ╔═╡ b99e7578-98a6-11f1-9b83-7fef937e3d49
begin
	cfg = MScB2B.SimulationConfig()
	
	epoched_cases = MScB2B.simulate_cases(
	    cfg_debug;
	    continuous = false,
	)
end

# ╔═╡ a366f720-7487-4098-8839-0580e8ac80a5
dat_clean, evts_clean = epoched_cases.clean

# ╔═╡ 13b4ab8f-d155-49d8-a90d-171ea1397285
(
    data_size = size(dat_clean),
    n_events = nrow(evts_clean),
    columns = names(evts_clean),
)

# ╔═╡ a52933c2-d729-45eb-aee1-e1c64c329061
begin
	fo_b2b = @formula(0 ~ 1 + condition + continuous)
	
	n_timepoints = size(dat_clean, 2)
	
	times_b2b =
	    (0:n_timepoints-1) ./ cfg.sfreq
	
	des_b2b_plain = [
	    Any => (
	        fo_b2b,
	        times_b2b,
	    )
	]
end

# ╔═╡ 317247f0-959e-4b30-9ff3-bdfb53826757
length(times_b2b) == size(dat_clean, 2)

# ╔═╡ 4a36a7c4-be6f-4737-818e-dab8db3fc5aa
uf_b2b_clean = Unfold.fit(
    UnfoldDecode.UnfoldModel,
    des_b2b_plain,
    evts_clean,
    dat_clean;
    solver = b2b_solver,
)

# ╔═╡ f099e1fb-d010-46e5-918a-6241ff8c929d
begin
	b2b_clean = coeftable(uf_b2b_clean)

b2b_clean.estimate =
    abs.(b2b_clean.estimate)

b2b_clean =
    b2b_clean[
        b2b_clean.coefname .!= "(Intercept)",
        :
    ]

plot_erp(
    b2b_clean;
    mapping = (; color = :coefname),
)
end

# ╔═╡ bdb1a270-3d60-4cd8-8067-f253b4f24b97
unique(b2b_clean.coefname)

# ╔═╡ a36dbbc5-c734-484d-831b-2d0e868ca798
first(b2b_clean, 100)

# ╔═╡ 70614002-5858-4807-8c82-ab7f79ef43e2
b2b_clean_nointercept = filter(
    :coefname => !=("(Intercept)"),
    b2b_clean,
)

# ╔═╡ 365ac1db-590c-43ff-8e4e-ab28a25b3d1f
unique(b2b_clean_nointercept.channel)

# ╔═╡ 50012a5e-76bb-43a9-a0c5-641c407d4cb7
begin
	fig = Figure()
	
	ax = Axis(
	    fig[1, 1],
	    xlabel = "Time (s)",
	    ylabel = "B2B estimate",
	    title = "Plain B2B — clean",
	)
	
	for coef in ["condition: face", "continuous"]
	    tmp = filter(
	        :coefname => ==(coef),
	        b2b_clean_nointercept,
	    )
	
	    lines!(
	        ax,
	        tmp.time,
	        abs.(tmp.estimate),
	        label = coef,
	    )
	end
	
	axislegend(ax)
	
	fig
end

# ╔═╡ c7884355-6cfb-4292-9bc1-eaed388526c8
begin
	condition_rows =
	    b2b_clean[b2b_clean.coefname .== "condition: face", :]
	
	continuous_rows =
	    b2b_clean[b2b_clean.coefname .== "continuous", :]
	
	condition_peak_time =
	    condition_rows.time[argmax(condition_rows.estimate)]
	
	continuous_peak_time =
	    continuous_rows.time[argmax(continuous_rows.estimate)]
	
	(condition_peak_time, continuous_peak_time)
end

# ╔═╡ 55cced88-76e4-4174-8604-aa0414fc71f1
begin
    n170 = UnfoldSim.n170(; sfreq = cfg.sfreq)
    p300 = UnfoldSim.p300(; sfreq = cfg.sfreq)

    t_n170 = (0:length(n170)-1) ./ cfg.sfreq
    t_p300 = (0:length(p300)-1) ./ cfg.sfreq

    (
        n170_peak = t_n170[argmin(n170)],
        p300_peak = t_p300[argmax(p300)],
        n170_length = length(n170),
        p300_length = length(p300),
    )
end

# ╔═╡ Cell order:
# ╠═8342e65e-a91b-4124-9f4f-9ebcb501fbe2
# ╠═f54f1b3a-44ac-4224-94d5-f4014a87d22a
# ╠═50b45037-f786-445d-a739-150f5bf58d53
# ╠═9fae2e81-09cf-479e-9b5a-c58f8f415cd1
# ╠═67cfe216-3b20-4bee-a063-7b607f7626f9
# ╠═2bd9a33e-cdb4-4110-905d-d7eead0e64e2
# ╠═454cce81-787c-42e5-bafb-d0cb791c8ced
# ╠═2ab62963-8107-4e3c-b995-80486df369c5
# ╠═c270273c-96fb-464a-94ee-4484a879a478
# ╠═ec67fe10-3d5d-4900-a0b2-11c783f15912
# ╠═b99e7578-98a6-11f1-9b83-7fef937e3d49
# ╠═a366f720-7487-4098-8839-0580e8ac80a5
# ╠═13b4ab8f-d155-49d8-a90d-171ea1397285
# ╠═a52933c2-d729-45eb-aee1-e1c64c329061
# ╠═11750070-aa2a-4467-bc97-45f832515637
# ╠═317247f0-959e-4b30-9ff3-bdfb53826757
# ╠═4a36a7c4-be6f-4737-818e-dab8db3fc5aa
# ╠═bdb1a270-3d60-4cd8-8067-f253b4f24b97
# ╠═a36dbbc5-c734-484d-831b-2d0e868ca798
# ╠═70614002-5858-4807-8c82-ab7f79ef43e2
# ╠═365ac1db-590c-43ff-8e4e-ab28a25b3d1f
# ╠═50012a5e-76bb-43a9-a0c5-641c407d4cb7
# ╠═f099e1fb-d010-46e5-918a-6241ff8c929d
# ╠═c7884355-6cfb-4292-9bc1-eaed388526c8
# ╠═55cced88-76e4-4174-8604-aa0414fc71f1
# ╠═0779b840-29c1-42b5-837f-d22b30aa70e6
