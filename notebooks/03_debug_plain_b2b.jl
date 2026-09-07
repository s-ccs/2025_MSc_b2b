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
	Pkg.instantiate()
end

# ╔═╡ 7b55f4b4-c299-43eb-b125-7d4219cef384
using PlutoLinks: @revise, @ingredients

# ╔═╡ c270273c-96fb-464a-94ee-4484a879a478
@revise using MScB2B

# ╔═╡ ec67fe10-3d5d-4900-a0b2-11c783f15912
begin
	using DataFrames
	using StatsModels
	using Unfold
	using UnfoldDecode
	using CairoMakie
	using Serialization
	using UnfoldMakie
	using UnfoldSim
	using PlutoUI
	using Statistics
	
end

# ╔═╡ 52de7d3b-a5ac-42b0-86dc-eaffdcc1877b
begin
	project_root = normpath(joinpath(@__DIR__, ".."))
	
	results_dir = joinpath(project_root, "results")
end

# ╔═╡ 15b3af3d-4ae7-4838-8f12-6bd40dafa6b6
Controls = @ingredients(joinpath(@__DIR__, "simulation_controls.jl"))

# ╔═╡ 9d9c41e2-56cd-4742-b166-68af0a1b3ea1
@bind sim Controls.simulation_controls()

# ╔═╡ 765521e2-bf87-4dba-9dfd-b986819e1961
cfg = MScB2B.SimulationConfig(
    n_trials = sim.n_trials,

    noiselevel = sim.noiselevel,
    channel_noise_sd = sim.channel_noise_sd,

    β0_n170 = sim.β0_n170,
    β0_p300 = sim.β0_p300,

    β_condition = sim.β_condition,
    β_continuous = sim.β_continuous,

    rho = sim.rho,

    overlap_interval_ms = sim.overlap_interval_ms,
    onset_condition_bias = sim.onset_condition_bias,
    shift_onset = sim.shift_onset,
)

# ╔═╡ 72658d97-4c75-4d0e-bd15-e1329457de3b
cases = MScB2B.simulate_cases(cfg)

# ╔═╡ c1da4a39-547c-4a07-bac0-1d85e403aacf
# ╠═╡ disabled = true
#=╠═╡
begin
	save(
	    "figures/fig_plain_b2b_condition_no_p300.svg",
	     fig_plain_b2b_condition_no_p300,
	)
	
	save(
		"figures/clean_condition_signed_b2b_no_p300.svg",
		clean_condition_signed_b2b_no_p300,
	)

	save(
		"figures/fig_plain_b2b_continuous_no_p300.svg",
		fig_plain_b2b_continuous_no_p300,
	)

end
  ╠═╡ =#

# ╔═╡ 58c8a219-a62f-41db-98d9-3933348d9584


# ╔═╡ 60fb6ecf-3eb0-47dd-a81f-9ca8ac88433f


# ╔═╡ f1b1ca38-f028-4bcf-af18-cf70eef13192


# ╔═╡ f914bad6-000f-4370-8da9-da9faf688a00


# ╔═╡ 7b5be90f-1611-4de4-9e96-d67216d7224c


# ╔═╡ ed63ae9b-292f-4e4e-96fd-5dee0f86d5b8


# ╔═╡ d5adc048-c72b-4fd5-a97f-87111964f4fd


# ╔═╡ 6695d158-6cda-47f3-a3fd-a89782b5d9a2


# ╔═╡ 71c19822-b1fd-4572-ac7d-d949bc0263ca


# ╔═╡ c4da28c0-77cf-4dfd-bc0b-efbf7b44065e


# ╔═╡ a8023e32-4fac-451f-9c98-f9063a7caceb


# ╔═╡ 36a589c0-969e-4803-9e6f-4f1bc3b27c84


# ╔═╡ 7c28bd0a-084d-4f64-804c-51749d8dcd71


# ╔═╡ b309cfa1-596f-450b-8ae7-6e6ead10ba9a


# ╔═╡ a6dc5ccc-c992-414f-a7b0-d76f909377e6
uf_clean = MScB2B.fit_plain_b2b_case(
    cases.clean;
    cross_val_reps = 3,
)

# ╔═╡ d31c92d0-09f9-483e-a6e4-8a9c155ee52d
# ╠═╡ disabled = true
#=╠═╡
begin
	res_clean = coeftable(uf_clean)

	UnfoldMakie.plot_erp(
	    res_clean;
	    mapping = (; color = :coefname),
	)
end
  ╠═╡ =#

# ╔═╡ 8bce183b-421b-4eea-85ec-d7e0afcb5cc4


# ╔═╡ 9187fc47-924a-4438-8b00-4415f30d8f7b
# ╠═╡ disabled = true
#=╠═╡
begin
	md"""### Null control
	
	β_condition = 0  
	β_continuous = 0  
	β0_n170 = 0  
	β0_p300 = 0  
	ρ = 0
	
	Expected:
	No systematic B2B recovery.
	
	Result:
	Condition and continuous estimates fluctuate around zero.
	✓ Passed
	"""
	
	
	res_clean = coeftable(uf_clean)

	UnfoldMakie.plot_erp(
		res_clean;
		mapping = (; color = :coefname),
	)
	
end
  ╠═╡ =#

# ╔═╡ 6f21484d-e2e7-49f9-bdba-7c6ab6361c26
# ╠═╡ disabled = true
#=╠═╡
begin
	md"""
	β_condition = 3
	Shift onset by one = false
	Onset condition bias = 0
	"""


	res_clean = coeftable(uf_clean)
	#res_clean.estimate = abs.(res_clean.estimate)
	res_clean = filter(:coefname => !=("(Intercept)"), res_clean)
	
	UnfoldMakie.plot_erp(
	    res_clean;
	    mapping = (; color = :coefname),
	    axis = (; ylabel = "B2B estimate"),
	)
end
  ╠═╡ =#

# ╔═╡ 84e1c5c4-3bc7-456e-b0f1-6a9e10ef2435
# ╠═╡ disabled = true
#=╠═╡
begin
	md"""
	β_condition = 0
	β_continuous = 3   
	
	β0_n170 = 0
	β0_p300 = 0
	
	rho = 0
	onset_condition_bias = 0
	shift_onset = false
	"""
	res_clean  = coeftable(uf_clean)
	res_clean.estimate = abs.(res_clean.estimate)
	res_clean = filter(:coefname => !=("(Intercept)"), res_clean)
	
	UnfoldMakie.plot_erp(
	    res_clean;
	    mapping = (; color = :coefname),
	    axis = (; ylabel = "B2B estimate"),
	)
end
	
  ╠═╡ =#

# ╔═╡ 2b5be1b0-9dbe-4b11-9585-c22f5fc2cb98
# ╠═╡ disabled = true
#=╠═╡
begin
	md"""
	β_condition = 3
	β_continuous = 3
	rho = 0
	
	β0_n170 = 0
	β0_p300 = 0
	onset_condition_bias = 0
	shift_onset = false
	"""
	res_clean  = coeftable(uf_clean)
	res_clean.estimate = abs.(res_clean.estimate)
	res_clean = filter(:coefname => !=("(Intercept)"), res_clean)
	
	UnfoldMakie.plot_erp(
	    res_clean;
	    mapping = (; color = :coefname),
	    axis = (; ylabel = "B2B estimate"),
	)
end
	
  ╠═╡ =#

# ╔═╡ 9713e635-4dc0-4156-a799-184f37de4063
# ╠═╡ disabled = true
#=╠═╡
begin
	md"""
	β_condition = 3
	β_continuous = 3
	rho = 0.8
	
	β0_n170 = 0
	β0_p300 = 0
	onset_condition_bias = 0
	shift_onset = false
	"""
	res_clean  = coeftable(uf_clean)
	res_clean.estimate = abs.(res_clean.estimate)
	res_clean = filter(:coefname => !=("(Intercept)"), res_clean)
	
	UnfoldMakie.plot_erp(
	    res_clean;
	    mapping = (; color = :coefname),
	    axis = (; ylabel = "B2B estimate"),
	)
end
  ╠═╡ =#

# ╔═╡ d43594a5-af65-460d-bacf-287398e1340b
begin
	md"""
	β_condition = 3
	β_continuous = 3
	rho = 0.8
	
	β0_n170 = 0
	β0_p300 = 0
	onset_condition_bias = 0
	shift_onset = false
	"""
	res_clean  = coeftable(uf_clean)
	res_clean.estimate = abs.(res_clean.estimate)
	res_clean = filter(:coefname => !=("(Intercept)"), res_clean)
	
	UnfoldMakie.plot_erp(
	    res_clean;
	    mapping = (; color = :coefname),
	    axis = (; ylabel = "B2B estimate"),
	)
end

# ╔═╡ c33d32db-f97a-4a76-bc5b-ace5f9635f55
extrema(res_clean.time)

# ╔═╡ 8a4cee87-23ba-4670-8e3b-d0a15de7a178


# ╔═╡ ad54eaec-a509-4695-b4ff-04aa22075321


# ╔═╡ f9ec1294-72f8-4756-a3f9-11bdb1204189


# ╔═╡ b8319361-15cc-45fc-9b03-8c5593e78d54


# ╔═╡ 79a4c4cc-9de9-4b3f-a091-88d28cdc22cc
# ╠═╡ disabled = true
#=╠═╡
let
	p300_original = UnfoldSim.p300(; sfreq = cfg.sfreq)

	p300_hanning = UnfoldSim.hanning(
		0.30,
		0.30,
		cfg.sfreq
	)
	maximum(abs.(p300_original .- p300_hanning))
end
  ╠═╡ =#

# ╔═╡ 70cc5cde-5242-42e8-93b6-c2b573d66282
plain_b2b =
    MScB2B.run_plain_b2b(
        cases;
        cross_val_reps = 3,
    )

# ╔═╡ 3cbbb808-c581-4b63-9715-d3c16668a464
fig_plain_b2b_condition = MScB2B.plot_b2b_grid(
	plain_b2b.score_tables,
	cfg;
	target = :condition,
	x_window = (-0.1, 0.6),
)

# ╔═╡ 48ab8ae0-7ad8-4ae2-b439-a933f2b61538
fig_plain_b2b_continuous = MScB2B.plot_b2b_grid(
    plain_b2b.score_tables,
    cfg;
    target = :continuous,
	x_window = (-0.1, 0.6),
)


# ╔═╡ 5ebd807f-1640-4b8a-96ee-8faada6a2f17
# ╠═╡ disabled = true
#=╠═╡
begin
	save(
		"03_plain_b2b_condition_hanning_1500trials_long_windows.svg",
		fig_plain_b2b_condition
	)

	save(
		"03_plain_b2b_continuous_hanning_1500trials_long_windows.svg",
		fig_plain_b2b_continuous
	)
end
  ╠═╡ =#

# ╔═╡ 88aef0ba-9559-49e3-b4b2-6f12fa93d2ca
# ╠═╡ disabled = true
#=╠═╡
begin
    dat = cases.clean.epoched
    evts = cases.clean.events_epoched
    times = collect(cases.clean.times)

    face_ix = evts.condition .== "face"
    car_ix  = evts.condition .== "car"

    # average over channels and trials
    erp_face = vec(StatsModels.mean(dat[:, :, face_ix]; dims = (1, 3)))
    erp_car  = vec(StatsModels.mean(dat[:, :, car_ix];  dims = (1, 3)))

    empirical_condition_effect = erp_face .- erp_car

    lines(
        times,
        empirical_condition_effect;
        axis = (
            xlabel = "Time [s]",
            ylabel = "Face - car",
            title = "Empirical condition effect in clean EEG",
        ),
    )
end
  ╠═╡ =#

# ╔═╡ b96998a9-1385-4bb1-a793-fef9b312dc8b
begin
    b2b_condition_clean =
        plain_b2b.score_tables[
            (string.(plain_b2b.score_tables.case) .== "clean") .&
            occursin.(
                "condition",
                string.(plain_b2b.score_tables.coefname),
            ),
            :
        ]

    sort!(b2b_condition_clean, :time)

    peak_ix =
        argmax(abs.(b2b_condition_clean.estimate))

    b2b_condition_clean[
        max(1, peak_ix - 5):min(nrow(b2b_condition_clean), peak_ix + 5),
        [:time, :estimate],
    ]
end

# ╔═╡ c4905ea9-4f24-43d4-bdda-7d3d0ec2f1e9
clean_condition_signed_b2b_no_p300 = lines(
    b2b_condition_clean.time,
    b2b_condition_clean.estimate;
    axis = (
        xlabel = "Time [s]",
        ylabel = "B2B estimate",
        title = "Clean condition B2B — signed",
    ),
)

# ╔═╡ 847dc5b8-6edc-4ba9-a50a-57e05d9f8ee7
clean_1ch = merge(
    cases.clean,
    (
        epoched = cases.clean.epoched[1:1, :, :],
    ),
)

# ╔═╡ 3fc3091b-108e-4ce6-92ab-c54197f3e415
# ╠═╡ disabled = true
#=╠═╡
b2b_clean_1ch =
    MScB2B.fit_plain_b2b_case(
        cfg,
        clean_1ch;
        cross_val_reps = 3,
    )
  ╠═╡ =#

# ╔═╡ 17b8912f-d170-436a-84dc-dc0774ed17b0
#=╠═╡
b2b_clean_1ch_scores =
    DataFrame(
        Unfold.coeftable(b2b_clean_1ch)
    )
  ╠═╡ =#

# ╔═╡ 9bc89901-09d0-4550-88ee-4c85dde242ec
# ╠═╡ disabled = true
#=╠═╡
b2b_clean_1ch_condition =
    b2b_clean_1ch_scores[
        occursin.(
            "condition",
            string.(b2b_clean_1ch_scores.coefname),
        ),
        :
    ]
  ╠═╡ =#

# ╔═╡ f7415676-6be7-4aa9-a5de-7cfb9981da12
# ╠═╡ disabled = true
#=╠═╡
lines(
    b2b_clean_1ch_condition.time,
    abs.(b2b_clean_1ch_condition.estimate);
    axis = (
        xlabel = "Time [s]",
        ylabel = "|B2B estimate|",
        title = "Clean condition B2B — one channel",
    ),
)
  ╠═╡ =#

# ╔═╡ 670671ac-fabd-460b-9f5a-d39ccd677676
# ╠═╡ disabled = true
#=╠═╡
begin
    reps_to_test = [3, 10, 50]

    fig = Figure(size = (800, 500))
    ax = Axis(
        fig[1, 1];
        xlabel = "Time [s]",
        ylabel = "|B2B estimate|",
        title = "Effect of cross_val_reps",
    )

    for reps in reps_to_test

        model = MScB2B.fit_plain_b2b_case(
            cfg,
            cases.clean;
            cross_val_reps = reps,
        )

        tbl = DataFrame(
            Unfold.coeftable(model)
        )

        tbl = tbl[
            occursin.(
                "condition",
                string.(tbl.coefname),
            ),
            :
        ]

        lines!(
            ax,
            tbl.time,
            abs.(tbl.estimate);
            label = "reps = $reps",
        )
    end

    axislegend(ax)

    fig
end
  ╠═╡ =#

# ╔═╡ 4ba29ad3-419a-42c4-b44c-220893ea48ff
# ╠═╡ disabled = true
#=╠═╡
begin
    formula_test =
        @formula(0 ~ 1 + condition)

    design_test = [
        Any => (
            formula_test,
            cases.clean.times,
        ),
    ]

    solver_test = (X, y) ->
        UnfoldDecode.solver_b2b(
            X,
            y;
            cross_val_reps = 50,
        )

    b2b_condition_only =
        Unfold.fit(
            UnfoldDecode.UnfoldModel,
            design_test,
            cases.clean.events_epoched,
            cases.clean.epoched;
            solver = solver_test,
        )

    tbl_condition_only =
        DataFrame(
            Unfold.coeftable(
                b2b_condition_only
            )
        )

    tbl_condition_only =
        tbl_condition_only[
            occursin.(
                "condition",
                string.(tbl_condition_only.coefname),
            ),
            :
        ]

    lines(
        tbl_condition_only.time,
        abs.(tbl_condition_only.estimate);
        axis = (
            xlabel = "Time [s]",
            ylabel = "|B2B estimate|",
            title = "B2B: condition-only design",
        ),
    )
end
  ╠═╡ =#

# ╔═╡ 63526444-05f0-42d4-ad11-09fdfb846a77
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
begin
    dat_check = cases.clean.epoched
    evts_check = cases.clean.events_epoched
    times_check = cases.clean.times

    solver = (X, y) ->
        UnfoldDecode.solver_b2b(
            X,
            y;
            cross_val_reps = 50,
        )

    designs = Dict(
        "current: condition + continuous" =>
            @formula(0 ~ 1 + condition + continuous),

        "condition only" =>
            @formula(0 ~ 1 + condition),

        "condition_num only" =>
            @formula(0 ~ 1 + condition_num),
    )

    fig_check = Figure(size = (900, 550))

    ax_check = Axis(
        fig_check[1, 1];
        xlabel = "Time [s]",
        ylabel = "|B2B estimate|",
        title = "B2B design comparison — clean condition",
    )

    for (label, formula) in designs

        design = [
            Any => (
                formula,
                times_check,
            ),
        ]

        model = Unfold.fit(
            UnfoldDecode.UnfoldModel,
            design,
            evts_check,
            dat_check;
            solver = solver,
        )

        tbl = DataFrame(
            Unfold.coeftable(model)
        )

        mask =
            occursin.(
                "condition",
                string.(tbl.coefname),
            ) .&
            .!occursin.(
                "Intercept",
                string.(tbl.coefname),
            )

        sub = tbl[mask, :]

        lines!(
            ax_check,
            sub.time,
            abs.(sub.estimate);
            label = label,
        )
    end

    axislegend(ax_check)

    fig_check
end
  ╠═╡ =#

# ╔═╡ c30bbd00-a849-446b-a6ad-a841e984fe7b
let
    dat = cases.clean.epoched
    evts = cases.clean.events_epoched
    times = collect(cases.clean.times)

    face_ix = evts.condition .== "face"
    car_ix = evts.condition .== "car"

    grand_mean =
        vec(StatsModels.mean(dat; dims = (1, 3)))

    condition_effect =
        vec(StatsModels.mean(dat[:, :, face_ix]; dims = (1, 3))) .-
        vec(StatsModels.mean(dat[:, :, car_ix]; dims = (1, 3)))

    fig_signal = Figure(size = (850, 550))

    ax1 = Axis(
        fig_signal[1, 1];
        xlabel = "Time [s]",
        ylabel = "EEG amplitude",
        title = "Clean EEG: common signal vs condition effect",
    )

    lines!(
        ax1,
        times,
        grand_mean;
        label = "Grand mean",
    )

    lines!(
        ax1,
        times,
        condition_effect;
        label = "Face - car",
    )

    hlines!(
        ax1,
        [0.0];
        linestyle = :dash,
    )

    axislegend(ax1)
	
    fig_signal

	save(
		"clean_common_vs_condition_no_p300.svg", 
		fig_signal,)
end

# ╔═╡ c85cc3f8-d4cb-432a-bb7a-30530cd4d0ef
let
    n170 = UnfoldSim.n170(; sfreq = cfg.sfreq)
    p300 = UnfoldSim.p300(; sfreq = cfg.sfreq)

    n170_common = cfg.β0_n170 .* n170
    p300_common = cfg.β0_p300 .* p300

    t_n170 = (0:length(n170)-1) ./ cfg.sfreq
    t_p300 = (0:length(p300)-1) ./ cfg.sfreq

    fig = Figure(size = (800, 500))

    ax = Axis(
        fig[1, 1];
        xlabel = "Time [s]",
        ylabel = "Amplitude",
        title = "Simulated common ERP components",
    )

    lines!(
        ax,
        t_n170,
        n170_common;
        label = "N170 intercept",
    )

    lines!(
        ax,
        t_p300,
        p300_common;
        label = "P300 intercept",
    )

    hlines!(ax, [0.0]; linestyle = :dash)

    axislegend(ax)

    fig
end

# ╔═╡ 8d931c76-2cf0-4ba6-b165-3a2201e0b4db
let
    dat = cases.clean.epoched

    dat_centered = dat .- StatsModels.mean(dat; dims = 3)

    (
        original_size = size(dat),
        centered_size = size(dat_centered),
        original_grand_mean = maximum(abs.(StatsModels.mean(dat; dims = 3))),
        centered_grand_mean = maximum(abs.(StatsModels.mean(dat_centered; dims = 3))),
    )
end

# ╔═╡ 6ff06cb2-391c-4778-8c07-72b1c3a64fc5
clean_centered = let
    dat = cases.clean.epoched
    dat_centered = dat .- StatsModels.mean(dat; dims = 3)

    merge(
        cases.clean,
        (; epoched = dat_centered),
    )
end

# ╔═╡ e0495416-0894-43bf-a03d-cfa076e03c8d
uf_clean_centered = MScB2B.fit_plain_b2b_case(
    cfg,
    clean_centered;
    cross_val_reps = 3,
)

# ╔═╡ 38d6762c-fabd-4763-b72e-7ac90936321e
tbl_clean_centered = Unfold.coeftable(uf_clean_centered)

# ╔═╡ 970b86a0-01cb-4e6d-a490-f82371a8062d
let
    # -------------------------
    # Original clean B2B
    # -------------------------
    uf_original = MScB2B.fit_plain_b2b_case(
        cfg,
        cases.clean;
        cross_val_reps = 3,
    )

    tbl_original = DataFrame(Unfold.coeftable(uf_original))
    tbl_original[!, :data] .= "Original"

    # -------------------------
    # Trial-centered clean B2B
    # -------------------------
    dat = cases.clean.epoched
    dat_centered = dat .- StatsModels.mean(dat; dims = 3)

    clean_centered = merge(
        cases.clean,
        (; epoched = dat_centered),
    )

    uf_centered = MScB2B.fit_plain_b2b_case(
        cfg,
        clean_centered;
        cross_val_reps = 3,
    )

    tbl_centered = DataFrame(Unfold.coeftable(uf_centered))
    tbl_centered[!, :data] .= "Trial-centered"

    # -------------------------
    # Combine + condition only
    # -------------------------
    tbl = vcat(tbl_original, tbl_centered)

    tbl_condition = filter(
        row -> occursin("condition", String(row.coefname)),
        tbl,
    )

    # -------------------------
    # Plot
    # -------------------------
    fig = Figure(size = (800, 450))
    ax = Axis(
        fig[1, 1],
        title = "Clean condition B2B: original vs trial-centered",
        xlabel = "Time [s]",
        ylabel = "B2B estimate",
    )

    for label in ["Original", "Trial-centered"]
        sub = filter(row -> row.data == label, tbl_condition)

        lines!(
            ax,
            sub.time,
            sub.estimate;
            label = label,
            linewidth = 2,
        )
    end

    axislegend(ax)

    fig
end

# ╔═╡ 679b95bd-8417-4c02-97af-81fb56874142


# ╔═╡ Cell order:
# ╠═8342e65e-a91b-4124-9f4f-9ebcb501fbe2
# ╠═7b55f4b4-c299-43eb-b125-7d4219cef384
# ╠═52de7d3b-a5ac-42b0-86dc-eaffdcc1877b
# ╠═c270273c-96fb-464a-94ee-4484a879a478
# ╠═ec67fe10-3d5d-4900-a0b2-11c783f15912
# ╠═15b3af3d-4ae7-4838-8f12-6bd40dafa6b6
# ╠═765521e2-bf87-4dba-9dfd-b986819e1961
# ╠═72658d97-4c75-4d0e-bd15-e1329457de3b
# ╠═9d9c41e2-56cd-4742-b166-68af0a1b3ea1
# ╠═3cbbb808-c581-4b63-9715-d3c16668a464
# ╠═48ab8ae0-7ad8-4ae2-b439-a933f2b61538
# ╠═5ebd807f-1640-4b8a-96ee-8faada6a2f17
# ╠═c1da4a39-547c-4a07-bac0-1d85e403aacf
# ╠═58c8a219-a62f-41db-98d9-3933348d9584
# ╠═60fb6ecf-3eb0-47dd-a81f-9ca8ac88433f
# ╠═f1b1ca38-f028-4bcf-af18-cf70eef13192
# ╠═f914bad6-000f-4370-8da9-da9faf688a00
# ╠═7b5be90f-1611-4de4-9e96-d67216d7224c
# ╠═ed63ae9b-292f-4e4e-96fd-5dee0f86d5b8
# ╠═d5adc048-c72b-4fd5-a97f-87111964f4fd
# ╠═6695d158-6cda-47f3-a3fd-a89782b5d9a2
# ╠═71c19822-b1fd-4572-ac7d-d949bc0263ca
# ╠═c4da28c0-77cf-4dfd-bc0b-efbf7b44065e
# ╠═a8023e32-4fac-451f-9c98-f9063a7caceb
# ╠═36a589c0-969e-4803-9e6f-4f1bc3b27c84
# ╠═7c28bd0a-084d-4f64-804c-51749d8dcd71
# ╠═b309cfa1-596f-450b-8ae7-6e6ead10ba9a
# ╠═a6dc5ccc-c992-414f-a7b0-d76f909377e6
# ╠═d31c92d0-09f9-483e-a6e4-8a9c155ee52d
# ╠═c33d32db-f97a-4a76-bc5b-ace5f9635f55
# ╠═8bce183b-421b-4eea-85ec-d7e0afcb5cc4
# ╠═9187fc47-924a-4438-8b00-4415f30d8f7b
# ╠═6f21484d-e2e7-49f9-bdba-7c6ab6361c26
# ╠═84e1c5c4-3bc7-456e-b0f1-6a9e10ef2435
# ╠═2b5be1b0-9dbe-4b11-9585-c22f5fc2cb98
# ╠═9713e635-4dc0-4156-a799-184f37de4063
# ╠═d43594a5-af65-460d-bacf-287398e1340b
# ╠═8a4cee87-23ba-4670-8e3b-d0a15de7a178
# ╠═ad54eaec-a509-4695-b4ff-04aa22075321
# ╠═f9ec1294-72f8-4756-a3f9-11bdb1204189
# ╠═b8319361-15cc-45fc-9b03-8c5593e78d54
# ╠═79a4c4cc-9de9-4b3f-a091-88d28cdc22cc
# ╠═70cc5cde-5242-42e8-93b6-c2b573d66282
# ╠═88aef0ba-9559-49e3-b4b2-6f12fa93d2ca
# ╠═b96998a9-1385-4bb1-a793-fef9b312dc8b
# ╠═c4905ea9-4f24-43d4-bdda-7d3d0ec2f1e9
# ╠═847dc5b8-6edc-4ba9-a50a-57e05d9f8ee7
# ╠═3fc3091b-108e-4ce6-92ab-c54197f3e415
# ╠═17b8912f-d170-436a-84dc-dc0774ed17b0
# ╠═9bc89901-09d0-4550-88ee-4c85dde242ec
# ╠═f7415676-6be7-4aa9-a5de-7cfb9981da12
# ╠═670671ac-fabd-460b-9f5a-d39ccd677676
# ╠═4ba29ad3-419a-42c4-b44c-220893ea48ff
# ╠═63526444-05f0-42d4-ad11-09fdfb846a77
# ╠═c30bbd00-a849-446b-a6ad-a841e984fe7b
# ╠═c85cc3f8-d4cb-432a-bb7a-30530cd4d0ef
# ╠═8d931c76-2cf0-4ba6-b165-3a2201e0b4db
# ╠═6ff06cb2-391c-4778-8c07-72b1c3a64fc5
# ╠═e0495416-0894-43bf-a03d-cfa076e03c8d
# ╠═38d6762c-fabd-4763-b72e-7ac90936321e
# ╠═970b86a0-01cb-4e6d-a490-f82371a8062d
# ╠═679b95bd-8417-4c02-97af-81fb56874142
