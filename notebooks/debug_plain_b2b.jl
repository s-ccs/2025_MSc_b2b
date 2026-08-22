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
	using UnfoldMakie
	using UnfoldSim
	using PlutoUI
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

# ╔═╡ 70cc5cde-5242-42e8-93b6-c2b573d66282
plain_b2b =
    MScB2B.run_plain_b2b(
        cfg,
        cases;
        cross_val_reps = 3,
    )

# ╔═╡ 9407e5ce-d6af-44f2-a253-847e51be6c1c
fig_plain_b2b_condition_no_p300 = MScB2B.plot_b2b_grid(
    plain_b2b.score_tables,
    cfg;
    target = :condition,
)

# ╔═╡ 4eb6bb0e-5dad-4948-85ab-1724428be360
fig_plain_b2b_continuous_no_p300 = MScB2B.plot_b2b_grid(
    plain_b2b.score_tables,
    cfg;
    target = :continuous,
)

# ╔═╡ 88aef0ba-9559-49e3-b4b2-6f12fa93d2ca
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

# ╔═╡ c1da4a39-547c-4a07-bac0-1d85e403aacf
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

# ╔═╡ 847dc5b8-6edc-4ba9-a50a-57e05d9f8ee7
clean_1ch = merge(
    cases.clean,
    (
        epoched = cases.clean.epoched[1:1, :, :],
    ),
)

# ╔═╡ 3fc3091b-108e-4ce6-92ab-c54197f3e415
b2b_clean_1ch =
    MScB2B.fit_plain_b2b_case(
        cfg,
        clean_1ch;
        cross_val_reps = 3,
    )

# ╔═╡ 17b8912f-d170-436a-84dc-dc0774ed17b0
b2b_clean_1ch_scores =
    DataFrame(
        Unfold.coeftable(b2b_clean_1ch)
    )

# ╔═╡ 9bc89901-09d0-4550-88ee-4c85dde242ec
b2b_clean_1ch_condition =
    b2b_clean_1ch_scores[
        occursin.(
            "condition",
            string.(b2b_clean_1ch_scores.coefname),
        ),
        :
    ]

# ╔═╡ f7415676-6be7-4aa9-a5de-7cfb9981da12
lines(
    b2b_clean_1ch_condition.time,
    abs.(b2b_clean_1ch_condition.estimate);
    axis = (
        xlabel = "Time [s]",
        ylabel = "|B2B estimate|",
        title = "Clean condition B2B — one channel",
    ),
)

# ╔═╡ 670671ac-fabd-460b-9f5a-d39ccd677676
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

# ╔═╡ 4ba29ad3-419a-42c4-b44c-220893ea48ff
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

# ╔═╡ 63526444-05f0-42d4-ad11-09fdfb846a77
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

# ╔═╡ Cell order:
# ╠═8342e65e-a91b-4124-9f4f-9ebcb501fbe2
# ╠═7b55f4b4-c299-43eb-b125-7d4219cef384
# ╠═c270273c-96fb-464a-94ee-4484a879a478
# ╠═ec67fe10-3d5d-4900-a0b2-11c783f15912
# ╠═15b3af3d-4ae7-4838-8f12-6bd40dafa6b6
# ╠═765521e2-bf87-4dba-9dfd-b986819e1961
# ╠═72658d97-4c75-4d0e-bd15-e1329457de3b
# ╠═70cc5cde-5242-42e8-93b6-c2b573d66282
# ╠═9407e5ce-d6af-44f2-a253-847e51be6c1c
# ╠═c1da4a39-547c-4a07-bac0-1d85e403aacf
# ╠═9d9c41e2-56cd-4742-b166-68af0a1b3ea1
# ╠═4eb6bb0e-5dad-4948-85ab-1724428be360
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
