### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 8925770f-fa3b-440b-8ff1-9a29ee70d93e
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 1019212a-9a2b-11f1-a030-175a0c1065d8
begin
	using Serialization
	using DataFrames
	using CairoMakie
	
	project_root = dirname(Base.active_project())
	
	one_step_scores = deserialize(
		joinpath(project_root, "results", "one_step_scores.jls")
	)
end

# ╔═╡ 49fc436f-194d-407a-93e0-b794a4c6d00b
first(one_step_scores, 10)

# ╔═╡ 2c8ab646-0667-47e2-9dda-375c5345b583
names(one_step_scores)

# ╔═╡ 827d6f10-9153-400f-8cb4-86b4ec23b053
unique(string.(one_step_scores.coefname))

# ╔═╡ 06b0526b-e462-4da2-8afd-c791438ec763
unique(string.(one_step_scores.case))

# ╔═╡ ef0a4d30-326b-4ce8-a95e-9e9aa25132a1
extrema(one_step_scores.time)

# ╔═╡ 28523aa3-9e20-4dd3-9f9c-f38668f8709a
function plot_one_step_b2b(
    scores;
    target::Symbol = :condition,
    use_abs::Bool = false,
)

    coefnames = string.(scores.coefname)

    mask =
        if target == :condition
            occursin.("condition", coefnames)
        elseif target == :continuous
            occursin.("continuous", coefnames)
        else
            error("target must be :condition or :continuous")
        end

    df = scores[mask, :]

    case_order = ["clean", "overlap", "confound", "both"]

    colors = Dict(
        "clean" => :dodgerblue,
        "overlap" => :darkorange,
        "confound" => :seagreen,
        "both" => :deeppink,
    )

    titles = Dict(
        "clean" => "Clean",
        "overlap" => "Overlap",
        "confound" => "Confound",
        "both" => "Both",
    )

    fig = Figure(size = (1000, 700))

    Label(
        fig[0, 1:2],
        target == :condition ?
            "One-step FIR + B2B — Condition" :
            "One-step FIR + B2B — Continuous",
        fontsize = 24,
    )

    for (i, case_name) in enumerate(case_order)

        row = i <= 2 ? 1 : 2
        col = i <= 2 ? i : i - 2

        ax = Axis(
            fig[row, col],
            title = titles[case_name],
            xlabel = row == 2 ? "Time [s]" : "",
            ylabel = col == 1 ? "B2B estimate" : "",
            topspinevisible = false,
            rightspinevisible = false,
        )

        subdf = df[string.(df.case) .== case_name, :]

        if nrow(subdf) == 0
            @warn "No data found" case_name target
            continue
        end

        time = Float64.(subdf.time)
        estimate = Float64.(subdf.estimate)

        if use_abs
            estimate = abs.(estimate)
        end

        order = sortperm(time)

        lines!(
            ax,
            time[order],
            estimate[order],
            color = colors[case_name],
            linewidth = 2.5,
        )

        hlines!(
            ax,
            [0.0],
            color = (:gray, 0.5),
            linestyle = :dash,
        )

        xlims!(ax, -0.1, 1.0)
    end

    return fig
end

# ╔═╡ 048db3a4-6f3d-4bd3-ba5b-abaad06a72cc
fig_one_step_condition = plot_one_step_b2b(
    one_step_scores;
    target = :condition,
    use_abs = false,
)


# ╔═╡ f86115d0-e0ed-48a1-9881-db4d22d389b4
fig_one_step_continuous = plot_one_step_b2b(
    one_step_scores;
    target = :continuous,
    use_abs = false,
)

# ╔═╡ c5355c7b-628b-4524-adc1-243dd2a5eee7
begin
	mkpath("figures")
	
	save(
	    "figures/one_step_condition.svg",
	    fig_one_step_condition,
	)
	
	save(
	    "figures/one_step_continuous.svg",
	    fig_one_step_continuous,
	)
end

# ╔═╡ Cell order:
# ╠═8925770f-fa3b-440b-8ff1-9a29ee70d93e
# ╠═1019212a-9a2b-11f1-a030-175a0c1065d8
# ╠═49fc436f-194d-407a-93e0-b794a4c6d00b
# ╠═2c8ab646-0667-47e2-9dda-375c5345b583
# ╠═827d6f10-9153-400f-8cb4-86b4ec23b053
# ╠═06b0526b-e462-4da2-8afd-c791438ec763
# ╠═ef0a4d30-326b-4ce8-a95e-9e9aa25132a1
# ╠═28523aa3-9e20-4dd3-9f9c-f38668f8709a
# ╠═048db3a4-6f3d-4bd3-ba5b-abaad06a72cc
# ╠═f86115d0-e0ed-48a1-9881-db4d22d389b4
# ╠═c5355c7b-628b-4524-adc1-243dd2a5eee7
