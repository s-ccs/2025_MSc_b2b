function plot_b2b_grid(
    scores::AbstractDataFrame,
    cfg;
    target::Symbol,
    score_col::Symbol = :estimate,
    use_abs::Bool = true,
    guide_peak::Real = 0.45,
)

    columns = propertynames(scores)

    :case in columns ||
        error("`scores` must contain a :case column.")

    :coefname in columns ||
        error("`scores` must contain a :coefname column.")

    score_col in columns ||
        error(
            "`scores` has no column $score_col. " *
            "Available columns: $columns"
        )

    time_col =
        if :time in columns
            :time
        elseif :timepoint in columns
            :timepoint
        else
            error("`scores` must contain :time or :timepoint.")
        end

    # ------------------------------------------
    # choose the coefficient rows for this target
    # ------------------------------------------

    coefnames = string.(scores.coefname)

    target_mask =
        if target == :continuous
            occursin.("continuous", coefnames)

        elseif target in (:condition, :condition_num)
            occursin.("condition", coefnames)

        else
            error("`target` must be :continuous or :condition")
        end

    target_scores = scores[target_mask, :]

    isempty(target_scores) &&
        error(
            "No rows found for target=$target. " *
            "Available coefnames: $(unique(coefnames))"
        )

    # ------------------------------------------
    # helper: extract one case curve
    # ------------------------------------------

    function get_curve(case_name::String)
        mask = string.(target_scores.case) .== case_name

        any(mask) ||
            error(
                "No case named \"$case_name\". " *
                "Available cases: $(unique(string.(target_scores.case)))"
            )

        time =
            if time_col == :time
                Float64.(target_scores[mask, :time])
            else
                (Float64.(target_scores[mask, :timepoint]) .- 1) ./ cfg.sfreq
            end

        values = Float64.(target_scores[mask, score_col])

        if use_abs
            values = abs.(values)
        end

        order = sortperm(time)
        return time[order], values[order]
    end

    # ------------------------------------------
    # timing guides
    # ------------------------------------------

    function scale_guide(effect)
        effect_abs = abs.(Float64.(effect))
        peak = maximum(effect_abs)

        return peak == 0 ?
               zeros(length(effect_abs)) :
               guide_peak .* effect_abs ./ peak
    end

    n170_effect =
        cfg.β_condition .*
        UnfoldSim.n170(; sfreq = cfg.sfreq)

    p300_effect =
        cfg.β_continuous .*
        UnfoldSim.p300(; sfreq = cfg.sfreq)

    n170_guide = scale_guide(n170_effect)
    p300_guide = scale_guide(p300_effect)

    n170_time = (0:length(n170_guide)-1) ./ cfg.sfreq
    p300_time = (0:length(p300_guide)-1) ./ cfg.sfreq

    # ------------------------------------------
    # labels and colors
    # ------------------------------------------

    condition_target = target in (:condition, :condition_num)

    figure_title =
        condition_target ?
        "Condition B2B" :
        "Continuous B2B"

    ylabel_text =
        if use_abs
            "B2B estimate magnitude"
        else
            "B2B estimate"
        end

    case_order = ("clean", "overlap", "confound", "both")

    case_titles = Dict(
        "clean" => "Clean",
        "overlap" => "Overlap",
        "confound" => "Confound",
        "both" => "Both",
    )

    colors = Dict(
        "clean" => :dodgerblue,
        "overlap" => :darkorange,
        "confound" => :seagreen,
        "both" => :deeppink,
    )

    # ------------------------------------------
    # figure
    # ------------------------------------------

    fig = Figure(
        size = (1000, 720),
        backgroundcolor = :white,
    )

    axes = Axis[]

    for (index, case_name) in enumerate(case_order)

        row = index <= 2 ? 1 : 2
        col = isodd(index) ? 1 : 2

        show_x = row == 2
        show_y = col == 1

        ax = Axis(
            fig[row, col];
            title = case_titles[case_name],
            xlabel = show_x ? "Time [s]" : "",
            ylabel = show_y ? ylabel_text : "",
            xticklabelsvisible = show_x,
            xticksvisible = show_x,
            yticklabelsvisible = show_y,
            yticksvisible = show_y,
            xgridvisible = false,
            ygridvisible = false,
            topspinevisible = false,
            rightspinevisible = false,
        )

        time, values = get_curve(case_name)

        lines!(
            ax,
            time,
            values;
            color = colors[case_name],
            linewidth = 2.5,
        )

        # --------------------------------------
        # decide which guides to show
        # --------------------------------------

        confounded_case = case_name in ("confound", "both")

        show_n170 =
            if condition_target
                true
            else
                confounded_case
            end

        show_p300 =
            if condition_target
                confounded_case
            else
                true
            end

        if show_n170
            lines!(
                ax,
                n170_time,
                n170_guide;
                color = :black,
                linestyle = :dash,
                linewidth = 2.5,
            )
        end

        if show_p300
            lines!(
                ax,
                p300_time,
                p300_guide;
                color = :gray40,
                linestyle = :dot,
                linewidth = 2.5,
            )
        end

        hlines!(
            ax,
            [0.0];
            color = (:gray, 0.35),
            linewidth = 1,
        )

        push!(axes, ax)
    end

    linkxaxes!(axes...)
    linkyaxes!(axes...)

    Label(
        fig[0, 1:2],
        figure_title;
        fontsize = 24,
        font = :bold,
    )

    legend_elements = [
        LineElement(color = colors["clean"], linewidth = 2.5),
        LineElement(color = colors["overlap"], linewidth = 2.5),
        LineElement(color = colors["confound"], linewidth = 2.5),
        LineElement(color = colors["both"], linewidth = 2.5),
        LineElement(color = :black, linestyle = :dash, linewidth = 2.5),
        LineElement(color = :gray40, linestyle = :dot, linewidth = 2.5),
    ]

    legend_labels = [
        "Clean",
        "Overlap",
        "Confound",
        "Both",
        "N170 timing guide",
        "P300 timing guide",
    ]

    Legend(
        fig[3, 1:2],
        legend_elements,
        legend_labels;
        orientation = :horizontal,
        framevisible = false,
        nbanks = 2,
    )

    colgap!(fig.layout, 24)
    rowgap!(fig.layout, 16)

    return fig
end