""" Plain B2B on raw epoched EEG"""
function fit_plain_b2b_case(
    case_data;
    cross_val_reps::Int
)
    dat = case_data.epoched
    evts = case_data.events_epoched
    times = case_data.times

    @assert ndims(dat) == 3
    @assert size(dat, 3) == nrow(evts)
    @assert size(dat, 2) == length(times)

    formula_b2b = @formula(0 ~ 1 + condition + continuous)

    design_plain_b2b = ["stimulus" => (formula_b2b, times)]
    solver_b2b = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps)

    return Unfold.fit(UnfoldDecode.UnfoldModel, design_plain_b2b, evts, dat; solver = solver_b2b)
end


""" Fit all four simulation cases with plain B2B. """
function run_plain_b2b(
    cases;
    cross_val_reps::Int
)
    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for case_name in (:clean, :overlap, :confound, :both)
        case_data = getproperty(cases, case_name)

        uf_plain_b2b = fit_plain_b2b_case(
            cfg,
            case_data;
            cross_val_reps = cross_val_reps
        )

        result = DataFrame(Unfold.coeftable(uf_plain_b2b))

        # keep raw B2B estimate for debugging
		result[!, :estimate_signed] = copy(result.estimate)

		# B2B has no sign
		result.estimate .= abs.(result.estimate)

		# remove intercept
        result = result[result.coefname .!= "(Intercept)", :]
        
        result[!, :case] = fill(String(case_name), nrow(result))

        push!(score_tables, result)
        fitted_models[case_name] = uf_plain_b2b
    end

    return (
        score_tables = vcat(score_tables...),
        fitted_models = fitted_models
    )
end


