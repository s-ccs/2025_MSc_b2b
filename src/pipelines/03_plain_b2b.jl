""" Plain B2B on raw epoched EEG"""
function _fit_plain_b2b_case(
    case_data,
    formula_b2b;
    cross_val_reps::Int
)
    dat = case_data.epoched
    evts = case_data.events_epoched
    times = case_data.times

    @assert ndims(dat) == 3
    @assert size(dat, 3) == nrow(evts)
    @assert size(dat, 2) == length(times)

    design_plain_b2b = ["stimulus" => (formula_b2b, times)]
    solver_b2b = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps)

    return Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design_plain_b2b,
        evts,
        dat;
        solver = solver_b2b
        )
end


""" Fit all four simulation cases with plain B2B. """
function _run_plain_b2b(
    simulation,
    cases,
    formula_b2b;
    cross_val_reps::Int
)
    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for case_name in cases
        case_data = getproperty(simulation, case_name)

        uf_plain_b2b = _fit_plain_b2b_case(
            case_data,
            formula_b2b;
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


function run_plain_b2b(
    cfg::ConditionContinuousConfig,
    simulation;
    cross_val_reps::Int = 3
)
    formula_b2b = @formula(0 ~ 1 + condition + continuous)
    return _run_plain_b2b(
        simulation,
        (:clean, :overlap, :confound, :both),
        formula_b2b;
        cross_val_reps = cross_val_reps
    )
end

function run_plain_b2b(
    cfg::CorrelatedContinuousConfig,
    simulation;
    cross_val_reps ::Int = 3
)
    formula_b2b = @formula(0 ~ 1 + continuous1 + continuous2 + continuous3)

    return _run_plain_b2b(
        simulation,
        (:no_overlap, :overlap),
        formula_b2b;
        cross_val_reps = cross_val_reps
    )
end

function fit_plain_b2b_roamm(
    dat,
    evts,
    times;
    cross_val_reps::Int = 5,
    scaling::Symbol = :center
)

    evts = copy(evts)
    predictors = [:word_length, :frequency_zipf, :word_surprisal]
    model_predictors = _prepare_predictors!(evts, predictors; scaling = scaling)

    formula_b2b = _make_formula(model_predictors)
    design = ["fixation" => (formula_b2b, times)]

    solver = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps, show_progress = false)
    return Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design,
        evts,
        dat;
        solver = solver
    )
end 

function predictive_evaluate_plain_b2b_roamm(
    dat,
    evts,
    times;
    nfolds::Int = 5,
    cross_val_reps::Int = 5,
    scaling:: Symbol = :center,
    seed::Int = 12
)
    predictors = [:word_length, :frequency_zipf, :word_surprisal]
    X_raw = Matrix{Float64}(evts[:, predictors])

    return calculate_b2b_R(
        X_raw,
        dat,
        times;
        nfolds = nfolds,
        cross_val_reps = cross_val_reps,
        scaling = scaling,
        seed = seed
    )
end 