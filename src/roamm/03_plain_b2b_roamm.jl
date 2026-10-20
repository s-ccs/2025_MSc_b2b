function fit_plain_b2b_roamm(
    dat,
    evts,
    times;
    cross_val_reps::Int = 5,
    scaling::Symbol = :center
)

    evts = copy(evts)
    predictors = [word_length, :frequency_zipf, :word_surprisal]
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