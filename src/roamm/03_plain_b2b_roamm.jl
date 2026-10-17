function fit_plain_b2b_roamm(
    dat,
    evts,
    times;
    cross_val_reps::Int = 5
)
    @assert ndims(dat) == 3
    @assert size(dat, 3) == nrow(evts)
    @assert size(dat, 2) == length(times)

    formula_b2b = @formula(0 ~ 1 + word_length + frequency_zipf + word_surprisal)
    design = ["fixation" => (formula_b2b, times)]

    solver = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps)
    return Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design,
        evts,
        dat;
        solver = solver
    )
end 