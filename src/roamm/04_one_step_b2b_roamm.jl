function fit_one_step_b2b_roamm(
    eeg,
    evts;
    cross_val_reps::Int = 5
)

    formula_b2b = @formula(0 ~ 1 + word_length + frequency_zipf + word_surprisal)
    design = [
        "fixation" => (
            formula_b2b,
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = 256.0
            ))]
    solver = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps, multithreading = false, show_progress = false)
    return Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design,
        evts,
        eeg;
        solver = solver
    )
end