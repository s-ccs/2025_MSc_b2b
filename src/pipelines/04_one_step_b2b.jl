function UnfoldDecode.solver_b2b(                     # when T is a number,
    X, # design matrix
    data::AbstractArray{T,2}; 
    cross_val_reps = 10,
    multithreading = true,
    show_progress = true,
    solver_G = UnfoldDecode.model_ridge,
    solver_H = UnfoldDecode.model_ridge,
) where {T<:Number}

    #@show model_fun

    E = zeros(size(X, 2), size(X, 2))
    #W = Array{Float64}(undef, size(X, 2), size(data, 1))
    prog = UnfoldDecode.Progress(cross_val_reps; dt = 0.1, enabled = show_progress)
    #@maybe_threads multithreading for t = 1:size(data, 2)

	UnfoldDecode.@maybe_threads multithreading for m = 1:cross_val_reps
		k_ix = collect(UnfoldDecode.Kfold(size(data, 2), 2))
		Y1 = data[:,  k_ix[1]] # views here made it much slower
		Y2 = data[:,  k_ix[2]]
		X1 = X[k_ix[1], :]
		X2 = X[k_ix[2], :]

		Y2G = solver_G(Y1', X1, Y2')
		H = solver_H(X2, Y2G)

		E[:, :] = E[:, :] + UnfoldDecode.Diagonal(H[UnfoldDecode.diagind(H)])
		UnfoldDecode.ProgressMeter.next!(prog; showvalues = [(:cross_val_rep, m)])
	end
	E[:, :] = E[:, :] ./ cross_val_reps
	W = (X * E[ :, :])' / data[:, :]

    # extract diagonal
    beta = UnfoldDecode.diag(E)
    # reshape to conform to ch x time x pred
    #beta = permutedims(beta, [3 1 2])
	beta  = reshape(beta, 1,:)
    modelinfo = Dict("W" => W, "E" => E, "cross_val_reps" => cross_val_reps) # no history implemented (yet?)

    # LinearModelFit{T, 2} requires a 2D array for beta
    return UnfoldDecode.Unfold.LinearModelFit{eltype(beta),2}(estimate = beta, 
    info = modelinfo, standarderror = similar(beta, 0, 0))
end


function _fit_one_step_b2b_case(
    cfg,
    case_data,
    formula_b2b;
    cross_val_reps::Int
)

    dat_cont = case_data.continuous
    evts = case_data.events_continuous

    design_one_step_b2b = [
        "stimulus" => (
            formula_b2b,
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = cfg.sfreq
            ))]
    solver_b2b = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps,)
    return Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design_one_step_b2b,
        evts,
        dat_cont;
        solver = solver_b2b
    )
end

function _run_one_step_b2b(
    cfg,
    simulation,
    cases,
    formula_b2b;
    cross_val_reps::Int
)   
    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for case_name in cases

        println("Starting one-step B2B: $case_name")
        flush(stdout)

        case_data = getproperty(simulation, case_name)

        uf_one_step_b2b = _fit_one_step_b2b_case(
            cfg,
            case_data,
            formula_b2b;
            cross_val_reps = cross_val_reps
        )

        result = DataFrame(Unfold.coeftable(uf_one_step_b2b))

        # keep raw B2B estimate for debugging
		result[!, :estimate_signed] = copy(result.estimate)

		# B2B has no sign
		result.estimate .= abs.(result.estimate)

		# remove intercept
        result = result[result.coefname .!= "(Intercept)", :]
        
        result[!, :case] = fill(String(case_name), nrow(result))

        push!(score_tables, result)
        fitted_models[case_name] = uf_one_step_b2b

        println("Finished one-step B2B: $case_name")
        flush(stdout)
    end
    return (
        score_tables = vcat(score_tables...),
        fitted_models = fitted_models
    )
end

function run_one_step_b2b(
    cfg::ConditionContinuousConfig,
    simulation;
    cross_val_reps::Int = 3
)
    formula_b2b = @formula(0 ~ 1 + condition + continuous)

    return _run_one_step_b2b(
        cfg,
        simulation,
        (:clean, :overlap, :confound, :both),
        formula_b2b;
        cross_val_reps = cross_val_reps
    )
end

function run_one_step_b2b(
    cfg::CorrelatedContinuousConfig,
    simulation;
    cross_val_reps::Int = 3
)
    formula_b2b = @formula(0 ~ 1 + continuous1 + continuous2 + continuous3)

    return _run_one_step_b2b(
        cfg,
        simulation,
        (:no_overlap, :overlap),
        formula_b2b;
        cross_val_reps = cross_val_reps
    )
end

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

function fit_one_step_b2b_roamm(
    eeg,
    evts;
    predictos::Vector{Symbol} = [:word_length, :frequency_zipf, :word_surprisal],
    sfreq::Float64 = 256.0,
    cross_val_reps::Int = 5,
    scaling::Symbol = :center
) 

    evts = copy(evts)
    sort!(evts, :latency)

    # Unfold event type
    evts[!, :event] = fill("fixation", nrow(evts))

    # Same perdictor preprocessing as two-step B2B
    model_predictors = _prepare_predictors!(evts, predictos; scaling = scaling)

    formula_b2b = _make_formula(model_predictors)

    design = [
        "fixation" => (
            formula_b2b,
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = sfreq
            ))]

    solver = (X, y) -> UnfoldDecode.solver_b2b(
        X,
        y;
        cross_val_reps = cross_val_reps,
        multithreading = false,
        show_progress = false,
        solver_G = UnfoldDecode.model_ridge,
        solver_H = UnfoldDecode.model_ridge
    )

    t0 = time()

    uf_b2b = Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design,
        evts,
        eeg;
        solver = solver,
        eventcolumn = :event
    )

    println("One-step B2B time: ", round(time() - t0; digits=2), " s")

    return(
        predictors = predictos,
        model_predictors = model_predictors,
        fitted_model = uf_b2b
    )
end

function calculate_one_step_b2b_R(
    X_fir,
    eeg;
    nfolds::Int = 5,
    cross_val_reps::Int = 5,
    sfreq::Float64 = 256.0,
    seed::Int = 12
)
    # X_fir: samples x FIR-expanded predictors
    # eeg: channels x samples

    X = Matrix{Float64}(X_fir)
    Y = Matrix{Float64}(permutedims(eeg))

    n_samples, n_predictors = size(X)

    # Outer CV: contiguous blocks
    edges = round.(Int, range(0, n_samples; length = nfolds + 1))

    fold_R = fill(NaN, nfolds, n_predictors)
    fold_R_knockout = fill(NaN, nfolds, n_predictors)
    fold_delta_R = fill(NaN, nfolds, n_predictors)
    fold_R_XH_true = fill(NaN, nfolds, n_predictors)

    for fold_idx in 1:nfolds
        println("Outer fold $fold_idx / $nfolds")  
        
        test_start = edges[fold_idx] + 1
        test_end = edges[fold_idx + 1]
        test_idx = collect(test_start:test_end)

        train_idx = [
            collect(1:test_start - 1);
            collect(test_end + 1:n_samples)
        ]

        X_train = X[train_idx, :]  
        Y_train = Y[train_idx, :]

        splits = _make_internal_splits(
            1:length(train_idx);
            cross_val_reps = cross_val_reps,
            seed = seed + 1000 * fold_idx
        )

        models = _fit_b2b(X_train, Y_train, splits)

        # Evaluate  on held-out samples
        scores = _r_scores(
            X[test_idx, :],
            Y[test_idx, :],
            models
        )

        fold_R[fold_idx, :] .= scores.R
        fold_R_knockout[fold_idx, :] .= scores.R_knockout
        fold_delta_R[fold_idx, :] .= scores.delta_R
        fold_R_XH_true[fold_idx, :] .= scores.R_XH_true
    end

    function mean_finite(scores)
        return[
            let vals = filter(isfinite, scores[;, i])
                isempty(vals) ? NaN : mean(vals)
            end
            for i in 1:n_predictors
        ]
    end

    return(
        R = mean_finite(fold_R),
        R_knockout = mean_finite(fold_R_knockout),
        delta_R = mean_finite(fold_delta_R),
        R_XH_true = mean_finite(fold_R_XH_true),
        fold_R = fold_R,
        fold_R_knockout = fold_R_knockout,
        fold_delta_R = fold_delta_R,
        fold_R_XH_true = fold_R_XH_true
    )
end