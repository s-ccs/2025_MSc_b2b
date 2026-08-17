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
    return UnfoldDecode.Unfold.LinearModelFit{eltype(beta),2}(beta, modelinfo)
end


function fit_one_step_b2b_case(
    cfg,
    case_data;
    cross_val_reps::Int = 3,
)

    dat_cont = case_data.continuous
    evts = case_data.events_continuous

    formula_b2b = @formula(0 ~ 1 + condition + continuous)

    design_one_step_b2b = [
        "stimulus" => (
            @formula(0 ~ 1 + condition + continuous),
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = cfg.sfreq,
                name = "stimulus",
            ))]
    solver_b2b = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps,)
    return Unfold.fit(
        UnfoldDecode.UnfoldModel,
        design,
        evts,
        dat_cont;
        solver = b2b_solver
    )

function run_one_step_b2b(
    cfg,
    cases;
    cross_val_reps::Int = 3,
)   
    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for case_name in (:clean, :overlap, :confound, :both,)
        case_data = getproperty(cases, case_name)

        uf_one_step_b2b = fit_one_step_b2b_case(
            cfg,
            case_data;
            cross_val_reps = cross_val_reps,
        )

        result = DataFrame(Unfold.coeftable(uf_one_step_b2b))
        result = result[result.coefname .!= "(intercept)", :]
        result[!, :case] = fill(String(case_name), nrow(result))

        push!(score_tables, result)
        fitted_models[case_name] = uf_b2b
    end
