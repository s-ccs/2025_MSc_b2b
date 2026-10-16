MLJ.machine(
    model::RidgeRegressor,       # Adapter: allow MLJ RidgeRegressor to accept matrix-valued inputs
    X::AbstractMatrix{Float64},  # by converting trials x channels matrices to an MLJ table
    y::AbstractVector{Float64};
    kwargs...,
) = MLJ.machine(model,
    MLJ.table(X),
    y;
    kwargs...,
)

function make_ridge_tuned_model(;
    inner_nfolds::Int = 3, 
    resolution::Int = 20
)
    ridge_model = RidgeRegressor()

    lambda_range = range(
        ridge_model,
        :lambda;
        lower = 1e-6,
        upper = 1e2,
        scale = :log
    )

    return TunedModel(
        model = ridge_model,
        resampling = CV(nfolds = inner_nfolds),
        tuning = Grid(resolution = resolution),
        range = lambda_range,
        measure = RSquared()
    )
end


function decode_standard_ridge(
    dat,
    evts,
    times,
    model;
    target::Symbol,
    nfolds::Int,
    seed::Int
)
    # dat: channels x timepoints x trials
    @assert length(times) == size(dat, 2)
    @assert ndims(dat) == 3
    @assert size(dat, 3) == nrow(evts)

    y = Float64.(evts[!, target])
    n_trials = size(dat, 3)
    n_timepoints = size(dat, 2)

    folds = 
        UnfoldDecode.MLJBase.train_test_pairs(
            UnfoldDecode.MLJBase.CV(
                nfolds = nfolds,
                shuffle = true,
                rng = seed
            ),
            1:n_trials
        )
    
    # trials x timepoints
    yhat = fill(NaN, n_trials, n_timepoints)

    for (train_indices, test_indices) in folds
        machines = UnfoldDecode.fit_timepoints(
            model,
            @view(dat[:, :, train_indices]),
            @view(y[train_indices])
        )

        predictions = UnfoldDecode.predict_timepoints(
            machines,
            @view(dat[:, :, test_indices])
        )

        for t in eachindex(predictions)
            yhat[test_indices, t] .= Float64.(predictions[t])
        end
    end

    scores = DataFrame(
        timepoint = collect(1:n_timepoints),
        time = Float64.(times),
        # time = (collect(1:n_timepoints) .- 1) ./ cfg.sfreq,
        r = [cor(yhat[:, t], y) for t in 1:n_timepoints],
        r_squared = [RSquared()(yhat[:, t], y) for t in 1:n_timepoints]
    )

    return scores, yhat
end

function fit_standard_decoding_roamm(
    dat,
    evts,
    times;
    model = make_ridge_tuned_model(),
    targets = (:word_length, :frequency_zipf, :word_surprisal),
    nfolds:Int = 5, 
    seed::Int = 12
)
    score_tables = DataFrame[]
    yhats = Dict{Symbol, Matrix{Float64}}()

    for target in targets
        scores, yhat = decode_standard_ridge(
            dat,
            evts,
            times,
            model;
            target = target,
            nfolds = nfolds,
            seed = seed
        )
        scores[!, :target] = fill(string(target), nrow(scores))
        push!(score_tables, scores)
        yhats[target] = yhat
    end

    return(
        scores = vcat(score_tables...),
        yhats = yhats
    )
end