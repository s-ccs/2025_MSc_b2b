function fit_standard_decoding_roamm(
    dat,
    evts,
    times;
    model = make_ridge_tuned_model(),
    targets = (:word_length, :frequency_zipf, :word_surprisal),
    nfolds::Int = 5, 
    seed::Int = 12,
    scaling::Symbol = :center
)

    evts = copy(evts)
    model_predictors = _prepare_predictors!(evts, collect(targets); scaling = scaling)

    score_tables = DataFrame[]
    yhats = Dict{Symbol, Matrix{Float64}}()

    for (target, model_target) in zip(targets, model_predictors)
        scores, yhat = decode_standard_ridge(
            dat,
            evts,
            times,
            model;
            target = model_target,
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