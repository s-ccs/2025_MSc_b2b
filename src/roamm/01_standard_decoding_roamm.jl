function fit_standard_decoding_roamm(
    dat,
    evts,
    times;
    model = make_ridge_tuned_model(),
    targets = (:word_length, :frequency_zipf, :word_surprisal),
    nfolds::Int = 5, 
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