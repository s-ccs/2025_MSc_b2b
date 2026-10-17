function fit_rerp_decoding_roamm(
    eeg,
    evts;
    model = make_ridge_tuned_model(),
    targets = (:word_length, :frequency_zipf, :word_surprisal),
    nfolds::Int = 5, 
    sfreq::Float64 = 256.0
)

    design = [
        "fixation" => (
            @formula(0 ~ 1 + word_length + frequency_zipf + word_surprisal),
            Unfold.firbasis(τ = [-0.1, 1.0], sfreq = sfreq)
        )
    ]

    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for target in targets
        uf_rerp = Unfold.fit(
            UnfoldDecodingModel,
            design,
            evts,
            eeg,
            model,
            "fixation" => target;
            nfolds = nfolds,
            predict_type = Continuous,
            eventcolumn = :event,
            multithreading = false
        )

        scores = DataFrame(
            Unfold.coeftable(
                uf_rerp;
                measure = Statistics.cor,
                averaged = true
            )
        )

        # rename the estimate column to r for correlation
        rename!(scores, :estimate => :r)
        scores[!, :target] = fill(String(target), nrow(scores))
        push!(score_tables, scores)
        fitted_models[target] = uf_rerp
    end

    return(
        scores = vcat(score_tables...),
        models = fitted_models
    )
end