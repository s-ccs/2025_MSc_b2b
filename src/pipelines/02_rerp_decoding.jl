function _fit_rerp_case(
    cfg,
    case_data,
    model,
    formula;
    target::Symbol,
    nfolds::Int
)

    # continuous EEG + correspoding events
    dat_cont = case_data.continuous
    evts = case_data.events_continuous

    # rERP / FIR analysis design
    design = ["stimulus" => (formula, Unfold.firbasis(τ = [-0.1, 1.0], sfreq = cfg.sfreq))]
    
    return Unfold.fit(
        UnfoldDecodingModel,
        design,
        evts,
        dat_cont,
        model,
        "stimulus" => target;
        nfolds = nfolds,
        predict_type = Continuous,
        eventcolumn = :event,
        multithreading = false
    )
end


function _run_rerp_decoding(
    cfg,
    simulation,
    cases,
    formula;
    model,
    targets,
    nfolds::Int 
)
    score_tables = DataFrame[]
    fitted_models = Dict{Tuple{Symbol, Symbol}, Any}()

    for case_name in cases
        case_data = getproperty(simulation, case_name)

        for target in targets
            uf_rerp = _fit_rerp_case(
                cfg,
                case_data,
                model,
                formula;
                target = target,
                nfolds = nfolds
            )

            score_case = DataFrame(
                Unfold.coeftable(
                    uf_rerp;
                    measure = Statistics.cor,
                    averaged = true
                )
            )

            # rename the estimate column to r for correlation
            rename!(score_case, :estimate => :r)


            score_case[!, :case] = fill(String(case_name), nrow(score_case))
            score_case[!, :target] = fill(String(target), nrow(score_case))
            push!(score_tables, score_case)
            fitted_models[(case_name, target)] = uf_rerp
        end
    end

    return(
        scores = vcat(score_tables...),
        models = fitted_models
    )
end


function run_rerp_decoding(
    cfg::ConditionContinuousConfig,
    simulation;
    model = UnfoldDecode.tunemodel(
        RidgeRegressor(fit_intercept = true);
        nfolds = 5,
        resolution = 10,
        measure = MLJ.rms
    ),
    target::Symbol = :continuous,
    nfolds::Int = 3
)
    formula = @formula(0 ~ 1 + condition + continuous)
    
    return _run_rerp_decoding(
        cfg,
        simulation,
        (:clean, :overlap, :confound, :both),
        formula;
        model = model,
        targets = (target,),
        nfolds = nfolds
    )
end

function run_rerp_decoding(
    cfg::CorrelatedContinuousConfig,
    simulation;
    model = UnfoldDecode.tunemodel(
        RidgeRegressor(fit_intercept = true);
        nfolds = 5,
        resolution = 10,
        measure = MLJ.rms
    ),
    targets = (
        :continuous1,
        :continuous2,
        :continuous3
    ),
    nfolds::Int = 3
)
    formula = @formula(0 ~ 1 + continuous1 + continuous2 + continuous3)
    
    return _run_rerp_decoding(
        cfg,
        simulation,
        (:no_overlap, :overlap),
        formula;
        model = model,
        targets = targets,
        nfolds = nfolds
    )
end

function fit_rerp_decoding_roamm(
    eeg,
    evts;
    # model = make_ridge_tuned_model(),
    model = UnfoldDecode.tunemodel(
        RidgeRegressor(fit_intercept = true);
        nfolds = 5,
        resolution = 10,
        measure = MLJ.rms
    ),
    targets = (:word_length, :frequency_zipf, :word_surprisal),
    nfolds::Int = 5, 
    sfreq::Float64 = 256.0,
    scaling::Symbol = :center
)
    evts = copy(evts)
    model_predictors = _prepare_predictors!(evts, collect(targets); scaling = scaling)
    design = [
        "fixation" => (
            _make_formula(model_predictors),
            Unfold.firbasis(τ = [-0.1, 1.0], sfreq = sfreq)
        )
    ]

    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for (target, model_target) in zip(targets, model_predictors)
        println("START decoding: $target")
        flush(stdout)
        uf_rerp = Unfold.fit(
            UnfoldDecodingModel,
            design,
            evts,
            eeg,
            model,
            "fixation" => model_target;
            nfolds = nfolds,
            predict_type = Continuous,
            eventcolumn = :event,
            multithreading = false,
            unfold_fit_options = (; show_progress = false)
        )
        println("FINISHED decoding $target")

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