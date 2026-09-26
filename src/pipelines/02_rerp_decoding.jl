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
    design = [
        "stimulus" => (
            formula,
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = cfg.sfreq
            ))]
    
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
    model = make_ridge_tuned_model(),
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
    model = make_ridge_tuned_model(),
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