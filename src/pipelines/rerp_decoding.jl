function fit_rerp_case(
    cfg,
    case_data,
    model;
    target::Symbol = :continuous,
    nfolds::Int = 3,
)

    # continuous EEG + correspoding events
    dat_cont = case_data.continuous
    evts = case_data.events_continuous

    # rERP / FIR analysis design
    design = [
        "stimulus" => (
            @formula(0 ~ 1 + condition + continuous),
            Unfold.firbasis(
                τ = [-0.1, 1.0],
                sfreq = cfg.sfreq,
                name = "stimulus",
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
        multithreading = false,
    )
end

function run_rerp_decoding(
    cfg,
    cases;
    model = make_ridge_tuned_model(),
    target::Symbol = :continuous,
    nfolds::Int = 3,
)
    score_tables = DataFrame[]
    fitted_models = Dict{Symbol, Any}()

    for case_name in (:clean, :overlap, :confound, :both)
        case_data = getproperty(cases, case_name)
        
        uf_rerp = fit_rerp_case(
            cfg,
            case_data,
            model;
            target = target,
            nfolds = nfolds,
        )

        score_case = DataFrame(
            Unfold.coeftable(
                uf_rerp;
                measure = Statistics.cor,
                averaged = true,
            )
        )

        # rename the estimate column to r for correlation
        rename!(score_case, :estimate => :r)

        # plotting helper currently expects :timepoint
        score_case[!, :timepoint] = collect(1:nrow(score_case))

        score_case[!, :case] = fill(String(case_name), nrow(score_case))
        push!(score_tables, score_case)
        fitted_models[case_name] = uf_rerp
    end

    return(
        scores = vcat(score_tables...),
        models = fitted_models,
    )
end