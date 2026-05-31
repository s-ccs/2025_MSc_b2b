predictor_formula() = @formula(0 ~ 1 + condition + continuous)

function make_components(cfg::SimConfig)
    n1 = LinearModelComponeny(;
        basis = n170(; sfreq = cfg.sfreq),
        β = [5.0, 3.0],
        formula = @formula(0 ~ 1 + condition),
    )

    p3 = LinearModelComponent(;
        basis = p300(; sfreq = cfg.sfreq),
        formula = @formula(0 ~ 1 + condinuous),
        β = [5.0, 1.0],        
    )

    return [n1, p3]
end

function ground_truth_kernels(cfg::SimConfig)
    n170_basis = n170(; sfreq = cfg.sfreq)
    p300_basis = p300(; sfreq = cfg.sfreq)

    t_n170 = collect(0:length(n170_basis) - 1) ./ cfg.sfreq
    t_p300 = collect(0:length(p300_basis) - 1) ./ cfg.sfreq

    return(
        t_n170 = t_n170,
        t_p300 = t_p300,
        n170_car = 5.0 .* n170_basis,
        n170_face = 8.0 .* n170_basis,
        n170_condition_effect = 3.0 .* n170_basis,
        p300_continous_effect = 1.0 .* p300_basis,
    )
end 