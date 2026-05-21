# src/Components.jl
# Ground-truth ERP component utilities.
# This file is included inside the B2BSim module.

export make_components,
       ground_truth_waveforms

"""
    make_components(; sfreq = 20,
                     use_positive_kernels = true,
                     condition_intercept = 5.0,
                     condition_effect = 3.0,
                     continuous_intercept = 5.0,
                     continuous_effect = 1.0)

Create two LinearModelComponents:

1. An early component modulated by `condition`
2. A late component modulated by `continuous`

By default, `use_positive_kernels = true` applies `abs.` to the N170/P300 basis.
This avoids interpreting B2B signs in the first simulation.
"""
function make_components(;
    sfreq::Real = 20,
    use_positive_kernels::Bool = true,
    condition_intercept::Real = 5.0,
    condition_effect::Real = 3.0,
    continuous_intercept::Real = 5.0,
    continuous_effect::Real = 1.0,
)
    early_basis = n170(; sfreq = sfreq)
    late_basis = p300(; sfreq = sfreq)

    if use_positive_kernels
        early_basis = abs.(early_basis)
        late_basis = abs.(late_basis)
    end

    condition_component = LinearModelComponent(;
        basis = early_basis,
        β = [condition_intercept, condition_effect],
        formula = @formula(0 ~ 1 + condition),
    )

    continuous_component = LinearModelComponent(;
        basis = late_basis,
        β = [continuous_intercept, continuous_effect],
        formula = @formula(0 ~ 1 + continuous),
    )

    return [condition_component, continuous_component]
end

"""
    ground_truth_waveforms(; sfreq = 20,
                            use_positive_kernels = true,
                            condition_intercept = 5.0,
                            condition_effect = 3.0,
                            continuous_intercept = 5.0,
                            continuous_effect = 1.0,
                            low_continuous = 1.0,
                            high_continuous = 4.0)

Return a NamedTuple with simple ground-truth waveforms for sanity-check plots.
This is useful in Pluto notebooks or plotting scripts.
"""
function ground_truth_waveforms(;
    sfreq::Real = 20,
    use_positive_kernels::Bool = true,
    condition_intercept::Real = 5.0,
    condition_effect::Real = 3.0,
    continuous_intercept::Real = 5.0,
    continuous_effect::Real = 1.0,
    low_continuous::Real = 1.0,
    high_continuous::Real = 4.0,
)
    early_basis = n170(; sfreq = sfreq)
    late_basis = p300(; sfreq = sfreq)

    if use_positive_kernels
        early_basis = abs.(early_basis)
        late_basis = abs.(late_basis)
    end

    t_early = collect(0:length(early_basis)-1) ./ sfreq
    t_late = collect(0:length(late_basis)-1) ./ sfreq

    low_condition = condition_intercept .* early_basis
    high_condition = (condition_intercept + condition_effect) .* early_basis

    low_continuous_wave =
        (continuous_intercept + continuous_effect * low_continuous) .* late_basis
    high_continuous_wave =
        (continuous_intercept + continuous_effect * high_continuous) .* late_basis

    return (
        t_early = t_early,
        t_late = t_late,
        low_condition = low_condition,
        high_condition = high_condition,
        low_continuous = low_continuous_wave,
        high_continuous = high_continuous_wave,
    )
end
