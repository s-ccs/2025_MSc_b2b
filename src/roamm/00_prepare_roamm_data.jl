function concatenate_runs(
    subject, 
    all_events; 
    eeg_dir
)
    eeg_parts = Matrix{Float64}[]
    events_parts = DataFrame[]

    sample_offset = 0
    for run in 1:5
        # Load continuous EEG for this run 
        eeg_path = joinpath(eeg_dir, "$(subject)_run$(run)_eeg.npy")
        eeg = Float64.(NPZ.npzread(eeg_path)) .* 1e6

        # Selecyt events belonging to this run
        evts = filter(row -> row.subject_id == subject && row.run_num == run, all_events)

        # Check run-local event latencies (alrerady 1-based)
        @assert all((evts.latency .>= 1) .& (evts.latency .<= size(eeg, 2)))

        # Convert run-local latencies to subject-level latencies
        evts.latency .+= sample_offset 
        push!(eeg_parts, eeg)
        push!(events_parts, evts)
        sample_offset += size(eeg, 2)
    end

    # Concatenate after loading all five runs    
    eeg_combined = hcat(eeg_parts...)
    events_combined = vcat(events_parts...)
    sort!(events_combined, :latency)
    @assert size(eeg_combined, 2) == sample_offset
    
    return (eeg_combined, events_combined)
end

function prepare_roamm_epochs(
    eeg,
    evts;
    sfreq = 256.0,
    epoch_window = (-0.1, 1.0)
)
    dat_epoched, times = Unfold.epoch(eeg, evts, epoch_window, sfreq)
    evts_epoched, dat_epoched = Unfold.drop_missing_epochs(evts, dat_epoched)
    return dat_epoched, evts_epoched, times
end

function _center_predictors!(evts, predictors)
	# prepare the centered numerical predictor values
	centered = Symbol[]

	for predictor in predictors
		predictor in propertynames(evts) || 
			error("Predictor $(predictor) not found in event table.")
		x = evts[!, predictor]
		any(ismissing, x) && error("Predictor $(predictor) contains missing values.")
		x = Float64.(x)
		new_name = Symbol(predictor, "_c")
		evts[!, new_name] = x .- mean(x)
		push!(centered, new_name)
	end
	return centered
end 

function _zscore_predictors!(evts, predictors)
	# prepare the z-scored numerical predictor values
	zscored = Symbol[]
	for predictor in predictors
		predictor in propertynames(evts) || 
			error("Predictor $(predictor) not found in event table.")
		x = evts[!, predictor]
		any(ismissing, x) && error("Predictor $(predictor) contains missing values.")
		x = Float64.(x)
		new_name = Symbol(predictor, "_z")
		evts[!, new_name] = (x .- mean(x)) ./ std(x)
		push!(zscored, new_name)
	end
	return zscored
end

function _prepare_predictors!(evts, predictors; scaling = :center)
    prepared = Symbol[]
    for predictor in predictors
        x = Float64.(evts[!, predictor])

        if scaling == :center
            x_new = x .- mean(x)
            suffix = "_c"
        elseif scaling == :zscore
            std(x) > 0 || error("Zero variance: $predictor")
            x_new = (x .- mean(x)) ./ std(x)
            suffix = "_z"
        else
            error("Unknown scaling: $scaling")
        end
        new_name = Symbol(predictor, suffix)
        evts[!, new_name] = x_new
        push!(prepared, new_name)
    end
    return prepared
end

function _make_formula(predictors)
	# construct the regression formula
	rhs = foldl(
		+,
		StatsModels.term.(predictors);
		init = StatsModels.term(1)
	)
	return StatsModels.term(0) ~ rhs
end