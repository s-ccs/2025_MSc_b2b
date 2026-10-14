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


function fit_two_step_b2b_roamm(
	dat_cont,
	evts;
	predictors::Vector{Symbol} = [:word_length],
	sfreq::Float64 = 256.0,
	cross_val_reps::Int = 3,
	scaling::Symbol = :center
)
	evts = copy(evts)
	sort!(evts, :latency)

	# Unfold needs an event typr macthing "stimulus"
	evts[!, :event] = fill("stimulus", nrow(evts))

	model_predictors = _prepare_predictors!(evts, predictors; scaling = scaling)
	formula = _make_formula(model_predictors)

	# Standarize surprisal for the smoke test
	#evts[!, :surprisal_z] = (evts.surprisal_bits .- mean(evts.surprisal_bits)) ./ std(evts.surprisal_bits)

	# Step 1: rERP 
	design_rerp = [
		"stimulus" => (
			formula,
			Unfold.firbasis(
				τ = [-0.1, 1.0],
				sfreq = sfreq
			))]
	
	uf_rerp = Unfold.fit(
		UnfoldLinearModelContinuousTime,
		design_rerp,
		evts,
		dat_cont;
		eventcolumn = :event
	)

	# Step2: recover single_trial corrected responses
	X_corrected = UnfoldDecode.singletrials(
		dat_cont,
		uf_rerp,
		evts,
		"stimulus",
		:event
	)
	println("Corrected trials size: ", size(X_corrected))

	# Step 3: trial_level B2B design, already epoched data, use a time vector
	times = Unfold.times(uf_rerp)[1]

	@assert length(times) == size(X_corrected, 2)
	
	design_two_step_b2b = [
		"stimulus" => (
			formula,
			times,
		)
	]	

	b2b_solver = (X, y) -> begin
		
		X_clean, y_clean = Unfold.drop_missing_epochs(X, y)
		X_clean = Float64.(X_clean)

		@assert ndims(X_clean) == 2
		@assert ndims(y_clean) == 3
		@assert size(X_clean, 1) == size(y_clean, 3)

		#UnfoldDecode.solver_b2b(X_clean, y_clean; cross_val_reps = cross_val_reps, show_progress = false)
		UnfoldDecode.solver_b2b(
			X_clean,
			y_clean;
			cross_val_reps = cross_val_reps,
			solver_G = UnfoldDecode.model_ridge,
			solver_H = UnfoldDecode.model_ridge,
			show_progress = false,
		)
	end

	uf_two_step_b2b = Unfold.fit(
		UnfoldModel,
		design_two_step_b2b,
		evts,
		X_corrected;
		solver = b2b_solver,
	)

	return (
		predictors = predictors,
		model_predictors = model_predictors,
		events = evts,
		rerp_model = uf_rerp,
		corrected_trials = X_corrected,
		b2b_fit_corrected = uf_two_step_b2b,
	)
end