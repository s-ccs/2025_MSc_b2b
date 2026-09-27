function _fit_two_step_b2b_case(
	cfg,
	case_data,
	formula;
	cross_val_reps::Int 
) 
	
	dat_cont = case_data.continuous
	evts = case_data.events_continuous


	# Step 1: Fit rERP model to continuous EEG data
	design_rerp = [
		"stimulus" => (
			formula,
			Unfold.firbasis(
				τ = [-0.1, 1.0],
				sfreq = cfg.sfreq
			))]

	uf_rerp = Unfold.fit(
		UnfoldLinearModelContinuousTime,
		design_rerp,
		evts,
		dat_cont;
		eventcolumn = :event
	)

	# Step 2: singletrials()
	X_corrected = UnfoldDecode.singletrials(
		dat_cont,
		uf_rerp,
		evts,
		"stimulus",
		:event
	)

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

		UnfoldDecode.solver_b2b(X_clean, y_clean; cross_val_reps = cross_val_reps, show_progress = false)
	end

	uf_two_step_b2b = Unfold.fit(
		UnfoldModel,
		design_two_step_b2b,
		evts,
		X_corrected;
		solver = b2b_solver,
	)

	return (
		rerp_model = uf_rerp,
		corrected_trials = X_corrected,
		b2b_fit_corrected = uf_two_step_b2b,
	)
end 


function _run_two_step_b2b(
	cfg,
	simulation,
	cases,
	formula;
	cross_val_reps::Int
)
	score_tables = DataFrame[]
	fitted_models = Dict{Symbol, Any}()

	for case_name in cases
		case_data = getproperty(simulation, case_name)

		two_step_b2b_fit = _fit_two_step_b2b_case(
			cfg,
			case_data,
			formula;
			cross_val_reps = cross_val_reps,
		)

		result = DataFrame(Unfold.coeftable(two_step_b2b_fit.b2b_fit_corrected))

		# keep raw B2B estimate for debugging
		result[!, :estimate_signed] = copy(result.estimate)

		# B2B has no sign
		result.estimate .= abs.(result.estimate)

		# remove intercept
		result = result[result.coefname .!= "(Intercept)", :]

		result[!, :case] = fill(String(case_name), nrow(result))

		push!(score_tables, result)
		fitted_models[case_name] = two_step_b2b_fit
	end

	return (
		score_tables = vcat(score_tables...),
		fitted_models = fitted_models
	)
end


function run_two_step_b2b(
	cfg::ConditionContinuousConfig,
	simulation;
	cross_val_reps::Int = 3
)
	formula = @formula(0 ~ 1 + condition + continuous)

	return _run_two_step_b2b(
		cfg,
		simulation,
		(:clean, :overlap, :confound, :both),
		formula;
		cross_val_reps = cross_val_reps
	)
end

function run_two_step_b2b(
	cfg::CorrelatedContinuousConfig,
	simulation;
	cross_val_reps::Int = 3
)
	formula = @formula(0 ~ 1 + continuous1 + continuous2 + continuous3)

	return _run_two_step_b2b(
		cfg,
		simulation,
		(:no_overlap, :overlap),
		formula;
		cross_val_reps = cross_val_reps
	)
end