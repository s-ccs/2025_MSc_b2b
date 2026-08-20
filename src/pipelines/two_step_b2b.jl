function fit_two_step_b2b_case(
	cfg,
	case_data;
	cross_val_reps::Int = 3,
) 
	
	dat_cont = case_data.continuous
	evts = case_data.events_continuous


	# Step 1: Fit rERP model to continuous EEG data
	design_rerp = [
		"stimulus" => (
			@formula(0 ~ 1 + condition + continuous),
			Unfold.firbasis(
				τ = [-0.1, 1.0],
				sfreq = cfg.sfreq,
				name = "stimulus",
			))]

	uf_rerp = Unfold.fit(
		UnfoldLinearModelContinuousTime,
		design_rerp,
		evts,
		dat_cont;
		eventcolumn = :event,
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
	
	design_b2b = [
		"stimulus" => (
			@formula(0 ~ 1 + condition + continuous),
			times,
		)
	]

	b2b_solver = (X, y) -> UnfoldDecode.solver_b2b(X, y; cross_val_reps = cross_val_reps)

	uf_b2b = Unfold.fit(
		UnfoldModel,
		design_b2b,
		evts,
		X_corrected;
		solver = b2b_solver,
	)

	return (
		rerp_model = uf_rerp,
		corrected_trials = X_corrected,
		b2b_fit_correted = uf_b2b,
	)

function run_two_step_b2b(
	cfg,
	cases;
	cross_val_reps::Int = 3,
)
	score_tables = DataFrame[]
	fitted_models = Dict{Symbol, Any}()

	for case_name in (:clean, :overlap, :confound, :both,)
		case_data = getproperty(cases, case_name)

		two_step_b2b_fit = fit_two_step_b2b_case(
			cfg,
			case_data;
			cross_val_reps = cross_val_reps,
		)

		result = DataFrame(Unfold.coeftable(two_step_b2b_fit.b2b_fit_correted))

		# keep raw B2B estimate for debugging
		result[!, :estimate_singed] = copy(result.estimate)

		# B2B has no sign
		result.estimate .= abs.(result.estimate)

		# remove intercept
		result = result[result.coefname .!= "(intercept)", :]

		result[!, :case] = fill(String(case_name), nrow(result))

		push!(score_tables, result)
		fitted_models[case_name] = two_step_b2b_fit
	end

	return (
		score_tables = vcat(score_tables...),
		fitted_models = fitted_models,
	)
end