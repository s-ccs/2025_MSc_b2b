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
	
	t0 = time()
	uf_rerp = Unfold.fit(
		UnfoldLinearModelContinuousTime,
		design_rerp,
		evts,
		dat_cont;
		eventcolumn = :event,
		show_progress = false
	)
	println("rERP time: ", round(time() - t0; digits=2), " s")

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

	t1 = time()
	uf_two_step_b2b = Unfold.fit(
		UnfoldModel,
		design_two_step_b2b,
		evts,
		X_corrected;
		solver = b2b_solver,
	)
	println("B2B time: ", round(time() - t1; digits=2), " s")

	return (
		predictors = predictors,
		model_predictors = model_predictors,
		corrected_trials_shape = size(X_corrected),
		rerp_model = uf_rerp,
		b2b_fit_corrected = uf_two_step_b2b
	)
end