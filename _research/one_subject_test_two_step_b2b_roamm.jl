function _center_predictors!(evts, predictors)
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

function _make_formula(predictors)
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
	cross_val_reps::Int = 3
)
	evts = copy(evts)
	sort!(evts, :latency)

	# Unfold needs an event typr macthing "stimulus"
	evts[!, :event] = fill("stimulus", nrow(evts))

	model_predictors = _center_predictors!(evts, predictors)
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


function fit_two_step_b2b_roamm_subject(
	subject,
	all_events;
	data_dir,
	predictors::Vector{Symbol} = [:word_length],
	sfreq::Float64 = 256.0,
	cross_val_reps::Int = 20
)

	evts_subject = copy(
		all_events[all_events.subject_id .== subject, :]
	)

	nrow(evts_subject) > 0 ||
		error("No events found for $(subject).")

	sort!(evts_subject, [:run_num, :latency])

	# Center predictors ONCE across the whole subject
	model_predictors = _center_predictors!(
		evts_subject,
		predictors
	)

	formula = _make_formula(model_predictors)

	runs = sort(unique(evts_subject.run_num))

	println("Subject: ", subject)
	println("Runs: ", runs)
	println("Total fixations: ", nrow(evts_subject))
	println("Formula: ", formula)


	event_runs = DataFrame[]

	X_subject = nothing
	times_ref = nothing

	trial_offset = 0


	for run in runs

		run_int = Int(run)

		println()
		println("------------------------------")
		println("Run ", run_int)
		println("------------------------------")

		evts_run = copy(
			evts_subject[
				evts_subject.run_num .== run,
				:
			]
		)

		evts_run[!, :event] =
			fill("stimulus", nrow(evts_run))


		# Load EEG
		eeg_file = joinpath(
			data_dir,
			"$(subject)_run$(run_int)_eeg.npy"
		)

		isfile(eeg_file) ||
			error("EEG file not found: $(eeg_file)")

		dat_cont = NPZ.npzread(eeg_file)

		# volts -> μV
		dat_cont = dat_cont .* 1e6


		# Basic alignment check
		maximum(evts_run.latency) <= size(dat_cont, 2) ||
			error(
				"Event latency exceeds EEG length " *
				"for $(subject), run $(run_int)"
			)

		println("EEG size: ", size(dat_cont))
		println("Fixations: ", nrow(evts_run))


		design_rerp = [
			"stimulus" => (
				formula,
				Unfold.firbasis(
					τ = [-0.1, 1.0],
					sfreq = sfreq
				)
			)
		]

		uf_rerp = Unfold.fit(
			UnfoldLinearModelContinuousTime,
			design_rerp,
			evts_run,
			dat_cont;
			eventcolumn = :event
		)


		X_corrected = UnfoldDecode.singletrials(
			dat_cont,
			uf_rerp,
			evts_run,
			"stimulus",
			:event
		)

		println(
			"Corrected trials size: ",
			size(X_corrected)
		)

		@assert size(X_corrected, 3) == nrow(evts_run)


		times = collect(
			Unfold.times(uf_rerp)[1]
		)

		if isnothing(times_ref)
			times_ref = times

			# Allocate subject-level array once
			X_subject = similar(
				X_corrected,
				size(X_corrected, 1),
				size(X_corrected, 2),
				nrow(evts_subject)
			)

		else
			@assert times == times_ref
		end


		# Copy this run into subject-level array
		n_trials = size(X_corrected, 3)

		X_subject[
			:,
			:,
			trial_offset + 1 :
			trial_offset + n_trials
		] .= X_corrected

		trial_offset += n_trials

		push!(event_runs, evts_run)


		# Release run-level objects
		dat_cont = nothing
		X_corrected = nothing
		uf_rerp = nothing

		GC.gc()
	end



	evts_b2b = vcat(event_runs...)

	@assert trial_offset == nrow(evts_b2b)
	@assert size(X_subject, 3) == nrow(evts_b2b)
	@assert length(times_ref) == size(X_subject, 2)


	println()
	println("==============================")
	println("Subject-level B2B")
	println("==============================")
	println("EEG size: ", size(X_subject))
	println("Events: ", nrow(evts_b2b))




	design_two_step_b2b = [
		"stimulus" => (
			formula,
			times_ref
		)
	]


	b2b_solver = (X, y) -> begin

		X_clean, y_clean =
			Unfold.drop_missing_epochs(X, y)

		X_clean = Float64.(X_clean)

		@assert ndims(X_clean) == 2
		@assert ndims(y_clean) == 3
		@assert size(X_clean, 1) == size(y_clean, 3)

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
		evts_b2b,
		X_subject;
		solver = b2b_solver,
	)


	return (
		subject = subject,
		runs = runs,
		predictors = predictors,
		model_predictors = model_predictors,
		n_trials = nrow(evts_b2b),
		b2b_fit_corrected = uf_two_step_b2b,
	)
end