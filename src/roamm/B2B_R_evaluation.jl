function _prepare_fold_predictors(X_raw, train_idx; scaling = :center)
	μ = mean(X_raw[train_idx, :], dims = 1)

	if scaling == :center
		X_scaled = X_raw .- μ
	else scaling == :zscore
		σ = std(X_raw[train_idx, :], dims = 1)
		X_scaled = (X_raw .- μ) ./ σ
	end
	return hcat(ones(size(X_raw, 1)), X_scaled)
end


function _make_internal_splits(train_idx; cross_val_reps = 5, seed = 12)
	splits = []

	for rep in 1: cross_val_reps
		rng = MersenneTwister(seed + rep)
		shuffled = shuffle(rng, collect(train_idx))
		mid = fld(length(shuffled), 2)

		idx_G = shuffled[1:mid]
		idx_H = shuffled[mid+1:end]

		push!(splits, (idx_G, idx_H))
	end 
	return splits
end

function _fit_b2b(X, Y, internal_splits)
	n_preds = size(X, 2)
	n_channels = size(Y, 2)
	n_reps = length(internal_splits)

	# ensembel coefficients
	G_mean = zeros(n_channels, n_preds)
	H_mean = zeros(n_preds, n_preds)
	# Knockout H for every predictor
	H_ko_mean = zeros(n_preds, n_preds, n_preds)

	# repeated b2b estimation
	for (idx_G, idx_H) in internal_splits
		G = UnfoldDecode.model_ridge(Y[idx_G, :], X[idx_G, :])
		Y2G = Y[idx_H, :] * G
		H = UnfoldDecode.model_ridge(X[idx_H, :], Y2G)

		G_mean .+= G
		H_mean .+= H

		# Knockout models
		for i in 2:n_preds
			X_H_ko = copy(X[idx_H, :])
			X_H_ko[:, i] .= 0

			H_ko = UnfoldDecod.model_ridge(X_H_ko, Y2G)
			H_ko_mean[:, :, i] .+= H_ko
		end
	end
	return (
		G = G_mean ./ n_reps,
		G = H_mean ./ n_reps,
		H_knockout = H_ko_mean ./ n_reps
	)
end

function _r_scores(X_test, Y_test, models)
	n_preds = size(X_test, 2)

	R = fill(NaN, n_preds)
	R_knockout = fill(NaN, n_preds)
	delta_R = fill(NaN, n_preds)

	YG = Y_test * models.G
	XH = X_test * models.H
	
	for i in 2:n_preds
		y_decoded = YG[:, i]
		x_reconstructed = XH[:, i]

		R[i] = cor(x_reconstructed, y_decoded)

		# evaluate knockout model
		X_ko = copy(X_test)
		X_ko[:, i] .= 0

		H_ko = models.H_knockout[:, :, i]
		x_ko = (X_ko * H_ko)[:, i]

		R_knockout[i] = cor(x_ko, y_decoded)
		delta_R[i] = R[i] - R_knockout[i]
	end
	return (
		R = R, 
		R_knockout = R_knockout,
		delta_R = delta_R
	)
end

function calculate_b2b_R(
	X_raw,
	dat,
	times;
	nfolds::Int = 5,
	cross_val_reps::Int = 5,
	scaling::Symbol = :center,
	seed = 12
)
	# X_raw: trials x predictors, no intercept
	# dat: channels x timepoints x trials

	n_trials, n_preds = size(X_raw)
	n_channels, n_times, n_epochs = size(dat)
	n_preds = n_predictors + 1

	@assert n_trials == n_epochs
	@assert length(times) == n_times
	@assert nfolds >= 2
	@assert cross_val_reps >= 1
	@assert all(isfinite, X_raw)
	
	# 1. Outer cross-validation
	
	folds = UnfoldDecode.MLJBase.train_test_pairs(
		UnfoldDecode.MLJBase.CV(
			nfolds = nfolds,
			shuffle = true,
			rng = seed
		),
		1:n_trials
	)
	
	n_outer = length(folds)

	# fold x time x predictor
	fold_R = fill(NaN, n_outer, n_times, n_preds)
	fold_delta_R = fill(NaN, n_outer, n_times, n_preds)
	fold_R_knockout = fill(NaN, n_outer, n_times, n_preds)
	
	# 2. Iterate over outer folds
	for (fold_idx, (train_idx, test_idx)) in enumerate(folds)

		println("Outer fold $fold_idx / $n_outer")

		X = _prepare_fold_predictors(X_raw, train_idx; scaling = scaling)

		@assert isempty(intersect(train_idx, test_idx))
		
		splits = _make_internal_splits(
			train_idx; 
			cross_val_reps = cross_val_reps,
			seed = seed + 1000 * fold_idx
		)
		
		for t in 1:n_times
			Y = Matrix(permutedims(dat[:, t, :]))

			# Fit models using only outer training trials
			models = _fit_b2b(X, Y, splits)

			# Evaluate on outer test trials
			scores = _r_scores(
				X[test_idx, :],
				Y[test_idx, :],
				models
			)

			fold_R[fold_idx, t, :] .= scores.R
			fold_R_knockout[fold_idx, t, :] .= scores.R_knockout
			fold_delta_R[fold_idx, t, :] .= scores.delta_R
		end
	end
		

	#Average acrosss outer CV folds
	function mean_finite_folds(scores)
		result = fill(NaN, n_times, n_preds)

		for t in 1:n_times, i in 2:n_preds
			vals = filter(isfinite, scores[:, t, i])
			if !isempty(vals)
				result[t, i] = mean(vals)
			end
		end
		
		return result
	end

	R = mean_finite_folds(fold_R)
	R_knockout = mean_finite_folds(fold_R_knockout)
	delta_R = mean_finite_folds(fold_delta_R)

	return(
		R = R,
		R_knockout = R_knockout,
		delta_R = delta_R,
		fold_R = fold_R,
		fold_R_knockout = fold_R_knockout,
		fold_delta_R = fold_delta_R,
		times = times
	)
end