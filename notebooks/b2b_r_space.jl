### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 604b38cc-d4b4-11f1-813e-6b8ca59cde12
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 80655def-a280-4ab4-896d-b44a1ebe778f
begin
	using Serialization, Unfold, DataFrames, UnfoldDecode
	using CairoMakie, Statistics, CSV, MScB2B, Random
end

# ╔═╡ e74535a3-82e9-4421-a7a8-fb874ff558c6
begin
	#base_dir = "/home/st/st_us-051900/st_st181570/2025_MSc_b2b/slurm_results/models/roamm"
	base_dir = "/store/users/xu/2025_MSc_b2b/slurm_results/models/roamm"
	
	readdir(base_dir)
end


# ╔═╡ 62fb99dd-d987-4645-8850-6c49f0d9367d
begin
	file = joinpath(base_dir, "two_step_roamm", "sub-10014_center_cv5.jls")

	result = deserialize(file)
	typeof(result)
end

# ╔═╡ 5b9a5b06-b8ae-42d5-b3ef-66dadd0629ff
propertynames(result)

# ╔═╡ 3e6f2b5c-2d34-422d-802d-7ace1dca4333
typeof(result.b2b_fit_corrected)

# ╔═╡ 4f7c9cb7-b5e4-4781-8228-9bdd5637d7cb
propertynames(result.b2b_fit_corrected)

# ╔═╡ 9bf14bb6-114a-4497-a224-13fc267636d2


# ╔═╡ 3fde22d9-0ed2-41c7-841d-57a15dede348
begin
	m = result.b2b_fit_corrected

typeof(m.modelfit)
propertynames(m.modelfit)
end

# ╔═╡ e6a6568b-252a-42d5-826f-c2b8003c1465
# ╠═╡ disabled = true
#=╠═╡
begin
	df = DataFrame(Unfold.coeftable(m))
first(df, 10)
end
  ╠═╡ =#

# ╔═╡ d4fc567d-eaf1-4fd4-a780-6975c62e41f8
begin
	df = DataFrame(Unfold.coeftable(m))
	
	filter(
	    :coefname => !=("(Intercept)"),
	    df
	)
end

# ╔═╡ b58604be-d5d2-474d-8c95-3f048227c763
function calculate_b2b_R(
	X,
	dat,
	times;
	nfolds=5,
	cross_val_reps = 5,
	seed = 12
)
	# X: trials x predictors, including intercept
	# dat: channels x timepoints x trials

	n_trials, n_preds = size(X)
	n_channels, n_times, n_epochs = size(dat)

	@assert n_trials == n_epochs
	@assert length(times) == n_times
	@assert nfolds >= 2
	@assert cross_val_reps >= 1
	
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

		@assert isempty(intersect(train_idx, test_idx))

		# Generate internal splits once per outer fold
		internal_splits = []

		for rep in 1:cross_val_reps
			rng = MersenneTwister(seed + 1000 * fold_idx + rep)

			shuffled = shuffle(rng, collect(train_idx))
			mid = fld(length(shuffled), 2)

			idx_G = shuffled[1:mid]
			idx_H = shuffled[mid+1:end]

			@assert !isempty(idx_G) && !isempty(idx_H)
			@assert isempty(intersect(idx_G, idx_H))

			push!(internal_splits, (idx_G, idx_H))
		end



		
		# Iterate over EEG time points
		for t in 1:n_times

			# trials x channels
			# : select everything along this dimension
			Y = Matrix(permutedims(dat[:, t, :]))

			# Ensemble coefficients
			G_mean = zeros(n_channels, n_preds)
			H_mean = zeros(n_preds, n_preds)

			# Knockout H for every predictor
			H_ko_mean = zeros(n_preds, n_preds, n_preds)

			# 4. Repeated B2B estimation
			for (idx_G, idx_H) in internal_splits
	
				# Stage 1: EEG -> predictors
				G = UnfoldDecode.model_ridge(
					Y[idx_G, :],
					X[idx_G, :]
				)

				# Decode the H-training trials
				Y2G = Y[idx_H, :] * G

				# Stage 2: full predictor set -> Y2G
				H = UnfoldDecode.model_ridge(
					X[idx_H, :],
					Y2G
				)

				# Evaluate on independent outer test trials
				G_mean .+= G
				H_mean .+= H

				# 5. Fit knockout H models using the same split
				for i in 2:n_preds
					X_H_ko = copy(X[idx_H, :])
					X_H_ko[:, i] .= 0

					H_ko = UnfoldDecode.model_ridge(
						X_H_ko, 
						Y2G
					)
					H_ko_mean[:, :, i] .+= H_ko
				end
			end
				
			# 6. Average coefficient across B2B repetitions
			G_mean ./= cross_val_reps
			H_mean ./= cross_val_reps
			H_ko_mean ./= cross_val_reps

			# 7. Evaluate full B2B on outer test set
			YG_test = Y[test_idx, :] * G_mean
			XH_test = X[test_idx, :] * H_mean

			for i in 2:n_preds
				y_decoded = YG_test[:, i]
				x_reconstructed = XH_test[:, i]

				if std(y_decoded) <= 1e-12 || 
					std(x_reconstructed) <= 1e-12
					continue
				end
				R_full = cor(x_reconstructed, y_decoded)

				fold_R[fold_idx, t, i] = R_full

				# 8. Evaluate predictor-i knockout
				X_test_ko = copy(X[test_idx, :])
				X_test_ko[:, i] .= 0

				X_ko = H_ko_mean[:, :, i]
				XH_ko = X_test_ko * H_ko

				if std(XH_ko[:, i]) <= 1e-12
					continue 
				end

				R_ko = cor(XH_ko[:, i], y_decoded)

				fold_R_knockout[fold_idx, t, i] = R_ko
				fold_delta_R[fold_idx, t, i] = R_full - R_ko
			end
		end
	end

	# 9. Average acrosss outer CV folds
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
	delta_R = mean_finite_folds(fold_delt_R)

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

# ╔═╡ 81b670aa-be58-4d35-bdef-4f3bb7c19e35
begin 
	case_data = sim_cases.no_overlap

	dat = case_data.epoched
	times = case_data.times
	evts = case_data. events_epoched

	X = hcat(
		ones(nrow(evts)),
		Matrix(select(
			evt,
			[:continuous1, :continuous2, :continuous3]
		))
	)

	Y = Matrix(permiutedims(dat[:, 1, :]))

	G_test = UnfoldDecode.model_ridge(
		Y[1:100, :],
		X[1:100, :]
	)

	@show typeof(G_test)
	@show size(G_test)

	@assert size(G_test) == (size(Y, 2), size(X, 2))
end
	

# ╔═╡ 3ffe9d32-25b8-4cc0-828f-5fcf2048d209
begin
	result_R = calculate_b2b_R(
		X,
		dat,
		times;
		nfolds = 5, 
		cross_val_reps = 5,
		seed = 12
	)

	@show size(result_R.R)
	@show size(result_R.delta_R)
end

# ╔═╡ f47d48d1-c6b6-4c49-aa31-fac5279d3968

begin
    fig = Figure(size = (950, 700))

    names = ["c1", "c2", "c3"]
    colors = [:steelblue, :darkorange, :seagreen]

    ax1 = Axis(
        fig[1, 1],
        xlabel = "Time (s)",
        ylabel = "Correlation (r)",
        title = "B2B — Held-out reconstruction correlation"
    )

    ax2 = Axis(
        fig[2, 1],
        xlabel = "Time (s)",
        ylabel = "ΔR",
        title = "B2B — Predictor knockout contribution"
    )

    for i in 1:3
        lines!(
            ax1,
            result_R.times,
            result_R.R[:, i+1],
            color = colors[i],
            label = names[i]
        )

        lines!(
            ax2,
            result_R.times,
            result_R.delta_R[:, i+1],
            color = colors[i],
            label = names[i]
        )
    end

    for ax in (ax1, ax2)
        vlines!(ax, [0.0], color = :gray, linestyle = :dash)
        hlines!(ax, [0.0], color = :gray, linestyle = :dash)
    end

    axislegend(ax1)
    fig
end


# ╔═╡ Cell order:
# ╠═604b38cc-d4b4-11f1-813e-6b8ca59cde12
# ╠═80655def-a280-4ab4-896d-b44a1ebe778f
# ╠═e74535a3-82e9-4421-a7a8-fb874ff558c6
# ╠═62fb99dd-d987-4645-8850-6c49f0d9367d
# ╠═5b9a5b06-b8ae-42d5-b3ef-66dadd0629ff
# ╠═3e6f2b5c-2d34-422d-802d-7ace1dca4333
# ╠═4f7c9cb7-b5e4-4781-8228-9bdd5637d7cb
# ╠═9bf14bb6-114a-4497-a224-13fc267636d2
# ╠═3fde22d9-0ed2-41c7-841d-57a15dede348
# ╠═e6a6568b-252a-42d5-826f-c2b8003c1465
# ╠═d4fc567d-eaf1-4fd4-a780-6975c62e41f8
# ╠═b58604be-d5d2-474d-8c95-3f048227c763
# ╠═81b670aa-be58-4d35-bdef-4f3bb7c19e35
# ╠═3ffe9d32-25b8-4cc0-828f-5fcf2048d209
# ╠═f47d48d1-c6b6-4c49-aa31-fac5279d3968
