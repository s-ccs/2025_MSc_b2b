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
	n_times = size(dat, 2)

	folds = UnfoldDecode.MLJBase.train_test_pairs(
		UnfoldDecode.MLJBase.CV(
			nfolds = nfolds,
			shuffle = true,
			rng = seed
		),
		1:n_trials
	)

	# fold x time x predictor
	fold_R = fill(NaN, length(folds), n_times, n_preds)

	for (fold_idx, (train_idx, test_idx)) in enumerate(folds)

		println("Outer fold $fold_idx / $(length(folds))")
		rep_R = fill(NaN, cross_val_reps, n_times, n_preds)

		for rep in 1:cross_val_reps
			# Internal B2B split: training trials only
			rng = MersenneTwister(seed + 1000 * fold_idx + rep)
			shuffled = shuffle(rng, collect(train_idx))
			mid = fld(length(shuffled, 2))

			idx_G = shuffled[1:mid]
			idx_H = shuffled[mid+1:end]

			for t in 1:n_times
				
				# EEG: trials x channels
				Y = Matrix(permutedims(dat[:, t, :]))

				# Stage 1: EEG -> predictors
				G = UnfoldDecode.model_ridge(
					Y[idx_G, :],
					X[idx_G, :]
				)

				# Decode the H-training trials
				Y2G = Y[idx_H, :] * G

				# Stage 2: predictors -> decoded predictors
				H = UnfoldDecode.model_ridge(
					X[idx_H, :],
					Y2G
				)

				# Evaluate on independent outer test trials
				XH_test = X[test_idx, :] * H
				YG_test = Y[test_idx, :] * G

				# Skip the intercept: column 1
				for i in 2:n_preds
					a = XH_test[:, i]
					b = YG_test[:, i]
					if std(a) > 1e-12 && std(b) > 1e-12
						rep_R[rep, t, i] = cor(a, b)
					end
				end
			end
		end

		# Mean across outer folds
		R = fill(NaN, n_times, n_preds)

		for t in 1:n_times, i in 2:n_preds
			vals = filter(isfinite, fold_R[:, t, i])
			if !isempty(vals)
				fold_R[fold_idx, t, i] = mean(vals)
			end
		end
	end

	# Mean across outer folds
	R = fill(NaN, n_times, n_preds)
	for t in 1:n_times, i in 2:n_preds
		vals = filter(isfinite, fold_R[:, t, i])
		if !isempty(vals)
			R[t, i] = mean(vals)
		end
	end
	return(
		R = R,
		fold_R = fold_R,
		times = times
	)
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
