function load_roamm_decoding_scores(results_dir)
	files = filter(f -> endswith(f, ".jls"), readdir(results_dir; join = true))
	isempty(files) && error("No .jls files found in $results_dir")
	tables = DataFrame[]
	for file in files
		result = deserialize(file)
		df = DataFrame(result.scores)
		df[!, :subject] = fill(String(result.subject), nrow(df))
		push!(tables, df)
	end

	return vcat(tables...)
end


function load_b2b_results(results_dir, model_field::Symbol)
	files = filter(f -> endswith(f, ".jls"), readdir(results_dir; join = true))
	isempty(files) && error("No .jls files found in $results_dir")

	tables = DataFrame[]
	for file in files
		result = deserialize(file)
		df = DataFrame(Unfold.coeftable(getproperty(result, model_field)))
		df.subject = fill(String(result.subject), nrow(df))
		df.estimate_signed = copy(df.estimate)
		df.estimate = abs.(df.estimate_signed)
		filter!(:coefname => !=("(Intercept)"), df)
		push!(tables, df)
	end

	return vcat(tables...)
end

function winsorized_mean(x; proportion=0.1)
	v = sort(collect(x))
	k = floor(Int, proportion * length(v))
	k == 0 && return mean(v)
	v[1:k] .= v[k+1]
	v[end-k+1:end] .= v[end-k]
	return mean(v)
end

function summarize_decoding(scores)
	# One value per participant, predictor and time point
	subject_scores = combine(groupby(scores, [:subject, :target, :time]), :r => mean => :r)
	summary = combine(
		groupby(subject_scores, [:target, :time]),
		:r => mean => :r_mean,
		:r => median => :r_median,
		:r => std => :r_std,
		:r => length => :n,
		:r => (x -> winsorized_mean(x; proportion=0.1)) => :winsor10,
		:r => (x -> winsorized_mean(x; proportion=0.2)) => :winsor20	
	)
	summary.r_sem = summary.r_std ./ sqrt.(summary.n)
	return summary
end

function summarize_b2b(scores)
	summary = combine(
		groupby(scores, [:coefname, :time]),
		:estimate => mean => :estimate_mean,
		:estimate => median => :estimate_median,
		:estimate => std => :estimate_std,
		:estimate => length => :n,
		:estimate => (x -> winsorized_mean(x; proportion=0.1)) => :winsor10,
		:estimate => (x -> winsorized_mean(x; proportion=0.2)) => :winsor20	
	)
	summary.estimate_sem = summary.estimate_std ./ sqrt.(summary.n)
	return summary
end

function load_b2b_evaluation_scores(results_dir; metric=:R)
    files = filter(
        f -> endswith(f, ".jls"),
        readdir(results_dir; join=true)
    )

    isempty(files) && error("No evaluation files in $results_dir")

    targets = ["word_length", "frequency_zipf", "word_surprisal"]

    rows = DataFrame(
        subject = String[],
        target = String[],
        time = Float64[],
        r = Float64[]
    )

    for file in files
        data = deserialize(file)
        ev = data.evaluation

        values = getproperty(ev, metric)

        @assert size(values, 1) == length(ev.times)
        @assert size(values, 2) == length(targets) + 1

        for (i, target) in enumerate(targets)
            for t in eachindex(ev.times)
                push!(rows, (
                    string(data.subject),
                    target,
                    Float64(ev.times[t]),
                    Float64(values[t, i + 1])
                ))
            end
        end
    end

    return rows
end