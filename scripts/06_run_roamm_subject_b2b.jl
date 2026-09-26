#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MScB2B
using CSV
using DataFrames
using Unfold
using Random

Random.seed!(42)


const DATA_DIR =
	"/scratch/data/ROAMM/outputs/eeg"

const PREPROC_DIR =
	"/scratch/data/ROAMM/outputs/preprocessing"

const RESULTS_DIR =
	joinpath(
		@__DIR__,
		"..",
		"results",
		"roamm_b2b"
	)

const PREDICTORS = [:word_length]

const CROSS_VAL_REPS = 20

mkpath(RESULTS_DIR)


# Load event table
all_events = CSV.read(
	joinpath(
		PREPROC_DIR,
		"roamm_b2b_events.csv"
	),
	DataFrame
)

subjects = ["sub-10177"]
#subjects = sort(unique(all_events.subject_id))

println(
	"Number of subjects: ",
	length(subjects)
)


for subject in subjects

	println()
	println("==============================")
	println("Processing ", subject)
	println("==============================")


	output_file = joinpath(
		RESULTS_DIR,
		"$(subject)_word_length_b2b.csv"
	)


	# Allows restarting the batch
	if isfile(output_file)
		println("Already finished, skipping.")
		continue
	end


	try

		fit =
			MScB2B.fit_two_step_b2b_roamm_subject(
				subject,
				all_events;
				data_dir = DATA_DIR,
				predictors = PREDICTORS,
				sfreq = 256.0,
				cross_val_reps = CROSS_VAL_REPS,
			)


		result = DataFrame(
			Unfold.coeftable(
				fit.b2b_fit_corrected
			)
		)

        # remove unused columns containing `nothing`
        select!(result, Not([:group, :stderror]))


		# Store signed values
		result[!, :estimate_signed] =
			copy(result.estimate)


		# B2B magnitude
		result.estimate .=
			abs.(result.estimate)


		# Remove intercept
		result = result[
			result.coefname .!= "(Intercept)",
			:
		]


		result[!, :subject] =
			fill(
				subject,
				nrow(result)
			)

		result[!, :n_trials] =
			fill(
				fit.n_trials,
				nrow(result)
			)


		CSV.write(
			output_file,
			result
		)

		println(
			"Saved: ",
			output_file
		)


	catch err

        if isfile(output_file)
            rm(output_file; force=true)
        end

		println(
			"FAILED: ",
			subject
		)

		showerror(
			stdout,
			err,
			catch_backtrace()
		)

		println()

	end


	GC.gc()
end


println()
println("Finished batch.")