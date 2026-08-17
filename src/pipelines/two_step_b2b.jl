function fit_two_step_b2b_cases(
    cases,
    des_rerp;
    eventname = "stimulus",
    eventcolumn = :event,
)
    result_tables = DataFrame[]

    for (case_name, (dat_cont, evts)) in pairs(cases)

		# ==========================================
		# Step 1: rERP
		# ==========================================
		
		uf_rerp = Unfold.fit(
		    UnfoldLinearModelContinuousTime,
		    des_rerp,
		    evts,
		    dat_cont;
		    eventcolumn = eventcolumn,
		)
		
		
		# ==========================================
		# Step 2: corrected single trials
		# ==========================================
		
		dat_corrected = UnfoldDecode.singletrials(
		    dat_cont,
		    uf_rerp,
		    evts,
		    eventname,
		    eventcolumn,
		)
		
		
		# ==========================================
		# Step 2.5: remove incomplete boundary epochs
		# ==========================================
		
		has_missing = dropdims(
		    any(
		        ismissing.(dat_corrected),
		        dims = (1, 2),
		    ),
		    dims = (1, 2),
		)
		
		good_trials = .!has_missing
		
		dat_corrected =
		    Float64.(
		        dat_corrected[:, :, good_trials]
		    )
		
		evts_corrected =
		    evts[good_trials, :]
		
		@assert size(dat_corrected, 3) == nrow(evts_corrected)
		@assert count(ismissing, dat_corrected) == 0
		
		
		# ==========================================
		# Step 3: epoched B2B design
		# ==========================================
				times_corrected =
		    Unfold.times(uf_rerp)[1]
		
		@assert length(times_corrected) ==
		        size(dat_corrected, 2)
		
		des_b2b_corrected = [
		    Any => (
		        @formula(0 ~ 1 + condition + continuous),
		        times_corrected,
		    )
		]
		
		
		# ==========================================
		# Step 4: B2B
		# ==========================================
		
		uf_b2b = Unfold.fit(
    UnfoldDecode.UnfoldModel,
    des_b2b_corrected,
    evts_corrected,
    dat_corrected;
    solver = b2b_solver,
)

        # ------------------------------------------
        # 6. collect coefficients
        # ------------------------------------------

        result = coeftable(uf_b2b)

        result = result[
            result.coefname .!= "(Intercept)",
            :
        ]

        result[!, :case] =
            fill(String(case_name), nrow(result))

        push!(result_tables, result)
    end

    return vcat(result_tables...)
end