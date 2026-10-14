function concatenate_runs(subject, all_events; eeg_dir)
    eeg_parts = Matrix{Float64}[]
    events_parts = DataFrame[]

    sample_offset = 0
    for run in 1:5
        # Load continuous EEG for this run 
        eeg_path = joinpath(eeg_dir, "$(subject)_run$(run)_eeg.npy")
        eeg = Float64.(NPZ.npzread(eeg_path)) .* 1e6

        # Selecyt events belonging to this run
        evts = filter(row -> row.subject_id == subject && row.run_num == run, all_events)

        # Check run-local event latencies (alrerady 1-based)
        @assert all((evts.latency .>= 1) .& (evts.latency .<= size(eeg, 2)))

        # Convert run-local latencies to subject-level latencies
        evts.latency .+= sample_offset 
        push!(eeg_parts, eeg)
        push!(events_parts, evts)
        sample_offset += size(eeg, 2)
    end

    # Concatenate after loading all five runs    
    eeg_combined = hcat(eeg_parts...)
    events_combined = vcat(events_parts...)
    sort!(events_combined, :latency)
    @assert size(eeg_combined, 2) == sample_offset
    
    return (eeg_combined, events_combined)
end