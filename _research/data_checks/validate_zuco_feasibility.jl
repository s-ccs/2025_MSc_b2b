### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ddf1f343-818e-460b-b2c2-6a1f60f81311

begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 24be1c80-ee34-4afa-ae7d-a289971171a7
begin
	
	using MAT
	using DataFrames
	using Statistics
end

# ╔═╡ 7be98ca1-f507-4499-b196-ffc7c64d9ea8
begin
	using CairoMakie
	
	fig = Figure()
	
	ax = Axis(
	    fig[1, 1],
	    xlabel = "Fixation",
	    ylabel = "Global - local sample offset",
	)
	
	scatter!(ax, 1:length(offsets), offsets)
	
	fig
end

# ╔═╡ d56938ec-8fdb-4f06-8d6f-9414b053f37d

md"""
# ZuCo 2.0 feasibility gate

This notebook gives a **go / conditional-go / stop** decision for using ZuCo in the continuous-EEG thesis pipeline.

It tests the complete chain:

1. raw ET ↔ raw EEG synchronization;
2. processed fixation → word mapping via `fixPositions`;
3. processed fixation ↔ raw ET fixation matching;
4. optional recovery of sentence-local onset from processed `rawEEG`;
5. construction of a fixation-level table containing word identity and EEG onset.

Recommended interpretation:

- **GO (processed EEG route)**: all major checks pass, including processed segment onset recovery.
- **GO (raw EEG route)**: sync + raw fixation matching + word mapping pass, but processed segment onset recovery is weak.
- **STOP / investigate**: the fixation linkage itself is unreliable.
"""

# ╔═╡ 0c2a70df-8c3d-4a05-a410-8a4acda94f0a

begin
    const SUBJECT = "YRP"
    const PROCESSED_PATH = "/scratch/data/ZuCo2/resultsYRP_NR.mat"
    const ET_PATH = "/scratch/data/ZuCo2/raw/YRP_NR1_ET.mat"
    const EEG_PATH = "/scratch/data/ZuCo2/raw/YRP_NR1_EEG.mat"
    const SFREQ = 500.0
end

# ╔═╡ dd188df4-d0ed-4031-97d4-d61168f76c9e

begin
    sentenceData = matopen(PROCESSED_PATH) do f
        read(f, "sentenceData")
    end

    et = matread(ET_PATH)
    eeg_file = matread(EEG_PATH)
    E = eeg_file["EEG"]
end

# ╔═╡ 823a827b-db50-4bab-b87f-fdf44b2ffab8

md"""
## Shared helpers
"""

# ╔═╡ 4d2d0861-bf08-48bc-b13b-43f3fc22f87b

begin
    normalize_word(w) = lowercase(replace(string(w), r"[^\p{L}']" => ""))

    function get_fixation_ids(x)
        if x === missing || x === nothing
            return Int[]
        elseif x isa Number
            return isnan(x) ? Int[] : [Int(round(x))]
        elseif x isa AbstractArray
            ids = Int[]
            for v in x
                append!(ids, get_fixation_ids(v))
            end
            return ids
        else
            return Int[]
        end
    end

    function invert_fixpositions(words, nfix::Int)
        mapping = Vector{Union{Missing, Int}}(missing, nfix)
        duplicates = Int[]
        invalid = Int[]

        for (word_id, pos) in enumerate(vec(words["fixPositions"]))
            for fix_id in get_fixation_ids(pos)
                if !(1 <= fix_id <= nfix)
                    push!(invalid, fix_id)
                elseif !ismissing(mapping[fix_id])
                    push!(duplicates, fix_id)
                else
                    mapping[fix_id] = word_id
                end
            end
        end

        return (
            mapping = mapping,
            duplicates = unique(duplicates),
            invalid = unique(invalid),
        )
    end

    function find_segment_start(eeg, seg)
        size(eeg, 1) == size(seg, 1) || return missing
        n = size(seg, 2)
        T = size(eeg, 2)
        n <= T || return missing

        anchor = findfirst(isfinite, seg)
        anchor === nothing && return missing
        ch, k = Tuple(anchor)
        value = seg[anchor]

        matches = Int[]
        for start in 1:(T - n + 1)
            isequal(eeg[ch, start + k - 1], value) || continue
            candidate = @view eeg[:, start:start+n-1]
            isequal(candidate, seg) && push!(matches, start)
        end

        isempty(matches) && return missing
        length(matches) == 1 && return only(matches)
        return missing
    end
end

# ╔═╡ 64f59fcf-4d9e-465d-b91d-bcbbf4ab8dfe

md"""
## 1. Raw EEG–ET synchronization
"""

# ╔═╡ 1664b49b-aea4-4555-b8c1-c63a5dfe3386
begin
	begin
		
		function synchronize_et_eeg(et, E)
		    et_ev = et["event"]
		    eeg_ev = E["event"]
		
		    et_types = round.(Int, et_ev[:, 2])
		    eeg_types = parse.(Int, strip.(string.(vec(eeg_ev["type"]))))
		    eeg_lat = Float64.(vec(eeg_ev["latency"]))
		
		    shared = intersect(unique(et_types), unique(eeg_types))
		    et_keep = findall(in(shared), et_types)
		    eeg_keep = findall(in(shared), eeg_types)
		
		    et_seq = et_types[et_keep]
		    eeg_seq = eeg_types[eeg_keep]
		
		    length(et_seq) == length(eeg_seq) || error("ET/EEG shared-event counts differ.")
		    et_seq == eeg_seq || error("ET/EEG trigger sequences differ.")
		
		    et_time = Float64.(et_ev[et_keep, 1])
		    eeg_sample = eeg_lat[eeg_keep]
		
		    X = hcat(ones(length(et_time)), et_time)
		    β = X \ eeg_sample
		
		    pred = X * β
		    sfreq = Float64(E["srate"])
		    err_ms = (eeg_sample .- pred) ./ sfreq .* 1000
		
		    return (
		        shared = shared,
		        et_seq = et_seq,
		        eeg_seq = eeg_seq,
		        β = β,
		        err_ms = err_ms,
		    )
		end
		
		sync = synchronize_et_eeg(et, E)
	end
end

# ╔═╡ dd41afae-9bb9-4d3a-b8e9-0ce409ba171c
begin
	
	sync_pass =
	    sync.et_seq == sync.eeg_seq &&
	    abs(sync.β[2] - Float64(E["srate"]) / 1000) <= 0.001 &&
	    mean(abs.(sync.err_ms)) <= 2.0
	
	sync_summary = (
	    pass = sync_pass,
	    n_triggers = length(sync.et_seq),
	    slope = sync.β[2],
	    expected_slope = Float64(E["srate"]) / 1000,
	    mean_abs_error_ms = mean(abs.(sync.err_ms)),
	    median_abs_error_ms = median(abs.(sync.err_ms)),
	    max_abs_error_ms = maximum(abs.(sync.err_ms)),
	)
end

# ╔═╡ fea44612-0ec1-49e8-9452-ceb1d6e70202

md"""
## 2. Match processed `allFixations` back to raw ET

Important observation to test explicitly:

`sentenceData["allFixations"]["duration"]` behaves like **500-Hz samples**, while raw ET uses milliseconds.
Therefore processed duration is converted as:

`duration_ms = duration_samples × 1000 / 500`.

The matcher uses duration + x-position + pupil size. `y` is deliberately not required because the processed structure may use a transformed/text-line y coordinate.
"""

# ╔═╡ 9d099b5f-f52a-443f-b2e6-96e478df1221
begin
	begin
		
		function processed_fixations(sentenceData, sentence_id::Int; sfreq::Real = SFREQ)
		    p = sentenceData["allFixations"][sentence_id]
		
		    dur_samples = Float64.(vec(p["duration"]))
		
		    DataFrame(
		        fixation_id = collect(1:length(dur_samples)),
		        duration_samples = dur_samples,
		        duration_ms = dur_samples .* (1000 / sfreq),
		        x = Float64.(vec(p["x"])),
		        y = Float64.(vec(p["y"])),
		        pupil = Float64.(vec(p["pupilsize"])),
		    )
		end
		
		function raw_fixations(et)
		    f = et["eyeevent"]["fixations"]["data"]
		
		    DataFrame(
		        raw_fixation_id = collect(1:size(f, 1)),
		        et_latency = Float64.(f[:, 1]),
		        et_endtime = Float64.(f[:, 2]),
		        duration_ms = Float64.(f[:, 2] .- f[:, 1]),
		        x = Float64.(f[:, 4]),
		        y = Float64.(f[:, 5]),
		        pupil = Float64.(f[:, 6]),
		    )
		end
		
		rawfix = raw_fixations(et)
	end
end

# ╔═╡ 26718990-7719-4df6-8e1f-31ebdaf92518
begin
	
	function same_fixation(
	    p,
	    r;
	    dur_tol_ms::Real = 4.0,
	    x_tol_px::Real = 0.6,
	    pupil_tol::Real = 2.0,
	)
	    abs(p.duration_ms - r.duration_ms) <= dur_tol_ms &&
	    abs(p.x - r.x) <= x_tol_px &&
	    abs(p.pupil - r.pupil) <= pupil_tol
	end
	
	function greedy_match_from(
	    proc::DataFrame,
	    raw::DataFrame,
	    start_idx::Int;
	    max_gap::Int = 6,
	)
	    matched = Vector{Union{Missing, Int}}(missing, nrow(proc))
	    raw_i = start_idx
	
	    for proc_i in 1:nrow(proc)
	        found = nothing
	
	        for ri in raw_i:min(nrow(raw), raw_i + max_gap)
	            if same_fixation(proc[proc_i, :], raw[ri, :])
	                found = ri
	                break
	            end
	        end
	
	        isnothing(found) && break
	
	        matched[proc_i] = found
	        raw_i = found + 1
	    end
	
	    return matched
	end
	
	function match_processed_sentence_to_raw(
	    sentenceData,
	    sentence_id::Int,
	    rawfix::DataFrame;
	    sfreq::Real = SFREQ,
	)
	    proc = processed_fixations(sentenceData, sentence_id; sfreq = sfreq)
	
	    starts = findall(
	        i -> same_fixation(proc[1, :], rawfix[i, :]),
	        1:nrow(rawfix),
	    )
	
	    isempty(starts) && return (
	        proc = proc,
	        matched = Vector{Union{Missing, Int}}(missing, nrow(proc)),
	        match_rate = 0.0,
	        raw_start = missing,
	        raw_end = missing,
	    )
	
	    candidates = [greedy_match_from(proc, rawfix, s) for s in starts]
	    counts = [count(!ismissing, c) for c in candidates]
	    best = candidates[argmax(counts)]
	
	    ids = collect(skipmissing(best))
	
	    return (
	        proc = proc,
	        matched = best,
	        match_rate = count(!ismissing, best) / length(best),
	        raw_start = isempty(ids) ? missing : first(ids),
	        raw_end = isempty(ids) ? missing : last(ids),
	    )
	end
end

# ╔═╡ 8d0b61b3-ad73-4a45-bd3d-2bf258f32911

md"""
## 3. Test the first few processed sentences against raw NR1
"""

# ╔═╡ e0a1d9c3-e96c-4a73-bfed-2f04bc275752
begin
	
	test_sentences = collect(1:min(5, length(sentenceData["allFixations"])))
	
	raw_matches = [
	    match_processed_sentence_to_raw(sentenceData, sid, rawfix)
	    for sid in test_sentences
	]
	
	raw_match_summary = DataFrame([
	    (
	        sentence_id = sid,
	        n_processed_fixations = nrow(m.proc),
	        match_rate = m.match_rate,
	        raw_start = m.raw_start,
	        raw_end = m.raw_end,
	    )
	    for (sid, m) in zip(test_sentences, raw_matches)
	])
end

# ╔═╡ 7d49ae2e-d5c3-4d26-ad84-4e48d5a96f7b

md"""
## 4. Validate fixation → word mapping and processed sentence-local onset recovery
"""

# ╔═╡ 162bd68c-28c3-4604-ac32-a51eed992b87
begin
	
	function processed_sentence_checks(sentenceData, sentence_id::Int; sfreq::Real = SFREQ)
	    words = sentenceData["word"][sentence_id]
	    allfix = sentenceData["allFixations"][sentence_id]
	    eeg = sentenceData["rawData"][sentence_id]
	
	    nfix = length(vec(allfix["duration"]))
	    inv = invert_fixpositions(words, nfix)
	
	    assigned_ids = findall(!ismissing, inv.mapping)
	
	    recovered = 0
	    checked = 0
	
	    raw_eeg_words = vec(words["rawEEG"])
	    fix_positions = vec(words["fixPositions"])
	
	    for word_id in eachindex(fix_positions)
	        ids = get_fixation_ids(fix_positions[word_id])
	        isempty(ids) && continue
	
	        segments = vec(raw_eeg_words[word_id])
	        length(ids) == length(segments) || continue
	
	        for seg in segments
	            checked += 1
	            !ismissing(find_segment_start(eeg, seg)) && (recovered += 1)
	        end
	    end
	
	    return (
	        sentence_id = sentence_id,
	        n_all = nfix,
	        n_assigned = length(assigned_ids),
	        n_unassigned = nfix - length(assigned_ids),
	        n_duplicates = length(inv.duplicates),
	        n_invalid = length(inv.invalid),
	        segment_checked = checked,
	        segment_recovered = recovered,
	        segment_recovery_rate = checked == 0 ? NaN : recovered / checked,
	    )
	end
	
	processed_checks = [
	    processed_sentence_checks(sentenceData, sid)
	    for sid in test_sentences
	]
	
	processed_summary = DataFrame(processed_checks)
end

# ╔═╡ 3efff367-0bbc-4132-b9a8-2987a4eca62a

md"""
## 5. Construct one final fixation-level event table

This joins:

`processed fixation ID → word`  
`processed fixation → raw ET fixation → ET onset`  
`ET onset → raw EEG sample`
"""

# ╔═╡ bc068290-64ad-41d7-8379-191ea49cf1d0
begin
	begin
		
		function build_final_events_for_sentence(
		    sentenceData,
		    sentence_id::Int,
		    rawfix::DataFrame,
		    sync;
		    subject::String = SUBJECT,
		)
		    m = match_processed_sentence_to_raw(sentenceData, sentence_id, rawfix)
		    proc = m.proc
		
		    words = sentenceData["word"][sentence_id]
		    nfix = nrow(proc)
		    inv = invert_fixpositions(words, nfix)
		    word_content = vec(words["content"])
		
		    rows = NamedTuple[]
		
		    for fix_id in 1:nfix
		        raw_id = m.matched[fix_id]
		        word_id = inv.mapping[fix_id]
		
		        (ismissing(raw_id) || ismissing(word_id)) && continue
		
		        et_onset = rawfix[raw_id, :et_latency]
		        eeg_sample = round(Int, sync.β[1] + sync.β[2] * et_onset)
		        word = string(word_content[word_id])
		
		        push!(
		            rows,
		            (
		                subject = subject,
		                sentence_id = sentence_id,
		                fixation_id = fix_id,
		                raw_fixation_id = raw_id,
		                word_id = word_id,
		                word = word,
		                word_lookup = normalize_word(word),
		                et_onset = et_onset,
		                eeg_sample_raw = eeg_sample,
		                fixation_duration_ms = rawfix[raw_id, :duration_ms],
		            ),
		        )
		    end
		
		    events = DataFrame(rows)
		    nrow(events) > 0 && sort!(events, :eeg_sample_raw)
		    return events
		end
		
		final_events_sentence1 =
		    build_final_events_for_sentence(sentenceData, 1, rawfix, sync)
	end
end

# ╔═╡ 31c2cdc4-7a8f-4f2b-9ca6-758323a84f51

md"""
## 6. Final feasibility decision
"""

# ╔═╡ 580ec354-4e97-4f30-9315-bd83255b9c1f
begin
	
	word_mapping_pass =
	    all(processed_summary.n_duplicates .== 0) &&
	    all(processed_summary.n_invalid .== 0)
	
	raw_fixation_link_pass =
	    all(raw_match_summary.match_rate .>= 0.95)
	
	processed_segment_pass =
	    all(processed_summary.segment_recovery_rate .>= 0.95)
	
	all_raw_samples_inside =
	    nrow(final_events_sentence1) > 0 &&
	    all(
	        (final_events_sentence1.eeg_sample_raw .>= 1) .&
	        (final_events_sentence1.eeg_sample_raw .<= size(E["data"], 2))
	    )
	
	verdict =
	    if sync_pass && word_mapping_pass && raw_fixation_link_pass && processed_segment_pass
	        "GO — processed sentence-continuous EEG route is feasible."
	    elseif sync_pass && word_mapping_pass && raw_fixation_link_pass && all_raw_samples_inside
	        "GO, WITH CONDITION — fixation/word/onset linkage works, but use/preprocess raw block EEG if processed sentence-local onset recovery remains weak."
	    else
	        "STOP / INVESTIGATE — ZuCo linkage is not yet reliable enough; consider ROAMM if this cannot be fixed quickly."
	    end
	
	(
	    verdict = verdict,
	    sync_pass = sync_pass,
	    word_mapping_pass = word_mapping_pass,
	    raw_fixation_link_pass = raw_fixation_link_pass,
	    processed_segment_pass = processed_segment_pass,
	    all_raw_samples_inside = all_raw_samples_inside,
	    raw_match_summary = raw_match_summary,
	    processed_summary = processed_summary,
	)
end

# ╔═╡ db706510-37f0-4f82-a1ee-471f31b2a5b8
function local_events_for_sentence(sentenceData, sentence_id::Int)

    eeg = sentenceData["rawData"][sentence_id]
    words = sentenceData["word"][sentence_id]

    word_content = vec(words["content"])
    fix_positions = vec(words["fixPositions"])
    raw_eeg_words = vec(words["rawEEG"])

    rows = NamedTuple[]

    for word_id in eachindex(fix_positions)

        fix_ids = get_fixation_ids(fix_positions[word_id])
        isempty(fix_ids) && continue

        segments = vec(raw_eeg_words[word_id])

        length(fix_ids) == length(segments) ||
            error(
                "word $word_id: $(length(fix_ids)) fixation IDs " *
                "but $(length(segments)) EEG segments"
            )

        for (fix_id, seg) in zip(fix_ids, segments)

            latency = find_segment_start(eeg, seg)

            push!(
                rows,
                (
                    fixation_id = fix_id,
                    word_id = word_id,
                    word = string(word_content[word_id]),
                    latency_local = latency,
                )
            )
        end
    end

    events = DataFrame(rows)
    sort!(events, :fixation_id)

    events
end

# ╔═╡ 085f5e33-6af0-4012-a224-01d2632da1e6
ev_ = local_events_for_sentence(sentenceData, 1)

# ╔═╡ aa6400aa-a4a9-4db0-8d08-c60a0edf4fcb
ev_

# ╔═╡ dfa95711-b13c-4bff-ba81-bc24ddc4b3d3
# ╠═╡ disabled = true
#=╠═╡
begin
    sid = 1

    ev = local_events_for_sentence(sentenceData, sid)
    raw_match = raw_matches[sid]

    offsets = Int[]

    for row in eachrow(ev)

        # sentence-final truncated segment etc.
        ismissing(row.latency_local) && continue

        fix_id = row.fixation_id

        raw_id = raw_match.matched[fix_id]
        ismissing(raw_id) && continue

        et_onset = rawfix[raw_id, :et_latency]

        global_sample =
            round(Int, sync.β[1] + sync.β[2] * et_onset)

        push!(
            offsets,
            global_sample - row.latency_local
        )
    end

    (
        n = length(offsets),
        unique_offsets = sort(unique(offsets)),
        min = minimum(offsets),
        max = maximum(offsets),
        range = maximum(offsets) - minimum(offsets),
        median = median(offsets),
    )
end
  ╠═╡ =#

# ╔═╡ 30ea467b-af7f-43cf-a1ef-05ffed99014f
begin
    sid = 1

    ev = local_events_for_sentence(sentenceData, sid)
    raw_match = raw_matches[sid]

    rows = NamedTuple[]

    for row in eachrow(ev)
        ismissing(row.latency_local) && continue

        fix_id = row.fixation_id
        raw_id = raw_match.matched[fix_id]
        ismissing(raw_id) && continue

        et_onset = rawfix[raw_id, :et_latency]

        global_sample =
            round(Int, sync.β[1] + sync.β[2] * et_onset)

        push!(rows, (
            fixation_id = fix_id,
            local_sample = row.latency_local,
            global_sample = global_sample,
            offset = global_sample - row.latency_local
        ))
    end

    offset_df = DataFrame(rows)
end

# ╔═╡ ad52dc91-9214-4735-92af-a031db4a4e52
let
    fig = Figure()

    ax = Axis(
        fig[1,1],
        xlabel = "Fixation ID",
        ylabel = "Global - local sample",
        title = "Sentence 1 offset"
    )

    scatterlines!(
        ax,
        offset_df.fixation_id,
        offset_df.offset
    )

    fig
end

# ╔═╡ 39cfaf55-13da-452a-b007-538caee37d59
let
	sid = 1
words = sentenceData["word"][sid]

rows = NamedTuple[]

for word_id in eachindex(vec(words["fixPositions"]))
    ids = get_fixation_ids(words["fixPositions"][word_id])
    isempty(ids) && continue

    segs = vec(words["rawEEG"][word_id])

    for (fix_id, seg) in zip(ids, segs)
        push!(rows, (
            fixation_id = fix_id,
            fixation_duration = 
                sentenceData["allFixations"][sid]["duration"][fix_id],
            segment_samples = size(seg, 2)
        ))
    end
end

sort(DataFrame(rows), :fixation_id)
end

# ╔═╡ Cell order:
# ╠═ddf1f343-818e-460b-b2c2-6a1f60f81311
# ╠═24be1c80-ee34-4afa-ae7d-a289971171a7
# ╟─d56938ec-8fdb-4f06-8d6f-9414b053f37d
# ╠═0c2a70df-8c3d-4a05-a410-8a4acda94f0a
# ╠═dd188df4-d0ed-4031-97d4-d61168f76c9e
# ╟─823a827b-db50-4bab-b87f-fdf44b2ffab8
# ╠═4d2d0861-bf08-48bc-b13b-43f3fc22f87b
# ╟─64f59fcf-4d9e-465d-b91d-bcbbf4ab8dfe
# ╠═1664b49b-aea4-4555-b8c1-c63a5dfe3386
# ╠═dd41afae-9bb9-4d3a-b8e9-0ce409ba171c
# ╠═fea44612-0ec1-49e8-9452-ceb1d6e70202
# ╠═9d099b5f-f52a-443f-b2e6-96e478df1221
# ╠═26718990-7719-4df6-8e1f-31ebdaf92518
# ╟─8d0b61b3-ad73-4a45-bd3d-2bf258f32911
# ╠═e0a1d9c3-e96c-4a73-bfed-2f04bc275752
# ╟─7d49ae2e-d5c3-4d26-ad84-4e48d5a96f7b
# ╠═162bd68c-28c3-4604-ac32-a51eed992b87
# ╟─3efff367-0bbc-4132-b9a8-2987a4eca62a
# ╠═bc068290-64ad-41d7-8379-191ea49cf1d0
# ╟─31c2cdc4-7a8f-4f2b-9ca6-758323a84f51
# ╠═580ec354-4e97-4f30-9315-bd83255b9c1f
# ╠═db706510-37f0-4f82-a1ee-471f31b2a5b8
# ╠═085f5e33-6af0-4012-a224-01d2632da1e6
# ╠═aa6400aa-a4a9-4db0-8d08-c60a0edf4fcb
# ╠═dfa95711-b13c-4bff-ba81-bc24ddc4b3d3
# ╠═7be98ca1-f507-4499-b196-ffc7c64d9ea8
# ╠═30ea467b-af7f-43cf-a1ef-05ffed99014f
# ╠═ad52dc91-9214-4735-92af-a031db4a4e52
# ╠═39cfaf55-13da-452a-b007-538caee37d59
