### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 18c10890-867d-43cf-a6ee-757556eabc47

begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ eb9fe0ce-c302-4a5e-bd10-8af82c4bce34
begin
	
	using MAT
	using DataFrames
	using Statistics
end

# ╔═╡ b00f28dc-162e-44d8-83f9-0eb6a531ee97

md"""
# ZuCo 2.0 — processed-data validation

Purpose:

1. verify that `fixPositions` gives a consistent fixation → word mapping;
2. recover fixation onset **within each processed sentence EEG** by locating the word-level `rawEEG` segment inside `sentenceData["rawData"]`;
3. quantify whether the processed ZuCo files are sufficient for continuous sentence-level Unfold analyses.

This notebook deliberately excludes the raw EEG/ET synchronization work. That is handled in `08_debug_zuco_raw_sync_clean.jl`.
"""

# ╔═╡ 7555975f-df2f-45a9-a187-f63830eb9868

begin
    const SUBJECT = "YRP"
    const PROCESSED_PATH = "/scratch/data/ZuCo2/resultsYRP_NR.mat"
    const SFREQ = 500.0
end

# ╔═╡ 36ddcca7-4f31-4c2e-9720-8cd7d253c5a9

sentenceData = matopen(PROCESSED_PATH) do f
    read(f, "sentenceData")
end

# ╔═╡ c48645ff-1cd7-428f-9b46-8012024e3acf

n_sentences = length(sentenceData["rawData"])

# ╔═╡ c98348f4-e226-4413-bdff-d4e297d5da13

md"""
## Helpers
"""

# ╔═╡ 5277344d-c1aa-4f35-8fd1-9f486c07ae53

begin
    normalize_word(w) = lowercase(replace(string(w), r"[^\p{L}\p{N}']" => ""))

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

    """
    Exact array equality that treats NaN == NaN for validation purposes.
    """
    same_segment(a, b) = size(a) == size(b) && isequal(a, b)

    """
    Find the unique sentence-local sample at which `seg` occurs in `eeg`.

    Uses one finite anchor value to avoid comparing the complete multichannel
    segment at every candidate sample. Returns `missing` if no exact match exists.
    """
    function find_segment_start(eeg, seg)
        size(eeg, 1) == size(seg, 1) ||
            error("EEG and segment have different channel counts.")

        n = size(seg, 2)
        T = size(eeg, 2)
        n <= T || return missing

        anchor = findfirst(isfinite, seg)
        anchor === nothing && return missing

        ch, k = Tuple(anchor)
        value = seg[anchor]
        matches = Int[]

        for start in 1:(T - n + 1)
            eeg_pos = start + k - 1
            isequal(eeg[ch, eeg_pos], value) || continue

            candidate = @view eeg[:, start:start+n-1]
            same_segment(candidate, seg) && push!(matches, start)
        end

        isempty(matches) && return missing
        length(matches) == 1 && return only(matches)

        error("EEG segment occurs at multiple positions: $(matches)")
    end

    """
    Invert ZuCo's word -> fixation IDs representation into fixation -> word ID.
    Returns the mapping plus diagnostics.
    """
    function invert_fixpositions(words, nfix::Int)
        mapping = Vector{Union{Missing, Int}}(missing, nfix)
        duplicates = Int[]
        invalid = Int[]

        for (word_id, pos) in enumerate(vec(words["fixPositions"]))
            for fix_id in get_fixation_ids(pos)
                if !(1 <= fix_id <= nfix)
                    push!(invalid, fix_id)
                    continue
                end
                if !ismissing(mapping[fix_id])
                    push!(duplicates, fix_id)
                else
                    mapping[fix_id] = word_id
                end
            end
        end

        return (
            word_id_for_fix = mapping,
            duplicate_fixations = unique(duplicates),
            invalid_fixations = unique(invalid),
        )
    end
end

# ╔═╡ c6a3b1b9-4fd4-45dc-be9d-e6d875fec1bd

md"""
## Build one fixation-level event table

One output row = one **word-assigned fixation**.

`fixation_id` is the processed `allFixations` index / reading-order index.
`latency_local` is the recovered onset sample within that sentence's processed continuous EEG.
"""

# ╔═╡ 2b855e74-cd71-4d27-97e3-ef669c45f163

function build_events_for_sentence(
    sentenceData,
    sentence_id::Int;
    subject::String = SUBJECT,
    sfreq::Real = SFREQ,
)
    eeg = sentenceData["rawData"][sentence_id]
    words = sentenceData["word"][sentence_id]
    all_fixations = sentenceData["allFixations"][sentence_id]

    word_content = vec(words["content"])
    fix_positions = vec(words["fixPositions"])
    raw_eeg_words = vec(words["rawEEG"])
    durations = vec(all_fixations["duration"])

    rows = NamedTuple[]

    for word_id in eachindex(word_content)
        word = string(word_content[word_id])
        word_lookup = normalize_word(word)

        fixation_ids = get_fixation_ids(fix_positions[word_id])
        isempty(fixation_ids) && continue

        raw_segments = raw_eeg_words[word_id]
        raw_segments isa AbstractArray ||
            error("Unexpected rawEEG structure for word $word_id: $(typeof(raw_segments))")

        segments = vec(raw_segments)

        length(segments) == length(fixation_ids) ||
            error(
                "Word $word_id ($word): $(length(fixation_ids)) fixation IDs but " *
                "$(length(segments)) EEG segments."
            )

        for (j, fixation_id) in enumerate(fixation_ids)
            seg = segments[j]
            seg isa AbstractMatrix ||
                error("rawEEG for word $word_id / fixation $fixation_id is not a matrix.")

            latency = find_segment_start(eeg, seg)

            push!(
                rows,
                (
                    subject = subject,
                    sentence_id = sentence_id,
                    fixation_id = fixation_id,
                    word_id = word_id,
                    word = word,
                    word_lookup = word_lookup,
                    word_length = length(word_lookup),
                    latency_local = latency,
                    fixation_duration_samples = durations[fixation_id],
                    segment_samples = size(seg, 2),
                    segment_duration_ms = size(seg, 2) / sfreq * 1000,
                ),
            )
        end
    end

    events = DataFrame(rows)

    if nrow(events) > 0
        sort!(
            events,
            :latency_local,
            by = x -> ismissing(x) ? typemax(Int) : x,
        )
    end

    return (eeg = eeg, events = events)
end

# ╔═╡ f7394b82-0fef-4104-bffd-696477635a0b

function validate_processed_sentence(
    sentenceData,
    sentence_id::Int;
    subject::String = SUBJECT,
    sfreq::Real = SFREQ,
)
    words = sentenceData["word"][sentence_id]
    allfix = sentenceData["allFixations"][sentence_id]
    n_all = length(vec(allfix["duration"]))

    inv = invert_fixpositions(words, n_all)
    assigned_ids = findall(!ismissing, inv.word_id_for_fix)

    result = build_events_for_sentence(
        sentenceData,
        sentence_id;
        subject = subject,
        sfreq = sfreq,
    )
    events = result.events

    n_assigned = length(assigned_ids)
    n_recovered = nrow(events) == 0 ? 0 : count(!ismissing, events.latency_local)
    recovery_rate = n_assigned == 0 ? NaN : n_recovered / n_assigned

    recovered = nrow(events) == 0 ? events : events[.!ismissing.(events.latency_local), :]
    temporal_order_ok =
        nrow(recovered) <= 1 ? true : issorted(recovered.fixation_id)

    duration_error_ms =
        nrow(events) == 0 ? Float64[] :
        abs.(
            events.segment_duration_ms .-
            events.fixation_duration_samples .* (1000 / sfreq)
        )

    return (
        sentence_id = sentence_id,
        n_all_fixations = n_all,
        n_word_assigned = n_assigned,
        n_unassigned = n_all - n_assigned,
        n_duplicate_assignments = length(inv.duplicate_fixations),
        n_invalid_fixation_ids = length(inv.invalid_fixations),
        n_latency_recovered = n_recovered,
        recovery_rate = recovery_rate,
        temporal_order_ok = temporal_order_ok,
        median_duration_error_ms =
            isempty(duration_error_ms) ? NaN : median(duration_error_ms),
        events = events,
    )
end

# ╔═╡ d61a57cb-cc58-4cff-a4bf-8fb7a4b050a3

md"""
## Quick validation on several sentences

Interpretation:

- `recovery_rate ≥ 0.95`: strong evidence that processed `rawEEG` can recover sentence-local onsets.
- duplicate/invalid fixation assignments should be zero.
- `temporal_order_ok == true`: recovered EEG onset order agrees with the fixation IDs stored by ZuCo.
"""

# ╔═╡ 1fc97a3d-0656-4bcc-9fb5-f21077fecda8

test_sentences = filter(<=(n_sentences), [1, 2, 3, 4, 5, 10, 20, 30])

# ╔═╡ 1e1eb279-5656-424b-be0d-f2c47b0e04b5
begin
	begin
		
		processed_checks = [validate_processed_sentence(sentenceData, sid) for sid in test_sentences]
		
		processed_summary = DataFrame([
		    (
		        sentence_id = x.sentence_id,
		        n_all_fixations = x.n_all_fixations,
		        n_word_assigned = x.n_word_assigned,
		        n_unassigned = x.n_unassigned,
		        n_duplicate_assignments = x.n_duplicate_assignments,
		        n_invalid_fixation_ids = x.n_invalid_fixation_ids,
		        n_latency_recovered = x.n_latency_recovered,
		        recovery_rate = x.recovery_rate,
		        temporal_order_ok = x.temporal_order_ok,
		        median_duration_error_ms = x.median_duration_error_ms,
		    )
		    for x in processed_checks
		])
	end
end

# ╔═╡ 6b4b2c0e-a05e-4b26-8156-28369eaeafcc

first_sentence_events = processed_checks[1].events

# ╔═╡ 942cfc61-c49e-436d-b202-d09c9288aeb5
begin
	
	processed_pass =
	    all(processed_summary.n_duplicate_assignments .== 0) &&
	    all(processed_summary.n_invalid_fixation_ids .== 0) &&
	    all(processed_summary.temporal_order_ok) &&
	    all(processed_summary.recovery_rate .>= 0.95)
	
	(
	    processed_pass = processed_pass,
	    median_recovery_rate = median(processed_summary.recovery_rate),
	)
end

# ╔═╡ b0df3385-8d94-4061-b1e4-a967d7513171
filter(
    row ->
        row.n_duplicate_assignments != 0 ||
        row.n_invalid_fixation_ids != 0 ||
        !row.temporal_order_ok ||
        row.recovery_rate < 0.95,
    processed_summary
)

# ╔═╡ 13d0cb29-cf2c-4849-acc0-b42fa111daf6
processed_summary

# ╔═╡ 8cbc297b-f451-4453-8eb9-fa0d42c48364
sort(processed_summary, :recovery_rate)

# ╔═╡ 995d86b0-cd5d-44c0-b214-34f269d8f690

md"""
### Decision

If `processed_pass == true`, the processed ZuCo files are sufficient to construct:

`processed continuous sentence EEG + fixation onset + word identity`

without rebuilding fixation-to-word assignment from `wordbounds`.
"""

# ╔═╡ a4cb4717-d363-4b44-a3e0-3235f7ae874c
for check in processed_checks
    if check.recovery_rate < 1
        bad = filter(
            row -> ismissing(row.latency_local),
            check.events
        )

        println("sentence ", check.sentence_id)
        display(bad)
    end
end

# ╔═╡ c6e494c5-a78d-4dcd-9b99-52859f2cbb24
for check in processed_checks
    bad = filter(
        row -> ismissing(row.latency_local),
        check.events
    )

    isempty(bad) && continue

    println(
        "sentence ", check.sentence_id,
        " | bad fixation = ", bad.fixation_id,
        " | max fixation = ", maximum(check.events.fixation_id)
    )
end

# ╔═╡ Cell order:
# ╠═18c10890-867d-43cf-a6ee-757556eabc47
# ╠═eb9fe0ce-c302-4a5e-bd10-8af82c4bce34
# ╟─b00f28dc-162e-44d8-83f9-0eb6a531ee97
# ╠═7555975f-df2f-45a9-a187-f63830eb9868
# ╠═36ddcca7-4f31-4c2e-9720-8cd7d253c5a9
# ╠═c48645ff-1cd7-428f-9b46-8012024e3acf
# ╟─c98348f4-e226-4413-bdff-d4e297d5da13
# ╠═5277344d-c1aa-4f35-8fd1-9f486c07ae53
# ╟─c6a3b1b9-4fd4-45dc-be9d-e6d875fec1bd
# ╠═2b855e74-cd71-4d27-97e3-ef669c45f163
# ╠═f7394b82-0fef-4104-bffd-696477635a0b
# ╟─d61a57cb-cc58-4cff-a4bf-8fb7a4b050a3
# ╠═1fc97a3d-0656-4bcc-9fb5-f21077fecda8
# ╠═1e1eb279-5656-424b-be0d-f2c47b0e04b5
# ╠═6b4b2c0e-a05e-4b26-8156-28369eaeafcc
# ╠═942cfc61-c49e-436d-b202-d09c9288aeb5
# ╠═b0df3385-8d94-4061-b1e4-a967d7513171
# ╠═13d0cb29-cf2c-4849-acc0-b42fa111daf6
# ╠═8cbc297b-f451-4453-8eb9-fa0d42c48364
# ╠═995d86b0-cd5d-44c0-b214-34f269d8f690
# ╠═a4cb4717-d363-4b44-a3e0-3235f7ae874c
# ╠═c6e494c5-a78d-4dcd-9b99-52859f2cbb24
