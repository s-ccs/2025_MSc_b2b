### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ ce27129c-a21b-11f1-b1fe-8155391cd42a
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.instantiate()
end

# ╔═╡ 6cda1d52-6a6d-47e4-92d4-c22fc490a782
using MAT

# ╔═╡ 5feef9ca-905e-4d6b-a879-fc79055da196
using DataFrames

# ╔═╡ 024db723-fe98-4a00-bee1-7af8207dd3d3
path = "/scratch/data/ZuCo2/resultsYRP_NR.mat"

# ╔═╡ 2c284ce2-312b-4634-b7da-2ed8bcbf47d3
f = matopen(path)

# ╔═╡ 7255394c-7636-4370-83e3-427045b7d16a
keys(f)

# ╔═╡ a83b5885-0587-414e-8746-8636b32d652c
sentenceData = read(f, "sentenceData")

# ╔═╡ 89742f49-c517-43db-a1be-b57a3aaa4787
close(f)

# ╔═╡ 0b429793-9a7f-4bbe-bae8-119947ed4768
# ╠═╡ disabled = true
#=╠═╡
events = DataFrame(
	latency = ...,
	word = ...,
	fixation_duration = ...,
	surprisal = ...
)
  ╠═╡ =#

# ╔═╡ 711e7021-5a26-47a7-8661-005c7968d1f4
"""
content
word
allFixations
rawData
wordbounds
"""

# ╔═╡ 40011078-01e6-4a15-89f7-8b35631538ef
sentenceData["content"][1]

# ╔═╡ a4d45681-fb15-4898-ad8e-47df1082d12e
sentenceData["content"][2]

# ╔═╡ 416eb0c5-c27e-4865-93a8-3468bda3aa2a
typeof(sentenceData["rawData"][1])

# ╔═╡ 49d6643d-b399-406f-9811-3fe9b8ba8ad7
size(sentenceData["rawData"][1])

# ╔═╡ 6a2d1d43-fa16-492a-b70c-128d89145de9
typeof(sentenceData["word"][1])

# ╔═╡ d49c911d-1d6d-4177-ac4b-e8164762e8f9
sentenceData["word"][1]

# ╔═╡ 40010f08-daee-4331-9518-628b5260994b
sentenceData["word"][2]

# ╔═╡ 32819606-fad5-4945-8b4f-12e96e027c0f
typeof(sentenceData["allFixations"][1])

# ╔═╡ 1bd55dc4-e23b-4b6d-86d9-56b2b0c85e03
w1 = sentenceData["word"][1]

# ╔═╡ e0f886af-fc56-4750-a6c6-e5176f9d803b
vec(w1["content"])

# ╔═╡ 3283af8b-abfb-4bd8-bf32-710d82321d1f
# ╠═╡ disabled = true
#=╠═╡
fix1 = sentenceData["allFixations"][1]
  ╠═╡ =#

# ╔═╡ a7058890-ffc3-4075-aa45-838b3cec09f4
eeg1 = sentenceData["rawData"][1]

# ╔═╡ d661eb50-e3c1-40cd-8e57-e43440e8bd31
typeof(eeg1)

# ╔═╡ 1e94e0e8-34f1-4748-b7c3-8102ed96f36c
size(eeg1)

# ╔═╡ c679bf2f-76ab-4d87-958c-b4f2b9035943
vec(w1["fixPositions"])

# ╔═╡ 46734394-1ebc-40a2-9894-c093b0012f0c
sentenceData["wordbounds"][1]

# ╔═╡ 568deede-3a8e-4325-a5f4-55ca91a637e4
typeof(sentenceData["wordbounds"][1])

# ╔═╡ 9e566430-d299-477d-a7b9-5fe4201c0ec1
size(sentenceData["wordbounds"][1])

# ╔═╡ 653824cf-4b3d-412c-9087-2d95516c3f53
filter(
    x -> occursin("fix", lowercase(String(x))),
    collect(keys(w1))
)

# ╔═╡ 31f006fc-cd6c-4771-9ffb-cd89ec563345
fix1 = sentenceData["allFixations"][1]

# ╔═╡ aaee6823-5f37-42d2-9f69-9f735356166e
keys(fix1)

# ╔═╡ 8feed164-f3e2-44f1-98eb-31511fc7640a
for (k, v) in fix1
    if v isa AbstractArray
        println(k, " => ", typeof(v), "  ", size(v))
    else
        println(k, " => ", typeof(v))
    end
end

# ╔═╡ 99c7e410-1bdf-4844-baf6-4d55a34ddd08
keys(fix1)

# ╔═╡ b614716f-30c5-44b1-a45f-d278d32d7d3b
for k in keys(fix1)
    v = fix1[k]
    println(k, " => ", typeof(v),
            v isa AbstractArray ? "  size=$(size(v))" : "")
end

# ╔═╡ 28266464-deea-4949-b0c5-7cffc6c36322
filter(
    k -> any(s -> occursin(s, lowercase(String(k))),
             ["time", "onset", "latency", "position", "fix"]),
    collect(keys(w1))
)

# ╔═╡ d58650a1-a819-4b3a-82bf-bf1aa66698c2
keys(w1)

# ╔═╡ 4203d4e3-9338-468b-8c1d-15d4c32ecb4a
begin
	raw_eeg_words = w1["rawEEG"]
	raw_et_words  = w1["rawET"]

	typeof(raw_eeg_words)
	typeof(raw_et_words)
end

# ╔═╡ 019f2f86-b2e5-424e-add9-152e79d0f1e2
for i in 1:5
    println("word $i: ", w1["content"][i])
    println("  rawEEG: ", typeof(raw_eeg_words[i]))
    println("  rawET:  ", typeof(raw_et_words[i]))

    if raw_eeg_words[i] isa AbstractArray
        println("  EEG size: ", size(raw_eeg_words[i]))
    end

    if raw_et_words[i] isa AbstractArray
        println("  ET size:  ", size(raw_et_words[i]))
    end
end

# ╔═╡ bf8676c0-8f35-4297-bb4e-886c08ed428a
raw_eeg_words[1]

# ╔═╡ bedfd2ef-89d1-4cfa-a785-a973f5892538
raw_eeg_words[2]

# ╔═╡ fa0b7857-2927-4cd9-9b88-08ce238326a6
typeof(raw_eeg_words[1][1])

# ╔═╡ fbfc916c-b716-456f-9e15-dae0294dae64
size(raw_eeg_words[1][1])

# ╔═╡ 32b0a93d-4fe4-44e7-9c4a-a74254fe3d18
size(raw_eeg_words[1][2])

# ╔═╡ 830bcd33-55a9-4e80-a4f7-a73d2d7b45f6
size(raw_eeg_words[2][1])

# ╔═╡ 739065d7-1066-475e-b8db-e13719562c08
fix1["duration"][fixpos]

# ╔═╡ 25041e5d-e8d8-4734-960d-c6d1b6ddbfa7
[size(raw_eeg_words[1][j], 2) / 500 * 1000
 for j in eachindex(raw_eeg_words[1])]

# ╔═╡ 1a316450-a2f5-408c-b20b-dab01fb0b43e
eeg1

# ╔═╡ 0f848587-c427-4993-8d4b-b73fbf1259a4
begin
	seg = raw_eeg_words[1][1]
n = size(seg, 2)

best_start = 0
best_error = Inf

for s in 1:(size(eeg1, 2) - n + 1)
    err = maximum(abs.(view(eeg1, :, s:s+n-1) .- seg))

    if err < best_error
        best_error = err
        best_start = s
    end
end

(best_start, best_error)
end


# ╔═╡ b05be546-0aa7-4fd0-ae97-59c3479d4fb1
eeg1[:, 188:248] == raw_eeg_words[1][1]

# ╔═╡ 4a59b563-823e-4e8d-8cb1-2928ec1dd474
# ╠═╡ disabled = true
#=╠═╡
function find_segment_start(eeg, seg)
    n = size(seg, 2)

    best_start = 0
    best_error = Inf

    for s in 1:(size(eeg, 2) - n + 1)
        err = maximum(abs.(view(eeg, :, s:s+n-1) .- seg))

        if err < best_error
            best_error = err
            best_start = s
        end
    end

    return best_start, best_error
end
  ╠═╡ =#

# ╔═╡ 95807ea0-14f0-46d1-a344-4976064994fb
begin
	
	# ---------------------------------------------------------
	# Helper 1: normalize a word for lexical predictors later 
	# e.g. "Ford" -> "ford"
	# ---------------------------------------------------------
	
	normalize_word(w) = lowercase(replace(string(w), r"[^\p{L}']" => ""))
	
	
	# ----------------------------------------------------------
	# Helper 2: convert ZuCo fixPositions into Int fixation IDs
	# -----------------------------------------------------------
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
	
	
	# -------------------------------------------------
	# Helper 3: compare two EEG matrices exactly
	# --------------------------------------------------
	function same_segment(a, b)
	
		size(a) == size(b) || return false
	
	
		return all(a .== b)
	
	end
	
	
	# ----------------------------------------------------------
	# Helper 4: Find where a fixation EEG segment occurs 
	# inside the sentence-level continuous EEG.
	#
	# Returns the sentence-local sample index of fixation onset.
	# ------------------------------------------------------------
	
	function find_segment_start(eeg, seg)
	
		size(eeg, 1) == size(seg, 1) ||
			error("EEG and segment have different numbers of channels.")
	
		n = size(seg, 2)
		T = size(eeg, 2)
		n <= T || return missing
	
		# Pick one non-zero finitive value from the segment as an anchor.
		# This avoids comparing the complete 105-channel segment at every possible sample
		#anchor = nothing
	
#		for ch in axes(seg, 1), k in axes(seg, 2)
#			value = seg[ch, k]
#	
#			if isfinite(value) && value != 0
#				anchor = (ch, k, value)
#				break
#			end
#		end
#	
#		anchor === nothing && error("Could not find a usable anchor in EEG segment.")
#	
#		ch, k, value = anchor
#	
		matches = Int[]
	
		for start in 1:(T - n + 1)
	
			# Position of the anchor within sentence EEG
			# anchor_position = start + k - 1
	
			#eeg[ch, anchor_position] == value || continue
	
			candidate = @view eeg[:, start:start+n-1]
	
			if same_segment(candidate, seg)
				push!(matches, start)
			end 
		end
	
		if isempty(matches)
			return missing
			
		elseif length(matches) == 1
			return only(matches)
	
		else
			error("EEG segment occurs at multiple positions: $(matches)")
		end

	end
end


# ╔═╡ bcec45d6-75b1-462a-b283-31a7b766efcc
find_segment_start(eeg1, raw_eeg_words[1][1])

# ╔═╡ 67d5625e-827b-4e87-9872-402d2147fff9
find_segment_start(eeg1, raw_eeg_words[1][2])

# ╔═╡ 23ba4c47-1ad0-4ef0-8e68-4da1fddc01b2
find_segment_start(eeg1, raw_eeg_words[2][1])

# ╔═╡ 6c613bb1-542d-403d-ac78-3c7b409aae8a
# -----------------------------------------------------
# MAIN FUNCTION
#
# Build one fixation-level event table for one sentence
# One row = one fixation
# ---------------------------------------------

function build_events_for_sentence(
	sentenceData,
	sentence_id::Int;
	subject::String = "YRP",
)

	# ------------------------
	# 1. Sentence-level EEG
	# ------------------------
	eeg = sentenceData["rawData"][sentence_id]

	# ------------------------
	# 2. Word information
	# ------------------------
	words = sentenceData["word"][sentence_id]
	word_content = vec(words["content"])
	fix_positions = vec(words["fixPositions"])
	raw_eeg_words = vec(words["rawEEG"])


	# -------------------------
	# 3. Fixation information
	# -------------------------
	all_fixations = sentenceData["allFixations"][sentence_id]
	durations = vec(all_fixations["duration"])

	rows = NamedTuple[]

	# --------------------------
	# 4. Walk through words
	# -----------------------------
	for word_id in eachindex(word_content)

		word = string(word_content[word_id])
		word_lookup = normalize_word(word)

		# Which fixation(s) belong to this word?
		fixation_ids = get_fixation_ids(fix_positions[word_id])

		isempty(fixation_ids) && continue

		# EEG segment(s) belonging to those fixation(s)
		raw_segments = raw_eeg_words[word_id]

		raw_segments isa AbstractArray ||
			error(
				"Unexpected rawEEG structure for word $word_id: " *
				"$(typeof(raw_segments))"
			)

		segments = vec(raw_segments)

		length(segments) == length(fixation_ids) ||
			error(
				"Word $word_id ($word): " *
				"$(length(fixation_ids)) fixation IDs but " *
				"$(length(segments)) EEG segments."
			)


		# -------------------------
		# 5. One row per fixation
		# -------------------------
		for (j, fixation_id) in enumerate(fixation_ids)

			seg = segments[j]

			seg isa AbstractMatrix ||
				error(
					"EEG segment for word $word_id / fixation " *
					"$fixation_id is not a matrix."
				)

			latency = find_segment_start(eeg, seg)

			if latency === missing
				@warn "Could not locate fixation in sentence EEG" subject sentence_id word_id word fixation_id
			end

			#latency === missing &&
			#	error(
			#		"Could not locate word $word_id / fixation " *
			#		"$fixation_id inside sentence EEG."
			#	)

			push!(
				rows,
				(
					subject = subject,
					#fixation_id = fixation_id,
					sentence_id = sentence_id,

					word_id = word_id,
					word = word,
					#word_lookup = word_lookup,
					word_length = length(word_lookup),

					# onset with THIS sentence
					latency_local = latency,

					# keep original ZuCo value for now;
					# don't assume units yet
					fixation_duration = durations[fixation_id],

					# useful sanity check
					#segment_samples = size(seg, 2),
					
				)
				
			)
		end
	end

	events = DataFrame(rows)

	# IMPORTANT;
	# word order != fixation order because people regress/re-read.
	sort!(events, :latency_local)

	return(
		eeg = eeg,
		events = events,
	)
end

# ╔═╡ 15d3af92-f6a5-436a-a1d8-48152f7ac091
s1 = build_events_for_sentence(
	sentenceData,
	1;
	subject = "YRP",
)

# ╔═╡ 0b3ea4c4-5d7e-44a3-9122-fd2599483e6b
s1.events

# ╔═╡ d0e6443f-89ba-4d77-98f8-ade32255f297
begin
	seg30 = segments[1]
	
	size(seg30)
end

# ╔═╡ 547f7b12-cff8-4110-b4d0-6391e528572d
function find_best_segment_match(eeg, seg)

    n = size(seg, 2)

    best_start = missing
    best_error = Inf

    for start in 1:(size(eeg, 2) - n + 1)

        candidate = @view eeg[:, start:start+n-1]

        # compare only finite values
        mask = isfinite.(candidate) .& isfinite.(seg)

        any(mask) || continue

        err = maximum(
            abs.(candidate[mask] .- seg[mask])
        )

        if err < best_error
            best_error = err
            best_start = start
        end
    end

    return best_start, best_error
end

# ╔═╡ 4d642db3-8d93-4030-badd-cd16a6fd1b0f
begin
	eeg = sentenceData["rawData"][1]

find_best_segment_match(eeg, seg30)
end


# ╔═╡ 111f8cd7-96f9-4cfd-b4cf-b92eed9cab43
find_segment_start(eeg, raw_eeg_words[1][1])

# ╔═╡ dacbd322-df7c-4cb0-b870-dd4e3a32469b
find_segment_start(eeg, raw_eeg_words[1][2])

# ╔═╡ d17d1ae5-cf73-4d41-8a3d-a73dde1a4d78
find_segment_start(eeg, seg30)

# ╔═╡ bd63431b-7794-47f5-9a30-c8ecf5faf91f
begin
	candidate30 = @view eeg[:, 3314:3314+size(seg30, 2)-1]

println("seg nonfinite:       ", count(!isfinite, seg30))
println("candidate nonfinite: ", count(!isfinite, candidate30))

finite_seg = isfinite.(seg30)
finite_candidate = isfinite.(candidate30)

println(
    "different finite masks: ",
    count(finite_seg .!= finite_candidate)
)
end

# ╔═╡ 11d23fcd-2c6e-4e32-99a8-67e302a0d4ab
begin
	mismatch = findall(
    finite_seg .!= finite_candidate
)

first(mismatch, min(10, length(mismatch)))
end

# ╔═╡ bbd4e059-978a-4d3d-b710-c15b135d6270
candidate30 == seg30

# ╔═╡ b30a0a40-130c-471e-982a-80b69d64f93e
begin
	bad = [
    I for I in eachindex(seg30)
    if !isequal(candidate30[I], seg30[I])
]

length(bad)
end

# ╔═╡ 7d27b479-1779-4fbb-8892-795c4d5f38b6
[(I, candidate30[I], seg30[I])
 for I in first(bad, min(10, length(bad)))]

# ╔═╡ 1d05b255-1c21-4054-8b9c-efe705ad4f73
begin
	words = string.(vec(w1["content"]))
fixpos = vec(w1["fixPositions"])

# 先看看 allFixations 有多少个
allfix = sentenceData["allFixations"][1]
nfix = length(vec(allfix["duration"]))

word_idx_for_fix = Vector{Union{Missing, Int}}(missing, nfix)

for (word_idx, pos) in enumerate(fixpos)

    # 空 matrix = 这个词没有 fixation
    isempty(pos) && continue

    inds = Int.(round.(vec(pos)))

    for fix_idx in inds
        word_idx_for_fix[fix_idx] = word_idx
    end
end
end

# ╔═╡ b486ccb6-42ae-4eb1-9f4e-4cc69c8eb6f9
fixpos = Int.(vec(w1["fixPositions"][1]))

# ╔═╡ 8a9f5d5b-e05f-4d5d-975a-d61de84aaa2d
begin
    word_id = 10

    words = sentenceData["word"][1]

    fix_ids = get_fixation_ids(
        words["fixPositions"][word_id]
    )

    segments = vec(
        words["rawEEG"][word_id]
    )

    fix_ids
end

# ╔═╡ Cell order:
# ╠═ce27129c-a21b-11f1-b1fe-8155391cd42a
# ╠═6cda1d52-6a6d-47e4-92d4-c22fc490a782
# ╠═024db723-fe98-4a00-bee1-7af8207dd3d3
# ╠═2c284ce2-312b-4634-b7da-2ed8bcbf47d3
# ╠═7255394c-7636-4370-83e3-427045b7d16a
# ╠═a83b5885-0587-414e-8746-8636b32d652c
# ╠═89742f49-c517-43db-a1be-b57a3aaa4787
# ╠═5feef9ca-905e-4d6b-a879-fc79055da196
# ╠═0b429793-9a7f-4bbe-bae8-119947ed4768
# ╠═711e7021-5a26-47a7-8661-005c7968d1f4
# ╠═40011078-01e6-4a15-89f7-8b35631538ef
# ╠═a4d45681-fb15-4898-ad8e-47df1082d12e
# ╠═416eb0c5-c27e-4865-93a8-3468bda3aa2a
# ╠═49d6643d-b399-406f-9811-3fe9b8ba8ad7
# ╠═6a2d1d43-fa16-492a-b70c-128d89145de9
# ╠═d49c911d-1d6d-4177-ac4b-e8164762e8f9
# ╠═40010f08-daee-4331-9518-628b5260994b
# ╠═32819606-fad5-4945-8b4f-12e96e027c0f
# ╠═1bd55dc4-e23b-4b6d-86d9-56b2b0c85e03
# ╠═e0f886af-fc56-4750-a6c6-e5176f9d803b
# ╠═3283af8b-abfb-4bd8-bf32-710d82321d1f
# ╠═aaee6823-5f37-42d2-9f69-9f735356166e
# ╠═8feed164-f3e2-44f1-98eb-31511fc7640a
# ╠═a7058890-ffc3-4075-aa45-838b3cec09f4
# ╠═d661eb50-e3c1-40cd-8e57-e43440e8bd31
# ╠═1e94e0e8-34f1-4748-b7c3-8102ed96f36c
# ╠═c679bf2f-76ab-4d87-958c-b4f2b9035943
# ╠═46734394-1ebc-40a2-9894-c093b0012f0c
# ╠═568deede-3a8e-4325-a5f4-55ca91a637e4
# ╠═9e566430-d299-477d-a7b9-5fe4201c0ec1
# ╠═653824cf-4b3d-412c-9087-2d95516c3f53
# ╠═31f006fc-cd6c-4771-9ffb-cd89ec563345
# ╠═99c7e410-1bdf-4844-baf6-4d55a34ddd08
# ╠═b614716f-30c5-44b1-a45f-d278d32d7d3b
# ╠═28266464-deea-4949-b0c5-7cffc6c36322
# ╠═d58650a1-a819-4b3a-82bf-bf1aa66698c2
# ╠═4203d4e3-9338-468b-8c1d-15d4c32ecb4a
# ╠═019f2f86-b2e5-424e-add9-152e79d0f1e2
# ╠═1d05b255-1c21-4054-8b9c-efe705ad4f73
# ╠═bf8676c0-8f35-4297-bb4e-886c08ed428a
# ╠═bedfd2ef-89d1-4cfa-a785-a973f5892538
# ╠═fa0b7857-2927-4cd9-9b88-08ce238326a6
# ╠═fbfc916c-b716-456f-9e15-dae0294dae64
# ╠═32b0a93d-4fe4-44e7-9c4a-a74254fe3d18
# ╠═830bcd33-55a9-4e80-a4f7-a73d2d7b45f6
# ╠═b486ccb6-42ae-4eb1-9f4e-4cc69c8eb6f9
# ╠═739065d7-1066-475e-b8db-e13719562c08
# ╠═25041e5d-e8d8-4734-960d-c6d1b6ddbfa7
# ╠═1a316450-a2f5-408c-b20b-dab01fb0b43e
# ╠═0f848587-c427-4993-8d4b-b73fbf1259a4
# ╠═b05be546-0aa7-4fd0-ae97-59c3479d4fb1
# ╠═4a59b563-823e-4e8d-8cb1-2928ec1dd474
# ╠═bcec45d6-75b1-462a-b283-31a7b766efcc
# ╠═67d5625e-827b-4e87-9872-402d2147fff9
# ╠═23ba4c47-1ad0-4ef0-8e68-4da1fddc01b2
# ╠═95807ea0-14f0-46d1-a344-4976064994fb
# ╠═6c613bb1-542d-403d-ac78-3c7b409aae8a
# ╠═15d3af92-f6a5-436a-a1d8-48152f7ac091
# ╠═111f8cd7-96f9-4cfd-b4cf-b92eed9cab43
# ╠═dacbd322-df7c-4cb0-b870-dd4e3a32469b
# ╠═d17d1ae5-cf73-4d41-8a3d-a73dde1a4d78
# ╠═0b3ea4c4-5d7e-44a3-9122-fd2599483e6b
# ╠═8a9f5d5b-e05f-4d5d-975a-d61de84aaa2d
# ╠═d0e6443f-89ba-4d77-98f8-ade32255f297
# ╠═547f7b12-cff8-4110-b4d0-6391e528572d
# ╠═4d642db3-8d93-4030-badd-cd16a6fd1b0f
# ╠═bd63431b-7794-47f5-9a30-c8ecf5faf91f
# ╠═11d23fcd-2c6e-4e32-99a8-67e302a0d4ab
# ╠═bbd4e059-978a-4d3d-b710-c15b135d6270
# ╠═b30a0a40-130c-471e-982a-80b69d64f93e
# ╠═7d27b479-1779-4fbb-8892-795c4d5f38b6
