### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ f267e967-3eef-4893-897f-108115017d67

begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
end

# ╔═╡ 00ca98ba-7c31-4968-b00d-88dc028cd78c
Pkg.add("GLM")

# ╔═╡ c97fd78f-a998-464d-931b-5e79294b2a7c
begin
	
	using MAT
	using Statistics
	using DataFrames
	using CairoMakie
end

# ╔═╡ 09b3136e-9c0b-4e7a-a2e4-0e1553d93e2d
using GLM

# ╔═╡ d48579ee-8a7e-4be7-9759-dd610f22a151

md"""
# ZuCo 2.0 — raw EEG / eye-tracking synchronization

"""

# ╔═╡ c4803cb3-27ed-4736-bdf5-b4d96ac4964b

begin
    const ET_PATH = "/scratch/data/ZuCo2/raw/YRP_NR1_ET.mat"
    const EEG_PATH = "/scratch/data/ZuCo2/raw/YRP_NR1_EEG.mat"
end

# ╔═╡ bf46bb90-bdc0-4323-ba30-eb9440ed567f

begin
    et = matread(ET_PATH)
    eeg_file = matread(EEG_PATH)

    E = eeg_file["EEG"]
    fix = et["eyeevent"]["fixations"]
    ev = E["event"]

    sfreq = Float64(E["srate"])
end

# ╔═╡ ea8d0da6-3b74-4f64-b8db-f54f41649697
et_ev = et["event"]

# ╔═╡ aa7dcccc-894b-47a1-a16b-9cf257281159
eeg_ev = E["event"]

# ╔═╡ 037aac89-2824-4e7d-9e00-0df3c4b6edde
size(et_ev)

# ╔═╡ fcac5306-a807-4120-83af-dddc2972ac72
begin
	et_time = et["event"][:, 1]
	et_type = Int.(round.(et["event"][:, 2]))

	eeg_type = parse.(Int, strip.(string.(vec(ev["type"]))))
	eeg_latency = Float64.(vec(ev["latency"]))
end

# ╔═╡ 74de8a2e-d713-46b0-89e5-ab33144ccb78
length(et_type), length(eeg_type)

# ╔═╡ 0a5fbd88-f916-4b81-b31e-5d35dabb17dc
et_type == eeg_type

# ╔═╡ 0a31ccf1-aabe-4fcd-b62e-2eb880197135
[
    (
        i = i,
        et_time = et_time[i],
        et_type = et_type[i],
        eeg_type = eeg_type[i],
        eeg_sample = eeg_latency[i]
    )
    for i in 1:10
]

# ╔═╡ f363c3b0-32c6-4ed3-940a-99a83612b9e0
begin
	delta_et = et_time[3] - et_time[2]

delta_eeg_samples =
    eeg_latency[3] - eeg_latency[2]

delta_eeg_ms =
    delta_eeg_samples / sfreq * 1000

(
    delta_et_ms = delta_et,
    delta_eeg_samples = delta_eeg_samples,
    delta_eeg_ms = delta_eeg_ms,
)
end


# ╔═╡ e4253d8b-e9e9-43a4-b8af-a6edc1d34b93
begin
	slope12 =
    (eeg_latency[2] - eeg_latency[1]) /
    (et_time[2] - et_time[1])

slope12
end

# ╔═╡ bc9730bb-020a-4e34-a440-2d8663cfe3a2
sync_df = DataFrame(
	et_time = et_time,
	eeg_sample = eeg_latency
)

# ╔═╡ 90af6572-c8f0-40c3-9122-0efd42cceb2c
model = lm(@formula(eeg_sample ~ et_time), sync_df)

# ╔═╡ 5a39914e-9e1b-4fa9-8d0c-00b880ea8908
intercept = coef(model)[1]

# ╔═╡ cbb79bf0-0d68-41ff-b7f1-17b84781cbff
slope = coef(model)[2]

# ╔═╡ 3ec6e5b2-96b0-4106-92ec-7f0802c563a5
predicted = predict(model)

# ╔═╡ ede7763e-7670-410a-b537-f3ee8fc50ae1
(
    actual = eeg_latency[1],
    predicted = predicted[1]
)

# ╔═╡ 9f75679e-084c-4048-a7b8-9456d93e0f55
residual_samples = eeg_latency .- predicted

# ╔═╡ bf91f0e0-75a6-40d7-a50a-d38579a0ae69
residual_ms =
    residual_samples ./ sfreq .* 1000

# ╔═╡ f5cd4253-d2f2-4471-bd6e-3a7860c74e1e
(
    mean_abs_error_ms = mean(abs.(residual_ms)),
    median_abs_error_ms = median(abs.(residual_ms)),
    max_abs_error_ms = maximum(abs.(residual_ms))
)

# ╔═╡ c8a036a6-2fbd-44c0-a22f-7fb6de35838d
fix_latency = Float64.(fix["data"][:, 1]) # fix onset

# ╔═╡ ae512b0e-fab5-4745-8283-7cc61c8a2223
fix_eeg_sample =
    round.(Int, intercept .+ slope .* fix_latency)

# ╔═╡ 37c4dc4a-53eb-4b88-9cb6-f4bd687b8cd6
sync_fixation = DataFrame(
	fixation_onset = fix_latency,
	eeg_sample = fix_eeg_sample
)

# ╔═╡ e06f1305-c3c7-4da3-a8e6-faf477eeb9f0


# ╔═╡ 5b0e36e2-d006-4c65-9b30-931bb72b5bc3
begin
	PROCESSED_PATH = "/scratch/data/ZuCo2/resultsYRP_NR.mat"
	
	sentenceData = matopen(PROCESSED_PATH) do f
	    read(f, "sentenceData")
	end
end

# ╔═╡ 15f43d79-8b10-4c11-b571-24853d59e386
words = sentenceData["word"][2]

# ╔═╡ 4602af5a-9854-43ad-837a-97a06a2f40f9
keys(words)

# ╔═╡ 440a4373-e6c8-4caa-a9b4-4cdd68db9e63
content = vec(words["content"])

# ╔═╡ d683dce5-39c6-4672-8c45-f65748687121
fixpos = vec(words["fixPositions"])

# ╔═╡ f96028e1-b454-45c9-bde1-984b3c651cc0
content[1] 

# ╔═╡ b25aee36-64d3-4a34-965c-d7ddda2cd279
fixpos[1]

# ╔═╡ 6a07f09f-3aeb-4dbb-834d-1ea8f067a5c1
sentenceData

# ╔═╡ 7c140490-7c94-4b09-9a49-7413a15e6b73
allfix = sentenceData["allFixations"][2]

# ╔═╡ eb94e881-a8f7-4891-aa86-ed260dba6421
begin
	fixation_table = DataFrame(
	    fixation_id = 1:length(vec(allfix["duration"])),
	    duration = vec(allfix["duration"]),
	    x = vec(allfix["x"]),
	    y = vec(allfix["y"]),
	    pupilsize = vec(allfix["pupilsize"]),
	)
	
	fixation_table
end

# ╔═╡ b8228c3c-8a88-4966-845c-c0f3e738a18a
begin
    nfix = length(vec(allfix["duration"]))

    word_for_fix = Vector{Union{Missing, String}}(missing, nfix)
    word_id_for_fix = Vector{Union{Missing, Int}}(missing, nfix)

    for word_id in eachindex(content)
        pos = fixpos[word_id]

        isempty(pos) && continue

        fixation_ids =
            pos isa Number ? [Int(round(pos))] : Int.(round.(vec(pos)))

        for fix_id in fixation_ids
            word_for_fix[fix_id] = string(content[word_id])
            word_id_for_fix[fix_id] = word_id
        end
    end

    fixation_table.word_id = word_id_for_fix
    fixation_table.word = word_for_fix

    fixation_table
end

# ╔═╡ 2a8a698b-f110-4168-ac52-c630300e5ce8
begin
    for sentence_id in 1:length(sentenceData["word"])
        words = sentenceData["word"][sentence_id]
        allfix = sentenceData["allFixations"][sentence_id]

        content = vec(words["content"])
        fixpos = vec(words["fixPositions"])

        nfix = length(vec(allfix["duration"]))

        for word_id in eachindex(content)
            pos = fixpos[word_id]
            isempty(pos) && continue

            fixation_ids =
                pos isa Number ?
                [Int(round(pos))] :
                Int.(round.(vec(pos)))

            for fix_id in fixation_ids
                if fix_id > nfix
                    println(
                        "sentence = ", sentence_id,
                        " | word_id = ", word_id,
                        " | word = ", content[word_id],
                        " | fix_id = ", fix_id,
                        " | n_allFixations = ", nfix
                    )
                end
            end
        end
    end
end

# ╔═╡ 1a07fd28-a412-407a-b086-c454cfbb9881
begin
    words253 = sentenceData["word"][253]
    allfix253 = sentenceData["allFixations"][253]

    (
        n_fixations = length(vec(allfix253["duration"])),
        words = vec(words253["content"]),
        fixPositions = vec(words253["fixPositions"])
    )
end

# ╔═╡ 7f1c16c2-11bf-42f7-a946-57a82b81b47a
vec(words253["content"])[4]

# ╔═╡ 616e53a9-c70b-415f-bf90-8b49da1dbd74
vec(words253["fixPositions"])[4]

# ╔═╡ 18bd6bb1-4eac-4e3c-82ba-2d0518b1e544
begin
    fixpos253 = vec(words253["fixPositions"])

    ids253 = Int[]

    for pos in fixpos253
        isempty(pos) && continue

        ids = pos isa Number ?
            [Int(round(pos))] :
            Int.(round.(vec(pos)))

        append!(ids253, ids)
    end

    (
        fixation_ids = sort(unique(ids253)),
        n_unique_ids = length(unique(ids253)),
        n_allFixations = length(vec(allfix253["duration"]))
    )
end

# ╔═╡ a02c0ed5-96ef-43eb-987a-fddf093d6e36
begin
    rows = NamedTuple[]

    for sentence_id in 1:length(sentenceData["word"])

        words = sentenceData["word"][sentence_id]
        allfix = sentenceData["allFixations"][sentence_id]

        content = vec(words["content"])
        fixpos = vec(words["fixPositions"])

        duration = vec(allfix["duration"])
        x = vec(allfix["x"])
        y = vec(allfix["y"])
        pupil = vec(allfix["pupilsize"])

        for word_id in eachindex(content)

            pos = fixpos[word_id]
            isempty(pos) && continue

            fixation_ids =
                pos isa Number ?
                [Int(round(pos))] :
                Int.(round.(vec(pos)))

            for fix_id in fixation_ids
                push!(rows, (
                    sentence_id = sentence_id,
                    fixation_id = fix_id,
                    word_id = word_id,
                    word = string(content[word_id]),
                    duration = duration[fix_id],
                    x = x[fix_id],
                    y = y[fix_id],
                    pupilsize = pupil[fix_id],
                ))
            end
        end
    end

    fixation_words = DataFrame(rows)
    sort!(fixation_words, [:sentence_id, :fixation_id])

    fixation_words
end

# ╔═╡ 087ff87b-1cf1-4759-8057-06100baf08cd


# ╔═╡ d3f374b3-e6f1-47a3-b4d9-2619ed76662f
md"""
# Processed sentence and raw ET, how to match their fixations?
"""

# ╔═╡ 6f9d8d33-b735-491c-9e80-cfc9795ef314
begin
    rawfix = fix["data"]

    raw_duration = rawfix[:, 2] .- rawfix[:, 1]
    raw_x = rawfix[:, 4]
    raw_y = rawfix[:, 5]

    mapping_rows = NamedTuple[]
    summary_rows = NamedTuple[]

    # 从 raw ET 开头开始往后找
    raw_pointer = 1

    for sentence_id in 1:length(sentenceData["allFixations"])

        allfix = sentenceData["allFixations"][sentence_id]

        proc_duration =
            (vec(allfix["duration"]) .- 1) .* (1000 / sfreq)

        proc_x = vec(allfix["x"])
        proc_y = vec(allfix["y"])

        matched_raw = Int[]

        search_from = raw_pointer

        for fix_row in eachindex(proc_duration)

            found = 0

            # 从当前位置向后搜索
            for raw_id in search_from:size(rawfix, 1)

                duration_ok =
                    abs(raw_duration[raw_id] - proc_duration[fix_row]) <= 1

                x_ok =
                    abs(raw_x[raw_id] - proc_x[fix_row]) <= 0.5


                if duration_ok && x_ok
                    found = raw_id
                    break
                end
            end

            push!(matched_raw, found)

            if found > 0
                push!(mapping_rows, (
                    sentence_id = sentence_id,
                    processed_row = fix_row,
                    raw_fixation_id = found,
                    duration_ms = proc_duration[fix_row],
                    x = proc_x[fix_row],
                    y = proc_y[fix_row],
                ))

                search_from = found + 1
            end
        end

        good = matched_raw[matched_raw .> 0]

        raw_start = isempty(good) ? missing : minimum(good)
        raw_end   = isempty(good) ? missing : maximum(good)

        push!(summary_rows, (
            sentence_id = sentence_id,
            n_processed = length(proc_duration),
            n_matched = length(good),
            match_rate = length(good) / length(proc_duration),
            raw_start = raw_start,
            raw_end = raw_end,
        ))

        # 下一句话从上一句话最后一个 match 后继续找
        if !isempty(good)
            raw_pointer = maximum(good) + 1
        end
    end

    fixation_mapping = DataFrame(mapping_rows)
    mapping_summary = DataFrame(summary_rows)

    mapping_summary
end

# ╔═╡ 98700dc9-c6ba-47df-8439-49e0fa9c3a5e
# ╠═╡ disabled = true
#=╠═╡
begin
    proc_duration =
        (vec(sentenceData["allFixations"][1]["duration"]) .- 1) .* (1000 / sfreq)

    proc_x =
        vec(sentenceData["allFixations"][1]["x"])

    DataFrame(
        processed_fix = 1:10,
        processed_ms = proc_duration[1:10],
        raw_fix = 5:14,
        raw_ms = raw_duration[5:14],
        processed_x = proc_x[1:10],
        raw_x = raw_x[5:14],
    )
end
  ╠═╡ =#

# ╔═╡ 26d77af1-7a62-46b1-81e1-d0304f3a560c
begin
    nr1_summary = mapping_summary[1:50, :]

    (
        n_sentences = nrow(nr1_summary),
        n_perfect_sentences = count(nr1_summary.match_rate .== 1.0),
        total_processed_fixations = sum(nr1_summary.n_processed),
        total_matched_fixations = sum(nr1_summary.n_matched),
        overall_match_rate =
            sum(nr1_summary.n_matched) / sum(nr1_summary.n_processed),
        median_sentence_match_rate = median(nr1_summary.match_rate),
        worst_sentence_match_rate = minimum(nr1_summary.match_rate),
    )
end

# ╔═╡ 8570975b-9f13-42ff-bc1e-4beced56a426
filter(:match_rate => <(1.0), mapping_summary[1:50, :])

# ╔═╡ 0228bd78-ee86-422e-8489-7713b7534ca7
begin
    nr1_map = filter(:sentence_id => <=(50), fixation_mapping)
    sort!(nr1_map, [:sentence_id, :processed_row])

    mapped_eeg_samples =
        fix_eeg_sample[nr1_map.raw_fixation_id]

    (
        raw_fixations_temporally_ordered =
            issorted(nr1_map.raw_fixation_id),

        no_duplicate_raw_assignments =
            length(unique(nr1_map.raw_fixation_id)) == nrow(nr1_map),

        all_mapped_fixations_inside_eeg =
            all(
                (mapped_eeg_samples .>= 1) .&
                (mapped_eeg_samples .<= size(E["data"], 2))
            )
    )
end

# ╔═╡ 988c4ab1-5cf6-42cc-8edb-46867cba1b06
begin
	raw_et_fixations = DataFrame(
	    raw_fixation_id = 1:size(rawfix, 1),
	    onset = rawfix[:, 1],
	    endtime = rawfix[:, 2],
	    duration = rawfix[:, 2] .- rawfix[:, 1],
	    x = rawfix[:, 4],
	    y = rawfix[:, 5],
	    pupilsize = rawfix[:, 6],
	)
	
	first(raw_et_fixations, 10)
end

# ╔═╡ 0c8ec583-7ca0-4007-880c-74cff50f5a69
begin
	md"""raw ET"""
	et["eyeevent"]["fixations"]["data"]
end

# ╔═╡ f0bc5575-9da5-47c0-9daf-9d0e06eaa72d
sentenceData["allFixations"][sentence_id]

# ╔═╡ 3336b8ee-7ce3-40ff-965c-5655efc58628
begin
    _fix = et["eyeevent"]["fixations"]

    println("fixation fields:")
    println(keys(_fix))

    println("\ncolumns:")
    println(vec(_fix["colheader"]))

    _rawfix = _fix["data"]

    _rawfix[1:10, :]
end

# ╔═╡ Cell order:
# ╠═f267e967-3eef-4893-897f-108115017d67
# ╠═c97fd78f-a998-464d-931b-5e79294b2a7c
# ╟─d48579ee-8a7e-4be7-9759-dd610f22a151
# ╠═c4803cb3-27ed-4736-bdf5-b4d96ac4964b
# ╠═bf46bb90-bdc0-4323-ba30-eb9440ed567f
# ╠═ea8d0da6-3b74-4f64-b8db-f54f41649697
# ╠═aa7dcccc-894b-47a1-a16b-9cf257281159
# ╠═037aac89-2824-4e7d-9e00-0df3c4b6edde
# ╠═fcac5306-a807-4120-83af-dddc2972ac72
# ╠═74de8a2e-d713-46b0-89e5-ab33144ccb78
# ╠═0a5fbd88-f916-4b81-b31e-5d35dabb17dc
# ╠═0a31ccf1-aabe-4fcd-b62e-2eb880197135
# ╠═f363c3b0-32c6-4ed3-940a-99a83612b9e0
# ╠═e4253d8b-e9e9-43a4-b8af-a6edc1d34b93
# ╠═00ca98ba-7c31-4968-b00d-88dc028cd78c
# ╠═09b3136e-9c0b-4e7a-a2e4-0e1553d93e2d
# ╠═bc9730bb-020a-4e34-a440-2d8663cfe3a2
# ╠═90af6572-c8f0-40c3-9122-0efd42cceb2c
# ╠═5a39914e-9e1b-4fa9-8d0c-00b880ea8908
# ╠═cbb79bf0-0d68-41ff-b7f1-17b84781cbff
# ╠═3ec6e5b2-96b0-4106-92ec-7f0802c563a5
# ╠═ede7763e-7670-410a-b537-f3ee8fc50ae1
# ╠═9f75679e-084c-4048-a7b8-9456d93e0f55
# ╠═bf91f0e0-75a6-40d7-a50a-d38579a0ae69
# ╠═f5cd4253-d2f2-4471-bd6e-3a7860c74e1e
# ╠═c8a036a6-2fbd-44c0-a22f-7fb6de35838d
# ╠═ae512b0e-fab5-4745-8283-7cc61c8a2223
# ╠═37c4dc4a-53eb-4b88-9cb6-f4bd687b8cd6
# ╠═e06f1305-c3c7-4da3-a8e6-faf477eeb9f0
# ╠═5b0e36e2-d006-4c65-9b30-931bb72b5bc3
# ╠═15f43d79-8b10-4c11-b571-24853d59e386
# ╠═4602af5a-9854-43ad-837a-97a06a2f40f9
# ╠═440a4373-e6c8-4caa-a9b4-4cdd68db9e63
# ╠═d683dce5-39c6-4672-8c45-f65748687121
# ╠═f96028e1-b454-45c9-bde1-984b3c651cc0
# ╠═b25aee36-64d3-4a34-965c-d7ddda2cd279
# ╠═6a07f09f-3aeb-4dbb-834d-1ea8f067a5c1
# ╠═7c140490-7c94-4b09-9a49-7413a15e6b73
# ╠═eb94e881-a8f7-4891-aa86-ed260dba6421
# ╠═b8228c3c-8a88-4966-845c-c0f3e738a18a
# ╠═2a8a698b-f110-4168-ac52-c630300e5ce8
# ╠═1a07fd28-a412-407a-b086-c454cfbb9881
# ╠═7f1c16c2-11bf-42f7-a946-57a82b81b47a
# ╠═616e53a9-c70b-415f-bf90-8b49da1dbd74
# ╠═18bd6bb1-4eac-4e3c-82ba-2d0518b1e544
# ╠═a02c0ed5-96ef-43eb-987a-fddf093d6e36
# ╠═087ff87b-1cf1-4759-8057-06100baf08cd
# ╠═98700dc9-c6ba-47df-8439-49e0fa9c3a5e
# ╟─d3f374b3-e6f1-47a3-b4d9-2619ed76662f
# ╠═6f9d8d33-b735-491c-9e80-cfc9795ef314
# ╠═26d77af1-7a62-46b1-81e1-d0304f3a560c
# ╠═8570975b-9f13-42ff-bc1e-4beced56a426
# ╠═0228bd78-ee86-422e-8489-7713b7534ca7
# ╠═988c4ab1-5cf6-42cc-8edb-46867cba1b06
# ╠═0c8ec583-7ca0-4007-880c-74cff50f5a69
# ╠═f0bc5575-9da5-47c0-9daf-9d0e06eaa72d
# ╠═3336b8ee-7ce3-40ff-965c-5655efc58628
