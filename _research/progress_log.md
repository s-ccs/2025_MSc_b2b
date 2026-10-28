# Progress Log

## May 2026 — Simulation setup and initial B2B checks

### Progress
- Built the first simulation prototype for testing:
  - temporal overlap,
  - condition–continuous confounding,
  - B2B predictor recovery.
- Added controllable event timing / overlap and condition–continuous correlation.
- Used simulated ERP components as known ground truth:
  - `condition → N170`
  - `continuous → P300`
- Added basic sanity checks for the simulated predictor relationships.


### Key findings 
- Investigated the sign of B2B coefficients.
- Current interpretation: B2B coefficients mainly reflect **predictor recovery strength**, rather than ERP-like signed amplitudes.
- For the current analyses, recovery magnitude is therefore more informative than directly interpreting the sign.

### Next
- Began separating the exploratory Pluto prototype from reusable simulation code.
- Initial modular structure included event design, onset models, ERP components, pipeline code, metrics, and scripts.

```
src/
├── B2BSim.jl          # main entry/ module files
├── SimConfig.jl       # simulation configuration
├── EventDesign.jl     # condition / continuous / collinearity design
├── OnsetModels.jl     # overlap / onset distribution
├── Components.jl      # N170 / P300 ground-truth components
├── Pipelines.jl       # simulate -> epoch -> fit B2B piplines
├── Metrics.jl         # recovery, magnitude leakage, correlation 
└── IOUtils.jl         # save csv / IO / path 

scripts/
├── run_single_case.jl # single setting，debug use
├── run_grid.jl        # multiple overlap × covariate × seed
├── summarize_grid.jl  # summarize results
└── make_figures.jl    # plotting from results
```

## June 2026 - Pipeline design and single-trial reconstruction

### Planned pipeline structure

gradually building up: 

- standard decoding on epoched EEG,
- rERP decoding on continuous firbasis designed EEG,
- plain B2B on epoched EEG,
- rERP reconstruction followed by B2B,
- one-step FIR-based B2B on continuous EEG. 




## July 2026 Simplification of the simulation and decoding benchmark

### Simulation
- two-predictor setting: categorical `condition`, `continuous` predictor
- Defined 4 main cases: `clean`, `overlap`, `confound`, `both` 
- Standard ridge decoding 
- rERP decoding

### Key findings / decisions
- Realised that the earlier simulation/pipeline structure was too complicated.
- Shifted toward a smaller benchmark where each bias could be tested separately.
- Spent substantial time understanding the `UnfoldDecode` source code, especially:
  - how time-point decoding is fitted,
  - how `MLJ.machine` is used internally,
  - how `coeftable` and decoding measures are returned.
- By the end of July, the overall pipeline design was clearer, but the implementation was still being debugged.


## Reduced activity
- progress during end of June and part of July was limited due to some circumstances. 
- After this period, the implementation strategy was simplified and rebuilt from basic, and started again from Pluto.



## August 2026 — Rebuilding and validating the decoding pipelines

### Progress
- Refactored the project into the local `MScB2B` package.
- Reworked the simulation so that each case provides both continuous and epoched data from the same simulated EEG.
- Implemented and debugged:
  - standard ridge decoding,
  - rERP decoding(partially debugged),
  - plain B2B.
- Added reusable plotting and interactive simulation controls.
- Pearson correlation r for decoding performance 
- Ridge hyperparameter tuning continues to use cross-validated RSquared()
- The simulation and the main decoding pipelines are now largely implemented.


## 17 Aug 2026 — Plain B2B debugging

### Progress
- Investigated an unexpected condition-B2B peak around ~200 ms.
- Checked cross-validation repetitions, channel count, and predictor coding; the peak remained.
- Removing the shared P300 intercept (`β₀,P300 = 0`) removed the peak while keeping the N170 intercept unchanged.

## 25 Aug 2026 - Five pipelines, numerical debugging

### Progress
- Implement the two-step overlap-corrected B2B pipeline: continuous EEG → rERP overlap correction → reconstructed single trials → B2B
- All five pipelines are now implemented and producing results.
- Generated comparable plots for `condition` and `continuous` in all conditions, see [../notebooks/06_compare_pipelines.jl](../notebooks/06_compare_pipelines.jl)
- Current focus is numerical validation/debugging of the B2B results.


## 01 Sep 2026 - ZuCo2 data explorarion

### Progress
- Added Hanning-window simulation with 1500 trials as a controlled temporal-support sanity check.❗ with current sampling rate and chosen size of basis function, there are aliasing/sampling artefacts.
- ❗ Why hanning-window shows up as much better decodable? 
- Explored ZuCo 2.0 NR with pilot subject YRP.
- Started building fixation-level events from processed Matlab data.
- Fixation onset is not stored directly; tested reconstructing it from fixation EEG segments.
- Proposed predictors: word length, SUBTLEX-US frequency, GPT-2 surprisal.
- Compared ZuCo2 with ROAMM, and decided to switch to [ROAMM](https://openneuro.org/datasets/ds007629/versions/1.3.0), it provides longer natural reading texts and synchronized EEG-eyetracking data.
- Reference: [ROAMM tutorials](https://data-brain-mind.github.io/tutorials/reading-observed-at-mindless-moments-roamm-a-simultaneous-eeg-and-eye-tracking-dataset-of-natural-reading-with-attention-annotations/)

## 08 Sep 2026 - ROAMM preprocessing

### Progress

#### ROAMM preprocessing
- Build a ROAMM preprocessing pipeline using Python in a separate repo: [Natural-Reading-Lexical-Features](https://github.com/xuyg16/Natural-Reading-Lexical-Features)
- Current scripts:
  > - `01_prepare_text.py`
  > - `02_compute_length_frequency.py`
  > - `03_compute_surprisal.py`
  > - `04_prepare_fixation_events.py`
  > - `05_merge_lexical_features.py`

- Generated word-level lexical predictors and fixation-level event tables aligned to EEG samples.
- Current predictors: word length, Zipf frequency, and GPT-2 surprisal.
- Next: test the ROAMM event table with B2B and decide whether to merge this preprocessing code back into the main repo.

#### Simulation / B2B debugging
- Compared all five pipelines under `clean`, `overlap`, `confound`, and `both`.
- Plain B2B does not correct temporal overlap
- One-step B2B:
  - Produces temporally narrow estimates compared with Plain and Two-step B2B.
  - There are two bumps when decoding continuous predictor in confounding condition.
  - It seems like deconvolution corrects temporal overlap but not predictor correlation.
- Two-step B2B:
  - Currently shows the most stable recovery.
  - Needs quantiative metrics.
 

## 15 Sep 2026 - ROAMM real-data pipeline

### Progress
- Added a sampling-rate check in the simulation; increasing the sampling rate does not remove the spike.
- Prepared the full ROAMM dataset for the real-data B2B analysis.
  - Converted all synced EEG runs to `.npy` format.
  - Current dataset: 44 subjects, 220 runs, 6000 trials per subject. 
- Implemented a provisional Two-step B2B pipeline for ROAMM:
  - run-wise Unfold overlap correction
  - extraction of corrected single trials
  - subject-level B2B across runs
  - successfully tested on `sub-10014` using all 5 runs.
  - noisy result, very samll estimates values.
- Compared LSQ and Ridge on a small real-data example.
  - LSQ is relatively stable under global EEG rescaling.
  - Ridge is much more sensitive to feature scaling.
- ❗ EEG data scaling, V or uV?

### Meeting note
- Do the sanity check on EEG data, what is the quality of data.
- Check the intercept, check channels, roughly go through some single channels.
- centrlize word_length, frenquency_zipf
- check saccade amplitude when page changes.


## 22 Sep 2026 - ROAMM sanity check & additional correlated predictors simulation

### Progess
- Ran sanity check on EEG data
  - all-channel FRP butterfly plot
  - raw FRP vs. deconvolved FRP
  - FRPs at posterior electrodes O1, O2, Oz, POz
- Applied b2b on one subject:
  - compared 1 run vs. 5 runs
  - using 5 runs reduced noise
- Add an additional simulation with 3 correlated continuous predictors:
  - p100, n170, p300
  - no bised overlap
  - plain b2b and two-step b2b show the unexplainable peaks

### Meeting note
- run sanity check scripts on all EEG data, to see if there is any unsual and suspicious data
- Raw FRP from a singe run looked better than deconvolved FRPs (this is not the central to the current RQ, no need to investigate further)
- why the unexpected peaks? change UnfoldSim default shape to Hanning window
- set predictor-biased overlap.


## 29 Sep 2026 - Midterm talk 

### Progress
- Updated the 3 correlated continuous predictors simulation, with hanning window and continuous-1 biased overlap
- Gave the midterm talk

### Feeback / questions to investigate


#### Simulation
- Show the simulation ground truth more explicitly:
  - visualize the correlation between condition and the continuous predictor
  - show how the predictors are correlated with each other

- One-step B2B does not separate the three correlated predictors as cleanly as
  plain B2B or two-step B2B.
  - Hypothesis: this may be related to strong regularization.
  - With the much longer continuous/FIR design matrix, regularization may shrink
    the signal strongly, leaving little decodable information.
  - Need to verify this rather than assuming it is the mechanism.

- Understand the effect of `ShiftOnsetByOne` more precisely:
  - why claim the first bump is from the overlap and second bump is the overlap from next fixation?
  - Why do the overlap-related bumps appear where they do?
  - Which neighboring event contributes to each bump?
  - How does shifting the onset-distance assignment change these bumps?
  - Verify this from the UnfoldSim implementation rather than inferring it only
    from the decoding curves.

#### ROAMM predictors
- Surprisal:
  - Why use GPT-2 rather than a newer language model?
  - Why compute surprisal using sentence-level context?
  - Be able to justify the model and context-window choice.

#### Methods / code
- Understand the feature-importance implementation:
  - trace through the code
  - understand exactly what quantity is being estimated and how it should be interpreted

## 06 Oct 2026 

### Progress
- Checked the original code, the code used a more sufficient way to define a ridge function, and also applied SVD decompose matrix.
- Sanity check and quality check of ROAMM data again, especially first-pass reading vs rereading, left and right eye.
- Due to the time limit, consider just evaluate S as main analysis, delta R as complementary analysis.
- Limitation found: a fixation contaion multiple information of multiple words.

- Why use GPT2:
  - The goal is not to maximize next-word prediction accuracy, but to model human language processing. Previous work has shown that larger and better-performing language models do not necessarily provide a better fit to human reading behaviour, and GPT-2 has been shown to provide a good predictor of human reading performance. It is also autoregressive, so its surprisal is based only on preceding context, which is appropriate for incremental reading.
  - Recalculate the word surprisal, because calculating surprisal naively and subtle tokenization issues can distort the results.

### Meeting notes
- run the model first, and then decide whether evaluate delta R or not. statistical test not decided yet.



## 13 Oct 2026

### Progress
- Extend the real-data analysis from the intial single-subject tests to the full dataset for all the pipelines.
- Insepected page/fixation boundaries in ROAMM event table, there is no overlapped boundaries.
- Set up cluster, using HPC, slurm bash scripts
- concatenated 5 runs in one subjects,  used centerised data
- Debug pipelines
- Ran pipelines on cluster

### Meeting
- DAG, suppessor effect, winsorized mean, logscale
- Do some research about the statistical test. what should we test, and how shoud we test? 


## 20 Oct 2026
### Progress
- The one-step B2B has not yet been included because it is computationally expensive.
- Subject-level outputs are now availble for group-level aggregation and statistical analysis.
- Inspected and summairzed the subject-level estimates across participants.
- Started writing the thesis.
- Drafted analysis decisions: 

### Meeting
- Meeting moved from evevry Tuesday to Thursday. 

## 30 Oct 2026

## Current status

## 20 Oct 2026 — ROAMM full-subject analysis

### Progress

#### ROAMM

- Extended the real-data analysis from the initial single-subject tests to the full dataset.
- Ran the main analysis pipelines across **44 ROAMM subjects**.
- Current full-subject results are available for **four pipelines**.
- The one-step FIR + B2B pipeline has not yet been included in the full ROAMM batch because the continuous FIR design is substantially more computationally expensive.
  - Need to decide whether it is necessary for the final thesis analysis.
- Subject-level outputs are now available for group-level aggregation and statistical analysis.

#### Predictors

- Continued preparing the lexical predictors for the ROAMM analyses:
  - word length
  - Zipf word frequency
  - GPT-2 surprisal
- Surprisal is now part of the real-data analysis workflow.
- Continued checking how lexical information should be associated with fixation events, especially given that neural activity during one fixation may also contain parafoveal processing of neighbouring words.

### Current interpretation

- The main goal of the ROAMM analysis is not to interpret an FRP as the response to exactly one isolated word.
- Instead, EEG is time-locked to fixation onset, while the lexical predictors describe properties of the currently fixated word.
- Processing of neighbouring/parafoveal words remains a limitation of this event definition and should be discussed in the thesis.

## 27 Oct 2026 — Final analysis planning and start of thesis writing

### Simulation

- The main simulation analyses are largely complete.

- Results are available for the main overlap × collinearity conditions and for the correlated-continuous-predictor simulation.

- The remaining simulation work is mainly:

  - final sanity checks,
  - choosing quantitative evaluation measures,
  - generating final figures,
  - deciding which exploratory results should be included in the thesis.

- The one-step B2B results remain difficult to interpret:

  - in the correlated-predictor simulation it does not separate the predictors as cleanly as plain B2B and two-step B2B,
  - it also produces additional bumps whose mechanism is not yet fully understood.

- Since one-step B2B is not currently planned as the main final pipeline, these effects may be treated as secondary/exploratory rather than investigated extensively.

### ROAMM

- Completed the computationally feasible full-subject analyses for **44 subjects**.
- The next step is no longer basic pipeline implementation, but **group-level evaluation and statistical inference**.
- Need to inspect and summarize the subject-level estimates across participants.
- Main unresolved analysis decisions:
  - statistical test for group-level B2B estimates,
  - whether to calculate **ΔR** as a complementary measure to S,
  - whether a separate train/test or held-out evaluation is necessary,
  - which comparisons between pipelines should form the main thesis result.

### Writing

- Started writing the thesis.
- Began drafting the **simulation methods** section.
- Current writing focus:
  - describe the simulated data and ground truth,
  - explain temporal overlap and predictor collinearity manipulations,
  - describe the decoding/B2B pipelines clearly,
  - avoid adding unnecessary implementation details before the final analysis decisions are fixed.

# Current status — 28 Oct 2026

## ✅ Implemented / completed

### Simulation

- [x] Main condition + continuous simulation
- [x] `clean` / `overlap` / `confound` / `both` conditions
- [x] Temporal-overlap manipulation
- [x] Predictor-collinearity manipulation
- [x] Standard ridge decoding
- [x] rERP decoding
- [x] Plain B2B
- [x] One-step FIR + B2B
- [x] Two-step rERP → reconstructed single trials → B2B
- [x] Comparison plots across pipelines
- [x] Sampling-rate sanity check
- [x] Three-correlated-continuous-predictor simulation
- [x] Hanning-window components
- [x] Predictor-biased overlap
- [x] Main simulation results generated

### ROAMM preprocessing

- [x] Prepare text
- [x] Compute word length
- [x] Compute Zipf frequency
- [x] Compute GPT-2 surprisal
- [x] Prepare fixation-level events
- [x] Merge lexical predictors with fixation events
- [x] Align fixation events to EEG samples
- [x] Export synchronized EEG data
- [x] Initial EEG / FRP sanity checks

### ROAMM analysis

- [x] Single-run analysis
- [x] Multi-run subject-level analysis
- [x] Two-step B2B tested successfully on real EEG
- [x] Full-subject analysis across **44 subjects**
- [x] Four main pipelines run across subjects
- [x] Subject-level results saved for group analysis

### Method understanding

- [x] Understand the basic B2B procedure:
  - forward model G
  - backward model H
  - predictor-recovery matrix S
- [x] Clarified that B2B S should be interpreted primarily as predictor-recovery / feature-importance strength rather than ERP amplitude.
- [x] Investigated the role of ridge regularization.
- [x] Investigated the SVD-based implementation used in the original code.
- [x] GPT-2 choice justified for modelling incremental human reading.
- [x] Identified parafoveal / multiple-word processing as a limitation of fixation-locked lexical analyses.

## 🔍 Remaining questions

### Simulation

- [ ] Finalize numerical evaluation of recovery.
- [ ] Decide which quantitative metrics are necessary for the final thesis.
- [ ] Decide how much of the one-step B2B behaviour needs to be explained.
- [ ] Finalize figures and figure layout.
- [ ] Make sure the simulation ground truth and predictor correlations are clearly visualized.

### ROAMM

- [ ] Aggregate subject-level results.
- [ ] Produce group-level figures for the main predictors.
- [ ] Inspect between-subject variability.
- [ ] Decide whether one-step B2B needs to be run on the complete ROAMM dataset.
- [ ] Select the final subset of real-data analyses to present in the thesis.

### Evaluation / statistics

- [ ] Decide the main statistical test for group-level S.
- [ ] Decide whether **ΔR** should be included as a complementary analysis.
- [ ] Decide whether held-out / train-test evaluation is necessary.
- [ ] If required, implement the minimum train/test analysis.
- [ ] Define the final comparisons between plain B2B and overlap-corrected B2B.

### Writing

- [x] Thesis writing started.
- [ ] Finish simulation Methods.
- [ ] Write ROAMM dataset and preprocessing Methods.
- [ ] Write B2B / overlap-correction Methods.
- [ ] Write simulation Results.
- [ ] Write ROAMM Results.
- [ ] Literature review / Introduction.
- [ ] Discussion.

## 🚩 Immediate priorities

1. **Finalize the statistical-analysis plan**

   - Ask about the group-level statistical test.
   - Decide whether ΔR is necessary.
   - Decide whether train/test evaluation is necessary.

2. **Finish the minimum remaining analyses**

   - implement the agreed statistical test,
   - implement ΔR / held-out evaluation only if needed,
   - avoid adding large exploratory analyses unless they directly contribute to the thesis.

3. **Aggregate the 44-subject ROAMM results**

   - group mean,
   - subject variability,
   - main predictor effects,
   - pipeline comparison.

4. **Freeze the simulation analysis**

   - final checks,
   - final metrics,
   - final figures.
   - Avoid further simulation expansion unless a result cannot otherwise be justified.

5. **Continue writing in parallel**

   - finish the simulation Methods first,
   - then ROAMM Methods,
   - write Results as soon as each analysis is frozen.





> [!IMPORTANT]
>
> ### How to get started
>
> Open the debugging notebook corresponding to the pipeline you want to inspect in Pluto:
>
> - [`00_debug_simulation.jl`](../notebooks/00_debug_simulation.jl) — simulation setup and generated data
> - [`01_debug_standard_decoding.jl`](../notebooks/01_debug_standard_decoding.jl) — standard decoding
> - [`02_debug_rerp_decoding.jl`](../notebooks/02_debug_rerp_decoding.jl) — rERP decoding
> - [`03_debug_plain_b2b.jl`](../notebooks/03_debug_plain_b2b.jl) — plain B2B
> - [`04_debug_one_step_b2b.jl`](../notebooks/04_debug_one_step_b2b.jl) — one-step B2B
> - [`05_debug_two_step_b2b.jl`](../notebooks/05_debug_two_step_b2b.jl) — two-step B2B
>
> Each notebook activates and instantiates the repository environment automatically.


> **Current structure:**
>```
> ├── notebooks/                                     # interactive debugging notebooks with sliders for different configurations
> │   ├── roamm                                      # debugging for ROAMM real data
> │   │   ├── roamm_sanity_check.jl
> │   │   ├── debug_roamm_two_step.jl
> │   │   ├── roamm.jl
> │   ├── 00_debug_simulation.jl 
> │   ├── 01_debug_standard_decoding.jl
> │   ├── 02_debug_rerp_decoding.jl
> │   ├── 03_debug_plain_b2b.jl           
> │   ├── 04_debug_one_step_b2b.jl
> │   ├── 05_debug_two_step_b2b.jl
> │   ├── 06_compare_pipelines.jl                    # overview of all five pipelines; no sliders
> │   ├── 07_correlated_continuous_sim_pipelines.jl
> │   └── simulation_controls.jl                     # shared PlutoUI sliders / simulation controls
>
> ├── plots/                                         # date based plots dirs
> ├── results/
> │   ├── condition_continuous_Sim/
> │   ├── correlated_continuous_Sim/
> │   ├── roamm/
>
> ├── ROAMM_preprocessing/
> │   ├── logs/
> │   ├── notebooks/
> │   │   ├── 00_inspect_roamm.ipynb
> │   │   ├── 01_check_roamm_preprocessing.ipynb
> │   ├── scripts/
> │   │   ├── 01_prepare_text.py
> │   │   ├── 02_compute_length_frequency.py
> │   │   ├── 03_compute_surprisal.py
> │   │   ├── 04_prepare_fixation_events.py
> │   │   ├── 05_merge_lexical_features.py
> │   │   ├── 06_export_all_eeg.py
>
>  
> ├── scripts/                                       # scripts for running and saving the five pipelines in folder `results/`
>
> ├── src/
> │   ├── pipelines                                 # implementations of the five pipelines
> │   │   ├── 01_standard_decoding.jl
> │   │   ├── 02_rerp_decoding.jl
> │   │   ├── 03_plain_b2b.jl
> │   │   ├── 04_one_step_b2b.jl
> │   │   └── 05_two_step_b2b.jl
> │   ├── plotting                                 # plotting functions
> │   │   ├── plot_b2b.jl
> │   │   └── plot_decoding.jl
> │   ├── roamm
> │   │   ├── two_step_b2b_roamm.jl
> │   ├── simulations                               # simulation setup and data generation for Cond_cont & 3 correlated continuous
> │   │   ├── ConditionContinuousSim.jl
> │   │   ├── CorrelatedContinuousSim.jl  
> │   ├── MScB2B.jl                                # main module of the custom MScB2B package
>```







