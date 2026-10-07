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

## Progress
- Checked the original data, the code used a more sufficient way to define a ridge function, and also applied SVD decompose matrix.
- Due to the time limit, consider just evaluate S as main analysis, delta R as complementary analysis.
- Limitation found: a fixation contaion multiple information of multiple words.

- Why use GPT2:
  - The goal is not to maximize next-word prediction accuracy, but to model human language processing. Previous work has shown that larger and better-performing language models do not necessarily provide a better fit to human reading behaviour, and GPT-2 has been shown to provide a good predictor of human reading performance. It is also autoregressive, so its surprisal is based only on preceding context, which is appropriate for incremental reading.
  - Recalculate the word surprisal, because calculating surprisal naively and subtle tokenization issues can distort the results.

### Meeting notes
- run the model first, and then decide whether evaluate delta R or not. statistical test not decided yet.






## Current status

### ✅ Implemented

#### Simulation / pipelines
- [x] Simulation: `clean` / `overlap` / `confound` / `both`
- [x] Standard ridge decoding
- [x] rERP decoding
- [x] Plain B2B
- [x] One-step FIR + B2B
- [x] Two-step rERP → reconstructed single trials → B2B
- [x] Plotting and saving results for all five pipelines
- [x] Sampling-rate control for the default P300 spike
- [x] Three-correlated-predictor simulation:
  - Hanning-window components
  - continuous-1-biased overlap

#### ROAMM
- [x] Text → lexical predictors → fixation onset → EEG latency preprocessing
- [x] Initial EEG sanity checks
- [x] ROAMM data successfully used as input to two-step B2B
- [x] Initial subject-level B2B analysis

#### Methodological decisions / understanding
- [x] GPT-2 choice justified:
  - goal is modelling human language processing rather than maximizing
    next-token prediction accuracy
  - autoregressive prediction is appropriate for incremental reading
- [x] Identified potential problems with naive word-level surprisal calculation
  due to GPT-2 subword tokenization



## Current status
### 🔍 Open questions

#### Simulation / B2B
- [ ] Why does the Hanning-window simulation appear more decodable?
  - Is this important enough to investigate further?

- [ ] Understand `ShiftOnsetByOne`
  - Why do the overlap-related bumps appear where they do?
  - Which neighboring event contributes to each bump?
  - How does shifting the onset-distance assignment move/change the bumps?

- [ ] Why does one-step B2B separate the three correlated predictors less
      cleanly than plain/two-step B2B?
  - Current hypothesis: strong regularization in the continuous FIR model
  - Needs verification

- [ ] How should B2B recovery be numerically evaluated?
  - B2B estimates and the simulated ERP waveform do not necessarily have
    the same shape
  - quantify target recovery and cross-talk across pipelines

#### ROAMM
- [ ] What is the appropriate context for surprisal?
  - sentence-level
  - story-level / longer preceding context

- [ ] How should GPT-2 subword probabilities be converted into word
      probabilities?
  - account for tokenization / whitespace issues

- [ ] Should fixation events use left-eye, right-eye, or binocular events?

- [ ] EEG scaling:
  - BIDS metadata reports `µV`
  - synced `.pkl` / exported `.npy` values are around `1e-5`
  - determine where scaling occurred upstream

- [ ] Run EEG quality/sanity checks across all subjects and flag suspicious data

#### Evaluation / statistics
- [ ] Decide whether ΔR is needed as a complementary analysis to S
- [ ] Decide on statistical testing for subject-level/group-level S estimates
- [ ] Decide whether held-out evaluation is necessary
- [ ] If used, compare plain B2B vs. two-step B2B on held-out ROAMM data

#### Methods / code understanding
- [ ] Understand the feature-importance implementation
  - trace through the code
  - understand exactly what quantity is estimated, SVD
  - determine how it should be interpreted


### 🚩 Immediate priorities

1. **Recalculate word surprisal**
   - decide context definition
   - implement tokenization-aware word probability

2. **Run two-step B2B on ROAMM**
   - word length
   - word frequency
   - word surprisal

3. **Obtain and inspect subject-level S estimates**

4. **Run EEG sanity checks across subjects**

5. **Then decide on evaluation/statistics**
   - ΔR?
   - statistical test?





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



