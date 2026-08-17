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

```
├── notebooks/
│   ├── debug_one_step_b2b.jl 
│   ├── debug_plain_b2b.jl
│   ├── debug_rerp_decoding.jl
│   ├── debug_simulation.jl                
│   ├── debug_standard_decoding.jl
│   └──  simulation_controls.jl  #PlutoUI. Slider

├── src/
│   ├── pipelines
│   │   ├── one_step_b2b.jl
│   │   ├── plain_b2b.jl
│   │   ├── rerp_decoding.jl
│   │   ├── standard_decoding.jl
│   │   └── two_step_b2b.jl
│   ├── plotting
│   │   ├── plot_b2b.jl
│   │   └── plot_decoding.jl
│   ├── MScB2B.jl
│   ├── Pipelines.jl
│   ├── plotting.jl
│   └── simulation.jl
```

## 17 Aug 2026 — Plain B2B debugging

### Progress
- Investigated an unexpected condition-B2B peak around ~200 ms.
- Checked cross-validation repetitions, channel count, and predictor coding; the peak remained.
- Removing the shared P300 intercept (`β₀,P300 = 0`) removed the peak while keeping the N170 intercept unchanged.



## Current status

### Working
- [x] Simulation: `clean` / `overlap` / `confound` / `both`
- [x] Standard ridge decoding
- [x] rERP decoding
- [x] Plain B2B
- [x] One-step FIR+B2B
- [x] Plotting and simulation controls

### In progress
- [ ] One-step FIR+B2B debug
- [ ] Finalise two-step rERP → reconstructed single trials → B2B
- [ ] Compare the decoding/B2B pipelines
- [ ] Run repeated seeds / parameter sweeps

### Next

**🚩 Implement the two-step overlap-corrected B2B pipeline:**

`continuous EEG → rERP overlap correction → reconstructed single trials → B2B`

** ❗Main open question:**  
How should the rERP model be used to reconstruct corrected single-trial epochs that can then be passed into the existing plain B2B pipeline?