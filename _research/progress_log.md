# Progress Log

## 2026-05-12 (week 1)

### Added / Updated
- Created `simulation_design/` for current simulation prototypes.
- Implemented a collinearity design prototype with `CollinearCovariateDesign`.
- Added a `collinearity_strength` control.
- Added sanity checks showing that `face` trials have larger continuous covariate values under high collinearity.
- Added N170/P300 simulation checks:
  - `condition → N170`
  - `continuous → P300`
  

- Started refactoring the exploratory Pluto prototype into a more reproducible simulation workflow.
- Began moving reusable parts of the simulation code into `src`.
- Added initial source files for:
  - `EventDesign.jl`: condition-continuous design and collinearity manipulation.
  - `OnsetModels.jl`: onset model for manipulating temporal overlap.
  - `Components.jl`: simulated ERP-like components.
  - `B2BSim.jl`: main entry point for loading the simulation code.
- Added `scripts/run_single_case.jl` as a minimal sanity-check script.



### Meeting
- ❓: B2B has no sign? flip the generator to negative to see whether the b2b estimate coefficents turn negative also.
- ❗: refer to 'b2b_sign_interpretation/' : My current interpretation is that B2B estimates are numerically signed, so I would not simply say that B2B has no sign. 
 

 ### Notes
- The current proposal demo remains in `proposal_demo/` and is treated as exploratory code.
- The main goal so far is to separate reusable simulation logic from interactive Pluto visualization.



## 2026-05-19 (week 2)

### Added / Updated
- reproducible simulation workflow (WIP)
- overview
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

### Meeting
- Skipped, stella;s sick



## 2026-05-26 (week 3)

### Added / Updated
- Continued working on script-based simulation framework.
- Worked on B2B sign interpretation and related sanity checks.


### Meeting 
- Skipped, supervisor's vacation.
- short meeting with Jevri.

### Notes
- ❗: Current interpretation: B2B estimates mainly reflect predictor recovery strength / decodability rather than signed ERP-like waveforms.
- ❓: This motivates using absolute estimates, target-window recovery, and cross-talk / leakage for evaluation?




## 2026-06-02 (week 4)

### Added / Updated
- Building and refactoring reusable source files and scripts.
- Applied single-trial reconstruction method into b2b.
- Planned the first small simulation grid:
- Planned first-pass outputs:
  - coefficient table
  - evaluation summary
  - recovery and leakage plots.

### Meeting 
- Skipped, supervisor vacation

### Notes 
- The main focus is still to make the simulation workflow reproducible outside Pluto.
- Next Step: connect the pipeline outputs to the evaluation code.


## 2026-06-09 (week 5)

### Added / Updated
- Continued building and optimize the larger simulation/evaluation framework. 

### Meeting
- Skipped, Group retreat

### Notes
- ❗: Implementation time was limited, not much was done this week.
- Priority for next week: prepare the evaluation pipeline and run the grid. 





## 2026-06-16 (week 6)

### Added

### Meeting

- Skipped, Supervisor retreat 

### Notes 
- To do: 
  - [ ] : Standardise coefficient output format across methods.
  - [ ] : Implement first-pass recovery metrics in `Metrics.jl`.
  - [ ] : Implement leakage / cross-talk metrics.
  - [ ] : Save coefficents tables and evaluation summaries as CSV.
  - [ ] : Make first-pass recovery and leakage plots.



