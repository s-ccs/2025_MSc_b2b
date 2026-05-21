# Progress Log

## 2026-05-12

### Added
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


### Notes
- The current proposal demo remains in `proposal_demo/` and is treated as exploratory code.
- The main goal so far is to separate reusable simulation logic from interactive Pluto visualization.


### Meeting
- ❓: B2B has no sign? flip the generator to negative to see whether the b2b estimate coefficents turn negative also.
- refer to 'b2b_sign_interpretation/' : My current interpretation is that B2B estimates are numerically signed, so I would not simply say that B2B has no sign. 
 

