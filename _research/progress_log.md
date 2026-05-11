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

### Notes
- The current proposal demo remains in `proposal_demo/` and is treated as exploratory code.
- The current focus is to combine the new collinearity manipulation with the existing overlap simulation.
