# **MSc-Thesis:** Combining overlap correction and back-to-back regression methods to analyse eye-tracking-EEG co-registration data during natural reading

**Author:** *Yanhong Xu*

**Supervisor(s):** *Supervisor 1*, *Supervisor 2*

**Year:** *2026*



> [!IMPORTANT]
> Current WIP scripts and simulation prototypes are mainly stored in [`_research/`](./_research/).
> For recent updates, see [`progress_log.md`](./_research/progress_log.md).  




## Project Description
### Research Questions
1. In simulations mimicking co-registration data with known ground truth, does combining overlap correction (deconvolution rERP) with B2B improve the disentanglement and recovery of collinear predictors when fixation-related responses overlap in time?

2. When applied to ZuCo, does the best-performing pipeline yield stable and interpretable estimates for lexical predictors (surprisal, frequency, length) while controlling for oculomotor covariates?

### Main goals
1. A simulation framework with controllable temporal overlap and predictor collinearity is implemented using Unfold toolbox.

2. Analysis pipelines spanning the overlap × collinearity conceptual method map, including baseline decoding, deconvolution rERP decoding, B2B, and rERP + B2B variants, are systematically benchmarked against known ground truth.

3. A fixation-aware surprisal predictor is implemented, distinguishing first fixations and refixations.

4. The selected best-performing pipeline(s) will be applied to the ZuCo dataset, including lexical predictors such as word surprisal (from a pretrained language model), log word frequency, and word length, while controlling for oculomotor covariates.

## Zotero Library Path
>Please provide the link to the Zotero group here or include a `Bib`-File in the `report` folder

## Instruction for a new student
>If a fellow student wants to reproduce all your results. What scripts, in which order, with which data need to be run?
>
>Be as specific as possible. Plan to spend **at least 1h** on this.
>
>Optional: Add a pipeline plot in which the different steps are displayed together with the corresponding scripts.

## Overview of Folder Structure 

```
│projectdir          <- Project's main folder. It is initialized as a Git
│                       repository with a reasonable .gitignore file.
│
├── report           <- **Immutable and add-only!**
│   ├── proposal     <- Proposal PDF
│   ├── thesis       <- Final Thesis PDF
│   ├── talks        <- PDFs (and optionally pptx etc) of the Intro,
|   |                   Midterm & Final-Talk
|
├── _research        <- WIP scripts, code, notes, comments,
│   |                   to-dos and anything in an alpha state.
│
├── plots            <- All exported plots go here, best in date folders.
|   |                   Note that to ensure reproducibility it is required that all plots can be
|   |                   recreated using the plotting scripts in the scripts folder.
|
├── notebooks        <- Pluto, Jupyter, Weave or any other mixed media notebooks.*
│
├── scripts          <- Various scripts, e.g. simulations, plotting, analysis,
│   │                   The scripts use the `src` folder for their base code.
│
├── src              <- Source code for use in this project. Contains functions,
│                       structures and modules that are used throughout
│                       the project and in multiple scripts.
│
├── test             <- Folder containing tests for `src`.
│   └── runtests.jl  <- Main test file
│   └── setup.jl     <- Setup test environment
│
├── README.md        <- Top-level README. A fellow student needs to be able to
|   |                   continue your project. Think about her!!
|
├── .gitignore       <- focused on Julia, but some Matlab things as well
│
├── (Manifest.toml)  <- Contains full list of exact package versions used currently.
|── (Project.toml)   <- Main project file, allows activation and installation.
└── (Requirements.txt)<- in case of python project - can also be an anaconda file, MakeFile etc.
                        
```

\*Instead of having a separate *notebooks* folder, you can also delete it and integrate your notebooks in the scripts folder. However, notebooks should always be marked by adding `nb_` in front of the file name.
