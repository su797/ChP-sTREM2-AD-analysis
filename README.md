# Choroid plexus, sTREM2, and cognition across the AD continuum

This repository contains the R analysis code used for a cross-sectional ADNI
study of choroid plexus volume, CSF sTREM2, tau-related biomarkers,
amyloid-beta, and cognition.

The code is limited to the analysis workflow. It does not generate submission
documents.

## Analyses accompanying the current manuscript

For **Choroid plexus volume as a candidate imaging biomarker in Alzheimer's
disease**, the materials from Additional file 3 are in
[`publication_analysis/`](publication_analysis/README.txt):

- Incremental explanatory models and their saved results are in the top level
  of that directory.
- [`SEM/`](publication_analysis/SEM/README.md) contains the three-indicator
  cognition models, model syntax, and bootstrap results.
- [`Matching/`](publication_analysis/Matching/README.txt) contains the paired
  diagnostic-group comparisons used for Supplementary Figure S11.

Use these versions for the manuscript's incremental, SEM, and S11 analyses.
They are run separately from `run_analysis.R`; the older SEM and matching
modules in the general pipeline below do not replace them. Each README gives
the input requirements and commands. The saved tables contain model estimates
and group summaries, not participant-level data.

## What is included

The original pipeline covers:

- construction of the baseline analytic dataset from local ADNI-derived files;
- descriptive statistics and diagnostic-stage group comparisons;
- primary covariate-adjusted regression models;
- partial-correlation and linear-model summaries for ChP/ICV, sTREM2, Tau,
  P-Tau, and amyloid-beta;
- false-discovery-rate adjustment;
- pairwise matching and overlap-weighting robustness checks;
- amyloid-beta ceiling/censoring sensitivity checks;
- alternative ChP definitions and structure-specificity models;
- diagnosis-adjusted, biologically defined, phase/site/protocol, time-window,
  and lateral-ventricle sensitivity analyses;
- latent-variable SEM/path-decomposition models and bootstrap confidence
  intervals;
- optional 80-participant manually corrected ChP segmentation subset checks.

## What is not included

ADNI participant-level data are not redistributed. Please obtain data through
the ADNI data-use process and place the required local files under `data/raw/`.

Submission materials, reviewer comments, and private project notes are
deliberately excluded from this public code release.

## Folder layout

```text
.
|-- R/                  # helper functions shared by the analysis scripts
|-- scripts/            # numbered analysis modules
|-- config/             # lightweight analysis settings
|-- publication_analysis/ # manuscript-specific code and aggregate results
|-- data/
|   |-- README.md
|   `-- raw/            # local ADNI-derived input files, not tracked by Git
|-- results/            # timestamped outputs, not tracked by Git
`-- run_analysis.R      # pipeline runner
```

## Input files

Place these files in `data/raw/`:

- `Data.csv`
- `Data_all.csv`

Optional:

- `manual_chp_qc_integrated_80.csv`

The expected variables are listed in `scripts/project_config.R`. In brief, the
main workflow expects diagnostic group, age, sex, education, APOE epsilon 4
status, ChP/ICV, bilateral ChP volume, estimated intracranial volume, sTREM2,
amyloid-beta, Tau, P-Tau, comparator brain volumes, and cognitive measures.

## Running the analysis

From the repository root:

```r
source("run_analysis.R")
```

Or from a shell:

```bash
Rscript run_analysis.R
```

Each run creates a timestamped folder under `results/`. The most recent run is
also recorded in `results/LATEST.txt`.

## Main outputs

Most outputs are plain CSV files under:

```text
results/<run_id>/summary/
```

Cleaned analysis data are written to:

```text
results/<run_id>/data_clean/
```

Figures used for analysis checks are written to:

```text
results/<run_id>/figures/
```

SEM diagnostics and small text diagnostics, when produced, are written to:

```text
results/<run_id>/diagnostics/
```

These diagnostic files are not article drafts; they simply keep model checks
near the CSV outputs.

## R packages

Core analyses use base R plus a small set of common packages:

- `readxl` for reading Excel files, if local inputs are provided as Excel files;
- `survival` for censored amyloid-beta sensitivity models;
- `lavaan` for latent-variable SEM models.

The scripts stop with a clear message if a required package is missing.

## Reproducibility notes

The random seed is set in `scripts/project_config.R`. The runner also writes
`sessionInfo.txt` into each result folder so package versions can be checked
later.

The overlap-weighting module is used as a demographic rebalancing sensitivity
check. The primary models remain full-sample covariate-adjusted models.

The SEM modules should be read as cross-sectional path-decomposition analyses,
not as causal mediation or temporal-ordering tests.

## Code availability statement

Analysis code is available in this repository. ADNI data are available to
qualified investigators through ADNI and are not included here.
