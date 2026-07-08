# Choroid plexus, sTREM2, and cognition across the AD continuum

This repository contains the R analysis code used for the cross-sectional ADNI
analysis of choroid plexus volume, CSF sTREM2, tau-related biomarkers, A・趣ｽｲ, and
cognition.

The code here is intentionally limited to the analysis workflow. It does not
generate submission documents.

## What is included

The pipeline covers the analyses reported in the study:

- construction of the baseline analytic dataset from the local ADNI-derived files;
- descriptive statistics and diagnostic-stage group comparisons;
- primary covariate-adjusted regression models;
- partial-correlation and linear-model summaries for ChP/ICV, sTREM2, Tau, P-Tau, and A・趣ｽｲ;
- false-discovery-rate adjustment;
- pairwise matching and overlap-weighting robustness checks;
- A・趣ｽｲ ceiling/censoring sensitivity checks;
- alternative ChP definitions and structure-specificity models;
- diagnosis-adjusted, biologically defined, phase/site/protocol, time-window, and lateral-ventricle sensitivity analyses;
- latent-variable SEM/path-decomposition models and bootstrap confidence intervals;
- optional 80-participant manually corrected ChP segmentation subset checks.

## What is not included

ADNI participant-level data are not redistributed. Please obtain data through
the ADNI data-use process and place the required local files under `data/raw/`.

Submission materials, reviewer comments, and private project notes are
deliberately excluded from this public code release.

## Folder layout

```text
.
隨乗㈹讌ｳ隨渉 R/                  # helper functions shared by the analysis scripts
隨乗㈹讌ｳ隨渉 scripts/            # numbered analysis modules
隨乗㈹讌ｳ隨渉 config/             # light-weight analysis settings
隨乗㈹讌ｳ隨渉 data/
隨上・  隨乗㈹讌ｳ隨渉 README.md
隨上・  隨乗喚讌ｳ隨渉 raw/            # local ADNI-derived input files, not tracked by Git
隨乗㈹讌ｳ隨渉 results/            # timestamped outputs, not tracked by Git
隨乗喚讌ｳ隨渉 run_analysis.R      # pipeline runner
```

## Input files

Place these files in `data/raw/`:

- `Data.csv`
- `Data_all.csv`

Optional:

- `manual_chp_qc_integrated_80.csv`

The expected variables are listed in `scripts/project_config.R`. In brief, the
main workflow expects diagnostic group, age, sex, education, APOE ・趣ｽｵ4 status,
ChP/ICV, bilateral ChP volume, estimated intracranial volume, sTREM2, A・趣ｽｲ, Tau,
P-Tau, comparator brain volumes, and cognitive measures.

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

These diagnostic files are not article drafts; they are just a convenient
way to keep model checks near the CSV outputs.

## R packages

Core analyses use base R plus a small set of common packages:

- `readxl` for reading Excel files, if local inputs are provided as Excel files;
- `survival` for censored A・趣ｽｲ sensitivity models;
- `lavaan` for latent-variable SEM models.

The scripts stop with a clear message if a required package is missing.

## Reproducibility notes

The random seed is set in `scripts/project_config.R`. The runner also writes
`sessionInfo.txt` into each result folder so that package versions are easy to
check later.

The overlap-weighting module is used as a demographic rebalancing sensitivity
check. The primary models remain full-sample covariate-adjusted models.

The SEM modules are best read as cross-sectional path-decomposition analyses,
not as causal mediation or temporal-ordering tests.

## Code availability statement

Analysis code is available in this repository. ADNI data are available to
qualified investigators through ADNI and are not included here.
