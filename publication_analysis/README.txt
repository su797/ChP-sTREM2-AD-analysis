ChP/ICV and CSF biomarker analyses

This folder contains the incremental-analysis script and its saved results.
SEM and matching analyses are in the SEM and Matching subfolders, each with
its own README.

Run the incremental models with:

    Rscript incremental_explanatory_models.R OUTPUT_ROOT INPUT_CSV

Use an existing, writable directory for OUTPUT_ROOT. Results are written to
its analysis_outputs subfolder. Pass INPUT_CSV explicitly; the default path
in the script refers to the original local analysis directory.

The analysis was run in R 4.6.1 using base R and stats, without additional
packages. INPUT_CSV is the harmonized 735-participant analysis dataset,
not a raw ADNI download. The script checks the required columns, participant
IDs, and the 701 complete cases used for the structural models. These checks
help identify a different input file before the models are fitted.

For each of sTREM2, Tau, and P-Tau, the script compares a base model with a
model that also includes ChP/ICV. Biomarkers and volume ratios are on the
natural-log scale. The base model includes age, sex, education in years,
APOE epsilon4 status, and ventricular, hippocampal, and whole-brain volumes
divided by ICV. Sensitivity models add diagnosis, use absolute comparator
volumes with ICV, or use log-transformed education. Holm correction is applied
to the three primary model comparisons.

Saved coefficients, model formulas, and inclusion summaries are in
aggregate_results. The archive does not contain participant-level data or
fitted model objects holding individual observations. Access to ADNI data
requires ADNI data-use approval. The association-analysis scripts are at:
https://github.com/su797/ChP-sTREM2-AD-analysis
