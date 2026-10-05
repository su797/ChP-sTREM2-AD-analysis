Matching sensitivity analysis (Supplementary Figure S11)

Run with base R:

    Rscript reproduce_matching.R DATA_CLEAN_DIR TRANSFORMATION_PLAN_CSV OUTPUT_DIR

DATA_CLEAN_DIR should contain:
    ChpSTREM2AD_analysis_dataset.csv
    matched_CN_vs_AD.csv
    matched_CN_vs_MCI.csv
    matched_MCI_vs_AD.csv

Use the transformation_plan.csv from the same preprocessing run for the
second argument. Participant-level inputs are not included in this archive;
they require authorized ADNI access. The preprocessing and matching scripts
are in the repository cited in the paper.

This script reads the existing matched samples; it does not redo matching.
For each diagnostic pair, it fits the same regression in the matched sample
and in all participants from those two groups before matching. The pair is
read from the filename, and both groups must be present. For example, the
unmatched CN_vs_AD comparison includes CN and AD, not CN alone or all three
diagnostic groups. AD is displayed as Dementia in the figure.

The outcomes are natural-log sTREM2, P-Tau, and Tau. The predictor is
natural-log ChP/ICV, adjusted for age, sex, natural-log education, and APOE
epsilon4 status. Outputs include standardized coefficients and 95% confidence
intervals. P values are nominal, two-sided t-test P values.

Saved results are in aggregate_results:
    matched_key_model_sensitivity.csv   regression estimates
    sample_and_model_validation.csv    sample counts and checks
    S11_plot_values.csv                coefficients and intervals used in S11

To redraw S11 from the saved estimates:

    Rscript plot_s11.R MATCHING_DIR

MATCHING_DIR is the directory containing aggregate_results. The plotting
script writes PNG and PDF files under its configured output subfolder
(see the out variable in plot_s11.R), and writes S11_plot_values.csv back
to aggregate_results. It does not refit the models.
