# SEM analyses

Run these scripts from the `SEM` directory with lavaan 0.6-21 installed:

```sh
Rscript reproduce_sem.R /path/ChpSTREM2AD_analysis_dataset.csv /path/output
```

This fits the models in `models` and writes point estimates, fit indices,
and session information. To run the 1000 bootstrap repetitions as well, add
`bootstrap` at the end of the command. That step can take several hours.
`bootstrap_fit.R` saves progress every 50 repetitions and reuses existing
output, so use a new output directory if the data or models have changed.

The input is the harmonized ADNI analysis dataset, which is not included.
It requires authorized ADNI access and the same preprocessing used for the
paper. In particular, the `L_` columns are already log-transformed; a different
assay field is not an interchangeable input. The script checks the required
columns and the 735-row sample.

## Model setup

Cognition is measured by MMSE, MoCA, and reversed ADAS13. These start as raw
scores and are standardized across the full cohort, not separately within
each diagnosis. Indicator residual covariances are not fitted. Point estimates
use MLR with FIML for missing data. The latent-variable analyses retain 735
participants; the complete-case observed-score sensitivity analysis has 718.

The multigroup models constrain factor loadings to be equal across CN, MCI,
and dementia. The cognition factor's residual variance is fixed to 1 in the
overall models and in the CN reference group, and is estimated in MCI and
dementia. The raw-score models are the ones reported in the paper. Other
score definitions, separate-group fits, and partial-MoCA-loading models in
the diagnostic tables are sensitivity analyses. No latent-mean comparisons
are reported. The factor specification was revised during manuscript
revision, rather than specified at study inception.

## Reading the results

`aggregate_results` holds the estimates and diagnostic summaries;
`bootstrap_replicates` holds the estimates from individual resamples.
These are model-level outputs, not participant records.

Bootstrap fits use ML; point estimates use MLR. The model specification,
transformations, FIML treatment, and factor scaling otherwise stay the same.
Reported MLR P values come from delta-method standard errors, not from the
bootstrap. Percentile intervals use the admissible bootstrap fits. Failed
fits, negative variances, and other inadmissible solutions are recorded but
are not used in those intervals. All overall models had 1000 valid fits;
multigroup counts are listed in `aggregate_results/bootstrap_status.csv`.

The q values belong to different test sets:

- Overall single-path results: eight indirect effects.
- `group_omnibus_tests.csv`: 16 tests, comprising eight single-path and eight
  parallel-effect tests.
- `raw_metric_group_wald.csv`: a separate diagnostic set of eight full-metric
  and eight partial-metric tests.

Use the q column for the relevant analysis; the three sets are not substitutes
for one another.

`serial_crossmodel_ratio_diagnostic.csv` compares estimates from separate
serial and single-path models. The paper uses the unstandardized
(`scale=est`) paired-bootstrap ratio. Names such as `proportion_mediated` and
`serial_share` in the syntax are parameter labels, not claims about causal
proportions. These are cross-sectional models. Fit statistics and bootstrap
stability also do not validate a clinical cognition scale.
