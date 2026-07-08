# Run the public analysis pipeline.
#
# This runner is intentionally boring in the best possible way: it sources the
# statistical modules in order and writes CSV/figure outputs under results/.
# Submission-document preparation lives outside this repository.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  repo_root <- dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

script_dir <- file.path(repo_root, "scripts")
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(script_dir)

Sys.unsetenv("CHP_STREM2_RUN_ID")
Sys.unsetenv("RESEARCHR_RESULT_VERSION")

pipeline <- c(
  "01_clean_project_data.R",
  "02_prepare_analysis_dataset.R",
  "02c_sample_flow_missingness_alignment.R",
  "02b_descriptive_statistics.R",
  "03_group_comparisons.R",
  "04_balance_and_matching.R",
  "05_chp_strem2_analysis.R",
  "06_tau_analysis.R",
  "06_ptau_analysis.R",
  "06b_advanced_biomarker_models.R",
  "06c_abeta_truncation_sensitivity.R",
  "06d_robustness_extensions.R",
  "06e_overlap_weighting_sensitivity.R",
  "06f_diagnosis_biologic_sensitivity.R",
  "06g_phase_site_protocol_sensitivity.R",
  "06h_timewindow_lv_sensitivity.R",
  "07_sem_mediation.R",
  "07b_serial_sem_followup.R",
  "07c_measurement_invariance.R",
  "07d_bootstrap_key_indirects.R"
)

for (script in pipeline) {
  message("Running ", script)
  source(file.path(script_dir, script))
}

manual_qc_file <- file.path(repo_root, "data", "raw", "manual_chp_qc_integrated_80.csv")
if (file.exists(manual_qc_file)) {
  message("Running 08_manual_chp_qc_subset.R")
  source(file.path(script_dir, "08_manual_chp_qc_subset.R"))
} else {
  message("Skipping manual ChP QC subset: data/raw/manual_chp_qc_integrated_80.csv was not found.")
}

source(file.path(script_dir, "00_setup.R"))
writeLines(capture.output(sessionInfo()), con = file.path(result_run_dir, "sessionInfo.txt"), useBytes = TRUE)
message("Done. Results are in: ", result_run_dir)
