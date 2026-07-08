source(file.path(getwd(), "00_setup.R"))

if (!file.exists(manual_qc_data_path)) {
  stop(
    "Manual ChP QC file not found. Put manual_chp_qc_integrated_80.csv in data/raw/ if you want to run this module.",
    call. = FALSE
  )
}

qc <- read_project_data(manual_qc_data_path)

manual_dir <- ensure_dir(file.path(result_summary_dir, "manual_chp_qc"))
manual_fig_dir <- ensure_dir(file.path(result_figures_dir, "manual_chp_qc"))

cor_row <- function(x, y, method) {
  ok <- stats::complete.cases(x, y)
  test <- stats::cor.test(x[ok], y[ok], method = method, exact = FALSE)
  data.frame(
    method = method,
    n = sum(ok),
    estimate = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

lin_concordance <- function(x, y) {
  ok <- stats::complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  2 * stats::cov(x, y) / (stats::var(x) + stats::var(y) + (mean(x) - mean(y))^2)
}

agreement_summary <- rbind(
  transform(
    cor_row(qc$manual_ChP_volume_mm3, qc$FreeSurfer_ChP_SUM_mm3, "pearson"),
    measure = "ChP volume"
  ),
  transform(
    cor_row(qc$manual_ChP_volume_mm3, qc$FreeSurfer_ChP_SUM_mm3, "spearman"),
    measure = "ChP volume"
  ),
  transform(
    cor_row(qc$manual_ChPICV, qc$FreeSurfer_ChPICV, "pearson"),
    measure = "ChP/ICV"
  ),
  transform(
    cor_row(qc$manual_ChPICV, qc$FreeSurfer_ChPICV, "spearman"),
    measure = "ChP/ICV"
  )
)

volume_diff <- qc$manual_ChP_volume_mm3 - qc$FreeSurfer_ChP_SUM_mm3
ratio_diff <- qc$manual_ChPICV - qc$FreeSurfer_ChPICV
ba_summary <- data.frame(
  measure = c("ChP volume", "ChP/ICV"),
  mean_difference = c(mean(volume_diff, na.rm = TRUE), mean(ratio_diff, na.rm = TRUE)),
  lower_limit = c(
    mean(volume_diff, na.rm = TRUE) - 1.96 * stats::sd(volume_diff, na.rm = TRUE),
    mean(ratio_diff, na.rm = TRUE) - 1.96 * stats::sd(ratio_diff, na.rm = TRUE)
  ),
  upper_limit = c(
    mean(volume_diff, na.rm = TRUE) + 1.96 * stats::sd(volume_diff, na.rm = TRUE),
    mean(ratio_diff, na.rm = TRUE) + 1.96 * stats::sd(ratio_diff, na.rm = TRUE)
  ),
  lin_concordance = c(
    lin_concordance(qc$manual_ChP_volume_mm3, qc$FreeSurfer_ChP_SUM_mm3),
    lin_concordance(qc$manual_ChPICV, qc$FreeSurfer_ChPICV)
  ),
  stringsAsFactors = FALSE
)

group_counts <- as.data.frame(table(qc$group_CN_MCI_Dementia), stringsAsFactors = FALSE)
names(group_counts) <- c("group", "n")

std_beta_lm <- function(data, outcome, exposure, covariates) {
  model_vars <- c(outcome, exposure, covariates)
  dd <- data[, model_vars, drop = FALSE]
  dd <- dd[stats::complete.cases(dd), , drop = FALSE]
  if (nrow(dd) < 10) {
    return(data.frame(n = nrow(dd), std_beta = NA_real_, p_value = NA_real_))
  }
  formula_text <- paste0("scale(", outcome, ") ~ scale(", exposure, ")")
  if (length(covariates) > 0) {
    formula_text <- paste(formula_text, paste(covariates, collapse = " + "), sep = " + ")
  }
  fit <- stats::lm(stats::as.formula(formula_text), data = dd)
  coefs <- summary(fit)$coefficients
  term <- paste0("scale(", exposure, ")")
  data.frame(
    n = nrow(dd),
    std_beta = unname(coefs[term, "Estimate"]),
    p_value = unname(coefs[term, "Pr(>|t|)"]),
    stringsAsFactors = FALSE
  )
}

# These names come from the manually corrected QC file. Keeping the mapping here
# makes the subset analysis easier to audit later.
covariates <- c("age", "female", "education_years", "APOE_e4_carrier")
biomarkers <- c(
  sTREM2 = "log_sTREM2",
  Tau = "log_Tau",
  P_Tau = "log_P_Tau",
  ABeta = "log_ABeta"
)

model_rows <- do.call(
  rbind,
  lapply(names(biomarkers), function(label) {
    exposure <- biomarkers[[label]]
    log_model <- std_beta_lm(qc, "manual_log_ChPICV", exposure, covariates)
    raw_model <- std_beta_lm(qc, "manual_ChPICV", exposure, covariates)
    rbind(
      data.frame(biomarker = label, outcome_scale = "log manual ChP/ICV", log_model),
      data.frame(biomarker = label, outcome_scale = "raw manual ChP/ICV", raw_model)
    )
  })
)
rownames(model_rows) <- NULL

write_csv_utf8(group_counts, file.path(manual_dir, "manual_chp_qc_group_counts.csv"), row.names = FALSE)
write_csv_utf8(agreement_summary, file.path(manual_dir, "manual_chp_agreement_correlations.csv"), row.names = FALSE)
write_csv_utf8(ba_summary, file.path(manual_dir, "manual_chp_bland_altman_summary.csv"), row.names = FALSE)
write_csv_utf8(model_rows, file.path(manual_dir, "manual_chp_replacement_models.csv"), row.names = FALSE)

scatter_path <- file.path(manual_fig_dir, "manual_vs_freesurfer_chpicv_scatter.png")
png(scatter_path, width = 1800, height = 1500, res = 220)
par(mar = c(5, 5, 2, 1), bty = "l")
plot(
  qc$FreeSurfer_ChPICV,
  qc$manual_ChPICV,
  pch = 16,
  col = "#1f5a85AA",
  xlab = "FreeSurfer ChP/ICV",
  ylab = "Manually corrected ChP/ICV"
)
abline(stats::lm(manual_ChPICV ~ FreeSurfer_ChPICV, data = qc), col = "#b31b34", lwd = 2)
abline(0, 1, col = "grey55", lty = 2)
dev.off()

ba_path <- file.path(manual_fig_dir, "manual_vs_freesurfer_chpicv_bland_altman.png")
ba_mean <- rowMeans(cbind(qc$manual_ChPICV, qc$FreeSurfer_ChPICV), na.rm = TRUE)
ba_diff <- qc$manual_ChPICV - qc$FreeSurfer_ChPICV
png(ba_path, width = 1800, height = 1500, res = 220)
par(mar = c(5, 5, 2, 1), bty = "l")
plot(
  ba_mean,
  ba_diff,
  pch = 16,
  col = "#1f5a85AA",
  xlab = "Mean ChP/ICV from both methods",
  ylab = "Manual - FreeSurfer ChP/ICV"
)
abline(h = mean(ba_diff, na.rm = TRUE), col = "#b31b34", lwd = 2)
abline(h = mean(ba_diff, na.rm = TRUE) + c(-1.96, 1.96) * stats::sd(ba_diff, na.rm = TRUE), col = "grey45", lty = 2)
dev.off()

append_analysis_log(
  project_root = repo_root,
  analysis_name = "08_manual_chp_qc_subset",
  output_files = c(
    file.path(manual_dir, "manual_chp_qc_group_counts.csv"),
    file.path(manual_dir, "manual_chp_agreement_correlations.csv"),
    file.path(manual_dir, "manual_chp_bland_altman_summary.csv"),
    file.path(manual_dir, "manual_chp_replacement_models.csv"),
    scatter_path,
    ba_path
  ),
  note = "Ran the manually corrected 80-participant ChP subset checks and exported agreement summaries, replacement models, and two simple QC figures.",
  summary_dir = result_summary_dir
)
