source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
raw_data_all <- read_project_data(raw_data_all_path)
transformation_plan <- read_project_data(transformation_table_path)

metadata_vars <- c("S_ID", "PTID", "EXAMDATE", "EXAMDATE.1")
metadata_df <- raw_data_all[, intersect(metadata_vars, names(raw_data_all)), drop = FALSE]
metadata_df <- metadata_df[!duplicated(metadata_df[, c("S_ID", "PTID"), drop = FALSE]), , drop = FALSE]

analysis_with_meta <- merge(
  analysis_data,
  metadata_df,
  by = c("S_ID", "PTID"),
  all.x = TRUE,
  sort = FALSE
)

safe_date <- function(x) as.Date(x)
date_mri <- safe_date(analysis_with_meta$EXAMDATE)
date_aux <- safe_date(analysis_with_meta$EXAMDATE.1)
analysis_with_meta$abs_days <- abs(as.numeric(date_aux - date_mri))

core_models <- list(
  list(name = "ChPICV_on_sTREM2", outcome = "ChPICV", exposure = "MSD_STREM2CORRECTED"),
  list(name = "ChPICV_on_TAU", outcome = "ChPICV", exposure = "TAU"),
  list(name = "ChPICV_on_PTAU", outcome = "ChPICV", exposure = "PTAU"),
  list(name = "ChPICV_on_ABETA", outcome = "ChPICV", exposure = "S_ABETA")
)

extract_exposure_row <- function(fit_df) {
  fit_df[fit_df$term == "exposure", , drop = FALSE]
}

run_time_window_models <- function(data_with_meta, windows) {
  rows <- list()
  for (window in windows) {
    if (is.infinite(window)) {
      sub_data <- data_with_meta
      label <- "Full sample"
    } else {
      sub_data <- data_with_meta[!is.na(data_with_meta$abs_days) & data_with_meta$abs_days <= window, , drop = FALSE]
      label <- paste0("<= ", window, " days")
    }
    for (cfg in core_models) {
      fit_rows <- fit_linear_model_base(
        data = sub_data,
        outcome = cfg$outcome,
        exposure = cfg$exposure,
        covariates = project_config$variables$covariates,
        factor_vars = project_config$variables$categorical_covariates,
        transformation_table = transformation_plan
      )
      fit_rows <- extract_exposure_row(fit_rows)
      if (nrow(fit_rows) == 0) next
      fit_rows$window_label <- label
      fit_rows$window_days <- ifelse(is.infinite(window), NA, window)
      fit_rows$n_with_dates <- sum(!is.na(sub_data$abs_days))
      fit_rows$median_abs_days <- if (all(is.na(sub_data$abs_days))) NA_real_ else stats::median(sub_data$abs_days, na.rm = TRUE)
      rows[[paste(label, cfg$name, sep = "_")]] <- fit_rows
    }
  }
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

run_lv_adjusted_models <- function(data) {
  rows <- list()
  for (cfg in core_models) {
    base_rows <- fit_linear_model_base(
      data = data,
      outcome = cfg$outcome,
      exposure = cfg$exposure,
      covariates = project_config$variables$covariates,
      factor_vars = project_config$variables$categorical_covariates,
      transformation_table = transformation_plan
    )
    base_rows <- extract_exposure_row(base_rows)
    if (nrow(base_rows) > 0) {
      base_rows$model_spec <- "Primary adjusted"
      rows[[paste(cfg$name, "base", sep = "_")]] <- base_rows
    }

    lv_rows <- fit_linear_model_base(
      data = data,
      outcome = cfg$outcome,
      exposure = cfg$exposure,
      covariates = c(project_config$variables$covariates, "LV_SUM"),
      factor_vars = project_config$variables$categorical_covariates,
      transformation_table = transformation_plan
    )
    lv_rows <- extract_exposure_row(lv_rows)
    if (nrow(lv_rows) > 0) {
      lv_rows$model_spec <- "Primary adjusted + lateral ventricle volume"
      rows[[paste(cfg$name, "lv", sep = "_")]] <- lv_rows
    }
  }
  if (length(rows) == 0) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

time_window_results <- run_time_window_models(analysis_with_meta, windows = c(Inf, 30, 90))
lv_adjusted_results <- run_lv_adjusted_models(analysis_data)

time_window_summary <- data.frame(
  subset = c("Full sample", "<= 30 days", "<= 90 days"),
  n = c(
    nrow(analysis_with_meta),
    sum(!is.na(analysis_with_meta$abs_days) & analysis_with_meta$abs_days <= 30),
    sum(!is.na(analysis_with_meta$abs_days) & analysis_with_meta$abs_days <= 90)
  ),
  median_abs_days = c(
    stats::median(analysis_with_meta$abs_days, na.rm = TRUE),
    stats::median(analysis_with_meta$abs_days[!is.na(analysis_with_meta$abs_days) & analysis_with_meta$abs_days <= 30], na.rm = TRUE),
    stats::median(analysis_with_meta$abs_days[!is.na(analysis_with_meta$abs_days) & analysis_with_meta$abs_days <= 90], na.rm = TRUE)
  ),
  max_abs_days = c(
    max(analysis_with_meta$abs_days, na.rm = TRUE),
    max(analysis_with_meta$abs_days[!is.na(analysis_with_meta$abs_days) & analysis_with_meta$abs_days <= 30], na.rm = TRUE),
    max(analysis_with_meta$abs_days[!is.na(analysis_with_meta$abs_days) & analysis_with_meta$abs_days <= 90], na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

time_window_path <- file.path(result_summary_dir, "time_window_sensitivity_models.csv")
time_window_summary_path <- file.path(result_summary_dir, "time_window_sensitivity_summary.csv")
lv_adjusted_path <- file.path(result_summary_dir, "lv_adjusted_models.csv")

write_csv_utf8(time_window_results, time_window_path, row.names = FALSE)
write_csv_utf8(time_window_summary, time_window_summary_path, row.names = FALSE)
write_csv_utf8(lv_adjusted_results, lv_adjusted_path, row.names = FALSE)

if (exists("log_analysis_step", mode = "function")) {
  log_analysis_step(
    script = "06h_timewindow_lv_sensitivity.R",
    outputs = c(time_window_path, time_window_summary_path, lv_adjusted_path),
    notes = "Time-window restricted sensitivity and lateral-ventricle-adjusted ChP models."
  )
}
