source(file.path(getwd(), "00_setup.R"))

selected_data <- read_project_data(selected_data_path)
analysis_data <- read_project_data(analysis_data_path)
raw_data_all <- read_project_data(raw_data_all_path)

bind_rows_fill <- function(items) {
  items <- items[vapply(items, function(x) !is.null(x) && nrow(x) > 0, FUN.VALUE = logical(1))]
  if (length(items) == 0) {
    return(data.frame())
  }
  all_names <- unique(unlist(lapply(items, names), use.names = FALSE))
  aligned <- lapply(items, function(df) {
    missing_names <- setdiff(all_names, names(df))
    if (length(missing_names) > 0) {
      for (nm in missing_names) {
        df[[nm]] <- NA
      }
    }
    df[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, aligned)
  rownames(out) <- NULL
  out
}

metadata_vars <- c(
  "S_ID", "PTID", "RID", "COLPROT", "ORIGPROT", "SITE", "VISCODE",
  "EXAMDATE", "EXAMDATE.1", "DX_bl", "DX", "T1WI", "ABETA_convert"
)
metadata_df <- raw_data_all[, intersect(metadata_vars, names(raw_data_all)), drop = FALSE]
metadata_df <- metadata_df[!duplicated(metadata_df[, c("S_ID", "PTID"), drop = FALSE]), , drop = FALSE]

selected_with_meta <- merge(
  selected_data,
  metadata_df,
  by = c("S_ID", "PTID"),
  all.x = TRUE,
  sort = FALSE
)

core_vars <- c(
  "ChPICV", "MSD_STREM2CORRECTED", "EstimatedTotalIntraCranialVol",
  "S_ABETA", "TAU", "PTAU", "MMSE", "MOCA"
)
core_vars_present <- intersect(core_vars, names(raw_data_all))
core_complete_idx <- rep(TRUE, nrow(raw_data_all))
for (nm in core_vars_present) {
  core_complete_idx <- core_complete_idx & !is.na(raw_data_all[[nm]])
}

sample_flow <- data.frame(
  source_dataset_rows = nrow(raw_data_all),
  source_unique_participants = length(unique(raw_data_all$PTID)),
  source_rows_complete_core_fields = sum(core_complete_idx),
  source_unique_participants_complete_core_fields = length(unique(raw_data_all$PTID[core_complete_idx])),
  selected_rows = nrow(selected_data),
  selected_unique_participants = length(unique(selected_data$PTID)),
  duplicate_ptid_after_selection = sum(duplicated(selected_data$PTID)),
  duplicate_sid_after_selection = sum(duplicated(selected_data$S_ID)),
  analysis_rows = nrow(analysis_data),
  analysis_unique_participants = length(unique(analysis_data$PTID)),
  stringsAsFactors = FALSE
)

summarise_top_counts <- function(x, top_n = 10) {
  x <- x[!is.na(x) & nzchar(trimws(as.character(x)))]
  if (length(x) == 0) {
    return("")
  }
  tab <- sort(table(as.character(x)), decreasing = TRUE)
  tab <- tab[seq_len(min(length(tab), top_n))]
  paste(names(tab), as.integer(tab), collapse = "; ")
}

platform_summary <- data.frame(
  total_n = nrow(selected_with_meta),
  adni_phase_distribution = summarise_top_counts(selected_with_meta$ORIGPROT, top_n = 10),
  collection_protocol_distribution = summarise_top_counts(selected_with_meta$COLPROT, top_n = 10),
  visitcode_distribution = summarise_top_counts(selected_with_meta$VISCODE, top_n = 10),
  baseline_dx_distribution = summarise_top_counts(selected_with_meta$DX_bl, top_n = 10),
  current_dx_distribution = summarise_top_counts(selected_with_meta$DX, top_n = 10),
  t1_sequence_distribution = summarise_top_counts(selected_with_meta$T1WI, top_n = 20),
  top_site_distribution = summarise_top_counts(selected_with_meta$SITE, top_n = 10),
  stringsAsFactors = FALSE
)

safe_date <- function(x) {
  out <- as.Date(x)
  out
}

date_mri <- safe_date(selected_with_meta$EXAMDATE)
date_aux <- safe_date(selected_with_meta$EXAMDATE.1)
date_diff <- abs(as.numeric(date_aux - date_mri))
date_diff_non_missing <- date_diff[!is.na(date_diff)]
visit_alignment <- data.frame(
  n_total = nrow(selected_with_meta),
  n_with_both_dates = length(date_diff_non_missing),
  proportion_with_both_dates = round(length(date_diff_non_missing) / max(1, nrow(selected_with_meta)), 4),
  same_day_n = sum(date_diff_non_missing == 0),
  same_day_pct = round(mean(date_diff_non_missing == 0), 4),
  median_abs_days = if (length(date_diff_non_missing) > 0) stats::median(date_diff_non_missing) else NA_real_,
  iqr_abs_days = if (length(date_diff_non_missing) > 0) stats::IQR(date_diff_non_missing) else NA_real_,
  p90_abs_days = if (length(date_diff_non_missing) > 0) stats::quantile(date_diff_non_missing, probs = 0.9, names = FALSE) else NA_real_,
  max_abs_days = if (length(date_diff_non_missing) > 0) max(date_diff_non_missing) else NA_real_,
  stringsAsFactors = FALSE
)

missingness_vars <- c(
  "S_AGE", "S_PTGENDER", "PTEDUCAT", "APOE401",
  "ChPICV", "ChP_SUM", "MSD_STREM2CORRECTED", "S_ABETA", "TAU", "PTAU",
  "MMSE", "MOCA", "mPACCdigit", "mPACCtrailsB", "NegaADAS13", "ADAS13", "CDRSB"
)
missingness_vars <- intersect(missingness_vars, names(selected_data))
missingness_summary <- bind_rows_fill(lapply(missingness_vars, function(var_name) {
  non_missing_n <- sum(!is.na(selected_data[[var_name]]))
  data.frame(
    variable = var_name,
    n_total = nrow(selected_data),
    non_missing_n = non_missing_n,
    missing_n = nrow(selected_data) - non_missing_n,
    missing_pct = round((nrow(selected_data) - non_missing_n) / max(1, nrow(selected_data)) * 100, 2),
    stringsAsFactors = FALSE
  )
}))

regression_files <- list(
  chp_strem2 = file.path(result_summary_dir, "chp_strem2_linear_overall.csv"),
  tau = file.path(result_summary_dir, "tau_linear_overall.csv"),
  ptau = file.path(result_summary_dir, "ptau_linear_overall.csv"),
  advanced = file.path(result_summary_dir, "advanced_biomarker_adjusted_overall.csv")
)

read_if_exists <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  read_project_data(path)
}

regression_n_rows <- bind_rows_fill(list(
  local({
    df <- read_if_exists(regression_files$chp_strem2)
    if (nrow(df) == 0) return(data.frame())
    df <- df[df$term == "exposure", , drop = FALSE]
    data.frame(
      analysis_block = "Primary adjusted regression",
      model_name = df$model_name,
      n = df$n,
      estimation = "Complete-case linear regression",
      stringsAsFactors = FALSE
    )
  }),
  local({
    df <- read_if_exists(regression_files$tau)
    if (nrow(df) == 0) return(data.frame())
    df <- df[df$term == "exposure", , drop = FALSE]
    data.frame(
      analysis_block = "Primary adjusted regression",
      model_name = df$model_name,
      n = df$n,
      estimation = "Complete-case linear regression",
      stringsAsFactors = FALSE
    )
  }),
  local({
    df <- read_if_exists(regression_files$ptau)
    if (nrow(df) == 0) return(data.frame())
    df <- df[df$term == "exposure", , drop = FALSE]
    data.frame(
      analysis_block = "Primary adjusted regression",
      model_name = df$model_name,
      n = df$n,
      estimation = "Complete-case linear regression",
      stringsAsFactors = FALSE
    )
  })
))

sem_fit <- read_if_exists(file.path(result_summary_dir, "sem_model_fit.csv"))
sem_n_rows <- data.frame()
if (nrow(sem_fit) > 0) {
  sem_fit <- sem_fit[sem_fit$group %in% c("Overall", "CN", "MCI", "AD"), , drop = FALSE]
  sem_fit <- sem_fit[!duplicated(sem_fit[, c("group")]), c("group", "n"), drop = FALSE]
  sem_fit$analysis_block <- "Formal SEM"
  sem_fit$model_name <- paste0("MMSE+MoCA SEM (", sem_fit$group, ")")
  sem_fit$estimation <- "FIML latent SEM"
  sem_n_rows <- sem_fit[, c("analysis_block", "model_name", "n", "estimation")]
}

fiml_core_vars <- c("ChPICV", "MSD_STREM2CORRECTED", "S_ABETA", "TAU", "PTAU", "MMSE", "MOCA", "S_PTGENDER", "S_AGE", "PTEDUCAT", "APOE401")
fiml_core_vars <- intersect(fiml_core_vars, names(selected_data))
sem_complete_cases <- stats::complete.cases(selected_data[, fiml_core_vars, drop = FALSE])
fiml_boundary <- bind_rows_fill(list(
  data.frame(
    analysis_block = "SEM variable set",
    model_name = "Complete cases if listwise deletion were used",
    n = sum(sem_complete_cases),
    estimation = "Hypothetical listwise count",
    stringsAsFactors = FALSE
  ),
  sem_n_rows,
  regression_n_rows
))

sample_flow_path <- file.path(result_summary_dir, "sample_flow_summary.csv")
platform_summary_path <- file.path(result_summary_dir, "sample_platform_summary.csv")
visit_alignment_path <- file.path(result_summary_dir, "visit_alignment_summary.csv")
missingness_path <- file.path(result_summary_dir, "sample_missingness_summary.csv")
fiml_boundary_path <- file.path(result_summary_dir, "analysis_n_summary.csv")

write_csv_utf8(sample_flow, sample_flow_path, row.names = FALSE)
write_csv_utf8(platform_summary, platform_summary_path, row.names = FALSE)
write_csv_utf8(visit_alignment, visit_alignment_path, row.names = FALSE)
write_csv_utf8(missingness_summary, missingness_path, row.names = FALSE)
write_csv_utf8(fiml_boundary, fiml_boundary_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "02c_sample_flow_missingness_alignment",
  output_files = c(
    sample_flow_path,
    platform_summary_path,
    visit_alignment_path,
    missingness_path,
    fiml_boundary_path
  ),
  note = "Added sample-flow, phase/site/sequence provenance, visit-alignment interval, variable missingness, and analysis-level sample-size/FIML summary tables.",
  summary_dir = result_summary_dir
)
