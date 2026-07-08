source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
group_var <- project_config$variables$group_label_var
desc_cfg <- runtime_settings$descriptive %||% list()
desc_digits <- desc_cfg$digits %||% 2

keep_existing_vars <- function(vars, data) {
  vars <- unlist(vars, use.names = FALSE)
  vars[vars %in% names(data)]
}

format_mean_sd <- function(mean_value, sd_value, digits = 2) {
  sprintf(paste0("%.", digits, "f +- %.", digits, "f"), mean_value, sd_value)
}

format_n_pct <- function(n_value, pct_value, digits = 1) {
  sprintf("%d (%.1f%%)", n_value, round(pct_value, digits))
}

overall_cont_vars <- keep_existing_vars(desc_cfg$continuous_overall, analysis_data)
overall_cat_vars <- keep_existing_vars(desc_cfg$categorical_overall, analysis_data)
group_cont_vars <- keep_existing_vars(desc_cfg$continuous_by_group, analysis_data)
group_cat_vars <- keep_existing_vars(desc_cfg$categorical_by_group, analysis_data)

overall_continuous <- if (length(overall_cont_vars) > 0) {
  out <- summarize_continuous(analysis_data, vars = overall_cont_vars, digits = desc_digits)
  out$summary_text <- mapply(format_mean_sd, out$mean, out$sd, MoreArgs = list(digits = desc_digits))
  out
} else {
  data.frame()
}

overall_categorical <- if (length(overall_cat_vars) > 0) {
  out <- summarize_categorical(analysis_data, vars = overall_cat_vars)
  out$summary_text <- mapply(format_n_pct, out$n, out$pct, MoreArgs = list(digits = 1))
  out
} else {
  data.frame()
}

group_continuous <- if (length(group_cont_vars) > 0) {
  out <- summarize_continuous(analysis_data, vars = group_cont_vars, group_var = group_var, digits = desc_digits)
  out$summary_text <- mapply(format_mean_sd, out$mean, out$sd, MoreArgs = list(digits = desc_digits))
  out
} else {
  data.frame()
}

group_categorical <- if (length(group_cat_vars) > 0) {
  out <- summarize_categorical(analysis_data, vars = group_cat_vars, group_var = group_var)
  out$summary_text <- mapply(format_n_pct, out$n, out$pct, MoreArgs = list(digits = 1))
  out
} else {
  data.frame()
}

overall_cont_path <- file.path(result_summary_dir, "descriptive_overall_continuous.csv")
overall_cat_path <- file.path(result_summary_dir, "descriptive_overall_categorical.csv")
group_cont_path <- file.path(result_summary_dir, "descriptive_by_group_continuous.csv")
group_cat_path <- file.path(result_summary_dir, "descriptive_by_group_categorical.csv")

write_csv_utf8(overall_continuous, overall_cont_path, row.names = FALSE)
write_csv_utf8(overall_categorical, overall_cat_path, row.names = FALSE)
write_csv_utf8(group_continuous, group_cont_path, row.names = FALSE)
write_csv_utf8(group_categorical, group_cat_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "02b_descriptive_statistics",
  output_files = c(overall_cont_path, overall_cat_path, group_cont_path, group_cat_path),
  note = "Generated overall and group-wise descriptive statistics for configured continuous and categorical variables.",
  summary_dir = result_summary_dir
)
