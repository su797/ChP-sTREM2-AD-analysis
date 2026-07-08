source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)

group_var <- project_config$variables$group_label_var

continuous_vars <- unique(c(
  project_config$group_comparison$target_vars,
  project_config$group_comparison$balance_continuous
))

continuous_results <- lapply(continuous_vars, function(var_name) {
  run_continuous_group_comparison(
    data = analysis_data,
    variable = var_name,
    group_var = group_var,
    transformation_table = transformation_plan,
    alpha = project_config$normality$alpha
  )
})

continuous_overall <- do.call(rbind, lapply(continuous_results, `[[`, "overall"))
continuous_summary <- do.call(rbind, lapply(continuous_results, `[[`, "summary"))
continuous_pairwise <- do.call(rbind, lapply(continuous_results, `[[`, "pairwise"))

categorical_vars <- project_config$group_comparison$balance_categorical
categorical_results <- lapply(categorical_vars, function(var_name) {
  run_categorical_group_comparison(
    data = analysis_data,
    variable = var_name,
    group_var = group_var
  )
})

categorical_overall <- do.call(rbind, lapply(categorical_results, `[[`, "overall"))
categorical_summary <- do.call(rbind, lapply(categorical_results, `[[`, "summary"))
categorical_pairwise <- do.call(rbind, lapply(categorical_results, `[[`, "pairwise"))

adjusted_models <- do.call(
  rbind,
  lapply(project_config$group_comparison$target_vars, function(var_name) {
    run_adjusted_linear_group_model(
      data = analysis_data,
      outcome = var_name,
      group_var = group_var,
      covariates = project_config$variables$covariates,
      factor_vars = project_config$variables$categorical_covariates,
      transformation_table = transformation_plan
    )
  })
)

group_overall_path <- file.path(result_summary_dir, "group_comparisons_overall.csv")
group_summary_path <- file.path(result_summary_dir, "group_comparisons_summary_by_group.csv")
group_pairwise_path <- file.path(result_summary_dir, "group_comparisons_pairwise.csv")
group_adjusted_path <- file.path(result_summary_dir, "group_comparisons_adjusted.csv")
group_overall_html <- file.path(result_tables_dir, "group_comparisons_overall.html")

group_overall <- rbind(continuous_overall, categorical_overall)
continuous_summary$level <- NA_character_
continuous_summary$pct <- NA_real_
continuous_summary <- continuous_summary[, c("variable", "analysis_var", "group", "level", "n", "pct", "summary")]

categorical_summary$analysis_var <- categorical_summary$variable
categorical_summary$summary <- NA_character_
categorical_summary <- categorical_summary[, c("variable", "analysis_var", "group", "level", "n", "pct", "summary")]

group_summary <- rbind(continuous_summary, categorical_summary)
group_pairwise <- rbind(continuous_pairwise, categorical_pairwise)

write_csv_utf8(group_overall, group_overall_path, row.names = FALSE)
write_csv_utf8(group_summary, group_summary_path, row.names = FALSE)
write_csv_utf8(group_pairwise, group_pairwise_path, row.names = FALSE)
write_csv_utf8(adjusted_models, group_adjusted_path, row.names = FALSE)

export_three_line_table(
  data = group_overall,
  csv_path = file.path(result_tables_dir, "group_comparisons_overall.csv"),
  html_path = group_overall_html,
  title = "Group Comparison Overview"
)

append_analysis_log(
  project_root = project_root,
  analysis_name = "03_group_comparisons",
  output_files = c(group_overall_path, group_summary_path, group_pairwise_path, group_adjusted_path),
  note = "Generated group comparisons for sTREM2, ChPICV, biomarkers, age, and sex; also exported adjusted group models into the current result version.",
  summary_dir = result_summary_dir
)
