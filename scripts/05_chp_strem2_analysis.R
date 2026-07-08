source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
group_var <- project_config$variables$group_label_var

linear_overall_list <- list()
linear_by_group_list <- list()
partial_overall_list <- list()
partial_by_group_list <- list()
plot_files <- c()

for (pair_cfg in project_config$chp_strem2$linear_pairs) {
  fit_overall <- fit_linear_model_base(
    data = analysis_data,
    outcome = pair_cfg$outcome,
    exposure = pair_cfg$exposure,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  fit_overall$model_name <- pair_cfg$name
  linear_overall_list[[pair_cfg$name]] <- fit_overall

  fit_group <- fit_linear_models_by_group(
    data = analysis_data,
    group_var = group_var,
    outcome = pair_cfg$outcome,
    exposure = pair_cfg$exposure,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  fit_group$model_name <- pair_cfg$name
  linear_by_group_list[[pair_cfg$name]] <- fit_group

  x_plot <- resolve_analysis_var(pair_cfg$exposure, transformation_plan)
  y_plot <- resolve_analysis_var(pair_cfg$outcome, transformation_plan)

  overall_plot_path <- file.path(result_figures_dir, paste0("scatter_", pair_cfg$name, "_overall.png"))
  save_scatter_plot_base(
    data = analysis_data,
    x = x_plot,
    y = y_plot,
    path = overall_plot_path,
    group_var = group_var,
    title = paste(pair_cfg$exposure, "vs", pair_cfg$outcome, "(Overall)")
  )
  plot_files <- c(plot_files, overall_plot_path)

  for (group_name in unique(as.character(analysis_data[[group_var]]))) {
    sub_data <- analysis_data[as.character(analysis_data[[group_var]]) == group_name, , drop = FALSE]
    group_plot_path <- file.path(result_figures_dir, paste0("scatter_", pair_cfg$name, "_", group_name, ".png"))
    save_scatter_plot_base(
      data = sub_data,
      x = x_plot,
      y = y_plot,
      path = group_plot_path,
      title = paste(pair_cfg$exposure, "vs", pair_cfg$outcome, "(", group_name, ")", sep = "")
    )
    plot_files <- c(plot_files, group_plot_path)
  }
}

for (partial_pair in project_config$chp_strem2$partial_pairs) {
  pair_name <- paste(partial_pair, collapse = "_with_")
  partial_overall <- run_partial_correlation(
    data = analysis_data,
    x_var = partial_pair[[1]],
    y_var = partial_pair[[2]],
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  partial_overall$model_name <- pair_name
  partial_overall_list[[pair_name]] <- partial_overall

  partial_group <- run_partial_correlation_by_group(
    data = analysis_data,
    group_var = group_var,
    x_var = partial_pair[[1]],
    y_var = partial_pair[[2]],
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  partial_group$model_name <- pair_name
  partial_by_group_list[[pair_name]] <- partial_group
}

linear_overall <- do.call(rbind, linear_overall_list)
linear_by_group <- do.call(rbind, linear_by_group_list)
partial_overall <- do.call(rbind, partial_overall_list)
partial_by_group <- do.call(rbind, partial_by_group_list)

linear_primary <- linear_overall[linear_overall$term == "exposure", , drop = FALSE]
linear_primary <- prepare_regression_table(linear_primary)

overall_path <- file.path(result_summary_dir, "chp_strem2_linear_overall.csv")
by_group_path <- file.path(result_summary_dir, "chp_strem2_linear_by_group.csv")
partial_overall_path <- file.path(result_summary_dir, "chp_strem2_partial_overall.csv")
partial_by_group_path <- file.path(result_summary_dir, "chp_strem2_partial_by_group.csv")
table_csv <- file.path(result_tables_dir, "table_chp_strem2_linear_primary.csv")
table_html <- file.path(result_tables_dir, "table_chp_strem2_linear_primary.html")

write_csv_utf8(linear_overall, overall_path, row.names = FALSE)
write_csv_utf8(linear_by_group, by_group_path, row.names = FALSE)
write_csv_utf8(partial_overall, partial_overall_path, row.names = FALSE)
write_csv_utf8(partial_by_group, partial_by_group_path, row.names = FALSE)

export_three_line_table(
  data = linear_primary,
  csv_path = table_csv,
  html_path = table_html,
  title = "ChP and sTREM2 Linear Regression"
)

append_analysis_log(
  project_root = project_root,
  analysis_name = "05_chp_strem2_analysis",
  output_files = c(overall_path, by_group_path, partial_overall_path, partial_by_group_path, table_csv, table_html, plot_files),
  note = "Completed adjusted linear models, partial correlations, and scatter plots for ChPICV and sTREM2 in the current result version.",
  summary_dir = result_summary_dir
)
