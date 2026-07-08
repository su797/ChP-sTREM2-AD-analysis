source(file.path(getwd(), "00_setup.R"))

selected_data <- read_project_data(selected_data_path)

transformation_plan <- build_transformation_plan(
  data = selected_data,
  continuous_vars = project_config$variables$continuous_screen,
  alpha = project_config$normality$alpha,
  prefix = project_config$normality$prefix
)

transformation_result <- apply_transformation_plan(selected_data, transformation_plan)
analysis_data <- transformation_result$data
transformation_plan <- transformation_result$transformation_table

cognition_map_list <- list()
for (model_name in names(project_config$cognition_models)) {
  model_config <- project_config$cognition_models[[model_name]]
  composite_result <- create_cognition_composite(
    data = analysis_data,
    composite_name = model_name,
    components = model_config$components,
    directions = model_config$directions,
    transformation_table = transformation_plan,
    min_non_missing = model_config$min_non_missing
  )
  analysis_data <- composite_result$data
  cognition_map_list[[model_name]] <- composite_result$component_map
}

cognition_map <- do.call(rbind, cognition_map_list)
rownames(cognition_map) <- NULL

write_csv_utf8(transformation_plan, normality_table_path, row.names = FALSE)
write_csv_utf8(transformation_plan, transformation_table_path, row.names = FALSE)
write_csv_utf8(cognition_map, cognition_map_path, row.names = FALSE)
write_csv_utf8(analysis_data, analysis_data_path, row.names = FALSE)
latest_outputs <- character(0)
if (isTRUE(project_config$results$write_latest_clean_copy)) {
  write_csv_utf8(analysis_data, analysis_data_latest_path, row.names = FALSE)
  latest_outputs <- analysis_data_latest_path
}

append_analysis_log(
  project_root = project_root,
  analysis_name = "02_prepare_analysis_dataset",
  output_files = c(normality_table_path, transformation_table_path, cognition_map_path, analysis_data_path, latest_outputs),
  note = "Screened normality, created transformed analysis variables, built cognition composite scores, and archived the analysis dataset into the current result version.",
  summary_dir = result_summary_dir
)
