source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
group_overall <- read_project_data(file.path(result_summary_dir, "group_comparisons_overall.csv"))

significant_continuous <- group_overall$variable[
  group_overall$variable %in% project_config$group_comparison$balance_continuous &
    group_overall$p_value < project_config$matching$alpha
]

significant_categorical <- group_overall$variable[
  group_overall$variable %in% project_config$group_comparison$balance_categorical &
    group_overall$p_value < project_config$matching$alpha
]

matching_required <- length(significant_continuous) > 0 || length(significant_categorical) > 0
matching_decision <- data.frame(
  matching_required = matching_required,
  significant_continuous = paste(significant_continuous, collapse = "; "),
  significant_categorical = paste(significant_categorical, collapse = "; "),
  stringsAsFactors = FALSE
)

decision_path <- file.path(result_summary_dir, "matching_decision.csv")
write_csv_utf8(matching_decision, decision_path, row.names = FALSE)

if (!matching_required) {
  append_analysis_log(
    project_root = project_root,
    analysis_name = "04_balance_and_matching",
    output_files = decision_path,
    note = "Age and sex were not significantly imbalanced across diagnosis groups; pairwise matching was skipped for this result version.",
    summary_dir = result_summary_dir
  )
} else {
  group_var <- project_config$variables$group_label_var
  matching_vars <- unique(c(significant_continuous, significant_categorical))
  exact_vars <- significant_categorical

  pair_summary_list <- list()
  pair_balance_list <- list()
  pair_target_list <- list()
  matched_files <- c(decision_path)

  for (pair_levels in project_config$matching$pairs) {
    match_result <- greedy_match_pair(
      data = analysis_data,
      group_var = group_var,
      pair_levels = pair_levels,
      match_vars = matching_vars,
      exact_vars = exact_vars,
      transformation_table = transformation_plan,
      seed = project_config$seed
    )

    pair_name <- paste(pair_levels, collapse = "_vs_")
    matched_path <- file.path(result_data_clean_dir, paste0("matched_", pair_name, ".csv"))
    write_csv_utf8(match_result$matched_data, matched_path, row.names = FALSE)
    matched_files <- c(matched_files, matched_path)

    pair_summary_list[[pair_name]] <- match_result$summary

    for (var_name in significant_continuous) {
      bal_out <- run_continuous_group_comparison(
        data = match_result$matched_data,
        variable = var_name,
        group_var = group_var,
        transformation_table = transformation_plan,
        alpha = project_config$normality$alpha
      )$overall
      bal_out$pair_name <- pair_name
      pair_balance_list[[paste(pair_name, var_name, sep = "_")]] <- bal_out
    }

    for (var_name in significant_categorical) {
      bal_out <- run_categorical_group_comparison(
        data = match_result$matched_data,
        variable = var_name,
        group_var = group_var
      )$overall
      bal_out$pair_name <- pair_name
      pair_balance_list[[paste(pair_name, var_name, sep = "_")]] <- bal_out
    }

    for (target_var in project_config$group_comparison$target_vars) {
      target_out <- run_continuous_group_comparison(
        data = match_result$matched_data,
        variable = target_var,
        group_var = group_var,
        transformation_table = transformation_plan,
        alpha = project_config$normality$alpha
      )$overall
      target_out$pair_name <- pair_name
      pair_target_list[[paste(pair_name, target_var, sep = "_")]] <- target_out
    }
  }

  matching_summary_path <- file.path(result_summary_dir, "matching_pair_summary.csv")
  matched_balance_path <- file.path(result_summary_dir, "matching_balance_after_matching.csv")
  matched_target_path <- file.path(result_summary_dir, "matched_group_comparisons_targets.csv")

  write_csv_utf8(do.call(rbind, pair_summary_list), matching_summary_path, row.names = FALSE)
  write_csv_utf8(do.call(rbind, pair_balance_list), matched_balance_path, row.names = FALSE)
  write_csv_utf8(do.call(rbind, pair_target_list), matched_target_path, row.names = FALSE)

  append_analysis_log(
    project_root = project_root,
    analysis_name = "04_balance_and_matching",
    output_files = c(matched_files, matching_summary_path, matched_balance_path, matched_target_path),
    note = "Performed pairwise matching for diagnosis groups when age and/or sex imbalance was detected and stored matched datasets inside the current result version.",
    summary_dir = result_summary_dir
  )
}
