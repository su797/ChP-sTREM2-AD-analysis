source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
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

build_alias_df <- function(data, alias_map, factor_aliases = character(0), transformation_table = NULL) {
  out <- list()
  for (alias_name in names(alias_map)) {
    source_var <- alias_map[[alias_name]]
    if (alias_name %in% factor_aliases) {
      out[[alias_name]] <- factor(as.character(data[[source_var]]))
    } else {
      analysis_var <- resolve_analysis_var(source_var, transformation_table)
      out[[alias_name]] <- as.numeric(data[[analysis_var]])
    }
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out[stats::complete.cases(out), , drop = FALSE]
}

metadata_cfg <- project_config$robustness$biologic_sensitivity
metadata_keep <- unique(c(metadata_cfg$metadata_merge_keys, metadata_cfg$metadata_vars))
metadata_df <- raw_data_all[, intersect(metadata_keep, names(raw_data_all)), drop = FALSE]
metadata_df <- metadata_df[!duplicated(metadata_df[, metadata_cfg$metadata_merge_keys, drop = FALSE]), , drop = FALSE]

analysis_with_meta <- merge(
  analysis_data,
  metadata_df,
  by = metadata_cfg$metadata_merge_keys,
  all.x = TRUE,
  sort = FALSE
)

fit_diagnosis_hierarchy <- function(data, cfg) {
  base_rows <- fit_linear_model_base(
    data = data,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  base_rows <- base_rows[base_rows$term == "exposure", , drop = FALSE]
  base_rows$model_scope <- "primary_adjusted"
  base_rows$group <- "Overall"
  base_rows$model_name <- cfg$name

  dx_rows <- fit_linear_model_base(
    data = data,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    covariates = c(project_config$variables$covariates, project_config$variables$group_label_var),
    factor_vars = c(project_config$variables$categorical_covariates, project_config$variables$group_label_var),
    transformation_table = transformation_plan
  )
  dx_rows <- dx_rows[dx_rows$term == "exposure", , drop = FALSE]
  dx_rows$model_scope <- "diagnosis_adjusted"
  dx_rows$group <- "Overall"
  dx_rows$model_name <- cfg$name

  by_group_rows <- fit_linear_models_by_group(
    data = data,
    group_var = project_config$variables$group_label_var,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  by_group_rows <- by_group_rows[by_group_rows$term == "exposure", , drop = FALSE]
  by_group_rows$model_scope <- "within_group"
  by_group_rows$model_name <- cfg$name

  alias_map <- c(
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    diagnosis = project_config$variables$group_label_var
  )
  for (cov_name in project_config$variables$covariates) {
    alias_map[[cov_name]] <- cov_name
  }
  comparison_df <- build_alias_df(
    data = data,
    alias_map = alias_map,
    factor_aliases = c(project_config$variables$categorical_covariates, "diagnosis"),
    transformation_table = transformation_plan
  )
  comparison_df$diagnosis <- stats::relevel(comparison_df$diagnosis, ref = "CN")
  reduced_formula <- build_formula_text(
    "outcome",
    c("exposure", project_config$variables$covariates),
    factor_vars = project_config$variables$categorical_covariates
  )
  full_formula <- build_formula_text(
    "outcome",
    c("exposure", project_config$variables$covariates, "diagnosis"),
    factor_vars = c(project_config$variables$categorical_covariates, "diagnosis")
  )
  reduced_fit <- stats::lm(stats::as.formula(reduced_formula), data = comparison_df)
  full_fit <- stats::lm(stats::as.formula(full_formula), data = comparison_df)
  compare_tab <- stats::anova(reduced_fit, full_fit)
  comparison_out <- data.frame(
    model_name = cfg$name,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    reduced_adj_r_squared = summary(reduced_fit)$adj.r.squared,
    full_adj_r_squared = summary(full_fit)$adj.r.squared,
    delta_adj_r_squared = summary(full_fit)$adj.r.squared - summary(reduced_fit)$adj.r.squared,
    compare_df = compare_tab$Df[[2]],
    compare_f = compare_tab$F[[2]],
    compare_p = compare_tab$`Pr(>F)`[[2]],
    n = nrow(comparison_df),
    stringsAsFactors = FALSE
  )

  list(
    hierarchy_rows = rbind(base_rows, dx_rows, by_group_rows),
    comparison_out = comparison_out
  )
}

run_biologic_sensitivity <- function(data, subset_cfg, model_cfg) {
  status_var <- metadata_cfg$stable_amyloid_var
  sub_data <- data[
    !is.na(data[[status_var]]) &
      as.character(data[[status_var]]) == subset_cfg$amyloid_status &
      as.character(data[[project_config$variables$group_label_var]]) %in% subset_cfg$groups,
    ,
    drop = FALSE
  ]
  if (nrow(sub_data) < 25) {
    return(data.frame())
  }

  fit_rows <- fit_linear_model_base(
    data = sub_data,
    outcome = model_cfg$outcome,
    exposure = model_cfg$exposure,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  fit_rows <- fit_rows[fit_rows$term == "exposure", , drop = FALSE]
  fit_rows$subset_name <- subset_cfg$name
  fit_rows$subset_label <- subset_cfg$label
  fit_rows$amyloid_status <- subset_cfg$amyloid_status
  fit_rows$groups_included <- paste(subset_cfg$groups, collapse = ", ")
  fit_rows$model_name <- model_cfg$name
  fit_rows$n_subset <- nrow(sub_data)
  fit_rows
}

hierarchy_list <- lapply(
  project_config$robustness$diagnosis_hierarchical_models,
  function(cfg) fit_diagnosis_hierarchy(analysis_with_meta, cfg)
)
hierarchy_rows <- bind_rows_fill(lapply(hierarchy_list, `[[`, "hierarchy_rows"))
hierarchy_compare <- bind_rows_fill(lapply(hierarchy_list, `[[`, "comparison_out"))
if (nrow(hierarchy_rows) > 0 && "p.value" %in% names(hierarchy_rows)) {
  hierarchy_rows$q_value <- stats::p.adjust(hierarchy_rows$p.value, method = "fdr")
}
if (nrow(hierarchy_compare) > 0 && "compare_p" %in% names(hierarchy_compare)) {
  hierarchy_compare$q_value <- stats::p.adjust(hierarchy_compare$compare_p, method = "fdr")
}

biologic_rows <- bind_rows_fill(unlist(
  lapply(metadata_cfg$stable_subsets, function(subset_cfg) {
    lapply(metadata_cfg$core_models, function(model_cfg) {
      run_biologic_sensitivity(analysis_with_meta, subset_cfg, model_cfg)
    })
  }),
  recursive = FALSE
))
if (nrow(biologic_rows) > 0 && "p.value" %in% names(biologic_rows)) {
  biologic_rows$q_value <- stats::p.adjust(biologic_rows$p.value, method = "fdr")
}

biologic_counts <- bind_rows_fill(lapply(metadata_cfg$stable_subsets, function(subset_cfg) {
  status_var <- metadata_cfg$stable_amyloid_var
  sub_data <- analysis_with_meta[
    !is.na(analysis_with_meta[[status_var]]) &
      as.character(analysis_with_meta[[status_var]]) == subset_cfg$amyloid_status &
      as.character(analysis_with_meta[[project_config$variables$group_label_var]]) %in% subset_cfg$groups,
    ,
    drop = FALSE
  ]
  data.frame(
    subset_name = subset_cfg$name,
    subset_label = subset_cfg$label,
    amyloid_status = subset_cfg$amyloid_status,
    groups_included = paste(subset_cfg$groups, collapse = ", "),
    n = nrow(sub_data),
    stringsAsFactors = FALSE
  )
}))

source_summary <- data.frame(
  total_n = nrow(analysis_with_meta),
  adni_phase_distribution = paste(
    names(sort(table(analysis_with_meta$COLPROT), decreasing = TRUE)),
    as.integer(sort(table(analysis_with_meta$COLPROT), decreasing = TRUE)),
    collapse = "; "
  ),
  dx_bl_distribution = paste(
    names(sort(table(analysis_with_meta$DX_bl), decreasing = TRUE)),
    as.integer(sort(table(analysis_with_meta$DX_bl), decreasing = TRUE)),
    collapse = "; "
  ),
  t1_sequence_distribution = paste(
    names(sort(table(analysis_with_meta$T1WI), decreasing = TRUE)),
    as.integer(sort(table(analysis_with_meta$T1WI), decreasing = TRUE)),
    collapse = "; "
  ),
  stringsAsFactors = FALSE
)

hierarchy_path <- file.path(result_summary_dir, "diagnosis_hierarchical_models.csv")
hierarchy_compare_path <- file.path(result_summary_dir, "diagnosis_hierarchical_model_comparison.csv")
biologic_path <- file.path(result_summary_dir, "biologic_sensitivity_models.csv")
biologic_counts_path <- file.path(result_summary_dir, "biologic_sensitivity_counts.csv")
source_summary_path <- file.path(result_summary_dir, "sample_source_summary.csv")

write_csv_utf8(hierarchy_rows, hierarchy_path, row.names = FALSE)
write_csv_utf8(hierarchy_compare, hierarchy_compare_path, row.names = FALSE)
write_csv_utf8(biologic_rows, biologic_path, row.names = FALSE)
write_csv_utf8(biologic_counts, biologic_counts_path, row.names = FALSE)
write_csv_utf8(source_summary, source_summary_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "06f_diagnosis_biologic_sensitivity",
  output_files = c(
    hierarchy_path,
    hierarchy_compare_path,
    biologic_path,
    biologic_counts_path,
    source_summary_path
  ),
  note = "Added diagnosis-adjusted hierarchical regression models, biologically defined amyloid-stable sensitivity analyses, and ADNI source/sequence summary tables.",
  summary_dir = result_summary_dir
)
