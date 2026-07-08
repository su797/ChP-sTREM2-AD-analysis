source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
raw_data_all <- read_project_data(raw_data_all_path)

cfg <- project_config$robustness$phase_site_protocol_sensitivity

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

clean_chr <- function(x) {
  trimws(as.character(x))
}

classify_t1_family <- function(x) {
  x <- tolower(clean_chr(x))
  out <- ifelse(
    grepl("mprage|mp-rage|mp_rage", x),
    "MPRAGE family",
    ifelse(
      grepl("spgr|fspgr", x),
      "IR-SPGR/FSPGR family",
      "Other T1 family"
    )
  )
  out[is.na(x) | x == ""] <- NA_character_
  out
}

metadata_keep <- unique(c(cfg$metadata_merge_keys, cfg$metadata_vars))
metadata_df <- raw_data_all[, intersect(metadata_keep, names(raw_data_all)), drop = FALSE]
metadata_df <- metadata_df[!duplicated(metadata_df[, cfg$metadata_merge_keys, drop = FALSE]), , drop = FALSE]

analysis_with_meta <- merge(
  analysis_data,
  metadata_df,
  by = cfg$metadata_merge_keys,
  all.x = TRUE,
  sort = FALSE
)

analysis_with_meta$ADNI_phase <- factor(clean_chr(analysis_with_meta[[cfg$phase_var]]), levels = cfg$phase_levels)
analysis_with_meta$T1_family <- factor(classify_t1_family(analysis_with_meta[[cfg$t1_var]]))

site_values <- clean_chr(analysis_with_meta[[cfg$site_var]])
site_counts <- table(site_values)
analysis_with_meta$SITE_group <- ifelse(
  site_values == "" | is.na(site_values),
  NA_character_,
  ifelse(
    site_counts[site_values] < cfg$sparse_site_min_n,
    "Other sites (<5 each)",
    paste0("Site ", site_values)
  )
)
analysis_with_meta$SITE_group <- factor(analysis_with_meta$SITE_group)

run_sensitivity_ladder <- function(data, model_cfg) {
  scopes <- list(
    list(
      scope = "primary_adjusted",
      covariates = project_config$variables$covariates,
      factor_vars = project_config$variables$categorical_covariates
    ),
    list(
      scope = "diagnosis_adjusted",
      covariates = c(project_config$variables$covariates, project_config$variables$group_label_var),
      factor_vars = c(project_config$variables$categorical_covariates, project_config$variables$group_label_var)
    ),
    list(
      scope = "diagnosis_phase_adjusted",
      covariates = c(project_config$variables$covariates, project_config$variables$group_label_var, "ADNI_phase"),
      factor_vars = c(project_config$variables$categorical_covariates, project_config$variables$group_label_var, "ADNI_phase")
    ),
    list(
      scope = "diagnosis_protocol_adjusted",
      covariates = c(project_config$variables$covariates, project_config$variables$group_label_var, "T1_family"),
      factor_vars = c(project_config$variables$categorical_covariates, project_config$variables$group_label_var, "T1_family")
    ),
    list(
      scope = "diagnosis_site_adjusted",
      covariates = c(project_config$variables$covariates, project_config$variables$group_label_var, "SITE_group"),
      factor_vars = c(project_config$variables$categorical_covariates, project_config$variables$group_label_var, "SITE_group")
    ),
    list(
      scope = "diagnosis_phase_protocol_site_adjusted",
      covariates = c(project_config$variables$covariates, project_config$variables$group_label_var, "ADNI_phase", "T1_family", "SITE_group"),
      factor_vars = c(project_config$variables$categorical_covariates, project_config$variables$group_label_var, "ADNI_phase", "T1_family", "SITE_group")
    )
  )

  rows <- lapply(scopes, function(scope_cfg) {
    fit_rows <- fit_linear_model_base(
      data = data,
      outcome = model_cfg$outcome,
      exposure = model_cfg$exposure,
      covariates = scope_cfg$covariates,
      factor_vars = scope_cfg$factor_vars,
      transformation_table = transformation_plan
    )
    fit_rows <- fit_rows[fit_rows$term == "exposure", , drop = FALSE]
    fit_rows$model_name <- model_cfg$name
    fit_rows$sensitivity_scope <- scope_cfg$scope
    fit_rows
  })

  bind_rows_fill(rows)
}

sensitivity_rows <- bind_rows_fill(lapply(cfg$core_models, function(model_cfg) {
  run_sensitivity_ladder(analysis_with_meta, model_cfg)
}))
if (nrow(sensitivity_rows) > 0 && "p.value" %in% names(sensitivity_rows)) {
  sensitivity_rows$q_value <- stats::p.adjust(sensitivity_rows$p.value, method = "fdr")
}

counts_rows <- rbind(
  data.frame(
    dimension = "ADNI phase",
    level = names(sort(table(analysis_with_meta$ADNI_phase), decreasing = TRUE)),
    n = as.integer(sort(table(analysis_with_meta$ADNI_phase), decreasing = TRUE)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dimension = "T1 family",
    level = names(sort(table(analysis_with_meta$T1_family), decreasing = TRUE)),
    n = as.integer(sort(table(analysis_with_meta$T1_family), decreasing = TRUE)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    dimension = "Grouped study site",
    level = names(sort(table(analysis_with_meta$SITE_group), decreasing = TRUE)),
    n = as.integer(sort(table(analysis_with_meta$SITE_group), decreasing = TRUE)),
    stringsAsFactors = FALSE
  )
)
counts_rows <- counts_rows[!is.na(counts_rows$level) & counts_rows$level != "", , drop = FALSE]

source_summary <- data.frame(
  n_total = nrow(analysis_with_meta),
  unique_sites = length(unique(na.omit(clean_chr(analysis_with_meta[[cfg$site_var]])))),
  grouped_sites = length(unique(na.omit(as.character(analysis_with_meta$SITE_group)))),
  sparse_site_min_n = cfg$sparse_site_min_n,
  phase_distribution = paste(
    names(sort(table(analysis_with_meta$ADNI_phase), decreasing = TRUE)),
    as.integer(sort(table(analysis_with_meta$ADNI_phase), decreasing = TRUE)),
    collapse = "; "
  ),
  t1_family_distribution = paste(
    names(sort(table(analysis_with_meta$T1_family), decreasing = TRUE)),
    as.integer(sort(table(analysis_with_meta$T1_family), decreasing = TRUE)),
    collapse = "; "
  ),
  stringsAsFactors = FALSE
)

sensitivity_path <- file.path(result_summary_dir, "phase_site_protocol_sensitivity.csv")
counts_path <- file.path(result_summary_dir, "phase_site_protocol_counts.csv")
source_summary_path <- file.path(result_summary_dir, "phase_site_protocol_summary.csv")

write_csv_utf8(sensitivity_rows, sensitivity_path, row.names = FALSE)
write_csv_utf8(counts_rows, counts_path, row.names = FALSE)
write_csv_utf8(source_summary, source_summary_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "06g_phase_site_protocol_sensitivity",
  output_files = c(sensitivity_path, counts_path, source_summary_path),
  note = "Added phase/site/T1-protocol sensitivity analyses for the four core ChP-biomarker association models.",
  summary_dir = result_summary_dir
)
