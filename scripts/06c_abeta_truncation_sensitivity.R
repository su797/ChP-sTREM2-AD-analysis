source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
group_var <- project_config$variables$group_label_var
group_order <- unlist(runtime_settings$group_order %||% c("CN", "MCI", "AD"), use.names = FALSE)
abeta_var <- project_config$advanced_models$amyloid_var %||% "S_ABETA"
abeta_threshold <- project_config$advanced_models$amyloid_ceiling_threshold %||% 1700
ceiling_flag_var <- "ABETA_ceiling_flag"

analysis_data[[ceiling_flag_var]] <- !is.na(analysis_data[[abeta_var]]) & analysis_data[[abeta_var]] >= abeta_threshold
restricted_data <- analysis_data[!analysis_data[[ceiling_flag_var]], , drop = FALSE]

bind_rows_fill <- function(...) {
  items <- list(...)
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
      values <- data[[source_var]]
      if (source_var == group_var) {
        out[[alias_name]] <- factor(as.character(values), levels = group_order)
      } else {
        out[[alias_name]] <- factor(values)
      }
    } else {
      analysis_var <- resolve_analysis_var(source_var, transformation_table)
      out[[alias_name]] <- as.numeric(data[[analysis_var]])
    }
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out[stats::complete.cases(out), , drop = FALSE]
}

tidy_custom_lm <- function(fit, model_df, model_name, outcome_name, exposure_name, formula_text, analysis_family, sample_set, group_name = NA_character_, term_label_map = NULL) {
  out <- tidy_lm_base(
    fit = fit,
    model_type = "linear",
    outcome = outcome_name,
    exposure = exposure_name,
    formula_text = formula_text
  )

  if (!is.null(term_label_map)) {
    out$term_label <- vapply(out$term, function(term_name) {
      if (term_name %in% names(term_label_map) && !is.null(term_label_map[[term_name]])) {
        return(term_label_map[[term_name]])
      }
      term_name
    }, FUN.VALUE = character(1))
  } else {
    out$term_label <- out$term
  }

  fit_summary <- summary(fit)
  out$model_name <- model_name
  out$analysis_family <- analysis_family
  out$sample_set <- sample_set
  out$group <- group_name
  out$adj_r_squared <- fit_summary$adj.r.squared
  out$r_squared <- fit_summary$r.squared
  out$aic <- stats::AIC(fit)
  out$bic <- stats::BIC(fit)
  out$n <- nrow(model_df)
  out
}

fit_biomarker_adjusted_model <- function(data, cfg, sample_set, group_name = NA_character_) {
  alias_map <- c(
    outcome = cfg$outcome,
    exposure = cfg$exposure
  )
  for (bio_var in cfg$biomarkers) {
    alias_map[[bio_var]] <- bio_var
  }
  for (cov_name in project_config$variables$covariates) {
    alias_map[[cov_name]] <- cov_name
  }

  factor_aliases <- project_config$variables$categorical_covariates
  model_df <- build_alias_df(
    data = data,
    alias_map = alias_map,
    factor_aliases = factor_aliases,
    transformation_table = transformation_plan
  )

  predictors <- c("exposure", cfg$biomarkers, project_config$variables$covariates)
  formula_text <- paste("outcome ~", paste(predictors, collapse = " + "))
  fit <- stats::lm(stats::as.formula(formula_text), data = model_df)

  term_label_map <- c(exposure = cfg$exposure)
  for (bio_var in cfg$biomarkers) {
    term_label_map[[bio_var]] <- bio_var
  }

  tidy_custom_lm(
    fit = fit,
    model_df = model_df,
    model_name = cfg$name,
    outcome_name = cfg$outcome,
    exposure_name = cfg$exposure,
    formula_text = formula_text,
    analysis_family = "abeta_truncation_biomarker_adjusted",
    sample_set = sample_set,
    group_name = group_name,
    term_label_map = term_label_map
  )
}

summarise_ceiling <- function(data, sample_set) {
  overall <- data.frame(
    sample_set = sample_set,
    group = "Overall",
    n_total = nrow(data),
    n_ceiling = sum(data[[ceiling_flag_var]], na.rm = TRUE),
    proportion_ceiling = mean(data[[ceiling_flag_var]], na.rm = TRUE),
    threshold = abeta_threshold,
    stringsAsFactors = FALSE
  )

  by_group <- do.call(
    rbind,
    lapply(group_order, function(group_name) {
      sub_data <- data[as.character(data[[group_var]]) == group_name, , drop = FALSE]
      data.frame(
        sample_set = sample_set,
        group = group_name,
        n_total = nrow(sub_data),
        n_ceiling = sum(sub_data[[ceiling_flag_var]], na.rm = TRUE),
        proportion_ceiling = mean(sub_data[[ceiling_flag_var]], na.rm = TRUE),
        threshold = abeta_threshold,
        stringsAsFactors = FALSE
      )
    })
  )

  rbind(overall, by_group)
}

simple_model_pairs <- list(
  list(name = "ChPICV_on_ABETA", outcome = "ChPICV", exposure = abeta_var),
  list(name = "ABETA_on_ChPICV", outcome = abeta_var, exposure = "ChPICV")
)

fit_simple_models <- function(data, sample_set) {
  overall_rows <- do.call(
    rbind,
    lapply(simple_model_pairs, function(pair_cfg) {
      out <- fit_linear_model_base(
        data = data,
        outcome = pair_cfg$outcome,
        exposure = pair_cfg$exposure,
        covariates = project_config$variables$covariates,
        factor_vars = project_config$variables$categorical_covariates,
        transformation_table = transformation_plan
      )
      out$model_name <- pair_cfg$name
      out$sample_set <- sample_set
      out$group <- "Overall"
      out$analysis_family <- "abeta_truncation_primary"
      out
    })
  )

  group_rows <- do.call(
    rbind,
    lapply(simple_model_pairs, function(pair_cfg) {
      out <- fit_linear_models_by_group(
        data = data,
        group_var = group_var,
        outcome = pair_cfg$outcome,
        exposure = pair_cfg$exposure,
        covariates = project_config$variables$covariates,
        factor_vars = project_config$variables$categorical_covariates,
        transformation_table = transformation_plan
      )
      out$model_name <- pair_cfg$name
      out$sample_set <- sample_set
      out$analysis_family <- "abeta_truncation_primary"
      out
    })
  )

  rbind(overall_rows, group_rows)
}

fit_partial_models <- function(data, sample_set) {
  overall <- run_partial_correlation(
    data = data,
    x_var = "ChPICV",
    y_var = abeta_var,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  overall$sample_set <- sample_set
  overall$group <- "Overall"

  by_group <- run_partial_correlation_by_group(
    data = data,
    group_var = group_var,
    x_var = "ChPICV",
    y_var = abeta_var,
    covariates = project_config$variables$covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  by_group$sample_set <- sample_set

  rbind(overall, by_group)
}

fit_adjusted_models <- function(data, sample_set) {
  rows <- do.call(
    rbind,
    lapply(project_config$advanced_models$biomarker_adjusted_models, function(cfg) {
      fit_biomarker_adjusted_model(data = data, cfg = cfg, sample_set = sample_set)
    })
  )
  rows
}

ceiling_summary <- rbind(
  summarise_ceiling(analysis_data, "full"),
  summarise_ceiling(restricted_data, "restricted_no_ge_1700")
)

simple_regression_rows <- rbind(
  fit_simple_models(analysis_data, "full"),
  fit_simple_models(restricted_data, "restricted_no_ge_1700")
)

partial_rows <- rbind(
  fit_partial_models(analysis_data, "full"),
  fit_partial_models(restricted_data, "restricted_no_ge_1700")
)

adjusted_rows <- rbind(
  fit_adjusted_models(analysis_data, "full"),
  fit_adjusted_models(restricted_data, "restricted_no_ge_1700")
)

regression_summary <- bind_rows_fill(
  simple_regression_rows,
  adjusted_rows
)

summary_path <- file.path(result_summary_dir, "abeta_truncation_summary.csv")
regression_path <- file.path(result_summary_dir, "abeta_truncation_regression_sensitivity.csv")
partial_path <- file.path(result_summary_dir, "abeta_truncation_partial_sensitivity.csv")
table_csv <- file.path(result_tables_dir, "table_abeta_truncation_sensitivity.csv")
table_html <- file.path(result_tables_dir, "table_abeta_truncation_sensitivity.html")

write_csv_utf8(ceiling_summary, summary_path, row.names = FALSE)
write_csv_utf8(regression_summary, regression_path, row.names = FALSE)
write_csv_utf8(partial_rows, partial_path, row.names = FALSE)

table_data <- regression_summary[
  regression_summary$term %in% c("exposure", "S_ABETA", "PTAU", "TAU") &
    regression_summary$group == "Overall",
  ,
  drop = FALSE
]
export_three_line_table(
  data = table_data,
  csv_path = table_csv,
  html_path = table_html,
  title = "amyloid-beta Ceiling Sensitivity Analysis"
)

append_analysis_log(
  project_root = project_root,
  analysis_name = "06c_abeta_truncation_sensitivity",
  output_files = c(summary_path, regression_path, partial_path, table_csv, table_html),
  note = "Completed amyloid-beta ceiling sensitivity analyses by comparing full-sample results with analyses excluding observations at the 1700 pg/mL assay ceiling.",
  summary_dir = result_summary_dir
)
