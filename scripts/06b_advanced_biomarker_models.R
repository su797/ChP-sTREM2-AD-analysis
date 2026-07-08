source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
group_var <- project_config$variables$group_label_var
group_order <- unlist(runtime_settings$group_order %||% c("CN", "MCI", "AD"), use.names = FALSE)

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

tidy_custom_lm <- function(fit, model_df, model_name, outcome_name, exposure_name, formula_text, analysis_family, group_name = NA_character_, term_label_map = NULL) {
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
  out$group <- group_name
  out$adj_r_squared <- fit_summary$adj.r.squared
  out$r_squared <- fit_summary$r.squared
  out$aic <- stats::AIC(fit)
  out$bic <- stats::BIC(fit)
  out$n <- nrow(model_df)
  out
}

fit_biomarker_adjusted_model <- function(data, cfg, group_name = NA_character_) {
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
    analysis_family = "biomarker_adjusted",
    group_name = group_name,
    term_label_map = term_label_map
  )
}

fit_interaction_model <- function(data, cfg) {
  alias_map <- c(
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    diagnosis = cfg$moderator
  )
  for (cov_name in project_config$variables$covariates) {
    alias_map[[cov_name]] <- cov_name
  }

  factor_aliases <- unique(c(project_config$variables$categorical_covariates, "diagnosis"))
  model_df <- build_alias_df(
    data = data,
    alias_map = alias_map,
    factor_aliases = factor_aliases,
    transformation_table = transformation_plan
  )
  model_df$diagnosis <- stats::relevel(model_df$diagnosis, ref = group_order[[1]])

  reduced_formula <- paste("outcome ~ exposure + diagnosis +", paste(project_config$variables$covariates, collapse = " + "))
  full_formula <- paste("outcome ~ exposure * diagnosis +", paste(project_config$variables$covariates, collapse = " + "))

  reduced_fit <- stats::lm(stats::as.formula(reduced_formula), data = model_df)
  full_fit <- stats::lm(stats::as.formula(full_formula), data = model_df)

  term_label_map <- list(
    exposure = cfg$exposure,
    diagnosisMCI = paste0(cfg$moderator, ": MCI vs CN"),
    diagnosisAD = paste0(cfg$moderator, ": AD vs CN"),
    "exposure:diagnosisMCI" = paste0(cfg$exposure, " x MCI"),
    "exposure:diagnosisAD" = paste0(cfg$exposure, " x AD")
  )

  terms_out <- tidy_custom_lm(
    fit = full_fit,
    model_df = model_df,
    model_name = cfg$name,
    outcome_name = cfg$outcome,
    exposure_name = cfg$exposure,
    formula_text = full_formula,
    analysis_family = "interaction",
    group_name = "overall",
    term_label_map = term_label_map
  )

  compare_tab <- stats::anova(reduced_fit, full_fit)
  comparison_out <- data.frame(
    model_name = cfg$name,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    moderator = cfg$moderator,
    n = nrow(model_df),
    reduced_formula = reduced_formula,
    full_formula = full_formula,
    reduced_aic = stats::AIC(reduced_fit),
    full_aic = stats::AIC(full_fit),
    reduced_adj_r_squared = summary(reduced_fit)$adj.r.squared,
    full_adj_r_squared = summary(full_fit)$adj.r.squared,
    interaction_df = compare_tab$Df[[2]],
    interaction_f = compare_tab$F[[2]],
    interaction_p_value = compare_tab$`Pr(>F)`[[2]],
    stringsAsFactors = FALSE
  )

  list(terms = terms_out, comparison = comparison_out)
}

fit_nonlinear_model_family <- function(data, cfg) {
  alias_map <- c(
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    group = group_var
  )
  for (cov_name in project_config$variables$covariates) {
    alias_map[[cov_name]] <- cov_name
  }

  factor_aliases <- unique(c(project_config$variables$categorical_covariates, "group"))
  model_df <- build_alias_df(
    data = data,
    alias_map = alias_map,
    factor_aliases = factor_aliases,
    transformation_table = transformation_plan
  )

  covariate_rhs <- paste(project_config$variables$covariates, collapse = " + ")
  linear_formula <- paste("outcome ~ exposure +", covariate_rhs)
  quadratic_formula <- paste("outcome ~ exposure + I(exposure^2) +", covariate_rhs)
  spline_formula <- paste("outcome ~ splines::ns(exposure, df =", cfg$spline_df, ") +", covariate_rhs)

  linear_fit <- stats::lm(stats::as.formula(linear_formula), data = model_df)
  quadratic_fit <- stats::lm(stats::as.formula(quadratic_formula), data = model_df)
  spline_fit <- stats::lm(stats::as.formula(spline_formula), data = model_df)

  model_summary <- data.frame(
    model_name = cfg$name,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    model_type = c("linear", "quadratic", "spline"),
    formula = c(linear_formula, quadratic_formula, spline_formula),
    n = nrow(model_df),
    aic = c(stats::AIC(linear_fit), stats::AIC(quadratic_fit), stats::AIC(spline_fit)),
    bic = c(stats::BIC(linear_fit), stats::BIC(quadratic_fit), stats::BIC(spline_fit)),
    adj_r_squared = c(summary(linear_fit)$adj.r.squared, summary(quadratic_fit)$adj.r.squared, summary(spline_fit)$adj.r.squared),
    stringsAsFactors = FALSE
  )
  model_summary$best_by_aic <- model_summary$aic == min(model_summary$aic, na.rm = TRUE)

  quadratic_terms <- tidy_custom_lm(
    fit = quadratic_fit,
    model_df = model_df,
    model_name = cfg$name,
    outcome_name = cfg$outcome,
    exposure_name = cfg$exposure,
    formula_text = quadratic_formula,
    analysis_family = "nonlinear_quadratic",
    group_name = "overall",
    term_label_map = list(exposure = cfg$exposure, "I(exposure^2)" = paste0(cfg$exposure, "^2"))
  )

  quad_row <- quadratic_terms[quadratic_terms$term == "I(exposure^2)", , drop = FALSE]
  spline_cmp <- stats::anova(linear_fit, spline_fit)
  linear_quad_cmp <- stats::anova(linear_fit, quadratic_fit)
  nonlinear_tests <- data.frame(
    model_name = cfg$name,
    outcome = cfg$outcome,
    exposure = cfg$exposure,
    test_name = c("quadratic_term", "linear_vs_quadratic", "linear_vs_spline"),
    statistic = c(
      if (nrow(quad_row) > 0) quad_row$statistic[[1]] else NA_real_,
      linear_quad_cmp$F[[2]],
      spline_cmp$F[[2]]
    ),
    p_value = c(
      if (nrow(quad_row) > 0) quad_row$p.value[[1]] else NA_real_,
      linear_quad_cmp$`Pr(>F)`[[2]],
      spline_cmp$`Pr(>F)`[[2]]
    ),
    stringsAsFactors = FALSE
  )

  list(
    model_df = model_df,
    model_summary = model_summary,
    nonlinear_tests = nonlinear_tests,
    linear_fit = linear_fit,
    quadratic_fit = quadratic_fit,
    spline_fit = spline_fit
  )
}

most_common_level <- function(x) {
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[[1]]
}

save_nonlinearity_plot <- function(model_info, cfg, path) {
  model_df <- model_info$model_df
  exposure_grid <- seq(min(model_df$exposure, na.rm = TRUE), max(model_df$exposure, na.rm = TRUE), length.out = 200)

  newdata <- data.frame(exposure = exposure_grid)
  for (cov_name in project_config$variables$covariates) {
    if (is.factor(model_df[[cov_name]])) {
      newdata[[cov_name]] <- factor(
        rep(most_common_level(model_df[[cov_name]]), length(exposure_grid)),
        levels = levels(model_df[[cov_name]])
      )
    } else {
      newdata[[cov_name]] <- rep(mean(model_df[[cov_name]], na.rm = TRUE), length(exposure_grid))
    }
  }

  linear_pred <- stats::predict(model_info$linear_fit, newdata = newdata)
  quadratic_pred <- stats::predict(model_info$quadratic_fit, newdata = newdata)
  spline_pred <- stats::predict(model_info$spline_fit, newdata = newdata)

  png(path, width = 1500, height = 1100, res = 180)
  graphics::par(mar = c(5, 5, 4, 2))
  point_cols <- c(CN = "#2b8cbe", MCI = "#fdae61", AD = "#d7191c")
  plot_cols <- point_cols[as.character(model_df$group)]
  plot_cols[is.na(plot_cols)] <- "#666666"
  graphics::plot(
    model_df$exposure,
    model_df$outcome,
    pch = 16,
    col = grDevices::adjustcolor(plot_cols, alpha.f = 0.55),
    xlab = cfg$exposure,
    ylab = cfg$outcome,
    main = paste(cfg$outcome, "vs", cfg$exposure, "- linearity check")
  )
  graphics::lines(exposure_grid, linear_pred, lwd = 3, col = "#1b9e77")
  graphics::lines(exposure_grid, quadratic_pred, lwd = 3, col = "#d95f02")
  graphics::lines(exposure_grid, spline_pred, lwd = 3, col = "#7570b3")
  graphics::legend(
    "topright",
    legend = c("Observed", "Linear", "Quadratic", "Spline", "CN", "MCI", "AD"),
    col = c("#666666", "#1b9e77", "#d95f02", "#7570b3", point_cols[["CN"]], point_cols[["MCI"]], point_cols[["AD"]]),
    pch = c(16, NA, NA, NA, 16, 16, 16),
    lty = c(NA, 1, 1, 1, NA, NA, NA),
    lwd = c(NA, 3, 3, 3, NA, NA, NA),
    bty = "n"
  )
  grDevices::dev.off()
}

biomarker_overall_list <- list()
biomarker_by_group_list <- list()
for (cfg in project_config$advanced_models$biomarker_adjusted_models) {
  biomarker_overall_list[[cfg$name]] <- fit_biomarker_adjusted_model(analysis_data, cfg)
  biomarker_by_group_list[[cfg$name]] <- do.call(
    rbind,
    lapply(group_order, function(group_name) {
      sub_data <- analysis_data[as.character(analysis_data[[group_var]]) == group_name, , drop = FALSE]
      fit_biomarker_adjusted_model(sub_data, cfg, group_name = group_name)
    })
  )
}

interaction_results <- lapply(project_config$advanced_models$interaction_models, function(cfg) fit_interaction_model(analysis_data, cfg))
names(interaction_results) <- vapply(project_config$advanced_models$interaction_models, `[[`, FUN.VALUE = character(1), "name")

nonlinear_results <- lapply(project_config$advanced_models$nonlinear_models, function(cfg) fit_nonlinear_model_family(analysis_data, cfg))
names(nonlinear_results) <- vapply(project_config$advanced_models$nonlinear_models, `[[`, FUN.VALUE = character(1), "name")

nonlinear_plot_files <- character(0)
for (cfg in project_config$advanced_models$nonlinear_models) {
  plot_path <- file.path(result_figures_dir, paste0("nonlinear_", cfg$name, "_overall.png"))
  save_nonlinearity_plot(nonlinear_results[[cfg$name]], cfg, plot_path)
  nonlinear_plot_files <- c(nonlinear_plot_files, plot_path)
}

biomarker_overall <- do.call(rbind, biomarker_overall_list)
biomarker_by_group <- do.call(rbind, biomarker_by_group_list)
interaction_terms <- do.call(rbind, lapply(interaction_results, `[[`, "terms"))
interaction_comparison <- do.call(rbind, lapply(interaction_results, `[[`, "comparison"))
nonlinear_model_summary <- do.call(rbind, lapply(nonlinear_results, `[[`, "model_summary"))
nonlinear_tests <- do.call(rbind, lapply(nonlinear_results, `[[`, "nonlinear_tests"))

write_csv_utf8(biomarker_overall, file.path(result_summary_dir, "advanced_biomarker_adjusted_overall.csv"), row.names = FALSE)
write_csv_utf8(biomarker_by_group, file.path(result_summary_dir, "advanced_biomarker_adjusted_by_group.csv"), row.names = FALSE)
write_csv_utf8(interaction_terms, file.path(result_summary_dir, "advanced_interaction_terms.csv"), row.names = FALSE)
write_csv_utf8(interaction_comparison, file.path(result_summary_dir, "advanced_interaction_comparison.csv"), row.names = FALSE)
write_csv_utf8(nonlinear_model_summary, file.path(result_summary_dir, "advanced_nonlinear_model_summary.csv"), row.names = FALSE)
write_csv_utf8(nonlinear_tests, file.path(result_summary_dir, "advanced_nonlinear_tests.csv"), row.names = FALSE)

export_three_line_table(
  data = biomarker_overall[biomarker_overall$term %in% c("exposure", "PTAU", "TAU", "S_ABETA"), , drop = FALSE],
  csv_path = file.path(result_tables_dir, "table_advanced_biomarker_adjusted.csv"),
  html_path = file.path(result_tables_dir, "table_advanced_biomarker_adjusted.html"),
  title = "Biomarker-adjusted ChPICV Models"
)

export_three_line_table(
  data = interaction_comparison,
  csv_path = file.path(result_tables_dir, "table_advanced_interaction_comparison.csv"),
  html_path = file.path(result_tables_dir, "table_advanced_interaction_comparison.html"),
  title = "Interaction Model Comparison"
)

append_analysis_log(
  project_root = project_root,
  analysis_name = "06b_advanced_biomarker_models",
  output_files = c(
    file.path(result_summary_dir, "advanced_biomarker_adjusted_overall.csv"),
    file.path(result_summary_dir, "advanced_biomarker_adjusted_by_group.csv"),
    file.path(result_summary_dir, "advanced_interaction_terms.csv"),
    file.path(result_summary_dir, "advanced_interaction_comparison.csv"),
    file.path(result_summary_dir, "advanced_nonlinear_model_summary.csv"),
    file.path(result_summary_dir, "advanced_nonlinear_tests.csv"),
    file.path(result_tables_dir, "table_advanced_biomarker_adjusted.csv"),
    file.path(result_tables_dir, "table_advanced_biomarker_adjusted.html"),
    file.path(result_tables_dir, "table_advanced_interaction_comparison.csv"),
    file.path(result_tables_dir, "table_advanced_interaction_comparison.html"),
    nonlinear_plot_files
  ),
  note = "Completed biomarker-adjusted regression, diagnosis interaction models, and nonlinearity checks for the ChPICV-sTREM2 relationship.",
  summary_dir = result_summary_dir
)
