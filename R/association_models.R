build_formula_text <- function(outcome, predictors, factor_vars = NULL) {
  predictor_terms <- vapply(
    predictors,
    function(var_name) {
      if (!is.null(factor_vars) && var_name %in% factor_vars) {
        sprintf("factor(%s)", var_name)
      } else {
        var_name
      }
    },
    FUN.VALUE = character(1)
  )

  sprintf("%s ~ %s", outcome, paste(predictor_terms, collapse = " + "))
}

tidy_lm_base <- function(fit, model_type = "linear", outcome = NULL, exposure = NULL, formula_text = NULL) {
  coef_table <- summary(fit)$coefficients
  ci_table <- suppressWarnings(stats::confint(fit))

  out <- data.frame(
    term = rownames(coef_table),
    estimate = coef_table[, 1],
    std.error = coef_table[, 2],
    statistic = coef_table[, 3],
    p.value = coef_table[, 4],
    conf.low = ci_table[, 1],
    conf.high = ci_table[, 2],
    stringsAsFactors = FALSE
  )

  out$model_type <- model_type
  out$outcome <- outcome
  out$exposure <- exposure
  out$formula <- formula_text
  out
}

fit_linear_model_base <- function(data, outcome, exposure, covariates = NULL, factor_vars = NULL, transformation_table = NULL) {
  outcome_var <- resolve_analysis_var(outcome, transformation_table)
  exposure_var <- resolve_analysis_var(exposure, transformation_table)
  covariate_vars <- resolve_analysis_vars(covariates, transformation_table)

  model_df <- data.frame(
    outcome = as.numeric(data[[outcome_var]]),
    exposure = as.numeric(data[[exposure_var]])
  )

  predictors <- "exposure"
  if (!is.null(covariates) && length(covariates) > 0) {
    for (i in seq_along(covariates)) {
      cov_name <- covariates[i]
      analysis_name <- covariate_vars[i]
      if (!is.null(factor_vars) && cov_name %in% factor_vars) {
        model_df[[cov_name]] <- factor(data[[cov_name]])
      } else {
        model_df[[cov_name]] <- as.numeric(data[[analysis_name]])
      }
      predictors <- c(predictors, cov_name)
    }
  }

  model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
  formula_text <- build_formula_text("outcome", predictors, factor_vars = factor_vars)
  fit <- stats::lm(stats::as.formula(formula_text), data = model_df)
  out <- tidy_lm_base(
    fit = fit,
    model_type = "linear",
    outcome = outcome,
    exposure = exposure,
    formula_text = formula_text
  )

  exposure_beta <- out$estimate[out$term == "exposure"]
  out$std_beta <- NA_real_
  if (length(exposure_beta) == 1 && stats::sd(model_df$outcome) > 0) {
    out$std_beta[out$term == "exposure"] <- exposure_beta * stats::sd(model_df$exposure) / stats::sd(model_df$outcome)
  }

  out$n <- nrow(model_df)
  out
}

fit_linear_models_by_group <- function(data, group_var, outcome, exposure, covariates = NULL, factor_vars = NULL, transformation_table = NULL) {
  group_levels <- unique(as.character(data[[group_var]]))
  out <- do.call(
    rbind,
    lapply(group_levels, function(group_name) {
      sub_data <- data[as.character(data[[group_var]]) == group_name, , drop = FALSE]
      fit_out <- fit_linear_model_base(
        data = sub_data,
        outcome = outcome,
        exposure = exposure,
        covariates = covariates,
        factor_vars = factor_vars,
        transformation_table = transformation_table
      )
      fit_out$group <- group_name
      fit_out
    })
  )
  rownames(out) <- NULL
  out
}

run_partial_correlation <- function(data, x_var, y_var, covariates = NULL, factor_vars = NULL, transformation_table = NULL) {
  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  y_analysis <- resolve_analysis_var(y_var, transformation_table)
  cov_analysis <- resolve_analysis_vars(covariates, transformation_table)

  model_df <- data.frame(
    x = as.numeric(data[[x_analysis]]),
    y = as.numeric(data[[y_analysis]])
  )

  if (!is.null(covariates) && length(covariates) > 0) {
    for (i in seq_along(covariates)) {
      cov_name <- covariates[i]
      analysis_name <- cov_analysis[i]
      if (!is.null(factor_vars) && cov_name %in% factor_vars) {
        model_df[[cov_name]] <- factor(data[[cov_name]])
      } else {
        model_df[[cov_name]] <- as.numeric(data[[analysis_name]])
      }
    }
  }

  model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
  rhs <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    "1"
  }

  fit_x <- stats::lm(stats::as.formula(paste("x ~", rhs)), data = model_df)
  fit_y <- stats::lm(stats::as.formula(paste("y ~", rhs)), data = model_df)
  res_x <- stats::residuals(fit_x)
  res_y <- stats::residuals(fit_y)

  estimate <- stats::cor(res_x, res_y)
  k <- ifelse(is.null(covariates), 0, length(covariates))
  n <- length(res_x)
  df <- n - k - 2
  statistic <- estimate * sqrt(df / (1 - estimate^2))
  p_value <- 2 * stats::pt(-abs(statistic), df = df)

  fisher_z <- 0.5 * log((1 + estimate) / (1 - estimate))
  fisher_se <- 1 / sqrt(n - k - 3)
  z_low <- fisher_z - 1.96 * fisher_se
  z_high <- fisher_z + 1.96 * fisher_se

  data.frame(
    x = x_var,
    x_analysis = x_analysis,
    y = y_var,
    y_analysis = y_analysis,
    n = n,
    estimate = estimate,
    statistic = statistic,
    p_value = p_value,
    conf_low = (exp(2 * z_low) - 1) / (exp(2 * z_low) + 1),
    conf_high = (exp(2 * z_high) - 1) / (exp(2 * z_high) + 1),
    method = "partial_pearson",
    stringsAsFactors = FALSE
  )
}

run_partial_correlation_by_group <- function(data, group_var, x_var, y_var, covariates = NULL, factor_vars = NULL, transformation_table = NULL) {
  group_levels <- unique(as.character(data[[group_var]]))
  out <- do.call(
    rbind,
    lapply(group_levels, function(group_name) {
      sub_data <- data[as.character(data[[group_var]]) == group_name, , drop = FALSE]
      cor_out <- run_partial_correlation(
        data = sub_data,
        x_var = x_var,
        y_var = y_var,
        covariates = covariates,
        factor_vars = factor_vars,
        transformation_table = transformation_table
      )
      cor_out$group <- group_name
      cor_out
    })
  )
  rownames(out) <- NULL
  out
}
