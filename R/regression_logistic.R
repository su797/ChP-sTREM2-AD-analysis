fit_logistic_model <- function(data, outcome, exposure, covariates = NULL, conf_level = 0.95) {
  rhs <- c(exposure, covariates)
  formula_text <- sprintf("%s ~ %s", outcome, paste(rhs, collapse = " + "))
  fit <- stats::glm(
    stats::as.formula(formula_text),
    data = data,
    family = stats::binomial()
  )

  coef_table <- summary(fit)$coefficients
  ci_table <- suppressWarnings(stats::confint.default(fit))
  result <- data.frame(
    term = rownames(coef_table),
    estimate = coef_table[, 1],
    std.error = coef_table[, 2],
    statistic = coef_table[, 3],
    p.value = coef_table[, 4],
    conf.low = ci_table[, 1],
    conf.high = ci_table[, 2],
    stringsAsFactors = FALSE
  )
  result$model_type <- "logistic"
  result$outcome <- outcome
  result$exposure <- exposure
  result$formula <- formula_text
  result$or <- exp(result$estimate)
  result$or_ci_low <- exp(result$conf.low)
  result$or_ci_high <- exp(result$conf.high)
  result
}

fit_logistic_models_batch <- function(data, outcomes, exposure, covariates = NULL, conf_level = 0.95) {
  out <- do.call(
    rbind,
    lapply(outcomes, function(y) {
      fit_logistic_model(
        data = data,
        outcome = y,
        exposure = exposure,
        covariates = covariates,
        conf_level = conf_level
      )
    })
  )
  rownames(out) <- NULL
  out
}
