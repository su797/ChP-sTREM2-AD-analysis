fit_linear_model <- function(data, outcome, exposure, covariates = NULL, conf_level = 0.95) {
  rhs <- c(exposure, covariates)
  formula_text <- sprintf("%s ~ %s", outcome, paste(rhs, collapse = " + "))
  fit <- stats::lm(stats::as.formula(formula_text), data = data)
  tidy_lm_base(
    fit = fit,
    model_type = "linear",
    outcome = outcome,
    exposure = exposure,
    formula_text = formula_text
  )
}

fit_linear_models_batch <- function(data, outcomes, exposure, covariates = NULL, conf_level = 0.95) {
  out <- do.call(
    rbind,
    lapply(outcomes, function(y) {
      fit_linear_model(
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

extract_primary_terms <- function(model_results, remove_intercept = TRUE) {
  out <- model_results
  if (remove_intercept) {
    out <- out[out$term != "(Intercept)", , drop = FALSE]
  }
  out
}
