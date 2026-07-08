source(file.path(getwd(), "00_setup.R"))

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)

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

weighted_mean_base <- function(x, w) {
  sum(w * x) / sum(w)
}

weighted_var_base <- function(x, w) {
  mu <- weighted_mean_base(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

weighted_sd_base <- function(x, w) {
  sqrt(weighted_var_base(x, w))
}

compute_smd_numeric_weighted <- function(x, g, w = NULL) {
  keep <- stats::complete.cases(x, g, w)
  x <- as.numeric(x[keep])
  g <- as.character(g[keep])
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w[keep])
  }
  if (length(unique(g)) != 2) return(NA_real_)
  lev <- unique(g)
  x1 <- x[g == lev[1]]
  x2 <- x[g == lev[2]]
  w1 <- w[g == lev[1]]
  w2 <- w[g == lev[2]]
  if (length(x1) < 2 || length(x2) < 2 || sum(w1) <= 0 || sum(w2) <= 0) return(NA_real_)
  m1 <- weighted_mean_base(x1, w1)
  m2 <- weighted_mean_base(x2, w2)
  s1 <- weighted_sd_base(x1, w1)
  s2 <- weighted_sd_base(x2, w2)
  sp <- sqrt((s1^2 + s2^2) / 2)
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (m1 - m2) / sp
}

compute_smd_binary_weighted <- function(x, g, w = NULL) {
  keep <- stats::complete.cases(x, g, w)
  x <- as.numeric(x[keep])
  g <- as.character(g[keep])
  if (is.null(w)) {
    w <- rep(1, length(x))
  } else {
    w <- as.numeric(w[keep])
  }
  if (length(unique(g)) != 2) return(NA_real_)
  lev <- unique(g)
  x1 <- x[g == lev[1]]
  x2 <- x[g == lev[2]]
  w1 <- w[g == lev[1]]
  w2 <- w[g == lev[2]]
  if (sum(w1) <= 0 || sum(w2) <= 0) return(NA_real_)
  p1 <- weighted_mean_base(x1, w1)
  p2 <- weighted_mean_base(x2, w2)
  p <- (p1 + p2) / 2
  denom <- sqrt(p * (1 - p))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (p1 - p2) / denom
}

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

tidy_weighted_lm <- function(fit, outcome = NULL, exposure = NULL, formula_text = NULL, analysis_type = "overlap_weighted") {
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
    outcome = outcome,
    exposure = exposure,
    formula = formula_text,
    analysis_type = analysis_type,
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

compute_model_scale_factor <- function(data, outcome, exposure, covariates = NULL, factor_vars = NULL, transformation_table = NULL) {
  outcome_var <- resolve_analysis_var(outcome, transformation_table)
  exposure_var <- resolve_analysis_var(exposure, transformation_table)
  covariate_vars <- resolve_analysis_vars(covariates, transformation_table)

  model_df <- data.frame(
    outcome = as.numeric(data[[outcome_var]]),
    exposure = as.numeric(data[[exposure_var]])
  )

  if (!is.null(covariates) && length(covariates) > 0) {
    for (i in seq_along(covariates)) {
      cov_name <- covariates[i]
      analysis_name <- covariate_vars[i]
      if (!is.null(factor_vars) && cov_name %in% factor_vars) {
        model_df[[cov_name]] <- factor(data[[cov_name]])
      } else {
        model_df[[cov_name]] <- as.numeric(data[[analysis_name]])
      }
    }
  }

  model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
  if (nrow(model_df) < 5) return(NA_real_)
  stats::sd(model_df$exposure) / stats::sd(model_df$outcome)
}

fit_weighted_linear_model_base <- function(data, outcome, exposure, covariates = NULL, factor_vars = NULL, transformation_table = NULL, weights_var) {
  outcome_var <- resolve_analysis_var(outcome, transformation_table)
  exposure_var <- resolve_analysis_var(exposure, transformation_table)
  covariate_vars <- resolve_analysis_vars(covariates, transformation_table)

  model_df <- data.frame(
    outcome = as.numeric(data[[outcome_var]]),
    exposure = as.numeric(data[[exposure_var]]),
    weights = as.numeric(data[[weights_var]])
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
  model_df <- model_df[is.finite(model_df$weights) & model_df$weights > 0, , drop = FALSE]
  formula_text <- build_formula_text("outcome", predictors, factor_vars = factor_vars)
  fit <- stats::lm(stats::as.formula(formula_text), data = model_df, weights = weights)
  out <- tidy_weighted_lm(
    fit = fit,
    outcome = outcome,
    exposure = exposure,
    formula_text = formula_text
  )

  exposure_beta <- out$estimate[out$term == "exposure"]
  out$std_beta <- NA_real_
  out$std_conf.low <- NA_real_
  out$std_conf.high <- NA_real_
  if (length(exposure_beta) == 1 && weighted_sd_base(model_df$outcome, model_df$weights) > 0) {
    scale_factor <- weighted_sd_base(model_df$exposure, model_df$weights) / weighted_sd_base(model_df$outcome, model_df$weights)
    out$std_beta[out$term == "exposure"] <- exposure_beta * scale_factor
    out$std_conf.low[out$term == "exposure"] <- out$conf.low[out$term == "exposure"] * scale_factor
    out$std_conf.high[out$term == "exposure"] <- out$conf.high[out$term == "exposure"] * scale_factor
  }
  out$n <- nrow(model_df)
  out$sum_weights <- sum(model_df$weights)
  out
}

cfg <- project_config$robustness$overlap_weighting
treat_var <- cfg$treatment_var
balance_vars <- cfg$balance_vars

if (!requireNamespace("nnet", quietly = TRUE)) {
  stop("Package 'nnet' is required for overlap weighting sensitivity analysis.", call. = FALSE)
}

ow_df <- data.frame(
  S_ID = analysis_data$S_ID,
  treatment = factor(analysis_data[[treat_var]])
)
ow_df$S_AGE <- as.numeric(analysis_data$S_AGE)
ow_df$S_PTGENDER <- factor(analysis_data$S_PTGENDER)
ow_df <- ow_df[stats::complete.cases(ow_df), , drop = FALSE]

ps_formula <- stats::as.formula("treatment ~ S_AGE + factor(S_PTGENDER)")
ps_fit <- nnet::multinom(ps_formula, data = ow_df, trace = FALSE)
ps_hat <- stats::predict(ps_fit, type = "probs")
if (is.null(dim(ps_hat))) {
  ps_hat <- matrix(ps_hat, ncol = length(levels(ow_df$treatment)))
}
colnames(ps_hat) <- levels(ow_df$treatment)
ps_hat <- pmax(ps_hat, 1e-6)
obs_group <- as.character(ow_df$treatment)
obs_prob <- ps_hat[cbind(seq_len(nrow(ps_hat)), match(obs_group, colnames(ps_hat)))]
harmonic_overlap <- 1 / rowSums(1 / ps_hat)
overlap_weight <- harmonic_overlap / obs_prob
overlap_weight <- overlap_weight / mean(overlap_weight)

weight_df <- data.frame(
  S_ID = ow_df$S_ID,
  treatment = obs_group,
  overlap_weight = overlap_weight,
  stringsAsFactors = FALSE
)

analysis_weighted <- merge(analysis_data, weight_df[, c("S_ID", "overlap_weight")], by = "S_ID", all.x = TRUE, sort = FALSE)

pair_map <- list(
  "CN vs MCI" = c("CN", "MCI"),
  "MCI vs AD" = c("MCI", "AD"),
  "CN vs AD" = c("CN", "AD")
)

balance_rows <- list()
for (pair_name in names(pair_map)) {
  groups <- pair_map[[pair_name]]
  pair_df <- analysis_weighted[as.character(analysis_weighted[[treat_var]]) %in% groups, , drop = FALSE]
  if (nrow(pair_df) == 0) next
  pair_df$treat_pair <- factor(as.character(pair_df[[treat_var]]), levels = groups)
  pair_df$sex_female <- ifelse(pair_df$S_PTGENDER == 2, 1, 0)

  balance_rows[[length(balance_rows) + 1]] <- data.frame(
    pair = pair_name,
    sample = "Unweighted",
    variable = "Age",
    smd = compute_smd_numeric_weighted(pair_df$S_AGE, pair_df$treat_pair, rep(1, nrow(pair_df))),
    stringsAsFactors = FALSE
  )
  balance_rows[[length(balance_rows) + 1]] <- data.frame(
    pair = pair_name,
    sample = "Unweighted",
    variable = "Female sex",
    smd = compute_smd_binary_weighted(pair_df$sex_female, pair_df$treat_pair, rep(1, nrow(pair_df))),
    stringsAsFactors = FALSE
  )
  balance_rows[[length(balance_rows) + 1]] <- data.frame(
    pair = pair_name,
    sample = "Overlap weighted",
    variable = "Age",
    smd = compute_smd_numeric_weighted(pair_df$S_AGE, pair_df$treat_pair, pair_df$overlap_weight),
    stringsAsFactors = FALSE
  )
  balance_rows[[length(balance_rows) + 1]] <- data.frame(
    pair = pair_name,
    sample = "Overlap weighted",
    variable = "Female sex",
    smd = compute_smd_binary_weighted(pair_df$sex_female, pair_df$treat_pair, pair_df$overlap_weight),
    stringsAsFactors = FALSE
  )
}

balance_df <- bind_rows_fill(balance_rows)
balance_path <- file.path(result_summary_dir, "overlap_weighting_balance.csv")
write_csv_utf8(balance_df, balance_path, row.names = FALSE)

weight_summary <- do.call(
  rbind,
  lapply(split(weight_df$overlap_weight, weight_df$treatment), function(x) {
    data.frame(
      n = length(x),
      min = min(x),
      p05 = unname(stats::quantile(x, 0.05)),
      median = stats::median(x),
      mean = mean(x),
      p95 = unname(stats::quantile(x, 0.95)),
      max = max(x),
      stringsAsFactors = FALSE
    )
  })
)
weight_summary$treatment <- rownames(weight_summary)
rownames(weight_summary) <- NULL
weight_summary <- weight_summary[, c("treatment", "n", "min", "p05", "median", "mean", "p95", "max")]
weight_summary_path <- file.path(result_summary_dir, "overlap_weighting_weight_summary.csv")
write_csv_utf8(weight_summary, weight_summary_path, row.names = FALSE)

weighted_model_rows <- list()
comparison_rows <- list()
for (model_cfg in cfg$models) {
  weighted_fit <- fit_weighted_linear_model_base(
    data = analysis_weighted,
    outcome = model_cfg$outcome,
    exposure = model_cfg$exposure,
    covariates = cfg$adjustment_covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan,
    weights_var = "overlap_weight"
  )
  weighted_fit$model_name <- model_cfg$name
  weighted_model_rows[[model_cfg$name]] <- weighted_fit

  unadjusted_fit <- fit_linear_model_base(
    data = analysis_data,
    outcome = model_cfg$outcome,
    exposure = model_cfg$exposure,
    covariates = NULL,
    factor_vars = NULL,
    transformation_table = transformation_plan
  )

  unweighted_fit <- fit_linear_model_base(
    data = analysis_data,
    outcome = model_cfg$outcome,
    exposure = model_cfg$exposure,
    covariates = cfg$adjustment_covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )

  unweighted_row <- unweighted_fit[unweighted_fit$term == "exposure", , drop = FALSE]
  unadjusted_row <- unadjusted_fit[unadjusted_fit$term == "exposure", , drop = FALSE]
  weighted_row <- weighted_fit[weighted_fit$term == "exposure", , drop = FALSE]
  scale_factor_unadjusted <- compute_model_scale_factor(
    data = analysis_data,
    outcome = model_cfg$outcome,
    exposure = model_cfg$exposure,
    covariates = NULL,
    factor_vars = NULL,
    transformation_table = transformation_plan
  )
  scale_factor_unweighted <- compute_model_scale_factor(
    data = analysis_data,
    outcome = model_cfg$outcome,
    exposure = model_cfg$exposure,
    covariates = cfg$adjustment_covariates,
    factor_vars = project_config$variables$categorical_covariates,
    transformation_table = transformation_plan
  )
  if (nrow(unadjusted_row) == 1) {
    unadjusted_row$std_conf.low <- if (is.finite(scale_factor_unadjusted)) unadjusted_row$conf.low * scale_factor_unadjusted else NA_real_
    unadjusted_row$std_conf.high <- if (is.finite(scale_factor_unadjusted)) unadjusted_row$conf.high * scale_factor_unadjusted else NA_real_
    unadjusted_row$model_name <- model_cfg$name
    unadjusted_row$analysis_type <- "Unadjusted"
    comparison_rows[[paste0(model_cfg$name, "_unadjusted")]] <- unadjusted_row
  }
  if (nrow(unweighted_row) == 1) {
    unweighted_row$std_conf.low <- if (is.finite(scale_factor_unweighted)) unweighted_row$conf.low * scale_factor_unweighted else NA_real_
    unweighted_row$std_conf.high <- if (is.finite(scale_factor_unweighted)) unweighted_row$conf.high * scale_factor_unweighted else NA_real_
    unweighted_row$model_name <- model_cfg$name
    unweighted_row$analysis_type <- "Primary adjusted"
    comparison_rows[[paste0(model_cfg$name, "_unweighted")]] <- unweighted_row
  }
  if (nrow(weighted_row) == 1) {
    weighted_row$model_name <- model_cfg$name
    weighted_row$analysis_type <- "Overlap weighted"
    comparison_rows[[paste0(model_cfg$name, "_weighted")]] <- weighted_row
  }
}

weighted_models_df <- bind_rows_fill(weighted_model_rows)
weighted_models_path <- file.path(result_summary_dir, "overlap_weighting_key_models.csv")
write_csv_utf8(weighted_models_df, weighted_models_path, row.names = FALSE)

comparison_df <- bind_rows_fill(comparison_rows)
comparison_path <- file.path(result_summary_dir, "overlap_weighting_compare_unweighted.csv")
write_csv_utf8(comparison_df, comparison_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "06e_overlap_weighting_sensitivity",
  output_files = c(balance_path, weight_summary_path, weighted_models_path, comparison_path),
  note = "Completed overlap weighting sensitivity analysis to rebalance age and sex across diagnosis groups in the full sample and refit key main-effect models.",
  summary_dir = result_summary_dir
)
