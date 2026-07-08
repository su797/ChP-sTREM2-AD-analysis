format_continuous_summary <- function(x, method = c("parametric", "nonparametric"), digits = 2) {
  method <- match.arg(method)
  x <- as.numeric(x)
  x <- x[!is.na(x)]

  if (length(x) == 0) {
    return(NA_character_)
  }

  if (method == "parametric") {
    sprintf("%0.*f +/- %0.*f", digits, mean(x), digits, stats::sd(x))
  } else {
    q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    sprintf("%0.*f [%0.*f, %0.*f]", digits, q[[2]], digits, q[[1]], digits, q[[3]])
  }
}

run_continuous_group_comparison <- function(data, variable, group_var, transformation_table = NULL, alpha = 0.05, digits = 2) {
  analysis_var <- resolve_analysis_var(variable, transformation_table)
  model_data <- data.frame(
    x = as.numeric(data[[analysis_var]]),
    group = factor(data[[group_var]])
  )
  model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]

  split_x <- split(model_data$x, model_data$group)
  group_normality <- vapply(split_x, function(x) {
    if (length(x) < 3) {
      return(FALSE)
    }
    test_p <- tryCatch(stats::shapiro.test(x)[["p.value"]], error = function(e) NA_real_)
    !is.na(test_p) && test_p >= alpha
  }, FUN.VALUE = logical(1))

  use_parametric <- all(group_normality)

  overall_test <- if (use_parametric) {
    fit <- stats::aov(x ~ group, data = model_data)
    data.frame(
      variable = variable,
      analysis_var = analysis_var,
      method = "anova",
      statistic = summary(fit)[[1]][["F value"]][1],
      p_value = summary(fit)[[1]][["Pr(>F)"]][1],
      stringsAsFactors = FALSE
    )
  } else {
    test_obj <- stats::kruskal.test(x ~ group, data = model_data)
    data.frame(
      variable = variable,
      analysis_var = analysis_var,
      method = "kruskal.test",
      statistic = unname(test_obj$statistic),
      p_value = test_obj$p.value,
      stringsAsFactors = FALSE
    )
  }

  summary_table <- do.call(
    rbind,
    lapply(levels(model_data$group), function(group_name) {
      group_x <- model_data$x[model_data$group == group_name]
      data.frame(
        variable = variable,
        analysis_var = analysis_var,
        group = group_name,
        n = sum(!is.na(group_x)),
        summary = format_continuous_summary(
          group_x,
          method = if (use_parametric) "parametric" else "nonparametric",
          digits = digits
        ),
        stringsAsFactors = FALSE
      )
    })
  )

  pairwise_table <- if (use_parametric) {
    fit <- stats::aov(x ~ group, data = model_data)
    pair_res <- stats::TukeyHSD(fit)$group
    data.frame(
      variable = variable,
      analysis_var = analysis_var,
      comparison = rownames(pair_res),
      estimate = pair_res[, "diff"],
      conf_low = pair_res[, "lwr"],
      conf_high = pair_res[, "upr"],
      p_value = pair_res[, "p adj"],
      method = "TukeyHSD",
      stringsAsFactors = FALSE
    )
  } else {
    pair_res <- stats::pairwise.wilcox.test(model_data$x, model_data$group, p.adjust.method = "BH")
    mat <- as.data.frame(as.table(pair_res$p.value), stringsAsFactors = FALSE)
    mat <- mat[!is.na(mat$Freq), , drop = FALSE]
    data.frame(
      variable = variable,
      analysis_var = analysis_var,
      comparison = paste(mat$Var1, "vs", mat$Var2),
      estimate = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      p_value = mat$Freq,
      method = "pairwise.wilcox",
      stringsAsFactors = FALSE
    )
  }

  list(
    overall = overall_test,
    summary = summary_table,
    pairwise = pairwise_table
  )
}

run_categorical_group_comparison <- function(data, variable, group_var) {
  model_data <- data.frame(
    x = factor(data[[variable]]),
    group = factor(data[[group_var]])
  )
  model_data <- model_data[stats::complete.cases(model_data), , drop = FALSE]

  tab <- table(model_data$x, model_data$group)
  chisq_obj <- suppressWarnings(stats::chisq.test(tab))
  overall_test <- data.frame(
    variable = variable,
    analysis_var = variable,
    method = "chisq.test",
    statistic = unname(chisq_obj$statistic),
    p_value = chisq_obj$p.value,
    stringsAsFactors = FALSE
  )

  summary_table <- data.frame(
    variable = variable,
    level = rep(rownames(tab), times = ncol(tab)),
    group = rep(colnames(tab), each = nrow(tab)),
    n = as.integer(tab),
    pct = round(as.vector(prop.table(tab, margin = 2)) * 100, 2),
    stringsAsFactors = FALSE
  )

  group_levels <- levels(model_data$group)
  pairwise_rows <- list()
  idx <- 1
  if (length(group_levels) >= 2) {
    combs <- utils::combn(group_levels, 2, simplify = FALSE)
    for (pair in combs) {
      sub_data <- model_data[model_data$group %in% pair, , drop = FALSE]
      sub_tab <- table(sub_data$x, sub_data$group)
      sub_test <- suppressWarnings(stats::chisq.test(sub_tab))
      pairwise_rows[[idx]] <- data.frame(
        variable = variable,
        analysis_var = variable,
        comparison = paste(pair, collapse = " vs "),
        estimate = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_,
        p_value = sub_test$p.value,
        method = "pairwise.chisq",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  list(
    overall = overall_test,
    summary = summary_table,
    pairwise = do.call(rbind, pairwise_rows)
  )
}

run_adjusted_linear_group_model <- function(data, outcome, group_var, covariates = NULL, factor_vars = NULL, transformation_table = NULL) {
  outcome_var <- resolve_analysis_var(outcome, transformation_table)
  covariate_vars <- resolve_analysis_vars(covariates, transformation_table)

  model_df <- data.frame(outcome = as.numeric(data[[outcome_var]]))
  model_df$group <- factor(data[[group_var]])

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
  rhs <- c("group", covariates)
  fit <- stats::lm(stats::as.formula(paste("outcome ~", paste(rhs, collapse = " + "))), data = model_df)
  fit_anova <- stats::anova(fit)

  data.frame(
    outcome = outcome,
    analysis_var = outcome_var,
    method = "adjusted_lm",
    n = nrow(model_df),
    group_p_value = fit_anova[["Pr(>F)"]][1],
    stringsAsFactors = FALSE
  )
}
