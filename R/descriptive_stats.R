summarize_continuous <- function(data, vars, group_var = NULL, digits = 2) {
  calc_one <- function(x) {
    data.frame(
      n = sum(!is.na(x)),
      mean = round(mean(x, na.rm = TRUE), digits),
      sd = round(stats::sd(x, na.rm = TRUE), digits),
      median = round(stats::median(x, na.rm = TRUE), digits),
      q1 = round(stats::quantile(x, 0.25, na.rm = TRUE), digits),
      q3 = round(stats::quantile(x, 0.75, na.rm = TRUE), digits),
      min = round(min(x, na.rm = TRUE), digits),
      max = round(max(x, na.rm = TRUE), digits)
    )
  }

  if (is.null(group_var)) {
    out <- do.call(
      rbind,
      lapply(vars, function(v) {
        cbind(variable = v, group = "overall", calc_one(data[[v]]))
      })
    )
    rownames(out) <- NULL
    return(out)
  }

  split_data <- split(data, data[[group_var]])
  out <- do.call(
    rbind,
    lapply(names(split_data), function(g) {
      do.call(
        rbind,
        lapply(vars, function(v) {
          cbind(variable = v, group = g, calc_one(split_data[[g]][[v]]))
        })
      )
    })
  )
  rownames(out) <- NULL
  out
}

summarize_categorical <- function(data, vars, group_var = NULL) {
  summarize_one <- function(df, var_name, group_name) {
    tab <- table(df[[var_name]], useNA = "ifany")
    prop <- prop.table(tab)
    data.frame(
      variable = var_name,
      group = group_name,
      level = names(tab),
      n = as.integer(tab),
      pct = round(as.numeric(prop) * 100, 2),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(group_var)) {
    return(do.call(rbind, lapply(vars, summarize_one, df = data, group_name = "overall")))
  }

  split_data <- split(data, data[[group_var]])
  out <- do.call(
    rbind,
    lapply(names(split_data), function(g) {
      do.call(rbind, lapply(vars, summarize_one, df = split_data[[g]], group_name = g))
    })
  )
  rownames(out) <- NULL
  out
}

compare_groups_ttest <- function(data, continuous_vars, group_var) {
  do.call(
    rbind,
    lapply(continuous_vars, function(v) {
      valid_idx <- !is.na(data[[v]]) & !is.na(data[[group_var]])
      sub_data <- data[valid_idx, , drop = FALSE]
      n_groups <- length(unique(sub_data[[group_var]]))
      formula_obj <- stats::as.formula(sprintf("%s ~ %s", v, group_var))

      if (n_groups == 2) {
        test_obj <- stats::t.test(formula_obj, data = sub_data)
        method_name <- "t.test"
      } else {
        fit_obj <- stats::aov(formula_obj, data = sub_data)
        fit_sum <- summary(fit_obj)[[1]]
        test_obj <- list(
          statistic = fit_sum[["F value"]][1],
          p.value = fit_sum[["Pr(>F)"]][1]
        )
        method_name <- "anova"
      }

      data.frame(
        variable = v,
        method = method_name,
        statistic = unname(test_obj$statistic),
        p_value = test_obj$p.value,
        stringsAsFactors = FALSE
      )
    })
  )
}

compare_groups_chisq <- function(data, categorical_vars, group_var) {
  do.call(
    rbind,
    lapply(categorical_vars, function(v) {
      tab <- table(data[[v]], data[[group_var]])
      test_obj <- suppressWarnings(stats::chisq.test(tab))
      data.frame(
        variable = v,
        method = "chisq.test",
        statistic = unname(test_obj$statistic),
        p_value = test_obj$p.value,
        stringsAsFactors = FALSE
      )
    })
  )
}
