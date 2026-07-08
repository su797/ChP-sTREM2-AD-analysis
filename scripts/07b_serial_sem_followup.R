source(file.path(getwd(), "00_setup.R"))

local_r_lib <- normalizePath(file.path(repo_root, ".r_libs"), winslash = "/", mustWork = FALSE)
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop("Package 'lavaan' is required for serial SEM follow-up analysis.", call. = FALSE)
}

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
cognition_map <- read_project_data(cognition_map_path)

safe_scale_serial <- function(x) {
  x <- as.numeric(x)
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(x - x_mean)
  }
  (x - x_mean) / x_sd
}

resolve_adjustment_serial <- function(cognition_name, sem_config) {
  adjustment <- sem_config$cognition_adjustments[[cognition_name]]
  if (is.null(adjustment)) {
    return(list(
      residual_covariances = list(),
      note = "No additional measurement-model adjustment was configured."
    ))
  }

  if (is.null(adjustment$residual_covariances)) {
    adjustment$residual_covariances <- list()
  }
  if (is.null(adjustment$note) || !nzchar(adjustment$note)) {
    adjustment$note <- "Model-specific measurement adjustment was applied."
  }
  adjustment
}

format_residual_covariances_serial <- function(residual_pairs) {
  if (length(residual_pairs) == 0) {
    return("None")
  }

  pieces <- vapply(residual_pairs, function(pair) {
    paste(pair, collapse = " ~~ ")
  }, FUN.VALUE = character(1))
  paste(pieces, collapse = "; ")
}

evaluate_sem_fit_serial <- function(fit_measures) {
  cfi <- unname(fit_measures[["cfi"]])
  tli <- unname(fit_measures[["tli"]])
  rmsea <- unname(fit_measures[["rmsea"]])
  srmr <- unname(fit_measures[["srmr"]])

  available <- all(!is.na(c(cfi, tli, rmsea, srmr)))
  if (!available) {
    return(list(
      fit_indices_available = FALSE,
      model_reasonable = "Model estimated successfully, but one or more global fit indices are unavailable or not informative.",
      fit_note = "Fit indices were not stably available for this serial mediation model."
    ))
  }

  if (cfi >= 0.95 && tli >= 0.95 && rmsea <= 0.06 && srmr <= 0.08) {
    return(list(
      fit_indices_available = TRUE,
      model_reasonable = "Global fit is good.",
      fit_note = "CFI/TLI >= 0.95 and RMSEA/SRMR are within common recommended thresholds."
    ))
  }

  if (cfi >= 0.90 && tli >= 0.90 && rmsea <= 0.08 && srmr <= 0.10) {
    return(list(
      fit_indices_available = TRUE,
      model_reasonable = "Global fit is acceptable but not ideal.",
      fit_note = "The serial model passes relaxed practical thresholds, but interpretation should remain cautious."
    ))
  }

  list(
    fit_indices_available = TRUE,
    model_reasonable = "Global fit is weak.",
    fit_note = "At least one of CFI/TLI/RMSEA/SRMR falls outside common acceptable thresholds."
  )
}

extract_defined_parameter_serial <- function(pe, lhs_name) {
  row <- pe[pe$lhs == lhs_name & pe$op == ":=", , drop = FALSE]
  if (nrow(row) == 0) {
    return(data.frame(
      estimate = NA_real_,
      std_estimate = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  row <- row[1, , drop = FALSE]
  data.frame(
    estimate = row$est[[1]],
    std_estimate = row$std.all[[1]],
    p_value = row$pvalue[[1]],
    conf_low = row$ci.lower[[1]],
    conf_high = row$ci.upper[[1]],
    stringsAsFactors = FALSE
  )
}

create_serial_sem_dataset <- function(data, cognition_name, cognition_map, x_var, mediator_1_var, mediator_2_var, covariates, sem_config, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  if (nrow(component_rows) == 0) {
    stop(paste("No cognition indicators found for", cognition_name), call. = FALSE)
  }

  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_1_analysis <- resolve_analysis_var(mediator_1_var, transformation_table)
  mediator_2_analysis <- resolve_analysis_var(mediator_2_var, transformation_table)

  sem_data <- data.frame(
    x = as.numeric(data[[x_analysis]]),
    mediator_1 = as.numeric(data[[mediator_1_analysis]]),
    mediator_2 = as.numeric(data[[mediator_2_analysis]])
  )

  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_source <- component_rows$analysis_var[[i]]
    indicator_direction <- as.numeric(component_rows$direction[[i]])
    indicator_values <- indicator_direction * as.numeric(data[[indicator_source]])
    if (isTRUE(sem_config$indicator_standardization)) {
      indicator_values <- safe_scale_serial(indicator_values)
    }
    sem_data[[indicator_name]] <- indicator_values
    indicator_names <- c(indicator_names, indicator_name)
  }

  for (cov_name in covariates) {
    cov_analysis <- resolve_analysis_var(cov_name, transformation_table)
    sem_data[[cov_name]] <- as.numeric(data[[cov_analysis]])
  }

  list(
    data = sem_data,
    indicator_names = indicator_names,
    component_rows = component_rows
  )
}

build_serial_mediation_model <- function(indicator_names, component_rows, covariates, adjustment) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL

  mediator_1_line <- if (is.null(covariate_rhs)) {
    "mediator_1 ~ a1*x"
  } else {
    paste("mediator_1 ~ a1*x +", covariate_rhs)
  }

  mediator_2_line <- if (is.null(covariate_rhs)) {
    "mediator_2 ~ d*mediator_1 + a2*x"
  } else {
    paste("mediator_2 ~ d*mediator_1 + a2*x +", covariate_rhs)
  }

  outcome_line <- if (is.null(covariate_rhs)) {
    "Cog ~ b1*mediator_1 + b2*mediator_2 + c_prime*x"
  } else {
    paste("Cog ~ b1*mediator_1 + b2*mediator_2 + c_prime*x +", covariate_rhs)
  }

  indicator_lookup <- setNames(indicator_names, component_rows$component)
  residual_lines <- character(0)
  for (pair in adjustment$residual_covariances) {
    if (length(pair) != 2) {
      next
    }
    lhs <- indicator_lookup[[pair[[1]]]]
    rhs <- indicator_lookup[[pair[[2]]]]
    if (is.null(lhs) || is.null(rhs)) {
      next
    }
    residual_lines <- c(residual_lines, paste(lhs, "~~", rhs))
  }

  paste(
    c(
      measurement_line,
      mediator_1_line,
      mediator_2_line,
      outcome_line,
      residual_lines,
      "indirect_m1 := a1*b1",
      "indirect_m2 := a2*b2",
      "serial_indirect := a1*d*b2",
      "total_indirect := indirect_m1 + indirect_m2 + serial_indirect",
      "c_total := c_prime + total_indirect",
      "serial_share := serial_indirect / c_total"
    ),
    collapse = "\n"
  )
}

extract_regression_row_serial <- function(pe, lhs, rhs, path_name) {
  row <- pe[pe$lhs == lhs & pe$op == "~" & pe$rhs == rhs, , drop = FALSE]
  if (nrow(row) == 0) {
    return(data.frame(
      path = path_name,
      estimate = NA_real_,
      std_estimate = NA_real_,
      p_value = NA_real_,
      conf_low = NA_real_,
      conf_high = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  row <- row[1, , drop = FALSE]
  data.frame(
    path = path_name,
    estimate = row$est[[1]],
    std_estimate = row$std.all[[1]],
    p_value = row$pvalue[[1]],
    conf_low = row$ci.lower[[1]],
    conf_high = row$ci.upper[[1]],
    stringsAsFactors = FALSE
  )
}

fit_serial_sem_result <- function(data, group_name, cognition_name, serial_model, sem_config, cognition_map, transformation_plan, covariates, covariate_text) {
  adjustment <- resolve_adjustment_serial(cognition_name, sem_config)
  sem_input <- create_serial_sem_dataset(
    data = data,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = serial_model$x,
    mediator_1_var = serial_model$mediator_1,
    mediator_2_var = serial_model$mediator_2,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_plan
  )

  lavaan_model <- build_serial_mediation_model(
    indicator_names = sem_input$indicator_names,
    component_rows = sem_input$component_rows,
    covariates = covariates,
    adjustment = adjustment
  )

  fit_object <- lavaan::sem(
    model = lavaan_model,
    data = sem_input$data,
    estimator = "MLR",
    missing = "fiml",
    fixed.x = FALSE,
    meanstructure = TRUE,
    std.lv = TRUE
  )

  fit_measures <- lavaan::fitMeasures(fit_object, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
  fit_eval <- evaluate_sem_fit_serial(fit_measures)
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)

  path_table <- do.call(
    rbind,
    list(
      extract_regression_row_serial(pe, "mediator_1", "x", "a1"),
      extract_regression_row_serial(pe, "mediator_2", "mediator_1", "d"),
      extract_regression_row_serial(pe, "mediator_2", "x", "a2"),
      extract_regression_row_serial(pe, "Cog", "mediator_1", "b1"),
      extract_regression_row_serial(pe, "Cog", "mediator_2", "b2"),
      extract_regression_row_serial(pe, "Cog", "x", "c_prime"),
      cbind(path = "indirect_m1", extract_defined_parameter_serial(pe, "indirect_m1"), stringsAsFactors = FALSE),
      cbind(path = "indirect_m2", extract_defined_parameter_serial(pe, "indirect_m2"), stringsAsFactors = FALSE),
      cbind(path = "serial_indirect", extract_defined_parameter_serial(pe, "serial_indirect"), stringsAsFactors = FALSE),
      cbind(path = "total_indirect", extract_defined_parameter_serial(pe, "total_indirect"), stringsAsFactors = FALSE),
      cbind(path = "c_total", extract_defined_parameter_serial(pe, "c_total"), stringsAsFactors = FALSE),
      cbind(path = "serial_share", extract_defined_parameter_serial(pe, "serial_share"), stringsAsFactors = FALSE)
    )
  )

  get_metric <- function(path_name, metric = "estimate") {
    row <- path_table[path_table$path == path_name, , drop = FALSE]
    if (nrow(row) == 0) {
      return(NA_real_)
    }
    row[[metric]][[1]]
  }

  serial_est <- get_metric("serial_indirect")
  direct_est <- get_metric("c_prime")
  total_est <- get_metric("c_total")
  serial_p <- get_metric("serial_indirect", "p_value")

  summary_row <- data.frame(
    group = group_name,
    sem_model = serial_model$name,
    x_var = serial_model$x,
    mediator_1_var = serial_model$mediator_1,
    mediator_2_var = serial_model$mediator_2,
    cognition_model = cognition_name,
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    a1 = get_metric("a1"),
    a1_p = get_metric("a1", "p_value"),
    d = get_metric("d"),
    d_p = get_metric("d", "p_value"),
    a2 = get_metric("a2"),
    a2_p = get_metric("a2", "p_value"),
    b1 = get_metric("b1"),
    b1_p = get_metric("b1", "p_value"),
    b2 = get_metric("b2"),
    b2_p = get_metric("b2", "p_value"),
    direct = direct_est,
    direct_p = get_metric("c_prime", "p_value"),
    indirect_m1 = get_metric("indirect_m1"),
    indirect_m1_p = get_metric("indirect_m1", "p_value"),
    indirect_m2 = get_metric("indirect_m2"),
    indirect_m2_p = get_metric("indirect_m2", "p_value"),
    serial_indirect = serial_est,
    serial_indirect_p = serial_p,
    total_indirect = get_metric("total_indirect"),
    total_indirect_p = get_metric("total_indirect", "p_value"),
    total = total_est,
    total_p = get_metric("c_total", "p_value"),
    serial_share_pct = if (is.na(get_metric("serial_share"))) NA_real_ else get_metric("serial_share") * 100,
    serial_opposite_to_direct = !is.na(serial_est) && !is.na(direct_est) && serial_est * direct_est < 0,
    serial_supported = !is.na(serial_p) && serial_p < 0.05,
    stringsAsFactors = FALSE
  )

  fit_row <- data.frame(
    group = group_name,
    sem_model = serial_model$name,
    cognition_model = cognition_name,
    x_var = serial_model$x,
    mediator_1_var = serial_model$mediator_1,
    mediator_2_var = serial_model$mediator_2,
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    cfi = fit_measures[["cfi"]],
    tli = fit_measures[["tli"]],
    rmsea = fit_measures[["rmsea"]],
    srmr = fit_measures[["srmr"]],
    chisq = fit_measures[["chisq"]],
    df = fit_measures[["df"]],
    fit_indices_available = fit_eval$fit_indices_available,
    model_reasonable = fit_eval$model_reasonable,
    fit_note = fit_eval$fit_note,
    residual_covariances = format_residual_covariances_serial(adjustment$residual_covariances),
    stringsAsFactors = FALSE
  )

  list(
    summary_row = summary_row,
    fit_row = fit_row,
    path_table = cbind(
      group = group_name,
      sem_model = serial_model$name,
      cognition_model = cognition_name,
      x_var = serial_model$x,
      mediator_1_var = serial_model$mediator_1,
      mediator_2_var = serial_model$mediator_2,
      covariates = covariate_text,
      path_table,
      stringsAsFactors = FALSE
    )
  )
}

create_multigroup_serial_sem_dataset <- function(data, group_var, group_levels, cognition_name, cognition_map, x_var, mediator_1_var, mediator_2_var, covariates, sem_config, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  if (nrow(component_rows) == 0) {
    stop(paste("No cognition indicators found for", cognition_name), call. = FALSE)
  }

  data <- data[data[[group_var]] %in% group_levels, , drop = FALSE]
  data[[group_var]] <- factor(data[[group_var]], levels = group_levels)

  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_1_analysis <- resolve_analysis_var(mediator_1_var, transformation_table)
  mediator_2_analysis <- resolve_analysis_var(mediator_2_var, transformation_table)

  sem_data <- data.frame(
    .group = data[[group_var]],
    x = as.numeric(data[[x_analysis]]),
    mediator_1 = as.numeric(data[[mediator_1_analysis]]),
    mediator_2 = as.numeric(data[[mediator_2_analysis]])
  )

  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_source <- component_rows$analysis_var[[i]]
    indicator_direction <- as.numeric(component_rows$direction[[i]])
    indicator_values <- indicator_direction * as.numeric(data[[indicator_source]])
    if (isTRUE(sem_config$indicator_standardization)) {
      indicator_values <- safe_scale_serial(indicator_values)
    }
    sem_data[[indicator_name]] <- indicator_values
    indicator_names <- c(indicator_names, indicator_name)
  }

  for (cov_name in covariates) {
    cov_analysis <- resolve_analysis_var(cov_name, transformation_table)
    sem_data[[cov_name]] <- as.numeric(data[[cov_analysis]])
  }

  list(
    data = sem_data,
    indicator_names = indicator_names,
    component_rows = component_rows
  )
}

build_multigroup_serial_lavaan_model <- function(indicator_names, component_rows, covariates, adjustment, group_levels) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  group_keys <- make.names(group_levels)

  label_vector <- function(prefix) {
    paste0("c(", paste0(prefix, "_", group_keys, collapse = ", "), ")*")
  }

  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL

  mediator_1_line <- if (is.null(covariate_rhs)) {
    paste("mediator_1 ~", paste0(label_vector("a1"), "x"))
  } else {
    paste("mediator_1 ~", paste0(label_vector("a1"), "x + ", covariate_rhs))
  }

  mediator_2_line <- if (is.null(covariate_rhs)) {
    paste("mediator_2 ~", paste0(label_vector("d"), "mediator_1 + ", label_vector("a2"), "x"))
  } else {
    paste("mediator_2 ~", paste0(label_vector("d"), "mediator_1 + ", label_vector("a2"), "x + ", covariate_rhs))
  }

  outcome_line <- if (is.null(covariate_rhs)) {
    paste("Cog ~", paste0(label_vector("b1"), "mediator_1 + ", label_vector("b2"), "mediator_2 + ", label_vector("cprime"), "x"))
  } else {
    paste("Cog ~", paste0(label_vector("b1"), "mediator_1 + ", label_vector("b2"), "mediator_2 + ", label_vector("cprime"), "x + ", covariate_rhs))
  }

  indicator_lookup <- setNames(indicator_names, component_rows$component)
  residual_lines <- character(0)
  for (pair in adjustment$residual_covariances) {
    if (length(pair) != 2) {
      next
    }
    lhs <- indicator_lookup[[pair[[1]]]]
    rhs <- indicator_lookup[[pair[[2]]]]
    if (is.null(lhs) || is.null(rhs)) {
      next
    }
    residual_lines <- c(residual_lines, paste(lhs, "~~", rhs))
  }

  defined_lines <- unlist(lapply(group_keys, function(key) {
    c(
      paste0("indirect_m1_", key, " := a1_", key, "*b1_", key),
      paste0("indirect_m2_", key, " := a2_", key, "*b2_", key),
      paste0("serial_indirect_", key, " := a1_", key, "*d_", key, "*b2_", key),
      paste0("total_indirect_", key, " := indirect_m1_", key, " + indirect_m2_", key, " + serial_indirect_", key)
    )
  }))

  paste(
    c(
      measurement_line,
      mediator_1_line,
      mediator_2_line,
      outcome_line,
      residual_lines,
      defined_lines
    ),
    collapse = "\n"
  )
}

build_serial_multigroup_constraints <- function(group_levels, path_name = c("serial_indirect", "d"), scope = c("omnibus", "pairwise")) {
  path_name <- match.arg(path_name)
  scope <- match.arg(scope)
  group_keys <- make.names(group_levels)

  path_expr <- switch(
    path_name,
    serial_indirect = paste0("a1_", group_keys, "*d_", group_keys, "*b2_", group_keys),
    d = paste0("d_", group_keys)
  )

  if (scope == "omnibus") {
    if (length(path_expr) < 2) {
      return(character(0))
    }
    return(vapply(seq(2, length(path_expr)), function(idx) {
      paste(path_expr[[1]], "==", path_expr[[idx]])
    }, FUN.VALUE = character(1)))
  }

  pair_indices <- utils::combn(seq_along(path_expr), 2)
  apply(pair_indices, 2, function(idx) {
    paste(path_expr[[idx[[1]]]], "==", path_expr[[idx[[2]]]])
  })
}

run_serial_wald_test <- function(fit_object, constraints, sem_model_name, cognition_name, path_name, scope, group_pair = NA_character_) {
  if (length(constraints) == 0) {
    return(data.frame())
  }

  test_result <- tryCatch(
    lavaan::lavTestWald(fit_object, constraints = constraints),
    error = function(e) NULL
  )

  if (is.null(test_result)) {
    return(data.frame())
  }

  data.frame(
    sem_model = sem_model_name,
    cognition_model = cognition_name,
    path = path_name,
    scope = scope,
    group_pair = group_pair,
    statistic = unname(test_result$stat[[1]]),
    df = unname(test_result$df[[1]]),
    p_value = unname(test_result$p.value[[1]]),
    significant = !is.na(unname(test_result$p.value[[1]])) && unname(test_result$p.value[[1]]) < 0.05,
    stringsAsFactors = FALSE
  )
}

extract_serial_multigroup_group_estimates <- function(fit_object, group_levels, serial_model, cognition_name, covariate_text) {
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)
  group_keys <- make.names(group_levels)
  out <- list()
  idx <- 1

  for (g in seq_along(group_levels)) {
    serial_row <- pe[pe$lhs == paste0("serial_indirect_", group_keys[[g]]) & pe$op == ":=", , drop = FALSE]
    total_row <- pe[pe$lhs == paste0("total_indirect_", group_keys[[g]]) & pe$op == ":=", , drop = FALSE]
    if (nrow(serial_row) == 0) {
      next
    }
    serial_row <- serial_row[1, , drop = FALSE]
    if (nrow(total_row) > 0) {
      total_row <- total_row[1, , drop = FALSE]
    } else {
      total_row <- NULL
    }

    out[[idx]] <- data.frame(
      sem_model = serial_model$name,
      cognition_model = cognition_name,
      group = group_levels[[g]],
      x_var = serial_model$x,
      mediator_1_var = serial_model$mediator_1,
      mediator_2_var = serial_model$mediator_2,
      effect_type = "serial_indirect",
      estimate = serial_row$est[[1]],
      std_estimate = serial_row$std.all[[1]],
      p_value = serial_row$pvalue[[1]],
      conf_low = serial_row$ci.lower[[1]],
      conf_high = serial_row$ci.upper[[1]],
      total_indirect = if (is.null(total_row)) NA_real_ else total_row$est[[1]],
      total_indirect_p = if (is.null(total_row)) NA_real_ else total_row$pvalue[[1]],
      covariates = covariate_text,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }

  if (length(out) == 0) {
    data.frame()
  } else {
    do.call(rbind, out)
  }
}

fit_serial_multigroup_sem_result <- function(data, group_var, group_levels, cognition_name, serial_model, sem_config, cognition_map, transformation_plan, covariates, covariate_text) {
  adjustment <- resolve_adjustment_serial(cognition_name, sem_config)
  sem_input <- create_multigroup_serial_sem_dataset(
    data = data,
    group_var = group_var,
    group_levels = group_levels,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = serial_model$x,
    mediator_1_var = serial_model$mediator_1,
    mediator_2_var = serial_model$mediator_2,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_plan
  )

  lavaan_model <- build_multigroup_serial_lavaan_model(
    indicator_names = sem_input$indicator_names,
    component_rows = sem_input$component_rows,
    covariates = covariates,
    adjustment = adjustment,
    group_levels = group_levels
  )

  fit_object <- lavaan::sem(
    model = lavaan_model,
    data = sem_input$data,
    group = ".group",
    estimator = "MLR",
    missing = "fiml",
    fixed.x = FALSE,
    meanstructure = TRUE,
    std.lv = TRUE
  )

  fit_measures <- lavaan::fitMeasures(fit_object, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
  fit_eval <- evaluate_sem_fit_serial(fit_measures)
  group_estimates <- extract_serial_multigroup_group_estimates(
    fit_object = fit_object,
    group_levels = group_levels,
    serial_model = serial_model,
    cognition_name = cognition_name,
    covariate_text = covariate_text
  )

  test_rows <- list()
  idx <- 1
  for (path_name in c("serial_indirect", "d")) {
    omnibus_constraints <- build_serial_multigroup_constraints(group_levels, path_name = path_name, scope = "omnibus")
    test_rows[[idx]] <- run_serial_wald_test(
      fit_object = fit_object,
      constraints = omnibus_constraints,
      sem_model_name = serial_model$name,
      cognition_name = cognition_name,
      path_name = path_name,
      scope = "omnibus"
    )
    idx <- idx + 1

    pair_constraints <- build_serial_multigroup_constraints(group_levels, path_name = path_name, scope = "pairwise")
    pair_groups <- utils::combn(group_levels, 2, simplify = FALSE)
    for (pair_idx in seq_along(pair_groups)) {
      test_rows[[idx]] <- run_serial_wald_test(
        fit_object = fit_object,
        constraints = pair_constraints[[pair_idx]],
        sem_model_name = serial_model$name,
        cognition_name = cognition_name,
        path_name = path_name,
        scope = "pairwise",
        group_pair = paste(pair_groups[[pair_idx]], collapse = "_vs_")
      )
      idx <- idx + 1
    }
  }

  fit_row <- data.frame(
    sem_model = serial_model$name,
    cognition_model = cognition_name,
    x_var = serial_model$x,
    mediator_1_var = serial_model$mediator_1,
    mediator_2_var = serial_model$mediator_2,
    groups = paste(group_levels, collapse = ", "),
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    cfi = fit_measures[["cfi"]],
    tli = fit_measures[["tli"]],
    rmsea = fit_measures[["rmsea"]],
    srmr = fit_measures[["srmr"]],
    chisq = fit_measures[["chisq"]],
    df = fit_measures[["df"]],
    fit_indices_available = fit_eval$fit_indices_available,
    model_reasonable = fit_eval$model_reasonable,
    fit_note = fit_eval$fit_note,
    residual_covariances = format_residual_covariances_serial(adjustment$residual_covariances),
    stringsAsFactors = FALSE
  )

  list(
    fit_row = fit_row,
    group_estimates = group_estimates,
    test_rows = do.call(rbind, Filter(function(x) nrow(x) > 0, test_rows))
  )
}

serial_summary_rows <- list()
serial_fit_rows <- list()
serial_path_rows <- list()
serial_idx <- 1

sem_group_var <- project_config$variables$group_label_var
sem_group_order <- unique(c("Overall", unlist(runtime_settings$group_order %||% c("CN", "MCI", "AD"), use.names = FALSE)))
sem_group_data <- setNames(vector("list", length(sem_group_order)), sem_group_order)
sem_group_data[["Overall"]] <- analysis_data
for (group_name in setdiff(sem_group_order, "Overall")) {
  sem_group_data[[group_name]] <- analysis_data[analysis_data[[sem_group_var]] == group_name, , drop = FALSE]
}

for (group_name in names(sem_group_data)) {
  group_data <- sem_group_data[[group_name]]
  if (nrow(group_data) == 0) {
    next
  }

  for (cognition_name in names(project_config$cognition_models)) {
    for (serial_model in project_config$sem$serial_models %||% list()) {
      model_covariates <- unique(c(project_config$variables$covariates, serial_model$extra_covariates %||% character(0)))
      covariate_text <- if (length(model_covariates) == 0) "None" else paste(model_covariates, collapse = " + ")

      serial_result <- fit_serial_sem_result(
        data = group_data,
        group_name = group_name,
        cognition_name = cognition_name,
        serial_model = serial_model,
        sem_config = project_config$sem,
        cognition_map = cognition_map,
        transformation_plan = transformation_plan,
        covariates = model_covariates,
        covariate_text = covariate_text
      )

      serial_summary_rows[[serial_idx]] <- serial_result$summary_row
      serial_fit_rows[[serial_idx]] <- serial_result$fit_row
      serial_path_rows[[serial_idx]] <- serial_result$path_table
      serial_idx <- serial_idx + 1
    }
  }
}

serial_summary <- if (length(serial_summary_rows) > 0) do.call(rbind, serial_summary_rows) else data.frame()
serial_fit <- if (length(serial_fit_rows) > 0) do.call(rbind, serial_fit_rows) else data.frame()
serial_paths <- if (length(serial_path_rows) > 0) do.call(rbind, serial_path_rows) else data.frame()

serial_multigroup_fit_rows <- list()
serial_multigroup_group_rows <- list()
serial_multigroup_test_rows <- list()
serial_mg_idx <- 1
serial_multigroup_groups <- c("CN", "MCI", "AD")

for (cognition_name in names(project_config$cognition_models)) {
  for (serial_model in project_config$sem$serial_models %||% list()) {
    model_covariates <- unique(c(project_config$variables$covariates, serial_model$extra_covariates %||% character(0)))
    covariate_text <- if (length(model_covariates) == 0) "None" else paste(model_covariates, collapse = " + ")

    serial_multigroup_result <- fit_serial_multigroup_sem_result(
      data = analysis_data,
      group_var = sem_group_var,
      group_levels = serial_multigroup_groups,
      cognition_name = cognition_name,
      serial_model = serial_model,
      sem_config = project_config$sem,
      cognition_map = cognition_map,
      transformation_plan = transformation_plan,
      covariates = model_covariates,
      covariate_text = covariate_text
    )

    serial_multigroup_fit_rows[[serial_mg_idx]] <- serial_multigroup_result$fit_row
    if (nrow(serial_multigroup_result$group_estimates) > 0) {
      serial_multigroup_group_rows[[serial_mg_idx]] <- serial_multigroup_result$group_estimates
    }
    if (nrow(serial_multigroup_result$test_rows) > 0) {
      serial_multigroup_test_rows[[serial_mg_idx]] <- serial_multigroup_result$test_rows
    }
    serial_mg_idx <- serial_mg_idx + 1
  }
}

serial_multigroup_fit <- if (length(serial_multigroup_fit_rows) > 0) do.call(rbind, serial_multigroup_fit_rows) else data.frame()
serial_multigroup_groups_df <- if (length(serial_multigroup_group_rows) > 0) do.call(rbind, serial_multigroup_group_rows) else data.frame()
serial_multigroup_tests <- if (length(serial_multigroup_test_rows) > 0) do.call(rbind, serial_multigroup_test_rows) else data.frame()

single_sem_summary_path <- file.path(result_summary_dir, "sem_mediation_summary.csv")
single_sem_summary <- if (file.exists(single_sem_summary_path)) read_project_data(single_sem_summary_path) else data.frame()

serial_coverage <- data.frame()
if (nrow(serial_summary) > 0 && nrow(single_sem_summary) > 0) {
  single_map <- data.frame(
    serial_model = c("PTAU_to_sTREM2_to_ChP_to_Cognition", "TAU_to_sTREM2_to_ChP_to_Cognition"),
    single_sem_model = c("PTAU_to_ChP_to_Cognition", "TAU_to_ChP_to_Cognition"),
    stringsAsFactors = FALSE
  )
  serial_coverage <- merge(serial_summary, single_map, by.x = "sem_model", by.y = "serial_model", all.x = TRUE)
  serial_coverage <- merge(
    serial_coverage,
    single_sem_summary[, c("sem_model", "group", "cognition_model", "indirect", "indirect_p", "direct", "direct_p", "total", "total_p")],
    by.x = c("single_sem_model", "group", "cognition_model"),
    by.y = c("sem_model", "group", "cognition_model"),
    all.x = TRUE,
    suffixes = c("", "_single")
  )
  serial_coverage$serial_share_of_single_indirect_pct <- ifelse(
    is.na(serial_coverage$indirect) | serial_coverage$indirect == 0,
    NA_real_,
    (serial_coverage$serial_indirect / serial_coverage$indirect) * 100
  )
  serial_coverage$supports_sTREM2_explaining_ptau_suppression <- !is.na(serial_coverage$serial_indirect_p) &
    serial_coverage$serial_indirect_p < 0.05 &
    !is.na(serial_coverage$serial_share_of_single_indirect_pct)
}

serial_summary_path <- file.path(result_summary_dir, "sem_serial_mediation_summary.csv")
serial_fit_path <- file.path(result_summary_dir, "sem_serial_model_fit.csv")
serial_paths_path <- file.path(result_summary_dir, "sem_serial_path_coefficients.csv")
serial_multigroup_fit_path <- file.path(result_summary_dir, "sem_serial_multigroup_model_fit.csv")
serial_multigroup_groups_path <- file.path(result_summary_dir, "sem_serial_multigroup_group_estimates.csv")
serial_multigroup_tests_path <- file.path(result_summary_dir, "sem_serial_multigroup_tests.csv")
serial_coverage_path <- file.path(result_summary_dir, "sem_serial_suppression_coverage.csv")
serial_report_path <- file.path(result_report_dir, "SEM_serial_followup_interpretation.md")

append_analysis_log(
  project_root = project_root,
  analysis_name = "07b_serial_sem_followup",
  output_files = c(
    serial_summary_path,
    serial_fit_path,
    serial_paths_path,
    serial_multigroup_fit_path,
    serial_multigroup_groups_path,
    serial_multigroup_tests_path,
    serial_coverage_path
  ),
  note = "Completed serial SEM follow-up models and exported CSV summaries. No submission text is generated in the public release.",
  summary_dir = result_summary_dir
)
