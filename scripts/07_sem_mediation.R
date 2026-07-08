source(file.path(getwd(), "00_setup.R"))

local_r_lib <- normalizePath(file.path(repo_root, ".r_libs"), winslash = "/", mustWork = FALSE)
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop("Package 'lavaan' is required for latent-variable SEM. Please install it before running 07_sem_mediation.R.", call. = FALSE)
}

find_latest_complete_result_dir_sem <- function(project_root) {
  result_root <- file.path(project_root, "results")
  dirs <- list.dirs(result_root, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[grepl("^\\d{8}_\\d{6}$", basename(dirs))]
  dirs <- sort(dirs, decreasing = TRUE)
  for (d in dirs) {
    required <- c(
      file.path(d, "data_clean", "ChpSTREM2AD_analysis_dataset.csv"),
      file.path(d, "summary", "sem_mediation_summary.csv")
    )
    if (all(file.exists(required))) return(d)
  }
  stop("No complete result directory found for SEM analysis.", call. = FALSE)
}

result_dir_use_sem <- dirname(result_summary_dir)
if (!file.exists(analysis_data_path)) {
  result_dir_use_sem <- find_latest_complete_result_dir_sem(project_root)
  result_summary_dir <- file.path(result_dir_use_sem, "summary")
  result_figures_dir <- file.path(result_dir_use_sem, "figures")
  result_tables_dir <- file.path(result_dir_use_sem, "tables")
  result_report_dir <- file.path(result_dir_use_sem, "diagnostics")
  analysis_data_path <- file.path(result_dir_use_sem, "data_clean", "ChpSTREM2AD_analysis_dataset.csv")
  transformation_table_path <- file.path(result_summary_dir, "transformation_plan.csv")
  cognition_map_path <- file.path(result_summary_dir, "cognition_component_map.csv")
}

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
cognition_map <- read_project_data(cognition_map_path)

safe_scale <- function(x) {
  x <- as.numeric(x)
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(x - x_mean)
  }
  (x - x_mean) / x_sd
}

format_residual_covariances <- function(residual_pairs) {
  if (length(residual_pairs) == 0) {
    return("")
  }

  pieces <- vapply(residual_pairs, function(pair) {
    paste(pair, collapse = " ~~ ")
  }, FUN.VALUE = character(1))
  paste(pieces, collapse = "; ")
}

resolve_adjustment <- function(cognition_name, sem_config) {
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

create_sem_dataset <- function(data, cognition_name, cognition_map, x_var, mediator_var, covariates, sem_config, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  if (nrow(component_rows) == 0) {
    stop(paste("No cognition indicators found for", cognition_name), call. = FALSE)
  }

  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_analysis <- resolve_analysis_var(mediator_var, transformation_table)

  sem_data <- data.frame(
    x = as.numeric(data[[x_analysis]]),
    mediator = as.numeric(data[[mediator_analysis]])
  )

  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_source <- component_rows$analysis_var[[i]]
    indicator_direction <- as.numeric(component_rows$direction[[i]])
    indicator_values <- indicator_direction * as.numeric(data[[indicator_source]])
    if (isTRUE(sem_config$indicator_standardization)) {
      indicator_values <- safe_scale(indicator_values)
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

build_cognition_component_map <- function(cognition_models, transformation_table = NULL) {
  do.call(
    rbind,
    lapply(names(cognition_models), function(composite_name) {
      cfg <- cognition_models[[composite_name]]
      do.call(
        rbind,
        lapply(seq_along(cfg$components), function(i) {
          data.frame(
            composite = composite_name,
            component = cfg$components[[i]],
            analysis_var = resolve_analysis_var(cfg$components[[i]], transformation_table),
            direction = cfg$directions[[i]],
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )
}

build_lavaan_model <- function(indicator_names, component_rows, covariates, adjustment) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  mediator_line <- if (is.null(covariate_rhs)) {
    "mediator ~ a*x"
  } else {
    paste("mediator ~ a*x +", covariate_rhs)
  }
  outcome_line <- if (is.null(covariate_rhs)) {
    "Cog ~ b*mediator + c_prime*x"
  } else {
    paste("Cog ~ b*mediator + c_prime*x +", covariate_rhs)
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
      mediator_line,
      outcome_line,
      residual_lines,
      "indirect := a*b",
      "c_total := c_prime + (a*b)",
      "proportion_mediated := indirect / c_total"
    ),
    collapse = "\n"
  )
}

evaluate_sem_fit <- function(fit_measures) {
  cfi <- unname(fit_measures[["cfi"]])
  tli <- unname(fit_measures[["tli"]])
  rmsea <- unname(fit_measures[["rmsea"]])
  srmr <- unname(fit_measures[["srmr"]])

  available <- all(!is.na(c(cfi, tli, rmsea, srmr)))
  if (!available) {
    return(list(
      fit_indices_available = FALSE,
      model_reasonable = "Model estimated successfully, but one or more global fit indices are unavailable or not informative.",
      fit_note = "This usually occurs when the model is close to just-identified or when fit indices cannot be stably computed from the available indicator structure."
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
      fit_note = "The model passes relaxed practical thresholds, but the measurement structure should still be interpreted cautiously."
    ))
  }

  list(
    fit_indices_available = TRUE,
    model_reasonable = "Global fit is weak.",
    fit_note = "At least one of CFI/TLI/RMSEA/SRMR falls outside common acceptable thresholds; consider revising indicators or model structure."
  )
}

extract_sem_path_table <- function(fit_object) {
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)

  extract_row <- function(lhs, op, rhs = NULL, path_name) {
    row <- pe[pe$lhs == lhs & pe$op == op, , drop = FALSE]
    if (!is.null(rhs)) {
      row <- row[row$rhs == rhs, , drop = FALSE]
    }
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

  do.call(
    rbind,
    list(
      extract_row("mediator", "~", "x", "a"),
      extract_row("Cog", "~", "mediator", "b"),
      extract_row("Cog", "~", "x", "c_prime"),
      extract_row("c_total", ":=", NULL, "c_total"),
      extract_row("indirect", ":=", NULL, "indirect"),
      extract_row("proportion_mediated", ":=", NULL, "proportion_mediated")
    )
  )
}

extract_modification_indices <- function(fit_object, sem_model_name, cognition_name, top_n = 10) {
  mi <- tryCatch(
    lavaan::modificationIndices(fit_object, sort. = TRUE),
    error = function(e) NULL
  )
  if (is.null(mi) || nrow(mi) == 0) {
    return(data.frame())
  }

  mi <- mi[mi$op %in% c("~~", "=~", "~"), c("lhs", "op", "rhs", "mi", "epc", "sepc.lv", "sepc.all"), drop = FALSE]
  mi <- mi[order(-mi$mi), , drop = FALSE]
  mi <- head(mi, top_n)
  if (nrow(mi) == 0) {
    return(data.frame())
  }

  mi$sem_model <- sem_model_name
  mi$cognition_model <- cognition_name
  mi
}

build_parallel_mediation_model <- function(indicator_names, component_rows, covariates, adjustment) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  mediator_1_line <- if (is.null(covariate_rhs)) {
    "mediator_1 ~ a1*x"
  } else {
    paste("mediator_1 ~ a1*x +", covariate_rhs)
  }
  mediator_2_line <- if (is.null(covariate_rhs)) {
    "mediator_2 ~ a2*x"
  } else {
    paste("mediator_2 ~ a2*x +", covariate_rhs)
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
      "mediator_1 ~~ mediator_2",
      outcome_line,
      residual_lines,
      "indirect_1 := a1*b1",
      "indirect_2 := a2*b2",
      "indirect_diff := indirect_1 - indirect_2"
    ),
    collapse = "\n"
  )
}

create_parallel_sem_dataset <- function(data, cognition_name, cognition_map, x_var, mediator_a, mediator_b, covariates, sem_config, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  if (nrow(component_rows) == 0) {
    stop(paste("No cognition indicators found for", cognition_name), call. = FALSE)
  }

  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_a_analysis <- resolve_analysis_var(mediator_a, transformation_table)
  mediator_b_analysis <- resolve_analysis_var(mediator_b, transformation_table)

  sem_data <- data.frame(
    x = as.numeric(data[[x_analysis]]),
    mediator_1 = as.numeric(data[[mediator_a_analysis]]),
    mediator_2 = as.numeric(data[[mediator_b_analysis]])
  )

  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_source <- component_rows$analysis_var[[i]]
    indicator_direction <- as.numeric(component_rows$direction[[i]])
    indicator_values <- indicator_direction * as.numeric(data[[indicator_source]])
    if (isTRUE(sem_config$indicator_standardization)) {
      indicator_values <- safe_scale(indicator_values)
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

extract_defined_parameter <- function(fit_object, lhs_name) {
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)
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

build_parallel_multi_mediation_model <- function(indicator_names, component_rows, covariates, adjustment, mediator_count) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL

  mediator_lines <- vapply(seq_len(mediator_count), function(i) {
    lhs <- paste0("mediator_", i)
    if (is.null(covariate_rhs)) {
      paste(lhs, "~", paste0("a", i, "*x"))
    } else {
      paste(lhs, "~", paste0("a", i, "*x + ", covariate_rhs))
    }
  }, FUN.VALUE = character(1))

  outcome_terms <- paste0("b", seq_len(mediator_count), "*mediator_", seq_len(mediator_count))
  outcome_rhs <- paste(c(outcome_terms, "c_prime*x", covariate_rhs), collapse = " + ")
  outcome_line <- paste("Cog ~", outcome_rhs)

  mediator_cov_lines <- character(0)
  if (mediator_count >= 2) {
    mediator_pairs <- utils::combn(seq_len(mediator_count), 2, simplify = FALSE)
    mediator_cov_lines <- vapply(mediator_pairs, function(idx) {
      paste0("mediator_", idx[[1]], " ~~ mediator_", idx[[2]])
    }, FUN.VALUE = character(1))
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

  indirect_names <- paste0("indirect_", seq_len(mediator_count))
  indirect_lines <- vapply(seq_len(mediator_count), function(i) {
    paste0(indirect_names[[i]], " := a", i, "*b", i)
  }, FUN.VALUE = character(1))
  total_indirect_line <- paste0("total_indirect := ", paste(indirect_names, collapse = " + "))
  c_total_line <- "c_total := c_prime + total_indirect"

  paste(
    c(
      measurement_line,
      mediator_lines,
      mediator_cov_lines,
      outcome_line,
      residual_lines,
      indirect_lines,
      total_indirect_line,
      c_total_line
    ),
    collapse = "\n"
  )
}

create_parallel_multi_sem_dataset <- function(data, cognition_name, cognition_map, x_var, mediators, covariates, sem_config, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  if (nrow(component_rows) == 0) {
    stop(paste("No cognition indicators found for", cognition_name), call. = FALSE)
  }

  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_analysis <- resolve_analysis_vars(mediators, transformation_table)

  sem_data <- data.frame(
    x = as.numeric(data[[x_analysis]])
  )

  for (i in seq_along(mediators)) {
    sem_data[[paste0("mediator_", i)]] <- as.numeric(data[[mediator_analysis[[i]]]])
  }

  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_source <- component_rows$analysis_var[[i]]
    indicator_direction <- as.numeric(component_rows$direction[[i]])
    indicator_values <- indicator_direction * as.numeric(data[[indicator_source]])
    if (isTRUE(sem_config$indicator_standardization)) {
      indicator_values <- safe_scale(indicator_values)
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
    component_rows = component_rows,
    mediator_vars = mediators
  )
}

fit_parallel_sem_result <- function(
  data,
  group_name,
  cognition_name,
  parallel_model,
  sem_config,
  cognition_map,
  transformation_plan,
  covariates,
  covariate_text
) {
  adjustment <- resolve_adjustment(cognition_name, sem_config)
  residual_text <- format_residual_covariances(adjustment$residual_covariances)

  sem_input <- create_parallel_multi_sem_dataset(
    data = data,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = parallel_model$x,
    mediators = parallel_model$mediators,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_plan
  )

  lavaan_model <- build_parallel_multi_mediation_model(
    indicator_names = sem_input$indicator_names,
    component_rows = sem_input$component_rows,
    covariates = covariates,
    adjustment = adjustment,
    mediator_count = length(parallel_model$mediators)
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
  fit_eval <- evaluate_sem_fit(fit_measures)
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)

  get_regression_row <- function(lhs, rhs) {
    row <- pe[pe$lhs == lhs & pe$op == "~" & pe$rhs == rhs, , drop = FALSE]
    if (nrow(row) == 0) {
      return(NULL)
    }
    row[1, , drop = FALSE]
  }

  direct_row <- get_regression_row("Cog", "x")
  total_indirect_row <- pe[pe$lhs == "total_indirect" & pe$op == ":=", , drop = FALSE]
  if (nrow(total_indirect_row) == 0) {
    total_indirect_row <- NULL
  } else {
    total_indirect_row <- total_indirect_row[1, , drop = FALSE]
  }
  total_row <- pe[pe$lhs == "c_total" & pe$op == ":=", , drop = FALSE]
  if (nrow(total_row) == 0) {
    total_row <- NULL
  } else {
    total_row <- total_row[1, , drop = FALSE]
  }

  mediator_rows <- lapply(seq_along(parallel_model$mediators), function(i) {
    a_row <- get_regression_row(paste0("mediator_", i), "x")
    b_row <- get_regression_row("Cog", paste0("mediator_", i))
    indirect_row <- pe[pe$lhs == paste0("indirect_", i) & pe$op == ":=", , drop = FALSE]
    if (nrow(indirect_row) == 0) {
      indirect_row <- NULL
    } else {
      indirect_row <- indirect_row[1, , drop = FALSE]
    }

    direct_est <- if (!is.null(direct_row)) direct_row$est[[1]] else NA_real_
    direct_p <- if (!is.null(direct_row)) direct_row$pvalue[[1]] else NA_real_
    indirect_est <- if (!is.null(indirect_row)) indirect_row$est[[1]] else NA_real_
    indirect_p <- if (!is.null(indirect_row)) indirect_row$pvalue[[1]] else NA_real_
    opposite_direction <- !is.na(direct_est) && !is.na(indirect_est) && direct_est * indirect_est < 0

    data.frame(
      parallel_model = parallel_model$name,
      x_var = parallel_model$x,
      mediator_var = parallel_model$mediators[[i]],
      cognition_model = cognition_name,
      group = group_name,
      covariates = covariate_text,
      n = lavaan::nobs(fit_object),
      a = if (!is.null(a_row)) a_row$est[[1]] else NA_real_,
      a_std = if (!is.null(a_row)) a_row$std.all[[1]] else NA_real_,
      a_p = if (!is.null(a_row)) a_row$pvalue[[1]] else NA_real_,
      b = if (!is.null(b_row)) b_row$est[[1]] else NA_real_,
      b_std = if (!is.null(b_row)) b_row$std.all[[1]] else NA_real_,
      b_p = if (!is.null(b_row)) b_row$pvalue[[1]] else NA_real_,
      indirect = indirect_est,
      indirect_std = if (!is.null(indirect_row)) indirect_row$std.all[[1]] else NA_real_,
      indirect_p = indirect_p,
      indirect_conf_low = if (!is.null(indirect_row)) indirect_row$ci.lower[[1]] else NA_real_,
      indirect_conf_high = if (!is.null(indirect_row)) indirect_row$ci.upper[[1]] else NA_real_,
      direct = direct_est,
      direct_std = if (!is.null(direct_row)) direct_row$std.all[[1]] else NA_real_,
      direct_p = direct_p,
      total_indirect = if (!is.null(total_indirect_row)) total_indirect_row$est[[1]] else NA_real_,
      total_indirect_std = if (!is.null(total_indirect_row)) total_indirect_row$std.all[[1]] else NA_real_,
      total_indirect_p = if (!is.null(total_indirect_row)) total_indirect_row$pvalue[[1]] else NA_real_,
      total = if (!is.null(total_row)) total_row$est[[1]] else NA_real_,
      total_std = if (!is.null(total_row)) total_row$std.all[[1]] else NA_real_,
      total_p = if (!is.null(total_row)) total_row$pvalue[[1]] else NA_real_,
      opposite_direction = opposite_direction,
      mediation_type = ifelse(opposite_direction, "inconsistent_mediation", "consistent_mediation"),
      stringsAsFactors = FALSE
    )
  })

  fit_row <- data.frame(
    parallel_model = parallel_model$name,
    x_var = parallel_model$x,
    mediators = paste(parallel_model$mediators, collapse = " + "),
    cognition_model = cognition_name,
    group = group_name,
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    cfi = unname(fit_measures[["cfi"]]),
    tli = unname(fit_measures[["tli"]]),
    rmsea = unname(fit_measures[["rmsea"]]),
    srmr = unname(fit_measures[["srmr"]]),
    chisq = unname(fit_measures[["chisq"]]),
    df = unname(fit_measures[["df"]]),
    fit_indices_available = fit_eval$fit_indices_available,
    model_reasonable = fit_eval$model_reasonable,
    fit_note = fit_eval$fit_note,
    indicator_standardized = isTRUE(sem_config$indicator_standardization),
    residual_covariances = residual_text,
    optimization_note = adjustment$note,
    stringsAsFactors = FALSE
  )

  list(
    fit_object = fit_object,
    fit_row = fit_row,
    mediator_rows = do.call(rbind, mediator_rows)
  )
}

create_multigroup_parallel_sem_dataset <- function(data, group_var, group_levels, cognition_name, cognition_map, x_var, mediators, covariates, sem_config, transformation_table = NULL) {
  filtered_data <- data[
    !is.na(data[[group_var]]) & data[[group_var]] %in% group_levels,
    ,
    drop = FALSE
  ]

  sem_input <- create_parallel_multi_sem_dataset(
    data = filtered_data,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = x_var,
    mediators = mediators,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_table
  )

  sem_input$data$sem_group <- factor(filtered_data[[group_var]], levels = group_levels)
  sem_input$group_levels <- group_levels
  sem_input$group_keys <- make.names(group_levels)
  sem_input
}

build_multigroup_parallel_lavaan_model <- function(indicator_names, component_rows, covariates, adjustment, group_levels, mediator_count) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  group_keys <- make.names(group_levels)
  multigroup_term <- function(base_name) {
    paste0("c(", paste(paste0(base_name, "_", group_keys), collapse = ", "), ")*")
  }

  mediator_lines <- vapply(seq_len(mediator_count), function(i) {
    lhs <- paste0("mediator_", i)
    if (is.null(covariate_rhs)) {
      paste(lhs, "~", multigroup_term(paste0("a", i)), "x")
    } else {
      paste(lhs, "~", multigroup_term(paste0("a", i)), "x +", covariate_rhs)
    }
  }, FUN.VALUE = character(1))

  outcome_terms <- vapply(seq_len(mediator_count), function(i) {
    paste(multigroup_term(paste0("b", i)), paste0("mediator_", i))
  }, FUN.VALUE = character(1))
  outcome_rhs <- paste(c(outcome_terms, paste(multigroup_term("c_prime"), "x"), covariate_rhs), collapse = " + ")
  outcome_line <- paste("Cog ~", outcome_rhs)

  mediator_cov_lines <- character(0)
  if (mediator_count >= 2) {
    mediator_pairs <- utils::combn(seq_len(mediator_count), 2, simplify = FALSE)
    mediator_cov_lines <- vapply(mediator_pairs, function(idx) {
      paste0("mediator_", idx[[1]], " ~~ mediator_", idx[[2]])
    }, FUN.VALUE = character(1))
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

  defined_lines <- unlist(lapply(seq_along(group_keys), function(g) {
    group_key <- group_keys[[g]]
    indirect_lines <- vapply(seq_len(mediator_count), function(i) {
      paste0("indirect_", i, "_", group_key, " := a", i, "_", group_key, "*b", i, "_", group_key)
    }, FUN.VALUE = character(1))
    total_indirect_terms <- paste0("indirect_", seq_len(mediator_count), "_", group_key)
    c(
      indirect_lines,
      paste0("total_indirect_", group_key, " := ", paste(total_indirect_terms, collapse = " + ")),
      paste0("c_total_", group_key, " := c_prime_", group_key, " + total_indirect_", group_key)
    )
  }))

  paste(
    c(
      measurement_line,
      mediator_lines,
      mediator_cov_lines,
      outcome_line,
      residual_lines,
      defined_lines
    ),
    collapse = "\n"
  )
}

extract_parallel_multigroup_indirects <- function(fit_object, group_levels, parallel_model, cognition_name, covariate_text, mediators) {
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)
  group_keys <- make.names(group_levels)
  out <- list()
  idx <- 1

  for (g in seq_along(group_levels)) {
    for (m in seq_along(mediators)) {
      lhs_name <- paste0("indirect_", m, "_", group_keys[[g]])
      row <- pe[pe$lhs == lhs_name & pe$op == ":=", , drop = FALSE]
      if (nrow(row) == 0) {
        next
      }
      row <- row[1, , drop = FALSE]
      out[[idx]] <- data.frame(
        parallel_model = parallel_model$name,
        x_var = parallel_model$x,
        mediator_var = mediators[[m]],
        cognition_model = cognition_name,
        group = group_levels[[g]],
        effect_type = "indirect",
        estimate = row$est[[1]],
        std_estimate = row$std.all[[1]],
        p_value = row$pvalue[[1]],
        conf_low = row$ci.lower[[1]],
        conf_high = row$ci.upper[[1]],
        covariates = covariate_text,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }

    total_lhs <- paste0("total_indirect_", group_keys[[g]])
    total_row <- pe[pe$lhs == total_lhs & pe$op == ":=", , drop = FALSE]
    if (nrow(total_row) > 0) {
      total_row <- total_row[1, , drop = FALSE]
      out[[idx]] <- data.frame(
        parallel_model = parallel_model$name,
        x_var = parallel_model$x,
        mediator_var = "TOTAL_INDIRECT",
        cognition_model = cognition_name,
        group = group_levels[[g]],
        effect_type = "total_indirect",
        estimate = total_row$est[[1]],
        std_estimate = total_row$std.all[[1]],
        p_value = total_row$pvalue[[1]],
        conf_low = total_row$ci.lower[[1]],
        conf_high = total_row$ci.upper[[1]],
        covariates = covariate_text,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  if (length(out) == 0) {
    return(data.frame())
  }
  do.call(rbind, out)
}

build_parallel_multigroup_constraints <- function(group_levels, mediator_index = NULL, effect_type = c("indirect", "total_indirect"), scope = c("omnibus", "pairwise")) {
  effect_type <- match.arg(effect_type)
  scope <- match.arg(scope)
  group_keys <- make.names(group_levels)
  path_expr <- if (effect_type == "indirect") {
    stopifnot(!is.null(mediator_index))
    paste0("a", mediator_index, "_", group_keys, "*b", mediator_index, "_", group_keys)
  } else {
    paste0("total_indirect_", group_keys)
  }

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

fit_parallel_multigroup_sem_result <- function(
  data,
  group_var,
  group_levels,
  cognition_name,
  parallel_model,
  sem_config,
  cognition_map,
  transformation_plan,
  covariates,
  covariate_text
) {
  adjustment <- resolve_adjustment(cognition_name, sem_config)
  residual_text <- format_residual_covariances(adjustment$residual_covariances)

  sem_input <- create_multigroup_parallel_sem_dataset(
    data = data,
    group_var = group_var,
    group_levels = group_levels,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = parallel_model$x,
    mediators = parallel_model$mediators,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_plan
  )

  lavaan_model <- build_multigroup_parallel_lavaan_model(
    indicator_names = sem_input$indicator_names,
    component_rows = sem_input$component_rows,
    covariates = covariates,
    adjustment = adjustment,
    group_levels = group_levels,
    mediator_count = length(parallel_model$mediators)
  )

  fit_object <- lavaan::sem(
    model = lavaan_model,
    data = sem_input$data,
    group = "sem_group",
    estimator = "MLR",
    missing = "fiml",
    fixed.x = FALSE,
    meanstructure = TRUE,
    std.lv = TRUE
  )

  fit_measures <- lavaan::fitMeasures(fit_object, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
  fit_eval <- evaluate_sem_fit(fit_measures)

  test_rows <- list()
  idx <- 1
  for (effect_type in sem_config$parallel_multigroup_tests$paths %||% c("indirect", "total_indirect")) {
    if (effect_type == "indirect") {
      for (m in seq_along(parallel_model$mediators)) {
        omnibus_constraints <- build_parallel_multigroup_constraints(group_levels, mediator_index = m, effect_type = "indirect", scope = "omnibus")
        test_rows[[idx]] <- run_multigroup_wald_test(
          fit_object = fit_object,
          constraints = omnibus_constraints,
          sem_model_name = parallel_model$name,
          cognition_name = cognition_name,
          path_name = paste0("indirect_", parallel_model$mediators[[m]]),
          scope = "omnibus"
        )
        idx <- idx + 1
        if (!isTRUE(sem_config$parallel_multigroup_tests$omnibus_only)) {
          pair_indices <- utils::combn(group_levels, 2, simplify = FALSE)
          pair_constraints <- build_parallel_multigroup_constraints(group_levels, mediator_index = m, effect_type = "indirect", scope = "pairwise")
          for (pair_idx in seq_along(pair_indices)) {
            pair_name <- paste(pair_indices[[pair_idx]], collapse = "_vs_")
            test_rows[[idx]] <- run_multigroup_wald_test(
              fit_object = fit_object,
              constraints = pair_constraints[[pair_idx]],
              sem_model_name = parallel_model$name,
              cognition_name = cognition_name,
              path_name = paste0("indirect_", parallel_model$mediators[[m]]),
              scope = "pairwise",
              group_pair = pair_name
            )
            idx <- idx + 1
          }
        }
      }
    } else {
      omnibus_constraints <- build_parallel_multigroup_constraints(group_levels, effect_type = "total_indirect", scope = "omnibus")
      test_rows[[idx]] <- run_multigroup_wald_test(
        fit_object = fit_object,
        constraints = omnibus_constraints,
        sem_model_name = parallel_model$name,
        cognition_name = cognition_name,
        path_name = "total_indirect",
        scope = "omnibus"
      )
      idx <- idx + 1
      if (!isTRUE(sem_config$parallel_multigroup_tests$omnibus_only)) {
        pair_indices <- utils::combn(group_levels, 2, simplify = FALSE)
        pair_constraints <- build_parallel_multigroup_constraints(group_levels, effect_type = "total_indirect", scope = "pairwise")
        for (pair_idx in seq_along(pair_indices)) {
          pair_name <- paste(pair_indices[[pair_idx]], collapse = "_vs_")
          test_rows[[idx]] <- run_multigroup_wald_test(
            fit_object = fit_object,
            constraints = pair_constraints[[pair_idx]],
            sem_model_name = parallel_model$name,
            cognition_name = cognition_name,
            path_name = "total_indirect",
            scope = "pairwise",
            group_pair = pair_name
          )
          idx <- idx + 1
        }
      }
    }
  }

  fit_row <- data.frame(
    parallel_model = parallel_model$name,
    x_var = parallel_model$x,
    mediators = paste(parallel_model$mediators, collapse = " + "),
    cognition_model = cognition_name,
    groups = paste(group_levels, collapse = ", "),
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    cfi = unname(fit_measures[["cfi"]]),
    tli = unname(fit_measures[["tli"]]),
    rmsea = unname(fit_measures[["rmsea"]]),
    srmr = unname(fit_measures[["srmr"]]),
    chisq = unname(fit_measures[["chisq"]]),
    df = unname(fit_measures[["df"]]),
    fit_indices_available = fit_eval$fit_indices_available,
    model_reasonable = fit_eval$model_reasonable,
    fit_note = fit_eval$fit_note,
    residual_covariances = residual_text,
    note = sem_config$parallel_multigroup_tests$note %||% "Parallel-mediator multigroup latent SEM was used to test indirect-effect differences across groups.",
    stringsAsFactors = FALSE
  )

  list(
    fit_row = fit_row,
    group_indirects = extract_parallel_multigroup_indirects(
      fit_object = fit_object,
      group_levels = group_levels,
      parallel_model = parallel_model,
      cognition_name = cognition_name,
      covariate_text = covariate_text,
      mediators = parallel_model$mediators
    ),
    test_rows = if (length(test_rows) > 0) do.call(rbind, test_rows) else data.frame()
  )
}

create_multigroup_sem_dataset <- function(data, group_var, group_levels, cognition_name, cognition_map, x_var, mediator_var, covariates, sem_config, transformation_table = NULL) {
  filtered_data <- data[
    !is.na(data[[group_var]]) & data[[group_var]] %in% group_levels,
    ,
    drop = FALSE
  ]

  sem_input <- create_sem_dataset(
    data = filtered_data,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = x_var,
    mediator_var = mediator_var,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_table
  )

  sem_input$data$sem_group <- factor(filtered_data[[group_var]], levels = group_levels)
  sem_input$group_levels <- group_levels
  sem_input$group_keys <- make.names(group_levels)
  sem_input
}

build_multigroup_lavaan_model <- function(indicator_names, component_rows, covariates, adjustment, group_levels) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  covariate_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  group_keys <- make.names(group_levels)

  multigroup_term <- function(base_name) {
    paste0("c(", paste(paste0(base_name, "_", group_keys), collapse = ", "), ")*")
  }

  mediator_line <- if (is.null(covariate_rhs)) {
    paste("mediator ~", multigroup_term("a"), "x")
  } else {
    paste("mediator ~", multigroup_term("a"), "x +", covariate_rhs)
  }

  outcome_line <- if (is.null(covariate_rhs)) {
    paste("Cog ~", multigroup_term("b"), "mediator +", multigroup_term("c_prime"), "x")
  } else {
    paste("Cog ~", multigroup_term("b"), "mediator +", multigroup_term("c_prime"), "x +", covariate_rhs)
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

  defined_lines <- unlist(lapply(group_keys, function(group_key) {
    c(
      paste0("indirect_", group_key, " := a_", group_key, "*b_", group_key),
      paste0("c_total_", group_key, " := c_prime_", group_key, " + (a_", group_key, "*b_", group_key, ")")
    )
  }))

  paste(
    c(
      measurement_line,
      mediator_line,
      outcome_line,
      residual_lines,
      defined_lines
    ),
    collapse = "\n"
  )
}

extract_multigroup_path_table <- function(fit_object, group_levels, sem_model, cognition_name, covariate_text) {
  pe <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)
  group_keys <- make.names(group_levels)
  group_lookup <- setNames(group_levels, seq_along(group_levels))

  extract_free_path <- function(lhs, rhs, path_name) {
    rows <- pe[pe$lhs == lhs & pe$op == "~" & pe$rhs == rhs, , drop = FALSE]
    if (nrow(rows) == 0) {
      return(data.frame())
    }

    data.frame(
      sem_model = sem_model$name,
      x_var = sem_model$x,
      mediator_var = sem_model$mediator,
      cognition_model = cognition_name,
      group = unname(group_lookup[as.character(rows$group)]),
      path = path_name,
      estimate = rows$est,
      std_estimate = rows$std.all,
      p_value = rows$pvalue,
      conf_low = rows$ci.lower,
      conf_high = rows$ci.upper,
      covariates = covariate_text,
      stringsAsFactors = FALSE
    )
  }

  extract_defined_by_group <- function(prefix, path_name) {
    rows <- do.call(
      rbind,
      lapply(seq_along(group_keys), function(i) {
        lhs_name <- paste0(prefix, "_", group_keys[[i]])
        row <- pe[pe$lhs == lhs_name & pe$op == ":=", , drop = FALSE]
        if (nrow(row) == 0) {
          return(NULL)
        }
        row <- row[1, , drop = FALSE]
        data.frame(
          sem_model = sem_model$name,
          x_var = sem_model$x,
          mediator_var = sem_model$mediator,
          cognition_model = cognition_name,
          group = group_levels[[i]],
          path = path_name,
          estimate = row$est[[1]],
          std_estimate = row$std.all[[1]],
          p_value = row$pvalue[[1]],
          conf_low = row$ci.lower[[1]],
          conf_high = row$ci.upper[[1]],
          covariates = covariate_text,
          stringsAsFactors = FALSE
        )
      })
    )
    if (is.null(rows)) {
      return(data.frame())
    }
    rows
  }

  do.call(
    rbind,
    list(
      extract_free_path("mediator", "x", "a"),
      extract_free_path("Cog", "mediator", "b"),
      extract_free_path("Cog", "x", "c_prime"),
      extract_defined_by_group("indirect", "indirect"),
      extract_defined_by_group("c_total", "c_total")
    )
  )
}

extract_wald_metrics <- function(test_object) {
  if (is.null(test_object)) {
    return(list(stat = NA_real_, df = NA_real_, p_value = NA_real_))
  }

  if (is.data.frame(test_object) || is.matrix(test_object)) {
    test_df <- as.data.frame(test_object)
    if (nrow(test_df) == 0) {
      return(list(stat = NA_real_, df = NA_real_, p_value = NA_real_))
    }
    stat_col <- grep("stat|chisq|wald", names(test_df), ignore.case = TRUE, value = TRUE)[1]
    df_col <- grep("^df$|df\\b", names(test_df), ignore.case = TRUE, value = TRUE)[1]
    p_col <- grep("p.value|pvalue|pr\\(", names(test_df), ignore.case = TRUE, value = TRUE)[1]
    return(list(
      stat = if (length(stat_col) > 0 && nzchar(stat_col)) test_df[[stat_col]][[1]] else NA_real_,
      df = if (length(df_col) > 0 && nzchar(df_col)) test_df[[df_col]][[1]] else NA_real_,
      p_value = if (length(p_col) > 0 && nzchar(p_col)) test_df[[p_col]][[1]] else NA_real_
    ))
  }

  if (is.list(test_object)) {
    return(list(
      stat = test_object$stat %||% NA_real_,
      df = test_object$df %||% NA_real_,
      p_value = test_object$p.value %||% NA_real_
    ))
  }

  list(stat = NA_real_, df = NA_real_, p_value = NA_real_)
}

build_multigroup_constraints <- function(group_levels, path_name, scope = c("omnibus", "pairwise")) {
  scope <- match.arg(scope)
  group_keys <- make.names(group_levels)
  path_expr <- switch(
    path_name,
    a = paste0("a_", group_keys),
    b = paste0("b_", group_keys),
    c_prime = paste0("c_prime_", group_keys),
    indirect = paste0("a_", group_keys, "*b_", group_keys),
    stop("Unsupported multigroup path name: ", path_name, call. = FALSE)
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

run_multigroup_wald_test <- function(fit_object, constraints, sem_model_name, cognition_name, path_name, scope, group_pair = NA_character_) {
  if (length(constraints) == 0) {
    return(data.frame())
  }

  test_result <- tryCatch(
    lavaan::lavTestWald(fit_object, constraints = paste(constraints, collapse = "\n")),
    error = function(e) NULL
  )
  metrics <- extract_wald_metrics(test_result)

  data.frame(
    sem_model = sem_model_name,
    cognition_model = cognition_name,
    path = path_name,
    scope = scope,
    group_pair = group_pair,
    statistic = metrics$stat,
    df = metrics$df,
    p_value = metrics$p_value,
    significant = !is.na(metrics$p_value) & metrics$p_value < 0.05,
    constraints = paste(constraints, collapse = " ; "),
    stringsAsFactors = FALSE
  )
}

fit_multigroup_sem_result <- function(
  data,
  group_var,
  group_levels,
  cognition_name,
  sem_model,
  sem_config,
  cognition_map,
  transformation_plan,
  covariates,
  covariate_text
) {
  adjustment <- resolve_adjustment(cognition_name, sem_config)
  residual_text <- format_residual_covariances(adjustment$residual_covariances)

  sem_input <- create_multigroup_sem_dataset(
    data = data,
    group_var = group_var,
    group_levels = group_levels,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = sem_model$x,
    mediator_var = sem_model$mediator,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_plan
  )

  lavaan_model <- build_multigroup_lavaan_model(
    indicator_names = sem_input$indicator_names,
    component_rows = sem_input$component_rows,
    covariates = covariates,
    adjustment = adjustment,
    group_levels = group_levels
  )

  fit_object <- lavaan::sem(
    model = lavaan_model,
    data = sem_input$data,
    group = "sem_group",
    estimator = "MLR",
    missing = "fiml",
    fixed.x = FALSE,
    meanstructure = TRUE,
    std.lv = TRUE
  )

  fit_measures <- lavaan::fitMeasures(fit_object, c("cfi", "tli", "rmsea", "srmr", "chisq", "df"))
  fit_eval <- evaluate_sem_fit(fit_measures)
  path_table <- extract_multigroup_path_table(
    fit_object = fit_object,
    group_levels = group_levels,
    sem_model = sem_model,
    cognition_name = cognition_name,
    covariate_text = covariate_text
  )

  test_rows <- list()
  row_idx <- 1
  for (path_name in sem_config$multigroup_tests$paths %||% c("a", "b", "c_prime", "indirect")) {
    omnibus_constraints <- build_multigroup_constraints(group_levels, path_name, scope = "omnibus")
    test_rows[[row_idx]] <- run_multigroup_wald_test(
      fit_object = fit_object,
      constraints = omnibus_constraints,
      sem_model_name = sem_model$name,
      cognition_name = cognition_name,
      path_name = path_name,
      scope = "omnibus"
    )
    row_idx <- row_idx + 1

    if (!isTRUE(sem_config$multigroup_tests$omnibus_only)) {
      pair_indices <- utils::combn(group_levels, 2, simplify = FALSE)
      pair_constraints <- build_multigroup_constraints(group_levels, path_name, scope = "pairwise")
      for (pair_idx in seq_along(pair_indices)) {
        pair_name <- paste(pair_indices[[pair_idx]], collapse = "_vs_")
        test_rows[[row_idx]] <- run_multigroup_wald_test(
          fit_object = fit_object,
          constraints = pair_constraints[[pair_idx]],
          sem_model_name = sem_model$name,
          cognition_name = cognition_name,
          path_name = path_name,
          scope = "pairwise",
          group_pair = pair_name
        )
        row_idx <- row_idx + 1
      }
    }
  }

  fit_row <- data.frame(
    sem_model = sem_model$name,
    x_var = sem_model$x,
    mediator_var = sem_model$mediator,
    cognition_model = cognition_name,
    groups = paste(group_levels, collapse = ", "),
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    cfi = unname(fit_measures[["cfi"]]),
    tli = unname(fit_measures[["tli"]]),
    rmsea = unname(fit_measures[["rmsea"]]),
    srmr = unname(fit_measures[["srmr"]]),
    chisq = unname(fit_measures[["chisq"]]),
    df = unname(fit_measures[["df"]]),
    fit_indices_available = fit_eval$fit_indices_available,
    model_reasonable = fit_eval$model_reasonable,
    fit_note = fit_eval$fit_note,
    residual_covariances = residual_text,
    note = sem_config$multigroup_tests$note %||% "Configural multigroup latent SEM was used to test structural-path differences.",
    stringsAsFactors = FALSE
  )

  list(
    fit_object = fit_object,
    fit_row = fit_row,
    path_table = path_table,
    test_rows = if (length(test_rows) > 0) do.call(rbind, test_rows) else data.frame()
  )
}

extract_path_metric <- function(path_table, path_name, metric = "estimate") {
  row <- path_table[path_table$path == path_name, , drop = FALSE]
  if (nrow(row) == 0 || !metric %in% names(row)) {
    return(NA_real_)
  }
  row[[metric]][[1]]
}

format_sem_effect_with_stars <- function(estimate, p_value, digits = report_digits) {
  if (is.na(estimate)) {
    return(NA_character_)
  }

  star_suffix <- ""
  if (!is.na(p_value) && p_value < 0.001) {
    star_suffix <- "**"
  } else if (!is.na(p_value) && p_value < 0.05) {
    star_suffix <- "*"
  }

  paste0(format_numeric_human(estimate, digits = digits), star_suffix)
}

format_sem_effect_plain <- function(estimate, digits = report_digits) {
  if (is.na(estimate)) {
    return(NA_character_)
  }
  format_numeric_human(estimate, digits = digits)
}

fit_sem_model_result <- function(
  data,
  group_name,
  cognition_name,
  sem_model,
  sem_config,
  cognition_map,
  transformation_plan,
  covariates,
  covariate_text
) {
  adjustment <- resolve_adjustment(cognition_name, sem_config)
  residual_text <- format_residual_covariances(adjustment$residual_covariances)

  sem_input <- create_sem_dataset(
    data = data,
    cognition_name = cognition_name,
    cognition_map = cognition_map,
    x_var = sem_model$x,
    mediator_var = sem_model$mediator,
    covariates = covariates,
    sem_config = sem_config,
    transformation_table = transformation_plan
  )

  lavaan_model <- build_lavaan_model(
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
  fit_eval <- evaluate_sem_fit(fit_measures)
  path_table <- extract_sem_path_table(fit_object)
  path_table$sem_model <- sem_model$name
  path_table$x_var <- sem_model$x
  path_table$mediator_var <- sem_model$mediator
  path_table$y_var <- cognition_name
  path_table$cognition_model <- cognition_name
  path_table$group <- group_name
  path_table$covariates <- covariate_text
  path_table$n <- lavaan::nobs(fit_object)

  summary_row <- data.frame(
    sem_model = sem_model$name,
    x_var = sem_model$x,
    mediator_var = sem_model$mediator,
    y_var = cognition_name,
    cognition_model = cognition_name,
    group = group_name,
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    indirect = path_table$estimate[path_table$path == "indirect"],
    indirect_std = path_table$std_estimate[path_table$path == "indirect"],
    indirect_p = path_table$p_value[path_table$path == "indirect"],
    total = path_table$estimate[path_table$path == "c_total"],
    total_std = path_table$std_estimate[path_table$path == "c_total"],
    total_p = path_table$p_value[path_table$path == "c_total"],
    direct = path_table$estimate[path_table$path == "c_prime"],
    direct_std = path_table$std_estimate[path_table$path == "c_prime"],
    direct_p = path_table$p_value[path_table$path == "c_prime"],
    proportion_mediated = path_table$estimate[path_table$path == "proportion_mediated"],
    proportion_mediated_pct = path_table$estimate[path_table$path == "proportion_mediated"] * 100,
    opposite_direction = (path_table$estimate[path_table$path == "c_prime"] * path_table$estimate[path_table$path == "indirect"]) < 0,
    mediation_type = ifelse(
      (path_table$estimate[path_table$path == "c_prime"] * path_table$estimate[path_table$path == "indirect"]) < 0,
      "inconsistent_mediation",
      "consistent_mediation"
    ),
    stringsAsFactors = FALSE
  )

  fit_row <- data.frame(
    sem_model = sem_model$name,
    x_var = sem_model$x,
    mediator_var = sem_model$mediator,
    y_var = cognition_name,
    cognition_model = cognition_name,
    group = group_name,
    covariates = covariate_text,
    n = lavaan::nobs(fit_object),
    engine = "lavaan_latent_sem",
    cfi = unname(fit_measures[["cfi"]]),
    tli = unname(fit_measures[["tli"]]),
    rmsea = unname(fit_measures[["rmsea"]]),
    srmr = unname(fit_measures[["srmr"]]),
    chisq = unname(fit_measures[["chisq"]]),
    df = unname(fit_measures[["df"]]),
    fit_indices_available = fit_eval$fit_indices_available,
    model_reasonable = fit_eval$model_reasonable,
    fit_note = fit_eval$fit_note,
    indicator_standardized = isTRUE(sem_config$indicator_standardization),
    residual_covariances = residual_text,
    optimization_note = adjustment$note,
    stringsAsFactors = FALSE
  )

  loadings <- lavaan::parameterEstimates(fit_object, standardized = TRUE, ci = TRUE)
  loadings <- loadings[loadings$lhs == "Cog" & loadings$op == "=~", c("rhs", "est", "std.all", "pvalue", "ci.lower", "ci.upper"), drop = FALSE]
  if (nrow(loadings) > 0) {
    names(loadings) <- c("indicator", "estimate", "std_estimate", "p_value", "conf_low", "conf_high")
    indicator_lookup <- setNames(sem_input$component_rows$component, sem_input$indicator_names)
    source_lookup <- setNames(sem_input$component_rows$analysis_var, sem_input$indicator_names)
    direction_lookup <- setNames(sem_input$component_rows$direction, sem_input$indicator_names)
    loadings$indicator_component <- unname(indicator_lookup[loadings$indicator])
    loadings$indicator_source <- unname(source_lookup[loadings$indicator])
    loadings$indicator_direction <- unname(direction_lookup[loadings$indicator])
    loadings$indicator_standardized <- isTRUE(sem_config$indicator_standardization)
    loadings$sem_model <- sem_model$name
    loadings$x_var <- sem_model$x
    loadings$mediator_var <- sem_model$mediator
    loadings$y_var <- cognition_name
    loadings$cognition_model <- cognition_name
    loadings$group <- group_name
  }

  mi_table <- extract_modification_indices(
    fit_object = fit_object,
    sem_model_name = sem_model$name,
    cognition_name = cognition_name,
    top_n = sem_config$modification_index_top_n %||% 10
  )
  if (nrow(mi_table) > 0) {
    mi_table$indicator_standardized <- isTRUE(sem_config$indicator_standardization)
    mi_table$residual_covariances <- residual_text
    mi_table$group <- group_name
  }

  list(
    fit_object = fit_object,
    path_table = path_table,
    summary_row = summary_row,
    fit_row = fit_row,
    loadings = loadings,
    mi_table = mi_table
  )
}

candidate_cognition_map <- build_cognition_component_map(
  cognition_models = project_config$cognition_model_candidates,
  transformation_table = transformation_plan
)

selection_cfg <- project_config$sem$cognition_model_selection
selection_sem_model <- Filter(function(x) identical(x$name, selection_cfg$reference_sem_model), project_config$sem$models)[[1]]
selection_group <- selection_cfg$reference_group %||% "Overall"
selection_data <- if (selection_group == "Overall") {
  analysis_data
} else {
  analysis_data[analysis_data[[project_config$variables$group_label_var]] == selection_group, , drop = FALSE]
}
selection_covariates <- unique(c(project_config$variables$covariates, selection_sem_model$extra_covariates %||% character(0)))
selection_covariate_text <- if (length(selection_covariates) == 0) "None" else paste(selection_covariates, collapse = " + ")

candidate_fit_rows <- lapply(names(project_config$cognition_model_candidates), function(cognition_name) {
  fit_result <- fit_sem_model_result(
    data = selection_data,
    group_name = selection_group,
    cognition_name = cognition_name,
    sem_model = selection_sem_model,
    sem_config = project_config$sem,
    cognition_map = candidate_cognition_map,
    transformation_plan = transformation_plan,
    covariates = selection_covariates,
    covariate_text = selection_covariate_text
  )

  cfg <- project_config$cognition_model_candidates[[cognition_name]]
  fit_row <- fit_result$fit_row
  data.frame(
    cognition_model = cognition_name,
    cognition_label_zh = resolve_variable_label(cognition_name, runtime_settings, "zh", fallback = cognition_name),
    cognition_label_en = resolve_variable_label(cognition_name, runtime_settings, "en", fallback = cognition_name),
    components = paste(cfg$components, collapse = " + "),
    n_indicators = length(cfg$components),
    reference_sem_model = selection_sem_model$name,
    reference_group = selection_group,
    n = fit_row$n[[1]],
    cfi = fit_row$cfi[[1]],
    tli = fit_row$tli[[1]],
    rmsea = fit_row$rmsea[[1]],
    srmr = fit_row$srmr[[1]],
    chisq = fit_row$chisq[[1]],
    df = fit_row$df[[1]],
    model_reasonable = fit_row$model_reasonable[[1]],
    fit_note = fit_row$fit_note[[1]],
    rmsea_good = !is.na(fit_row$rmsea[[1]]) && fit_row$rmsea[[1]] <= 0.06,
    srmr_good = !is.na(fit_row$srmr[[1]]) && fit_row$srmr[[1]] <= 0.08,
    selected_for_main_sem = cognition_name == (selection_cfg$preferred_model %||% "Cog_MMSE_MOCA"),
    stringsAsFactors = FALSE
  )
})
candidate_fit_table <- do.call(rbind, candidate_fit_rows)

path_rows <- list()
summary_rows <- list()
fit_rows <- list()
loading_rows <- list()
mi_rows <- list()
sensitivity_rows <- list()
difference_rows <- list()
parallel_summary_rows <- list()
parallel_fit_rows <- list()
parallel_multigroup_fit_rows <- list()
parallel_multigroup_group_rows <- list()
parallel_multigroup_test_rows <- list()
multigroup_fit_rows <- list()
multigroup_path_rows <- list()
multigroup_test_rows <- list()
plot_files <- c()
combined_plot_specs <- list()
idx <- 1
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
    for (sem_model in project_config$sem$models) {
      model_covariates <- unique(c(project_config$variables$covariates, sem_model$extra_covariates %||% character(0)))
      covariate_text <- if (length(model_covariates) == 0) "None" else paste(model_covariates, collapse = " + ")
      primary_result <- fit_sem_model_result(
        data = group_data,
        group_name = group_name,
        cognition_name = cognition_name,
        sem_model = sem_model,
        covariates = model_covariates,
        sem_config = project_config$sem,
        cognition_map = cognition_map,
        transformation_plan = transformation_plan,
        covariate_text = covariate_text
      )

      path_rows[[idx]] <- primary_result$path_table
      summary_rows[[idx]] <- primary_result$summary_row
      fit_rows[[idx]] <- primary_result$fit_row
      if (nrow(primary_result$loadings) > 0) {
        loading_rows[[idx]] <- primary_result$loadings
      }
      if (nrow(primary_result$mi_table) > 0) {
        mi_rows[[idx]] <- primary_result$mi_table
      }

      if (isTRUE(project_config$sem$indicator_standardization)) {
        sensitivity_sem_config <- project_config$sem
        sensitivity_sem_config$indicator_standardization <- FALSE
        sensitivity_result <- fit_sem_model_result(
          data = group_data,
          group_name = group_name,
          cognition_name = cognition_name,
          sem_model = sem_model,
          covariates = model_covariates,
          sem_config = sensitivity_sem_config,
          cognition_map = cognition_map,
          transformation_plan = transformation_plan,
          covariate_text = covariate_text
        )

        sensitivity_rows[[idx]] <- data.frame(
          group = group_name,
          sem_model = sem_model$name,
          x_var = sem_model$x,
          mediator_var = sem_model$mediator,
          cognition_model = cognition_name,
          covariates = covariate_text,
          n = primary_result$fit_row$n[[1]],
          primary_indicator_standardized = primary_result$fit_row$indicator_standardized[[1]],
          primary_cfi = primary_result$fit_row$cfi[[1]],
          primary_tli = primary_result$fit_row$tli[[1]],
          primary_rmsea = primary_result$fit_row$rmsea[[1]],
          primary_srmr = primary_result$fit_row$srmr[[1]],
          primary_model_reasonable = primary_result$fit_row$model_reasonable[[1]],
          primary_direct = primary_result$summary_row$direct[[1]],
          primary_direct_p = primary_result$summary_row$direct_p[[1]],
          primary_indirect = primary_result$summary_row$indirect[[1]],
          primary_indirect_p = primary_result$summary_row$indirect_p[[1]],
          sensitivity_indicator_standardized = sensitivity_result$fit_row$indicator_standardized[[1]],
          sensitivity_cfi = sensitivity_result$fit_row$cfi[[1]],
          sensitivity_tli = sensitivity_result$fit_row$tli[[1]],
          sensitivity_rmsea = sensitivity_result$fit_row$rmsea[[1]],
          sensitivity_srmr = sensitivity_result$fit_row$srmr[[1]],
          sensitivity_model_reasonable = sensitivity_result$fit_row$model_reasonable[[1]],
          sensitivity_direct = sensitivity_result$summary_row$direct[[1]],
          sensitivity_direct_p = sensitivity_result$summary_row$direct_p[[1]],
          sensitivity_indirect = sensitivity_result$summary_row$indirect[[1]],
          sensitivity_indirect_p = sensitivity_result$summary_row$indirect_p[[1]],
          direct_delta = primary_result$summary_row$direct[[1]] - sensitivity_result$summary_row$direct[[1]],
          indirect_delta = primary_result$summary_row$indirect[[1]] - sensitivity_result$summary_row$indirect[[1]],
          stringsAsFactors = FALSE
        )
      }

      x_label_en <- resolve_variable_label(sem_model$x, runtime_settings, "en", fallback = sem_model$x)
      mediator_label_en <- resolve_variable_label(sem_model$mediator, runtime_settings, "en", fallback = sem_model$mediator)
      y_label_en <- resolve_variable_label(cognition_name, runtime_settings, "en", fallback = cognition_name)
      model_label_en <- resolve_variable_label(sem_model$name, runtime_settings, "en", fallback = sem_model$name)
      group_label_en <- ifelse(group_name == "Overall", "Overall", group_name)

      plot_path <- file.path(result_figures_dir, paste0("sem_", group_name, "_", sem_model$name, "_", cognition_name, ".png"))
      plot_mediation_diagram(
        x_label = x_label_en,
        mediator_label = mediator_label_en,
        y_label = y_label_en,
        path_table = primary_result$path_table,
        path = plot_path,
        title = paste(group_label_en, model_label_en, sep = " | ")
      )
      plot_files <- c(plot_files, plot_path)

      panel_key <- cognition_name
      spec_key <- paste(sem_model$name, group_name, sep = "||")
      combined_plot_specs[[panel_key]][[spec_key]] <- list(
        x_label = x_label_en,
        mediator_label = mediator_label_en,
        y_label = y_label_en,
        path_table = primary_result$path_table,
        title = group_label_en,
        model_label = model_label_en
      )

      idx <- idx + 1
    }
  }
}

for (cognition_name in names(combined_plot_specs)) {
  combined_path <- file.path(result_figures_dir, paste0("sem_combined_", cognition_name, ".png"))
  cognition_label_en <- resolve_variable_label(cognition_name, runtime_settings, "en", fallback = cognition_name)
  ordered_specs <- list()
  model_names <- vapply(project_config$sem$models, function(x) x$name, FUN.VALUE = character(1))
  row_letters <- LETTERS[seq_along(model_names)]

  for (model_idx in seq_along(model_names)) {
    model_name <- model_names[[model_idx]]
    model_label_en <- resolve_variable_label(model_name, runtime_settings, "en", fallback = model_name)
    for (group_name in sem_group_order) {
      spec_key <- paste(model_name, group_name, sep = "||")
      spec <- combined_plot_specs[[cognition_name]][[spec_key]]
      if (is.null(spec)) {
        next
      }
      spec$row_label <- if (identical(group_name, sem_group_order[[1]])) {
        paste0(row_letters[[model_idx]], ". ", model_label_en)
      } else {
        ""
      }
      ordered_specs[[length(ordered_specs) + 1]] <- spec
    }
  }

  plot_mediation_panel(
    diagram_specs = ordered_specs,
    path = combined_path,
    title = paste("Combined SEM diagrams -", cognition_label_en),
    ncol = 4
  )
  plot_files <- c(plot_files, combined_path)
}

p_idx <- 1
for (group_name in names(sem_group_data)) {
  group_data <- sem_group_data[[group_name]]
  if (nrow(group_data) == 0) {
    next
  }

  for (cognition_name in names(project_config$cognition_models)) {
    for (parallel_model in project_config$sem$parallel_models %||% list()) {
      model_covariates <- unique(c(project_config$variables$covariates, parallel_model$extra_covariates %||% character(0)))
      covariate_text <- if (length(model_covariates) == 0) "None" else paste(model_covariates, collapse = " + ")

      parallel_result <- fit_parallel_sem_result(
        data = group_data,
        group_name = group_name,
        cognition_name = cognition_name,
        parallel_model = parallel_model,
        sem_config = project_config$sem,
        cognition_map = cognition_map,
        transformation_plan = transformation_plan,
        covariates = model_covariates,
        covariate_text = covariate_text
      )

      parallel_summary_rows[[p_idx]] <- parallel_result$mediator_rows
      parallel_fit_rows[[p_idx]] <- parallel_result$fit_row
      p_idx <- p_idx + 1
    }
  }
}

diff_idx <- 1
for (group_name in names(sem_group_data)) {
  group_data <- sem_group_data[[group_name]]
  if (nrow(group_data) == 0) {
    next
  }

  for (cognition_name in names(project_config$cognition_models)) {
    for (test_cfg in project_config$sem$indirect_difference_tests %||% list()) {
      adjustment <- resolve_adjustment(cognition_name, project_config$sem)
      comparison_covariates <- project_config$variables$covariates
      covariate_text <- if (length(comparison_covariates) == 0) "None" else paste(comparison_covariates, collapse = " + ")
      parallel_input <- create_parallel_sem_dataset(
        data = group_data,
        cognition_name = cognition_name,
        cognition_map = cognition_map,
        x_var = test_cfg$x,
        mediator_a = test_cfg$mediator_a,
        mediator_b = test_cfg$mediator_b,
        covariates = comparison_covariates,
        sem_config = project_config$sem,
        transformation_table = transformation_plan
      )

      parallel_model <- build_parallel_mediation_model(
        indicator_names = parallel_input$indicator_names,
        component_rows = parallel_input$component_rows,
        covariates = comparison_covariates,
        adjustment = adjustment
      )

      parallel_fit <- lavaan::sem(
        model = parallel_model,
        data = parallel_input$data,
        estimator = "MLR",
        missing = "fiml",
        fixed.x = FALSE,
        meanstructure = TRUE,
        std.lv = TRUE
      )

      ind_1 <- extract_defined_parameter(parallel_fit, "indirect_1")
      ind_2 <- extract_defined_parameter(parallel_fit, "indirect_2")
      ind_diff <- extract_defined_parameter(parallel_fit, "indirect_diff")

      difference_rows[[diff_idx]] <- data.frame(
        comparison_name = test_cfg$name,
        group = group_name,
        cognition_model = cognition_name,
        x_var = test_cfg$x,
        mediator_a = test_cfg$mediator_a,
        mediator_b = test_cfg$mediator_b,
        covariates = covariate_text,
        n = lavaan::nobs(parallel_fit),
        indirect_a = ind_1$estimate[[1]],
        indirect_a_std = ind_1$std_estimate[[1]],
        indirect_a_p = ind_1$p_value[[1]],
        indirect_b = ind_2$estimate[[1]],
        indirect_b_std = ind_2$std_estimate[[1]],
        indirect_b_p = ind_2$p_value[[1]],
        indirect_diff = ind_diff$estimate[[1]],
        indirect_diff_std = ind_diff$std_estimate[[1]],
        indirect_diff_p = ind_diff$p_value[[1]],
        indirect_diff_conf_low = ind_diff$conf_low[[1]],
        indirect_diff_conf_high = ind_diff$conf_high[[1]],
        difference_significant = !is.na(ind_diff$p_value[[1]]) & ind_diff$p_value[[1]] < 0.05,
        note = "Indirect-effect difference test was estimated in a parallel-mediator latent SEM with the same exposure, the same cognitive latent outcome, the same covariates, and both mediators entered simultaneously.",
        stringsAsFactors = FALSE
      )
      diff_idx <- diff_idx + 1
    }
  }
}

multigroup_cfg <- project_config$sem$multigroup_tests %||% list(enabled = FALSE)
if (isTRUE(multigroup_cfg$enabled)) {
  multigroup_groups <- multigroup_cfg$groups %||% c("CN", "MCI", "AD")
  mg_idx <- 1

  for (cognition_name in names(project_config$cognition_models)) {
    for (sem_model in project_config$sem$models) {
      model_covariates <- unique(c(project_config$variables$covariates, sem_model$extra_covariates %||% character(0)))
      covariate_text <- if (length(model_covariates) == 0) "None" else paste(model_covariates, collapse = " + ")

      multigroup_result <- fit_multigroup_sem_result(
        data = analysis_data,
        group_var = sem_group_var,
        group_levels = multigroup_groups,
        cognition_name = cognition_name,
        sem_model = sem_model,
        sem_config = project_config$sem,
        cognition_map = cognition_map,
        transformation_plan = transformation_plan,
        covariates = model_covariates,
        covariate_text = covariate_text
      )

      multigroup_fit_rows[[mg_idx]] <- multigroup_result$fit_row
      if (nrow(multigroup_result$path_table) > 0) {
        multigroup_path_rows[[mg_idx]] <- multigroup_result$path_table
      }
      if (nrow(multigroup_result$test_rows) > 0) {
        multigroup_test_rows[[mg_idx]] <- multigroup_result$test_rows
      }
      mg_idx <- mg_idx + 1
    }
  }
}

parallel_multigroup_cfg <- project_config$sem$parallel_multigroup_tests %||% list(enabled = FALSE)
if (isTRUE(parallel_multigroup_cfg$enabled)) {
  parallel_groups <- parallel_multigroup_cfg$groups %||% c("CN", "MCI", "AD")
  pmg_idx <- 1

  for (cognition_name in names(project_config$cognition_models)) {
    for (parallel_model in project_config$sem$parallel_models %||% list()) {
      model_covariates <- unique(c(project_config$variables$covariates, parallel_model$extra_covariates %||% character(0)))
      covariate_text <- if (length(model_covariates) == 0) "None" else paste(model_covariates, collapse = " + ")

      parallel_multigroup_result <- fit_parallel_multigroup_sem_result(
        data = analysis_data,
        group_var = sem_group_var,
        group_levels = parallel_groups,
        cognition_name = cognition_name,
        parallel_model = parallel_model,
        sem_config = project_config$sem,
        cognition_map = cognition_map,
        transformation_plan = transformation_plan,
        covariates = model_covariates,
        covariate_text = covariate_text
      )

      parallel_multigroup_fit_rows[[pmg_idx]] <- parallel_multigroup_result$fit_row
      if (nrow(parallel_multigroup_result$group_indirects) > 0) {
        parallel_multigroup_group_rows[[pmg_idx]] <- parallel_multigroup_result$group_indirects
      }
      if (nrow(parallel_multigroup_result$test_rows) > 0) {
        parallel_multigroup_test_rows[[pmg_idx]] <- parallel_multigroup_result$test_rows
      }
      pmg_idx <- pmg_idx + 1
    }
  }
}

sem_paths <- do.call(rbind, path_rows)
sem_summary <- do.call(rbind, summary_rows)
sem_fit <- do.call(rbind, fit_rows)
sem_loadings <- if (length(loading_rows) > 0) do.call(rbind, loading_rows) else data.frame()
sem_modification_indices <- if (length(mi_rows) > 0) do.call(rbind, mi_rows) else data.frame()
sem_standardization_sensitivity <- if (length(sensitivity_rows) > 0) do.call(rbind, sensitivity_rows) else data.frame()
sem_indirect_difference_tests <- if (length(difference_rows) > 0) do.call(rbind, difference_rows) else data.frame()
sem_parallel_summary <- if (length(parallel_summary_rows) > 0) do.call(rbind, parallel_summary_rows) else data.frame()
sem_parallel_fit <- if (length(parallel_fit_rows) > 0) do.call(rbind, parallel_fit_rows) else data.frame()
sem_parallel_multigroup_fit <- if (length(parallel_multigroup_fit_rows) > 0) do.call(rbind, parallel_multigroup_fit_rows) else data.frame()
sem_parallel_multigroup_groups <- if (length(parallel_multigroup_group_rows) > 0) do.call(rbind, parallel_multigroup_group_rows) else data.frame()
sem_parallel_multigroup_tests <- if (length(parallel_multigroup_test_rows) > 0) do.call(rbind, parallel_multigroup_test_rows) else data.frame()
sem_multigroup_fit <- if (length(multigroup_fit_rows) > 0) do.call(rbind, multigroup_fit_rows) else data.frame()
sem_multigroup_paths <- if (length(multigroup_path_rows) > 0) do.call(rbind, multigroup_path_rows) else data.frame()
sem_multigroup_tests <- if (length(multigroup_test_rows) > 0) do.call(rbind, multigroup_test_rows) else data.frame()

sem_stage_direction_summary <- if (nrow(sem_summary) > 0) {
  do.call(rbind, lapply(seq_len(nrow(sem_summary)), function(i) {
    row <- sem_summary[i, , drop = FALSE]
    row_paths <- sem_paths[
      sem_paths$sem_model == row$sem_model[[1]] &
        sem_paths$cognition_model == row$cognition_model[[1]] &
        sem_paths$group == row$group[[1]],
      ,
      drop = FALSE
    ]

    a_est <- extract_path_metric(row_paths, "a")
    a_p <- extract_path_metric(row_paths, "a", "p_value")
    b_est <- extract_path_metric(row_paths, "b")
    b_p <- extract_path_metric(row_paths, "b", "p_value")
    c_est <- extract_path_metric(row_paths, "c_prime")
    c_p <- extract_path_metric(row_paths, "c_prime", "p_value")

    effect_label <- if (isTRUE(row$opposite_direction[[1]])) {
      "opposite-direction path"
    } else {
      "same-direction path"
    }
    contribution_text <- if (is.na(row$proportion_mediated_pct[[1]])) {
      effect_label
    } else {
      paste0(effect_label, "; model-relative share = ", format_numeric_human(row$proportion_mediated_pct[[1]], digits = 2), "%")
    }
    data.frame(
      x = row$x_var[[1]],
      mediator = row$mediator_var[[1]],
      cognition = row$cognition_model[[1]],
      group = row$group[[1]],
      a = format_sem_effect_with_stars(a_est, a_p),
      b = format_sem_effect_with_stars(b_est, b_p),
      c = format_sem_effect_with_stars(c_est, c_p),
      indirect_effect = format_sem_effect_plain(row$indirect[[1]]),
      direct_effect = format_sem_effect_plain(row$direct[[1]]),
      total_effect = format_sem_effect_plain(row$total[[1]]),
      effect_pattern = contribution_text,
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame()
}

sem_mechanism_interpretation <- if (nrow(sem_summary) > 0) {
  do.call(rbind, lapply(seq_len(nrow(sem_summary)), function(i) {
    row <- sem_summary[i, , drop = FALSE]
    significant_indirect <- !is.na(row$indirect_p[[1]]) && row$indirect_p[[1]] < 0.05
    pattern_key <- if (!significant_indirect) {
      "no_significant_mediation"
    } else if (isTRUE(row$opposite_direction[[1]])) {
      "inconsistent_mediation"
    } else {
      "consistent_mediation"
    }

    pattern_label <- if (pattern_key == "consistent_mediation") {
      "same-direction indirect path"
    } else if (pattern_key == "inconsistent_mediation") {
      "opposite-direction indirect path"
    } else {
      "no statistically clear indirect path"
    }

    interpretation <- if (pattern_key == "consistent_mediation") {
      "The indirect and direct paths point in the same direction within this cross-sectional path-decomposition model."
    } else if (pattern_key == "inconsistent_mediation") {
      "The indirect and direct paths point in opposite directions; this is reported as a suppression-type pattern, not as evidence of protection."
    } else {
      "The indirect path did not meet the nominal significance threshold in this model."
    }
    data.frame(
      sem_model = row$sem_model[[1]],
      x_var = row$x_var[[1]],
      mediator_var = row$mediator_var[[1]],
      cognition_model = row$cognition_model[[1]],
      group = row$group[[1]],
      indirect = row$indirect[[1]],
      indirect_p = row$indirect_p[[1]],
      direct = row$direct[[1]],
      direct_p = row$direct_p[[1]],
      total = row$total[[1]],
      total_p = row$total_p[[1]],
      pattern = pattern_key,
      pattern_label = pattern_label,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame()
}

sem_moderated_mediation_summary <- {
  single_rows <- if (nrow(sem_multigroup_tests) > 0) {
    data.frame(
      analysis_type = "single_mediator",
      model_name = sem_multigroup_tests$sem_model,
      mediator_component = sem_multigroup_tests$path,
      cognition_model = sem_multigroup_tests$cognition_model,
      scope = sem_multigroup_tests$scope,
      group_pair = sem_multigroup_tests$group_pair,
      statistic = sem_multigroup_tests$statistic,
      df = sem_multigroup_tests$df,
      p_value = sem_multigroup_tests$p_value,
      significant = sem_multigroup_tests$significant,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  if (nrow(single_rows) > 0) {
    single_rows <- single_rows[single_rows$mediator_component %in% c("indirect"), , drop = FALSE]
  }

  parallel_rows <- if (nrow(sem_parallel_multigroup_tests) > 0) {
    data.frame(
      analysis_type = "parallel_mediator",
      model_name = sem_parallel_multigroup_tests$sem_model,
      mediator_component = sem_parallel_multigroup_tests$path,
      cognition_model = sem_parallel_multigroup_tests$cognition_model,
      scope = sem_parallel_multigroup_tests$scope,
      group_pair = sem_parallel_multigroup_tests$group_pair,
      statistic = sem_parallel_multigroup_tests$statistic,
      df = sem_parallel_multigroup_tests$df,
      p_value = sem_parallel_multigroup_tests$p_value,
      significant = sem_parallel_multigroup_tests$significant,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame()
  }

  out <- rbind(single_rows, parallel_rows)
  if (nrow(out) == 0) {
    data.frame()
  } else {
    out
  }
}

paths_path <- file.path(result_summary_dir, "sem_path_coefficients.csv")
summary_path <- file.path(result_summary_dir, "sem_mediation_summary.csv")
fit_path <- file.path(result_summary_dir, "sem_model_fit.csv")
loadings_path <- file.path(result_summary_dir, "sem_factor_loadings.csv")
mi_path <- file.path(result_summary_dir, "sem_modification_indices_top.csv")
sensitivity_path <- file.path(result_summary_dir, "sem_standardization_sensitivity.csv")
indirect_difference_path <- file.path(result_summary_dir, "sem_indirect_difference_tests.csv")
parallel_summary_path <- file.path(result_summary_dir, "sem_parallel_mediation_summary.csv")
parallel_fit_path <- file.path(result_summary_dir, "sem_parallel_model_fit.csv")
parallel_multigroup_fit_path <- file.path(result_summary_dir, "sem_parallel_multigroup_model_fit.csv")
parallel_multigroup_groups_path <- file.path(result_summary_dir, "sem_parallel_multigroup_group_estimates.csv")
parallel_multigroup_tests_path <- file.path(result_summary_dir, "sem_parallel_multigroup_tests.csv")
moderated_mediation_summary_path <- file.path(result_summary_dir, "sem_moderated_mediation_summary.csv")
stage_summary_path <- file.path(result_summary_dir, "sem_stage_direction_summary.csv")
multigroup_fit_path <- file.path(result_summary_dir, "sem_multigroup_model_fit.csv")
multigroup_paths_path <- file.path(result_summary_dir, "sem_multigroup_group_estimates.csv")
multigroup_tests_path <- file.path(result_summary_dir, "sem_multigroup_path_tests.csv")
mechanism_path <- file.path(result_summary_dir, "sem_mechanism_interpretation.csv")
candidate_fit_path <- file.path(result_summary_dir, "cognition_model_selection_fit.csv")
table_csv <- file.path(result_tables_dir, "table_sem_mediation_summary.csv")
table_html <- file.path(result_tables_dir, "table_sem_mediation_summary.html")
candidate_md_path <- file.path(result_report_dir, "Cognition_model_selection.md")

write_csv_utf8(sem_paths, paths_path, row.names = FALSE)
write_csv_utf8(sem_summary, summary_path, row.names = FALSE)
write_csv_utf8(sem_fit, fit_path, row.names = FALSE)
write_csv_utf8(sem_loadings, loadings_path, row.names = FALSE)
write_csv_utf8(sem_modification_indices, mi_path, row.names = FALSE)
write_csv_utf8(sem_standardization_sensitivity, sensitivity_path, row.names = FALSE)
write_csv_utf8(sem_indirect_difference_tests, indirect_difference_path, row.names = FALSE)
write_csv_utf8(sem_parallel_summary, parallel_summary_path, row.names = FALSE)
write_csv_utf8(sem_parallel_fit, parallel_fit_path, row.names = FALSE)
write_csv_utf8(sem_parallel_multigroup_fit, parallel_multigroup_fit_path, row.names = FALSE)
write_csv_utf8(sem_parallel_multigroup_groups, parallel_multigroup_groups_path, row.names = FALSE)
write_csv_utf8(sem_parallel_multigroup_tests, parallel_multigroup_tests_path, row.names = FALSE)
write_csv_utf8(sem_moderated_mediation_summary, moderated_mediation_summary_path, row.names = FALSE)
write_csv_utf8(sem_stage_direction_summary, stage_summary_path, row.names = FALSE)
write_csv_utf8(sem_multigroup_fit, multigroup_fit_path, row.names = FALSE)
write_csv_utf8(sem_multigroup_paths, multigroup_paths_path, row.names = FALSE)
write_csv_utf8(sem_multigroup_tests, multigroup_tests_path, row.names = FALSE)
write_csv_utf8(sem_mechanism_interpretation, mechanism_path, row.names = FALSE)
write_csv_utf8(candidate_fit_table, candidate_fit_path, row.names = FALSE)

export_three_line_table(
  data = sem_summary,
  csv_path = table_csv,
  html_path = table_html,
  title = "Latent-variable SEM Mediation Summary"
)


append_analysis_log(
  project_root = project_root,
  analysis_name = "07_sem_mediation",
  output_files = c(
    paths_path, summary_path, fit_path, loadings_path, mi_path, sensitivity_path,
    indirect_difference_path, parallel_summary_path, parallel_fit_path,
    parallel_multigroup_fit_path, parallel_multigroup_groups_path,
    parallel_multigroup_tests_path, moderated_mediation_summary_path,
    stage_summary_path, multigroup_fit_path, multigroup_paths_path,
    multigroup_tests_path, mechanism_path, candidate_fit_path, table_csv,
    table_html, plot_files
  ),
  note = "Completed latent-variable SEM and path-decomposition models and exported analysis tables/figures. No submission text is generated in the public release.",
  summary_dir = result_summary_dir
)
