source(file.path(getwd(), "00_setup.R"))

local_r_lib <- normalizePath(file.path(repo_root, ".r_libs"), winslash = "/", mustWork = FALSE)
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop("Package 'lavaan' is required for bootstrap SEM summaries.", call. = FALSE)
}

analysis_data <- read_project_data(analysis_data_path)
transformation_plan <- read_project_data(transformation_table_path)
cognition_map <- read_project_data(cognition_map_path)
preferred_cognition <- project_config$sem$cognition_model_selection$preferred_model %||% names(project_config$cognition_models)[1]
boot_n <- project_config$sem$bootstrap %||% 1000

create_sem_dataset_local <- function(data, cognition_name, cognition_map, x_var, mediator_var, covariates, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_analysis <- resolve_analysis_var(mediator_var, transformation_table)
  sem_data <- data.frame(
    x = as.numeric(data[[x_analysis]]),
    mediator = as.numeric(data[[mediator_analysis]])
  )
  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_names <- c(indicator_names, indicator_name)
    sem_data[[indicator_name]] <- as.numeric(data[[component_rows$analysis_var[[i]]]]) * as.numeric(component_rows$direction[[i]])
  }
  for (cov_name in covariates) {
    cov_analysis <- resolve_analysis_var(cov_name, transformation_table)
    if (cov_name %in% project_config$variables$categorical_covariates) {
      sem_data[[cov_name]] <- factor(data[[cov_name]])
    } else {
      sem_data[[cov_name]] <- as.numeric(data[[cov_analysis]])
    }
  }
  sem_data <- sem_data[stats::complete.cases(sem_data), , drop = FALSE]
  list(data = sem_data, indicator_names = indicator_names)
}

build_single_mediator_model <- function(indicator_names, covariates) {
  cov_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  mediator_line <- if (is.null(cov_rhs)) {
    "mediator ~ a*x"
  } else {
    paste("mediator ~ a*x +", cov_rhs)
  }
  outcome_line <- if (is.null(cov_rhs)) {
    "Cog ~ b*mediator + c_prime*x"
  } else {
    paste("Cog ~ b*mediator + c_prime*x +", cov_rhs)
  }
  paste(
    c(
      paste("Cog =~", paste(indicator_names, collapse = " + ")),
      mediator_line,
      outcome_line,
      "indirect := a*b",
      "c_total := c_prime + (a*b)"
    ),
    collapse = "\n"
  )
}

create_serial_sem_dataset_local <- function(data, cognition_name, cognition_map, x_var, mediator_1, mediator_2, covariates, transformation_table = NULL) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  sem_data <- data.frame(
    x = as.numeric(data[[resolve_analysis_var(x_var, transformation_table)]]),
    mediator_1 = as.numeric(data[[resolve_analysis_var(mediator_1, transformation_table)]]),
    mediator_2 = as.numeric(data[[resolve_analysis_var(mediator_2, transformation_table)]])
  )
  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_names <- c(indicator_names, indicator_name)
    sem_data[[indicator_name]] <- as.numeric(data[[component_rows$analysis_var[[i]]]]) * as.numeric(component_rows$direction[[i]])
  }
  for (cov_name in covariates) {
    cov_analysis <- resolve_analysis_var(cov_name, transformation_table)
    if (cov_name %in% project_config$variables$categorical_covariates) {
      sem_data[[cov_name]] <- factor(data[[cov_name]])
    } else {
      sem_data[[cov_name]] <- as.numeric(data[[cov_analysis]])
    }
  }
  sem_data <- sem_data[stats::complete.cases(sem_data), , drop = FALSE]
  list(data = sem_data, indicator_names = indicator_names)
}

build_serial_mediator_model <- function(indicator_names, covariates) {
  cov_rhs <- if (length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  m1_line <- if (is.null(cov_rhs)) {
    "mediator_1 ~ a1*x"
  } else {
    paste("mediator_1 ~ a1*x +", cov_rhs)
  }
  m2_line <- if (is.null(cov_rhs)) {
    "mediator_2 ~ d*mediator_1 + a2*x"
  } else {
    paste("mediator_2 ~ d*mediator_1 + a2*x +", cov_rhs)
  }
  y_line <- if (is.null(cov_rhs)) {
    "Cog ~ b*mediator_2 + b1*mediator_1 + c_prime*x"
  } else {
    paste("Cog ~ b*mediator_2 + b1*mediator_1 + c_prime*x +", cov_rhs)
  }
  paste(
    c(
      paste("Cog =~", paste(indicator_names, collapse = " + ")),
      m1_line,
      m2_line,
      y_line,
      "serial_indirect := a1*d*b",
      "m1_only_indirect := a1*b1",
      "m2_only_indirect := a2*b",
      "total_indirect := serial_indirect + m1_only_indirect + m2_only_indirect",
      "c_total := c_prime + total_indirect"
    ),
    collapse = "\n"
  )
}

extract_defined_effect <- function(pe, lhs_name, effect_name) {
  row <- pe[pe$lhs == lhs_name & pe$op == ":=", , drop = FALSE]
  if (nrow(row) == 0) {
    return(data.frame())
  }
  data.frame(
    effect = effect_name,
    estimate = row$est[[1]],
    conf_low = row$ci.lower[[1]],
    conf_high = row$ci.upper[[1]],
    p_value = row$pvalue[[1]],
    stringsAsFactors = FALSE
  )
}

run_single_bootstrap <- function(model_cfg) {
  sem_parts <- create_sem_dataset_local(
    data = analysis_data,
    cognition_name = preferred_cognition,
    cognition_map = cognition_map,
    x_var = model_cfg$x,
    mediator_var = model_cfg$mediator,
    covariates = project_config$variables$covariates,
    transformation_table = transformation_plan
  )
  model_syntax <- build_single_mediator_model(sem_parts$indicator_names, project_config$variables$covariates)
  fit <- lavaan::sem(model_syntax, data = sem_parts$data, se = "bootstrap", bootstrap = boot_n, missing = "fiml")
  pe <- lavaan::parameterEstimates(fit, ci = TRUE)
  out <- rbind(
    extract_defined_effect(pe, "indirect", "indirect"),
    extract_defined_effect(pe, "c_total", "total")
  )
  cprime_row <- pe[pe$lhs == "Cog" & pe$op == "~" & pe$rhs == "x", , drop = FALSE]
  if (nrow(cprime_row) > 0) {
    out <- rbind(
      out,
      data.frame(
        effect = "direct",
        estimate = cprime_row$est[[1]],
        conf_low = cprime_row$ci.lower[[1]],
        conf_high = cprime_row$ci.upper[[1]],
        p_value = cprime_row$pvalue[[1]],
        stringsAsFactors = FALSE
      )
    )
  }
  out$model_name <- model_cfg$name
  out$x <- model_cfg$x
  out$mediator <- model_cfg$mediator
  out$cognition_model <- preferred_cognition
  out$n <- nrow(sem_parts$data)
  out
}

run_serial_bootstrap <- function(model_cfg) {
  sem_parts <- create_serial_sem_dataset_local(
    data = analysis_data,
    cognition_name = preferred_cognition,
    cognition_map = cognition_map,
    x_var = model_cfg$x,
    mediator_1 = model_cfg$mediator_1,
    mediator_2 = model_cfg$mediator_2,
    covariates = project_config$variables$covariates,
    transformation_table = transformation_plan
  )
  model_syntax <- build_serial_mediator_model(sem_parts$indicator_names, project_config$variables$covariates)
  fit <- lavaan::sem(model_syntax, data = sem_parts$data, se = "bootstrap", bootstrap = boot_n, missing = "fiml")
  pe <- lavaan::parameterEstimates(fit, ci = TRUE)
  out <- rbind(
    extract_defined_effect(pe, "serial_indirect", "serial_indirect"),
    extract_defined_effect(pe, "total_indirect", "total_indirect"),
    extract_defined_effect(pe, "c_total", "total")
  )
  cprime_row <- pe[pe$lhs == "Cog" & pe$op == "~" & pe$rhs == "x", , drop = FALSE]
  if (nrow(cprime_row) > 0) {
    out <- rbind(
      out,
      data.frame(
        effect = "direct",
        estimate = cprime_row$est[[1]],
        conf_low = cprime_row$ci.lower[[1]],
        conf_high = cprime_row$ci.upper[[1]],
        p_value = cprime_row$pvalue[[1]],
        stringsAsFactors = FALSE
      )
    )
  }
  out$model_name <- model_cfg$name
  out$x <- model_cfg$x
  out$mediator <- paste(model_cfg$mediator_1, model_cfg$mediator_2, sep = " -> ")
  out$cognition_model <- preferred_cognition
  out$n <- nrow(sem_parts$data)
  out
}

single_rows <- do.call(rbind, lapply(project_config$sem$models, run_single_bootstrap))
serial_rows <- do.call(rbind, lapply(project_config$sem$serial_models, run_serial_bootstrap))
bootstrap_rows <- rbind(single_rows, serial_rows)
bootstrap_rows$analysis_type <- ifelse(grepl("serial", bootstrap_rows$effect), "serial", "single_or_total")

bootstrap_path <- file.path(result_summary_dir, "sem_bootstrap_key_indirects.csv")
write_csv_utf8(bootstrap_rows, bootstrap_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "07d_bootstrap_key_indirects",
  output_files = bootstrap_path,
  note = "Re-estimated key single-mediator and serial-mediator SEM paths using bootstrap confidence intervals for the preferred cognition model.",
  summary_dir = result_summary_dir
)
