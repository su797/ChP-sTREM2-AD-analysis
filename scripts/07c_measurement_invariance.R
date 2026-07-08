source(file.path(getwd(), "00_setup.R"))

local_r_lib <- normalizePath(file.path(repo_root, ".r_libs"), winslash = "/", mustWork = FALSE)
if (dir.exists(local_r_lib)) {
  .libPaths(c(local_r_lib, .libPaths()))
}

if (!requireNamespace("lavaan", quietly = TRUE)) {
  stop("Package 'lavaan' is required for measurement invariance analysis.", call. = FALSE)
}

analysis_data <- read_project_data(analysis_data_path)
cognition_map <- read_project_data(cognition_map_path)
transformation_plan <- read_project_data(transformation_table_path)

safe_scale_invariance <- function(x) {
  x <- as.numeric(x)
  x_mean <- mean(x, na.rm = TRUE)
  x_sd <- stats::sd(x, na.rm = TRUE)
  if (is.na(x_sd) || x_sd == 0) {
    return(x - x_mean)
  }
  (x - x_mean) / x_sd
}

prepare_invariance_dataset <- function(data, cognition_name, group_var, group_levels, sem_config, transformation_table) {
  component_rows <- cognition_map[cognition_map$composite == cognition_name, , drop = FALSE]
  data <- data[data[[group_var]] %in% group_levels, , drop = FALSE]
  data[[group_var]] <- factor(data[[group_var]], levels = group_levels)

  out <- data.frame(.group = data[[group_var]])
  indicator_names <- character(0)
  for (i in seq_len(nrow(component_rows))) {
    indicator_name <- paste0("indicator_", i)
    indicator_source <- component_rows$analysis_var[[i]]
    indicator_direction <- as.numeric(component_rows$direction[[i]])
    indicator_values <- indicator_direction * as.numeric(data[[indicator_source]])
    if (isTRUE(sem_config$indicator_standardization)) {
      indicator_values <- safe_scale_invariance(indicator_values)
    }
    out[[indicator_name]] <- indicator_values
    indicator_names <- c(indicator_names, indicator_name)
  }

  list(data = out, indicator_names = indicator_names, component_rows = component_rows)
}

resolve_adjustment_invariance <- function(cognition_name, sem_config) {
  adjustment <- sem_config$cognition_adjustments[[cognition_name]]
  if (is.null(adjustment)) {
    return(list(residual_covariances = list()))
  }
  if (is.null(adjustment$residual_covariances)) {
    adjustment$residual_covariances <- list()
  }
  adjustment
}

build_invariance_model <- function(indicator_names, component_rows, adjustment) {
  measurement_line <- paste("Cog =~", paste(indicator_names, collapse = " + "))
  indicator_lookup <- setNames(indicator_names, component_rows$component)
  residual_lines <- character(0)
  for (pair in adjustment$residual_covariances) {
    if (length(pair) != 2) next
    lhs <- indicator_lookup[[pair[[1]]]]
    rhs <- indicator_lookup[[pair[[2]]]]
    if (is.null(lhs) || is.null(rhs)) next
    residual_lines <- c(residual_lines, paste(lhs, "~~", rhs))
  }
  paste(c(measurement_line, residual_lines), collapse = "\n")
}

extract_fit_row_invariance <- function(model_name, fit_object) {
  if (is.null(fit_object) || !inherits(fit_object, "lavaan")) {
    return(data.frame(
      invariance_model = model_name,
      cfi = NA_real_,
      tli = NA_real_,
      rmsea = NA_real_,
      srmr = NA_real_,
      chisq = NA_real_,
      df = NA_real_,
      aic = NA_real_,
      bic = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  converged <- tryCatch(lavaan::lavInspect(fit_object, "converged"), error = function(e) FALSE)
  if (!isTRUE(converged)) {
    return(data.frame(
      invariance_model = model_name,
      cfi = NA_real_,
      tli = NA_real_,
      rmsea = NA_real_,
      srmr = NA_real_,
      chisq = NA_real_,
      df = NA_real_,
      aic = NA_real_,
      bic = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  fm <- lavaan::fitMeasures(fit_object, c("cfi", "tli", "rmsea", "srmr", "chisq", "df", "aic", "bic"))
  data.frame(
    invariance_model = model_name,
    cfi = fm[["cfi"]],
    tli = fm[["tli"]],
    rmsea = fm[["rmsea"]],
    srmr = fm[["srmr"]],
    chisq = fm[["chisq"]],
    df = fm[["df"]],
    aic = fm[["aic"]],
    bic = fm[["bic"]],
    stringsAsFactors = FALSE
  )
}

evaluate_invariance_step <- function(delta_cfi, delta_rmsea, delta_srmr) {
  if (is.na(delta_cfi) || is.na(delta_rmsea) || is.na(delta_srmr)) {
    return("Unable to judge.")
  }
  if (abs(delta_cfi) <= 0.010 && abs(delta_rmsea) <= 0.015 && abs(delta_srmr) <= 0.030) {
    return("Supported.")
  }
  "Not supported."
}

group_var <- project_config$variables$group_label_var
group_levels <- c("CN", "MCI", "AD")
cognition_name <- "Cog_MMSE_MOCA"

invariance_input <- prepare_invariance_dataset(
  data = analysis_data,
  cognition_name = cognition_name,
  group_var = group_var,
  group_levels = group_levels,
  sem_config = project_config$sem,
  transformation_table = transformation_plan
)

adjustment <- resolve_adjustment_invariance(cognition_name, project_config$sem)
lavaan_model <- build_invariance_model(
  indicator_names = invariance_input$indicator_names,
  component_rows = invariance_input$component_rows,
  adjustment = adjustment
)

configural_fit <- lavaan::cfa(
  model = lavaan_model,
  data = invariance_input$data,
  group = ".group",
  estimator = "MLR",
  missing = "fiml",
  std.lv = TRUE
)

metric_fit <- lavaan::cfa(
  model = lavaan_model,
  data = invariance_input$data,
  group = ".group",
  group.equal = c("loadings"),
  estimator = "MLR",
  missing = "fiml",
  std.lv = TRUE
)

scalar_fit <- lavaan::cfa(
  model = lavaan_model,
  data = invariance_input$data,
  group = ".group",
  group.equal = c("loadings", "intercepts"),
  estimator = "MLR",
  missing = "fiml",
  std.lv = TRUE
)

strict_fit <- tryCatch(
  lavaan::cfa(
    model = lavaan_model,
    data = invariance_input$data,
    group = ".group",
    group.equal = c("loadings", "intercepts", "residuals"),
    estimator = "MLR",
    missing = "fiml",
    std.lv = TRUE
  ),
  error = function(e) NULL
)

fit_table <- do.call(
  rbind,
  list(
    extract_fit_row_invariance("configural", configural_fit),
    extract_fit_row_invariance("metric", metric_fit),
    extract_fit_row_invariance("scalar", scalar_fit),
    extract_fit_row_invariance("strict", strict_fit)
  )
)

comparison_table <- fit_table
comparison_table$delta_cfi <- c(NA_real_, diff(comparison_table$cfi))
comparison_table$delta_rmsea <- c(NA_real_, diff(comparison_table$rmsea))
comparison_table$delta_srmr <- c(NA_real_, diff(comparison_table$srmr))
comparison_table$decision <- c(
  "Reference model.",
  evaluate_invariance_step(comparison_table$delta_cfi[[2]], comparison_table$delta_rmsea[[2]], comparison_table$delta_srmr[[2]]),
  evaluate_invariance_step(comparison_table$delta_cfi[[3]], comparison_table$delta_rmsea[[3]], comparison_table$delta_srmr[[3]]),
  if (all(is.na(comparison_table[4, c("cfi", "rmsea", "srmr")]))) {
    "Strict model did not converge."
  } else {
    evaluate_invariance_step(comparison_table$delta_cfi[[4]], comparison_table$delta_rmsea[[4]], comparison_table$delta_srmr[[4]])
  }
)

report_lines <- c(
  "# Measurement Invariance of MMSE + MoCA",
  "",
  "This follow-up analysis tested whether the latent cognition model based on MMSE + MoCA is comparable across CN, MCI, and AD.",
  "",
  "## Key rule",
  "",
  "A more constrained model was considered practically acceptable when delta CFI <= 0.010, delta RMSEA <= 0.015, and delta SRMR <= 0.030.",
  "",
  "## Interpretation",
  ""
)

for (i in seq_len(nrow(comparison_table))) {
  row <- comparison_table[i, , drop = FALSE]
  report_lines <- c(
    report_lines,
    paste0(
      "- `", row$invariance_model[[1]], "`: CFI=", format_numeric_human(row$cfi[[1]], digits = 3),
      ", TLI=", format_numeric_human(row$tli[[1]], digits = 3),
      ", RMSEA=", format_numeric_human(row$rmsea[[1]], digits = 3),
      ", SRMR=", format_numeric_human(row$srmr[[1]], digits = 3),
      if (!is.na(row$delta_cfi[[1]])) paste0("; delta CFI=", format_numeric_human(row$delta_cfi[[1]], digits = 3)) else "",
      if (!is.na(row$delta_rmsea[[1]])) paste0(", delta RMSEA=", format_numeric_human(row$delta_rmsea[[1]], digits = 3)) else "",
      if (!is.na(row$delta_srmr[[1]])) paste0(", delta SRMR=", format_numeric_human(row$delta_srmr[[1]], digits = 3)) else "",
      ". Decision: ", row$decision[[1]]
    )
  )
}

report_lines <- c(
  report_lines,
  "",
  "## Why it matters",
  "",
  "If metric/scalar invariance is reasonably supported, then CN/MCI/AD group differences in structural SEM paths can be interpreted with greater confidence at the latent-cognition level."
)

fit_path <- file.path(result_summary_dir, "sem_measurement_invariance_fit.csv")
comparison_path <- file.path(result_summary_dir, "sem_measurement_invariance_comparison.csv")
report_path <- file.path(result_report_dir, "SEM_measurement_invariance.md")

write_csv_utf8(fit_table, fit_path, row.names = FALSE)
write_csv_utf8(comparison_table, comparison_path, row.names = FALSE)
writeLines(enc2utf8(report_lines), con = report_path, useBytes = TRUE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "07c_measurement_invariance",
  note = "Completed MMSE+MoCA measurement invariance analysis across CN, MCI, and AD.",
  summary_dir = result_summary_dir
)