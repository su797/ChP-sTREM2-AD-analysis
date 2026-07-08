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

build_alias_df <- function(data, alias_map, factor_aliases = character(0), transformation_table = NULL) {
  out <- list()
  for (alias_name in names(alias_map)) {
    source_var <- alias_map[[alias_name]]
    if (alias_name %in% factor_aliases) {
      out[[alias_name]] <- factor(data[[source_var]])
    } else {
      analysis_var <- resolve_analysis_var(source_var, transformation_table)
      out[[alias_name]] <- as.numeric(data[[analysis_var]])
    }
  }
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out[stats::complete.cases(out), , drop = FALSE]
}

calculate_vif_from_matrix <- function(model_df, predictor_names) {
  design <- stats::model.matrix(
    stats::as.formula(paste("~", paste(predictor_names, collapse = " + "))),
    data = model_df
  )
  design <- design[, setdiff(colnames(design), "(Intercept)"), drop = FALSE]
  if (ncol(design) <= 1) {
    return(data.frame(term = colnames(design), vif = 1, tolerance = 1, stringsAsFactors = FALSE))
  }

  out <- do.call(
    rbind,
    lapply(seq_len(ncol(design)), function(i) {
      y <- design[, i]
      x <- design[, -i, drop = FALSE]
      fit <- stats::lm(y ~ x)
      r2 <- summary(fit)$r.squared
      vif <- 1 / (1 - r2)
      data.frame(
        term = colnames(design)[i],
        vif = vif,
        tolerance = 1 / vif,
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(out) <- NULL
  out
}

run_abeta_censored_models <- function(data, threshold, group_name = "Overall") {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for censored amyloid-beta sensitivity analysis.", call. = FALSE)
  }

  alias_map <- c(
    abeta = "S_ABETA",
    chp = "ChPICV",
    S_PTGENDER = "S_PTGENDER",
    S_AGE = "S_AGE",
    PTEDUCAT = "PTEDUCAT",
    APOE401 = "APOE401"
  )
  model_df <- build_alias_df(
    data = data,
    alias_map = alias_map,
    factor_aliases = c("S_PTGENDER", "APOE401"),
    transformation_table = NULL
  )
  if (nrow(model_df) < 20) {
    return(data.frame())
  }

  model_df$abeta_obs <- pmin(model_df$abeta, threshold)
  model_df$event_observed <- model_df$abeta < threshold
  surv_obj <- survival::Surv(model_df$abeta_obs, model_df$event_observed, type = "right")
  fit <- survival::survreg(
    surv_obj ~ chp + S_PTGENDER + S_AGE + PTEDUCAT + APOE401,
    data = model_df,
    dist = "gaussian"
  )
  coef_table <- summary(fit)$table
  out <- data.frame(
    group = group_name,
    term = rownames(coef_table),
    estimate = coef_table[, "Value"],
    std.error = coef_table[, "Std. Error"],
    statistic = coef_table[, "z"],
    p_value = coef_table[, "p"],
    threshold = threshold,
    n = nrow(model_df),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

run_matched_key_models <- function() {
  matched_files <- list.files(result_data_clean_dir, pattern = "^matched_.*\\.csv$", full.names = TRUE)
  if (length(matched_files) == 0) {
    return(data.frame())
  }

  rows <- list()
  for (matched_path in matched_files) {
    matched_data <- read_project_data(matched_path)
    pair_name <- sub("^matched_|\\.csv$", "", basename(matched_path))
    pair_levels <- strsplit(pair_name, "_vs_")[[1]]
    unmatched_pair <- analysis_data[as.character(analysis_data[[project_config$variables$group_label_var]]) %in% pair_levels, , drop = FALSE]

    for (sample_set in c("unmatched_pair", "matched_pair")) {
      sample_data <- if (sample_set == "matched_pair") matched_data else unmatched_pair
      for (cfg in project_config$robustness$matched_key_models) {
        fit_rows <- fit_linear_model_base(
          data = sample_data,
          outcome = cfg$outcome,
          exposure = cfg$exposure,
          covariates = project_config$variables$covariates,
          factor_vars = project_config$variables$categorical_covariates,
          transformation_table = transformation_plan
        )
        fit_rows <- fit_rows[fit_rows$term == "exposure", , drop = FALSE]
        fit_rows$pair_name <- pair_name
        fit_rows$sample_set <- sample_set
        fit_rows$model_name <- cfg$name
        rows[[paste(pair_name, sample_set, cfg$name, sep = "_")]] <- fit_rows
      }
    }
  }
  bind_rows_fill(rows)
}

run_collinearity_models <- function() {
  rows <- list()
  for (cfg in project_config$robustness$collinearity_models) {
    alias_map <- stats::setNames(cfg$predictors, cfg$predictors)
    alias_map <- c(outcome = cfg$outcome, alias_map)
    model_df <- build_alias_df(
      data = analysis_data,
      alias_map = alias_map,
      factor_aliases = project_config$variables$categorical_covariates,
      transformation_table = transformation_plan
    )
    if (nrow(model_df) == 0) {
      next
    }
    predictor_names <- cfg$predictors
    vif_rows <- calculate_vif_from_matrix(model_df = model_df, predictor_names = predictor_names)
    vif_rows$model_name <- cfg$name
    vif_rows$outcome <- cfg$outcome
    vif_rows$n <- nrow(model_df)
    rows[[cfg$name]] <- vif_rows
  }
  bind_rows_fill(rows)
}

run_alternative_chp_models <- function() {
  rows <- list()
  for (cfg in project_config$robustness$alternative_chp_models) {
    covariates <- unique(c(project_config$variables$covariates, cfg$extra_covariates))
    overall <- fit_linear_model_base(
      data = analysis_data,
      outcome = cfg$outcome,
      exposure = cfg$exposure,
      covariates = covariates,
      factor_vars = project_config$variables$categorical_covariates,
      transformation_table = transformation_plan
    )
    overall <- overall[overall$term == "exposure", , drop = FALSE]
    overall$group <- "Overall"
    overall$model_name <- cfg$name

    by_group <- fit_linear_models_by_group(
      data = analysis_data,
      group_var = project_config$variables$group_label_var,
      outcome = cfg$outcome,
      exposure = cfg$exposure,
      covariates = covariates,
      factor_vars = project_config$variables$categorical_covariates,
      transformation_table = transformation_plan
    )
    by_group <- by_group[by_group$term == "exposure", , drop = FALSE]
    by_group$model_name <- cfg$name
    rows[[cfg$name]] <- rbind(overall, by_group)
  }
  bind_rows_fill(rows)
}

run_structure_specificity_models <- function() {
  biomarker_vars <- c("MSD_STREM2CORRECTED", "TAU", "PTAU", "S_ABETA")
  rows <- list()
  for (cfg in project_config$robustness$structure_specificity_models) {
    for (bio in biomarker_vars) {
      fit_rows <- fit_linear_model_base(
        data = analysis_data,
        outcome = cfg$outcome,
        exposure = bio,
        covariates = unique(c(project_config$variables$covariates, cfg$add_covariates)),
        factor_vars = project_config$variables$categorical_covariates,
        transformation_table = transformation_plan
      )
      fit_rows <- fit_rows[fit_rows$term == "exposure", , drop = FALSE]
      fit_rows$structure <- cfg$name
      fit_rows$biomarker <- bio
      fit_rows$group <- "Overall"
      rows[[paste(cfg$name, bio, sep = "_")]] <- fit_rows
    }
  }
  bind_rows_fill(rows)
}

apply_fdr_to_summary_files <- function(summary_dir, family_map) {
  manifest_rows <- list()
  for (family_name in names(family_map)) {
    files <- family_map[[family_name]]
    for (file_name in files) {
      input_path <- file.path(summary_dir, file_name)
      if (!file.exists(input_path)) {
        next
      }
      df <- read_project_data(input_path)
      changed <- FALSE
      if ("p.value" %in% names(df)) {
        df$p_fdr <- stats::p.adjust(df$p.value, method = "fdr")
        changed <- TRUE
      }
      if ("p_value" %in% names(df)) {
        df$p_fdr <- stats::p.adjust(df$p_value, method = "fdr")
        changed <- TRUE
      }
      if ("omnibus_p" %in% names(df)) {
        df$omnibus_p_fdr <- stats::p.adjust(df$omnibus_p, method = "fdr")
        changed <- TRUE
      }
      if ("pairwise_p" %in% names(df)) {
        df$pairwise_p_fdr <- stats::p.adjust(df$pairwise_p, method = "fdr")
        changed <- TRUE
      }
      if (!changed) {
        next
      }
      output_name <- sub("\\.csv$", "_fdr.csv", file_name)
      output_path <- file.path(summary_dir, output_name)
      write_csv_utf8(df, output_path, row.names = FALSE)
      manifest_rows[[paste(family_name, output_name, sep = "_")]] <- data.frame(
        family = family_name,
        input_file = file_name,
        output_file = output_name,
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows_fill(manifest_rows)
}

abeta_threshold <- project_config$robustness$abeta_censoring$threshold %||% project_config$advanced_models$amyloid_ceiling_threshold

abeta_censored_rows <- bind_rows_fill(c(
  list(run_abeta_censored_models(analysis_data, threshold = abeta_threshold, group_name = "Overall")),
  lapply(c("CN", "MCI", "AD"), function(group_name) {
    sub_data <- analysis_data[as.character(analysis_data[[project_config$variables$group_label_var]]) == group_name, , drop = FALSE]
    run_abeta_censored_models(sub_data, threshold = abeta_threshold, group_name = group_name)
  })
))

matched_key_rows <- run_matched_key_models()
collinearity_rows <- run_collinearity_models()
alternative_chp_rows <- run_alternative_chp_models()
structure_specificity_rows <- run_structure_specificity_models()
fdr_manifest <- apply_fdr_to_summary_files(
  summary_dir = result_summary_dir,
  family_map = project_config$robustness$fdr_families
)

abeta_censored_path <- file.path(result_summary_dir, "abeta_censored_sensitivity.csv")
matched_key_path <- file.path(result_summary_dir, "matched_key_model_sensitivity.csv")
collinearity_path <- file.path(result_summary_dir, "collinearity_diagnostics.csv")
alternative_chp_path <- file.path(result_summary_dir, "alternative_chp_definition_models.csv")
structure_specificity_path <- file.path(result_summary_dir, "structure_specificity_models.csv")
fdr_manifest_path <- file.path(result_summary_dir, "fdr_adjustment_manifest.csv")

write_csv_utf8(abeta_censored_rows, abeta_censored_path, row.names = FALSE)
write_csv_utf8(matched_key_rows, matched_key_path, row.names = FALSE)
write_csv_utf8(collinearity_rows, collinearity_path, row.names = FALSE)
write_csv_utf8(alternative_chp_rows, alternative_chp_path, row.names = FALSE)
write_csv_utf8(structure_specificity_rows, structure_specificity_path, row.names = FALSE)
write_csv_utf8(fdr_manifest, fdr_manifest_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "06d_robustness_extensions",
  output_files = c(
    abeta_censored_path,
    matched_key_path,
    collinearity_path,
    alternative_chp_path,
    structure_specificity_path,
    fdr_manifest_path
  ),
  note = "Added FDR-adjusted summary outputs, right-censored amyloid-beta sensitivity models, matched-sample reruns of key main findings, collinearity diagnostics, alternative ChP definition analyses, and structure-specificity regression models.",
  summary_dir = result_summary_dir
)
