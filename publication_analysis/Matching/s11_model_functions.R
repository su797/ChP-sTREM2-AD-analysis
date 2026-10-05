# CSV 表头可能带 BOM，读取后先去掉，避免同名变量匹配失败。
strip_bom <-
function (x) 
{
    sub("^﻿", "", x)
}
read_project_data <-
function (path, sheet = NULL) 
{
    ext <- tolower(tools::file_ext(path))
    if (ext == "csv") {
        data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, 
            na.strings = c("", "NA", "#NULL!", "#N/A", "NULL"), 
            fileEncoding = "UTF-8-BOM")
        names(data) <- strip_bom(names(data))
        return(data)
    }
    if (ext %in% c("xlsx", "xls")) {
        ensure_packages("readxl")
        data <- readxl::read_excel(path, sheet = sheet)
        names(data) <- strip_bom(names(data))
        return(data)
    }
    if (ext == "rds") {
        return(readRDS(path))
    }
    stop("Unsupported file format: ", ext, call. = FALSE)
}
ensure_packages <-
function (pkgs) 
{
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, 
        FUN.VALUE = logical(1))]
    if (length(missing_pkgs) > 0) {
        stop(sprintf("Missing required packages: %s. Please install them first.", 
            paste(missing_pkgs, collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}
# 解析用の変数名は変換表から取得する。ここで新たに変換を推測しない。
resolve_analysis_var <-
function (var_name, transformation_table = NULL) 
{
    if (is.null(transformation_table)) {
        return(var_name)
    }
    hit <- transformation_table[transformation_table$variable == 
        var_name, , drop = FALSE]
    if (nrow(hit) == 0) {
        return(var_name)
    }
    hit$analysis_var[[1]]
}
resolve_analysis_vars <-
function (var_names, transformation_table = NULL) 
{
    vapply(var_names, resolve_analysis_var, transformation_table = transformation_table, 
        FUN.VALUE = character(1))
}
build_formula_text <-
function (outcome, predictors, factor_vars = NULL) 
{
    predictor_terms <- vapply(predictors, function(var_name) {
        if (!is.null(factor_vars) && var_name %in% factor_vars) {
            sprintf("factor(%s)", var_name)
        }
        else {
            var_name
        }
    }, FUN.VALUE = character(1))
    sprintf("%s ~ %s", outcome, paste(predictor_terms, collapse = " + "))
}
tidy_lm_base <-
function (fit, model_type = "linear", outcome = NULL, exposure = NULL, 
    formula_text = NULL) 
{
    coef_table <- summary(fit)$coefficients
    ci_table <- suppressWarnings(stats::confint(fit))
    out <- data.frame(term = rownames(coef_table), estimate = coef_table[, 
        1], std.error = coef_table[, 2], statistic = coef_table[, 
        3], p.value = coef_table[, 4], conf.low = ci_table[, 
        1], conf.high = ci_table[, 2], stringsAsFactors = FALSE)
    out$model_type <- model_type
    out$outcome <- outcome
    out$exposure <- exposure
    out$formula <- formula_text
    out
}
fit_linear_model_base <-
function (data, outcome, exposure, covariates = NULL, factor_vars = NULL, 
    transformation_table = NULL) 
{
    outcome_var <- resolve_analysis_var(outcome, transformation_table)
    exposure_var <- resolve_analysis_var(exposure, transformation_table)
    covariate_vars <- resolve_analysis_vars(covariates, transformation_table)
    model_df <- data.frame(outcome = as.numeric(data[[outcome_var]]), 
        exposure = as.numeric(data[[exposure_var]]))
    predictors <- "exposure"
    if (!is.null(covariates) && length(covariates) > 0) {
        for (i in seq_along(covariates)) {
            cov_name <- covariates[i]
            analysis_name <- covariate_vars[i]
            if (!is.null(factor_vars) && cov_name %in% factor_vars) {
                model_df[[cov_name]] <- factor(data[[cov_name]])
            }
            else {
                model_df[[cov_name]] <- as.numeric(data[[analysis_name]])
            }
            predictors <- c(predictors, cov_name)
        }
    }
# 完整病例按当前模型的实际变量确定，下面的 n 也以此为准。
    model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
    formula_text <- build_formula_text("outcome", predictors, 
        factor_vars = factor_vars)
    fit <- stats::lm(stats::as.formula(formula_text), data = model_df)
    out <- tidy_lm_base(fit = fit, model_type = "linear", outcome = outcome, 
        exposure = exposure, formula_text = formula_text)
# 標準化係数は、このモデルに入った対象者の標準偏差で計算する。
    exposure_beta <- out$estimate[out$term == "exposure"]
    out$std_beta <- NA_real_
    if (length(exposure_beta) == 1 && stats::sd(model_df$outcome) > 
        0) {
        out$std_beta[out$term == "exposure"] <- exposure_beta * 
            stats::sd(model_df$exposure)/stats::sd(model_df$outcome)
    }
    out$n <- nrow(model_df)
    out
}
bind_rows_fill <-
function (items) 
{
    items <- items[vapply(items, function(x) !is.null(x) && nrow(x) > 
        0, FUN.VALUE = logical(1))]
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
# 从 matched_CN_vs_AD.csv 这类文件名读取两个组名，而不是从单个组名猜比较对象。
parse_matched_pair <-
function (matched_path, allowed_groups) 
{
    filename <- basename(matched_path)
    if (!grepl("^matched_.+\\.csv$", filename)) {
        stop("Invalid matched-data filename: ", filename, call. = FALSE)
    }
    pair_name <- sub("^matched_", "", tools::file_path_sans_ext(filename))
    pair_levels <- strsplit(pair_name, "_vs_", fixed = TRUE)[[1]]
    if (length(pair_levels) != 2L || anyDuplicated(pair_levels) || 
        any(!pair_levels %in% allowed_groups)) {
        stop("Matched filename must identify two distinct known groups: ", 
            filename, call. = FALSE)
    }
    list(name = pair_name, levels = pair_levels)
}
run_matched_key_models <-
function () 
{
    matched_files <- list.files(result_data_clean_dir, pattern = "^matched_.*\\.csv$", 
        full.names = TRUE)
    if (length(matched_files) == 0) {
        return(data.frame())
    }
    rows <- list()
    for (matched_path in matched_files) {
        matched_data <- read_project_data(matched_path)
        group_var <- project_config$variables$group_label_var
        pair <- parse_matched_pair(matched_path, unique(as.character(analysis_data[[group_var]])))
        pair_name <- pair$name
        pair_levels <- pair$levels
# 未匹配参照必须同时包含这一对诊断组；不能只取其中一组，也不能混入第三组。
        unmatched_pair <- analysis_data[as.character(analysis_data[[project_config$variables$group_label_var]]) %in% 
            pair_levels, , drop = FALSE]
        for (sample_set in c("unmatched_pair", "matched_pair")) {
            sample_data <- if (sample_set == "matched_pair") 
                matched_data
            else unmatched_pair
# 両群が含まれていることを確認してから回帰に進む。
            observed_groups <- unique(as.character(sample_data[[group_var]]))
            if (anyNA(observed_groups) || !setequal(observed_groups, 
                pair_levels)) {
                stop("Both requested diagnostic groups must be present in ", 
                  pair_name, " / ", sample_set, call. = FALSE)
            }
            for (cfg in project_config$robustness$matched_key_models) {
                fit_rows <- fit_linear_model_base(data = sample_data, 
                  outcome = cfg$outcome, exposure = cfg$exposure, 
                  covariates = project_config$variables$covariates, 
                  factor_vars = project_config$variables$categorical_covariates, 
                  transformation_table = transformation_plan)
                fit_rows <- fit_rows[fit_rows$term == "exposure", 
                  , drop = FALSE]
                fit_rows$pair_name <- pair_name
                fit_rows$sample_set <- sample_set
                fit_rows$model_name <- cfg$name
                rows[[paste(pair_name, sample_set, cfg$name, 
                  sep = "_")]] <- fit_rows
            }
        }
    }
    bind_rows_fill(rows)
}
