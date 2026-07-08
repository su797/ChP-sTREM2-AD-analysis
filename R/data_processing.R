safe_character <- function(x) {
  trimws(strip_bom(as.character(x)))
}

safe_numeric <- function(x) {
  x <- safe_character(x)
  x[x %in% c("", "NA", "#NULL!", "#N/A", "NULL", "NaN")] <- NA
  suppressWarnings(as.numeric(x))
}

recode_from_map <- function(x, value_map, ordered = FALSE) {
  x_chr <- safe_character(x)
  out <- unname(value_map[x_chr])
  factor(out, levels = unique(unname(value_map)), ordered = ordered)
}

convert_selected_to_numeric <- function(data, numeric_vars) {
  conversion_log <- do.call(
    rbind,
    lapply(numeric_vars, function(var_name) {
      raw_value <- data[[var_name]]
      converted <- safe_numeric(raw_value)
      bad_count <- sum(
        !is.na(safe_character(raw_value)) &
          safe_character(raw_value) != "" &
          is.na(converted)
      )
      data[[var_name]] <<- converted
      data.frame(
        variable = var_name,
        introduced_na = bad_count,
        stringsAsFactors = FALSE
      )
    })
  )

  list(data = data, conversion_log = conversion_log)
}

evaluate_normality <- function(x, alpha = 0.05) {
  x <- x[!is.na(x)]
  n <- length(x)

  if (n < 3) {
    return(data.frame(
      n = n,
      shapiro_p = NA_real_,
      skewness = NA_real_,
      is_normal = NA,
      can_log = FALSE,
      can_log1p = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  shapiro_p <- tryCatch(stats::shapiro.test(x)[["p.value"]], error = function(e) NA_real_)
  skewness <- if (stats::sd(x) > 0) {
    mean((x - mean(x))^3) / stats::sd(x)^3
  } else {
    0
  }

  data.frame(
    n = n,
    shapiro_p = shapiro_p,
    skewness = skewness,
    is_normal = !is.na(shapiro_p) && shapiro_p >= alpha,
    can_log = all(x > 0),
    can_log1p = all(x >= 0),
    stringsAsFactors = FALSE
  )
}

build_transformation_plan <- function(data, continuous_vars, alpha = 0.05, prefix = "L_") {
  out <- do.call(
    rbind,
    lapply(continuous_vars, function(var_name) {
      norm_info <- evaluate_normality(data[[var_name]], alpha = alpha)
      transform_method <- "none"
      analysis_var <- var_name
      note <- ""

      if (!isTRUE(norm_info$is_normal)) {
        if (isTRUE(norm_info$can_log)) {
          transform_method <- "log"
          analysis_var <- paste0(prefix, var_name)
        } else if (isTRUE(norm_info$can_log1p)) {
          transform_method <- "log1p"
          analysis_var <- paste0(prefix, var_name)
        } else {
          note <- "Non-normal but cannot log-transform because values are negative."
        }
      }

      cbind(
        data.frame(
          variable = var_name,
          analysis_var = analysis_var,
          transform_method = transform_method,
          note = note,
          stringsAsFactors = FALSE
        ),
        norm_info
      )
    })
  )

  rownames(out) <- NULL
  out
}

apply_transformation_plan <- function(data, plan_table) {
  transformed_table <- plan_table
  transformed_table$transformed_p <- NA_real_

  for (i in seq_len(nrow(plan_table))) {
    row_info <- plan_table[i, , drop = FALSE]
    var_name <- row_info$variable
    analysis_var <- row_info$analysis_var
    transform_method <- row_info$transform_method

    if (transform_method == "log") {
      data[[analysis_var]] <- log(data[[var_name]])
      transformed_table$transformed_p[i] <- evaluate_normality(data[[analysis_var]])$shapiro_p
    } else if (transform_method == "log1p") {
      data[[analysis_var]] <- log1p(data[[var_name]])
      transformed_table$transformed_p[i] <- evaluate_normality(data[[analysis_var]])$shapiro_p
    } else {
      transformed_table$transformed_p[i] <- transformed_table$shapiro_p[i]
    }
  }

  list(data = data, transformation_table = transformed_table)
}

resolve_analysis_var <- function(var_name, transformation_table = NULL) {
  if (is.null(transformation_table)) {
    return(var_name)
  }

  hit <- transformation_table[transformation_table$variable == var_name, , drop = FALSE]
  if (nrow(hit) == 0) {
    return(var_name)
  }

  hit$analysis_var[[1]]
}

resolve_analysis_vars <- function(var_names, transformation_table = NULL) {
  vapply(var_names, resolve_analysis_var, transformation_table = transformation_table, FUN.VALUE = character(1))
}

z_score <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

create_cognition_composite <- function(data, composite_name, components, directions, transformation_table = NULL, min_non_missing = NULL) {
  if (length(components) != length(directions)) {
    stop("components and directions must have the same length.", call. = FALSE)
  }

  if (is.null(min_non_missing)) {
    min_non_missing <- length(components)
  }

  component_map <- do.call(
    rbind,
    lapply(seq_along(components), function(i) {
      analysis_var <- resolve_analysis_var(components[i], transformation_table)
      data.frame(
        composite = composite_name,
        component = components[i],
        analysis_var = analysis_var,
        direction = directions[i],
        stringsAsFactors = FALSE
      )
    })
  )

  z_matrix <- do.call(
    cbind,
    lapply(seq_along(components), function(i) {
      analysis_var <- component_map$analysis_var[i]
      z_score(directions[i] * data[[analysis_var]])
    })
  )

  composite_score <- apply(z_matrix, 1, function(row_value) {
    if (sum(!is.na(row_value)) >= min_non_missing) {
      mean(row_value, na.rm = TRUE)
    } else {
      NA_real_
    }
  })

  data[[composite_name]] <- composite_score
  list(data = data, component_map = component_map)
}
