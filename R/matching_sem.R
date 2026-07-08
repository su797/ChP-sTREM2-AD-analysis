greedy_match_pair <- function(data, group_var, pair_levels, match_vars, exact_vars = NULL, transformation_table = NULL, seed = 20260414) {
  set.seed(seed)

  sub_data <- data[as.character(data[[group_var]]) %in% pair_levels, , drop = FALSE]
  sub_data$.group <- as.character(sub_data[[group_var]])

  group_sizes <- table(sub_data$.group)
  anchor_group <- names(group_sizes)[which.min(group_sizes)]
  candidate_group <- setdiff(pair_levels, anchor_group)

  anchor_data <- sub_data[sub_data$.group == anchor_group, , drop = FALSE]
  candidate_data <- sub_data[sub_data$.group == candidate_group, , drop = FALSE]

  match_analysis_vars <- resolve_analysis_vars(match_vars, transformation_table)
  exact_analysis_vars <- if (!is.null(exact_vars)) exact_vars else character(0)

  matched_rows <- list()
  match_id <- 1
  candidate_available <- rep(TRUE, nrow(candidate_data))

  numeric_sds <- vapply(match_analysis_vars, function(var_name) {
    sd_value <- stats::sd(as.numeric(sub_data[[var_name]]), na.rm = TRUE)
    ifelse(is.na(sd_value) || sd_value == 0, 1, sd_value)
  }, FUN.VALUE = numeric(1))

  for (i in seq_len(nrow(anchor_data))) {
    anchor_row <- anchor_data[i, , drop = FALSE]
    candidate_idx <- which(candidate_available)
    if (length(candidate_idx) == 0) {
      next
    }

    candidate_pool <- candidate_data[candidate_idx, , drop = FALSE]
    if (!is.null(exact_vars) && length(exact_vars) > 0) {
      exact_mask <- apply(
        sapply(exact_analysis_vars, function(var_name) {
          candidate_pool[[var_name]] == anchor_row[[var_name]]
        }),
        1,
        all
      )
      candidate_pool <- candidate_pool[exact_mask, , drop = FALSE]
      candidate_idx <- candidate_idx[exact_mask]
    }

    if (nrow(candidate_pool) == 0) {
      next
    }

    distance_values <- vapply(seq_len(nrow(candidate_pool)), function(j) {
      candidate_row <- candidate_pool[j, , drop = FALSE]
      sum(
        ((as.numeric(anchor_row[match_analysis_vars]) - as.numeric(candidate_row[match_analysis_vars])) / numeric_sds)^2,
        na.rm = TRUE
      )
    }, FUN.VALUE = numeric(1))

    best_idx <- candidate_idx[which.min(distance_values)]
    candidate_available[best_idx] <- FALSE

    anchor_row$match_id <- match_id
    candidate_data[best_idx, "match_id"] <- match_id

    matched_rows[[length(matched_rows) + 1]] <- anchor_row
    matched_rows[[length(matched_rows) + 1]] <- candidate_data[best_idx, , drop = FALSE]
    match_id <- match_id + 1
  }

  matched_data <- if (length(matched_rows) > 0) {
    do.call(rbind, matched_rows)
  } else {
    sub_data[0, , drop = FALSE]
  }

  matched_data$pair_name <- paste(pair_levels, collapse = "_vs_")
  rownames(matched_data) <- NULL

  list(
    matched_data = matched_data,
    summary = data.frame(
      pair_name = paste(pair_levels, collapse = "_vs_"),
      anchor_group = anchor_group,
      candidate_group = candidate_group,
      matched_pairs = ifelse(nrow(matched_data) == 0, 0, length(unique(matched_data$match_id))),
      matched_n = nrow(matched_data),
      stringsAsFactors = FALSE
    )
  )
}

wrap_plot_label <- function(label, width = 12) {
  label <- as.character(label)
  label <- gsub("_", "_ ", label, fixed = TRUE)
  paste(strwrap(label, width = width), collapse = "\n")
}

format_p_for_sem_plot <- function(p_value, digits = 3) {
  if (is.null(p_value) || length(p_value) == 0 || is.na(p_value)) {
    return("NA")
  }
  threshold <- 10^(-digits)
  if (p_value < threshold) {
    return(paste0("<", format(threshold, scientific = FALSE, trim = TRUE)))
  }
  sprintf(paste0("%.", digits, "f"), round(p_value, digits))
}

make_sem_path_label <- function(estimate, p_value, prefix = "\u03b2", digits = 3, include_name = NULL) {
  if (is.null(estimate) || length(estimate) == 0 || is.na(estimate)) {
    return("NA")
  }
  beta_text <- sprintf("%s=%s", prefix, format_numeric_human(estimate, digits = digits))
  p_text <- paste0("p=", format_p_for_sem_plot(p_value, digits = digits))
  if (!is.null(include_name) && nzchar(include_name)) {
    return(sprintf("%s: %s, %s", include_name, beta_text, p_text))
  }
  sprintf("%s, %s", beta_text, p_text)
}

is_sem_plot_significant <- function(p_value) {
  !is.null(p_value) && length(p_value) > 0 && !is.na(p_value) && p_value < 0.05
}

compute_segment_label_position <- function(x0, y0, x1, y1, t = 0.5, offset = 0.03) {
  x <- x0 + t * (x1 - x0)
  y <- y0 + t * (y1 - y0)
  dx <- x1 - x0
  dy <- y1 - y0
  seg_len <- sqrt(dx^2 + dy^2)
  if (is.na(seg_len) || seg_len == 0) {
    return(list(x = x, y = y, angle = 0))
  }
  nx <- -dy / seg_len
  ny <- dx / seg_len
  list(
    x = x + offset * nx,
    y = y + offset * ny,
    angle = atan2(dy, dx) * 180 / pi
  )
}

normalize_text_angle <- function(angle) {
  if (is.na(angle)) {
    return(0)
  }
  if (angle > 90) {
    return(angle - 180)
  }
  if (angle < -90) {
    return(angle + 180)
  }
  angle
}

make_mediation_center_label <- function(indirect_value, indirect_p, proportion_value, direct_value = NA_real_) {
  if (!is_sem_plot_significant(indirect_p) || is.null(proportion_value) || length(proportion_value) == 0 || is.na(proportion_value)) {
    return(NULL)
  }

  opposite_direction <- !is.na(direct_value) && sign(indirect_value) != 0 && sign(direct_value) != 0 && sign(indirect_value) != sign(direct_value)
  if (opposite_direction) {
    percent_value <- -abs(100 * proportion_value)
    return(list(
      label = sprintf("Suppression\n(%s%%)", format_numeric_human(percent_value, digits = 1)),
      font = 2
    ))
  }

  percent_value <- abs(100 * proportion_value)
  list(
    label = sprintf("Facilitative\n(%s%%)", format_numeric_human(percent_value, digits = 1)),
    font = 2
  )
}

draw_mediation_diagram <- function(x_label, mediator_label, y_label, path_table, title = NULL, panel_cex = 2, margin_cex = 2) {
  a_row <- path_table[path_table$path == "a", , drop = FALSE]
  b_row <- path_table[path_table$path == "b", , drop = FALSE]
  c_row <- path_table[path_table$path == "c_prime", , drop = FALSE]
  indirect_row <- path_table[path_table$path == "indirect", , drop = FALSE]
  total_row <- path_table[path_table$path == "c_total", , drop = FALSE]
  prop_row <- path_table[path_table$path == "proportion_mediated", , drop = FALSE]

  a_value <- a_row$estimate[[1]]
  a_p <- a_row$p_value[[1]]
  b_value <- b_row$estimate[[1]]
  b_p <- b_row$p_value[[1]]
  c_prime_value <- c_row$estimate[[1]]
  c_prime_p <- c_row$p_value[[1]]
  indirect_value <- indirect_row$estimate[[1]]
  indirect_p <- indirect_row$p_value[[1]]
  total_value <- total_row$estimate[[1]]
  total_p <- total_row$p_value[[1]]
  prop_value <- prop_row$estimate[[1]]

  x_box_label <- wrap_plot_label(x_label, width = 12)
  mediator_box_label <- wrap_plot_label(mediator_label, width = 12)
  y_box_label <- wrap_plot_label(y_label, width = 12)

  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
  graphics::par(xpd = NA)

  if (!is.null(title)) {
    graphics::title(main = title, font.main = 2, cex.main = 1.12 * margin_cex, line = 0.55)
  }

  x_box <- c(xleft = 0.04, ybottom = 0.22, xright = 0.26, ytop = 0.44)
  mediator_box <- c(xleft = 0.42, ybottom = 0.58, xright = 0.64, ytop = 0.80)
  y_box <- c(xleft = 0.75, ybottom = 0.20, xright = 0.99, ytop = 0.46)

  graphics::rect(x_box["xleft"], x_box["ybottom"], x_box["xright"], x_box["ytop"], col = "#ffffff", border = "#222222", lwd = 6.2)
  graphics::rect(mediator_box["xleft"], mediator_box["ybottom"], mediator_box["xright"], mediator_box["ytop"], col = "#ffffff", border = "#222222", lwd = 6.2)
  graphics::rect(y_box["xleft"], y_box["ybottom"], y_box["xright"], y_box["ytop"], col = "#ffffff", border = "#222222", lwd = 6.2)

  graphics::text(mean(c(x_box["xleft"], x_box["xright"])), mean(c(x_box["ybottom"], x_box["ytop"])), labels = x_box_label, cex = 1.08 * panel_cex, font = 1)
  graphics::text(mean(c(mediator_box["xleft"], mediator_box["xright"])), mean(c(mediator_box["ybottom"], mediator_box["ytop"])), labels = mediator_box_label, cex = 1.08 * panel_cex, font = 1)
  graphics::text(mean(c(y_box["xleft"], y_box["xright"])), mean(c(y_box["ybottom"], y_box["ytop"])), labels = y_box_label, cex = 1.02 * panel_cex, font = 1)

  a_x0 <- x_box["xright"]
  a_y0 <- x_box["ytop"] - 0.01
  a_x1 <- mediator_box["xleft"]
  a_y1 <- mediator_box["ybottom"] + 0.02
  b_x0 <- mediator_box["xright"]
  b_y0 <- mediator_box["ybottom"] + 0.02
  b_x1 <- y_box["xleft"]
  b_y1 <- y_box["ytop"] - 0.01

  graphics::arrows(a_x0, a_y0, a_x1, a_y1, length = 0.07, lwd = 6.0)
  graphics::arrows(b_x0, b_y0, b_x1, b_y1, length = 0.07, lwd = 6.0)
  direct_y <- 0.25
  total_y <- 0.19
  graphics::arrows(x_box["xright"], direct_y, y_box["xleft"], direct_y, length = 0.07, lwd = 6.3)
  graphics::arrows(x_box["xright"], total_y, y_box["xleft"], total_y, length = 0.07, lwd = 5.2, lty = 2)

  graphics::text(
    0.30, 0.655,
    labels = make_sem_path_label(a_value, a_p),
    cex = 0.82 * panel_cex,
    font = ifelse(is_sem_plot_significant(a_p), 2, 1),
    srt = 0,
    adj = c(1, 0.5)
  )
  graphics::text(
    0.71, 0.655,
    labels = make_sem_path_label(b_value, b_p),
    cex = 0.82 * panel_cex,
    font = ifelse(is_sem_plot_significant(b_p), 2, 1),
    srt = 0,
    adj = c(0, 0.5)
  )
  graphics::text(
    0.50, direct_y + 0.03,
    labels = sprintf("Direct (c'): %s", make_sem_path_label(c_prime_value, c_prime_p)),
    cex = 0.92 * panel_cex,
    font = ifelse(is_sem_plot_significant(c_prime_p), 2, 1)
  )
  graphics::text(
    0.50, total_y - 0.03,
    labels = sprintf("Total (c): %s", make_sem_path_label(total_value, total_p)),
    cex = 0.92 * panel_cex,
    font = ifelse(is_sem_plot_significant(total_p), 2, 1)
  )

  center_label <- make_mediation_center_label(indirect_value, indirect_p, prop_value, c_prime_value)
  if (!is.null(center_label)) {
    graphics::text(
      0.50, 0.48,
      labels = center_label$label,
      cex = 1.22 * panel_cex,
      font = center_label$font
    )
  }
}

run_mediation_model <- function(data, x_var, mediator_var, y_var, covariates = NULL, factor_vars = NULL, transformation_table = NULL, n_boot = 2000, seed = 20260414) {
  x_analysis <- resolve_analysis_var(x_var, transformation_table)
  mediator_analysis <- resolve_analysis_var(mediator_var, transformation_table)
  cov_analysis <- resolve_analysis_vars(covariates, transformation_table)

  model_df <- data.frame(
    x = z_score(data[[x_analysis]]),
    mediator = z_score(data[[mediator_analysis]]),
    y = z_score(data[[y_var]])
  )

  if (!is.null(covariates) && length(covariates) > 0) {
    for (i in seq_along(covariates)) {
      cov_name <- covariates[i]
      analysis_name <- cov_analysis[i]
      if (!is.null(factor_vars) && cov_name %in% factor_vars) {
        model_df[[cov_name]] <- factor(data[[cov_name]])
      } else {
        model_df[[cov_name]] <- as.numeric(data[[analysis_name]])
      }
    }
  }

  model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
  rhs_cov <- if (!is.null(covariates) && length(covariates) > 0) paste(covariates, collapse = " + ") else NULL
  mediator_formula <- if (is.null(rhs_cov)) "mediator ~ x" else paste("mediator ~ x +", rhs_cov)
  outcome_formula <- if (is.null(rhs_cov)) "y ~ x + mediator" else paste("y ~ x + mediator +", rhs_cov)
  total_formula <- if (is.null(rhs_cov)) "y ~ x" else paste("y ~ x +", rhs_cov)

  fit_a <- stats::lm(stats::as.formula(mediator_formula), data = model_df)
  fit_b <- stats::lm(stats::as.formula(outcome_formula), data = model_df)
  fit_total <- stats::lm(stats::as.formula(total_formula), data = model_df)

  a_path <- summary(fit_a)$coefficients["x", , drop = FALSE]
  b_path <- summary(fit_b)$coefficients["mediator", , drop = FALSE]
  c_prime <- summary(fit_b)$coefficients["x", , drop = FALSE]
  c_total <- summary(fit_total)$coefficients["x", , drop = FALSE]

  indirect_value <- a_path[1, 1] * b_path[1, 1]
  total_value <- c_total[1, 1]
  proportion_mediated <- ifelse(total_value == 0, NA_real_, indirect_value / total_value)

  set.seed(seed)
  boot_values <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(model_df)), replace = TRUE)
    boot_df <- model_df[idx, , drop = FALSE]
    boot_a <- stats::lm(stats::as.formula(mediator_formula), data = boot_df)
    boot_b <- stats::lm(stats::as.formula(outcome_formula), data = boot_df)
    stats::coef(boot_a)[["x"]] * stats::coef(boot_b)[["mediator"]]
  })

  indirect_ci <- stats::quantile(boot_values, probs = c(0.025, 0.975), na.rm = TRUE)
  indirect_p <- 2 * min(mean(boot_values <= 0), mean(boot_values >= 0))

  path_table <- data.frame(
    path = c("a", "b", "c_prime", "c_total", "indirect", "proportion_mediated"),
    estimate = c(a_path[1, 1], b_path[1, 1], c_prime[1, 1], c_total[1, 1], indirect_value, proportion_mediated),
    p_value = c(a_path[1, 4], b_path[1, 4], c_prime[1, 4], c_total[1, 4], indirect_p, NA_real_),
    conf_low = c(NA_real_, NA_real_, NA_real_, NA_real_, indirect_ci[[1]], NA_real_),
    conf_high = c(NA_real_, NA_real_, NA_real_, NA_real_, indirect_ci[[2]], NA_real_),
    stringsAsFactors = FALSE
  )

  list(
    model_df = model_df,
    path_table = path_table,
    fits = list(mediator = fit_a, outcome = fit_b, total = fit_total)
  )
}

plot_mediation_diagram <- function(x_label, mediator_label, y_label, path_table, path, title = NULL) {
  save_plot_file(function() {
    draw_mediation_diagram(
      x_label = x_label,
      mediator_label = mediator_label,
      y_label = y_label,
      path_table = path_table,
      title = title,
      panel_cex = 2.60,
      margin_cex = 2.55
    )
  }, path = path, width = 14.0, height = 9.0, dpi = 320)
}

plot_mediation_panel <- function(diagram_specs, path, title = NULL, ncol = 3) {
  if (length(diagram_specs) == 0) {
    return(invisible(NULL))
  }

  n_panels <- length(diagram_specs)
  nrow <- ceiling(n_panels / ncol)
  width <- max(26, ncol * 8.8)
  height <- max(18, nrow * 7.2 + ifelse(is.null(title), 0, 1.9))

  save_plot_file(function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(nrow, ncol), mar = c(2.0, 1.6, 3.0, 1.6), oma = c(0, 0, ifelse(is.null(title), 0, 2.4), 0))

    for (spec in diagram_specs) {
      panel_title <- spec$title
      if (!is.null(spec$row_label) && !isTRUE(spec$row_label == "")) {
        panel_title <- paste(spec$row_label, spec$title, sep = " | ")
      }

      draw_mediation_diagram(
        x_label = spec$x_label,
        mediator_label = spec$mediator_label,
        y_label = spec$y_label,
        path_table = spec$path_table,
        title = panel_title,
        panel_cex = 2.35,
        margin_cex = 2.30
      )
    }

    if (!is.null(title)) {
      graphics::mtext(title, outer = TRUE, side = 3, line = 0.35, cex = 2.45, font = 2)
    }
  }, path = path, width = width, height = height, dpi = 320)
}
