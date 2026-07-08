format_p_value <- function(x, digits = 3) {
  ifelse(is.na(x), NA, ifelse(x < 0.001, "<0.001", format(round(x, digits), nsmall = digits)))
}

prepare_regression_table <- function(model_results, digits = 3) {
  out <- model_results

  if ("estimate" %in% names(out)) {
    out$estimate <- round(out$estimate, digits)
  }
  if ("std.error" %in% names(out)) {
    out$std.error <- round(out$std.error, digits)
  }
  if ("statistic" %in% names(out)) {
    out$statistic <- round(out$statistic, digits)
  }
  if ("conf.low" %in% names(out)) {
    out$conf.low <- round(out$conf.low, digits)
  }
  if ("conf.high" %in% names(out)) {
    out$conf.high <- round(out$conf.high, digits)
  }
  if ("or" %in% names(out)) {
    out$or <- round(out$or, digits)
  }
  if ("or_ci_low" %in% names(out)) {
    out$or_ci_low <- round(out$or_ci_low, digits)
  }
  if ("or_ci_high" %in% names(out)) {
    out$or_ci_high <- round(out$or_ci_high, digits)
  }
  if ("p.value" %in% names(out)) {
    out$p_value_display <- format_p_value(out$p.value, digits)
  }

  out
}

export_three_line_table <- function(data, csv_path, html_path = NULL, title = NULL) {
  write_csv_utf8(data, csv_path, row.names = FALSE)

  if (!is.null(html_path)) {
    if (requireNamespace("gt", quietly = TRUE)) {
      gt_tbl <- gt::gt(data)
      if (!is.null(title)) {
        gt_tbl <- gt::tab_header(gt_tbl, title = title)
      }
      gt_tbl <- gt::opt_table_outline(gt_tbl)
      gt::gtsave(gt_tbl, html_path)
    } else {
      write_simple_html_table(data, html_path, title = title)
    }
  }

  invisible(
    list(
      csv = csv_path,
      html = html_path
    )
  )
}
