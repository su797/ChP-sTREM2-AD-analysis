plot_scatter_with_lm <- function(data, x, y, color_var = NULL, title = NULL) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    aes_args <- list(x = as.name(x), y = as.name(y))
    if (!is.null(color_var)) {
      aes_args$color <- as.name(color_var)
    }

    return(
      ggplot2::ggplot(data, do.call(ggplot2::aes, aes_args)) +
        ggplot2::geom_point(alpha = 0.8, size = 2.5) +
        ggplot2::geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = title, x = x, y = y)
    )
  }

  function() {
    x_value <- data[[x]]
    y_value <- data[[y]]
    plot(x_value, y_value, pch = 19, col = "#3b7ea1", xlab = x, ylab = y, main = title)
    abline(stats::lm(y_value ~ x_value), col = "#b54d3d", lwd = 2)
  }
}

plot_box_by_group <- function(data, x_group, y_value, title = NULL) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    return(
      ggplot2::ggplot(data, ggplot2::aes_string(x = x_group, y = y_value)) +
        ggplot2::geom_boxplot(outlier.alpha = 0.5) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(title = title, x = x_group, y = y_value)
    )
  }

  function() {
    graphics::boxplot(data[[y_value]] ~ data[[x_group]], xlab = x_group, ylab = y_value, main = title, col = "#d8e6ef")
  }
}

save_scatter_plot_base <- function(data, x, y, path, group_var = NULL, title = NULL, weights_var = NULL) {
  plot_fun <- function() {
    keep_cols <- c(x, y)
    if (!is.null(weights_var) && weights_var %in% names(data)) {
      keep_cols <- c(keep_cols, weights_var)
    }
    complete_idx <- stats::complete.cases(data[, keep_cols, drop = FALSE])
    plot_df <- data[complete_idx, , drop = FALSE]

    if (!is.null(group_var)) {
      group_factor <- factor(plot_df[[group_var]])
      cols <- c("#3b7ea1", "#b54d3d", "#5b8f3d", "#9467bd", "#8c564b")
      plot(plot_df[[x]], plot_df[[y]],
        pch = 19,
        col = cols[as.integer(group_factor)],
        xlab = x,
        ylab = y,
        main = title
      )
      graphics::legend("topright", legend = levels(group_factor), col = cols[seq_along(levels(group_factor))], pch = 19, bty = "n")
    } else {
      plot(plot_df[[x]], plot_df[[y]], pch = 19, col = "#3b7ea1", xlab = x, ylab = y, main = title)
    }

    if (nrow(plot_df) >= 3) {
      if (!is.null(weights_var) && weights_var %in% names(plot_df)) {
        plot_df <- plot_df[is.finite(plot_df[[weights_var]]) & plot_df[[weights_var]] > 0, , drop = FALSE]
      }
      if (nrow(plot_df) >= 3) {
        if (!is.null(weights_var) && weights_var %in% names(plot_df)) {
          graphics::abline(stats::lm(plot_df[[y]] ~ plot_df[[x]], weights = plot_df[[weights_var]]), col = "#b54d3d", lwd = 2)
        } else {
          graphics::abline(stats::lm(plot_df[[y]] ~ plot_df[[x]]), col = "#b54d3d", lwd = 2)
        }
      }
    }
  }

  save_plot_file(plot_fun, path = path, width = 7, height = 5, dpi = 300)
}
