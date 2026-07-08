escape_md_cell <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\|", "\\\\|", x)
  x
}

significance_mark_suffix <- function(p_value) {
  if (is.null(p_value) || length(p_value) == 0 || is.na(p_value)) {
    return("")
  }
  if (p_value < 0.001) {
    return("\\*\\*")
  }
  if (p_value < 0.05) {
    return("\\*")
  }
  ""
}

style_significant_markdown <- function(text_value, p_value) {
  if (is.null(p_value) || length(p_value) == 0 || is.na(p_value) || p_value >= 0.05) {
    return(text_value)
  }
  paste0("**", text_value, "**", significance_mark_suffix(p_value))
}

format_report_value <- function(x, digits = 4) {
  if (is.numeric(x)) {
    out <- ifelse(
      is.na(x),
      "",
      ifelse(abs(x) >= 1000, format(round(x, 2), big.mark = ",", trim = TRUE), format(round(x, digits), trim = TRUE))
    )
    return(out)
  }
  escape_md_cell(x)
}

is_p_value_column <- function(col_name) {
  grepl("(^p\\.value$|^p_value$|_p$|_p_value$|_p\\.value$)", col_name)
}

markdown_table <- function(data, digits = 4, empty_text = "_No data available._") {
  if (is.null(data) || nrow(data) == 0) {
    return(c(empty_text, ""))
  }

  formatted <- data
  for (col_name in names(formatted)) {
    if (is.numeric(formatted[[col_name]]) && is_p_value_column(col_name)) {
      formatted[[col_name]] <- vapply(
        formatted[[col_name]],
        function(p_value) format_p_value_human(p_value, digits = digits, markdown = TRUE),
        FUN.VALUE = character(1)
      )
    } else {
      formatted[[col_name]] <- format_report_value(formatted[[col_name]], digits = digits)
    }
  }

  header <- paste0("| ", paste(names(formatted), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(formatted)), collapse = " | "), " |")
  rows <- apply(formatted, 1, function(row_value) {
    paste0("| ", paste(row_value, collapse = " | "), " |")
  })

  c(header, separator, rows, "")
}

markdown_link <- function(label, target) {
  sprintf("[%s](%s)", label, target)
}

markdown_image <- function(path, alt = NULL) {
  if (is.null(alt)) {
    alt <- basename(path)
  }
  sprintf("![%s](%s)", alt, path)
}

markdown_bullet_links <- function(paths, relative_dir = ".", empty_text = "- _No files found._") {
  if (length(paths) == 0) {
    return(c(empty_text, ""))
  }

  c(
    vapply(paths, function(path) {
      rel_path <- file.path(relative_dir, basename(path))
      sprintf("- %s", markdown_link(basename(path), rel_path))
    }, FUN.VALUE = character(1)),
    ""
  )
}

pick_term_rows <- function(data, term_name = "exposure") {
  if (is.null(data) || !"term" %in% names(data)) {
    return(data.frame())
  }
  data[data$term == term_name, , drop = FALSE]
}

read_result_or_empty <- function(path) {
  if (!file.exists(path)) {
    return(data.frame())
  }
  file_info <- file.info(path)
  if (is.na(file_info$size) || file_info$size == 0) {
    return(data.frame())
  }
  raw_lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  stripped_lines <- trimws(gsub("[\",]", "", raw_lines))
  if (length(raw_lines) == 0 || !any(nzchar(stripped_lines))) {
    return(data.frame())
  }
  read_project_data(path)
}

collapse_values <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return("None")
  }
  paste(x, collapse = "; ")
}

format_p_value_human <- function(p_value, digits = 4, markdown = TRUE) {
  if (is.null(p_value) || length(p_value) == 0 || is.na(p_value)) {
    return("NA")
  }

  threshold <- 10^(-digits)
  base_text <- if (p_value < threshold) {
    paste0("< ", format(threshold, scientific = FALSE, trim = TRUE))
  } else {
    sprintf(paste0("%.", digits, "f"), round(p_value, digits))
  }

  if (!isTRUE(markdown)) {
    return(base_text)
  }

  style_significant_markdown(base_text, p_value)
}

format_numeric_human <- function(value, digits = 4) {
  if (is.null(value) || length(value) == 0 || is.na(value)) {
    return("NA")
  }
  sprintf(paste0("%.", digits, "f"), round(value, digits))
}

format_beta_human <- function(estimate, std_beta = NA, digits = 4) {
  if (is.null(std_beta) || length(std_beta) == 0 || is.na(std_beta)) {
    beta_value <- estimate
  } else {
    beta_value <- std_beta
  }
  format_numeric_human(beta_value, digits = digits)
}

order_by_group <- function(data, group_col, group_order = NULL) {
  if (is.null(data) || nrow(data) == 0 || is.null(group_order) || !group_col %in% names(data)) {
    return(data)
  }
  data[[group_col]] <- factor(as.character(data[[group_col]]), levels = group_order)
  data <- data[order(data[[group_col]]), , drop = FALSE]
  data[[group_col]] <- as.character(data[[group_col]])
  data
}

make_paragraph_block <- function(sentences, empty_text = NULL) {
  sentences <- sentences[!is.na(sentences) & sentences != ""]
  if (length(sentences) == 0) {
    if (is.null(empty_text)) {
      return(character(0))
    }
    return(c(empty_text, ""))
  }
  c(sentences, "")
}
