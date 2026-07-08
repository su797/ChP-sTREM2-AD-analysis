ensure_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing_pkgs) > 0) {
    stop(
      sprintf(
        "Missing required packages: %s. Please install them first.",
        paste(missing_pkgs, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

load_common_functions <- function(common_r_dir) {
  r_files <- list.files(common_r_dir, pattern = "\\.R$", full.names = TRUE)
  r_files <- r_files[basename(r_files) != "utils_io.R"]
  invisible(lapply(r_files, source))
}

strip_bom <- function(x) {
  sub("^\ufeff", "", x)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_latest_result_pointer <- function(result_context, pointer_name = "LATEST.txt") {
  pointer_path <- file.path(result_context$result_base_dir, pointer_name)
  pointer_lines <- c(
    paste0("run_id=", result_context$run_id),
    paste0("run_dir=", result_context$run_dir)
  )
  writeLines(pointer_lines, con = pointer_path, useBytes = TRUE)
  invisible(pointer_path)
}

create_result_dirs <- function(project_root, result_config = list()) {
  default_config <- list(
    dir_time_format = "%Y%m%d_%H%M%S",
    data_subdir = "data_clean",
    write_latest_pointer = TRUE,
    latest_pointer_name = "LATEST.txt",
    run_id = NULL
  )
  cfg <- utils::modifyList(default_config, result_config)

  result_base_dir <- ensure_dir(file.path(project_root, "result"))
  if (!is.null(cfg$run_id) && nzchar(cfg$run_id)) {
    run_id <- cfg$run_id
  } else {
    run_id_base <- format(Sys.time(), cfg$dir_time_format)
    run_id <- run_id_base
    suffix_idx <- 1
    while (dir.exists(file.path(result_base_dir, run_id))) {
      run_id <- sprintf("%s_%02d", run_id_base, suffix_idx)
      suffix_idx <- suffix_idx + 1
    }
  }
  run_dir <- ensure_dir(file.path(result_base_dir, run_id))

  result_context <- list(
    result_base_dir = result_base_dir,
    run_id = run_id,
    run_dir = run_dir,
    summary_dir = ensure_dir(file.path(run_dir, "summary")),
    tables_dir = ensure_dir(file.path(run_dir, "tables")),
    figures_dir = ensure_dir(file.path(run_dir, "figures")),
    report_dir = ensure_dir(file.path(run_dir, "report")),
    data_clean_dir = ensure_dir(file.path(run_dir, cfg$data_subdir)),
    latest_clean_dir = ensure_dir(file.path(project_root, "data", "clean"))
  )

  if (isTRUE(cfg$write_latest_pointer)) {
    write_latest_result_pointer(result_context, pointer_name = cfg$latest_pointer_name)
  }

  invisible(result_context)
}

read_project_data <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))

  if (ext == "csv") {
    data <- read.csv(
      path,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("", "NA", "#NULL!", "#N/A", "NULL"),
      fileEncoding = "UTF-8-BOM"
    )
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

write_csv_utf8 <- function(data, path, row.names = FALSE) {
  ensure_dir(dirname(path))
  utils::write.csv(data, file = path, row.names = row.names, fileEncoding = "UTF-8")
  invisible(path)
}

save_plot_file <- function(plot_obj, path, width = 7, height = 5, dpi = 300) {
  ensure_dir(dirname(path))

  if (requireNamespace("ggplot2", quietly = TRUE) && inherits(plot_obj, "ggplot")) {
    ggplot2::ggsave(filename = path, plot = plot_obj, width = width, height = height, dpi = dpi)
    return(invisible(path))
  }

  if (!is.function(plot_obj)) {
    stop("plot_obj must be a ggplot object or a plotting function.", call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext == "pdf") {
    grDevices::pdf(path, width = width, height = height)
  } else {
    grDevices::png(path, width = width, height = height, units = "in", res = dpi)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_obj()
  invisible(path)
}

append_analysis_log <- function(project_root = NULL, analysis_name, output_files, note = NULL, summary_dir = NULL) {
  if (is.null(summary_dir)) {
    if (is.null(project_root)) {
      stop("Either project_root or summary_dir must be provided.", call. = FALSE)
    }
    summary_dir <- file.path(project_root, "result", "summary")
  }

  ensure_dir(summary_dir)
  log_path <- file.path(summary_dir, "analysis_log.csv")
  log_row <- data.frame(
    run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    analysis_name = analysis_name,
    output_files = paste(normalizePath(output_files, winslash = "/", mustWork = FALSE), collapse = "; "),
    note = ifelse(is.null(note), "", note),
    stringsAsFactors = FALSE
  )

  if (file.exists(log_path)) {
    old_log <- read.csv(log_path, stringsAsFactors = FALSE, check.names = FALSE)
    new_log <- rbind(old_log, log_row)
  } else {
    new_log <- log_row
  }

  write_csv_utf8(new_log, log_path, row.names = FALSE)
}

standardize_model_output <- function(model_results, analysis_name, outcome = NULL) {
  model_results$analysis_name <- analysis_name
  if (!is.null(outcome)) {
    model_results$outcome <- outcome
  }
  model_results
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

write_simple_html_table <- function(data, path, title = NULL) {
  ensure_dir(dirname(path))

  header_html <- paste(sprintf("<th>%s</th>", html_escape(names(data))), collapse = "")
  body_html <- apply(data, 1, function(row) {
    paste0(
      "<tr>",
      paste(sprintf("<td>%s</td>", html_escape(as.character(row))), collapse = ""),
      "</tr>"
    )
  })

  lines <- c(
    "<!doctype html>",
    "<html>",
    "<head>",
    "<meta charset='utf-8'>",
    "<style>",
    "body { font-family: Arial, sans-serif; margin: 24px; }",
    "h2 { margin-bottom: 12px; }",
    "table { border-collapse: collapse; width: 100%; }",
    "thead tr { border-top: 2px solid #000; border-bottom: 1.5px solid #000; }",
    "tbody tr:last-child { border-bottom: 2px solid #000; }",
    "th, td { padding: 6px 10px; text-align: left; font-size: 13px; }",
    "</style>",
    "</head>",
    "<body>"
  )

  if (!is.null(title)) {
    lines <- c(lines, sprintf("<h2>%s</h2>", html_escape(title)))
  }

  lines <- c(
    lines,
    "<table>",
    sprintf("<thead><tr>%s</tr></thead>", header_html),
    "<tbody>",
    body_html,
    "</tbody>",
    "</table>",
    "</body>",
    "</html>"
  )

  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}
