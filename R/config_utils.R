`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    return(y)
  }
  x
}

read_json_config <- function(path, default = list()) {
  if (!file.exists(path)) {
    return(default)
  }

  json_skip_ws <- function(text, pos) {
    text_len <- nchar(text, type = "chars")
    while (pos <= text_len) {
      ch <- substr(text, pos, pos)
      if (!grepl("[[:space:]]", ch)) {
        break
      }
      pos <- pos + 1
    }
    pos
  }

  json_parse_string <- function(text, pos) {
    text_len <- nchar(text, type = "chars")
    pos <- pos + 1
    out <- character(0)

    while (pos <= text_len) {
      ch <- substr(text, pos, pos)
      if (identical(ch, "\"")) {
        return(list(value = paste(out, collapse = ""), pos = pos + 1))
      }
      if (identical(ch, "\\")) {
        pos <- pos + 1
        esc <- substr(text, pos, pos)
        esc_value <- switch(
          esc,
          "\"" = "\"",
          "\\" = "\\",
          "/" = "/",
          "b" = "\b",
          "f" = "\f",
          "n" = "\n",
          "r" = "\r",
          "t" = "\t",
          esc
        )
        out <- c(out, esc_value)
        pos <- pos + 1
      } else {
        out <- c(out, ch)
        pos <- pos + 1
      }
    }

    stop("Invalid JSON string.", call. = FALSE)
  }

  json_parse_number <- function(text, pos) {
    remaining <- substr(text, pos, nchar(text, type = "chars"))
    match <- regexpr("^-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?", remaining, perl = TRUE)
    if (match[[1]] != 1) {
      stop("Invalid JSON number.", call. = FALSE)
    }
    token <- regmatches(remaining, match)[[1]]
    list(value = as.numeric(token), pos = pos + nchar(token, type = "chars"))
  }

  json_parse_literal <- function(text, pos) {
    remaining <- substr(text, pos, nchar(text, type = "chars"))
    if (startsWith(remaining, "true")) {
      return(list(value = TRUE, pos = pos + 4))
    }
    if (startsWith(remaining, "false")) {
      return(list(value = FALSE, pos = pos + 5))
    }
    if (startsWith(remaining, "null")) {
      return(list(value = NULL, pos = pos + 4))
    }
    stop("Invalid JSON literal.", call. = FALSE)
  }

  json_parse_array <- function(text, pos) {
    items <- list()
    pos <- json_skip_ws(text, pos + 1)
    if (substr(text, pos, pos) == "]") {
      return(list(value = items, pos = pos + 1))
    }

    idx <- 1
    repeat {
      parsed <- json_parse_value(text, pos)
      items[[idx]] <- parsed$value
      idx <- idx + 1
      pos <- json_skip_ws(text, parsed$pos)
      next_char <- substr(text, pos, pos)
      if (next_char == "]") {
        return(list(value = items, pos = pos + 1))
      }
      if (next_char != ",") {
        stop("Invalid JSON array.", call. = FALSE)
      }
      pos <- json_skip_ws(text, pos + 1)
    }
  }

  json_parse_object <- function(text, pos) {
    object <- list()
    pos <- json_skip_ws(text, pos + 1)
    if (substr(text, pos, pos) == "}") {
      return(list(value = object, pos = pos + 1))
    }

    repeat {
      key_parsed <- json_parse_string(text, pos)
      key <- key_parsed$value
      pos <- json_skip_ws(text, key_parsed$pos)
      if (substr(text, pos, pos) != ":") {
        stop("Invalid JSON object.", call. = FALSE)
      }
      pos <- json_skip_ws(text, pos + 1)
      value_parsed <- json_parse_value(text, pos)
      object[[key]] <- value_parsed$value
      pos <- json_skip_ws(text, value_parsed$pos)
      next_char <- substr(text, pos, pos)
      if (next_char == "}") {
        return(list(value = object, pos = pos + 1))
      }
      if (next_char != ",") {
        stop("Invalid JSON object.", call. = FALSE)
      }
      pos <- json_skip_ws(text, pos + 1)
    }
  }

  json_parse_value <- function(text, pos) {
    pos <- json_skip_ws(text, pos)
    first_char <- substr(text, pos, pos)

    if (first_char == "{") {
      return(json_parse_object(text, pos))
    }
    if (first_char == "[") {
      return(json_parse_array(text, pos))
    }
    if (first_char == "\"") {
      return(json_parse_string(text, pos))
    }
    if (grepl("[-0-9]", first_char)) {
      return(json_parse_number(text, pos))
    }
    json_parse_literal(text, pos)
  }

  json_text <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  parsed <- json_parse_value(json_text, 1)$value
  if (is.null(parsed)) {
    return(default)
  }
  parsed
}

merge_lists_recursive <- function(base, override) {
  if (is.null(base) || length(base) == 0) {
    return(override)
  }
  if (is.null(override) || length(override) == 0) {
    return(base)
  }

  if (!is.list(base) || !is.list(override) || is.null(names(override))) {
    return(override)
  }

  merged <- base
  for (name in names(override)) {
    if (!is.null(merged[[name]]) && is.list(merged[[name]]) && is.list(override[[name]]) && !is.null(names(override[[name]]))) {
      merged[[name]] <- merge_lists_recursive(merged[[name]], override[[name]])
    } else {
      merged[[name]] <- override[[name]]
    }
  }
  merged
}

load_runtime_settings <- function(common_defaults_path, project_local_path = NULL) {
  base_settings <- read_json_config(common_defaults_path, default = list())
  project_settings <- read_json_config(project_local_path, default = list())
  merge_lists_recursive(base_settings, project_settings)
}

resolve_language <- function(settings, default = "zh", supported = c("zh", "en", "ja")) {
  language <- settings$language %||% default
  if (!language %in% supported) {
    return(default)
  }
  language
}

resolve_multilingual_value <- function(entry, language, fallback = NULL) {
  if (is.null(entry)) {
    return(fallback)
  }

  if (is.atomic(entry) && length(entry) == 1) {
    return(as.character(entry))
  }

  if (is.list(entry)) {
    if (!is.null(entry[[language]])) {
      return(as.character(entry[[language]]))
    }
    if (!is.null(entry$default)) {
      return(as.character(entry$default))
    }
    if (!is.null(entry$zh)) {
      return(as.character(entry$zh))
    }
    if (!is.null(entry$en)) {
      return(as.character(entry$en))
    }
    if (!is.null(entry$ja)) {
      return(as.character(entry$ja))
    }
  }

  fallback
}

interpolate_text <- function(template, values = list()) {
  out <- template
  if (length(values) == 0) {
    return(out)
  }

  for (name in names(values)) {
    placeholder <- paste0("{", name, "}")
    value <- ifelse(is.null(values[[name]]) || is.na(values[[name]]), "", as.character(values[[name]]))
    out <- gsub(placeholder, value, out, fixed = TRUE)
  }
  out
}

translate_text <- function(locale_bundle, language, key, values = list(), fallback_language = "zh") {
  lang_pack <- locale_bundle[[language]] %||% list()
  template <- lang_pack[[key]]

  if (is.null(template)) {
    template <- locale_bundle[[fallback_language]][[key]] %||% key
  }

  interpolate_text(as.character(template), values = values)
}

resolve_variable_label <- function(variable_name, settings, language, fallback = NULL) {
  fallback <- fallback %||% variable_name
  label_entry <- settings$display_labels$variables[[variable_name]]
  resolve_multilingual_value(label_entry, language, fallback = fallback)
}

resolve_level_label <- function(variable_name, level_value, settings, language, fallback = NULL) {
  fallback <- fallback %||% as.character(level_value)
  level_key <- as.character(level_value)
  label_entry <- settings$display_labels$levels[[variable_name]][[level_key]]
  resolve_multilingual_value(label_entry, language, fallback = fallback)
}
