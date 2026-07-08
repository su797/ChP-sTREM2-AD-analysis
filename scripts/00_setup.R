# Shared setup for the public analysis scripts.
#
# The original project had a larger reporting layer. For the GitHub release we
# keep this file deliberately plain: set paths, load helper functions, create a
# timestamped results folder, and leave submission-document generation out of the picture.

script_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(script_root, ".."), winslash = "/", mustWork = TRUE)
common_root <- file.path(repo_root, "R")

# Backward-compatible names used by the original analysis modules.
statistics_root <- script_root
project_root <- repo_root

source(file.path(common_root, "utils_io.R"))
load_common_functions(common_root)
source(file.path(script_root, "project_config.R"))

settings_path <- file.path(repo_root, "config", "analysis_settings.json")
runtime_settings <- read_json_config(settings_path, default = list())
runtime_settings$language <- runtime_settings$language %||% "en"
runtime_settings$group_order <- runtime_settings$group_order %||% c("CN", "MCI", "AD")
runtime_settings$descriptive <- runtime_settings$descriptive %||% list(
  digits = 2,
  continuous_overall = c("S_AGE"),
  categorical_overall = c("S_PTGENDER_label"),
  continuous_by_group = c("S_AGE"),
  categorical_by_group = c("S_PTGENDER_label")
)
runtime_settings$display_labels <- runtime_settings$display_labels %||% list(
  variables = list(),
  levels = list()
)
runtime_settings$report <- runtime_settings$report %||% list(digits = 4)
report_digits <- runtime_settings$report$digits %||% project_config$output$digits %||% 4
report_language <- runtime_settings$language
locale_bundle <- list()
project_report_title <- project_config$project_name

if (!is.null(runtime_settings$project_config_override)) {
  project_config <- merge_lists_recursive(project_config, runtime_settings$project_config_override)
}

set.seed(project_config$seed)

# Reuse the same run id when the pipeline sources this setup file repeatedly.
current_result_version <- Sys.getenv(
  "CHP_STREM2_RUN_ID",
  unset = Sys.getenv("RESEARCHR_RESULT_VERSION", unset = "")
)
if (!nzchar(current_result_version)) {
  current_result_version <- format(Sys.time(), project_config$results$dir_time_format %||% "%Y%m%d_%H%M%S")
}
Sys.setenv(CHP_STREM2_RUN_ID = current_result_version)
Sys.setenv(RESEARCHR_RESULT_VERSION = current_result_version)

ensure_release_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

result_base_dir <- ensure_release_dir(file.path(repo_root, "results"))
result_run_dir <- ensure_release_dir(file.path(result_base_dir, current_result_version))
result_summary_dir <- ensure_release_dir(file.path(result_run_dir, "summary"))
result_tables_dir <- ensure_release_dir(file.path(result_run_dir, "tables"))
result_figures_dir <- ensure_release_dir(file.path(result_run_dir, "figures"))
result_report_dir <- ensure_release_dir(file.path(result_run_dir, "diagnostics"))
result_data_clean_dir <- ensure_release_dir(file.path(result_run_dir, project_config$results$data_subdir %||% "data_clean"))
latest_clean_dir <- ensure_release_dir(file.path(repo_root, "data", "clean"))

result_version <- current_result_version
raw_data_path <- file.path(repo_root, "data", "raw", "Data.csv")
raw_data_all_path <- file.path(repo_root, "data", "raw", "Data_all.csv")
manual_qc_data_path <- file.path(repo_root, "data", "raw", "manual_chp_qc_integrated_80.csv")

selected_data_path <- file.path(result_data_clean_dir, "ChpSTREM2AD_selected_dataset.csv")
analysis_data_path <- file.path(result_data_clean_dir, "ChpSTREM2AD_analysis_dataset.csv")
selected_data_latest_path <- file.path(latest_clean_dir, "ChpSTREM2AD_selected_dataset.csv")
analysis_data_latest_path <- file.path(latest_clean_dir, "ChpSTREM2AD_analysis_dataset.csv")

normality_table_path <- file.path(result_summary_dir, "normality_screening.csv")
conversion_log_path <- file.path(result_summary_dir, "numeric_conversion_log.csv")
transformation_table_path <- file.path(result_summary_dir, "transformation_plan.csv")
cognition_map_path <- file.path(result_summary_dir, "cognition_component_map.csv")
variable_dictionary_path <- file.path(result_summary_dir, "selected_variable_dictionary.csv")
result_run_metadata_path <- file.path(result_run_dir, "run_metadata.csv")

writeLines(
  c(
    paste0("run_id=", result_version),
    paste0("run_dir=", result_run_dir)
  ),
  con = file.path(result_base_dir, "LATEST.txt"),
  useBytes = TRUE
)

write_csv_utf8(
  data.frame(
    project_name = project_config$project_name,
    result_version = result_version,
    generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    raw_data_path = raw_data_path,
    raw_data_all_path = raw_data_all_path,
    stringsAsFactors = FALSE
  ),
  result_run_metadata_path,
  row.names = FALSE
)
