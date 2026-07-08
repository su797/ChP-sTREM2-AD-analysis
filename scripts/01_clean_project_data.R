source(file.path(getwd(), "00_setup.R"))

raw_data <- read_project_data(raw_data_path)
raw_data_all <- read_project_data(raw_data_all_path)

supplementary_name_map <- c(
  "choroid-plexus_SUM" = "ChP_SUM",
  "Lateral-Ventricle" = "LV_SUM"
)
for (source_name in names(supplementary_name_map)) {
  target_name <- supplementary_name_map[[source_name]]
  if (source_name %in% names(raw_data_all) && !(target_name %in% names(raw_data_all))) {
    names(raw_data_all)[names(raw_data_all) == source_name] <- target_name
  }
}

merge_keys <- c("S_ID", "PTID")
supplementary_vars <- setdiff(project_config$variables$keep_vars, names(raw_data))
supplementary_data <- raw_data_all[, unique(c(merge_keys, supplementary_vars)), drop = FALSE]
selected_data <- merge(
  raw_data[, intersect(names(raw_data), setdiff(project_config$variables$keep_vars, supplementary_vars)), drop = FALSE],
  supplementary_data,
  by = merge_keys,
  all.x = TRUE,
  sort = FALSE
)
selected_data <- selected_data[, project_config$variables$keep_vars, drop = FALSE]

numeric_conversion <- convert_selected_to_numeric(
  data = selected_data,
  numeric_vars = setdiff(project_config$variables$keep_vars, c("PTID"))
)
selected_data <- numeric_conversion$data

selected_data$S_DX_label <- recode_from_map(selected_data$S_DX, project_config$labels$S_DX)
selected_data$S_PTGENDER_label <- recode_from_map(selected_data$S_PTGENDER, project_config$labels$S_PTGENDER)
selected_data$APOE401_label <- recode_from_map(selected_data$APOE401, project_config$labels$APOE401)

variable_dictionary <- data.frame(
  variable = c(
    "S_ID",
    "RID",
    "PTID",
    "S_DX",
    "S_DX_label",
    "S_PTGENDER",
    "S_PTGENDER_label",
    "S_AGE",
    "PTEDUCAT",
    "APOE401",
    "APOE401_label",
    "ChPICV",
    "RChPICV",
    "ChP_SUM",
    "MSD_STREM2CORRECTED",
    "EstimatedTotalIntraCranialVol",
    "S_ABETA",
    "TAU",
    "PTAU",
    "Hippocampus_SUM",
    "Amygdala_SUM",
    "LV_SUM",
    "WholeBrain",
    "MMSE",
    "MOCA",
    "mPACCdigit",
    "mPACCtrailsB",
    "ADAS13",
    "CDRSB"
  ),
  role = c(
    "ID",
    "ID",
    "ID",
    "Group",
    "DerivedLabel",
    "Covariate",
    "DerivedLabel",
    "Covariate",
    "Covariate",
    "Covariate",
    "DerivedLabel",
    "Primary",
    "Primary",
    "Primary",
    "Primary",
    "Primary",
    "Biomarker",
    "Biomarker",
    "Biomarker",
    "StructureControl",
    "StructureControl",
    "StructureControl",
    "StructureControl",
    "Cognition",
    "Cognition",
    "Cognition",
    "Cognition",
    "Cognition",
    "Cognition"
  ),
  stringsAsFactors = FALSE
)

merge_audit <- data.frame(
  source_rows_data = nrow(raw_data),
  source_rows_data_all = nrow(raw_data_all),
  merged_rows = nrow(selected_data),
  matched_on_S_ID_PTID = sum(!is.na(selected_data$RID)),
  chpicv_available = sum(!is.na(selected_data$ChPICV)),
  chp_sum_available = sum(!is.na(selected_data$ChP_SUM)),
  stringsAsFactors = FALSE
)
merge_audit_path <- file.path(result_summary_dir, "data_merge_audit.csv")

write_csv_utf8(selected_data, selected_data_path, row.names = FALSE)
latest_outputs <- character(0)
if (isTRUE(project_config$results$write_latest_clean_copy)) {
  write_csv_utf8(selected_data, selected_data_latest_path, row.names = FALSE)
  latest_outputs <- selected_data_latest_path
}
write_csv_utf8(numeric_conversion$conversion_log, conversion_log_path, row.names = FALSE)
write_csv_utf8(variable_dictionary, variable_dictionary_path, row.names = FALSE)
write_csv_utf8(merge_audit, merge_audit_path, row.names = FALSE)

append_analysis_log(
  project_root = project_root,
  analysis_name = "01_clean_project_data",
  output_files = c(selected_data_path, latest_outputs, conversion_log_path, variable_dictionary_path, merge_audit_path),
  note = "Merged Data.csv with Data_all.csv by S_ID and PTID, preserved the original analysis variables, added absolute choroid plexus and structure-control measures, converted numeric fields, and archived the cleaned dataset into the current result version.",
  summary_dir = result_summary_dir
)
