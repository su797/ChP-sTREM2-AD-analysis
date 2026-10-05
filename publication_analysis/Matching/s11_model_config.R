# S11 使用的协变量和三个模型。ChPICV 为解释变量，生物标志物为结局。
# 変数の変換は別の変換表で指定するため、ここは元の列名のまま。
list(variables = list(group_label_var = "S_DX_label", covariates = c("S_PTGENDER", 
"S_AGE", "PTEDUCAT", "APOE401"), categorical_covariates = c("S_PTGENDER", 
"APOE401")), robustness = list(matched_key_models = list(list(
    name = "sTREM2_on_ChPICV", outcome = "MSD_STREM2CORRECTED", 
    exposure = "ChPICV"), list(name = "PTAU_on_ChPICV", outcome = "PTAU", 
    exposure = "ChPICV"), list(name = "TAU_on_ChPICV", outcome = "TAU", 
    exposure = "ChPICV"))))
