# 这里读取本地分析数据；投稿附件只保留汇总输出。
# 実行時は第1引数に作業先、第2引数にデータCSVを指定する。省略時のパスは元の解析環境用。
options(stringsAsFactors = FALSE, digits = 17)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1, 2))
root <- args[1]
stopifnot(dir.exists(root))
out <- file.path(root, "analysis_outputs")
source_csv <- if (length(args) == 2) args[2] else "I:/researchR/project/ChpSTREM2AD/result/20260419_155221/data_clean/ChpSTREM2AD_analysis_dataset.csv"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
d <- read.csv(source_csv, check.names = FALSE)
# 先核对人数和 PTID，避免误用重复记录或不同版本的数据。
stopifnot(nrow(d) == 735, !anyDuplicated(d$PTID))
required <- c("S_AGE", "S_PTGENDER", "PTEDUCAT", "APOE401", "S_DX_label", "ChPICV", "MSD_STREM2CORRECTED", "TAU", "PTAU", "LV_SUM", "Hippocampus_SUM", "WholeBrain", "EstimatedTotalIntraCranialVol")
stopifnot(all(required %in% names(d)))
positive <- c("ChPICV", "MSD_STREM2CORRECTED", "TAU", "PTAU", "LV_SUM", "Hippocampus_SUM", "WholeBrain", "EstimatedTotalIntraCranialVol", "PTEDUCAT")
stopifnot(all(vapply(d[positive], function(x) all(is.finite(x[!is.na(x)]) & x[!is.na(x)] > 0), logical(1))))
d$age <- d$S_AGE
d$sex <- factor(d$S_PTGENDER)
# 主增量模型用教育年数原值，log 版本留给敏感性分析。
d$education <- d$PTEDUCAT
d$log_education <- log(d$PTEDUCAT)
d$apoe <- factor(d$APOE401)
d$diagnosis <- factor(d$S_DX_label, levels = c("CN", "MCI", "AD"))
d$chp <- log(d$ChPICV)
d$lv <- log(d$LV_SUM / d$EstimatedTotalIntraCranialVol)
d$hipp <- log(d$Hippocampus_SUM / d$EstimatedTotalIntraCranialVol)
d$brain <- log(d$WholeBrain / d$EstimatedTotalIntraCranialVol)
d$log_icv <- log(d$EstimatedTotalIntraCranialVol)
d$abs_lv <- log(d$LV_SUM)
d$abs_hipp <- log(d$Hippocampus_SUM)
d$abs_brain <- log(d$WholeBrain)
# Base と Extended は同じ701人で比較する。モデルごとに対象者を変えない。
cc <- complete.cases(d[required])
stopifnot(sum(cc) == 701)
a <- d[cc, ]
outcomes <- c(sTREM2 = "MSD_STREM2CORRECTED", Tau = "TAU", `P-Tau` = "PTAU")
specs <- list(
  primary = c("age", "sex", "education", "apoe", "lv", "hipp", "brain"),
  diagnosis_adjusted = c("age", "sex", "education", "apoe", "lv", "hipp", "brain", "diagnosis"),
  absolute_volumes_ICV = c("age", "sex", "education", "apoe", "abs_lv", "abs_hipp", "abs_brain", "log_icv"),
  log_education = c("age", "sex", "log_education", "apoe", "lv", "hipp", "brain")
)
calc_vif <- function(fit) {
  x <- model.matrix(fit)[, -1, drop = FALSE]
  sapply(seq_len(ncol(x)), function(j) 1 / (1 - summary(lm(x[, j] ~ x[, -j, drop = FALSE]))$r.squared))
}
# HC3 は標準誤差の感度確認用。下の通常の入れ子モデルF検定とは別に保存する。
hc3 <- function(fit, term) {
  x <- model.matrix(fit)
  bread <- solve(crossprod(x))
  w <- residuals(fit)^2 / (1 - hatvalues(fit))^2
  v <- bread %*% crossprod(x, x * w) %*% bread
  se <- sqrt(v[term, term]); b <- coef(fit)[term]; df <- df.residual(fit)
  c(se = unname(se), p = unname(2 * pt(-abs(b / se), df)), low = unname(b - qt(.975, df) * se), high = unname(b + qt(.975, df) * se))
}
rows <- list(); coefficients <- list(); formulas <- list(); nr <- 0
pdf(file.path(out, "incremental_model_diagnostics.pdf"), width = 10, height = 7)
for (spec in names(specs)) {
  for (label in names(outcomes)) {
    a$y <- log(a[[outcomes[[label]]]])
    f0 <- reformulate(specs[[spec]], response = "y")
    f1 <- reformulate(c(specs[[spec]], "chp"), response = "y")
    m0 <- lm(f0, data = a, na.action = na.fail)
    m1 <- lm(f1, data = a, na.action = na.fail)
    stopifnot(identical(rownames(model.frame(m0)), rownames(model.frame(m1))))
    stopifnot(m0$rank == ncol(model.matrix(m0)), m1$rank == ncol(model.matrix(m1)))
    cmp <- anova(m0, m1)
    s0 <- summary(m0); s1 <- summary(m1)
    delta <- s1$r.squared - s0$r.squared
# partial R2 的分母是基础模型的残差平方和，不是总平方和。
    partial <- (deviance(m0) - deviance(m1)) / deviance(m0)
    stopifnot(abs(partial - delta / (1 - s0$r.squared)) < 1e-10, delta >= -1e-12)
    ct <- coef(s1)["chp", ]; ci <- confint(m1, "chp")
    stopifnot(abs(cmp$F[2] - ct["t value"]^2) < 1e-8)
    stopifnot(abs(cmp$`Pr(>F)`[2] - ct["Pr(>|t|)"]) < 1e-10)
# 标准化使用当前分析样本中、对数尺度上的标准差。
    scale <- sd(a$chp) / sd(a$y)
    robust <- hc3(m1, "chp")
    nr <- nr + 1
    rows[[nr]] <- data.frame(specification = spec, biomarker = label, n = nobs(m1),
      base_r2 = s0$r.squared, extended_r2 = s1$r.squared,
      base_adj_r2 = s0$adj.r.squared, extended_adj_r2 = s1$adj.r.squared,
      delta_r2 = delta, delta_adj_r2 = s1$adj.r.squared - s0$adj.r.squared,
      partial_r2 = partial, F = cmp$F[2], df1 = cmp$Df[2], df2 = df.residual(m1),
      comparison_p = cmp$`Pr(>F)`[2], beta = ct["Estimate"], beta_low = ci[1], beta_high = ci[2],
      standardized_beta = ct["Estimate"] * scale, standardized_low = ci[1] * scale, standardized_high = ci[2] * scale,
      hc3_se = robust["se"], hc3_p = robust["p"], hc3_low = robust["low"], hc3_high = robust["high"],
      max_vif = max(calc_vif(m1)), max_cooks_d = max(cooks.distance(m1)),
      count_cooks_above_4_n = sum(cooks.distance(m1) > 4 / nobs(m1)),
      max_abs_studentized_residual = max(abs(rstudent(m1))),
      bp_p = pchisq(nobs(m1) * summary(lm(residuals(m1)^2 ~ model.matrix(m1)[, -1]))$r.squared,
                    df = m1$rank - 1, lower.tail = FALSE), row.names = NULL)
    for (model_name in c("base", "extended")) {
      m <- if (model_name == "base") m0 else m1
      tab <- coef(summary(m)); cit <- confint(m)
      coefficients[[length(coefficients) + 1]] <- data.frame(specification = spec, biomarker = label, model = model_name,
        term = rownames(tab), beta = tab[, 1], se = tab[, 2], t = tab[, 3], p = tab[, 4], low = cit[, 1], high = cit[, 2], row.names = NULL)
    }
    formulas[[nr]] <- data.frame(specification = spec, biomarker = label,
      base_formula = paste(deparse(f0), collapse = " "), extended_formula = paste(deparse(f1), collapse = " "))
    if (spec == "primary") {
      par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
      plot(m1, which = c(1, 2, 3, 5), id.n = 0)
      mtext(paste("Primary incremental model:", label, "(n = 701)"), outer = TRUE)
    }
  }
}
dev.off()
results <- do.call(rbind, rows)
results$holm_p <- NA_real_
main <- results$specification == "primary"
# Holm補正は主解析の3つのバイオマーカー間で行う。感度分析のP値はここには含めない。
results$holm_p[main] <- p.adjust(results$comparison_p[main], method = "holm")
write.csv(results, file.path(out, "incremental_results.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, coefficients), file.path(out, "all_model_coefficients.csv"), row.names = FALSE)
write.csv(do.call(rbind, formulas), file.path(out, "model_formulas.csv"), row.names = FALSE)

# Verify the standardized sTREM2 coefficient in the full-cohort association model.
# 这个核对模型以 ChP/ICV 为结局；不要和上面的生物标志物结局模型混淆。
legacy <- lm(log(ChPICV) ~ log(MSD_STREM2CORRECTED) + factor(S_PTGENDER) + S_AGE + log(PTEDUCAT) + factor(APOE401), data = d)
legacy_beta <- coef(legacy)[2] * sd(log(d$MSD_STREM2CORRECTED)) / sd(log(d$ChPICV))
stopifnot(abs(legacy_beta - (-0.172289629230331)) < 1e-10)
write.csv(data.frame(check = "legacy_primary_sTREM2_standardized_beta", computed = unname(legacy_beta), archived = -0.172289629230331, absolute_difference = abs(legacy_beta + 0.172289629230331)), file.path(out, "legacy_reproduction_check.csv"), row.names = FALSE)

# Aggregate inclusion summaries only: no participant identifiers in outputs.
counts <- as.data.frame(table(included = cc, diagnosis = d$diagnosis))
write.csv(counts, file.path(out, "inclusion_counts.csv"), row.names = FALSE)
continuous <- c("S_AGE", "PTEDUCAT", "ChPICV", "MSD_STREM2CORRECTED", "TAU", "PTAU", "LV_SUM", "Hippocampus_SUM", "EstimatedTotalIntraCranialVol")
summaries <- lapply(continuous, function(v) {
  x <- d[[v]][cc]; y <- d[[v]][!cc]
  data.frame(variable = v, included_n = length(x), included_mean = mean(x), included_sd = sd(x), excluded_n = length(y), excluded_mean = mean(y), excluded_sd = sd(y), smd = (mean(x) - mean(y)) / sqrt((var(x) + var(y)) / 2))
})
write.csv(do.call(rbind, summaries), file.path(out, "inclusion_continuous_summary.csv"), row.names = FALSE)
cats <- list()
for (v in c("S_PTGENDER", "APOE401", "S_DX_label")) {
  for (level in unique(d[[v]])) {
    p1 <- mean(d[[v]][cc] == level); p0 <- mean(d[[v]][!cc] == level)
    cats[[length(cats) + 1]] <- data.frame(variable = v, level = level, included_n = sum(d[[v]][cc] == level), excluded_n = sum(d[[v]][!cc] == level), included_percent = 100 * p1, excluded_percent = 100 * p0, smd = (p1 - p0) / sqrt((p1 * (1-p1) + p0 * (1-p0))/2))
  }
}
write.csv(do.call(rbind, cats), file.path(out, "inclusion_categorical_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "R_session_info.txt"))
saveRDS(list(results = results, formulas = formulas), file.path(out, "aggregate_analysis_objects.rds"))
print(results[, c("specification", "biomarker", "n", "base_r2", "extended_r2", "delta_r2", "partial_r2", "comparison_p", "holm_p", "standardized_beta", "hc3_p", "max_vif")], row.names = FALSE)
