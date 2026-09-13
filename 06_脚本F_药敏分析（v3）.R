# =====================================================================
# 脚本 F v3 · 药物敏感性分析（修复版）
# 对应手册第 12 章：基于 GDSC2 的 TCGA-THY 药物 IC50 预测
#                     → 高危/低危组间药物敏感性差异
# ------------------------------------------------------------
# 【v3 修复（2026-09-13，严格审查后，依据官方 vignette）】
#  ① calcPhenotype() 1.3.1 参数签名（官方 vignette 确认）：
#     trainingExprData / trainingPtype / testExprData / batchCorrect /
#     powerTransformPhenotype / removeLowVaryingGenes /
#     removeLowVaringGenesFrom / minNumSamples / selection / printOutput /
#     pcr / report_pc / cc / rsq / percent
#     → 没有 pValCutoff、没有 seed（v2 传这两个参数导致
#       "参数没有用(pValCutoff = 0.05, seed = 123)"）；seed 用 set.seed() 设定。
#  ② GDSC IC50 是 log 尺度：vignette 明确 "The GDSC IC50 values are already
#     log-transformed. Convert them back" → trainingPtype <- exp(GDSC2_Res)
#  ③ batchCorrect 选择：训练=GDSC 微阵列、测试=TCGA RNA-seq →
#     vignette 推荐 "standardize"（z-score 标准化），非 "eb"。
#  ④ 其余 v2 修复保留：dat_final→dat、内嵌 RData、分组对齐、历史对账。
# 【前置】oncoPredict 1.3.1 + sva 3.60.0 已装（v2 已确认）
# 【耗时】预测约 10-30 分钟
# 【运行方式】Ctrl+A 全选 → Run → 把全部输出发回
# =====================================================================

# ================= 0. 包与数据加载 =================
if (!requireNamespace("oncoPredict", quietly = TRUE)) {
  install.packages("oncoPredict", repos = "https://cloud.r-project.org")
}
library(oncoPredict)
library(sva)
cat("✓ oncoPredict 版本：", as.character(packageVersion("oncoPredict")),
    "；sva 版本：", as.character(packageVersion("sva")), "\n")

# ================= 1. 数据加载 =================
stopifnot("tcga_thca.RData 不存在" = file.exists("D:/gwas/tcga_thca.RData"))
stopifnot("risk_model.RData 不存在" = file.exists("D:/gwas/risk_model.RData"))
load("D:/gwas/tcga_thca.RData")      # GDSC2_Expr / GDSC2_Res / ic50_pred / expr_mat_use ...
load("D:/gwas/risk_model.RData")     # dat（409×30，含 riskGroup）
stopifnot("RData 中无 GDSC2_Expr" = exists("GDSC2_Expr"))
stopifnot("RData 中无 GDSC2_Res" = exists("GDSC2_Res"))
stopifnot("RData 中无 dat（应为 dat，不是 dat_final）" = exists("dat"))
cat("✓ GDSC2_Expr：", dim(GDSC2_Expr), "；GDSC2_Res：", dim(GDSC2_Res), "；dat：", dim(dat), "\n")

# ================= 2. 分组对齐（以 dat 顺序为准）=================
cat("\n========== 1. 样本对齐与分组 ==========\n")
if (exists("expr_mat_use")) {
  expr_for_test <- expr_mat_use
} else {
  expr_for_test <- get("expr_mat", envir = .GlobalEnv)
}
stopifnot("表达矩阵行名非患者 ID" = nchar(rownames(expr_for_test)[1]) == 12)
common <- intersect(rownames(expr_for_test), as.character(dat$patientId))
cat("共同患者：", length(common), "/", nrow(dat), "\n")
stopifnot(length(common) >= 100)
idx <- match(as.character(dat$patientId), rownames(expr_for_test))
idx <- idx[!is.na(idx)]
expr_test_aligned <- expr_for_test[idx, ]
grp <- dat$riskGroup[match(rownames(expr_test_aligned), as.character(dat$patientId))]
stopifnot(all(!is.na(grp)))
cat("测试表达矩阵（患者 × 基因）：", dim(expr_test_aligned), "；分组分布：\n")
print(table(grp, useNA = "ifany"))

# ================= 3. GDSC2 数据预处理（v3：log IC50 → exp 转回）=================
cat("\n========== 2. GDSC2 数据预处理 ==========\n")
cat("GDSC2_Res 数值范围（log 尺度）：",
    paste(round(range(GDSC2_Res, na.rm = TRUE), 4), collapse = " ~ "),
    "；NA 数：", sum(is.na(GDSC2_Res)), "\n")
GDSC2_Res_exp <- exp(as.matrix(GDSC2_Res))      # IC50 转回原始尺度
cat("exp() 后范围：", paste(round(range(GDSC2_Res_exp, na.rm = TRUE), 4), collapse = " ~ "), "\n")
if (exists("ic50_pred")) {
  cat("历史 ic50_pred 维度：", dim(ic50_pred),
      "；范围：", paste(round(range(ic50_pred, na.rm = TRUE), 4), collapse = " ~ "), "\n")
}

# ================= 4. 药物敏感性预测（oncoPredict，v3 修正参数）=================
cat("\n========== 3. oncoPredict 预测（约 10-30 分钟，请耐心等待）==========\n")
test_expr_t <- t(expr_test_aligned)                 # 基因 × 患者
cat("训练表达矩阵（基因×细胞系）：", dim(GDSC2_Expr), "\n")
cat("训练药物表（细胞系×药物，已 exp）：", dim(GDSC2_Res_exp), "\n")
cat("测试表达矩阵（基因×患者）：", dim(test_expr_t), "\n")

set.seed(123)                                       # v3：seed 在外部设定

pred_ic50 <- calcPhenotype(
  trainingExprData = GDSC2_Expr,
  trainingPtype    = GDSC2_Res_exp,                # v3：exp 转回
  testExprData     = test_expr_t,
  batchCorrect     = "standardize",                # v3：microarray→RNA-seq 推荐
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  removeLowVaringGenesFrom = "homogenizeData",
  minNumSamples = 10,
  selection = 1,
  printOutput = TRUE,
  pcr = FALSE,
  report_pc = FALSE,
  cc = FALSE,
  rsq = FALSE,
  percent = 80
)
cat("✓ 预测 IC50 矩阵维度（患者 × 药物）：", dim(pred_ic50), "\n")

# ================= 5. 与历史 ic50_pred 对账 =================
cat("\n========== 4. 与历史 ic50_pred 对账 ==========\n")
if (exists("ic50_pred") && is.matrix(ic50_pred)) {
  cat("历史 ic50_pred 维度：", dim(ic50_pred), "\n")
  cat("历史 ic50_pred 行名前 3：", paste(head(rownames(ic50_pred), 3), collapse = ", "), "\n")
  cat("历史 ic50_pred 列名前 5：", paste(head(colnames(ic50_pred), 5), collapse = ", "), "\n")
  common_drugs <- intersect(colnames(pred_ic50), colnames(ic50_pred))
  cat("共同药物数：", length(common_drugs), "/", ncol(ic50_pred), "\n")
  common_pats <- intersect(rownames(pred_ic50), rownames(ic50_pred))
  cat("共同患者数：", length(common_pats), "/", nrow(ic50_pred), "\n")
  if (length(common_drugs) > 10 && length(common_pats) > 100) {
    cors <- sapply(common_drugs, function(d) {
      cor(pred_ic50[common_pats, d], ic50_pred[common_pats, d],
          use = "complete.obs")
    })
    cat("逐药物 cor 汇总：中位数 =", round(median(cors, na.rm = TRUE), 4),
        "；范围 =", paste(round(range(cors, na.rm = TRUE), 4), collapse = " ~ "),
        "；NA 数 =", sum(is.na(cors)), "\n")
    cat("cor > 0.9 的药物数：", sum(cors > 0.9, na.rm = TRUE), "/", length(common_drugs), "\n")
    cat("cor > 0.5 的药物数：", sum(cors > 0.5, na.rm = TRUE), "/", length(common_drugs), "\n")
    cat("cor > 0.3 的药物数：", sum(cors > 0.3, na.rm = TRUE), "/", length(common_drugs), "\n")
  }
} else {
  cat("⚠ 当前会话无 ic50_pred（历史预测），跳过对账。\n")
}

# ================= 6. High vs Low 组间差异 =================
cat("\n========== 5. 高危 vs 低危组药物敏感性差异（Wilcoxon）==========\n")
hi <- grp == "High"; lo <- grp == "Low"
cat("High：", sum(hi), "；Low：", sum(lo), "\n")
drug_res <- data.frame(drug = colnames(pred_ic50),
                       mean_High = NA_real_, mean_Low = NA_real_,
                       P = NA_real_, P_BH = NA_real_,
                       stringsAsFactors = FALSE)
for (j in seq_len(ncol(pred_ic50))) {
  a <- pred_ic50[hi, j]; b <- pred_ic50[lo, j]
  w <- tryCatch(wilcox.test(a, b), error = function(e) NULL)
  drug_res$mean_High[j] <- mean(a, na.rm = TRUE)
  drug_res$mean_Low[j]  <- mean(b, na.rm = TRUE)
  drug_res$P[j] <- if (!is.null(w)) w$p.value else NA
}
drug_res$P_BH <- p.adjust(drug_res$P, method = "BH")
cat("P_BH < 0.05 的药物数：", sum(drug_res$P_BH < 0.05, na.rm = TRUE), "/", nrow(drug_res), "\n")
sig_drugs <- drug_res[order(drug_res$P, na.last = TRUE)[1:min(20, nrow(drug_res))], ]
cat("\nP 值最小的前 20 个药物：\n")
print(sig_drugs, row.names = FALSE, digits = 5)

# ================= 7. 保存（备份保护）=================
cat("\n========== 6. 保存 drug_sensitivity_res ==========\n")
if (file.exists("D:/gwas/drug_sensitivity_res.RData")) {
  if (!file.exists("D:/gwas/drug_sensitivity_res_backup.RData")) {
    file.copy("D:/gwas/drug_sensitivity_res.RData",
              "D:/gwas/drug_sensitivity_res_backup.RData")
    cat("✓ 历史 drug_sensitivity_res.RData 已备份\n")
  }
}
save(pred_ic50, drug_res, grp, file = "D:/gwas/drug_sensitivity_res.RData")
write.csv(drug_res, "D:/gwas/drug_sensitivity_res.csv", row.names = FALSE)
cat("✓ 已保存 drug_sensitivity_res.RData / .csv\n")

cat("\n========== 脚本 F v3 完成：请把全部输出发回 ==========\n")
