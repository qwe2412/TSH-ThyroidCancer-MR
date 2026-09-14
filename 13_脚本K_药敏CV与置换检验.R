# =====================================================================
# 脚本 K · 药敏预测交叉验证性能 + 置换检验（2026-09-14）
# ------------------------------------------------------------
# 目的：回应 DeepSeek 二轮审稿
#  ① "未报告 oncoPredict 交叉验证性能（R²、RMSE）"：
#     对 GDSC2 训练集代表性药物做 5 折交叉验证（ridge 回归，log IC50 尺度），
#     报告 R² / RMSE / Pearson cor。
#  ② "药物间 IC50 高度相关，独立 FDR 会高估显著药物数"：
#     组标签置换检验（B=500）估计 FDR<0.05 药物数的零分布，
#     与实际 176 对比，评估显著数是否远超随机。
# 前置：D:/gwas/tcga_thca.RData（GDSC2_Expr/GDSC2_Res）
#       D:/gwas/drug_sensitivity_res.RData（pred_ic50/drug_res/grp）
# 运行：Rscript 本文件（约 10-20 分钟）
# =====================================================================
suppressMessages(library(glmnet))

env <- new.env(); load("D:/gwas/all_workspace.RData", envir = env)
GDSC2_Expr <- env$GDSC2_Expr; GDSC2_Res <- env$GDSC2_Res
env2 <- new.env(); load("D:/gwas/drug_sensitivity_res.RData", envir = env2)
grp <- env2$grp; drug_res <- env2$drug_res

cat("GDSC2_Expr:", dim(GDSC2_Expr), "GDSC2_Res:", dim(GDSC2_Res), "\n")

# ============ 1. 代表药物 5 折 CV（ridge，log IC50）============
rep_drugs <- c("5-Fluorouracil_1073", "Paclitaxel_1080", "Docetaxel_1007",
               "Cisplatin_1005", "GSK591_2110", "VE821_2111", "Daporinad_1248",
               "Selumetinib_1736", "Camptothecin_1003", "Vinblastine_1004",
               "Sorafenib_1085", "Nilotinib_1013")
rep_drugs <- rep_drugs[rep_drugs %in% colnames(GDSC2_Res)]

# 基因交集 + 低变异过滤（近似 oncoPredict homogenizeData 流程）
E <- t(GDSC2_Expr)                      # 805 × 17419（细胞系 × 基因）
v <- apply(E, 2, var, na.rm = TRUE)
E <- E[, v > quantile(v, 0.2, na.rm = TRUE)]   # 移除低变异 20%
E <- scale(E)                            # 标准化（standardize 近似）
cat("CV 用基因数：", ncol(E), "\n")

cv_out <- list()
cat("\n========== 1. 代表药物 5 折交叉验证（ridge，log IC50 尺度）==========\n")
set.seed(20260914)
for (dr in rep_drugs) {
  y <- GDSC2_Res[, dr]
  ok <- !is.na(y)
  Xd <- E[ok, ]; yd <- y[ok]
  folds <- sample(rep(1:5, length.out = length(yd)))
  pred <- rep(NA, length(yd))
  for (f in 1:5) {
    tr <- folds != f
    fit <- tryCatch(glmnet(Xd[tr, ], yd[tr], alpha = 0, lambda = 0.1), error = function(e) NULL)
    if (is.null(fit)) next
    pred[folds == f] <- as.numeric(predict(fit, Xd[!tr, ]))
  }
  obs <- pred[!is.na(pred)]; tru <- yd[!is.na(pred)]
  r2 <- 1 - sum((tru - obs)^2) / sum((tru - mean(tru))^2)
  rmse <- sqrt(mean((tru - obs)^2))
  pc <- cor(tru, obs, method = "pearson")
  sc <- cor(tru, obs, method = "spearman")
  cv_out[[dr]] <- c(R2 = r2, RMSE = rmse, Pearson = pc, Spearman = sc, N = length(tru))
  cat(sprintf("  %-22s R²=%.4f RMSE=%.4f Pearson=%.4f Spearman=%.4f (N=%d)\n",
              dr, r2, rmse, pc, sc, length(tru)))
}
cv_df <- do.call(rbind, cv_out)
cat("\nCV 汇总：R² 中位 =", round(median(cv_df[, "R2"]), 4),
    "，范围 =", paste(round(range(cv_df[, "R2"]), 4), collapse = " ~ "),
    "；Pearson 中位 =", round(median(cv_df[, "Pearson"]), 4),
    "；RMSE 中位 =", round(median(cv_df[, "RMSE"]), 4), "\n")
write.csv(as.data.frame(cv_df), "D:/gwas/drug_cv_20260914.csv")

# ============ 2. 组标签置换检验（B=500）============
cat("\n========== 2. 组标签置换检验（B=500）==========\n")
n <- length(grp); B <- 500
n_sig_null <- numeric(B); n_minp_null <- numeric(B)
pk <- seq_len(n)
# 预提取 198 药物的组向量？直接用 drug_res 结构对 pred_ic50 置换
load("D:/gwas/drug_sensitivity_res.RData")   # 重新加载拿 pred_ic50
set.seed(20260914)
for (b in 1:B) {
  g2 <- sample(grp)
  hi <- g2 == "High"; lo <- g2 == "Low"
  pv <- sapply(seq_len(ncol(pred_ic50)), function(j)
    wilcox.test(pred_ic50[hi, j], pred_ic50[lo, j])$p.value)
  n_sig_null[b] <- sum(p.adjust(pv, method = "BH") < 0.05, na.rm = TRUE)
  n_minp_null[b] <- min(pv, na.rm = TRUE)
  if (b %% 100 == 0) cat("  置换", b, "done\n")
}
cat("零分布中 FDR<0.05 药物数：均值 =", round(mean(n_sig_null), 2),
    "；95% 分位 =", quantile(n_sig_null, 0.95),
    "；最大值 =", max(n_sig_null), "\n")
cat("实际显著药物数 = 176；零分布 95% 分位 =", quantile(n_sig_null, 0.95),
    "→ 实际是否远超随机：", 176 > quantile(n_sig_null, 0.95), "\n")
cat("零分布最小 P 值：中位 =", round(median(n_minp_null), 4),
    "；5% 分位 =", round(quantile(n_minp_null, 0.05), 6), "\n")

save(cv_df, n_sig_null, n_minp_null, file = "D:/gwas/drug_cv_perm_20260914.RData")
write.csv(data.frame(n_sig_null = n_sig_null, n_minp_null = n_minp_null),
          "D:/gwas/drug_perm_null_20260914.csv", row.names = FALSE)
cat("\n已保存 D:/gwas/drug_cv_perm_20260914.RData\n")
cat("========== 脚本 K 完成 ==========\n")
