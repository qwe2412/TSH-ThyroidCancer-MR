# =====================================================================
# 脚本 I · Lasso-Cox 与弹性网变量选择（严格审查补跑，2026-09-13）
# ------------------------------------------------------------
# 【背景】审稿指出"Lasso-Cox（lambda.min 与 lambda.1se）及弹性网
#   （alpha=0.5）三条正则化路径下均被保留"缺乏脚本证据。
# 【本次重跑结论（seed=20260913，10 折交叉验证）】
#   Lasso (alpha=1)   lambda.min → 43 变量保留（含 GNG7 b=-0.0145、
#                     NFIA b=-0.0008，均负向）；lambda.1se → 0 变量
#   弹性网 (alpha=0.5) lambda.min = lambda.1se → 0 变量
#   → "三条路径均保留"表述不成立；仅 lambda.min 路径保留 GNG7/NFIA。
# 【数据】lasso_cv.RData 内 dat（409×68：patient/sample + 62 候选基因
#   + 临床协变量 + DFS_time/DFS_event）；缺失值以中位数填补。
# 【运行方式】Ctrl+A 全选 → Run（约 1-2 分钟）
# =====================================================================
suppressMessages(library(glmnet))
suppressMessages(library(survival))

# ============ 1. 数据 ============
stopifnot(file.exists("D:/gwas/lasso_cv.RData"))
env <- new.env(); load("D:/gwas/lasso_cv.RData", envir = env)
dat <- env$dat
time_col <- "DFS_time"; event_col <- "DFS_event"
genes <- setdiff(colnames(dat), c("patient", "sample", time_col, event_col,
                                  "DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS"))
cat("候选基因列数：", length(genes), "\n")
X <- as.matrix(dat[, genes])
cat("X NA 数：", sum(is.na(X)), "→ 以中位数填补\n")
X[is.na(X)] <- median(X, na.rm = TRUE)
y <- Surv(dat[[time_col]], dat[[event_col]])
cat("事件数：", sum(dat[[event_col]]), "/", nrow(dat), "\n")

# ============ 2. Lasso-Cox（alpha=1）============
cat("\n========== 2. Lasso-Cox（alpha=1，10 折 CV，seed=20260913）==========\n")
set.seed(20260913)
cv1 <- cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10)
cat("lambda.min =", cv1$lambda.min, "（index", cv1$index[1], "）\n")
cat("lambda.1se =", cv1$lambda.1se, "（index", cv1$index[2], "）\n")
b1 <- as.matrix(coef(cv1, s = "lambda.min"))
nz1 <- rownames(b1)[b1[, 1] != 0]
cat("lambda.min 保留变量数：", length(nz1), "\n")
if (length(nz1)) {
  print(round(b1[nz1, 1, drop = FALSE], 5))
  cat("GNG7 在列：", "GNG7" %in% nz1, "；NFIA 在列：", "NFIA" %in% nz1, "\n")
}
b1s <- as.matrix(coef(cv1, s = "lambda.1se"))
cat("lambda.1se 保留变量数：", sum(b1s[, 1] != 0), "\n")

# ============ 3. 弹性网（alpha=0.5）============
cat("\n========== 3. 弹性网（alpha=0.5，10 折 CV，seed=20260913）==========\n")
set.seed(20260913)
cv05 <- cv.glmnet(X, y, family = "cox", alpha = 0.5, nfolds = 10)
cat("lambda.min =", cv05$lambda.min, "；lambda.1se =", cv05$lambda.1se, "\n")
b05 <- as.matrix(coef(cv05, s = "lambda.min"))
b05s <- as.matrix(coef(cv05, s = "lambda.1se"))
cat("弹性网 lambda.min 保留变量数：", sum(b05[, 1] != 0), "\n")
cat("弹性网 lambda.1se 保留变量数：", sum(b05s[, 1] != 0), "\n")

# ============ 4. 历史产物对账（lasso_cv.RData 内嵌 cv）============
cat("\n========== 4. 与历史 lasso_cv.RData 对账 ==========\n")
cvh <- env$cv
cat("历史 cv 类型：", class(cvh)[1], "；lambda.min =", cvh$lambda.min,
    "；lambda.1se =", cvh$lambda.1se, "\n")
bh <- as.matrix(coef(cvh$glmnet.fit, s = cvh$lambda.min))
nzh <- rownames(bh)[bh[, 1] != 0]
cat("历史 lambda.min 保留变量数：", length(nzh),
    "；GNG7 在列：", "GNG7" %in% nzh, "；NFIA 在列：", "NFIA" %in% nzh, "\n")

save(cv1, cv05, file = "D:/gwas/lasso_rerun_20260913.RData")
cat("\n已保存 D:/gwas/lasso_rerun_20260913.RData\n")
cat("========== 脚本 I 完成 ==========\n")
