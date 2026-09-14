# =====================================================================
# 脚本 J · Lasso/弹性网变量选择稳定性 + 重复交叉验证（2026-09-14）
# ------------------------------------------------------------
# 目的：回应 DeepSeek 二轮审稿 "变量选择极不稳定" 与
#   "报告变量选择稳定性：Bootstrap 选择频率、LASSO 路径稳定性、
#    重复交叉验证"。
#  ① Bootstrap 变量选择频率（B=200）：每次有放回重采样 → cv.glmnet
#     (alpha=1, 10 折, seed=20260914+i) → lambda.min 路径记录保留变量
#     → 输出 GNG7/NFIA/各变量被保留频率
#  ② 重复 5×10 折交叉验证：拟合 GNG7+NFIA Cox，评估测试折 C-index
#     → 输出 C-index 均值/范围
#  ③ 报告 GNG7/NFIA 在 bootstrap 中的系数方向一致性
# 前置：D:/gwas/lasso_cv.RData（dat 409×68）；glmnet、survival
# 运行：Rscript 本文件（约 8-15 分钟）
# =====================================================================
suppressMessages(library(glmnet))
suppressMessages(library(survival))

env <- new.env(); load("D:/gwas/lasso_cv.RData", envir = env)
dat <- env$dat
time_col <- "DFS_time"; event_col <- "DFS_event"
genes <- setdiff(colnames(dat), c("patient", "sample", time_col, event_col,
                                  "DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS"))
X <- as.matrix(dat[, genes])
X[is.na(X)] <- median(X, na.rm = TRUE)
y <- Surv(dat[[time_col]], dat[[event_col]])
cat("样本:", nrow(dat), "事件:", sum(dat[[event_col]]), "基因:", length(genes), "\n")

# ============ 1. Bootstrap 变量选择频率（B=200）============
cat("\n========== 1. Bootstrap 变量选择频率（B = 200）==========\n")
B <- 200
keep_freq <- numeric(length(genes)); names(keep_freq) <- genes
coef_dir <- matrix(0, length(genes), B, dimnames = list(genes, NULL))
gng7_freq <- nfia_freq <- 0
set.seed(20260914)
for (i in seq_len(B)) {
  idx <- sample(nrow(dat), replace = TRUE)
  Xb <- X[idx, , drop = FALSE]
  yb <- y[idx, , drop = FALSE]
  if (length(unique(yb[, 2])) < 2) next            # 单事件 bootstrap 样本跳过
  cvb <- tryCatch(cv.glmnet(Xb, yb, family = "cox", alpha = 1, nfolds = 10),
                  error = function(e) NULL)
  if (is.null(cvb)) next
  bb <- as.matrix(coef(cvb, s = "lambda.min"))
  nz <- rownames(bb)[bb[, 1] != 0]
  keep_freq[nz] <- keep_freq[nz] + 1
  coef_dir[, i] <- bb[, 1]
  if ("GNG7" %in% nz) gng7_freq <- gng7_freq + 1
  if ("NFIA" %in% nz) nfia_freq <- nfia_freq + 1
  if (i %% 50 == 0) cat("  bootstrap", i, "done\n")
}
cat("有效 bootstrap 迭代:", sum(coef_dir[1, ] != 0 | colSums(abs(coef_dir)) > 0), "\n")
top <- sort(keep_freq, decreasing = TRUE)[1:15]
cat("\n保留频率 Top 15 变量（次数 /", B, "）：\n")
print(top)
cat("\nGNG7 保留频率:", gng7_freq, "/", B, "（", round(100*gng7_freq/B,1), "%）\n")
cat("NFIA 保留频率:", nfia_freq, "/", B, "（", round(100*nfia_freq/B,1), "%）\n")
g7 <- coef_dir["GNG7", ]; n7 <- coef_dir["NFIA", ]
cat("GNG7 非零迭代中系数均值为负:", mean(g7[g7 != 0]) < 0, "（均值", round(mean(g7[g7 != 0]), 4), "）\n")
cat("NFIA 非零迭代中系数均值为负:", mean(n7[n7 != 0]) < 0, "（均值", round(mean(n7[n7 != 0]), 4), "）\n")

# ============ 2. 重复 5×10 折交叉验证 C-index（GNG7+NFIA Cox）============
cat("\n========== 2. 重复 5×10 折交叉验证 C-index（GNG7+NFIA）==========\n")
rep_cindex <- numeric(50)
k <- 0
set.seed(20260914)
for (r in 1:5) {
  folds <- sample(rep(1:10, length.out = nrow(dat)))
  for (f in 1:10) {
    te <- folds == f
    fit <- tryCatch(coxph(Surv(dat[[time_col]][!te], dat[[event_col]][!te]) ~
                            X[!te, "GNG7"] + X[!te, "NFIA"]),
                    error = function(e) NULL)
    if (is.null(fit)) next
    lp <- predict(fit, newdata = data.frame(GNG7 = X[te, "GNG7"], NFIA = X[te, "NFIA"]))
    # 基于风险分层的 C-index（无 concordance 依赖：用 survival::concordance）
    cc <- tryCatch(concordance(Surv(dat[[time_col]][te], dat[[event_col]][te]) ~ lp)$concordance,
                   error = function(e) NA)
    k <- k + 1; rep_cindex[k] <- cc
  }
}
rep_cindex <- rep_cindex[1:k]
cat("重复 CV C-index：均值 =", round(mean(rep_cindex, na.rm = TRUE), 4),
    "；范围 =", paste(round(range(rep_cindex, na.rm = TRUE), 4), collapse = " ~ "),
    "；NA 数 =", sum(is.na(rep_cindex)), "\n")

save(keep_freq, gng7_freq, nfia_freq, rep_cindex, coef_dir,
     file = "D:/gwas/lasso_boot_20260914.RData")
cat("\n已保存 D:/gwas/lasso_boot_20260914.RData\n")
cat("========== 脚本 J 完成 ==========\n")
