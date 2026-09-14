# =====================================================================
# 脚本 L1b · MR 二轮补跑-主结局纯离线兜底版（2026-09-14）
# 不重跑 PRESSO：直接复用 2026-09-13 脚本 A 封存的
#   Distortion Outliers Indices=[1,12,33,66] 与 Outlier 显著行（rownames 1/12/34/68），
#   剔除后重估 IVW（PRESSO corrected 等价结果）；再算 I²/tau²/RAPS/Steiger。
# 全部离线、秒级完成。
# =====================================================================
suppressMessages(library(TwoSampleMR))
dat_mr <- readRDS("D:/gwas/tmp/dat_mr_main.rds")
cat("主结局有效 SNP：", nrow(dat_mr), "\n")

# ---------- 1. PRESSO corrected IVW（封存 Distortion Indices）----------
cat("\n========== 1. PRESSO corrected IVW（按封存 Distortion Indices）==========\n")
dist_idx <- c(1, 12, 33, 66)   # 2026-09-13 脚本 A 实测 Distortion Outliers Indices
cat("Distortion Outliers Indices：", dist_idx, "\n")
cat("离群 SNP rsID：", dat_mr$SNP[dist_idx], "\n")
dat_clean <- dat_mr[-dist_idx, ]
cat("剔除后 SNP 数：", nrow(dat_clean), "\n")
cat("\n--- 剔除离群 SNP 后 IVW（PRESSO corrected）---\n")
print(generate_odds_ratios(mr(dat_clean, method_list = c("mr_ivw_mre", "mr_ivw_fe"))))

# ---------- 2. I² 与 tau² ----------
cat("\n========== 2. 异质性 I² 与 tau²（DL）==========\n")
het <- mr_heterogeneity(dat_mr, method_list = c("mr_ivw", "mr_egger_regression"))
print(het)
for (i in seq_len(nrow(het))) {
  Q <- het$Q[i]; df <- het$Q_df[i]
  cat(het$method[i], "：Q =", Q, "df =", df, "→ I² =", round(ifelse(Q > df, (Q - df) / Q * 100, 0), 2), "%\n")
}
w <- 1 / dat_mr$se.outcome^2
Q_ivw <- het$Q[1]; k <- nrow(dat_mr)
tau2 <- (Q_ivw - (k - 1)) / (sum(w) - sum(w^2) / sum(w))
cat("tau²（DL，IVW Q）：", round(tau2, 5), "；tau（SD）：", round(sqrt(tau2), 5), "\n")

# ---------- 3. MR-RAPS ----------
cat("\n========== 3. MR-RAPS ==========\n")
cat("mr.raps 已安装：", requireNamespace("mr.raps", quietly = TRUE), "\n")
if (requireNamespace("mr.raps", quietly = TRUE)) {
  rr <- mr.raps::mr.raps.all(dat_mr$beta.exposure, dat_mr$beta.outcome,
                             dat_mr$se.exposure, dat_mr$se.outcome)
  print(rr)
}

# ---------- 4. Steiger ----------
cat("\n========== 4. Steiger 方向性检验 ==========\n")
d1 <- dat_mr
d1$r.exposure <- get_r_from_pn(d1$pval.exposure, d1$samplesize.exposure)
d1$r.outcome  <- get_r_from_pn(d1$pval.outcome, 491974)
d1$samplesize.outcome <- 491974
ste <- tryCatch(directionality_test(d1), error = function(e) {
  cat("Steiger 失败：", conditionMessage(e), "\n"); NULL })
if (!is.null(ste)) print(ste)

cat("\n========== 脚本 L1b 完成 ==========\n")
