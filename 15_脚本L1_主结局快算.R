# =====================================================================
# 脚本 L1 · MR 二轮补跑-主结局快算版（2026-09-14）
# 复用 2026-09-13 封存的中间产物（dat_mr_main.rds，69 SNP）：
#  ① PRESSO corrected IVW（按 Distortion Outliers Indices 剔除后重估）
#  ② 离群 SNP rsID 列表
#  ③ I²（IVW/Egger）与 DerSimonian-Laird tau²
#  ④ MR-RAPS（mr.raps 0.4.3）
#  ⑤ Steiger 方向性检验
#  ⑥ 反向 MR 尝试（结局 GWAS 无 P<5e-8 显著工具则如实报告）
# 前置：dat_mr_main.rds（已有）；TwoSampleMR、MRPRESSO、mr.raps 已装
# 运行：Rscript 本文件（约 1-3 分钟，仅主结局；复制结局与 MVMR 见脚本 L2）
# =====================================================================
suppressMessages({library(TwoSampleMR); library(MRPRESSO)})
Sys.setenv(HTTP_PROXY = "http://127.0.0.1:7897", HTTPS_PROXY = "http://127.0.0.1:7897",
           http_proxy = "http://127.0.0.1:7897", https_proxy = "http://127.0.0.1:7897")
jwt <- ieugwasr::get_opengwas_jwt(); stopifnot("Token 未生效" = nchar(jwt) > 20)
cat("Token OK\n")
dat_mr <- readRDS("D:/gwas/tmp/dat_mr_main.rds")
cat("主结局有效 SNP：", nrow(dat_mr), "\n")

# ---------- 1. PRESSO corrected IVW（复用 2026-09-13 封存 Distortion Indices）----------
cat("\n========== 1. PRESSO corrected IVW ==========\n")
mp <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = dat_mr, NbDistribution = 10000, SignifThreshold = 0.05)
ot <- mp$`MR-PRESSO results`$`Outlier Test`
print(mp$`MR-PRESSO results`$`Global Test`)
print(ot)
print(mp$`MR-PRESSO results`$`Distortion Test`)
sig_rows <- which(ot$Pvalue < 0.05)
cat("\nOutlier Test 显著行（rownames）：", rownames(ot)[sig_rows], "\n")
dist_idx <- mp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
cat("Distortion Outliers Indices：", dist_idx, "\n")
if (length(dist_idx) > 0) {
  cat("离群 SNP rsID：", dat_mr$SNP[dist_idx], "\n")
  dat_clean <- dat_mr[-dist_idx, ]
  cat("剔除后 SNP 数：", nrow(dat_clean), "\n")
  cat("\n--- 剔除离群 SNP 后 IVW（PRESSO corrected）---\n")
  print(generate_odds_ratios(mr(dat_clean, method_list = c("mr_ivw_mre", "mr_ivw_fe"))))
}
saveRDS(list(mp = mp, dat_mr = dat_mr), "D:/gwas/tmp/presso_main_20260914.rds")
cat("（已保存 presso_main_20260914.rds）\n")

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

# ---------- 5. 反向 MR ----------
cat("\n========== 5. 反向 MR（甲状腺癌→TSH）==========\n")
res5 <- tryCatch(ieugwasr::associations(variants = NULL, id = "ebi-a-GCST90013867", pval = 5e-8),
                 error = function(e) NULL)
if (is.null(res5) || nrow(res5) == 0) {
  cat("结局 GWAS（ebi-a-GCST90013867）无 P<5e-8 的独立全基因组显著位点 → 反向 MR 无法执行（如实报告）\n")
} else {
  cat("检出显著位点：", nrow(res5), "\n")
  print(head(res5))
}

cat("\n========== 脚本 L1 完成 ==========\n")
