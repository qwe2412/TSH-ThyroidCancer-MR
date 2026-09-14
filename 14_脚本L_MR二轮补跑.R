# =====================================================================
# 脚本 L · MR 二轮补跑（回应 DeepSeek 二轮审稿 2026-09-14）
# ------------------------------------------------------------
# 目的：
#  ① MR-PRESSO outlier-corrected 估计（按 Distortion Outliers Indices 剔除）
#     + 离群 SNP rsID 列表 + I² / tau²（DerSimonian-Laird）
#  ② 复制结局①/② 完整敏感性：五法 MR + Q + Egger 截距 + PRESSO(2000)
#  ③ MR-RAPS（稳健调整轮廓评分；若包可装）
#  ④ Steiger 方向性检验（r 由 p、n 计算）
#  ⑤ 反向 MR（甲状腺癌→TSH）：若结局 GWAS 无显著工具则如实报告
#  ⑥ MVMR 全暴露 qhet_mvmr 估计（TSH/BMI/睾酮）
# 前置：Token 在 ~/.Renviron；xlsx 两个；MVMR 0.4.8 已装
# 运行：Rscript 本文件（后台，PRESSO 较慢）
# =====================================================================
suppressMessages({library(readxl); library(TwoSampleMR); library(MRPRESSO)})

jwt <- ieugwasr::get_opengwas_jwt()
stopifnot("Token 未生效" = nchar(jwt) > 20)
cat("Token OK\n")

# ---------- 复用脚本 A 1-6 步 ----------
num <- function(x) as.numeric(gsub(",", "", trimws(x)))
extract_ivs <- function(raw, hdr) {
  d <- raw[(hdr + 1):nrow(raw), ]
  data.frame(rsid = trimws(as.character(d[[2]])),
             effect_allele = trimws(as.character(d[[5]])),
             other_allele  = trimws(as.character(d[[4]])),
             eaf  = num(d[[8]]), N = num(d[[9]]),
             beta = num(d[[10]]), se = num(d[[11]]), pval = num(d[[12]]))
}
raw4 <- as.data.frame(read_excel("D:/gwas/41467_2020_17718_MOESM4_ESM.xlsx", col_names = FALSE, col_types = "text"))
raw5 <- as.data.frame(read_excel("D:/gwas/41467_2020_17718_MOESM5_ESM.xlsx", col_names = FALSE, col_types = "text"))
ivs <- rbind(extract_ivs(raw4, 3), extract_ivs(raw5, 3))
ivs <- ivs[!is.na(ivs$rsid) & !is.na(ivs$beta) & !is.na(ivs$se), ]
stopifnot(nrow(ivs) == 99)
exp <- format_data(ivs, type = "exposure",
                   snp_col = "rsid", beta_col = "beta", se_col = "se",
                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                   eaf_col = "eaf", pval_col = "pval", samplesize_col = "N")
exp <- clump_data(exp, clump_kb = 10000, clump_r2 = 0.001)
cat("clump 后 SNP 数：", nrow(exp), "\n")
saveRDS(exp, "D:/gwas/tmp/exp_70snps.rds")

# ---------- 主结局 ----------
out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST90018929")
dat <- harmonise_data(exp, out)
dat_mr <- dat[dat$mr_keep == TRUE, ]
cat("主结局有效 SNP：", nrow(dat_mr), "\n")
saveRDS(dat_mr, "D:/gwas/tmp/dat_mr_main.rds")

# ============ 1. PRESSO 10000 + corrected IVW ============
cat("\n========== 1. MR-PRESSO（主结局，10000 次）==========\n")
mp <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = dat_mr, NbDistribution = 10000, SignifThreshold = 0.05)
print(mp$`MR-PRESSO results`$`Global Test`)
print(mp$`MR-PRESSO results`$`Outlier Test`)
print(mp$`MR-PRESSO results`$`Distortion Test`)
ot <- mp$`MR-PRESSO results`$`Outlier Test`
sig_rows <- which(ot$Pvalue < 0.05)
cat("\nOutlier Test 显著行（rownames）：", rownames(ot)[sig_rows], "\n")
dist_idx <- mp$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
cat("Distortion Outliers Indices：", dist_idx, "\n")
if (length(dist_idx) > 0) {
  cat("离群 SNP rsID（dat_mr 位置）：", dat_mr$SNP[dist_idx], "\n")
  dat_clean <- dat_mr[-dist_idx, ]
  cat("\n--- 剔除离群 SNP 后 IVW（PRESSO corrected）---\n")
  print(generate_odds_ratios(mr(dat_clean, method_list = c("mr_ivw_mre", "mr_ivw_fe"))))
  cat("剔除后 SNP 数：", nrow(dat_clean), "\n")
}

# ============ 2. I² 与 tau² ============
cat("\n========== 2. 异质性 I² 与 tau²（DerSimonian-Laird）==========\n")
het <- mr_heterogeneity(dat_mr, method_list = c("mr_ivw", "mr_egger_regression"))
print(het)
for (i in seq_len(nrow(het))) {
  Q <- het$Q[i]; df <- het$Q_df[i]
  I2 <- ifelse(Q > df, (Q - df) / Q * 100, 0)
  cat(het$method[i], "：Q =", Q, "df =", df, "→ I² =", round(I2, 2), "%\n")
}
w <- 1 / dat_mr$se.outcome^2
Q_ivw <- het$Q[1]; k <- nrow(dat_mr)
tau2 <- (Q_ivw - (k - 1)) / (sum(w) - sum(w^2) / sum(w))
cat("tau²（DL，IVW Q）：", round(tau2, 5), "\n")
cat("tau（SD）：", round(sqrt(tau2), 5), "\n")

# ============ 3. MR-RAPS（若可装）===========
cat("\n========== 3. MR-RAPS ==========\n")
rap_ok <- requireNamespace("mr.raps", quietly = TRUE)
cat("mr.raps 已安装：", rap_ok, "\n")
if (rap_ok) {
  rr <- mr.raps::mr.raps.all(dat_mr$beta.exposure, dat_mr$beta.outcome,
                             dat_mr$se.exposure, dat_mr$se.outcome)
  print(rr)
}

# ============ 4. Steiger 方向性检验 ============
cat("\n========== 4. Steiger 方向性检验 ==========\n")
d1 <- dat_mr
d1$r.exposure <- get_r_from_pn(d1$pval.exposure, d1$samplesize.exposure)
d1$r.outcome  <- get_r_from_pn(d1$pval.outcome, 491974)
d1$samplesize.outcome <- 491974
ste <- tryCatch(directionality_test(d1), error = function(e) {
  cat("Steiger 失败：", conditionMessage(e), "\n"); NULL })
if (!is.null(ste)) print(ste)

# ============ 5. 反向 MR（甲状腺癌 → TSH）===========
cat("\n========== 5. 反向 MR ==========\n")
# 主结局病例仅 1054，先试复制①（n=407746）是否有 P<5e-8 的甲状腺癌显著 SNP
tryCatch({
  sig_snps <- ieugwasr::associations(variants = NULL, id = "ebi-a-GCST90013867",
                                     pval = 5e-8)  # NULL variants + pval 过滤
}, error = function(e) { cat("反向 MR 提取失败：", conditionMessage(e), "\n") })
# 若无显著位点，ieugwasr 返回空/报错 → 如实报告
cat("（反向 MR 需要结局 GWAS 的独立全基因组显著位点；若上一步无输出则无法执行）\n")

# ============ 6. 复制结局敏感性 ============
run_rep_sens <- function(id, label) {
  cat("\n========== 6.", label, id, "==========\n")
  out_i <- extract_outcome_data(snps = exp$SNP, outcomes = id)
  dat_i <- harmonise_data(exp, out_i)
  dat_i <- dat_i[dat_i$mr_keep == TRUE, ]
  cat("有效 SNP：", nrow(dat_i), "\n")
  cat("--- 五法 MR ---\n")
  print(generate_odds_ratios(mr(dat_i)))
  cat("--- 异质性 ---\n")
  het_i <- mr_heterogeneity(dat_i, method_list = c("mr_ivw", "mr_egger_regression"))
  print(het_i)
  for (j in seq_len(nrow(het_i))) {
    Q <- het_i$Q[j]; df <- het_i$Q_df[j]
    cat(het_i$method[j], "：I² =", round(ifelse(Q > df, (Q - df) / Q * 100, 0), 2), "%\n")
  }
  cat("--- Egger 截距 ---\n")
  print(mr_pleiotropy_test(dat_i))
  cat("--- PRESSO（2000 次）---\n")
  mp_i <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                    SdOutcome = "se.outcome", SdExposure = "se.exposure",
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                    data = dat_i, NbDistribution = 2000, SignifThreshold = 0.05)
  print(mp_i$`MR-PRESSO results`$`Global Test`)
  print(mp_i$`MR-PRESSO results`$`Distortion Test`)
  saveRDS(dat_i, sprintf("D:/gwas/tmp/dat_mr_rep_%s.rds", gsub("[:/]", "_", id)))
}
run_rep_sens("ebi-a-GCST90013867", "复制结局①")
run_rep_sens("finn-b-C3_THYROID_GLAND", "复制结局②")

# ============ 7. MVMR 全暴露 qhet ============
cat("\n========== 7. MVMR 全暴露 qhet_mvmr ==========\n")
if (requireNamespace("MVMR", quietly = TRUE)) {
  library(MVMR)
  # 重建 F.data（同脚本 11）：TSH 工具骨架 + BMI/睾酮/结局效应
  bmi_eff <- extract_outcome_data(snps = exp$SNP, outcomes = "ieu-b-40")
  tt_eff  <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST90014013")
  d_bmi <- harmonise_data(exp, bmi_eff); d_bmi <- d_bmi[d_bmi$mr_keep == TRUE, ]
  d_tt  <- harmonise_data(exp, tt_eff);  d_tt  <- d_tt[d_tt$mr_keep == TRUE, ]
  keep <- Reduce(intersect, list(dat_mr$SNP, d_bmi$SNP, d_tt$SNP))
  cat("三暴露+结局有效 SNP：", length(keep), "\n")
  m <- data.frame(
    SNP = keep,
    b_tsh = dat_mr$beta.exposure[match(keep, dat_mr$SNP)],  se_tsh = dat_mr$se.exposure[match(keep, dat_mr$SNP)],
    b_bmi = d_bmi$beta.outcome[match(keep, d_bmi$SNP)],     se_bmi = d_bmi$se.outcome[match(keep, d_bmi$SNP)],
    b_tt  = d_tt$beta.outcome[match(keep, d_tt$SNP)],       se_tt  = d_tt$se.outcome[match(keep, d_tt$SNP)],
    b_out = dat_mr$beta.outcome[match(keep, dat_mr$SNP)],   se_out = dat_mr$se.outcome[match(keep, dat_mr$SNP)])
  BX <- as.matrix(m[, c("b_tsh","b_bmi","b_tt")]); BY <- m$b_out
  BXse <- as.matrix(m[, c("se_tsh","se_bmi","se_tt")]); BYse <- m$se_out
  F.data <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = BXse, seBYG = BYse, RSID = m$SNP)
  pheno_cor <- matrix(c(1.00, 0.05, 0.00, 0.05, 1.00, -0.10, 0.00, -0.10, 1.00),
                      nrow = 3, byrow = TRUE)
  qh <- tryCatch(qhet_mvmr(F.data, pheno_cor, CI = FALSE, iterations = 1000),
                 error = function(e) { cat("失败：", conditionMessage(e), "\n"); NULL })
  if (!is.null(qh)) {
    cat("qhet_mvmr 对象 class：", paste(class(qh), collapse = ","), "\n")
    cat("结构：\n"); str(qh, max.level = 3)
    cat("打印：\n"); print(qh)
  }
} else { cat("MVMR 包未安装\n") }

cat("\n========== 脚本 L 完成 ==========\n")
