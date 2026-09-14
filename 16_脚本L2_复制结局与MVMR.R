# =====================================================================
# 脚本 L2 · MR 二轮补跑-复制结局与 MVMR 全暴露（2026-09-14）
#  ① 复制结局① ebi-a-GCST90013867：五法 MR + Q/I² + Egger 截距 + PRESSO(2000)
#  ② 复制结局② finn-b-C3_THYROID_GLAND：同上
#  ③ MVMR 重建 F.data + qhet_mvmr(iterations=1000) 全暴露估计（TSH/BMI/睾酮）
# 前置：exp_70snps.rds；Token；MVMR 0.4.8
# 运行：Rscript 本文件（约 10-20 分钟；每步自动保存 rds）
# =====================================================================
suppressMessages({library(TwoSampleMR); library(MRPRESSO); library(MVMR)})
Sys.setenv(HTTP_PROXY = "http://127.0.0.1:7897", HTTPS_PROXY = "http://127.0.0.1:7897",
           http_proxy = "http://127.0.0.1:7897", https_proxy = "http://127.0.0.1:7897")
jwt <- ieugwasr::get_opengwas_jwt(); stopifnot("Token 未生效" = nchar(jwt) > 20)
cat("Token OK\n")
exp <- readRDS("D:/gwas/tmp/exp_70snps.rds")

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
    cat(het_i$method[j], "：Q =", Q, "df =", df,
        "→ I² =", round(ifelse(Q > df, (Q - df) / Q * 100, 0), 2), "%\n")
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
  saveRDS(list(dat_i = dat_i, het_i = het_i, mp_i = mp_i),
          sprintf("D:/gwas/tmp/rep_%s.rds", gsub("[:/]", "_", id)))
  cat("（已保存 rep_", id, ".rds）\n", sep = "")
}
run_rep_sens("ebi-a-GCST90013867", "复制结局①")
run_rep_sens("finn-b-C3_THYROID_GLAND", "复制结局②")

# ============ 7. MVMR 全暴露 qhet ============
cat("\n========== 7. MVMR 全暴露 qhet_mvmr ==========\n")
dat_mr <- readRDS("D:/gwas/tmp/dat_mr_main.rds")
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
BX <- as.matrix(m[, c("b_tsh", "b_bmi", "b_tt")]); BY <- m$b_out
BXse <- as.matrix(m[, c("se_tsh", "se_bmi", "se_tt")]); BYse <- m$se_out
F.data <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = BXse, seBYG = BYse, RSID = m$SNP)
pheno_cor <- matrix(c(1.00, 0.05, 0.00, 0.05, 1.00, -0.10, 0.00, -0.10, 1.00),
                    nrow = 3, byrow = TRUE)
# MVMR-IVW 复核
fit_mvmr <- mr_mvivw(mr_mvinput(bx = BX, bxse = BXse, by = BY, byse = BYse))
print(fit_mvmr)
cat("--- qhet_mvmr（1000 次）---\n")
qh <- tryCatch(qhet_mvmr(F.data, pheno_cor, CI = FALSE, iterations = 1000),
               error = function(e) { cat("失败：", conditionMessage(e), "\n"); NULL })
if (!is.null(qh)) {
  cat("class:", paste(class(qh), collapse = ","), "\n")
  str(qh, max.level = 3)
  print(qh)
}
saveRDS(list(m = m, F.data = F.data, qh = qh, fit_mvmr = fit_mvmr),
        "D:/gwas/tmp/mvmr_qhet_20260914.rds")
cat("（已保存 mvmr_qhet_20260914.rds）\n")
cat("\n========== 脚本 L2 完成 ==========\n")
