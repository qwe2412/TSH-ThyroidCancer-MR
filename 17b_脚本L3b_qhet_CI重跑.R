# 脚本 L3b · qhet_mvmr CI=TRUE 重跑（2026-09-15）
suppressMessages({library(TwoSampleMR); library(MVMR); library(MendelianRandomization)})
Sys.setenv(HTTP_PROXY = "http://127.0.0.1:7897", HTTPS_PROXY = "http://127.0.0.1:7897",
           http_proxy = "http://127.0.0.1:7897", https_proxy = "http://127.0.0.1:7897")
jwt <- ieugwasr::get_opengwas_jwt(); stopifnot("Token 未生效" = nchar(jwt) > 20)
cat("Token OK\n")
exp <- readRDS("D:/gwas/tmp/exp_70snps.rds")
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

cat("\n--- MVMR-IVW ---\n")
mv_in <- mr_mvinput(bx = BX, bxse = BXse, by = BY, byse = BYse)
fit_mvmr <- mr_mvivw(mv_in)
print(fit_mvmr)
cat("TSH b =", fit_mvmr@Estimate[1], "; se =", fit_mvmr@StdError[1],
    "; p =", fit_mvmr@Pvalue[1], "\n")
cat("BMI b =", fit_mvmr@Estimate[2], "; se =", fit_mvmr@StdError[2],
    "; p =", fit_mvmr@Pvalue[2], "\n")
cat("睾酮 b =", fit_mvmr@Estimate[3], "; se =", fit_mvmr@StdError[3],
    "; p =", fit_mvmr@Pvalue[3], "\n")

F.data <- format_mvmr(BXGs = BX, BYG = BY, seBXGs = BXse, seBYG = BYse, RSID = m$SNP)
pheno_cor <- matrix(c(1.00, 0.05, 0.00, 0.05, 1.00, -0.10, 0.00, -0.10, 1.00),
                    nrow = 3, byrow = TRUE)
cat("\n--- qhet_mvmr（1000 次，CI=TRUE）---\n")
qh <- tryCatch(qhet_mvmr(F.data, pheno_cor, CI = TRUE, iterations = 1000),
               error = function(e) { cat("失败：", conditionMessage(e), "\n"); NULL })
if (!is.null(qh)) {
  cat("class:", paste(class(qh), collapse = ","), "\n")
  cat("结构 str:\n"); str(qh, max.level = 3)
  cat("print:\n"); print(qh)
}
saveRDS(list(m = m, fit_mvmr = fit_mvmr, qh = qh), "D:/gwas/tmp/mvmr_qhet_CI_20260915.rds")
cat("（已保存 mvmr_qhet_CI_20260915.rds）\n")
cat("\n========== 脚本 L3b 完成 ==========\n")
