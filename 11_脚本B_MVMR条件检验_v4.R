# =====================================================================
# 脚本 B v4 · MVMR 条件检验完整版（严格审查补跑 2026-09-13）
# ------------------------------------------------------------
# 目的：补 DeepSeek 审稿意见要求但缺失的 MVMR 报告项：
#   ① 各暴露（BMI/睾酮）本身的 MVMR 效应
#   ② 条件 F 统计量（弱工具检验，Sanderson & Windmeijer 2016 / sim.9133）
#   ③ 条件异质性 Q（pleiotropy_mvmr，水平多效性检验）
#   ④ Q 统计量最小化稳健估计（qhet_mvmr）
# 方法：TSH 工具 SNP（clump 后，与脚本 A 一致）为骨架，提取
#   BMI(ieu-b-40)/睾酮(ebi-a-GCST90014013)/结局(ebi-a-GCST90018929) 效应，
#   三暴露 MVMR。
# 说明：SNP 对暴露效应的协方差 gencov 用 phenocov_mvmr（表型相关近似）
#   与 0（独立假设）双版本报告；表型相关为文献近似值，已在输出中注明假设。
# 前置：~/.Renviron 有 OPENGWAS_JWT；MVMR 0.4.8 已装（GitHub WSpiller/mvmr）
# 运行：Rscript 本文件（独立拉数据，不依赖脚本 A 会话）
# =====================================================================
suppressMessages({library(TwoSampleMR); library(MVMR); library(readxl)})

jwt <- ieugwasr::get_opengwas_jwt()
stopifnot("Token 未生效" = nchar(jwt) > 20)
cat("Token 长度:", nchar(jwt), "\n")

# ---------- 1. TSH 工具变量（同脚本 A：MOESM4+MOESM5 → clump） ----------
cat("\n========== 1. TSH 工具变量 ==========\n")
num <- function(x) as.numeric(gsub(",", "", trimws(x)))
extract_ivs <- function(raw, hdr) {
  d <- raw[(hdr + 1):nrow(raw), ]
  data.frame(rsid=trimws(as.character(d[[2]])), chr_pos=trimws(as.character(d[[3]])),
             effect_allele=trimws(as.character(d[[5]])), other_allele=trimws(as.character(d[[4]])),
             eaf=num(d[[8]]), N=num(d[[9]]), beta=num(d[[10]]), se=num(d[[11]]), pval=num(d[[12]]))
}
raw4 <- as.data.frame(read_excel("D:/gwas/41467_2020_17718_MOESM4_ESM.xlsx", col_names=FALSE, col_types="text"), stringsAsFactors=FALSE)
raw5 <- as.data.frame(read_excel("D:/gwas/41467_2020_17718_MOESM5_ESM.xlsx", col_names=FALSE, col_types="text"), stringsAsFactors=FALSE)
ivs <- rbind(extract_ivs(raw4, 3), extract_ivs(raw5, 3))
ivs <- ivs[!is.na(ivs$rsid) & !is.na(ivs$beta) & !is.na(ivs$se) & !is.na(ivs$pval), ]
cat("TSH 原始 SNP:", nrow(ivs), "\n")
exp <- format_data(ivs, type="exposure", snp_col="rsid", beta_col="beta", se_col="se",
                   effect_allele_col="effect_allele", other_allele_col="other_allele",
                   eaf_col="eaf", pval_col="pval", samplesize_col="N")
exp <- clump_data(exp, clump_kb=10000, clump_r2=0.001)
cat("clump 后 SNP:", nrow(exp), "\n")

# ---------- 2. 效应提取 ----------
cat("\n========== 2. 效应提取 ==========\n")
bmi_eff <- extract_outcome_data(snps=exp$SNP, outcomes="ieu-b-40")
tt_eff  <- extract_outcome_data(snps=exp$SNP, outcomes="ebi-a-GCST90014013")
out     <- extract_outcome_data(snps=exp$SNP, outcomes="ebi-a-GCST90018929")
cat("BMI 行数:", nrow(bmi_eff), "；睾酮:", nrow(tt_eff), "；结局:", nrow(out), "\n")

d_out <- harmonise_data(exp, out);     d_out <- d_out[d_out$mr_keep == TRUE, ]
d_bmi <- harmonise_data(exp, bmi_eff); d_bmi <- d_bmi[d_bmi$mr_keep == TRUE, ]
d_tt  <- harmonise_data(exp, tt_eff);  d_tt  <- d_tt[d_tt$mr_keep == TRUE, ]
cat("结局协调:", nrow(d_out), "；BMI:", nrow(d_bmi), "；睾酮:", nrow(d_tt), "\n")

keep <- Reduce(intersect, list(d_out$SNP, d_bmi$SNP, d_tt$SNP))
cat("三暴露+结局均有效的 SNP 数:", length(keep), "\n")
m <- data.frame(
  SNP=keep,
  b_tsh=d_out$beta.exposure[match(keep, d_out$SNP)],  se_tsh=d_out$se.exposure[match(keep, d_out$SNP)],
  b_bmi=d_bmi$beta.outcome[match(keep, d_bmi$SNP)],   se_bmi=d_bmi$se.outcome[match(keep, d_bmi$SNP)],
  b_tt =d_tt$beta.outcome[match(keep, d_tt$SNP)],     se_tt =d_tt$se.outcome[match(keep, d_tt$SNP)],
  b_out=d_out$beta.outcome[match(keep, d_out$SNP)],   se_out=d_out$se.outcome[match(keep, d_out$SNP)])
BX <- as.matrix(m[, c("b_tsh","b_bmi","b_tt")])
BXse <- as.matrix(m[, c("se_tsh","se_bmi","se_tt")])
BY <- m$b_out; BYse <- m$se_out

# ---------- 3. format_mvmr ----------
F.data <- format_mvmr(BXGs=BX, BYG=BY, seBXGs=BXse, seBYG=BYse, RSID=m$SNP)
cat("format_mvmr 完成，维度:", dim(F.data), "\n")
# F.data 列：X1..Xk, Y, X1se..Xkse, Yse, SNP? —— 用 7:9 是 se 列（k=3 时列 7-9）

# ---------- 4. 表型相关协方差（phenocov_mvmr） ----------
# 近似表型相关矩阵（文献近似值，注明假设）：
#   TSH–BMI r≈0.05，TSH–睾酮 r≈0.00，BMI–睾酮 r≈-0.10（男性为主样本）
cat("\n========== 3. 表型相关协方差 ==========\n")
pheno_cor <- matrix(c(1.00, 0.05, 0.00,
                      0.05, 1.00, -0.10,
                      0.00, -0.10, 1.00), nrow=3, byrow=TRUE)
cat("表型相关矩阵（近似假设）:\n"); print(pheno_cor)
Xcovmat <- phenocov_mvmr(pheno_cor, F.data[, 7:9])
cat("phenocov_mvmr 输出矩阵数:", length(Xcovmat), "\n")

# ---------- 5. 条件 F（弱工具检验） ----------
cat("\n========== 4. 条件 F 统计量 ==========\n")
cat("--- 4a. gencov=0（独立假设）---\n")
cf0 <- strength_mvmr(r_input=F.data, gencov=0)
print(cf0)
cat("--- 4b. gencov=phenocov（表型相关校正）---\n")
cf1 <- tryCatch(strength_mvmr(r_input=F.data, gencov=Xcovmat),
                error=function(e) { cat("失败:", conditionMessage(e), "\n"); NULL })
if (!is.null(cf1)) print(cf1)
write.csv(as.data.frame(cf0), "D:/gwas/mvmr_condF_20260913.csv")

# ---------- 6. 条件异质性 Q（pleiotropy_mvmr） ----------
cat("\n========== 5. 条件异质性 Q（水平多效性检验）==========\n")
cat("--- 5a. gencov=0 ---\n")
q0 <- pleiotropy_mvmr(r_input=F.data, gencov=0)
print(q0)
cat("--- 5b. gencov=phenocov ---\n")
q1 <- tryCatch(pleiotropy_mvmr(r_input=F.data, gencov=Xcovmat),
               error=function(e) { cat("失败:", conditionMessage(e), "\n"); NULL })
if (!is.null(q1)) print(q1)
if (!is.null(q0)) write.csv(as.data.frame(q0), "D:/gwas/mvmr_hetQ_20260913.csv")

# ---------- 7. MVMR-IVW（主估计） ----------
cat("\n========== 6. MVMR-IVW 效应 ==========\n")
ivw <- ivw_mvmr(r_input=F.data)
print(ivw)
res_mv <- data.frame(Exposure=c("TSH","BMI","Testosterone"),
                     b=ivw[,1], se=ivw[,2])
res_mv$P <- 2 * pnorm(-abs(res_mv$b / res_mv$se))
res_mv$OR <- exp(res_mv$b)
res_mv$OR_lci <- exp(res_mv$b - 1.96 * res_mv$se)
res_mv$OR_uci <- exp(res_mv$b + 1.96 * res_mv$se)
print(res_mv, digits=5)
write.csv(res_mv, "D:/gwas/mvmr_ivw_20260913.csv", row.names=FALSE)

# ---------- 8. Q 统计量最小化稳健估计（qhet_mvmr） ----------
cat("\n========== 7. Q 统计量最小化估计（qhet_mvmr, CI=FALSE）==========\n")
qh <- tryCatch(qhet_mvmr(F.data, pheno_cor, CI=FALSE, iterations=1000),
               error=function(e) { cat("失败:", conditionMessage(e), "\n"); NULL })
if (!is.null(qh)) print(qh)

cat("\n========== 脚本 B v4 完成 ==========\n")
