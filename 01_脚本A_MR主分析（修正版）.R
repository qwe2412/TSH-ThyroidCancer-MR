# =====================================================================
# 全流程重跑 · 脚本 A：MR 主分析全流程（对应手册第 1-7 章）【修正版】
# ------------------------------------------------------------
# 覆盖：数据读取 → 99 SNP → LD clump → 结局提取 → 协调
#       → 五法 MR → 异质性/多效性/留一 → 随机效应 IVW → 出图
#       → MR-PRESSO（10000次）→ 剔除离群后 IVW → 复制结局 ×2
# ------------------------------------------------------------
# 【本次修正（2026-09-12）】
#   harmonise_data 不会删除无效行，而是标记 mr_keep=FALSE（本次运行中
#   rs7966590 等位基因不兼容、rs2979181 回文，2 行 mr_keep=FALSE）。
#   mr() 会自动过滤，但 mr_presso() 不过滤 → 导致"71 vs 69"不一致。
#   修正：协调后立即用 dat_mr <- dat[dat$mr_keep == TRUE, ]，
#   主分析、敏感性、MR-PRESSO 全部统一为 69 个有效 SNP。
# ------------------------------------------------------------
# 运行方式：RStudio 打开本文件 → Ctrl+A 全选 → 点 Run
# 前置：① Token 已在 ~/.Renviron 并重启过 R
#       ② D:/gwas/ 下两个 xlsx 存在
# 耗时：数据准备 1-2 分钟 + MR-PRESSO 10000 次约 10-30 分钟（耐心等）
# =====================================================================

# ---- 缺失包时先运行（只需一次）----
# install.packages(c("readxl", "TwoSampleMR", "MRPRESSO", "MendelianRandomization"))

library(readxl)
library(TwoSampleMR)
library(MRPRESSO)

# ================= 0. 运行前检查 =================
stopifnot("MOESM4 文件不存在" = file.exists("D:/gwas/41467_2020_17718_MOESM4_ESM.xlsx"))
stopifnot("MOESM5 文件不存在" = file.exists("D:/gwas/41467_2020_17718_MOESM5_ESM.xlsx"))
cat("✓ 数据文件存在\n")

jwt <- ieugwasr::get_opengwas_jwt()
stopifnot("Token 未生效" = nchar(jwt) > 20)
cat("✓ Token 已生效（长度", nchar(jwt), "）\n")

# ================= 1. 读取补充数据（手册 4.1） =================
raw4 <- as.data.frame(read_excel("D:/gwas/41467_2020_17718_MOESM4_ESM.xlsx",
                                 col_names = FALSE, col_types = "text"), stringsAsFactors = FALSE)
raw5 <- as.data.frame(read_excel("D:/gwas/41467_2020_17718_MOESM5_ESM.xlsx",
                                 col_names = FALSE, col_types = "text"), stringsAsFactors = FALSE)

# ================= 2. 提取 99 个 SNP（手册 4.2，表头行=第3行） =================
num <- function(x) as.numeric(gsub(",", "", trimws(x)))
extract_ivs <- function(raw, hdr) {
  d <- raw[(hdr + 1):nrow(raw), ]
  data.frame(
    rsid = trimws(as.character(d[[2]])),          # col2  = rsID
    chr_pos = trimws(as.character(d[[3]])),       # col3  = Position
    effect_allele = trimws(as.character(d[[5]])), # col5  = Alt
    other_allele  = trimws(as.character(d[[4]])), # col4  = Ref
    eaf  = num(d[[8]]),                           # col8  = Allele freq
    N    = num(d[[9]]),                           # col9  = N (total)
    beta = num(d[[10]]),                          # col10 = Effect（SD 单位）
    se   = num(d[[11]]),                          # col11 = SE
    pval = num(d[[12]])                           # col12 = P-value
  )
}
iv4 <- extract_ivs(raw4, 3)
iv5 <- extract_ivs(raw5, 3)
ivs <- rbind(iv4, iv5)
ivs <- ivs[!is.na(ivs$rsid) & !is.na(ivs$beta) & !is.na(ivs$se) & !is.na(ivs$pval), ]
stopifnot("SNP 提取数 ≠ 99" = nrow(ivs) == 99)
cat("✓ 提取 99 个 SNP\n")

# ================= 3. 格式化暴露（手册 4.3） =================
exp <- format_data(ivs, type = "exposure",
                   snp_col = "rsid", beta_col = "beta", se_col = "se",
                   effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                   eaf_col = "eaf", pval_col = "pval", samplesize_col = "N")

# ================= 4. LD clump（手册 5.1，EUR，10000kb，r²=0.001） =================
exp <- clump_data(exp, clump_kb = 10000, clump_r2 = 0.001)
cat("✓ clump 后独立 SNP 数：", nrow(exp), "\n")

# ================= 5. 结局提取（手册 5.2） =================
out <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST90018929")
cat("结局提取返回行数：", nrow(out), "（可能含重复/无效记录）\n")
stopifnot("结局匹配 SNP 数异常（<60）" = nrow(out) >= 60)

# ================= 6. 协调（手册 5.3）【关键修正】 =================
dat <- harmonise_data(exp, out)
# 只保留 mr_keep=TRUE 的有效协调 SNP（与主分析完全一致）
dat_mr <- dat[dat$mr_keep == TRUE, ]
cat("✓ 协调总行数：", nrow(dat), "；有效 SNP（mr_keep=TRUE）：", nrow(dat_mr), "\n")
if (nrow(dat_mr) < nrow(dat)) {
  cat("  已剔除无效协调行：", paste(dat$SNP[dat$mr_keep == FALSE], collapse = ", "), "\n")
}

# ================= 7. 五法 MR（手册 5.3 + 6.4，基于 dat_mr） =================
cat("\n========== 7. 五法 MR 主分析（nsnp =", nrow(dat_mr), "）==========\n")
res_main <- mr(dat_mr, method_list = c("mr_ivw", "mr_ivw_mre",
                                       "mr_egger_regression", "mr_weighted_median",
                                       "mr_simple_mode", "mr_weighted_mode"))
print(res_main)
cat("\n--- OR（95% CI）---\n")
print(generate_odds_ratios(res_main))

# ================= 8. 敏感性分析（手册 6.1-6.4，基于 dat_mr） =================
cat("\n========== 8.1 异质性（Cochran's Q）==========\n")
print(mr_heterogeneity(dat_mr, method_list = c("mr_ivw", "mr_ivw_mre", "mr_egger_regression")))

cat("\n========== 8.2 MR-Egger 截距（多效性）==========\n")
print(mr_pleiotropy_test(dat_mr))

cat("\n========== 8.3 留一法 ==========\n")
res_loo <- mr_leaveoneout(dat_mr)
cat("留一法 b 范围：", range(res_loo$b), "\n")

cat("\n========== 8.4 随机效应 IVW（对照）==========\n")
print(mr(dat_mr, method_list = "mr_ivw_mre"))

# ================= 9. 单 SNP F 统计量 =================
cat("\n========== 9. F 统计量（弱工具检验）==========\n")
dat_mr$F_stat <- dat_mr$beta.exposure^2 / dat_mr$se.exposure^2
cat("F：最小值 =", min(dat_mr$F_stat), "，中位数 =", median(dat_mr$F_stat), "\n")

# ================= 10. 出图（手册 6.5，基于 dat_mr） =================
dir.create("D:/gwas/plots", showWarnings = FALSE)
png("D:/gwas/plots/scatter.png", width = 1800, height = 1200, res = 200)
print(mr_scatter_plot(res_main, dat_mr)); dev.off()
png("D:/gwas/plots/forest.png", width = 1800, height = 1200, res = 200)
print(mr_forest_plot(mr_singlesnp(dat_mr))); dev.off()
png("D:/gwas/plots/funnel.png", width = 1800, height = 1200, res = 200)
print(mr_funnel_plot(mr_singlesnp(dat_mr))); dev.off()
png("D:/gwas/plots/leaveoneout.png", width = 1800, height = 2400, res = 200)
print(mr_leaveoneout_plot(res_loo)); dev.off()
cat("✓ 4 张图已保存到 D:/gwas/plots/\n")

# ================= 11. MR-PRESSO（手册 6.6，10000 次置换，基于 dat_mr = 69 SNP） =================
cat("\n========== 11. MR-PRESSO（10000 次置换，约 10-30 分钟）==========\n")
cat("使用 SNP 数：", nrow(dat_mr), "（与主分析一致）\n")
mp <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
                data = dat_mr, NbDistribution = 10000, SignifThreshold = 0.05)
cat("\n--- Global Test ---\n")
print(mp$`MR-PRESSO results`$`Global Test`)
cat("\n--- Outlier Test ---\n")
print(mp$`MR-PRESSO results`$`Outlier Test`)
cat("\n--- Distortion Test ---\n")
print(mp$`MR-PRESSO results`$`Distortion Test`)

# ================= 12. 剔除离群 SNP 后 IVW（基于 dat_mr） =================
cat("\n========== 12. 剔除离群 SNP 后 IVW ==========\n")
ot <- mp$`MR-PRESSO results`$`Outlier Test`
if (!is.null(ot) && length(ot$Index) > 0) {
  cat("离群 SNP 索引：", ot$Index, "\n")
  cat("离群 SNP：", dat_mr$SNP[ot$Index], "\n")
  dat_clean <- dat_mr[-ot$Index, ]
} else {
  cat("未检出离群 SNP\n")
  dat_clean <- dat_mr
}
mr(dat_clean, method_list = c("mr_ivw_mre", "mr_ivw_fe"))

# ================= 13. 复制结局（手册 6.7，同样过滤 mr_keep） =================
cat("\n========== 13.1 复制结局 1：ebi-a-GCST90013867 ==========\n")
out1 <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST90013867")
dat1 <- harmonise_data(exp, out1)
dat1 <- dat1[dat1$mr_keep == TRUE, ]
cat("协调后有效 SNP 数：", nrow(dat1), "\n")
print(generate_odds_ratios(mr(dat1)))

cat("\n========== 13.2 复制结局 2：finn-b-C3_THYROID_GLAND ==========\n")
out2 <- extract_outcome_data(snps = exp$SNP, outcomes = "finn-b-C3_THYROID_GLAND")
dat2 <- harmonise_data(exp, out2)
dat2 <- dat2[dat2$mr_keep == TRUE, ]
cat("协调后有效 SNP 数：", nrow(dat2), "\n")
print(generate_odds_ratios(mr(dat2)))

# ================= 14. 完成 =================
cat("\n========== 脚本 A 完成 ==========\n")
cat("请把 Console 全部输出 + D:/gwas/plots/ 4 张图发回。\n")
