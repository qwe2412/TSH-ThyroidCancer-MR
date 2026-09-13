# =====================================================================
# 全流程重跑 · 脚本 B：MV-MR 多变量孟德尔随机化（对应手册 5.4）【v3 重写】
# ------------------------------------------------------------
# 【方法学说明（2026-09-13 重写）】
#   MVMR 的 SNP 集 = TSH 的 70 个工具 SNP（与主分析一致），
#   对每个 SNP 从 BMI（ieu-b-40）/ 睾酮（ebi-a-GCST90014013）/ 结局
#   （ebi-a-GCST90018929）GWAS 提取效应，构建四组效应矩阵后做 MVMR-IVW。
#   ⚠ 不要用"三暴露工具并集"：TSH 暴露（Zhou 2020）只有 99 个位点有
#     效应数据，BMI/睾酮的工具 SNP 大多没有 TSH 效应 → 交集必为 0。
# ------------------------------------------------------------
# 运行方式：RStudio 中本文件 → Ctrl+A 全选 → 点 Run
# 前置：脚本 A 已在本会话运行过（exp、out 存在）
# 耗时：BMI/睾酮/结局效应提取约 1-2 分钟（在线）
# =====================================================================

library(TwoSampleMR)
library(MendelianRandomization)

# ================= 0. 检查脚本 A 的环境 =================
stopifnot("未检测到 exp（TSH 工具变量）。请先运行脚本 A，再运行本脚本。" = exists("exp"))
stopifnot("未检测到 out（结局提取）。请先运行脚本 A。" = exists("out"))
cat("✓ 复用脚本 A 的 TSH 工具变量（", nrow(exp), " SNP）\n")

# ================= 1. 在 TSH 工具 SNP 上提取 BMI 效应 =================
cat("\n========== 1. 提取 BMI 效应（ieu-b-40）==========\n")
bmi_eff <- extract_outcome_data(snps = exp$SNP, outcomes = "ieu-b-40")
cat("BMI 效应返回行数：", nrow(bmi_eff), "\n")

# ================= 2. 在 TSH 工具 SNP 上提取睾酮效应 =================
cat("\n========== 2. 提取睾酮效应（ebi-a-GCST90014013）==========\n")
tt_eff <- extract_outcome_data(snps = exp$SNP, outcomes = "ebi-a-GCST90014013")
cat("睾酮效应返回行数：", nrow(tt_eff), "\n")

# ================= 3. 三组协调（均与 TSH 暴露对齐，过滤 mr_keep） =================
cat("\n========== 3. 协调（结局/BMI/睾酮）==========\n")
d_out <- harmonise_data(exp, out);          d_out <- d_out[d_out$mr_keep == TRUE, ]
d_bmi <- harmonise_data(exp, bmi_eff);      d_bmi <- d_bmi[d_bmi$mr_keep == TRUE, ]
d_tt  <- harmonise_data(exp, tt_eff);       d_tt  <- d_tt[d_tt$mr_keep == TRUE, ]
cat("结局协调有效 SNP：", nrow(d_out), "；BMI：", nrow(d_bmi), "；睾酮：", nrow(d_tt), "\n")

# ================= 4. 按 SNP 对齐四组效应 =================
cat("\n========== 4. 四组效应对齐 ==========\n")
keep <- Reduce(intersect, list(d_out$SNP, d_bmi$SNP, d_tt$SNP))
cat("四组效应均存在、纳入 MVMR 的 SNP 数：", length(keep), "\n")
m <- data.frame(
  SNP    = keep,
  b_tsh  = d_out$beta.exposure[match(keep, d_out$SNP)],
  se_tsh = d_out$se.exposure[match(keep, d_out$SNP)],
  b_bmi  = d_bmi$beta.outcome[match(keep, d_bmi$SNP)],
  se_bmi = d_bmi$se.outcome[match(keep, d_bmi$SNP)],
  b_tt   = d_tt$beta.outcome[match(keep, d_tt$SNP)],
  se_tt  = d_tt$se.outcome[match(keep, d_tt$SNP)],
  b_out  = d_out$beta.outcome[match(keep, d_out$SNP)],
  se_out = d_out$se.outcome[match(keep, d_out$SNP)])

# ================= 5. MVMR-IVW（MendelianRandomization） =================
cat("\n========== 5. MVMR-IVW 结果 ==========\n")
mvmr <- mr_mvinput(bx   = as.matrix(m[, c("b_tsh", "b_bmi", "b_tt")]),
                   bxse = as.matrix(m[, c("se_tsh", "se_bmi", "se_tt")]),
                   by   = m$b_out, byse = m$se_out)
res_mv <- mr_mvivw(mvmr)
print(res_mv)
cat("\n--- TSH 的 MVMR 效应（校正 BMI 与睾酮）---\n")
cat("b =", res_mv@Estimate[1], "，SE =", res_mv@StdError[1], "，P =", res_mv@Pvalue[1], "\n")
cat("OR =", exp(res_mv@Estimate[1]),
    "（95% CI", exp(res_mv@Estimate[1] - 1.96 * res_mv@StdError[1]),
    "–", exp(res_mv@Estimate[1] + 1.96 * res_mv@StdError[1]), "）\n")

# ================= 6. 完成 =================
cat("\n========== 脚本 B 完成 ==========\n")
cat("请把第 4 步 SNP 数与第 5 步 MVMR-IVW 结果发回。\n")
