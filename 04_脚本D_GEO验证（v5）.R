# =====================================================================
# 脚本 D v5 · GSE33630 梯度验证（历史口径 riskScore_ext 对账版）
# ------------------------------------------------------------
# 【v4 运行确证（2026-09-13，审查结论，勿改动）】
#   df_ext 第 4 列 = riskScore_ext（0.0166 尺度、正值），历史
#   geo_cohens_d.csv 的 mean（0.01663/0.01659/0.01657）与 n（11/49/45）
#   与 riskScore_ext 完全吻合 → 历史 GEO 分析直接用 riskScore_ext 列。
#   v4 用论文系数（-0.0053257×gng7 - 0.0003788×nfia）重算的负数 score
#   与历史口径不一致（Cohen's d：4.13 vs 历史 3.871 等）→ 论文数值
#   基准 = 历史 csv，须以 riskScore_ext 复算。
# 本版任务：
#   ① 用 df_ext$riskScore_ext（历史口径）重算 KW + 两两 Wilcoxon + Cohen's d
#   ② 与 geo_cohens_d.csv 自动逐项比对（组均值 / n / d）
#   ③ 回归反推 riskScore_ext 的线性公式（带/不带截距），与论文系数
#      对照 → 锁定历史 GEO 的 riskScore 公式（决定论文方法学措辞）
#   ④ 若 riskScore_ext 的 d 与历史 csv 完全一致 → GEO 验证闭环
# 运行方式：Ctrl+A 全选 → Run → 把全部输出发回
# =====================================================================

stopifnot("D:/gwas/geo/gse33630_ext.RData 不存在" =
            file.exists("D:/gwas/geo/gse33630_ext.RData"))

cat("\n========== GSE33630 梯度验证（v5 · riskScore_ext 历史口径）==========\n")

# ================= 0. 隔离加载 =================
env_g <- new.env()
load("D:/gwas/geo/gse33630_ext.RData", envir = env_g)
df_ext <- get("df_ext", envir = env_g)
cat("df_ext 列名：", paste(names(df_ext), collapse = ", "), "\n")
stopifnot(all(c("grp", "gng7", "nfia", "riskScore_ext") %in% names(df_ext)))

grp    <- as.character(df_ext$grp)
g7     <- as.numeric(df_ext$gng7)
nfa    <- as.numeric(df_ext$nfia)
rs_ext <- as.numeric(df_ext$riskScore_ext)

cat("riskScore_ext 范围：", paste(range(rs_ext, na.rm = TRUE), collapse = " ~ "),
    "；NA 数：", sum(is.na(rs_ext)), "\n")
cat("grp 分布：\n"); print(table(grp, useNA = "ifany"))

# ================= 1. 历史口径分析（riskScore_ext）=================
d <- data.frame(grp = grp, score = rs_ext, stringsAsFactors = FALSE)
d <- d[!is.na(d$score), ]
tab <- table(d$grp)
d <- d[d$grp %in% names(tab)[tab >= 2], ]
cat("\n--- 各组 riskScore_ext 分布 ---\n")
print(aggregate(score ~ grp, data = d, function(x)
  round(c(n = length(x), mean = mean(x), median = median(x),
          min = min(x), max = max(x)), 6)), row.names = FALSE)

cat("\n--- Kruskal-Wallis（riskScore_ext）---\n")
kw <- kruskal.test(score ~ grp, data = d)
cat("χ² =", kw$statistic, "，df =", kw$parameter,
    "，P =", format(kw$p.value, scientific = TRUE), "\n")

cohens_d <- function(a, b) {
  na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * var(a) + (nb - 1) * var(b)) / (na + nb - 2))
  (mean(a) - mean(b)) / sp
}
cat("\n--- 两两 Wilcoxon + Cohen's d（riskScore_ext）---\n")
grps <- unique(d$grp)
res_new <- data.frame(compare = character(), mean1 = numeric(),
                      mean2 = numeric(), cohens_d = numeric(),
                      n1 = integer(), n2 = integer(), w_p = character(),
                      stringsAsFactors = FALSE)
for (i in 1:(length(grps) - 1)) for (j in (i + 1):length(grps)) {
  a <- d$score[d$grp == grps[i]]
  b <- d$score[d$grp == grps[j]]
  w <- wilcox.test(a, b)
  dd <- cohens_d(a, b)
  cat(grps[i], "vs", grps[j],
      "：Wilcoxon P =", format(w$p.value, scientific = TRUE),
      "；Cohen's d =", round(dd, 2),
      "（mean:", round(mean(a), 6), "vs", round(mean(b), 6), "）\n")
  res_new <- rbind(res_new, data.frame(
    compare = paste(grps[i], "vs", grps[j]),
    mean1 = mean(a), mean2 = mean(b), cohens_d = dd,
    n1 = length(a), n2 = length(b),
    w_p = format(w$p.value, scientific = TRUE)))
}

# ================= 2. 与历史 geo_cohens_d.csv 自动对账 =================
cat("\n========== 与历史 geo_cohens_d.csv 自动对账 ==========\n")
stopifnot("无历史 geo_cohens_d.csv 可对账" =
            file.exists("D:/gwas/geo_cohens_d.csv"))
geo_hist <- read.csv("D:/gwas/geo_cohens_d.csv", stringsAsFactors = FALSE)
cat("历史 csv 列名：", paste(names(geo_hist), collapse = ", "), "；行数：", nrow(geo_hist), "\n")
print(geo_hist)

ok_d   <- all(abs(res_new$cohens_d - geo_hist$cohens_d) < 1e-6)
ok_mean <- all(abs(res_new$mean1 - geo_hist$mean1) < 1e-6 &
                 abs(res_new$mean2 - geo_hist$mean2) < 1e-6)
ok_n   <- all(res_new$n1 == geo_hist$n1 & res_new$n2 == geo_hist$n2)
cat("\n对账结果：\n")
cat("  Cohen's d 与历史一致：", ok_d,
    "（新算", paste(round(res_new$cohens_d, 3), collapse = ", "),
    "vs 历史", paste(round(geo_hist$cohens_d, 3), collapse = ", "), "）\n")
cat("  组均值与历史一致：", ok_mean, "\n")
cat("  样本量与历史一致：", ok_n, "\n")
if (ok_d && ok_mean && ok_n) {
  cat("\n✓✓ GEO 梯度验证闭环：用 riskScore_ext 列复算与历史结果完全一致，\n")
  cat("   论文 GEO 数值 = 历史 geo_cohens_d.csv（Cohen's d ",
      paste(round(geo_hist$cohens_d, 3), collapse = "/"), "）。\n")
} else {
  cat("\n⚠ 存在差异：请把上方两表发回，逐项核对差异来源（公式/样本）。\n")
}

# ================= 3. 反推 riskScore_ext 公式（论文方法学依据）=================
cat("\n========== 反推 riskScore_ext 线性公式 ==========\n")
fit1 <- lm(riskScore_ext ~ gng7 + nfia, data = df_ext)
fit0 <- lm(riskScore_ext ~ gng7 + nfia - 1, data = df_ext)
cat("带截距：riskScore_ext = ", format(coef(fit1)[1], digits = 10),
    " + ", format(coef(fit1)[2], digits = 10), "×gng7 + ",
    format(coef(fit1)[3], digits = 10), "×nfia", "（R² =",
    round(summary(fit1)$r.squared, 6), "）\n")
cat("无截距：riskScore_ext = ", format(coef(fit0)[1], digits = 10),
    "×gng7 + ", format(coef(fit0)[2], digits = 10), "×nfia", "（R² =",
    round(summary(fit0)$r.squared, 6), "）\n")
b1 <- -0.0053257; b2 <- -0.0003788
cat("论文系数：gng7 = -0.0053257，nfia = -0.0003788\n")
cat("无截距反推系数 / 论文系数 比值：gng7 =",
    format(coef(fit0)[1] / b1, digits = 6), "；nfia =",
    format(coef(fit0)[2] / b2, digits = 6), "\n")
cat("（两比值若接近 → riskScore_ext ≈ 论文公式的常数缩放；\n")
cat("  若相差明显 → 历史 GEO 用的是另一套系数，论文方法学需按反推公式如实交代）\n")
score_paper <- -0.0053257 * g7 - 0.0003788 * nfa
cat("riskScore_ext 与论文系数 score 的 Pearson cor：",
    cor(rs_ext, score_paper), "\n")

cat("\n========== 脚本 D v5 完成：请把全部输出发回 ==========\n")
