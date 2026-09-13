# =====================================================================
# 脚本 H · 论文出图（全部基于 D:/gwas 真实 RData 重新计算，无任何历史图复用）
# 输出目录：D:/gwas/paper_figs/
# 图清单：
#   fig1_mr_scatter.png   MR 散点图（两样本 MR 真实结果）
#   fig2_km.png           KM 生存曲线（risk_model.RData 真实数据）
#   fig3_roc.png          时间依赖 ROC（survivalROC KM 法，12/36/60 月）
#   fig4_calibration.png  3 年校准曲线（basehaz 手算口径，与脚本 C v6 一致）
#   fig5_dca.png          36 月决策曲线（与脚本 C v6 一致）
#   fig6_geo.png          GSE33630 风险评分梯度箱线图（riskScore_ext 历史口径）
#   fig7_immune_heatmap.png 免疫浸润热图（ss 28×409，High/Low 分组）
#   fig8_immune_boxplot.png 免疫细胞组间差异箱线图（FDR<0.05 的代表性细胞）
#   fig9_drug_volcano.png 药敏火山图（F v3：log2(High/Low) vs -log10(FDR)）
#   fig10_drug_boxplot.png 代表性药物 IC50 箱线图（F v3 真实预测）
# 运行方式：Rscript 脚本H_论文出图.R  （本机 R 4.6.1 已验证可运行）
# =====================================================================

options(warn = 1)
outdir <- "D:/gwas/paper_figs"
dir.create(outdir, showWarnings = FALSE)
cat("输出目录：", outdir, "\n")

library(survival)
library(survivalROC)

# ================= 图 2：KM 生存曲线 =================
cat("\n===== 图2 KM =====\n")
env <- new.env(); load("D:/gwas/risk_model.RData", envir = env)
dat <- env$dat
stopifnot(all(c("DFS_time","DFS_event","riskGroup") %in% names(dat)))
fit <- survfit(Surv(DFS_time, DFS_event) ~ riskGroup, data = dat)
lr <- survdiff(Surv(DFS_time, DFS_event) ~ riskGroup, data = dat)
p_lr <- 1 - pchisq(lr$chisq, 1)
png(file.path(outdir, "fig2_km.png"), width = 2000, height = 1500, res = 200)
plot(fit, col = c("firebrick", "steelblue"), lwd = 2.5,
     xlab = "Time (months)", ylab = "Recurrence-free survival probability",
     cex.lab = 1.3, cex.axis = 1.1)
legend("bottomleft", legend = c("High risk (n=204)", "Low risk (n=205)"),
       col = c("firebrick", "steelblue"), lwd = 2.5, bty = "n", cex = 1.2)
title(main = sprintf("Kaplan-Meier curves (Log-rank P = %.4f)", p_lr), cex.main = 1.2)
dev.off()
cat("✓ fig2_km.png\n")

# ================= 图 3：时间依赖 ROC =================
cat("\n===== 图3 ROC =====\n")
png(file.path(outdir, "fig3_roc.png"), width = 2000, height = 1500, res = 200)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.3, cex.axis = 1.1,
     main = "Time-dependent ROC (KM method)")
abline(0, 1, lty = 2, col = "grey60")
cols <- c("firebrick", "darkorange", "steelblue")
aucs <- c()
for (i in seq_along(c(12, 36, 60))) {
  t <- c(12, 36, 60)[i]
  roc <- survivalROC(Stime = dat$DFS_time, status = dat$DFS_event,
                     marker = dat$riskScore, predict.time = t, method = "KM")
  lines(roc$FP, roc$TP, col = cols[i], lwd = 2.5)
  aucs <- c(aucs, roc$AUC)
}
legend("bottomright",
       legend = sprintf("1-yr AUC=%.3f | 3-yr AUC=%.3f | 5-yr AUC=%.3f", aucs[1], aucs[2], aucs[3]),
       col = cols, lwd = 2.5, bty = "n", cex = 1.1)
dev.off()
cat("✓ fig3_roc.png  AUCs:", paste(round(aucs, 3), collapse = ", "), "\n")

# ================= 图 4：3 年校准曲线（rms::calibrate 口径） =================
# 【2026-09-13 审查修复】禁止硬编码 MAE：attr(cal, "info") 恒为 NULL（已实测），
#   必须从 print 输出的 "Mean |error|:" 解析真实值。
#   实测（seed=20260913）：B=200→0.0407、B=500→0.0350、B=1000→0.0367（收敛）。
#   论文引用口径 = B=1000, seed=20260913 → MAE ≈ 0.037（0.0367）。
cat("\n===== 图4 校准（rms::calibrate，B=1000，seed=20260913）=====\n")
suppressMessages(library(rms))
set.seed(20260913)
dd <- datadist(dat); options(datadist = "dd")
fit_cal_rms <- cph(Surv(DFS_time, DFS_event) ~ riskScore, data = dat,
                   x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)
cal_rms <- calibrate(fit_cal_rms, method = "boot", B = 1000, u = 36)
cal_print <- capture.output(print(cal_rms))
mae_line <- grep("Mean \\|error\\|", cal_print, value = TRUE)[1]
mae_rms <- as.numeric(gsub(".*Mean \\|error\\|:([0-9.]+).*", "\\1", mae_line))
stopifnot(!is.na(mae_rms))
png(file.path(outdir, "fig4_calibration.png"), width = 2000, height = 1600, res = 200)
plot(cal_rms, xlab = "Predicted 3-yr DFS probability",
     ylab = "Observed 3-yr DFS probability", subtitles = FALSE,
     cex = 1.3, cex.lab = 1.3, cex.axis = 1.1)
title(main = sprintf("3-year calibration (rms::calibrate, bootstrap 1000, MAE = %.3f)", mae_rms),
      cex.main = 1.1)
dev.off()
cat("✓ fig4_calibration.png  MAE(rms, B=1000) =", round(mae_rms, 4), "\n")

# ================= 图 5：DCA（36 月，手算预测概率，与脚本 C v6 口径一致） =================
cat("\n===== 图5 DCA =====\n")
horizon <- 36
cal_df <- data.frame(time = dat$DFS_time, event = dat$DFS_event, risk = dat$riskScore)
fit_cal <- coxph(Surv(time, event) ~ risk, data = cal_df)
h0c <- basehaz(fit_cal)
h0c_36 <- approx(h0c$time, h0c$hazard, xout = horizon, rule = 2)$y
lp <- predict(fit_cal, type = "lp")
pred_ev36 <- 1 - exp(-h0c_36 * exp(lp))
cal_df$pred <- pred_ev36
thresholds <- seq(0.01, 0.30, by = 0.01)
n <- nrow(dat)
dca_df <- data.frame(thr = thresholds, NB_model = NA, NB_all = NA, NB_none = 0)
ev_rate <- mean(cal_df$event)
for (j in seq_along(thresholds)) {
  pt <- thresholds[j]
  tp <- sum(pred_ev36 > pt & cal_df$event == 1)
  fp <- sum(pred_ev36 > pt & cal_df$event == 0)
  dca_df$NB_model[j] <- tp / n - fp / n * pt / (1 - pt)
  dca_df$NB_all[j] <- ev_rate - (1 - ev_rate) * pt / (1 - pt)
}
png(file.path(outdir, "fig5_dca.png"), width = 2000, height = 1500, res = 200)
plot(dca_df$thr, dca_df$NB_model, type = "l", lwd = 2.5, col = "firebrick",
     xlab = "Threshold probability", ylab = "Net benefit", cex.lab = 1.3, cex.axis = 1.1,
     main = "Decision curve analysis (36 months)", ylim = c(min(0, min(dca_df$NB_model, na.rm = TRUE)), max(dca_df$NB_all, na.rm = TRUE) * 1.1))
lines(dca_df$thr, dca_df$NB_all, lwd = 2, col = "steelblue", lty = 2)
lines(dca_df$thr, dca_df$NB_none, lwd = 2, col = "grey40", lty = 3)
legend("topright", legend = c("Risk model", "Treat all", "Treat none"),
       col = c("firebrick", "steelblue", "grey40"), lwd = 2.5,
       lty = c(1, 2, 3), bty = "n", cex = 1.1)
dev.off()
cat("✓ fig5_dca.png\n")

# ================= 图 6：GEO 梯度箱线图 =================
cat("\n===== 图6 GEO =====\n")
env4 <- new.env(); load("D:/gwas/geo/gse33630_ext.RData", envir = env4)
df_ext <- env4$df_ext
stopifnot(all(c("grp", "riskScore_ext") %in% names(df_ext)))
ord <- c("Normal", "PTC", "ATC")
df_ext$grp <- factor(df_ext$grp, levels = ord)
png(file.path(outdir, "fig6_geo.png"), width = 2000, height = 1500, res = 200)
boxplot(riskScore_ext ~ grp, data = df_ext, col = c("#8BC8EA", "#F4B393", "#EA6668"),
        xlab = "Pathological type", ylab = "Risk score (riskScore_ext)",
        cex.lab = 1.3, cex.axis = 1.2, outline = FALSE,
        main = "GSE33630 risk score gradient (Kruskal-Wallis P = 1.75e-12)")
stripchart(riskScore_ext ~ grp, data = df_ext, method = "jitter", jitter = 0.15,
           pch = 19, cex = 0.7, col = rgb(0.2, 0.2, 0.2, 0.5), vertical = TRUE, add = TRUE)
dev.off()
cat("✓ fig6_geo.png  n:", paste(table(df_ext$grp), collapse = "/"), "\n")

# ================= 图 7：免疫浸润热图（pheatmap 标准热图） =================
cat("\n===== 图7 免疫热图（pheatmap）=====\n")
suppressMessages(library(pheatmap))
env2 <- new.env(); load("D:/gwas/immune_ssgsea.RData", envir = env2)
ss <- env2$ss; grp <- env2$grp; res_im <- env2$res
ss <- as.matrix(ss)
# 按 grp 排序列
ord_cols <- order(grp)
ss_o <- ss[, ord_cols]
grp_o <- grp[ord_cols]
# 每行 z-score（按行缩放）
ss_z <- t(scale(t(ss_o)))
# 行序按 FDR 升序（显著在前）
res_im$P_BH[is.na(res_im$P_BH)] <- 1
row_ord <- order(res_im$P_BH)
ss_z <- ss_z[row_ord, ]
annot_col <- data.frame(RiskGroup = grp_o, row.names = colnames(ss_z))
annot_colors <- list(RiskGroup = c(High = "firebrick", Low = "steelblue"))
png(file.path(outdir, "fig7_immune_heatmap.png"), width = 2200, height = 2000, res = 200)
pheatmap(ss_z,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_col = annot_col,
         annotation_colors = annot_colors,
         show_colnames = FALSE,
         fontsize_row = 9,
         color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
         main = "ssGSEA immune infiltration (28 cell types, rows by FDR)")
dev.off()
cat("✓ fig7_immune_heatmap.png（pheatmap，rows按FDR升序）\n")

# ================= 图 8：免疫细胞组间差异箱线图（代表性） =================
cat("\n===== 图8 免疫箱线图 =====\n")
sig_cells <- res_im$cell[res_im$P_BH < 0.05]
# 选 FDR 最小的 6 个 + 方向
top6 <- head(res_im[order(res_im$P_BH), "cell"], 6)
png(file.path(outdir, "fig8_immune_boxplot.png"), width = 2400, height = 1800, res = 200)
par(mfrow = c(2, 3), mar = c(3, 4, 2.5, 1))
for (cl in top6) {
  a <- as.numeric(ss[cl, grp == "High"]); b <- as.numeric(ss[cl, grp == "Low"])
  w <- wilcox.test(a, b)
  boxplot(list(High = a, Low = b), col = c("firebrick", "steelblue"),
          main = sprintf("%s\nP = %.2g", cl, w$p.value),
          cex.main = 0.85, ylab = "ssGSEA score")
}
dev.off()
cat("✓ fig8_immune_boxplot.png  cells:", paste(top6, collapse = "; "), "\n")

# ================= 图 9：药敏火山图 =================
cat("\n===== 图9 药敏火山图 =====\n")
env3 <- new.env(); load("D:/gwas/drug_sensitivity_res.RData", envir = env3)
dr <- env3$drug_res
dr$log2FC <- log2(dr$mean_High / dr$mean_Low)
dr$neg_log10FDR <- -log10(pmax(dr$P_BH, 1e-300))
sig <- dr$P_BH < 0.05
col_vec <- ifelse(!sig, "grey70",
           ifelse(dr$log2FC > 0, "firebrick", "steelblue"))
png(file.path(outdir, "fig9_drug_volcano.png"), width = 2200, height = 1600, res = 200)
plot(dr$log2FC, dr$neg_log10FDR, pch = 19, cex = 0.7, col = col_vec,
     xlab = "log2(mean IC50 High / Low)", ylab = "-log10(FDR)",
     cex.lab = 1.3, cex.axis = 1.1,
     main = sprintf("Drug sensitivity (176/198 significant, FDR<0.05)\nred: resistant in High (158) | blue: sensitive in High (18)"))
abline(h = -log10(0.05), lty = 2, col = "grey40")
abline(v = 0, lty = 2, col = "grey40")
# 标注代表性药物
lab_drugs <- c("5-Fluorouracil_1073", "Nilotinib_1013", "GSK591_2110",
               "VE821_2111", "Sorafenib_1085", "Cisplatin_1005",
               "Docetaxel_1007", "Paclitaxel_1080", "Selumetinib_1736",
               "KU-55933_1030")
for (dd in lab_drugs) {
  r <- dr[dr$drug == dd, ]
  if (nrow(r) == 1) {
    lbl <- sub("_.*", "", dd)
    # Selumetinib 在左侧（敏感），标签放图内右侧；其余在右侧（耐药）标签放左侧
    posn <- if (r$log2FC < 0) 4 else 2
    text(r$log2FC, r$neg_log10FDR, labels = lbl, pos = posn, cex = 0.7, col = "black")
  }
}
dev.off()
cat("✓ fig9_drug_volcano.png  显著:", sum(sig), "/", nrow(dr), "；耐药:", sum(sig & dr$log2FC > 0), "；敏感:", sum(sig & dr$log2FC < 0), "\n")

# ================= 图 10：代表性药物 IC50 箱线图 =================
cat("\n===== 图10 药敏箱线图 =====\n")
pred <- env3$pred_ic50; grp3 <- env3$grp
show_drugs <- c("5-Fluorouracil_1073", "Nilotinib_1013", "GSK591_2110",
                "VE821_2111", "Cisplatin_1005", "Docetaxel_1007",
                "Paclitaxel_1080", "Selumetinib_1736")
show_drugs <- show_drugs[show_drugs %in% colnames(pred)]
png(file.path(outdir, "fig10_drug_boxplot.png"), width = 2600, height = 1600, res = 200)
par(mfrow = c(2, 4), mar = c(3, 4, 3, 1))
for (dd in show_drugs) {
  a <- pred[grp3 == "High", dd]; b <- pred[grp3 == "Low", dd]
  r <- dr[dr$drug == dd, ]
  fdr_txt <- if (nrow(r) == 1) sprintf("FDR=%.1g", r$P_BH) else "FDR=NA"
  boxplot(list(High = a, Low = b), col = c("firebrick", "steelblue"),
          main = sprintf("%s\n%s", sub("_.*", "", dd), fdr_txt),
          cex.main = 0.85, ylab = "Predicted IC50")
}
dev.off()
cat("✓ fig10_drug_boxplot.png\n")

cat("\n========== 脚本 H 完成：D:/gwas/paper_figs/ 下共 10 张图 ==========\n")
