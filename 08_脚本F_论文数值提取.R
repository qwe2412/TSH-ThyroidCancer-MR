# =====================================================================
# 脚本 F 论文数值提取 · 从 F v3 预测结果提取论文 3.6/4.5 段所需全部数值
# ------------------------------------------------------------
# 目的：论文 F 段（药物敏感性）需更新为 F v3 真实预测数值。
#       本脚本从 D:/gwas/drug_sensitivity_res.RData（F v3 已保存）提取：
#        ① 显著药物总数/比例（FDR<0.05）
#        ② 显著药物中"High 组 IC50 更低(敏感)"vs"High 组更高(耐药)"的数量
#        ③ 论文 3.6/4.5 段引用的具体药物逐一输出（mean_High/mean_Low/Δ/P/FDR）
#        ④ 全部显著药物按 |Δ| 降序前 30（供论文选代表性药物）
#        ⑤ 全量结果存 csv 供引用
# 【耗时】< 1 分钟
# 【运行方式】Ctrl+A 全选 → Run → 把全部输出发回
# =====================================================================

# ================= 0. 加载 F v3 预测结果 =================
stopifnot(file.exists("D:/gwas/drug_sensitivity_res.RData"))
load("D:/gwas/drug_sensitivity_res.RData")   # pred_ic50(409×198), drug_res(198×6), grp
cat("pred_ic50 dim：", paste(dim(pred_ic50), collapse = " x "), "\n")
cat("drug_res 列名：", paste(colnames(drug_res), collapse = ", "), "\n")
cat("分组：High", sum(grp == "High"), "/ Low", sum(grp == "Low"), "\n")

hi <- grp == "High"; lo <- grp == "Low"

# ================= 1. 显著药物统计 =================
cat("\n========== 1. 显著药物统计（FDR < 0.05）==========\n")
sig <- drug_res$P_BH < 0.05
cat("FDR<0.05 药物数：", sum(sig, na.rm = TRUE), "/", nrow(drug_res),
    "（", round(100 * sum(sig, na.rm = TRUE) / nrow(drug_res), 1), "%）\n")

# ================= 2. 方向统计 =================
cat("\n========== 2. 显著药物中方向统计 ==========\n")
# Δ = mean_High - mean_Low；Δ<0 → High 组 IC50 更低（更敏感）
dlt <- drug_res$mean_High - drug_res$mean_Low
sig_dlt <- dlt[sig]
cat("显著药物中 High 组更敏感（mean_High < mean_Low）：",
    sum(sig_dlt < 0, na.rm = TRUE), "个\n")
cat("显著药物中 High 组更耐药（mean_High > mean_Low）：",
    sum(sig_dlt > 0, na.rm = TRUE), "个\n")
cat("→ 论文叙述方向判定依据（旧稿写“High 组对细胞毒药更敏感”，需按此核实）。\n")

# ================= 3. 论文引用的具体药物 =================
cat("\n========== 3. 论文 3.6/4.5 段引用药物逐一提取 ==========\n")
# 论文提及（GDSC 列名格式通常为"药物名_ID"，此处按关键词匹配）
paper_drugs <- c(
  "Docetaxel", "Paclitaxel", "Camptothecin", "Cisplatin", "Vinblastine",
  "Fluorouracil", "Doxorubicin", "YK-4-279", "KU-55933", "Daporinad",
  "GSK591", "VE821", "Cytarabine", "Sorafenib", "Nilotinib"
)
dr_cols <- colnames(pred_ic50)
cat("GDSC 列名示例（前 10）：", paste(head(dr_cols, 10), collapse = ", "), "\n")
for (kw in paper_drugs) {
  hit <- dr_cols[grepl(kw, dr_cols, ignore.case = TRUE)]
  if (length(hit) == 0) {
    cat(sprintf("  %-14s：未匹配到列\n", kw))
    next
  }
  for (h in hit) {
    r <- drug_res[drug_res$drug == h, ]
    if (nrow(r) == 1) {
      cat(sprintf("  %-14s -> %-22s | High=%.4f Low=%.4f Δ=%.4f | P=%.4g FDR=%.4g%s\n",
                  kw, h, r$mean_High, r$mean_Low, r$mean_High - r$mean_Low,
                  r$P, r$P_BH, ifelse(is.na(r$P_BH) || r$P_BH >= 0.05, "  [NS]", "")))
    }
  }
}

# ================= 4. 显著药物 |Δ| 前 30 =================
cat("\n========== 4. 显著药物按 |Δ| 降序前 30（供论文选代表性药物）==========\n")
tmp <- drug_res[sig, ]
tmp$abs_delta <- abs(tmp$mean_High - tmp$mean_Low)
tmp <- tmp[order(-tmp$abs_delta), ]
tmp$direction <- ifelse(tmp$mean_High < tmp$mean_Low, "High敏感", "High耐药")
print(head(tmp[, c("drug", "mean_High", "mean_Low", "P", "P_BH", "direction")], 30),
      row.names = FALSE, digits = 5)

# ================= 5. 全量结果存 csv =================
cat("\n========== 5. 保存 ==========\n")
drug_res_out <- drug_res
drug_res_out$mean_High_minus_Low <- dlt
drug_res_out$direction <- ifelse(dlt < 0, "High敏感", "High耐药")
drug_res_out$significant <- sig
write.csv(drug_res_out, "D:/gwas/drug_sensitivity_res_full.csv", row.names = FALSE)
cat("✓ 已保存 D:/gwas/drug_sensitivity_res_full.csv（", nrow(drug_res_out),
    " 行 × ", ncol(drug_res_out), " 列）\n", sep = "")

cat("\n========== 脚本 F 论文数值提取完成：请把全部输出发回 ==========\n")
