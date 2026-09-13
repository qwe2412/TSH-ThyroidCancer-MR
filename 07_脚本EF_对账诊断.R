# =====================================================================
# 脚本 EF 对账诊断 · E 免疫浸润 & F 药敏 历史数值对账（严格审查）
# ------------------------------------------------------------
# 目的：E v4 / F v3 已产出新数值。本脚本与历史产物逐项对账，
#       判定"可复现"或"口径差异/真错误"，供后续决策。
# 【E 部分】本次 ss(28×409)/res  vs 历史 immune_ssgsea_backup.RData
#           (ssgsea_score / res_immune) → 逐细胞 cor + 数值差
# 【F 部分】本次 pred_ic50(409×198) vs 历史 ic50_pred（tcga_thca 内嵌）
#           → 三套 cor（原始 Pearson / log 后 Pearson / Spearman）
#           + 历史异常值统计（判断 1.8e32 等爆炸值占比）
#           + 历史 drug_sensitivity_res_backup.RData 结构
# 【耗时】< 2 分钟
# 【运行方式】Ctrl+A 全选 → Run → 把全部输出发回
# =====================================================================

# =====================================================================
# E 部分：本次 ss/res vs 历史 ssgsea_score/res_immune
# =====================================================================
cat("\n############ E 免疫浸润对账 ############\n")
stopifnot(file.exists("D:/gwas/immune_ssgsea.RData"))
load("D:/gwas/immune_ssgsea.RData")          # 本次：ss(28×409), res(28×6), grp

if (file.exists("D:/gwas/immune_ssgsea_backup.RData")) {
  env_hist <- new.env()
  ok <- tryCatch({ load("D:/gwas/immune_ssgsea_backup.RData", envir = env_hist); TRUE },
                 error = function(e) { cat("历史加载失败：", conditionMessage(e), "\n"); FALSE })
  if (ok) {
    cat("历史对象：", paste(ls(env_hist), collapse = ", "), "\n")
    hist_ss  <- if (exists("ssgsea_score", envir = env_hist)) get("ssgsea_score", envir = env_hist) else NULL
    hist_res <- if (exists("res_immune", envir = env_hist)) get("res_immune", envir = env_hist) else NULL

    # ---- 1a. ss 得分矩阵对账 ----
    if (!is.null(hist_ss)) {
      hist_ss <- as.matrix(hist_ss)
      cat("本次 ss dim：", paste(dim(ss), collapse = " x "),
          "；历史 ssgsea_score dim：", paste(dim(hist_ss), collapse = " x "), "\n")
      ccells <- intersect(rownames(ss), rownames(hist_ss))
      csamp  <- intersect(colnames(ss), colnames(hist_ss))
      cat("共同细胞类型：", length(ccells), "/", nrow(ss),
          "；共同样本：", length(csamp), "/", ncol(ss), "\n")
      if (length(ccells) >= 10 && length(csamp) >= 100) {
        ss_new  <- ss[ccells, csamp]
        ss_hist <- hist_ss[ccells, csamp]
        cors <- sapply(ccells, function(cl) cor(as.numeric(ss_new[cl, ]),
                                                as.numeric(ss_hist[cl, ])))
        cat("\n逐细胞类型 cor（本次 vs 历史）：\n")
        print(round(cors, 4))
        cat("→ cor 中位数：", round(median(cors), 4),
            "；最小值：", round(min(cors), 4), "\n")
        diffs <- as.numeric(ss_new) - as.numeric(ss_hist)
        cat("→ 数值差绝对最大值：", round(max(abs(diffs)), 6),
            "；|diff|>0.01 占比：", round(mean(abs(diffs) > 0.01), 4), "\n")
      }
    }

    # ---- 1b. 组间差异表对账 ----
    if (!is.null(hist_res)) {
      cat("\n历史 res_immune dim：", paste(dim(hist_res), collapse = " x "),
          "；列名：", paste(colnames(hist_res), collapse = ", "), "\n")
      cat("历史 res_immune 全表：\n")
      print(as.data.frame(hist_res), digits = 5)
      cat("\n本次 res（前 6 列，cell/mean_High/mean_Low/W/P/P_BH）：\n")
      print(res, digits = 5)
      cat("→ 请人工核对：历史与本次的 P 值列是否同源、方向是否一致。\n")
    }
  }
} else {
  cat("⚠ 未找到 immune_ssgsea_backup.RData，跳过 E 对账。\n")
}

# =====================================================================
# F 部分：本次 pred_ic50 vs 历史 ic50_pred + 历史药敏备份结构
# =====================================================================
cat("\n############ F 药敏对账 ############\n")
stopifnot(file.exists("D:/gwas/drug_sensitivity_res.RData"))
load("D:/gwas/drug_sensitivity_res.RData")   # 本次：pred_ic50(409×198), drug_res(198×6), grp

if (exists("pred_ic50") && is.matrix(pred_ic50)) {
  cat("本次 pred_ic50 dim：", paste(dim(pred_ic50), collapse = " x "),
      "；范围：", paste(round(range(pred_ic50, na.rm = TRUE), 4), collapse = " ~ "),
      "；负值数：", sum(pred_ic50 < 0, na.rm = TRUE), "\n")
}

# ---- 2a. 与 tcga_thca 内嵌历史 ic50_pred 对账 ----
if (file.exists("D:/gwas/tcga_thca.RData")) {
  env_t <- new.env()
  load("D:/gwas/tcga_thca.RData", envir = env_t)
  if (exists("ic50_pred", envir = env_t)) {
    hist_pred <- as.matrix(get("ic50_pred", envir = env_t))
    cat("\n历史 ic50_pred（tcga_thca 内嵌）dim：", paste(dim(hist_pred), collapse = " x "), "\n")
    cat("历史范围：", paste(round(range(hist_pred, na.rm = TRUE), 4), collapse = " ~ "),
        "；负值数：", sum(hist_pred < 0, na.rm = TRUE),
        "；NA 数：", sum(is.na(hist_pred)), "\n")
    cat("历史 |x|>100 占比：", round(mean(abs(hist_pred) > 100, na.rm = TRUE), 4),
        "；|x|>1e6 占比：", round(mean(abs(hist_pred) > 1e6, na.rm = TRUE), 6), "\n")

    cdr <- intersect(colnames(pred_ic50), colnames(hist_pred))
    cpa <- intersect(rownames(pred_ic50), rownames(hist_pred))
    cat("共同药物：", length(cdr), "/", ncol(hist_pred),
        "；共同患者：", length(cpa), "/", nrow(hist_pred), "\n")

    if (length(cdr) > 10 && length(cpa) > 100) {
      # ① 原始 Pearson（尺度可能不匹配）
      p1 <- sapply(cdr, function(d) cor(pred_ic50[cpa, d], hist_pred[cpa, d], use = "complete.obs"))
      # ② 若本次为正尺度：log(本次) vs 历史（历史若为 log 尺度则应吻合）
      p2 <- tryCatch(sapply(cdr, function(d) cor(log(pred_ic50[cpa, d]), hist_pred[cpa, d], use = "complete.obs")),
                     error = function(e) rep(NA, length(cdr)))
      # ③ Spearman（不受单调变换影响，最能反映相对排序一致性）
      p3 <- sapply(cdr, function(d) cor(pred_ic50[cpa, d], hist_pred[cpa, d],
                                        method = "spearman", use = "complete.obs"))
      cat("\n三套逐药物 cor 汇总（中位数 / 范围 / >0.9 数 / >0.5 数）：\n")
      for (lab in c("原始Pearson", "log(本次)Pearson", "Spearman")) {
        v <- switch(lab, "原始Pearson" = p1, "log(本次)Pearson" = p2, "Spearman" = p3)
        cat(sprintf("  %-14s：中位 %.4f | 范围 %.4f ~ %.4f | >0.9: %d | >0.5: %d\n",
                    lab, median(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE),
                    sum(v > 0.9, na.rm = TRUE), sum(v > 0.5, na.rm = TRUE)))
      }
      cat("\n→ 判定参考：若 Spearman 高而 Pearson 低 = 纯尺度差异（变换一致）；
         若 Spearman 也低 = 预测流程/数据版本差异，需进一步排查。\n")
    }
  } else {
    cat("⚠ tcga_thca.RData 中无 ic50_pred。\n")
  }
}

# ---- 2b. 历史 drug_sensitivity_res_backup.RData 结构 ----
if (file.exists("D:/gwas/drug_sensitivity_res_backup.RData")) {
  cat("\n历史 drug_sensitivity_res_backup.RData 结构：\n")
  env_hf <- new.env()
  okf <- tryCatch({ load("D:/gwas/drug_sensitivity_res_backup.RData", envir = env_hf); TRUE },
                  error = function(e) { cat("历史加载失败：", conditionMessage(e), "\n"); FALSE })
  if (okf) {
    cat("历史对象：", paste(ls(env_hf), collapse = ", "), "\n")
    for (nm in ls(env_hf)) {
      obj <- get(nm, envir = env_hf)
      if (is.matrix(obj) || is.data.frame(obj)) {
        nums <- tryCatch(as.numeric(as.matrix(obj)), error = function(e) NULL)
        cat("  ", nm, ": dim", paste(dim(obj), collapse = " x "),
            "；行名前3", paste(head(rownames(obj), 3), collapse = ", "),
            "；列名前3", paste(head(colnames(obj), 3), collapse = ", "),
            if (!is.null(nums)) paste("；范围", paste(round(range(nums, na.rm = TRUE), 3), collapse = " ~ ")),
            "\n")
      } else {
        cat("  ", nm, ": class", class(obj)[1], "\n")
      }
    }
  }
} else {
  cat("⚠ 未找到 drug_sensitivity_res_backup.RData（历史药敏备份）。\n")
}

cat("\n========== 脚本 EF 对账诊断完成：请把全部输出发回 ==========\n")
