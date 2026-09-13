# =====================================================================
# 脚本 E v4 · 免疫浸润分析（修复版）
# 对应手册 11.2-11.4：ssGSEA（28 免疫细胞，Charoentong 2017）
#                     → 组间差异（Wilcoxon + BH）→ 免疫检查点差异
# ------------------------------------------------------------
# 【v4 修复（2026-09-13，严格审查后）】
#  ① GSVA 2.6 的 gsva() 是单参数 S4 泛型：正确调用 = gsva(ssgseaParam(...))，
#     v3 写成 gsva(expr, ssgseaParam(...)) 会把 expr(matrix) 当 param 分派，
#     报"函数'gsva'标签'param = "matrix"'找不到继承方法"。
#     → 本版改为 gsva(ssgseaParam(expr, gene_sets)) 单参数形式。
# 【v3 修复（2026-09-13）】
#  ① 检测 GSVA 是否导出 ssgseaParam 分流新旧 API。
# 【v2 修复（2026-09-13）】
#  ① dat_final → dat；② 分组对齐以 dat 顺序为准；③ gmt 路径自适应探测；
#  ④ 历史产物保护 + 隔离加载对账。
# 【前置】GSVA / GSEABase 已装；D:/gwas/Charoentong2017.gmt（已确认存在）
# 【运行方式】Ctrl+A 全选 → Run → 把全部输出发回
# =====================================================================

# ================= 0. 包与数据加载 =================
if (!requireNamespace("GSVA", quietly = TRUE)) {
  stop("GSVA 未安装。请先运行 BiocManager::install('GSVA') 后重试。")
}
if (!requireNamespace("GSEABase", quietly = TRUE)) {
  stop("GSEABase 未安装。请先运行 BiocManager::install('GSEABase') 后重试。")
}
library(GSVA)
library(GSEABase)

stopifnot("expr_immune.RData 不存在" = file.exists("D:/gwas/expr_immune.RData"))
stopifnot("risk_model.RData 不存在" = file.exists("D:/gwas/risk_model.RData"))
load("D:/gwas/expr_immune.RData")     # expr_mat_use（样本 × 基因，行名=12位患者ID）
load("D:/gwas/risk_model.RData")      # dat（409×30，含 riskGroup）
stopifnot("expr_immune.RData 中无 expr_mat_use" = exists("expr_mat_use"))
stopifnot("risk_model.RData 中无 dat（应为 dat，不是 dat_final）" = exists("dat"))
cat("✓ expr_mat_use 维度：", dim(expr_mat_use),
    "；dat 维度：", dim(dat), "\n")

# ================= 1. 分组对齐（修复：以 dat 顺序为准）=================
cat("\n========== 1. 样本对齐与分组 ==========\n")
common <- intersect(rownames(expr_mat_use), as.character(dat$patientId))
cat("共同患者：", length(common), "/", nrow(dat), "\n")
stopifnot("共同患者 < 100，请检查患者 ID 格式" = length(common) >= 100)

idx <- match(as.character(dat$patientId), rownames(expr_mat_use))
idx <- idx[!is.na(idx)]
expr_aligned <- expr_mat_use[idx, ]
grp <- dat$riskGroup[match(rownames(expr_aligned), as.character(dat$patientId))]
stopifnot("riskGroup 匹配失败" = all(!is.na(grp)))
cat("样本数：", nrow(expr_aligned), "；分组分布：\n")
print(table(grp, useNA = "ifany"))
cat("（High/Low 应为 204/205，与金标准一致）\n")

# ================= 2. gmt 路径探测 =================
cat("\n========== 2. Charoentong2017.gmt 探测 ==========\n")
gmt_candidates <- c(
  "D:/gwas/Charoentong2017.gmt",
  "D:/gwas/gmt/Charoentong2017.gmt",
  "D:/gwas/immune/Charoentong2017.gmt",
  file.path(getwd(), "Charoentong2017.gmt")
)
found <- gmt_candidates[file.exists(gmt_candidates)]
if (length(found) == 0) {
  cat("⚠ 未在以下位置找到 Charoentong2017.gmt：\n")
  for (p in gmt_candidates) cat("   ", p, "\n")
  stop("请告知 gmt 文件实际路径，或将其复制到 D:/gwas/ 后重跑。")
}
gmt_path <- found[1]
cat("✓ 使用 gmt：", gmt_path, "\n")

# ================= 3. 解析 gmt（MSigDB 格式：名称\t描述\t基因...）=================
lines <- readLines(gmt_path, warn = FALSE)
lines <- lines[nzchar(trimws(lines))]
gs_split <- strsplit(lines, "\t")
gene_sets <- lapply(gs_split, function(x) x[-(1:2)])
names(gene_sets) <- sapply(gs_split, function(x) x[1])
cat("gmt 基因集总数：", length(gene_sets), "\n")
cat("基因集名称（前 10）：", paste(head(names(gene_sets), 10), collapse = ", "), "\n")

# 与表达矩阵列名（基因符号）取交集，保留重叠 ≥3 的基因集
gene_sets <- lapply(gene_sets, function(genes) intersect(genes, colnames(expr_aligned)))
keep <- sapply(gene_sets, length) >= 3
gene_sets <- gene_sets[keep]
cat("与表达矩阵重叠 ≥3 基因的基因集数：", length(gene_sets), "\n")
stopifnot("有效基因集 < 5，请检查基因符号格式（大小写）" = length(gene_sets) >= 5)

# ================= 4. ssGSEA（v4：单参数 param 调用）=================
cat("\n========== 3. ssGSEA（GSVA 包）==========\n")
expr_gene_x_sample <- t(expr_aligned)          # 基因 × 样本
cat("GSVA 版本：", as.character(packageVersion("GSVA")), "\n")
cat("GSVA 是否导出 ssgseaParam：", "ssgseaParam" %in% getNamespaceExports("GSVA"), "\n")
if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
  # GSVA >= 1.50：gsva(param) 单参数泛型（2.x 起不再接受 gsva(expr, sets, method=)）
  cat("→ 使用新版 API：gsva(ssgseaParam(expr, gene_sets))\n")
  ss <- GSVA::gsva(GSVA::ssgseaParam(expr_gene_x_sample, gene_sets),
                   verbose = FALSE)
} else {
  # 旧版 GSVA（< 1.50）：method="ssgsea"
  cat("→ 使用旧版 API：gsva(expr, gene_sets, method='ssgsea')\n")
  ss <- GSVA::gsva(expr_gene_x_sample, gene_sets,
                   method = "ssgsea", verbose = FALSE)
}
cat("ssGSEA 得分矩阵维度（基因集 × 样本）：", dim(ss), "\n")
if (is(ss, "SummarizedExperiment")) ss <- assay(ss)
if (!is.matrix(ss)) ss <- as.matrix(ss)

# ================= 5. 组间差异（Wilcoxon + BH）=================
cat("\n========== 4. 免疫细胞浸润组间差异（High vs Low）==========\n")
hi <- grp == "High"; lo <- grp == "Low"
res <- data.frame(cell = rownames(ss),
                  mean_High = NA_real_, mean_Low = NA_real_,
                  W = NA_real_, P = NA_real_, P_BH = NA_real_,
                  stringsAsFactors = FALSE)
for (i in seq_len(nrow(ss))) {
  a <- as.numeric(ss[i, hi]); b <- as.numeric(ss[i, lo])
  w <- wilcox.test(a, b)
  res$mean_High[i] <- mean(a); res$mean_Low[i] <- mean(b)
  res$W[i] <- w$statistic; res$P[i] <- w$p.value
}
res$P_BH <- p.adjust(res$P, method = "BH")
print(res, row.names = FALSE, digits = 5)
cat("\nP_BH < 0.05 的细胞类型：",
    paste(res$cell[res$P_BH < 0.05], collapse = ", "), "\n")

# ================= 6. 免疫检查点 / 共抑制标志物差异 =================
cat("\n========== 5. 免疫检查点基因差异（High vs Low）==========\n")
checkpoint_genes <- c("CD274", "PDCD1", "PDCD1LG2", "CTLA4", "LAG3",
                      "HAVCR2", "TIGIT", "IDO1", "SIGLEC15", "CD80",
                      "CD86", "TNFRSF9", "IFNG", "GZMB", "PRF1")
ck <- intersect(checkpoint_genes, colnames(expr_aligned))
cat("表达矩阵中可用的检查点基因（", length(ck), " 个）：",
    paste(ck, collapse = ", "), "\n")
if (length(ck) > 0) {
  ck_res <- data.frame(gene = ck, mean_High = NA_real_, mean_Low = NA_real_,
                       P = NA_real_, P_BH = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(ck)) {
    g <- ck[i]
    a <- as.numeric(expr_aligned[hi, g]); b <- as.numeric(expr_aligned[lo, g])
    w <- wilcox.test(a, b)
    ck_res$mean_High[i] <- mean(a); ck_res$mean_Low[i] <- mean(b)
    ck_res$P[i] <- w$p.value
  }
  ck_res$P_BH <- p.adjust(ck_res$P, method = "BH")
  print(ck_res, row.names = FALSE, digits = 5)
}

# ================= 7. 保存（备份保护）=================
cat("\n========== 6. 保存 immune_ssgsea.RData ==========\n")
if (file.exists("D:/gwas/immune_ssgsea.RData")) {
  if (!file.exists("D:/gwas/immune_ssgsea_backup.RData")) {
    file.copy("D:/gwas/immune_ssgsea.RData", "D:/gwas/immune_ssgsea_backup.RData")
    cat("✓ 历史 immune_ssgsea.RData 已备份为 immune_ssgsea_backup.RData\n")
  }
}
save(ss, res, grp, expr_aligned, file = "D:/gwas/immune_ssgsea.RData")
cat("✓ 已保存 D:/gwas/immune_ssgsea.RData（ss：", dim(ss)[1], "×", dim(ss)[2], "）\n")

# ================= 8. 与历史 immune_ssgsea_backup.RData 对账 =================
if (file.exists("D:/gwas/immune_ssgsea_backup.RData")) {
  cat("\n========== 7. 历史 immune_ssgsea_backup.RData 对账 ==========\n")
  env_old <- new.env()
  ok <- tryCatch({ load("D:/gwas/immune_ssgsea_backup.RData", envir = env_old); TRUE },
                 error = function(e) { cat("历史产物加载失败：", conditionMessage(e), "\n"); FALSE })
  if (ok) {
    cat("历史对象：", paste(ls(env_old), collapse = ", "), "\n")
    for (nm in ls(env_old)) {
      obj <- get(nm, envir = env_old)
      if (is.matrix(obj) || is.data.frame(obj)) {
        cat("  ", nm, ": dim =", paste(dim(obj), collapse = " x "),
            "；rownames 前 5：", paste(utils::head(rownames(obj), 5), collapse = ", "), "\n")
      } else if (is.numeric(obj) && length(obj) <= 100) {
        cat("  ", nm, ": length", length(obj), "；", paste(utils::head(obj, 10), collapse = ", "), "\n")
      } else {
        cat("  ", nm, ": class", class(obj)[1], "\n")
      }
    }
    cat("→ 请核对：本次 ssGSEA 得分与历史是否同源（行名/样本集一致、数值可比）。\n")
  }
}

cat("\n========== 脚本 E v4 完成：请把全部输出发回 ==========\n")
