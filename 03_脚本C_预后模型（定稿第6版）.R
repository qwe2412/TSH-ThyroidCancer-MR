# =====================================================================
# 脚本 C 定稿（第 6 版）· 预后风险模型全流程重跑【建模集=金标准 409 例】
# ------------------------------------------------------------
# 【已锁定事实（勿改动）】
#   ① riskScore = expr_wide$symbol 第12行GNG7/第14行NFIA + 15位条形码列名
#      + 论文系数 −0.0053257×GNG7 − 0.0003788×NFIA → 与金标准 cor=1 差=0
#   ② 建模集 = 表达矩阵内 409 例；DFS 与金标准 409/409 一致
#   ③ 中位数分组 204/205，与金标准一致
#   ④ C-index 点估计 0.666、Log-rank P=0.0412、时间ROC用 KM 法（论文口径）
# 【本版修复（2026-09-13，v6）】
#   ① 校准：改用 basehaz(centered=TRUE) + predict(type="lp") 手算 36 月生存概率
#      （v5 的 survfit summary 对多曲线返回 matrix，赋给 data.frame 列报
#       "replacement has 1 row, data has 409" → 中断；手算方案无此问题）
#      并打印总体 KM 36 月生存作 sanity check；输出 3 组/5 组 MAE
#   ② Bootstrap optimism-corrected 修复方向 bug：
#      v5 用 concordance(Surv~lp) 公式接口，方向约定与 coxph 相反
#      → 假 optimism=0.331、假 optimism-C=0.335。
#      本版改用 concordance(fb, newdata=dat)（coxph 接口，方向一致）。
#      预期真实 optimism 很小、optimism-C ≈ 0.663 → 锁定论文 Bootstrap 口径
#   ③ 历史产物（cindex_opt/cindex_boot/cox_final/time_roc/dca_res）加载
#      提前到 Bootstrap 之后立即执行，确保拿到论文 CI 的真实口径
#   ④ DCA 在修复后的预测概率上重算
# 运行方式：Ctrl+A 全选 → Run → 把全部输出发回
# =====================================================================

library(survival); library(Hmisc); library(survivalROC); library(rms)

# ---- 0. 数据加载 ----
stopifnot(file.exists("D:/gwas/tcga_thca.RData"))
load("D:/gwas/tcga_thca.RData")
env_rm <- new.env(); load("D:/gwas/risk_model.RData", envir = env_rm)
gold <- env_rm$dat
cat("✓ 已加载数据；金标准 dim：", dim(gold), "\n")

# ---- 1. riskScore 重建 + 建模集 409 ----
i_g7  <- which(expr_wide$symbol == "GNG7")
i_nfa <- which(expr_wide$symbol == "NFIA")
b1 <- -0.0053257; b2 <- -0.0003788
bc <- colnames(expr_wide)[-1]
rs <- as.numeric(b1 * expr_wide[i_g7, -1] + b2 * expr_wide[i_nfa, -1])
names(rs) <- bc
mA <- match(as.character(gold$sampleId), bc)
cat("\n===== ① riskScore 对账（409）=====\n")
cat("匹配：", sum(!is.na(mA)), "/409；cor：", cor(rs[mA], gold$riskScore),
    "；最大差：", max(abs(rs[mA] - gold$riskScore)), "\n")
stopifnot(sum(!is.na(mA)) == 409, all(rs[mA] == gold$riskScore))

# ---- 2. DFS 重建并对账 ----
mcl <- match(gold$patientId, clin_wide$patientId)
DFS_time_new <- as.numeric(clin_wide$DFS_MONTHS[mcl])
DFS_status_new <- clin_wide$DFS_STATUS[mcl]
parse_dfs <- function(x) {
  s <- toupper(trimws(as.character(x)))
  ifelse(grepl("^1|DECEASED|RECURR|PROGRESS|METAST|EVENT|^YES", s), 1,
  ifelse(grepl("^0|LIVING|DISEASEFREE|FREE|^NO", s), 0, NA_real_))
}
DFS_event_new <- parse_dfs(DFS_status_new)
cat("\n===== ② DFS 对账（409）=====\n")
cat("DFS_time 一致：", sum(DFS_time_new == gold$DFS_time), "/409；",
    "DFS_event 一致：", sum(DFS_event_new == gold$DFS_event), "/409\n")
stopifnot(all(DFS_time_new == gold$DFS_time), all(DFS_event_new == gold$DFS_event))
cat("事件 27、删失 382\n")

dat <- data.frame(
  patientId = gold$patientId, sampleId = gold$sampleId, riskScore = rs[mA],
  AGE = gold$AGE, DARK = gold$DARK, ETHNICITY = gold$ETHNICITY,
  FUSION_OTHER = gold$FUSION_OTHER, HISTOLOGICAL_TYPE = gold$HISTOLOGICAL_TYPE,
  HISTORY_RADIATION_EXPOSURE = gold$HISTORY_RADIATION_EXPOSURE,
  MACIS_SCORE = gold$MACIS_SCORE, MB_COV = gold$MB_COV,
  PATH_M_STAGE = gold$PATH_M_STAGE, PATH_N_STAGE = gold$PATH_N_STAGE,
  PATH_T_STAGE = gold$PATH_T_STAGE, RACE = gold$RACE,
  RISK_GROUP = gold$RISK_GROUP, SAMPLE_COUNT = gold$SAMPLE_COUNT,
  SEX = gold$SEX, TUMOR_STATUS = gold$TUMOR_STATUS,
  DFS_MONTHS = gold$DFS_MONTHS, DFS_STATUS = gold$DFS_STATUS,
  OS_MONTHS = gold$OS_MONTHS, OS_STATUS = gold$OS_STATUS,
  DFS_time = DFS_time_new, DFS_event = DFS_event_new,
  stringsAsFactors = FALSE
)
dat$riskGroup <- ifelse(dat$riskScore > median(dat$riskScore), "High", "Low")
cat("\n===== ③ 分组 =====\n")
print(table(dat$riskGroup))
cat("与金标准一致：", sum(dat$riskGroup == gold$riskGroup), "/409\n")

# ---- 3. Cox + C-index（多口径）----
cat("\n===== ④ Cox（n=409）=====\n")
fit_cox <- coxph(Surv(DFS_time, DFS_event) ~ riskScore, data = dat)
print(summary(fit_cox))
c_all <- concordance(fit_cox)
app_C <- c_all$concordance
cat("C-index（concordance）：", app_C, "\n")
cat("CI-a（concordance 正态）：", app_C - 1.96*sqrt(c_all$var), "–",
    app_C + 1.96*sqrt(c_all$var), "\n")
cat("  → 半宽对应 SE =", sqrt(c_all$var), "；论文 CI 0.553–0.772 半宽/1.96 =",
    (0.772 - 0.553) / 2 / 1.96, "（对比）\n")
rc <- rcorr.cens(dat$riskScore, Surv(dat$DFS_time, dat$DFS_event))
cat("rcorr.cens C Index（方向相反）：", rc["C Index"], "；SE：", rc["S.D."], "\n")

# ---- 4. Bootstrap 双口径（optimism 方向修复）----
cat("\n===== ⑤ Bootstrap C-index（1000 次，seed=123）=====\n")
set.seed(123)
n <- nrow(dat)
c_boot_in <- numeric(1000); opt <- numeric(1000)
for (i in 1:1000) {
  idx <- sample(n, n, replace = TRUE)
  fb <- coxph(Surv(DFS_time, DFS_event) ~ riskScore, data = dat[idx, ])
  c_boot_in[i] <- concordance(fb)$concordance
  # 【修复】coxph 接口带 newdata：在原数据上评估 bootstrap 模型，方向与 coxph 一致
  c_app_b <- concordance(fb, newdata = dat)$concordance
  opt[i] <- c_boot_in[i] - c_app_b
}
cat("口径A（bootstrap 样本内 C-index 均值）：", mean(c_boot_in), "\n")
cat("口径A 95% 分位：", quantile(c_boot_in, c(0.025, 0.975)), "\n")
cat("口径B（optimism-corrected = 表观 - mean(optimism)）：",
    app_C - mean(opt), "\n")
cat("  其中表观 C =", app_C, "；mean(optimism) =", mean(opt), "\n")

# ---- 5. 历史统计量产物加载（锁定论文 CI / Bootstrap 口径）----
cat("\n===== ⑥ 历史 cindex_opt / cindex_boot / cox_final / time_roc / dca_res =====\n")
for (f in c("cindex_opt.RData", "cindex_boot.RData", "cox_final.RData",
            "time_roc.RData", "dca_res.RData")) {
  fp <- file.path("D:/gwas", f)
  if (!file.exists(fp)) { cat(f, "：不存在\n"); next }
  env_h <- new.env()
  ok <- tryCatch({ load(fp, envir = env_h); TRUE },
                 error = function(e) { cat(f, "加载失败：", conditionMessage(e), "\n"); FALSE })
  if (!ok) next
  cat("---", f, "对象：", paste(ls(env_h), collapse = ", "), "---\n")
  for (nm in ls(env_h)) {
    obj <- get(nm, envir = env_h)
    if (is.numeric(obj)) {
      cat("  ", nm, "（数值 length", length(obj), "）：",
          paste(utils::head(as.vector(obj), 20), collapse = ", "), "\n")
    } else if (is.data.frame(obj) || is.matrix(obj)) {
      cat("  ", nm, "：dim", paste(dim(obj), collapse = "x"), "\n")
      print(utils::head(obj, 5))
    } else {
      cat("  ", nm, "：class", class(obj)[1], "\n")
    }
  }
}

# ---- 6. Log-rank + KM ----
cat("\n===== ⑦ Log-rank =====\n")
lr <- survdiff(Surv(DFS_time, DFS_event) ~ riskGroup, data = dat)
cat("chisq：", lr$chisq, "，P =", 1 - pchisq(lr$chisq, 1), "\n")
print(summary(survfit(Surv(DFS_time, DFS_event) ~ riskGroup, data = dat))$table)

# ---- 7. 时间 ROC（KM 为主，论文口径）----
cat("\n===== ⑧ 时间依赖 ROC =====\n")
for (t in c(12, 36, 60)) {
  roc_km  <- survivalROC(Stime = dat$DFS_time, status = dat$DFS_event,
                         marker = dat$riskScore, predict.time = t, method = "KM")
  roc_nne <- survivalROC(Stime = dat$DFS_time, status = dat$DFS_event,
                         marker = dat$riskScore, predict.time = t,
                         method = "NNE", span = 0.25 * n^( -0.20))
  cat(t, "个月：AUC(KM) =", roc_km$AUC, "；AUC(NNE) =", roc_nne$AUC, "\n")
}

# ---- 8. 3 年校准【v6：basehaz(centered=TRUE) + lp 手算】----
cat("\n===== ⑨ 3 年校准（v6 手算）=====\n")
horizon <- 36
cal_df <- data.frame(time = dat$DFS_time, event = dat$DFS_event, risk = dat$riskScore)
fit_cal <- coxph(Surv(time, event) ~ risk, data = cal_df)
h0c <- basehaz(fit_cal)                    # centered=TRUE（默认），与 predict(lp) 匹配
h0c_36 <- approx(h0c$time, h0c$hazard, xout = horizon, rule = 2)$y
lp <- predict(fit_cal, type = "lp")        # 居中线性预测
pred_surv36 <- exp(-h0c_36 * exp(lp))
pred_ev36 <- 1 - pred_surv36
cal_df$pred <- pred_ev36
cat("36 月预测事件率：范围", range(pred_ev36), "；均值", mean(pred_ev36), "\n")
km_all <- survfit(Surv(time, event) ~ 1, data = cal_df)
km36 <- 1 - summary(km_all, times = horizon)$surv
cat("sanity check：总体 KM 36 月事件率 =", km36, "（预测均值应与其接近）\n")

cal_mae <- function(k) {
  cal_df$grp <- cut(cal_df$pred,
                    breaks = quantile(cal_df$pred, probs = seq(0, 1, 1/k), na.rm = TRUE),
                    include.lowest = TRUE, labels = FALSE)
  obs <- sapply(sort(unique(cal_df$grp)), function(g) {
    dg <- cal_df[cal_df$grp == g, ]
    kf <- survfit(Surv(time, event) ~ 1, data = dg)
    ss <- summary(kf, times = horizon)
    if (length(ss$surv) == 0) { ss2 <- summary(kf); 1 - ss2$surv[length(ss2$surv)] }
    else 1 - ss$surv[length(ss$surv)]
  })
  prd <- sapply(sort(unique(cal_df$grp)), function(g) mean(cal_df$pred[cal_df$grp == g]))
  list(obs = obs, prd = prd, mae = mean(abs(prd - obs)))
}
r3 <- cal_mae(3); r5 <- cal_mae(5)
cat("3 组：预测 =", paste(round(r3$prd, 3), collapse = ", "), "；KM 实际 =",
    paste(round(r3$obs, 3), collapse = ", "), "；MAE =", r3$mae, "\n")
cat("5 组：预测 =", paste(round(r5$prd, 3), collapse = ", "), "；KM 实际 =",
    paste(round(r5$obs, 3), collapse = ", "), "；MAE =", r5$mae, "\n")

# ---- 9. DCA（修复后的预测概率）----
cat("\n===== ⑩ DCA（36 月 net benefit）=====\n")
thresholds <- seq(0.05, 0.30, by = 0.05)
dca_df <- data.frame(thr = thresholds, NB_model = NA, NB_all = NA, NB_none = 0)
for (j in seq_along(thresholds)) {
  pt <- thresholds[j]
  tp <- sum(pred_ev36 > pt & cal_df$event == 1)
  fp <- sum(pred_ev36 > pt & cal_df$event == 0)
  dca_df$NB_model[j] <- tp / n - fp / n * pt / (1 - pt)
  ev_rate <- mean(cal_df$event)
  dca_df$NB_all[j] <- ev_rate - (1 - ev_rate) * pt / (1 - pt)
}
print(dca_df)
cat("模型净获益优于全治/不治的阈值：",
    paste(thresholds[dca_df$NB_model > pmax(dca_df$NB_all, dca_df$NB_none)], collapse = ", "), "\n")

# ---- 10. 保存 ----
cat("\n===== ⑪ 保存 risk_model.RData =====\n")
if (!file.exists("D:/gwas/risk_model_backup.RData")) {
  file.copy("D:/gwas/risk_model.RData", "D:/gwas/risk_model_backup.RData")
}
dat$t_os <- as.numeric(gold$OS_MONTHS)
dat$riskGroup2 <- gold$riskGroup2; dat$riskGroup3 <- gold$riskGroup3
dat$riskGroup5 <- gold$riskGroup5
save(dat, file = "D:/gwas/risk_model.RData")
cat("✓ 已保存（dat", nrow(dat), "×", ncol(dat), "）\n")

cat("\n========== 脚本 C v6 完成：请把全部输出发回 ==========\n")
