# 重复 5x5 折 CV C-index：手写 Harrell's C + coxph iter.max=100
suppressMessages({library(survival)})
env <- new.env(); load("D:/gwas/lasso_cv.RData", envir = env)
dat <- env$dat
genes <- setdiff(colnames(dat), c("patient", "sample", "DFS_time", "DFS_event",
                                  "DFS_MONTHS", "DFS_STATUS", "OS_MONTHS", "OS_STATUS"))
X <- as.matrix(dat[, genes]); X[is.na(X)] <- median(X, na.rm = TRUE)
cat("事件数:", sum(dat$DFS_event), "样本:", nrow(dat), "\n")

harrell_c <- function(time, status, lp) {
  n <- length(time); ev <- which(status == 1)
  ct <- 0; cc <- 0
  for (i in ev) {
    later <- which(time > time[i])
    if (length(later) == 0) next
    ct <- ct + length(later)
    cc <- cc + sum(lp[later] < lp[i])
  }
  if (ct == 0) NA else cc / ct
}

rep_cindex <- numeric(25); k <- 0
set.seed(20260914)
for (r in 1:5) {
  folds <- sample(rep(1:5, length.out = nrow(dat)))
  for (f in 1:5) {
    te <- folds == f
    fit <- tryCatch(
      coxph(Surv(dat$DFS_time[!te], dat$DFS_event[!te]) ~ X[!te, "GNG7"] + X[!te, "NFIA"],
            control = coxph.control(iter.max = 100)),
      error = function(e) NULL)
    if (is.null(fit)) { k <- k + 1; rep_cindex[k] <- NA; next }
    lp <- predict(fit, newdata = data.frame(GNG7 = X[te, "GNG7"], NFIA = X[te, "NFIA"]))
    k <- k + 1
    if (sum(dat$DFS_event[te]) == 0 || length(unique(lp)) < 2) { rep_cindex[k] <- NA; next }
    rep_cindex[k] <- harrell_c(dat$DFS_time[te], dat$DFS_event[te], lp)
  }
}
rep_cindex <- rep_cindex[1:k]
cat("\n重复 5x5 折 CV C-index：均值 =", round(mean(rep_cindex, na.rm = TRUE), 4),
    "；范围 =", paste(round(range(rep_cindex, na.rm = TRUE), 4), collapse = " ~ "),
    "；有效折 =", sum(!is.na(rep_cindex)), "/", k, "\n")
cat("各折 C-index:", paste(round(rep_cindex, 3), collapse = ", "), "\n")
