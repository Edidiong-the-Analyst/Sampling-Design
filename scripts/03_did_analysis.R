## 03_did_analysis.R
##
## Difference-in-Differences analysis of the (synthetic) contraceptive
## counselling programme's effect on modern contraceptive prevalence (mCPR).
##
## Demonstrates:
##   1. A standard two-period DiD regression (2024 pre vs 2025 post)
##   2. Cluster-robust standard errors computed manually (CR1 sandwich
##      estimator) -- no external packages required
##   3. A pre-trend / parallel-trends check using the extra 2023 wave
##   4. A permutation (randomization-inference) placebo test as a
##      robustness check that does not rely on asymptotic assumptions
##
## Written in base R only.

panel <- read.csv("../data/cluster_panel_synthetic.csv", stringsAsFactors = FALSE)

## ---- Manual cluster-robust SE helper (CR1 sandwich estimator) ----
cluster_robust_se <- function(model, cluster) {
  X <- model.matrix(model)
  u <- residuals(model)
  cluster <- factor(cluster)
  M <- nlevels(cluster)
  N <- nrow(X)
  K <- ncol(X)

  meat <- matrix(0, K, K)
  for (g in levels(cluster)) {
    idx <- which(cluster == g)
    Xg <- X[idx, , drop = FALSE]
    ug <- u[idx]
    meat <- meat + t(Xg) %*% (ug %*% t(ug)) %*% Xg
  }
  bread <- solve(t(X) %*% X)
  dfc <- (M / (M - 1)) * ((N - 1) / (N - K))
  vcov_cr <- dfc * bread %*% meat %*% bread
  se <- sqrt(diag(vcov_cr))
  se
}

report_did <- function(model, cluster, label) {
  se <- cluster_robust_se(model, cluster)
  coefs <- coef(model)
  common <- intersect(names(coefs), names(se))
  tvals <- coefs[common] / se[common]
  pvals <- 2 * pt(-abs(tvals), df = model$df.residual)
  out <- data.frame(
    term = common,
    estimate = round(coefs[common], 3),
    cluster_robust_se = round(se[common], 3),
    t_value = round(tvals, 2),
    p_value = signif(pvals, 3)
  )
  cat(sprintf("\n=== %s ===\n", label))
  print(out, row.names = FALSE)
  out
}

## ---- 1. Two-period DiD (2024 pre vs 2025 post) ----
two_period <- panel[panel$period %in% c("2024_wave1_pre", "2025_wave2_post"), ]
two_period$post <- ifelse(two_period$period == "2025_wave2_post", 1, 0)

did_model <- lm(mCPR_pct ~ treated * post, data = two_period,
                 weights = n_surveyed)  # weight by survey size, larger clusters more precise

did_results <- report_did(did_model, two_period$cluster_id, "DiD Regression: mCPR (%) ~ treated * post (weighted by n_surveyed)")

cat(sprintf(
  "\nInterpretation: the 'treated:post' coefficient is the estimated causal effect of the\nprogramme on modern contraceptive prevalence, in percentage points, net of the secular trend.\n"
))

## ---- 2. Pre-trend / parallel-trends check ----
# Compare the treated-vs-control gap in the two PRE periods only (2023 -> 2024).
# If the programme assignment is as-if random with respect to trend, this
# interaction should NOT be statistically distinguishable from zero.
pre_only <- panel[panel$period %in% c("2023_wave0_pre", "2024_wave1_pre"), ]
pre_only$wave1 <- ifelse(pre_only$period == "2024_wave1_pre", 1, 0)

pretrend_model <- lm(mCPR_pct ~ treated * wave1, data = pre_only, weights = n_surveyed)
pretrend_results <- report_did(pretrend_model, pre_only$cluster_id,
                                "Pre-trend check: mCPR (%) ~ treated * wave1, PRE-PROGRAMME PERIODS ONLY")

pretrend_p <- pretrend_results$p_value[pretrend_results$term == "treated:wave1"]
cat(sprintf(
  "\nPre-trend check result: p = %s for treated x wave1 interaction (pre-programme).\n%s\n",
  pretrend_p,
  ifelse(pretrend_p > 0.10,
         "No evidence of differential pre-trends (p > 0.10) -- supports the parallel-trends assumption.",
         "WARNING: evidence of differential pre-trends -- the DiD estimate above should be interpreted with caution.")
))

## ---- 3. Permutation / randomization-inference placebo test ----
# Re-randomise which clusters are "treated" (within state, preserving the
# true allocation ratio) many times, re-estimate the DiD coefficient each
# time, and see how extreme the REAL estimate is relative to this placebo
# distribution. This is a robustness check that does not rely on the
# normal-theory assumptions behind the cluster-robust SE above.
set.seed(2026)
n_perm <- 1000
real_effect <- did_results$estimate[did_results$term == "treated:post"]

perm_effects <- numeric(n_perm)
states <- unique(two_period$state)
base_clusters <- unique(two_period[, c("cluster_id", "state")])

for (p in seq_len(n_perm)) {
  perm_treat <- unlist(lapply(states, function(s) {
    ids <- base_clusters$cluster_id[base_clusters$state == s]
    n_t <- sum(two_period$treated[two_period$state == s & two_period$post == 0] == 1) / 1 # per period count already 1:1
    n_t <- length(ids) / 2
    sample(c(rep(1, n_t), rep(0, length(ids) - n_t)))
  }))
  names(perm_treat) <- unlist(lapply(states, function(s) base_clusters$cluster_id[base_clusters$state == s]))

  tmp <- two_period
  tmp$perm_treated <- perm_treat[tmp$cluster_id]
  m <- lm(mCPR_pct ~ perm_treated * post, data = tmp, weights = n_surveyed)
  cf <- coef(m)
  perm_effects[p] <- if ("perm_treated:post" %in% names(cf)) cf["perm_treated:post"] else NA
}

perm_effects <- perm_effects[!is.na(perm_effects)]
perm_p_value <- mean(abs(perm_effects) >= abs(real_effect))

cat(sprintf(
  "\n=== Permutation placebo test (%d re-randomisations) ===\nReal estimated effect: %.2f pp\nPermutation-based p-value: %.3f\n",
  length(perm_effects), real_effect, perm_p_value
))

## ---- Save outputs ----
if (!dir.exists("../output")) dir.create("../output", recursive = TRUE)
write.csv(did_results, "../output/did_regression_results.csv", row.names = FALSE)
write.csv(pretrend_results, "../output/pretrend_check_results.csv", row.names = FALSE)

png("../output/pretrend_and_effect_plot.png", width = 900, height = 600, res = 120)
agg <- aggregate(mCPR_pct ~ period + treated, data = panel, FUN = mean)
agg$period_num <- match(agg$period, c("2023_wave0_pre", "2024_wave1_pre", "2025_wave2_post"))
plot(NA, xlim = c(1, 3), ylim = range(agg$mCPR_pct) + c(-2, 2),
     xaxt = "n", xlab = "Survey wave", ylab = "Mean mCPR (%)",
     main = "Treated vs. control mCPR over time (synthetic data)")
axis(1, at = 1:3, labels = c("2023 (pre)", "2024 (pre)", "2025 (post)"))
abline(v = 2.5, lty = 2, col = "grey60")
text(2.5, max(agg$mCPR_pct) + 1.5, "Programme launch", col = "grey40", cex = 0.8, pos = 4)
lines(agg$period_num[agg$treated == 1], agg$mCPR_pct[agg$treated == 1], type = "o", col = "#0B6E4F", lwd = 2, pch = 16)
lines(agg$period_num[agg$treated == 0], agg$mCPR_pct[agg$treated == 0], type = "o", col = "#888888", lwd = 2, pch = 16)
legend("topleft", legend = c("Treated clusters", "Control clusters"),
       col = c("#0B6E4F", "#888888"), lwd = 2, pch = 16, bty = "n")
dev.off()

sink("../output/did_analysis_summary.txt")
cat("Difference-in-Differences Analysis Summary\n===========================================\n\n")
cat("1. Main DiD estimate (2024 pre vs 2025 post):\n")
print(did_results, row.names = FALSE)
cat(sprintf("\nEstimated programme effect: %.2f percentage points (cluster-robust SE reported above)\n", real_effect))
cat("\n2. Pre-trend check (2023 vs 2024, pre-programme):\n")
print(pretrend_results, row.names = FALSE)
cat(sprintf("\np-value on treated x wave1 (pre-trend): %s\n", pretrend_p))
cat(sprintf("\n3. Permutation placebo test: p = %.3f (%d re-randomisations)\n", perm_p_value, length(perm_effects)))
sink()

cat("\nSaved: ../output/did_regression_results.csv\n")
cat("Saved: ../output/pretrend_check_results.csv\n")
cat("Saved: ../output/pretrend_and_effect_plot.png\n")
cat("Saved: ../output/did_analysis_summary.txt\n")
