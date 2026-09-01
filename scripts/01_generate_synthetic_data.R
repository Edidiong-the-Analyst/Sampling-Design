## 01_generate_synthetic_data.R
##
## Generates a FULLY SYNTHETIC cluster-level panel dataset simulating a
## community-based contraceptive counselling programme evaluation across
## six Nigerian states, at baseline (pre) and endline (post).
##
## This is a portfolio/demonstration dataset only. It is not real programme
## data from any organisation. Structure and rough parameter values are
## loosely inspired by publicly known patterns in Nigerian mCPR (modern
## contraceptive prevalence rate) statistics, but all clusters, values and
## outcomes below are simulated.
##
## Output: ../data/cluster_panel_synthetic.csv

set.seed(2026)

states <- c("Lagos", "Kano", "Kaduna", "Rivers", "Anambra", "Oyo")
clusters_per_state <- 24            # e.g. enumeration areas / communities
n_states <- length(states)
n_clusters <- n_states * clusters_per_state

cluster_id   <- sprintf("CL-%03d", seq_len(n_clusters))
state        <- rep(states, each = clusters_per_state)

# Assign treatment (programme rollout) to half the clusters within each
# state, mimicking a stratified rollout so treatment/control are balanced
# within state from the start (relevant to the sampling design script).
treated <- unlist(lapply(seq_len(n_states), function(i) {
  rep(c(1, 0), each = clusters_per_state / 2)
}))

# Cluster-level baseline characteristics (fixed effects / confounders)
baseline_mCPR   <- pmin(pmax(rnorm(n_clusters, mean = 24, sd = 6), 5), 55)   # %
remoteness_index <- round(runif(n_clusters, 0, 10), 1)   # 0 = urban, 10 = very hard-to-reach
cluster_effect   <- rnorm(n_clusters, mean = 0, sd = 3)   # unobserved cluster quality

# TRUE programme effect used to simulate the data (the value our DiD
# analysis in script 03 should recover approximately, within noise):
TRUE_EFFECT_PP <- 8   # percentage-point increase in mCPR attributable to programme

build_period <- function(period_label, post, trend_step, apply_effect) {
  # secular trend common to both arms (contraceptive access improving nationally),
  # accumulated by trend_step periods since the first wave (2023 = 0)
  secular_trend <- 1.2 * trend_step
  # programme effect only applies to treated clusters once the programme has
  # actually launched (apply_effect == 1), i.e. the post period only
  programme_effect <- ifelse(treated == 1 & apply_effect == 1, TRUE_EFFECT_PP, 0)

  mCPR <- baseline_mCPR + secular_trend + programme_effect +
    cluster_effect + rnorm(n_clusters, 0, 3.5) - 0.3 * remoteness_index * apply_effect

  mCPR <- pmin(pmax(mCPR, 0), 100)

  # sample size per cluster survey round (varies a bit, as real fieldwork does)
  n_surveyed <- round(rnorm(n_clusters, mean = 180, sd = 25))
  n_surveyed <- pmax(n_surveyed, 60)
  n_using_modern_method <- round(n_surveyed * mCPR / 100)

  data.frame(
    cluster_id, state, treated,
    period = period_label, post,
    remoteness_index,
    n_surveyed, n_using_modern_method,
    mCPR_pct = round(100 * n_using_modern_method / n_surveyed, 2)
  )
}

# Three survey waves: two pre-programme waves (2023, 2024) let us check
# parallel pre-trends before trusting the DiD estimate off the third wave.
# Programme launches between the 2024 and 2025 waves.
panel <- rbind(
  build_period("2023_wave0_pre",  post = 0, trend_step = 0, apply_effect = 0),
  build_period("2024_wave1_pre",  post = 0, trend_step = 1, apply_effect = 0),
  build_period("2025_wave2_post", post = 1, trend_step = 2, apply_effect = 1)
)

panel <- panel[order(panel$cluster_id, panel$post), ]
rownames(panel) <- NULL

out_dir <- "../data"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
write.csv(panel, file.path(out_dir, "cluster_panel_synthetic.csv"), row.names = FALSE)

cat(sprintf(
  "Generated synthetic panel: %d clusters x 3 survey waves = %d rows\n",
  n_clusters, nrow(panel)
))
cat(sprintf("True simulated programme effect: +%.1f percentage points on mCPR\n", TRUE_EFFECT_PP))
cat("Saved to: ../data/cluster_panel_synthetic.csv\n")
