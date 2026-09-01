## 02_sampling_design.R
##
## Sampling strategy for a cluster-level evaluation of a contraceptive
## counselling programme's effect on modern contraceptive prevalence
## (mCPR). Demonstrates: minimum-detectable-effect-driven sample size
## calculation for a two-arm cluster design, design-effect adjustment for
## intra-cluster correlation, and proportional stratified allocation
## across states.
##
## Written in base R only (no external packages required) so it runs
## anywhere R is installed.

## ---- 1. Parameters (would normally come from formative research / prior data) ----
p_control  <- 0.25   # assumed baseline mCPR in control clusters (25%)
mde_pp     <- 0.08   # minimum detectable effect: 8 percentage points
p_treat    <- p_control + mde_pp
alpha      <- 0.05   # two-sided significance level
power      <- 0.80
icc        <- 0.02   # assumed intra-cluster correlation for mCPR at community level
m_per_cluster <- 180 # planned respondents surveyed per cluster (from prior fieldwork norms)

z_alpha <- qnorm(1 - alpha / 2)
z_power <- qnorm(power)

## ---- 2. Individual-level sample size per arm (standard two-proportion test) ----
p_bar <- (p_control + p_treat) / 2
n_individual_per_arm <- ((z_alpha * sqrt(2 * p_bar * (1 - p_bar)) +
                             z_power * sqrt(p_control * (1 - p_control) + p_treat * (1 - p_treat)))^2) /
  (mde_pp^2)

## ---- 3. Design effect adjustment for clustering ----
design_effect <- 1 + (m_per_cluster - 1) * icc
n_effective_per_arm <- n_individual_per_arm * design_effect

## ---- 4. Convert to number of clusters needed per arm ----
clusters_per_arm <- ceiling(n_effective_per_arm / m_per_cluster)
total_clusters   <- clusters_per_arm * 2
total_respondents <- total_clusters * m_per_cluster

cat("=== Sample Size / Power Calculation ===\n")
cat(sprintf("Baseline mCPR (control):        %.0f%%\n", p_control * 100))
cat(sprintf("Minimum detectable effect:      %.0f percentage points\n", mde_pp * 100))
cat(sprintf("Significance level (alpha):     %.2f (two-sided)\n", alpha))
cat(sprintf("Power:                           %.2f\n", power))
cat(sprintf("Assumed ICC:                     %.3f\n", icc))
cat(sprintf("Respondents per cluster:         %d\n", m_per_cluster))
cat(sprintf("Design effect:                   %.2f\n", design_effect))
cat(sprintf("Required respondents per arm:    %d (before clustering)\n", ceiling(n_individual_per_arm)))
cat(sprintf("Required respondents per arm:    %d (after design effect)\n", ceiling(n_effective_per_arm)))
cat(sprintf("Required clusters per arm:       %d\n", clusters_per_arm))
cat(sprintf("Total clusters (both arms):      %d\n", total_clusters))
cat(sprintf("Total planned respondents:       %d\n\n", total_respondents))

## ---- 5. Proportional stratified allocation across states ----
# Allocate the required clusters proportionally to each state's share of
# eligible communities, then round while preserving the total (largest
# remainder method) -- a standard, defensible allocation approach.
states <- c("Lagos", "Kano", "Kaduna", "Rivers", "Anambra", "Oyo")
eligible_communities <- c(Lagos = 340, Kano = 410, Kaduna = 300,
                           Rivers = 260, Anambra = 190, Oyo = 250)  # illustrative sampling frame sizes

raw_alloc <- total_clusters * eligible_communities / sum(eligible_communities)
floor_alloc <- floor(raw_alloc)
remainder <- total_clusters - sum(floor_alloc)
remainders <- raw_alloc - floor_alloc
add_to <- order(remainders, decreasing = TRUE)[seq_len(remainder)]
final_alloc <- floor_alloc
final_alloc[add_to] <- final_alloc[add_to] + 1

alloc_table <- data.frame(
  state = names(eligible_communities),
  eligible_communities = as.integer(eligible_communities),
  clusters_allocated = as.integer(final_alloc),
  clusters_treatment = as.integer(round(final_alloc / 2)),
  clusters_control = as.integer(final_alloc - round(final_alloc / 2))
)

cat("=== Stratified Allocation by State (proportional to sampling frame, largest-remainder rounding) ===\n")
print(alloc_table, row.names = FALSE)

if (!dir.exists("../output")) dir.create("../output", recursive = TRUE)
write.csv(alloc_table, "../output/sampling_allocation_by_state.csv", row.names = FALSE)

sink("../output/sampling_design_summary.txt")
cat("Sampling Design Summary\n========================\n\n")
cat(sprintf("Design: two-arm (treatment vs. control), cluster-level, single post-baseline comparison feeding into a DiD analysis (see script 03).\n"))
cat(sprintf("Baseline mCPR assumed:      %.0f%%\n", p_control * 100))
cat(sprintf("Minimum detectable effect:  %.0f percentage points\n", mde_pp * 100))
cat(sprintf("Alpha / Power:               %.2f / %.2f\n", alpha, power))
cat(sprintf("Assumed ICC:                 %.3f\n", icc))
cat(sprintf("Design effect:               %.2f\n", design_effect))
cat(sprintf("Clusters required per arm:   %d\n", clusters_per_arm))
cat(sprintf("Total clusters:              %d\n", total_clusters))
cat(sprintf("Total planned respondents:   %d\n\n", total_respondents))
cat("Allocation by state:\n")
print(alloc_table, row.names = FALSE)
sink()

cat("\nSaved: ../output/sampling_allocation_by_state.csv\n")
cat("Saved: ../output/sampling_design_summary.txt\n")
