# ==============================================================================
# TUDOR PIPELINE: STEP 05 — LANCET/BMJ-REQUIRED STATISTICS
# ==============================================================================
# PURPOSE: Compute all statistics required for Lancet/BMJ submission:
#          NRI, IDI, DCA, Calibration (slope/intercept/H-L), Brier score.
#
# REQUIRES: tudor_analysis_ready.rds and 02_validation_results.rds
#
# OUTPUTS:  - NRI (categorical & continuous)
#           - IDI with bootstrap CI
#           - Decision Curve Analysis
#           - Calibration metrics
#           - Brier scores
#           - 05_lancet_stats.rds
# ==============================================================================

set.seed(42)

library(data.table)
library(dplyr)
library(pROC)

DATA_DIR <- Sys.getenv("TUDOR_DATA_DIR", unset = "")
if (DATA_DIR == "") {
  if (file.exists(file.path(getwd(), "TUDOR_UKB_Features.csv"))) {
    DATA_DIR <- getwd()
  } else {
    DATA_DIR <- "C:/Users/nader/Downloads"
  }
}
OUTPUT_DIR <- file.path(DATA_DIR, "tudor_pipeline_output")

cat("=== TUDOR PIPELINE: 05_lancet_statistics.R ===\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
df <- readRDS(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))
if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
val_results <- readRDS(file.path(OUTPUT_DIR, "02_validation_results.rds"))

hr <- df[df$cohort_high_risk == TRUE, ]
cat("High-risk cohort:", nrow(hr), "participants\n")
cat("Genetic FH cases:", sum(hr$is_fh_genetic), "\n\n")

# ==============================================================================
# 2. CALIBRATED PROBABILITIES (via GLM)
# ==============================================================================
# TUDOR provides a logistic score; we need calibrated probabilities for NRI/DCA.
# Fit a recalibration GLM using the TUDOR score as the single predictor.
# This preserves the discrimination (AUC) while recalibrating to UKB prevalence.

cat("--- Calibrating Probabilities ---\n")

# TUDOR calibrated probability
glm_tudor <- glm(is_fh_genetic ~ tudor_score, data = hr, family = binomial)
hr$tudor_calib_prob <- predict(glm_tudor, type = "response")

# eDLCN calibrated probability
glm_edlcn <- glm(is_fh_genetic ~ edlcn_score, data = hr, family = binomial)
hr$edlcn_calib_prob <- predict(glm_edlcn, type = "response")

cat(sprintf("TUDOR calibrated prob range: [%.4f, %.4f]\n",
            min(hr$tudor_calib_prob), max(hr$tudor_calib_prob)))
cat(sprintf("eDLCN calibrated prob range: [%.4f, %.4f]\n\n",
            min(hr$edlcn_calib_prob), max(hr$edlcn_calib_prob)))

# ==============================================================================
# 3. NET RECLASSIFICATION INDEX (NRI)
# ==============================================================================
cat("==========================================================\n")
cat(" NET RECLASSIFICATION INDEX (NRI)\n")
cat("==========================================================\n")

# --- 3a. Categorical NRI (at clinical thresholds) ---
# Risk categories: Low (<1%), Moderate (1-5%), High (>5%)
categorize_risk <- function(prob) {
  cut(prob, breaks = c(-Inf, 0.01, 0.05, Inf),
      labels = c("Low", "Moderate", "High"))
}

hr$tudor_cat <- categorize_risk(hr$tudor_calib_prob)
hr$edlcn_cat <- categorize_risk(hr$edlcn_calib_prob)

cat("\n--- Reclassification Table ---\n")
reclass_table <- table(eDLCN = hr$edlcn_cat, TUDOR = hr$tudor_cat)
cat("All participants:\n")
print(reclass_table)

cat("\nFH Cases:\n")
print(table(eDLCN = hr$edlcn_cat[hr$is_fh_genetic == 1],
            TUDOR = hr$tudor_cat[hr$is_fh_genetic == 1]))

cat("\nNon-FH:\n")
print(table(eDLCN = hr$edlcn_cat[hr$is_fh_genetic == 0],
            TUDOR = hr$tudor_cat[hr$is_fh_genetic == 0]))

# Calculate categorical NRI manually
calculate_nri <- function(outcome, prob_new, prob_old, thresholds = c(0.01, 0.05)) {
  cat_new <- as.numeric(cut(prob_new, breaks = c(-Inf, thresholds, Inf)))
  cat_old <- as.numeric(cut(prob_old, breaks = c(-Inf, thresholds, Inf)))

  events <- outcome == 1
  nonevents <- outcome == 0

  # Events: proportion moving up minus proportion moving down
  up_events <- mean(cat_new[events] > cat_old[events])
  down_events <- mean(cat_new[events] < cat_old[events])
  nri_events <- up_events - down_events

  # Non-events: proportion moving down minus proportion moving up
  up_nonevents <- mean(cat_new[nonevents] > cat_old[nonevents])
  down_nonevents <- mean(cat_new[nonevents] < cat_old[nonevents])
  nri_nonevents <- down_nonevents - up_nonevents

  nri_total <- nri_events + nri_nonevents

  list(
    nri = nri_total,
    nri_events = nri_events,
    nri_nonevents = nri_nonevents,
    up_events = up_events, down_events = down_events,
    up_nonevents = up_nonevents, down_nonevents = down_nonevents
  )
}

nri_cat <- calculate_nri(hr$is_fh_genetic, hr$tudor_calib_prob, hr$edlcn_calib_prob)

cat(sprintf("\nCategorical NRI (TUDOR vs eDLCN): %.3f\n", nri_cat$nri))
cat(sprintf("  Events NRI:     %+.3f (up: %.1f%%, down: %.1f%%)\n",
            nri_cat$nri_events,
            nri_cat$up_events * 100, nri_cat$down_events * 100))
cat(sprintf("  Non-events NRI: %+.3f (up: %.1f%%, down: %.1f%%)\n",
            nri_cat$nri_nonevents,
            nri_cat$up_nonevents * 100, nri_cat$down_nonevents * 100))

# --- 3b. Bootstrap CI for NRI ---
n_boot <- 2000
nri_boot <- numeric(n_boot)

for (b in seq_len(n_boot)) {
  idx <- sample(nrow(hr), replace = TRUE)
  boot_nri <- calculate_nri(hr$is_fh_genetic[idx],
                             hr$tudor_calib_prob[idx],
                             hr$edlcn_calib_prob[idx])
  nri_boot[b] <- boot_nri$nri
}

nri_ci <- quantile(nri_boot, c(0.025, 0.975))
cat(sprintf("  Bootstrap 95%% CI: [%.3f, %.3f]\n", nri_ci[1], nri_ci[2]))

# --- 3c. Continuous NRI ---
cat("\n--- Continuous NRI ---\n")
events <- hr$is_fh_genetic == 1
nonevents <- hr$is_fh_genetic == 0

cnri_events <- mean(hr$tudor_calib_prob[events] > hr$edlcn_calib_prob[events]) -
               mean(hr$tudor_calib_prob[events] < hr$edlcn_calib_prob[events])
cnri_nonevents <- mean(hr$tudor_calib_prob[nonevents] < hr$edlcn_calib_prob[nonevents]) -
                  mean(hr$tudor_calib_prob[nonevents] > hr$edlcn_calib_prob[nonevents])
cnri_total <- cnri_events + cnri_nonevents

cat(sprintf("Continuous NRI: %.3f\n", cnri_total))
cat(sprintf("  Events:     %+.3f\n", cnri_events))
cat(sprintf("  Non-events: %+.3f\n\n", cnri_nonevents))

# ==============================================================================
# 4. INTEGRATED DISCRIMINATION IMPROVEMENT (IDI)
# ==============================================================================
cat("==========================================================\n")
cat(" INTEGRATED DISCRIMINATION IMPROVEMENT (IDI)\n")
cat("==========================================================\n")

# IDI = (mean predicted prob in events_new - mean predicted prob in events_old)
#     - (mean predicted prob in nonevents_new - mean predicted prob in nonevents_old)

idi_events <- mean(hr$tudor_calib_prob[events]) - mean(hr$edlcn_calib_prob[events])
idi_nonevents <- mean(hr$tudor_calib_prob[nonevents]) - mean(hr$edlcn_calib_prob[nonevents])
idi <- idi_events - idi_nonevents

cat(sprintf("IDI (TUDOR vs eDLCN): %.4f\n", idi))
cat(sprintf("  IS (events):     %+.4f\n", idi_events))
cat(sprintf("  IP (non-events): %+.4f\n", idi_nonevents))

# Bootstrap CI for IDI
idi_boot <- numeric(n_boot)
for (b in seq_len(n_boot)) {
  idx <- sample(nrow(hr), replace = TRUE)
  ev <- hr$is_fh_genetic[idx] == 1
  nev <- hr$is_fh_genetic[idx] == 0
  ie <- mean(hr$tudor_calib_prob[idx][ev]) - mean(hr$edlcn_calib_prob[idx][ev])
  ine <- mean(hr$tudor_calib_prob[idx][nev]) - mean(hr$edlcn_calib_prob[idx][nev])
  idi_boot[b] <- ie - ine
}

idi_ci <- quantile(idi_boot, c(0.025, 0.975))
cat(sprintf("  Bootstrap 95%% CI: [%.4f, %.4f]\n", idi_ci[1], idi_ci[2]))
cat(sprintf("  P-value: %.2e\n\n",
            2 * min(mean(idi_boot <= 0), mean(idi_boot >= 0))))

# ==============================================================================
# 5. DECISION CURVE ANALYSIS (DCA)
# ==============================================================================
cat("==========================================================\n")
cat(" DECISION CURVE ANALYSIS (DCA)\n")
cat("==========================================================\n")

# DCA: Net benefit at various threshold probabilities
dca_thresholds <- seq(0.001, 0.10, by = 0.001)
prevalence <- mean(hr$is_fh_genetic)

dca_results <- data.frame(
  threshold = dca_thresholds,
  nb_tudor = NA_real_,
  nb_edlcn = NA_real_,
  nb_treat_all = NA_real_,
  nb_treat_none = 0
)

for (i in seq_along(dca_thresholds)) {
  pt <- dca_thresholds[i]
  odds <- pt / (1 - pt)

  # TUDOR
  tudor_pos <- hr$tudor_calib_prob >= pt
  tp_tudor <- sum(tudor_pos & hr$is_fh_genetic == 1) / nrow(hr)
  fp_tudor <- sum(tudor_pos & hr$is_fh_genetic == 0) / nrow(hr)
  dca_results$nb_tudor[i] <- tp_tudor - fp_tudor * odds

  # eDLCN
  edlcn_pos <- hr$edlcn_calib_prob >= pt
  tp_edlcn <- sum(edlcn_pos & hr$is_fh_genetic == 1) / nrow(hr)
  fp_edlcn <- sum(edlcn_pos & hr$is_fh_genetic == 0) / nrow(hr)
  dca_results$nb_edlcn[i] <- tp_edlcn - fp_edlcn * odds

  # Treat all
  dca_results$nb_treat_all[i] <- prevalence - (1 - prevalence) * odds
}

# Summary of DCA
best_tudor <- sum(dca_results$nb_tudor > dca_results$nb_edlcn, na.rm = TRUE)
cat(sprintf("TUDOR has higher net benefit at %d/%d thresholds (%.0f%%)\n",
            best_tudor, length(dca_thresholds),
            best_tudor / length(dca_thresholds) * 100))
cat(sprintf("Threshold range evaluated: %.1f%% to %.1f%%\n",
            min(dca_thresholds) * 100, max(dca_thresholds) * 100))
cat("\nSample thresholds:\n")
for (pt in c(0.01, 0.02, 0.05)) {
  idx <- which.min(abs(dca_thresholds - pt))
  cat(sprintf("  At %.0f%%: TUDOR NB = %.4f, eDLCN NB = %.4f, Treat-All = %.4f\n",
              pt * 100, dca_results$nb_tudor[idx],
              dca_results$nb_edlcn[idx], dca_results$nb_treat_all[idx]))
}
cat("\n")

# ==============================================================================
# 6. CALIBRATION
# ==============================================================================
cat("==========================================================\n")
cat(" CALIBRATION ANALYSIS\n")
cat("==========================================================\n")

# --- 6a. Calibration Slope and Intercept ---
# Fit logistic regression of outcome on predicted log-odds
# Perfect calibration: intercept = 0, slope = 1

tudor_logodds <- log(hr$tudor_calib_prob / (1 - hr$tudor_calib_prob))
calib_fit <- glm(is_fh_genetic ~ tudor_logodds, data = hr, family = binomial)

calib_intercept <- coef(calib_fit)[1]
calib_slope <- coef(calib_fit)[2]

cat(sprintf("Calibration intercept: %.3f (ideal: 0)\n", calib_intercept))
cat(sprintf("Calibration slope:     %.3f (ideal: 1)\n", calib_slope))

# CIs for calibration slope
calib_ci <- confint.default(calib_fit)
cat(sprintf("  Slope 95%% CI: [%.3f, %.3f]\n",
            calib_ci["tudor_logodds", 1], calib_ci["tudor_logodds", 2]))

# --- 6b. Hosmer-Lemeshow Test ---
# Divide into deciles of predicted probability
hr$prob_decile <- cut(hr$tudor_calib_prob,
                       breaks = quantile(hr$tudor_calib_prob,
                                          probs = seq(0, 1, 0.1),
                                          na.rm = TRUE),
                       include.lowest = TRUE, labels = FALSE)

hl_table <- hr %>%
  group_by(prob_decile) %>%
  summarise(
    n = n(),
    observed = sum(is_fh_genetic),
    expected = sum(tudor_calib_prob),
    obs_rate = mean(is_fh_genetic),
    exp_rate = mean(tudor_calib_prob),
    .groups = "drop"
  )

# H-L chi-squared statistic
hl_chi2 <- sum((hl_table$observed - hl_table$expected)^2 /
                (hl_table$expected * (1 - hl_table$exp_rate) + 0.001))
hl_df <- nrow(hl_table) - 2
hl_p <- 1 - pchisq(hl_chi2, df = hl_df)

cat(sprintf("\nHosmer-Lemeshow chi2: %.2f (df=%d, p=%.3f)\n", hl_chi2, hl_df, hl_p))
cat("(p > 0.05 indicates acceptable calibration)\n\n")

cat("Calibration by decile:\n")
cat(sprintf("%-8s | %6s | %8s | %8s | %10s | %10s\n",
            "Decile", "N", "Observed", "Expected", "Obs Rate", "Exp Rate"))
cat(paste(rep("-", 65), collapse = ""), "\n")
for (i in seq_len(nrow(hl_table))) {
  r <- hl_table[i, ]
  cat(sprintf("%-8d | %6d | %8d | %8.1f | %9.4f | %9.4f\n",
              r$prob_decile, r$n, r$observed, r$expected,
              r$obs_rate, r$exp_rate))
}
cat(paste(rep("-", 65), collapse = ""), "\n\n")

# ==============================================================================
# 7. BRIER SCORE
# ==============================================================================
cat("==========================================================\n")
cat(" BRIER SCORE\n")
cat("==========================================================\n")

brier_tudor <- mean((hr$tudor_calib_prob - hr$is_fh_genetic)^2)
brier_edlcn <- mean((hr$edlcn_calib_prob - hr$is_fh_genetic)^2)

# Scaled Brier score (accounts for prevalence)
brier_max <- prevalence * (1 - prevalence)
scaled_brier_tudor <- 1 - brier_tudor / brier_max
scaled_brier_edlcn <- 1 - brier_edlcn / brier_max

cat(sprintf("TUDOR Brier score:  %.4f (scaled: %.3f)\n", brier_tudor, scaled_brier_tudor))
cat(sprintf("eDLCN Brier score:  %.4f (scaled: %.3f)\n", brier_edlcn, scaled_brier_edlcn))
cat(sprintf("Max Brier (null):   %.4f\n", brier_max))
cat("(Lower Brier = better; higher scaled Brier = better)\n\n")

# Bootstrap CI for Brier difference
brier_diff_boot <- numeric(n_boot)
for (b in seq_len(n_boot)) {
  idx <- sample(nrow(hr), replace = TRUE)
  bt <- mean((hr$tudor_calib_prob[idx] - hr$is_fh_genetic[idx])^2)
  be <- mean((hr$edlcn_calib_prob[idx] - hr$is_fh_genetic[idx])^2)
  brier_diff_boot[b] <- bt - be
}

brier_diff_ci <- quantile(brier_diff_boot, c(0.025, 0.975))
cat(sprintf("Brier difference (TUDOR - eDLCN): %.4f (95%% CI: [%.4f, %.4f])\n",
            brier_tudor - brier_edlcn, brier_diff_ci[1], brier_diff_ci[2]))
cat("(Negative = TUDOR better)\n\n")

# ==============================================================================
# 8. SAVE RESULTS
# ==============================================================================
lancet_results <- list(
  nri = list(
    categorical = nri_cat,
    categorical_ci = nri_ci,
    continuous = list(total = cnri_total, events = cnri_events,
                      nonevents = cnri_nonevents)
  ),
  idi = list(
    idi = idi, idi_events = idi_events, idi_nonevents = idi_nonevents,
    ci = idi_ci
  ),
  dca = dca_results,
  calibration = list(
    intercept = calib_intercept, slope = calib_slope,
    slope_ci = calib_ci["tudor_logodds", ],
    hl_chi2 = hl_chi2, hl_df = hl_df, hl_p = hl_p,
    decile_table = hl_table
  ),
  brier = list(
    tudor = brier_tudor, edlcn = brier_edlcn,
    scaled_tudor = scaled_brier_tudor, scaled_edlcn = scaled_brier_edlcn,
    diff_ci = brier_diff_ci
  )
)

saveRDS(lancet_results, file.path(OUTPUT_DIR, "05_lancet_stats.rds"))

cat("=== 05_lancet_statistics.R COMPLETE ===\n")
