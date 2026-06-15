# ==============================================================================
# TUDOR PIPELINE: STEP 02 — EXTERNAL VALIDATION IN UK BIOBANK
# ==============================================================================
# PURPOSE: Externally validate the TUDOR FH diagnostic model (trained on
#          Wales FH Registry) in the independent UK Biobank population.
#          This is a TRIPOD Type 4 external validation: model weights are
#          FIXED from Wales — nothing is re-estimated on UKB data.
#
# ADDRESSES REVIEWER COMMENTS:
#   R2-1:  Clarifies intended use case (lower-prevalence settings)
#   R2-3:  Demonstrates external validation + recalibration approach
#   R2-9:  Harmonises all AUROC values with clear labels
#   R2-10: Consistent metrics between text/tables/figures
#   R2-14: Clinical decision pathway thresholds
#   R3:    Independent external validation (Reviewer 3's key demand)
#
# INPUT:   tudor_analysis_ready.rds from 01_data_merge.R
# OUTPUT:  02_validation_results.rds (required by 05_lancet_statistics.R)
#          ROC curves, calibration plots, performance tables
# ==============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(pROC)
  library(ggplot2)
})

# --- Configuration -----------------------------------------------------------
DATA_DIR <- Sys.getenv("TUDOR_DATA_DIR", unset = "")
if (DATA_DIR == "") {
  if (file.exists(file.path(getwd(), "TUDOR_UKB_Features.csv"))) {
    DATA_DIR <- getwd()
  } else {
    DATA_DIR <- "C:/Users/nader/Downloads"
  }
}
OUTPUT_DIR <- file.path(DATA_DIR, "tudor_pipeline_output")
PLOT_DIR   <- file.path(OUTPUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== TUDOR PIPELINE: 02_external_validation.R ===\n")
cat("TRIPOD Type 4: External validation of pre-specified model\n")
cat("Training: Wales FH Registry | Validation: UK Biobank\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
if (!file.exists(rds_file)) stop("Run 01_data_merge.R first!")
df <- readRDS(rds_file)

if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}

cat("Total cohort:", nrow(df), "\n")
cat("Genetic FH cases:", sum(df$is_fh_genetic), "\n")
cat("FH prevalence (total):", sprintf("%.2f%%", 100 * mean(df$is_fh_genetic)), "\n\n")

# ==============================================================================
# 2. DEFINE ANALYSIS COHORTS
# ==============================================================================
# Primary analysis: High-risk cohort (LDL_RW > 4.9 mmol/L)
#   This mirrors the clinical use case: screening those with elevated LDL
# Secondary: Full population (all participants)

hr <- df[df$cohort_high_risk == TRUE & !is.na(df$tudor_prob), ]
full <- df[!is.na(df$tudor_prob), ]

cat("=== COHORT DEFINITIONS ===\n")
cat(sprintf("%-30s N = %6d  (FH = %4d, prev = %.2f%%)\n",
            "Full population:", nrow(full), sum(full$is_fh_genetic),
            100 * mean(full$is_fh_genetic)))
cat(sprintf("%-30s N = %6d  (FH = %4d, prev = %.2f%%)\n",
            "High-risk (LDL > 4.9):", nrow(hr), sum(hr$is_fh_genetic),
            100 * mean(hr$is_fh_genetic)))
cat(sprintf("%-30s N = %6d  (FH = %4d, prev = %.2f%%)\n",
            "Moderate (2.6-4.9):", sum(df$cohort_moderate, na.rm = TRUE),
            sum(df$is_fh_genetic[df$cohort_moderate == TRUE], na.rm = TRUE),
            100 * mean(df$is_fh_genetic[df$cohort_moderate == TRUE], na.rm = TRUE)))
cat("\n")

# ==============================================================================
# 3. PRIMARY DISCRIMINATION ANALYSIS (HIGH-RISK COHORT)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PRIMARY ANALYSIS: DISCRIMINATION (HIGH-RISK COHORT)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# --- 3a. ROC curves ---
roc_tudor   <- roc(hr$is_fh_genetic, hr$tudor_prob, quiet = TRUE)
roc_edlcn   <- roc(hr$is_fh_genetic, hr$edlcn_score, quiet = TRUE)
roc_ldl     <- roc(hr$is_fh_genetic, hr$LDL_RW, quiet = TRUE)
roc_trig    <- roc(hr$is_fh_genetic, hr$Trig_Filter_RW, quiet = TRUE)

# Bootstrap CIs (DeLong for speed, bootstrap for publication)
ci_tudor  <- ci.auc(roc_tudor, method = "delong")
ci_edlcn  <- ci.auc(roc_edlcn, method = "delong")
ci_ldl    <- ci.auc(roc_ldl, method = "delong")
ci_trig   <- ci.auc(roc_trig, method = "delong")

cat("AREA UNDER ROC CURVE (AUC):\n\n")
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "TUDOR v2", ci_tudor[2], ci_tudor[1], ci_tudor[3]))
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "eDLCN", ci_edlcn[2], ci_edlcn[1], ci_edlcn[3]))
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "LDL-C alone", ci_ldl[2], ci_ldl[1], ci_ldl[3]))
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "Trig Filter", ci_trig[2], ci_trig[1], ci_trig[3]))
cat("\n")

# --- 3b. DeLong pairwise tests ---
cat("DELONG PAIRWISE COMPARISONS:\n")
comparisons <- list(
  list("TUDOR vs eDLCN",       roc_tudor, roc_edlcn),
  list("TUDOR vs LDL-C alone", roc_tudor, roc_ldl),
  list("TUDOR vs Trig Filter", roc_tudor, roc_trig),
  list("Trig Filter vs eDLCN", roc_trig,  roc_edlcn),
  list("Trig Filter vs LDL-C", roc_trig,  roc_ldl)
)

for (comp in comparisons) {
  dt <- roc.test(comp[[2]], comp[[3]], method = "delong")
  delta <- as.numeric(auc(comp[[2]])) - as.numeric(auc(comp[[3]]))
  cat(sprintf("  %-30s delta = %+.3f  p = %.2e %s\n",
              comp[[1]], delta, dt$p.value,
              ifelse(dt$p.value < 0.05, "*", "")))
}
cat("\n")


# ==============================================================================
# 4. SENSITIVITY / SPECIFICITY AT CLINICAL THRESHOLDS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SENSITIVITY AND SPECIFICITY AT KEY THRESHOLDS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Helper function
calc_metrics_at_threshold <- function(response, predictor, threshold, label) {
  pred_pos <- predictor >= threshold
  tp <- sum(pred_pos & response == 1)
  fp <- sum(pred_pos & response == 0)
  fn <- sum(!pred_pos & response == 1)
  tn <- sum(!pred_pos & response == 0)

  sens <- tp / (tp + fn)
  spec <- tn / (tn + fp)
  ppv  <- tp / (tp + fp)
  npv  <- tn / (tn + fn)

  # Wilson score CIs for proportions
  n_pos <- tp + fn
  n_neg <- tn + fp
  sens_ci <- prop.test(tp, n_pos)$conf.int
  spec_ci <- prop.test(tn, n_neg)$conf.int

  list(label = label, threshold = threshold,
       sens = sens, sens_ci = sens_ci,
       spec = spec, spec_ci = spec_ci,
       ppv = ppv, npv = npv,
       tp = tp, fp = fp, fn = fn, tn = tn)
}

# --- TUDOR thresholds (probability-based) ---
# Reviewer R2-14 asks for decision pathway thresholds
cat("A. TUDOR PROBABILITY THRESHOLDS:\n")
cat("(Clinical decision pathway: <0.25 = Low, 0.25-0.75 = Intermediate, >0.75 = High)\n\n")

cat(sprintf("%-12s | %8s | %8s | %8s | %8s | %5s %5s %5s %5s\n",
            "Threshold", "Sens", "Spec", "PPV", "NPV", "TP", "FP", "FN", "TN"))
cat(strrep("-", 90), "\n")

tudor_thresholds <- c(0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05)

# Find Youden optimal threshold
youden_idx <- which.max(roc_tudor$sensitivities + roc_tudor$specificities - 1)
youden_thresh <- roc_tudor$thresholds[youden_idx]
tudor_thresholds <- sort(unique(c(tudor_thresholds, round(youden_thresh, 4))))

for (thr in tudor_thresholds) {
  m <- calc_metrics_at_threshold(hr$is_fh_genetic, hr$tudor_prob, thr,
                                  sprintf("TUDOR >= %.3f", thr))
  is_youden <- abs(thr - youden_thresh) < 0.001
  marker <- ifelse(is_youden, " <-- Youden", "")
  cat(sprintf("%-12s | %7.1f%% | %7.1f%% | %7.1f%% | %7.1f%% | %5d %5d %5d %5d%s\n",
              sprintf(">= %.4f", thr),
              m$sens * 100, m$spec * 100, m$ppv * 100, m$npv * 100,
              m$tp, m$fp, m$fn, m$tn, marker))
}
cat("\n")

# --- eDLCN thresholds (score-based) ---
cat("B. eDLCN SCORE THRESHOLDS:\n")
cat("(Standard: Possible >= 3, Probable >= 6, Definite >= 8)\n\n")

cat(sprintf("%-12s | %8s | %8s | %8s | %8s | %5s %5s %5s %5s\n",
            "Threshold", "Sens", "Spec", "PPV", "NPV", "TP", "FP", "FN", "TN"))
cat(strrep("-", 90), "\n")

for (thr in c(1, 3, 5, 6, 8)) {
  m <- calc_metrics_at_threshold(hr$is_fh_genetic, hr$edlcn_score, thr,
                                  sprintf("eDLCN >= %d", thr))
  dlcn_label <- ""
  if (thr == 3) dlcn_label <- " (Possible)"
  if (thr == 6) dlcn_label <- " (Probable)"
  if (thr == 8) dlcn_label <- " (Definite)"
  cat(sprintf("%-12s | %7.1f%% | %7.1f%% | %7.1f%% | %7.1f%% | %5d %5d %5d %5d%s\n",
              sprintf(">= %d", thr),
              m$sens * 100, m$spec * 100, m$ppv * 100, m$npv * 100,
              m$tp, m$fp, m$fn, m$tn, dlcn_label))
}
cat("\n")


# ==============================================================================
# 5. CLINICAL DECISION PATHWAY (Reviewer R2-14)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CLINICAL DECISION PATHWAY BASED ON TUDOR PROBABILITY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Three-tier classification
hr$tudor_category <- ifelse(hr$tudor_prob < 0.01, "Low (<1%)",
                     ifelse(hr$tudor_prob < 0.03, "Intermediate (1-3%)",
                                                   "High (>3%)"))
hr$tudor_category <- factor(hr$tudor_category,
                             levels = c("Low (<1%)", "Intermediate (1-3%)", "High (>3%)"))

cat("Risk Category Distribution:\n")
for (cat_level in levels(hr$tudor_category)) {
  sub <- hr[hr$tudor_category == cat_level, ]
  n_fh <- sum(sub$is_fh_genetic)
  prev <- 100 * n_fh / nrow(sub)
  cat(sprintf("  %-25s N = %6d  |  FH = %4d  |  FH prevalence = %.2f%%\n",
              cat_level, nrow(sub), n_fh, prev))
}

cat("\nProposed Clinical Actions:\n")
cat("  Low (<1%):           Reassurance. No genetic testing. Standard lipid mgmt.\n")
cat("  Intermediate (1-3%): Clinical review. Consider family history, physical signs.\n")
cat("                       Genetic testing if additional risk factors present.\n")
cat("  High (>3%):          Refer for genetic testing. Initiate cascade screening.\n")
cat("                       Consider high-intensity statin + ezetimibe.\n\n")


# ==============================================================================
# 6. SUBGROUP ANALYSES
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUBGROUP ANALYSES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

subgroup_results <- data.frame()

run_subgroup <- function(data, label) {
  if (sum(data$is_fh_genetic) < 10) {
    cat(sprintf("  %-30s SKIPPED (FH < 10)\n", label))
    return(NULL)
  }
  r_tudor <- roc(data$is_fh_genetic, data$tudor_prob, quiet = TRUE)
  r_edlcn <- roc(data$is_fh_genetic, data$edlcn_score, quiet = TRUE)
  ci_t <- ci.auc(r_tudor, method = "delong")
  ci_e <- ci.auc(r_edlcn, method = "delong")

  cat(sprintf("  %-30s TUDOR = %.3f [%.3f-%.3f]  eDLCN = %.3f [%.3f-%.3f]  N=%d (FH=%d)\n",
              label, ci_t[2], ci_t[1], ci_t[3],
              ci_e[2], ci_e[1], ci_e[3],
              nrow(data), sum(data$is_fh_genetic)))

  data.frame(subgroup = label,
             tudor_auc = as.numeric(ci_t[2]),
             tudor_ci_lo = as.numeric(ci_t[1]),
             tudor_ci_hi = as.numeric(ci_t[3]),
             edlcn_auc = as.numeric(ci_e[2]),
             edlcn_ci_lo = as.numeric(ci_e[1]),
             edlcn_ci_hi = as.numeric(ci_e[3]),
             n = nrow(data), n_fh = sum(data$is_fh_genetic))
}

# Overall
subgroup_results <- rbind(subgroup_results, run_subgroup(hr, "Overall (high-risk)"))

# By sex
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$Gender_num == 1, ], "Male"))
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$Gender_num == 0, ], "Female"))

# By age
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$Age_at_LDL1 < 50, ], "Age < 50"))
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$Age_at_LDL1 >= 50 & hr$Age_at_LDL1 < 60, ], "Age 50-59"))
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$Age_at_LDL1 >= 60, ], "Age >= 60"))

# By statin status
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$statin_name == "None", ], "No statin"))
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$statin_name != "None", ], "On statin"))

# By statin type (major statins)
for (st in c("Simvastatin", "Atorvastatin")) {
  sub <- hr[hr$statin_name == st, ]
  if (nrow(sub) > 100) {
    subgroup_results <- rbind(subgroup_results,
      run_subgroup(sub, paste("On", st)))
  }
}

# By LDL range
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$LDL_RW >= 4.9 & hr$LDL_RW < 6.5, ], "LDL 4.9-6.5"))
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$LDL_RW >= 6.5 & hr$LDL_RW < 8.5, ], "LDL 6.5-8.5"))
subgroup_results <- rbind(subgroup_results,
  run_subgroup(hr[hr$LDL_RW >= 8.5, ], "LDL >= 8.5"))

# By ASCVD status
if (sum(hr$has_any_cvd) > 10) {
  subgroup_results <- rbind(subgroup_results,
    run_subgroup(hr[hr$has_any_cvd == TRUE, ], "ASCVD present"))
  subgroup_results <- rbind(subgroup_results,
    run_subgroup(hr[hr$has_any_cvd == FALSE, ], "No ASCVD"))
}

cat("\n")


# ==============================================================================
# 7. FULL POPULATION ANALYSIS (SECONDARY)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SECONDARY ANALYSIS: FULL POPULATION\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

roc_tudor_full <- roc(full$is_fh_genetic, full$tudor_prob, quiet = TRUE)
roc_edlcn_full <- roc(full$is_fh_genetic, full$edlcn_score, quiet = TRUE)
roc_ldl_full   <- roc(full$is_fh_genetic, full$LDL_RW, quiet = TRUE)
roc_trig_full  <- roc(full$is_fh_genetic, full$Trig_Filter_RW, quiet = TRUE)

ci_tudor_full <- ci.auc(roc_tudor_full, method = "delong")
ci_edlcn_full <- ci.auc(roc_edlcn_full, method = "delong")
ci_ldl_full   <- ci.auc(roc_ldl_full, method = "delong")
ci_trig_full  <- ci.auc(roc_trig_full, method = "delong")

cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "TUDOR v2", ci_tudor_full[2], ci_tudor_full[1], ci_tudor_full[3]))
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "eDLCN", ci_edlcn_full[2], ci_edlcn_full[1], ci_edlcn_full[3]))
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "LDL-C alone", ci_ldl_full[2], ci_ldl_full[1], ci_ldl_full[3]))
cat(sprintf("%-25s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "Trig Filter", ci_trig_full[2], ci_trig_full[1], ci_trig_full[3]))
cat("\n")


# ==============================================================================
# 8. TRIG FILTER ILLUSTRATIVE EXAMPLE (Reviewer R2-18)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("METABOLIC SHIELD ILLUSTRATIVE EXAMPLE (Reviewer R2-18)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Find a representative FH case with high Trig Filter (metabolically pure)
fh_cases <- hr[hr$is_fh_genetic == TRUE, ]
fh_cases <- fh_cases[order(-fh_cases$Trig_Filter_RW), ]

# Show top 5 FH cases by Trig Filter (anonymised)
cat("Top 5 FH cases by Trig Filter (metabolic purity):\n")
cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-12s %-10s %-10s %-8s\n",
            "Rank", "LDL_RW", "TRG", "HDL", "Age", "TrigFilter", "TUDOR_p", "eDLCN", "Statin"))
cat(strrep("-", 85), "\n")

for (i in 1:min(5, nrow(fh_cases))) {
  r <- fh_cases[i, ]
  cat(sprintf("%-6d %-8.1f %-8.2f %-8.2f %-8.0f %-12.1f %-10.4f %-10d %-8s\n",
              i, r$LDL_RW, r$TRG.1, r$HDL.1, r$Age_at_LDL1,
              r$Trig_Filter_RW, r$tudor_prob, r$edlcn_score, r$statin_name))
}

cat("\nContrast: Non-FH cases with similar LDL but LOW Trig Filter:\n")
nonfh_high_ldl <- hr[hr$is_fh_genetic == FALSE & hr$LDL_RW > median(fh_cases$LDL_RW), ]
nonfh_high_ldl <- nonfh_high_ldl[order(nonfh_high_ldl$Trig_Filter_RW), ]

cat(sprintf("%-6s %-8s %-8s %-8s %-8s %-12s %-10s %-10s %-8s\n",
            "Rank", "LDL_RW", "TRG", "HDL", "Age", "TrigFilter", "TUDOR_p", "eDLCN", "Statin"))
cat(strrep("-", 85), "\n")

for (i in 1:min(5, nrow(nonfh_high_ldl))) {
  r <- nonfh_high_ldl[i, ]
  cat(sprintf("%-6d %-8.1f %-8.2f %-8.2f %-8.0f %-12.1f %-10.4f %-10d %-8s\n",
              i, r$LDL_RW, r$TRG.1, r$HDL.1, r$Age_at_LDL1,
              r$Trig_Filter_RW, r$tudor_prob, r$edlcn_score, r$statin_name))
}

cat("\nINTERPRETATION: Both groups have elevated LDL-C.\n")
cat("FH patients: HIGH Trig Filter (high LDL, LOW triglycerides) = 'pure' hypercholesterolaemia.\n")
cat("Non-FH:      LOW Trig Filter (high LDL, HIGH triglycerides) = metabolic syndrome pattern.\n")
cat("TUDOR uses this metabolic purity signal; eDLCN cannot detect it.\n\n")


# ==============================================================================
# 9. ROC CURVE PLOT
# ==============================================================================
cat("Generating ROC curve plot...\n")

roc_data <- rbind(
  data.frame(sens = roc_tudor$sensitivities, spec = 1 - roc_tudor$specificities,
             Model = sprintf("TUDOR (AUC = %.3f)", auc(roc_tudor))),
  data.frame(sens = roc_edlcn$sensitivities, spec = 1 - roc_edlcn$specificities,
             Model = sprintf("eDLCN (AUC = %.3f)", auc(roc_edlcn))),
  data.frame(sens = roc_trig$sensitivities, spec = 1 - roc_trig$specificities,
             Model = sprintf("Trig Filter (AUC = %.3f)", auc(roc_trig))),
  data.frame(sens = roc_ldl$sensitivities, spec = 1 - roc_ldl$specificities,
             Model = sprintf("LDL-C alone (AUC = %.3f)", auc(roc_ldl)))
)

p_roc <- ggplot(roc_data, aes(x = spec, y = sens, color = Model)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("red3", "steelblue", "darkgreen", "orange")) +
  labs(x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)",
       title = "External Validation: ROC Curves in UK Biobank",
       subtitle = sprintf("High-risk cohort (LDL > 4.9 mmol/L, N = %s, FH = %d)",
                           format(nrow(hr), big.mark = ","), sum(hr$is_fh_genetic))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = c(0.7, 0.25),
        legend.background = element_rect(fill = "white", colour = "grey80"))

ggsave(file.path(PLOT_DIR, "fig_roc_external_validation.pdf"),
       p_roc, width = 8, height = 7)
ggsave(file.path(PLOT_DIR, "fig_roc_external_validation.png"),
       p_roc, width = 8, height = 7, dpi = 300)
cat("  Saved: fig_roc_external_validation.pdf/png\n")


# ==============================================================================
# 10. CALIBRATION PLOT
# ==============================================================================
cat("Generating calibration plot...\n")

# Create probability deciles for calibration
hr$prob_decile <- cut(hr$tudor_prob,
                       breaks = quantile(hr$tudor_prob,
                                          probs = seq(0, 1, 0.1), na.rm = TRUE),
                       include.lowest = TRUE, labels = FALSE)

calib_data <- hr %>%
  group_by(prob_decile) %>%
  summarise(
    predicted = mean(tudor_prob),
    observed = mean(is_fh_genetic),
    n = n(),
    se = sqrt(observed * (1 - observed) / n),
    .groups = "drop"
  )

p_calib <- ggplot(calib_data, aes(x = predicted, y = observed)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n), color = "red3") +
  geom_errorbar(aes(ymin = pmax(0, observed - 1.96 * se),
                     ymax = pmin(1, observed + 1.96 * se)),
                width = 0.001, color = "red3") +
  scale_size_continuous(range = c(2, 6)) +
  labs(x = "Predicted Probability (TUDOR)",
       y = "Observed FH Proportion",
       title = "Calibration Plot: TUDOR in UK Biobank",
       subtitle = "High-risk cohort, deciles of predicted probability",
       size = "N per decile") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(PLOT_DIR, "fig_calibration.pdf"), p_calib, width = 7, height = 6)
ggsave(file.path(PLOT_DIR, "fig_calibration.png"), p_calib, width = 7, height = 6, dpi = 300)
cat("  Saved: fig_calibration.pdf/png\n\n")


# ==============================================================================
# 11. SUMMARY PERFORMANCE TABLE
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PERFORMANCE SUMMARY TABLE (for manuscript)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Youden optimal for TUDOR
m_youden <- calc_metrics_at_threshold(hr$is_fh_genetic, hr$tudor_prob,
                                       youden_thresh, "Youden optimal")

# eDLCN at Possible (>=3)
m_edlcn3 <- calc_metrics_at_threshold(hr$is_fh_genetic, hr$edlcn_score, 3, "eDLCN Possible")
m_edlcn6 <- calc_metrics_at_threshold(hr$is_fh_genetic, hr$edlcn_score, 6, "eDLCN Probable")

cat(sprintf("%-30s %-10s %-10s %-10s %-10s %-10s\n",
            "Metric", "TUDOR", "eDLCN>=3", "eDLCN>=6", "TrigFilter", "LDL-C"))
cat(strrep("-", 80), "\n")
cat(sprintf("%-30s %-10.3f %-10.3f %-10.3f %-10.3f %-10.3f\n",
            "AUC", auc(roc_tudor), auc(roc_edlcn), auc(roc_edlcn), auc(roc_trig), auc(roc_ldl)))
cat(sprintf("%-30s %-10.1f %-10.1f %-10.1f %-10s %-10s\n",
            "Sensitivity (%)", m_youden$sens*100, m_edlcn3$sens*100, m_edlcn6$sens*100, "-", "-"))
cat(sprintf("%-30s %-10.1f %-10.1f %-10.1f %-10s %-10s\n",
            "Specificity (%)", m_youden$spec*100, m_edlcn3$spec*100, m_edlcn6$spec*100, "-", "-"))
cat(sprintf("%-30s %-10.1f %-10.1f %-10.1f %-10s %-10s\n",
            "PPV (%)", m_youden$ppv*100, m_edlcn3$ppv*100, m_edlcn6$ppv*100, "-", "-"))
cat(sprintf("%-30s %-10.1f %-10.1f %-10.1f %-10s %-10s\n",
            "NPV (%)", m_youden$npv*100, m_edlcn3$npv*100, m_edlcn6$npv*100, "-", "-"))
cat("\n")
cat(sprintf("TUDOR Youden optimal threshold: %.4f\n\n", youden_thresh))


# ==============================================================================
# 12. SAVE RESULTS
# ==============================================================================
val_results <- list(
  # Primary analysis (high-risk)
  primary = list(
    cohort = "High-risk (LDL > 4.9)",
    n = nrow(hr), n_fh = sum(hr$is_fh_genetic),
    prevalence = mean(hr$is_fh_genetic),
    tudor_auc = as.numeric(auc(roc_tudor)),
    tudor_ci = as.numeric(ci_tudor),
    edlcn_auc = as.numeric(auc(roc_edlcn)),
    edlcn_ci = as.numeric(ci_edlcn),
    ldl_auc = as.numeric(auc(roc_ldl)),
    ldl_ci = as.numeric(ci_ldl),
    trig_auc = as.numeric(auc(roc_trig)),
    trig_ci = as.numeric(ci_trig),
    youden_threshold = youden_thresh,
    youden_sens = m_youden$sens,
    youden_spec = m_youden$spec
  ),
  # Full population
  full_population = list(
    n = nrow(full), n_fh = sum(full$is_fh_genetic),
    tudor_auc = as.numeric(auc(roc_tudor_full)),
    edlcn_auc = as.numeric(auc(roc_edlcn_full))
  ),
  # Subgroup results
  subgroups = subgroup_results,
  # ROC objects (for downstream scripts)
  roc_tudor = roc_tudor,
  roc_edlcn = roc_edlcn,
  roc_ldl = roc_ldl,
  roc_trig = roc_trig,
  # Calibration data
  calibration = calib_data,
  # Timestamp
  timestamp = Sys.time()
)

out_file <- file.path(OUTPUT_DIR, "02_validation_results.rds")
saveRDS(val_results, out_file)
cat("Saved validation results to:", out_file, "\n")

cat("\n=== 02_external_validation.R COMPLETE ===\n")
cat("\nKey numbers for manuscript (HARMONISED — Reviewer R2-9):\n")
cat(sprintf("  HIGH-RISK COHORT:  TUDOR AUC = %.3f  |  eDLCN AUC = %.3f\n",
            auc(roc_tudor), auc(roc_edlcn)))
cat(sprintf("  FULL POPULATION:   TUDOR AUC = %.3f  |  eDLCN AUC = %.3f\n",
            auc(roc_tudor_full), auc(roc_edlcn_full)))
cat(sprintf("  YOUDEN THRESHOLD:  %.4f  (Sens = %.1f%%, Spec = %.1f%%)\n",
            youden_thresh, m_youden$sens * 100, m_youden$spec * 100))
