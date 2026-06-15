# ==============================================================================
# TUDOR PIPELINE: STEP 10 — WALES INTERNAL VALIDATION (Aligned to UKB Format)
# ==============================================================================
# PURPOSE: Present Wales FH Registry TUDOR validation results in the SAME
#          format as the UKB external validation pipeline (scripts 02-07),
#          enabling side-by-side comparison for the manuscript.
#
# METHODOLOGY DIFFERENCES (Wales vs UKB):
#   Wales: Elastic Net (alpha=0.5), Index->Relatives (TRIPOD Type 2b)
#   UKB:   Fixed logistic regression weights (TRIPOD Type 4)
#   Wales: Literature-based dose-specific statin factors
#   UKB:   Real-world drug-class-level reduction factors
#   Wales: DLCN + Simon Broome (clinical only, no DNA)
#   UKB:   eDLCN only (no family history/physical signs)
#
# PREREQUISITE: Run TUDOR_v2_Clean.R first, then save workspace:
#   source("TUDOR_v2_Clean.R")
#   save.image("tudor_v2_workspace.RData")
#
# Or: source("TUDOR_v2_Clean.R") will be called if workspace not found.
#
# INPUT:  tudor_v2_workspace.RData (or TUDOR_v2_Clean.R + Publication.R)
# OUTPUT: wales_* tables and figures matching UKB pipeline naming convention
# ==============================================================================

set.seed(42)

# Null-coalescing operator (avoids rlang dependency)
`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
})

cat("\n")
cat("================================================================\n")
cat("  TUDOR PIPELINE: 10_wales_validation.R                        \n")
cat("  Wales Internal Validation — UKB-Aligned Output Format        \n")
cat("================================================================\n\n")

# ==============================================================================
# SECTION 0: SETUP & DATA LOADING
# ==============================================================================

cat("--- Section 0: Loading Wales Workspace ---\n")

# Try to load saved workspace first; fall back to sourcing Clean.R
# Detect script directory (works both via source() and interactive console)
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
wales_ws <- file.path(script_dir, "tudor_v2_workspace.RData")
if (!file.exists(wales_ws)) {
  wales_ws <- "C:/Users/nader/Downloads/tudor_v2_workspace.RData"
}

if (file.exists(wales_ws)) {
  load(wales_ws)
  cat("  Loaded workspace from:", wales_ws, "\n")
} else {
  cat("  Workspace not found. Sourcing TUDOR_v2_Clean.R...\n")
  source("C:/Users/nader/Downloads/TUDOR_v2_Clean.R")
}

# Verify critical objects exist
required_objects <- c("df", "model_df_v2", "cv_v2", "features_v2",
                      "en_pred_v2_te", "yte_v2", "roc_v2_ext",
                      "auc_v2_ext", "ci_v2_ext", "cen_v2", "sd_v2",
                      "features_base", "cv_v1", "en_pred_v1_te", "yte_v1",
                      "roc_v1_ext", "auc_v1_ext")
missing <- required_objects[!sapply(required_objects, exists)]
if (length(missing) > 0) {
  stop("Missing workspace objects: ", paste(missing, collapse = ", "),
       "\nPlease re-run TUDOR_v2_Clean.R and save workspace.")
}

# Output directories
OUTPUT_DIR <- file.path("C:/Users/nader/Downloads", "wales_pipeline_output")
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Build external test set dataframe (Relatives used for validation)
ext_test_rows_v2 <- as.integer(rownames(model_df_v2))[model_df_v2$I_Vs_R == 2]
te_m2_mask <- complete.cases(
  as.matrix(model_df_v2[model_df_v2$I_Vs_R == 2, features_v2]),
  model_df_v2$Positive1[model_df_v2$I_Vs_R == 2]
)
ext_df <- df[ext_test_rows_v2[te_m2_mask], ]

tudor_pred_ext <- as.numeric(en_pred_v2_te)
tudor_v1_pred_ext <- as.numeric(en_pred_v1_te)

cat(sprintf("  Wales Registry: %d total patients\n", nrow(df)))
cat(sprintf("  External test set (Relatives): %d patients\n", nrow(ext_df)))
cat(sprintf("  FH-positive: %d (%.1f%%)\n",
            sum(yte_v2 == 1), mean(yte_v2) * 100))
cat("\n")

# ==============================================================================
# SECTION 1: PRIMARY DISCRIMINATION (mirrors 02_validation_results.rds)
# ==============================================================================

cat("================================================================\n")
cat("SECTION 1: PRIMARY DISCRIMINATION ANALYSIS\n")
cat("  External Validation: Index -> Relatives (TRIPOD Type 2b)\n")
cat("================================================================\n\n")

# --- 1a. DLCN Scoring (clinical only, no DNA) ---
compute_dlcn_clinical <- function(ldl, tendon_xanth, corneal_lt40, ascvd_event,
                                   fam_hist_score, gender, ascvd_age) {
  score <- 0
  if (!is.na(ldl)) {
    if (ldl >= 8.5)      score <- score + 8
    else if (ldl >= 6.5) score <- score + 5
    else if (ldl >= 5.0) score <- score + 3
    else if (ldl >= 4.0) score <- score + 1
  }
  if (!is.na(tendon_xanth) && tendon_xanth == 1) score <- score + 6
  if (!is.na(corneal_lt40) && corneal_lt40 == 1) score <- score + 4
  if (!is.na(ascvd_event) && ascvd_event == 1) {
    if (!is.na(ascvd_age)) {
      if ((!is.na(gender) && gender == "M" && ascvd_age < 55) ||
          (!is.na(gender) && gender == "F" && ascvd_age < 60)) {
        score <- score + 2
      } else {
        score <- score + 1
      }
    } else {
      score <- score + 1
    }
  }
  if (!is.na(fam_hist_score) && fam_hist_score >= 1) {
    score <- score + min(as.integer(fam_hist_score), 2)
  }
  return(score)
}

# Simon Broome scoring
compute_sb_clinical <- function(ldl, tc, age, tendon_xanth, fam_hist_score) {
  has_high_chol <- FALSE
  if (!is.na(age) && !is.na(tc) && !is.na(ldl)) {
    if (age >= 16) has_high_chol <- (tc > 7.5 || ldl > 4.9)
    else           has_high_chol <- (tc > 6.7 || ldl > 4.0)
  }
  has_tx <- (!is.na(tendon_xanth) && tendon_xanth == 1)
  has_fh <- (!is.na(fam_hist_score) && fam_hist_score >= 1)
  if (has_high_chol && has_tx) return(3)
  if (has_high_chol && has_fh) return(2)
  if (has_high_chol) return(1)
  return(0)
}

# Apply DLCN to external test set
ext_df$DLCN_score <- mapply(compute_dlcn_clinical,
  ldl = ext_df$LDL_untreated,
  tendon_xanth = ifelse("tendon_xanth" %in% names(ext_df), ext_df$tendon_xanth,
                         ifelse("TendonXanthomata" %in% names(ext_df), ext_df$TendonXanthomata, NA)),
  corneal_lt40 = ext_df$corneal_less_40,
  ascvd_event  = ext_df$ascvd_combine,
  fam_hist_score = ext_df$DFamHist_score,
  gender       = ext_df$Gender,
  ascvd_age    = ext_df$min_ascvd_age
)

# Apply Simon Broome
ext_df$SB_score <- mapply(compute_sb_clinical,
  ldl = ext_df$LDL_untreated,
  tc  = ext_df$TC.1,
  age = ext_df$Age_at_LDL1,
  tendon_xanth = ifelse("tendon_xanth" %in% names(ext_df), ext_df$tendon_xanth,
                         ifelse("TendonXanthomata" %in% names(ext_df), ext_df$TendonXanthomata, NA)),
  fam_hist_score = ext_df$DFamHist_score
)

# --- 1b. ROC curves ---
roc_tudor   <- roc_v2_ext  # already computed
roc_v1      <- roc_v1_ext
roc_dlcn    <- roc(yte_v2, ext_df$DLCN_score, quiet = TRUE)
roc_sb      <- roc(yte_v2, ext_df$SB_score, quiet = TRUE)
roc_ldl     <- roc(yte_v2, ext_df$LDL_untreated, quiet = TRUE)
roc_trig    <- roc(yte_v2, ext_df$Trig_Filter, quiet = TRUE)

ci_tudor  <- ci_v2_ext
ci_v1     <- tryCatch(ci.auc(roc_v1, method = "delong"), error = function(e) c(NA, auc_v1_ext, NA))
ci_dlcn   <- ci.auc(roc_dlcn, method = "delong")
ci_sb     <- ci.auc(roc_sb, method = "delong")
ci_ldl    <- ci.auc(roc_ldl, method = "delong")
ci_trig   <- ci.auc(roc_trig, method = "delong")

cat("AREA UNDER ROC CURVE (AUC):\n\n")
cat(sprintf("%-30s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "TUDOR v2 (Elastic Net)", ci_tudor[2], ci_tudor[1], ci_tudor[3]))
cat(sprintf("%-30s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "TUDOR v1 (base features)", ci_v1[2], ci_v1[1], ci_v1[3]))
cat(sprintf("%-30s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "DLCN (clinical, no DNA)", ci_dlcn[2], ci_dlcn[1], ci_dlcn[3]))
cat(sprintf("%-30s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "Simon Broome (clinical)", ci_sb[2], ci_sb[1], ci_sb[3]))
cat(sprintf("%-30s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "LDL_untreated alone", ci_ldl[2], ci_ldl[1], ci_ldl[3]))
cat(sprintf("%-30s AUC = %.3f  (95%% CI: %.3f - %.3f)\n",
            "Trig Filter", ci_trig[2], ci_trig[1], ci_trig[3]))
cat("\n")

# --- 1c. DeLong pairwise tests ---
cat("DELONG PAIRWISE COMPARISONS:\n")
comparisons <- list(
  list("TUDOR v2 vs DLCN",        roc_tudor, roc_dlcn),
  list("TUDOR v2 vs Simon Broome", roc_tudor, roc_sb),
  list("TUDOR v2 vs TUDOR v1",    roc_tudor, roc_v1),
  list("TUDOR v2 vs LDL alone",   roc_tudor, roc_ldl),
  list("TUDOR v2 vs Trig Filter", roc_tudor, roc_trig)
)

delong_results <- list()
for (comp in comparisons) {
  dt <- tryCatch(roc.test(comp[[2]], comp[[3]], method = "delong"),
                 error = function(e) list(p.value = NA, statistic = NA))
  delta <- as.numeric(auc(comp[[2]])) - as.numeric(auc(comp[[3]]))
  cat(sprintf("  %-35s delta = %+.3f  p = %.2e %s\n",
              comp[[1]], delta, dt$p.value,
              ifelse(!is.na(dt$p.value) && dt$p.value < 0.05, "*", "")))
  delong_results[[comp[[1]]]] <- list(delta = delta, p = dt$p.value)
}
cat("\n")

# --- 1d. Youden optimal threshold ---
calc_metrics_at_threshold <- function(response, predictor, threshold, label) {
  pred_pos <- predictor >= threshold
  tp <- sum(pred_pos & response == 1)
  fp <- sum(pred_pos & response == 0)
  fn <- sum(!pred_pos & response == 1)
  tn <- sum(!pred_pos & response == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  ppv  <- tp / max(tp + fp, 1); npv  <- tn / max(tn + fn, 1)
  list(label = label, threshold = threshold,
       sens = sens, spec = spec, ppv = ppv, npv = npv,
       tp = tp, fp = fp, fn = fn, tn = tn)
}

youden_idx <- which.max(roc_tudor$sensitivities + roc_tudor$specificities - 1)
youden_thresh <- roc_tudor$thresholds[youden_idx]
m_youden <- calc_metrics_at_threshold(yte_v2, tudor_pred_ext, youden_thresh, "Youden")

cat(sprintf("YOUDEN OPTIMAL THRESHOLD: %.4f\n", youden_thresh))
cat(sprintf("  Sensitivity: %.1f%%\n", m_youden$sens * 100))
cat(sprintf("  Specificity: %.1f%%\n", m_youden$spec * 100))
cat(sprintf("  PPV:         %.1f%%\n", m_youden$ppv * 100))
cat(sprintf("  NPV:         %.1f%%\n\n", m_youden$npv * 100))


# ==============================================================================
# SECTION 2: SUBGROUP ANALYSES (mirrors val_results$subgroups)
# ==============================================================================

cat("================================================================\n")
cat("SECTION 2: SUBGROUP ANALYSES\n")
cat("================================================================\n\n")

subgroup_results <- data.frame()

run_subgroup <- function(mask, label) {
  sub_y <- yte_v2[mask]
  sub_p <- tudor_pred_ext[mask]
  sub_d <- ext_df$DLCN_score[mask]

  if (sum(sub_y == 1) < 10 || sum(sub_y == 0) < 10) {
    cat(sprintf("  %-30s SKIPPED (insufficient cases)\n", label))
    return(NULL)
  }

  r_tudor <- roc(sub_y, sub_p, quiet = TRUE)
  r_dlcn  <- roc(sub_y, sub_d, quiet = TRUE)
  ci_t <- ci.auc(r_tudor, method = "delong")
  ci_d <- ci.auc(r_dlcn, method = "delong")

  cat(sprintf("  %-30s TUDOR = %.3f [%.3f-%.3f]  DLCN = %.3f [%.3f-%.3f]  N=%d (FH=%d)\n",
              label, ci_t[2], ci_t[1], ci_t[3],
              ci_d[2], ci_d[1], ci_d[3],
              length(sub_y), sum(sub_y == 1)))

  data.frame(subgroup = label,
             tudor_auc = as.numeric(ci_t[2]),
             tudor_ci_lo = as.numeric(ci_t[1]),
             tudor_ci_hi = as.numeric(ci_t[3]),
             edlcn_auc = as.numeric(ci_d[2]),
             edlcn_ci_lo = as.numeric(ci_d[1]),
             edlcn_ci_hi = as.numeric(ci_d[3]),
             n = length(sub_y), n_fh = sum(sub_y == 1),
             stringsAsFactors = FALSE)
}

# Overall
subgroup_results <- rbind(subgroup_results,
  run_subgroup(rep(TRUE, length(yte_v2)), "Overall (Relatives)"))

# By Gender
subgroup_results <- rbind(subgroup_results,
  run_subgroup(ext_df$Gender == "M", "Male"),
  run_subgroup(ext_df$Gender == "F", "Female"))

# By Age
subgroup_results <- rbind(subgroup_results,
  run_subgroup(!is.na(ext_df$Age_at_LDL1) & ext_df$Age_at_LDL1 < 40, "Age < 40"),
  run_subgroup(!is.na(ext_df$Age_at_LDL1) & ext_df$Age_at_LDL1 >= 40 &
               ext_df$Age_at_LDL1 < 60, "Age 40-59"),
  run_subgroup(!is.na(ext_df$Age_at_LDL1) & ext_df$Age_at_LDL1 >= 60, "Age >= 60"))

# By Gene Type
subgroup_results <- rbind(subgroup_results,
  run_subgroup(!is.na(ext_df$Gene1) & ext_df$Gene1 == "LDLR", "LDLR mutation"),
  run_subgroup(!is.na(ext_df$Gene1) & ext_df$Gene1 == "APOB", "APOB mutation"))

# By Treatment Status
subgroup_results <- rbind(subgroup_results,
  run_subgroup(ext_df$treatment_status %in% c("untreated", "untreated_compliance"), "Untreated"),
  run_subgroup(ext_df$treatment_status == "treated", "On treatment"))

cat("\n")


# ==============================================================================
# SECTION 3: LANCET STATISTICS (mirrors 05_lancet_stats.rds)
# ==============================================================================

cat("================================================================\n")
cat("SECTION 3: LANCET / BMJ STATISTICS\n")
cat("================================================================\n\n")

n_boot <- 2000

# --- 3a. Calibrated probabilities ---
# TUDOR v2 predictions are already probability-scaled from Elastic Net
# For DLCN, calibrate via GLM
glm_dlcn <- glm(yte_v2 ~ ext_df$DLCN_score, family = binomial)
dlcn_calib_prob <- predict(glm_dlcn, type = "response")

cat("--- Calibration ---\n")
cat(sprintf("TUDOR prob range: [%.4f, %.4f]\n",
            min(tudor_pred_ext), max(tudor_pred_ext)))
cat(sprintf("DLCN calib prob range: [%.4f, %.4f]\n\n",
            min(dlcn_calib_prob), max(dlcn_calib_prob)))

# --- 3b. NRI (TUDOR v2 vs DLCN) ---
cat("--- Net Reclassification Index (NRI) ---\n")

# Risk categories adapted for Wales prevalence (~30-50%)
# Using 0.20, 0.50 as thresholds (Low <20%, Mod 20-50%, High >50%)
categorize_risk_wales <- function(prob) {
  cut(prob, breaks = c(-Inf, 0.20, 0.50, Inf),
      labels = c("Low", "Moderate", "High"))
}

calculate_nri <- function(outcome, prob_new, prob_old, thresholds = c(0.20, 0.50)) {
  cat_new <- as.numeric(cut(prob_new, breaks = c(-Inf, thresholds, Inf)))
  cat_old <- as.numeric(cut(prob_old, breaks = c(-Inf, thresholds, Inf)))
  events <- outcome == 1; nonevents <- outcome == 0
  up_events <- mean(cat_new[events] > cat_old[events])
  down_events <- mean(cat_new[events] < cat_old[events])
  nri_events <- up_events - down_events
  up_nonevents <- mean(cat_new[nonevents] > cat_old[nonevents])
  down_nonevents <- mean(cat_new[nonevents] < cat_old[nonevents])
  nri_nonevents <- down_nonevents - up_nonevents
  list(nri = nri_events + nri_nonevents,
       nri_events = nri_events, nri_nonevents = nri_nonevents,
       up_events = up_events, down_events = down_events,
       up_nonevents = up_nonevents, down_nonevents = down_nonevents)
}

nri_cat <- calculate_nri(yte_v2, tudor_pred_ext, dlcn_calib_prob)

cat(sprintf("Categorical NRI (TUDOR vs DLCN): %.3f\n", nri_cat$nri))
cat(sprintf("  Events NRI:     %+.3f\n", nri_cat$nri_events))
cat(sprintf("  Non-events NRI: %+.3f\n", nri_cat$nri_nonevents))

# Bootstrap CI
nri_boot <- numeric(n_boot)
for (b in seq_len(n_boot)) {
  idx <- sample(length(yte_v2), replace = TRUE)
  boot_nri <- calculate_nri(yte_v2[idx], tudor_pred_ext[idx], dlcn_calib_prob[idx])
  nri_boot[b] <- boot_nri$nri
}
nri_ci <- quantile(nri_boot, c(0.025, 0.975))
cat(sprintf("  Bootstrap 95%% CI: [%.3f, %.3f]\n", nri_ci[1], nri_ci[2]))

# Continuous NRI
events <- yte_v2 == 1; nonevents <- yte_v2 == 0
cnri_events <- mean(tudor_pred_ext[events] > dlcn_calib_prob[events]) -
               mean(tudor_pred_ext[events] < dlcn_calib_prob[events])
cnri_nonevents <- mean(tudor_pred_ext[nonevents] < dlcn_calib_prob[nonevents]) -
                  mean(tudor_pred_ext[nonevents] > dlcn_calib_prob[nonevents])
cnri_total <- cnri_events + cnri_nonevents
cat(sprintf("\nContinuous NRI: %.3f\n\n", cnri_total))

# --- 3c. IDI ---
cat("--- Integrated Discrimination Improvement (IDI) ---\n")
idi_events <- mean(tudor_pred_ext[events]) - mean(dlcn_calib_prob[events])
idi_nonevents <- mean(tudor_pred_ext[nonevents]) - mean(dlcn_calib_prob[nonevents])
idi <- idi_events - idi_nonevents

idi_boot <- numeric(n_boot)
for (b in seq_len(n_boot)) {
  idx <- sample(length(yte_v2), replace = TRUE)
  ev <- yte_v2[idx] == 1; nev <- yte_v2[idx] == 0
  ie <- mean(tudor_pred_ext[idx][ev]) - mean(dlcn_calib_prob[idx][ev])
  ine <- mean(tudor_pred_ext[idx][nev]) - mean(dlcn_calib_prob[idx][nev])
  idi_boot[b] <- ie - ine
}
idi_ci <- quantile(idi_boot, c(0.025, 0.975))

cat(sprintf("IDI (TUDOR vs DLCN): %.4f [%.4f, %.4f]\n\n", idi, idi_ci[1], idi_ci[2]))

# --- 3d. DCA ---
cat("--- Decision Curve Analysis (DCA) ---\n")
dca_thresholds <- seq(0.01, 0.50, by = 0.01)
prevalence <- mean(yte_v2)

dca_results <- data.frame(
  threshold = dca_thresholds,
  nb_tudor = NA_real_, nb_edlcn = NA_real_,
  nb_treat_all = NA_real_, nb_treat_none = 0
)

for (i in seq_along(dca_thresholds)) {
  pt <- dca_thresholds[i]; odds <- pt / (1 - pt)
  tudor_pos <- tudor_pred_ext >= pt
  tp_t <- sum(tudor_pos & yte_v2 == 1) / length(yte_v2)
  fp_t <- sum(tudor_pos & yte_v2 == 0) / length(yte_v2)
  dca_results$nb_tudor[i] <- tp_t - fp_t * odds

  dlcn_pos <- dlcn_calib_prob >= pt
  tp_d <- sum(dlcn_pos & yte_v2 == 1) / length(yte_v2)
  fp_d <- sum(dlcn_pos & yte_v2 == 0) / length(yte_v2)
  dca_results$nb_edlcn[i] <- tp_d - fp_d * odds

  dca_results$nb_treat_all[i] <- prevalence - (1 - prevalence) * odds
}

best_tudor <- sum(dca_results$nb_tudor > dca_results$nb_edlcn, na.rm = TRUE)
cat(sprintf("TUDOR has higher net benefit at %d/%d thresholds (%.0f%%)\n\n",
            best_tudor, length(dca_thresholds),
            best_tudor / length(dca_thresholds) * 100))

# --- 3e. Calibration ---
cat("--- Calibration Analysis ---\n")
tudor_logodds <- log(tudor_pred_ext / (1 - tudor_pred_ext + 1e-10))
calib_fit <- glm(yte_v2 ~ tudor_logodds, family = binomial)
calib_intercept <- coef(calib_fit)[1]
calib_slope <- coef(calib_fit)[2]
calib_ci_slope <- confint.default(calib_fit)

cat(sprintf("Calibration intercept: %.3f (ideal: 0)\n", calib_intercept))
cat(sprintf("Calibration slope:     %.3f (ideal: 1)\n", calib_slope))
cat(sprintf("  Slope 95%% CI: [%.3f, %.3f]\n",
            calib_ci_slope["tudor_logodds", 1], calib_ci_slope["tudor_logodds", 2]))

# H-L test
prob_decile <- cut(tudor_pred_ext,
                    breaks = quantile(tudor_pred_ext, probs = seq(0, 1, 0.1), na.rm = TRUE),
                    include.lowest = TRUE, labels = FALSE)

hl_df_data <- data.frame(decile = prob_decile, outcome = yte_v2, pred = tudor_pred_ext)
hl_table <- aggregate(cbind(n = outcome, observed = outcome, expected = pred) ~ decile,
                       data = hl_df_data,
                       FUN = function(x) c(length(x), sum(x), sum(x)))
# Recompute properly
hl_table <- do.call(rbind, lapply(sort(unique(prob_decile)), function(d) {
  mask <- prob_decile == d
  data.frame(prob_decile = d, n = sum(mask),
             observed = sum(yte_v2[mask]),
             expected = sum(tudor_pred_ext[mask]),
             obs_rate = mean(yte_v2[mask]),
             exp_rate = mean(tudor_pred_ext[mask]))
}))

hl_chi2 <- sum((hl_table$observed - hl_table$expected)^2 /
                (hl_table$expected * (1 - hl_table$exp_rate) + 0.001))
hl_df <- nrow(hl_table) - 2
hl_p <- 1 - pchisq(hl_chi2, df = hl_df)

cat(sprintf("Hosmer-Lemeshow chi2: %.2f (df=%d, p=%.3f)\n\n", hl_chi2, hl_df, hl_p))

# --- 3f. Brier scores ---
cat("--- Brier Score ---\n")
brier_tudor <- mean((tudor_pred_ext - yte_v2)^2)
brier_dlcn  <- mean((dlcn_calib_prob - yte_v2)^2)
brier_max   <- prevalence * (1 - prevalence)
scaled_brier_tudor <- 1 - brier_tudor / brier_max
scaled_brier_dlcn  <- 1 - brier_dlcn / brier_max

cat(sprintf("TUDOR Brier score:  %.4f (scaled: %.3f)\n", brier_tudor, scaled_brier_tudor))
cat(sprintf("DLCN Brier score:   %.4f (scaled: %.3f)\n\n", brier_dlcn, scaled_brier_dlcn))


# ==============================================================================
# SECTION 4: SENSITIVITY ANALYSES (mirrors 06_sensitivity_results.rds)
# ==============================================================================

cat("================================================================\n")
cat("SECTION 4: SENSITIVITY ANALYSES\n")
cat("================================================================\n\n")

sensitivity_results <- list()

auc_with_ci <- function(outcome, predictor, n_boot = 500) {
  valid <- !is.na(outcome) & !is.na(predictor)
  if (sum(outcome[valid] == 1) < 5) return(list(auc = NA, ci = c(NA, NA), n = sum(valid)))
  r <- roc(outcome[valid], predictor[valid], quiet = TRUE)
  ci <- ci.auc(r, conf.level = 0.95, method = "bootstrap",
               boot.n = n_boot, boot.stratified = TRUE)
  list(auc = as.numeric(ci[2]), ci = as.numeric(ci[c(1, 3)]),
       n = sum(valid), n_cases = sum(outcome[valid] == 1), roc = r)
}

# S_A: Without I_Vs_R (Wales-specific)
cat("--- S_A: Sensitivity WITHOUT I_Vs_R ---\n")
if (exists("auc_noIR_ext") && exists("ci_noIR_ext")) {
  cat(sprintf("  v2 WITH I_Vs_R:    %.4f [%.4f - %.4f]\n",
              auc_v2_ext, ci_v2_ext[1], ci_v2_ext[3]))
  cat(sprintf("  v2 WITHOUT I_Vs_R: %.4f [%.4f - %.4f]\n",
              auc_noIR_ext, ci_noIR_ext[1], ci_noIR_ext[3]))
  sensitivity_results$sa_without_ivsr <- list(
    with_ivsr = list(auc = as.numeric(auc_v2_ext), ci = as.numeric(ci_v2_ext[c(1,3)])),
    without_ivsr = list(auc = as.numeric(auc_noIR_ext), ci = as.numeric(ci_noIR_ext[c(1,3)]))
  )
} else {
  cat("  Objects not available (run Publication.R Task A first).\n")
}
cat("\n")

# S_B: v2 vs v1 comparison
cat("--- S_B: TUDOR v2 vs v1 ---\n")
cat(sprintf("  v1 External AUC: %.4f\n", as.numeric(auc_v1_ext)))
cat(sprintf("  v2 External AUC: %.4f\n", as.numeric(auc_v2_ext)))
cat(sprintf("  Delta v2 - v1: %+.4f\n\n", as.numeric(auc_v2_ext) - as.numeric(auc_v1_ext)))
sensitivity_results$sb_v2_vs_v1 <- list(
  v1_auc = as.numeric(auc_v1_ext), v2_auc = as.numeric(auc_v2_ext),
  delta = as.numeric(auc_v2_ext) - as.numeric(auc_v1_ext)
)

# S7: Prevalence-adjusted PPV/NPV
cat("--- S7: Prevalence-Adjusted PPV/NPV ---\n")
assumed_prevalences <- c(1/500, 1/250, 1/200, 1/100)
sens_at_youden <- m_youden$sens
spec_at_youden <- m_youden$spec

s7_results <- list()
cat(sprintf("%-15s | %8s | %8s\n", "Prevalence", "PPV", "NPV"))
cat(paste(rep("-", 40), collapse = ""), "\n")
for (prev in assumed_prevalences) {
  ppv <- (sens_at_youden * prev) / (sens_at_youden * prev + (1 - spec_at_youden) * (1 - prev))
  npv <- (spec_at_youden * (1 - prev)) / (spec_at_youden * (1 - prev) + (1 - sens_at_youden) * prev)
  label <- sprintf("1/%d", round(1/prev))
  cat(sprintf("%-15s | %7.1f%% | %7.1f%%\n", label, ppv * 100, npv * 100))
  s7_results[[label]] <- list(prevalence = prev, ppv = ppv, npv = npv)
}
sensitivity_results$s7_prevalence <- s7_results
cat("\n")

# S8: Sex interaction
cat("--- S8: Sex Interaction ---\n")
men_mask <- ext_df$Gender == "M"; women_mask <- ext_df$Gender == "F"
res_men <- auc_with_ci(yte_v2[men_mask], tudor_pred_ext[men_mask])
res_women <- auc_with_ci(yte_v2[women_mask], tudor_pred_ext[women_mask])
cat(sprintf("  Men:   AUC = %.3f [%.3f - %.3f] (N=%d)\n",
            res_men$auc, res_men$ci[1], res_men$ci[2], res_men$n))
cat(sprintf("  Women: AUC = %.3f [%.3f - %.3f] (N=%d)\n\n",
            res_women$auc, res_women$ci[1], res_women$ci[2], res_women$n))
sensitivity_results$s8_sex <- list(
  men = list(auc = res_men$auc, ci = res_men$ci),
  women = list(auc = res_women$auc, ci = res_women$ci)
)

# S10: Statin-free subgroup
cat("--- S10: Statin-Free Subgroup ---\n")
statin_free_mask <- ext_df$treatment_status %in% c("untreated", "untreated_compliance")
on_statin_mask <- ext_df$treatment_status == "treated"
res_free <- auc_with_ci(yte_v2[statin_free_mask], tudor_pred_ext[statin_free_mask])
res_statin <- auc_with_ci(yte_v2[on_statin_mask], tudor_pred_ext[on_statin_mask])

if (!is.na(res_free$auc)) {
  cat(sprintf("  Statin-Free: AUC = %.3f [%.3f - %.3f] (N=%d, Cases=%d)\n",
              res_free$auc, res_free$ci[1], res_free$ci[2], res_free$n, res_free$n_cases))
}
if (!is.na(res_statin$auc)) {
  cat(sprintf("  On Statin:   AUC = %.3f [%.3f - %.3f] (N=%d, Cases=%d)\n",
              res_statin$auc, res_statin$ci[1], res_statin$ci[2], res_statin$n, res_statin$n_cases))
}
sensitivity_results$s10_statin_free <- list(
  statin_free = list(auc = res_free$auc, ci = res_free$ci,
                     n = res_free$n, n_cases = res_free$n_cases),
  on_statin = list(auc = res_statin$auc, ci = res_statin$ci,
                   n = res_statin$n, n_cases = res_statin$n_cases)
)
cat("\n")


# ==============================================================================
# SECTION 5: TABLES & FIGURES (mirrors script 07 outputs)
# ==============================================================================

cat("================================================================\n")
cat("SECTION 5: PUBLICATION TABLES & FIGURES\n")
cat("================================================================\n\n")

# --- Helpers ---
median_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%.1f [%.1f - %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}
mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}
n_pct <- function(x, total) {
  if (length(x) == 0 || all(is.na(x))) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), mean(x, na.rm = TRUE) * 100)
}

# --- Table 1: Baseline Characteristics ---
cat("--- Table 1: Baseline Characteristics ---\n")

make_table1_row <- function(sub, label) {
  data.frame(
    Characteristic = label,
    N = nrow(sub),
    Age = mean_sd(sub$Age_at_LDL1),
    Female_pct = sprintf("%.1f%%", mean(sub$Gender == "F", na.rm = TRUE) * 100),
    LDL_mmol = median_iqr(sub$LDL.1),
    LDL_RW_mmol = median_iqr(sub$LDL_untreated),
    HDL_mmol = median_iqr(sub$HDL.1),
    Trig_mmol = median_iqr(sub$TRG.1),
    TC_mmol = median_iqr(sub$TC.1),
    BMI = mean_sd(sub$BMI_clean %||% sub$BMI),
    On_Statin_pct = sprintf("%.1f%%", mean(sub$treatment_status == "treated", na.rm = TRUE) * 100),
    FH_Cases = n_pct(sub$Positive1, nrow(sub)),
    Patient_Type = sprintf("Index: %d, Relative: %d",
                            sum(sub$I_Vs_R == 1, na.rm = TRUE),
                            sum(sub$I_Vs_R == 2, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

table1 <- rbind(
  make_table1_row(df, "All Patients"),
  make_table1_row(df[df$I_Vs_R == 1, ], "Index Patients"),
  make_table1_row(df[df$I_Vs_R == 2, ], "Relatives (External Test)"),
  make_table1_row(ext_df[yte_v2 == 1, ], "FH+ Relatives"),
  make_table1_row(ext_df[yte_v2 == 0, ], "FH- Relatives")
)

write.csv(table1, file.path(TABLE_DIR, "wales_table1_baseline.csv"), row.names = FALSE)
cat("  Saved: wales_table1_baseline.csv\n")

# --- Table 2: Subgroup AUCs ---
cat("--- Table 2: Subgroup AUCs ---\n")
if (nrow(subgroup_results) > 0) {
  write.csv(subgroup_results, file.path(TABLE_DIR, "wales_table2_subgroup_aucs.csv"),
            row.names = FALSE)
  cat("  Saved: wales_table2_subgroup_aucs.csv\n")
}

# --- Table 3: Head-to-Head ---
cat("--- Table 3: Head-to-Head Comparison ---\n")
table3 <- data.frame(
  Metric = c(
    "TUDOR v2 AUC", "TUDOR v1 AUC", "DLCN AUC (clinical)", "Simon Broome AUC",
    "LDL alone AUC", "Trig Filter AUC",
    "DeLong p: TUDOR vs DLCN", "DeLong p: TUDOR vs SB",
    "Categorical NRI (TUDOR vs DLCN)", "Continuous NRI",
    "IDI (TUDOR vs DLCN)",
    "Calibration Intercept", "Calibration Slope",
    "Hosmer-Lemeshow p",
    "Brier Score (TUDOR)", "Brier Score (DLCN)",
    "Scaled Brier (TUDOR)", "Scaled Brier (DLCN)"
  ),
  Value = c(
    sprintf("%.3f [%.3f - %.3f]", ci_tudor[2], ci_tudor[1], ci_tudor[3]),
    sprintf("%.3f [%.3f - %.3f]", ci_v1[2], ci_v1[1], ci_v1[3]),
    sprintf("%.3f [%.3f - %.3f]", ci_dlcn[2], ci_dlcn[1], ci_dlcn[3]),
    sprintf("%.3f [%.3f - %.3f]", ci_sb[2], ci_sb[1], ci_sb[3]),
    sprintf("%.3f [%.3f - %.3f]", ci_ldl[2], ci_ldl[1], ci_ldl[3]),
    sprintf("%.3f [%.3f - %.3f]", ci_trig[2], ci_trig[1], ci_trig[3]),
    sprintf("%.2e", delong_results[["TUDOR v2 vs DLCN"]]$p),
    sprintf("%.2e", delong_results[["TUDOR v2 vs Simon Broome"]]$p),
    sprintf("%.3f [%.3f, %.3f]", nri_cat$nri, nri_ci[1], nri_ci[2]),
    sprintf("%.3f", cnri_total),
    sprintf("%.4f [%.4f, %.4f]", idi, idi_ci[1], idi_ci[2]),
    sprintf("%.3f", calib_intercept),
    sprintf("%.3f [%.3f, %.3f]", calib_slope,
            calib_ci_slope["tudor_logodds", 1], calib_ci_slope["tudor_logodds", 2]),
    sprintf("%.3f", hl_p),
    sprintf("%.4f", brier_tudor),
    sprintf("%.4f", brier_dlcn),
    sprintf("%.3f", scaled_brier_tudor),
    sprintf("%.3f", scaled_brier_dlcn)
  ),
  stringsAsFactors = FALSE
)

write.csv(table3, file.path(TABLE_DIR, "wales_table3_lancet_stats.csv"), row.names = FALSE)
cat("  Saved: wales_table3_lancet_stats.csv\n")

# --- Table 4: Sensitivity Summary ---
cat("--- Table 4: Sensitivity Summary ---\n")
s_rows <- list()
if (!is.null(sensitivity_results$sa_without_ivsr)) {
  s <- sensitivity_results$sa_without_ivsr$without_ivsr
  s_rows <- c(s_rows, list(data.frame(Analysis = "S_A: Without I_Vs_R",
    AUC = sprintf("%.3f", s$auc),
    CI = sprintf("[%.3f - %.3f]", s$ci[1], s$ci[2]),
    N = NA, stringsAsFactors = FALSE)))
}
s_rows <- c(s_rows, list(data.frame(Analysis = "S_B: v2 vs v1 (delta)",
  AUC = sprintf("%+.4f", sensitivity_results$sb_v2_vs_v1$delta),
  CI = "-", N = NA, stringsAsFactors = FALSE)))
if (!is.null(sensitivity_results$s8_sex$men$auc) && !is.na(sensitivity_results$s8_sex$men$auc)) {
  s_rows <- c(s_rows, list(data.frame(Analysis = "S8: Men", AUC = sprintf("%.3f", sensitivity_results$s8_sex$men$auc),
    CI = sprintf("[%.3f - %.3f]", sensitivity_results$s8_sex$men$ci[1], sensitivity_results$s8_sex$men$ci[2]),
    N = NA, stringsAsFactors = FALSE)))
  s_rows <- c(s_rows, list(data.frame(Analysis = "S8: Women", AUC = sprintf("%.3f", sensitivity_results$s8_sex$women$auc),
    CI = sprintf("[%.3f - %.3f]", sensitivity_results$s8_sex$women$ci[1], sensitivity_results$s8_sex$women$ci[2]),
    N = NA, stringsAsFactors = FALSE)))
}
if (!is.null(sensitivity_results$s10_statin_free$statin_free$auc) &&
    !is.na(sensitivity_results$s10_statin_free$statin_free$auc)) {
  s <- sensitivity_results$s10_statin_free$statin_free
  s_rows <- c(s_rows, list(data.frame(Analysis = "S10: Statin-Free",
    AUC = sprintf("%.3f", s$auc),
    CI = sprintf("[%.3f - %.3f]", s$ci[1], s$ci[2]),
    N = s$n, stringsAsFactors = FALSE)))
}

if (length(s_rows) > 0) {
  table4 <- do.call(rbind, s_rows)
  write.csv(table4, file.path(TABLE_DIR, "wales_table4_sensitivity.csv"), row.names = FALSE)
  cat("  Saved: wales_table4_sensitivity.csv\n")
}

# --- Figure 1: ROC Curves ---
cat("\n--- Figure 1: ROC Curves ---\n")

roc_to_df <- function(roc_obj, label) {
  data.frame(Sensitivity = roc_obj$sensitivities,
             Specificity = 1 - roc_obj$specificities,
             Model = label, stringsAsFactors = FALSE)
}

plot_data <- rbind(
  roc_to_df(roc_tudor, sprintf("TUDOR v2 (AUC: %.3f)", auc(roc_tudor))),
  roc_to_df(roc_v1,    sprintf("TUDOR v1 (AUC: %.3f)", auc(roc_v1))),
  roc_to_df(roc_dlcn,  sprintf("DLCN (AUC: %.3f)", auc(roc_dlcn))),
  roc_to_df(roc_sb,    sprintf("Simon Broome (AUC: %.3f)", auc(roc_sb))),
  roc_to_df(roc_ldl,   sprintf("LDL alone (AUC: %.3f)", auc(roc_ldl)))
)

fig1 <- ggplot(plot_data, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#E41A1C", "#FF7F00", "#377EB8", "#984EA3", "#4DAF4A")) +
  labs(
    title = "Wales Internal Validation: TUDOR v2 vs Comparators",
    subtitle = sprintf("External Test Set (Relatives), N = %d, FH Cases = %d",
                        length(yte_v2), sum(yte_v2 == 1)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.65, 0.25),
        legend.background = element_rect(fill = "white", color = "grey80"),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank()) +
  coord_equal()

ggsave(file.path(FIG_DIR, "wales_figure1_roc_curves.pdf"), fig1, width = 7, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "wales_figure1_roc_curves.png"), fig1, width = 7, height = 7, dpi = 300)
cat("  Saved: wales_figure1_roc_curves.pdf/.png\n")

# --- Figure 2: DCA ---
cat("--- Figure 2: Decision Curve Analysis ---\n")

dca_long <- data.frame(
  Threshold = rep(dca_results$threshold, 4),
  Net_Benefit = c(dca_results$nb_tudor, dca_results$nb_edlcn,
                  dca_results$nb_treat_all, dca_results$nb_treat_none),
  Strategy = rep(c("TUDOR v2", "DLCN", "Test All", "Test None"), each = nrow(dca_results)),
  stringsAsFactors = FALSE
)

fig2 <- ggplot(dca_long, aes(x = Threshold * 100, y = Net_Benefit,
                              color = Strategy, linetype = Strategy)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "grey40", "grey60")) +
  scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted")) +
  labs(title = "Decision Curve Analysis: TUDOR v2 vs DLCN (Wales)",
       x = "Threshold Probability (%)", y = "Net Benefit",
       color = NULL, linetype = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill = "white", color = "grey80"),
        plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank()) +
  ylim(-0.05, NA)

ggsave(file.path(FIG_DIR, "wales_figure2_dca.pdf"), fig2, width = 8, height = 6, dpi = 300)
ggsave(file.path(FIG_DIR, "wales_figure2_dca.png"), fig2, width = 8, height = 6, dpi = 300)
cat("  Saved: wales_figure2_dca.pdf/.png\n")

# --- Figure 3: Calibration Plot ---
cat("--- Figure 3: Calibration Plot ---\n")

calib_plot_data <- hl_table

fig3 <- ggplot(calib_plot_data, aes(x = exp_rate, y = obs_rate)) +
  geom_point(aes(size = n), color = "#377EB8") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  geom_smooth(method = "loess", se = TRUE, color = "#E41A1C",
              fill = "#E41A1C", alpha = 0.2) +
  scale_size_continuous(range = c(3, 8), name = "N per decile") +
  labs(title = "Calibration Plot: TUDOR v2 (Wales Internal Validation)",
       subtitle = sprintf("Calibration slope = %.2f, H-L p = %.3f",
                            calib_slope, hl_p),
       x = "Mean Predicted Probability", y = "Observed Proportion") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "wales_figure3_calibration.pdf"), fig3, width = 7, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "wales_figure3_calibration.png"), fig3, width = 7, height = 7, dpi = 300)
cat("  Saved: wales_figure3_calibration.pdf/.png\n")

# --- Supplementary: Forest Plot ---
cat("--- Supplementary: Forest Plot ---\n")

if (nrow(subgroup_results) > 0 && any(!is.na(subgroup_results$tudor_auc))) {
  forest_data <- subgroup_results[!is.na(subgroup_results$tudor_auc), ]
  forest_data$subgroup <- factor(forest_data$subgroup,
                                  levels = rev(forest_data$subgroup))

  fig_forest <- ggplot(forest_data, aes(x = tudor_auc, y = subgroup)) +
    geom_point(size = 3, color = "#E41A1C") +
    geom_errorbarh(aes(xmin = tudor_ci_lo, xmax = tudor_ci_hi), height = 0.2,
                   color = "#E41A1C") +
    geom_vline(xintercept = as.numeric(auc_v2_ext), linetype = "dashed",
               color = "grey50") +
    labs(title = "TUDOR v2 AUC by Subgroup (Wales)",
         subtitle = sprintf("Dashed line = overall AUC (%.3f)", auc_v2_ext),
         x = "Area Under the ROC Curve (AUC)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    xlim(0.5, 1.0)

  ggsave(file.path(FIG_DIR, "wales_supp_forest_plot.pdf"), fig_forest,
         width = 8, height = 5, dpi = 300)
  ggsave(file.path(FIG_DIR, "wales_supp_forest_plot.png"), fig_forest,
         width = 8, height = 5, dpi = 300)
  cat("  Saved: wales_supp_forest_plot.pdf/.png\n")
}

# --- TRIPOD Checklist ---
cat("--- TRIPOD Checklist ---\n")

tripod <- data.frame(
  Item = c("1. Title", "2. Abstract", "3a. Background", "3b. Objectives",
           "4a. Study design", "4b. Setting",
           "5a. Participants", "5b. Eligibility", "6a. Outcome",
           "6b. Predictors", "7a. Sample size", "7b. Missing data",
           "8. Statistical methods", "9. Risk groups", "10a. Participant flow",
           "10b. Demographics", "11. Model development", "12. Model performance",
           "13a. Results", "13b. Sensitivity", "14. Discussion",
           "15. Limitations", "16. Implications"),
  Description = c(
    "Development and internal validation of TUDOR v2 FH diagnostic model",
    "AUC with 95% CI (bootstrap), DeLong comparisons vs DLCN & Simon Broome",
    "FH underdiagnosis; clinical scoring limitations; personalised calibration gap",
    "Develop TUDOR v2 with personalised pharmacogenomic calibration on Wales FH Registry",
    "Development + internal validation: Elastic Net, Index->Relatives (TRIPOD Type 2b)",
    "All Wales Familial Hypercholesterolaemia Registry",
    "Index patients (training) and Relatives (external test); genetic FH confirmed",
    "Patients with lipid panel data; genetic testing results available",
    "Genetic FH status (binary); defined by Mutation1 gene-level classification",
    "LDL_untreated (reverse-engineered), Trig_Filter, HDL, TRG, Age, Sex, pers_residual",
    "See Table 1 for Ns per group",
    "BMI: median imputation; treatment: coded from PASS drug codes",
    "Elastic Net (alpha=0.5), 5-fold CV, bootstrap CIs (2000), NRI, IDI, DCA",
    "FH+ vs FH- (genetic gold standard)",
    "See TUDOR_v2_Clean.R filtering log",
    "Table 1 from this script",
    "Two-stage hybrid: flat drug-specific + personalised calibration residual",
    "AUC, calibration slope/intercept, Brier score, NRI, IDI",
    "Table 2 subgroups, Table 3 Lancet stats",
    "Sensitivity analyses in Section 4",
    "Strengths: genetic gold standard, personalised calibration. Compare DLCN/SB.",
    "Literature-based statin factors, single-centre registry, compliance assumptions",
    "TUDOR v2 as pre-genetic screening tool; personalisation enables treated patient inclusion"
  ),
  Status = rep("Script-generated", 23),
  stringsAsFactors = FALSE
)

write.csv(tripod, file.path(TABLE_DIR, "wales_tripod_checklist.csv"), row.names = FALSE)
cat("  Saved: wales_tripod_checklist.csv\n")


# ==============================================================================
# SECTION 6: SIDE-BY-SIDE COMPARISON TABLE (Wales vs UKB)
# ==============================================================================

cat("\n================================================================\n")
cat("SECTION 6: SIDE-BY-SIDE COMPARISON (WALES vs UKB)\n")
cat("================================================================\n\n")

# Try to load UKB results if available
ukb_output_dir <- file.path("C:/Users/nader/Downloads", "tudor_pipeline_output")
ukb_val_file   <- file.path(ukb_output_dir, "02_validation_results.rds")
ukb_lancet_file <- file.path(ukb_output_dir, "05_lancet_stats.rds")

has_ukb_val    <- file.exists(ukb_val_file)
has_ukb_lancet <- file.exists(ukb_lancet_file)

if (has_ukb_val) {
  ukb_val <- readRDS(ukb_val_file)
  cat("  Loaded UKB validation results\n")
}
if (has_ukb_lancet) {
  ukb_lancet <- readRDS(ukb_lancet_file)
  cat("  Loaded UKB Lancet statistics\n")
}

comparison <- data.frame(
  Metric = c("TUDOR AUC", "Comparator AUC", "DeLong p",
             "Sensitivity (Youden)", "Specificity (Youden)",
             "PPV", "NPV",
             "Categorical NRI", "Continuous NRI", "IDI",
             "Calibration Slope", "Calibration Intercept",
             "H-L p-value", "Brier (TUDOR)", "Brier (Comparator)"),
  Wales_Value = c(
    sprintf("%.3f", as.numeric(auc_v2_ext)),
    sprintf("%.3f (DLCN)", as.numeric(auc(roc_dlcn))),
    sprintf("%.2e", delong_results[["TUDOR v2 vs DLCN"]]$p),
    sprintf("%.1f%%", m_youden$sens * 100),
    sprintf("%.1f%%", m_youden$spec * 100),
    sprintf("%.1f%%", m_youden$ppv * 100),
    sprintf("%.1f%%", m_youden$npv * 100),
    sprintf("%.3f", nri_cat$nri),
    sprintf("%.3f", cnri_total),
    sprintf("%.4f", idi),
    sprintf("%.3f", calib_slope),
    sprintf("%.3f", calib_intercept),
    sprintf("%.3f", hl_p),
    sprintf("%.4f", brier_tudor),
    sprintf("%.4f", brier_dlcn)
  ),
  Wales_CI = c(
    sprintf("[%.3f - %.3f]", ci_tudor[1], ci_tudor[3]),
    sprintf("[%.3f - %.3f]", ci_dlcn[1], ci_dlcn[3]),
    "-",
    "-", "-", "-", "-",
    sprintf("[%.3f, %.3f]", nri_ci[1], nri_ci[2]),
    "-",
    sprintf("[%.4f, %.4f]", idi_ci[1], idi_ci[2]),
    sprintf("[%.3f, %.3f]", calib_ci_slope["tudor_logodds", 1], calib_ci_slope["tudor_logodds", 2]),
    "-", "-", "-", "-"
  ),
  UKB_Value = rep("N/A", 15),
  UKB_CI = rep("N/A", 15),
  stringsAsFactors = FALSE
)

# Fill UKB values if available
if (has_ukb_val) {
  comparison$UKB_Value[1] <- sprintf("%.3f", ukb_val$primary$tudor_auc)
  comparison$UKB_CI[1]    <- sprintf("[%.3f - %.3f]", ukb_val$primary$tudor_ci[1], ukb_val$primary$tudor_ci[3])
  comparison$UKB_Value[2] <- sprintf("%.3f (eDLCN)", ukb_val$primary$edlcn_auc)
  comparison$UKB_CI[2]    <- sprintf("[%.3f - %.3f]", ukb_val$primary$edlcn_ci[1], ukb_val$primary$edlcn_ci[3])
  comparison$UKB_Value[4] <- sprintf("%.1f%%", ukb_val$primary$youden_sens * 100)
  comparison$UKB_Value[5] <- sprintf("%.1f%%", ukb_val$primary$youden_spec * 100)
}
if (has_ukb_lancet) {
  comparison$UKB_Value[8]  <- sprintf("%.3f", ukb_lancet$nri$categorical$nri)
  comparison$UKB_CI[8]     <- sprintf("[%.3f, %.3f]", ukb_lancet$nri$categorical_ci[1], ukb_lancet$nri$categorical_ci[2])
  comparison$UKB_Value[9]  <- sprintf("%.3f", ukb_lancet$nri$continuous$total)
  comparison$UKB_Value[10] <- sprintf("%.4f", ukb_lancet$idi$idi)
  comparison$UKB_CI[10]    <- sprintf("[%.4f, %.4f]", ukb_lancet$idi$ci[1], ukb_lancet$idi$ci[2])
  comparison$UKB_Value[11] <- sprintf("%.3f", ukb_lancet$calibration$slope)
  comparison$UKB_CI[11]    <- sprintf("[%.3f, %.3f]", ukb_lancet$calibration$slope_ci[1], ukb_lancet$calibration$slope_ci[2])
  comparison$UKB_Value[12] <- sprintf("%.3f", ukb_lancet$calibration$intercept)
  comparison$UKB_Value[13] <- sprintf("%.3f", ukb_lancet$calibration$hl_p)
  comparison$UKB_Value[14] <- sprintf("%.4f", ukb_lancet$brier$tudor)
  comparison$UKB_Value[15] <- sprintf("%.4f", ukb_lancet$brier$edlcn)
}

write.csv(comparison, file.path(TABLE_DIR, "wales_vs_ukb_comparison.csv"), row.names = FALSE)
cat("  Saved: wales_vs_ukb_comparison.csv\n")

# Print comparison
cat("\n  Side-by-Side Comparison:\n")
cat(sprintf("  %-25s | %-20s | %-20s\n", "Metric", "Wales", "UKB"))
cat(paste(rep("-", 70), collapse = ""), "\n")
for (i in seq_len(nrow(comparison))) {
  cat(sprintf("  %-25s | %-20s | %-20s\n",
              comparison$Metric[i], comparison$Wales_Value[i], comparison$UKB_Value[i]))
}
cat("\n")


# ==============================================================================
# SECTION 7: FH PATIENTS WITH LDL_UNTREATED < 4 mmol/L
# ==============================================================================

cat("================================================================\n")
cat("SECTION 7: FH PATIENTS WITH LDL_untreated < 4 mmol/L\n")
cat("  Clinically significant: missed by LDL-based screening (DLCN)\n")
cat("================================================================\n\n")

# --- Wales Analysis ---
fh_pos <- df[df$Positive1 == 1 & !is.na(df$LDL_untreated), ]
fh_low_ldl <- fh_pos[fh_pos$LDL_untreated < 4, ]

cat("=== WALES FH REGISTRY ===\n")
cat(sprintf("  Total FH+ patients with LDL_untreated data: %d\n", nrow(fh_pos)))
cat(sprintf("  FH+ with LDL_untreated < 4 mmol/L: %d (%.1f%%)\n",
            nrow(fh_low_ldl), nrow(fh_low_ldl) / nrow(fh_pos) * 100))

if (nrow(fh_low_ldl) > 0) {
  cat("\n  By Gene Type:\n")
  if ("Gene1" %in% names(fh_low_ldl)) {
    gene_tab <- table(fh_low_ldl$Gene1, useNA = "ifany")
    for (g in names(gene_tab)) {
      cat(sprintf("    %-15s: %d (%.1f%%)\n", ifelse(is.na(g), "Unknown", g),
                  gene_tab[g], gene_tab[g] / nrow(fh_low_ldl) * 100))
    }
  }

  cat("\n  By Treatment Status:\n")
  treat_tab <- table(fh_low_ldl$treatment_status, useNA = "ifany")
  for (t in names(treat_tab)) {
    cat(sprintf("    %-20s: %d (%.1f%%)\n", ifelse(is.na(t), "Unknown", t),
                treat_tab[t], treat_tab[t] / nrow(fh_low_ldl) * 100))
  }

  cat("\n  By Patient Type:\n")
  ivr_tab <- table(fh_low_ldl$I_Vs_R, useNA = "ifany")
  for (v in names(ivr_tab)) {
    label <- ifelse(v == "1", "Index", ifelse(v == "2", "Relative", "Unknown"))
    cat(sprintf("    %-15s: %d\n", label, ivr_tab[v]))
  }

  cat(sprintf("\n  Mean LDL_untreated in this group: %.2f mmol/L\n",
              mean(fh_low_ldl$LDL_untreated)))
  cat(sprintf("  Median LDL_untreated: %.2f mmol/L\n",
              median(fh_low_ldl$LDL_untreated)))
}

# --- TUDOR sensitivity for FH cases with LDL_untreated < 4 ---
# Among the external test set relatives
ext_fh_mask <- yte_v2 == 1
ext_low_ldl_mask <- ext_fh_mask & ext_df$LDL_untreated < 4

n_ext_fh_low <- sum(ext_low_ldl_mask, na.rm = TRUE)
cat(sprintf("\n  External test set FH+ with LDL_untreated < 4: %d\n", n_ext_fh_low))

if (n_ext_fh_low > 0) {
  tudor_pred_low <- tudor_pred_ext[ext_low_ldl_mask]
  tudor_correct_low <- sum(tudor_pred_low >= youden_thresh)
  cat(sprintf("  TUDOR correctly identifies at Youden threshold: %d / %d (%.1f%%)\n",
              tudor_correct_low, n_ext_fh_low,
              tudor_correct_low / n_ext_fh_low * 100))

  # DLCN for these patients
  dlcn_low <- ext_df$DLCN_score[ext_low_ldl_mask]
  dlcn_possible <- sum(dlcn_low >= 3)
  dlcn_probable <- sum(dlcn_low >= 6)
  cat(sprintf("  DLCN >= 3 (Possible): %d / %d (%.1f%%)\n",
              dlcn_possible, n_ext_fh_low, dlcn_possible / n_ext_fh_low * 100))
  cat(sprintf("  DLCN >= 6 (Probable): %d / %d (%.1f%%)\n",
              dlcn_probable, n_ext_fh_low, dlcn_probable / n_ext_fh_low * 100))
}
cat("\n")

# --- UKB Analysis (if data available) ---
ukb_rds <- file.path(ukb_output_dir, "tudor_analysis_ready.rds")
if (file.exists(ukb_rds)) {
  cat("=== UK BIOBANK ===\n")
  ukb_df <- readRDS(ukb_rds)
  if ("participant.eid" %in% names(ukb_df) && !"eid" %in% names(ukb_df)) {
    data.table::setnames(ukb_df, "participant.eid", "eid")
  }

  ukb_fh <- ukb_df[ukb_df$is_fh_genetic == 1 & !is.na(ukb_df$LDL_RW), ]
  ukb_fh_low <- ukb_fh[ukb_fh$LDL_RW < 4, ]

  cat(sprintf("  Total FH+ with LDL_RW data: %d\n", nrow(ukb_fh)))
  cat(sprintf("  FH+ with LDL_RW < 4 mmol/L: %d (%.1f%%)\n",
              nrow(ukb_fh_low), nrow(ukb_fh_low) / nrow(ukb_fh) * 100))

  if (nrow(ukb_fh_low) > 0) {
    cat("\n  By Statin Status:\n")
    statin_tab <- table(ukb_fh_low$statin_name, useNA = "ifany")
    for (s in names(statin_tab)) {
      cat(sprintf("    %-20s: %d\n", ifelse(is.na(s), "Unknown", s), statin_tab[s]))
    }

    cat(sprintf("\n  Mean LDL_RW: %.2f mmol/L\n", mean(ukb_fh_low$LDL_RW)))

    # TUDOR sensitivity for these patients
    if ("tudor_prob" %in% names(ukb_fh_low)) {
      ukb_hr_low <- ukb_fh_low[ukb_fh_low$cohort_high_risk == TRUE &
                                 !is.na(ukb_fh_low$tudor_prob), ]
      if (nrow(ukb_hr_low) > 0) {
        # Use UKB Youden threshold from val_results if available
        if (has_ukb_val) {
          ukb_youden <- ukb_val$primary$youden_threshold
          ukb_correct <- sum(ukb_hr_low$tudor_prob >= ukb_youden)
          cat(sprintf("  TUDOR correct (high-risk, Youden %.4f): %d / %d (%.1f%%)\n",
                      ukb_youden, ukb_correct, nrow(ukb_hr_low),
                      ukb_correct / nrow(ukb_hr_low) * 100))
        }
      }
    }
  }
  cat("\n")
  rm(ukb_df)  # Free memory
}

# --- Save FH LDL<4 summary ---
fh_low_summary <- data.frame(
  Dataset = c("Wales", "Wales"),
  Group = c("All FH+", "FH+ LDL_untreated < 4"),
  N = c(nrow(fh_pos), nrow(fh_low_ldl)),
  Pct = c(100, nrow(fh_low_ldl) / nrow(fh_pos) * 100),
  Mean_LDL_untreated = c(mean(fh_pos$LDL_untreated, na.rm = TRUE),
                          mean(fh_low_ldl$LDL_untreated, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

write.csv(fh_low_summary, file.path(TABLE_DIR, "wales_fh_low_ldl_analysis.csv"),
          row.names = FALSE)
cat("Saved: wales_fh_low_ldl_analysis.csv\n")


# ==============================================================================
# SAVE ALL RESULTS
# ==============================================================================

wales_val_results <- list(
  primary = list(
    cohort = "Wales FH Registry (Relatives - Internal Validation)",
    n = length(yte_v2), n_fh = sum(yte_v2 == 1),
    prevalence = mean(yte_v2),
    tudor_auc = as.numeric(auc_v2_ext),
    tudor_ci = as.numeric(ci_v2_ext),
    dlcn_auc = as.numeric(auc(roc_dlcn)),
    dlcn_ci = as.numeric(ci_dlcn),
    sb_auc = as.numeric(auc(roc_sb)),
    sb_ci = as.numeric(ci_sb),
    ldl_auc = as.numeric(auc(roc_ldl)),
    ldl_ci = as.numeric(ci_ldl),
    youden_threshold = youden_thresh,
    youden_sens = m_youden$sens,
    youden_spec = m_youden$spec
  ),
  subgroups = subgroup_results,
  delong = delong_results,
  roc_tudor = roc_tudor, roc_dlcn = roc_dlcn,
  roc_sb = roc_sb, roc_ldl = roc_ldl,
  calibration = hl_table
)

wales_lancet_results <- list(
  nri = list(categorical = nri_cat, categorical_ci = nri_ci,
             continuous = list(total = cnri_total, events = cnri_events,
                               nonevents = cnri_nonevents)),
  idi = list(idi = idi, idi_events = idi_events, idi_nonevents = idi_nonevents,
             ci = idi_ci),
  dca = dca_results,
  calibration = list(intercept = calib_intercept, slope = calib_slope,
                     slope_ci = calib_ci_slope["tudor_logodds", ],
                     hl_chi2 = hl_chi2, hl_df = hl_df, hl_p = hl_p,
                     decile_table = hl_table),
  brier = list(tudor = brier_tudor, dlcn = brier_dlcn,
               scaled_tudor = scaled_brier_tudor, scaled_dlcn = scaled_brier_dlcn)
)

saveRDS(wales_val_results, file.path(OUTPUT_DIR, "wales_02_validation_results.rds"))
saveRDS(wales_lancet_results, file.path(OUTPUT_DIR, "wales_05_lancet_stats.rds"))
saveRDS(sensitivity_results, file.path(OUTPUT_DIR, "wales_06_sensitivity_results.rds"))

cat("\n==========================================================\n")
cat(" OUTPUT SUMMARY\n")
cat("==========================================================\n")

tables <- list.files(TABLE_DIR, full.names = FALSE)
figures <- list.files(FIG_DIR, full.names = FALSE)
cat("Tables:\n"); for (t in tables) cat("  ", t, "\n")
cat("\nFigures:\n"); for (f in figures) cat("  ", f, "\n")
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))

cat("\n================================================================\n")
cat("  KEY NUMBERS FOR MANUSCRIPT (Wales Internal Validation):\n")
cat(sprintf("  TUDOR v2 AUC:  %.3f [%.3f - %.3f]\n",
            auc_v2_ext, ci_v2_ext[1], ci_v2_ext[3]))
cat(sprintf("  DLCN AUC:      %.3f [%.3f - %.3f]\n",
            auc(roc_dlcn), ci_dlcn[1], ci_dlcn[3]))
cat(sprintf("  NRI (cat):     %.3f [%.3f, %.3f]\n", nri_cat$nri, nri_ci[1], nri_ci[2]))
cat(sprintf("  IDI:           %.4f [%.4f, %.4f]\n", idi, idi_ci[1], idi_ci[2]))
cat(sprintf("  Calib slope:   %.3f\n", calib_slope))
cat(sprintf("  Brier (TUDOR): %.4f\n", brier_tudor))
cat(sprintf("  FH+ with LDL<4: %d / %d (%.1f%%)\n",
            nrow(fh_low_ldl), nrow(fh_pos), nrow(fh_low_ldl)/nrow(fh_pos)*100))
cat("================================================================\n")

cat("\n=== 10_wales_validation.R COMPLETE ===\n")
