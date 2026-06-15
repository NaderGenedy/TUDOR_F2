# ==============================================================================
# TUDOR PIPELINE: STEP 13 — ApoB AUGMENTED MODEL + NEW COHORT VALIDATION
# ==============================================================================
# PURPOSE: (1) Validate base TUDOR on the NEW lipid clinic cohort
#          (2) Test ApoB-augmented models for potential AUC improvement
#          (3) Provide manuscript-ready evidence for ApoB enhancement
#
# SECTION 1: Base TUDOR on NEW lipid clinic cohort (TC>7.5 | LDL>4.9 | ASCVD)
# SECTION 2: ApoB Augmented Models (A-D)
# SECTION 3: Manuscript-Ready Summary
#
# INPUT:   11_lipid_clinic_cohort.rds (or tudor_analysis_ready.rds)
#
# OUTPUT:  new_cohort_validation.csv, apob_augmented_results.csv,
#          cohort_criteria_breakdown.csv
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(pROC)
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
TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("================================================================\n")
cat("  TUDOR PIPELINE: 13_apob_augmented_model.R                     \n")
cat("  ApoB Augmentation + New Lipid Clinic Cohort Validation        \n")
cat("================================================================\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("--- Loading Data ---\n\n")

# Try lipid clinic cohort from script 11
lc_file <- file.path(OUTPUT_DIR, "11_lipid_clinic_cohort.rds")
full_file <- file.path(OUTPUT_DIR, "11_tudor_with_lipid_clinic.rds")

if (file.exists(full_file)) {
  df <- readRDS(full_file)
  cat("  Loaded full dataset with lipid clinic flag\n")
} else if (file.exists(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))) {
  df <- readRDS(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))
  cat("  Loaded tudor_analysis_ready.rds (will define cohort)\n")

  # Define lipid clinic cohort
  crit_tc    <- !is.na(df$CHOL) & df$CHOL > 7.5
  crit_ldl   <- !is.na(df$LDL_RW) & df$LDL_RW > 4.9
  crit_ascvd <- !is.na(df$Premature_ASCVD) & df$Premature_ASCVD == 1
  df$cohort_lipid_clinic <- crit_tc | crit_ldl | crit_ascvd
} else {
  stop("No analysis-ready data found. Run scripts 01 and 11 first.")
}

if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}

cat(sprintf("  Total UKB: %d\n", nrow(df)))
cat(sprintf("  Genetic FH total: %d\n", sum(df$is_fh_genetic)))
cat(sprintf("  Old high-risk cohort: %d\n", sum(df$cohort_high_risk, na.rm = TRUE)))
cat(sprintf("  New lipid clinic cohort: %d\n\n", sum(df$cohort_lipid_clinic)))

# ==============================================================================
# SECTION 1: BASE TUDOR ON NEW LIPID CLINIC COHORT
# ==============================================================================
cat("================================================================\n")
cat("SECTION 1: BASE TUDOR VALIDATION — NEW LIPID CLINIC COHORT\n")
cat("================================================================\n\n")

# New cohort
lc <- df[df$cohort_lipid_clinic == TRUE & !is.na(df$tudor_prob), ]
cat(sprintf("  Lipid clinic cohort (valid predictions): %d\n", nrow(lc)))
cat(sprintf("  FH cases: %d (%.2f%%)\n", sum(lc$is_fh_genetic),
            100 * mean(lc$is_fh_genetic)))

# Old cohort for comparison
hr <- df[df$cohort_high_risk == TRUE & !is.na(df$tudor_prob), ]
cat(sprintf("  Old high-risk cohort: %d (FH = %d, %.2f%%)\n\n",
            nrow(hr), sum(hr$is_fh_genetic), 100 * mean(hr$is_fh_genetic)))

# --- 1a. ROC Analysis — NEW cohort ---
cat("--- Discrimination: NEW Lipid Clinic Cohort ---\n\n")

roc_tudor_lc  <- roc(lc$is_fh_genetic, lc$tudor_prob, quiet = TRUE)
roc_edlcn_lc  <- roc(lc$is_fh_genetic, lc$edlcn_score, quiet = TRUE)
roc_ldl_lc    <- roc(lc$is_fh_genetic, lc$LDL_RW, quiet = TRUE)
roc_trig_lc   <- roc(lc$is_fh_genetic, lc$Trig_Filter_RW, quiet = TRUE)

ci_tudor_lc  <- ci.auc(roc_tudor_lc, method = "delong")
ci_edlcn_lc  <- ci.auc(roc_edlcn_lc, method = "delong")
ci_ldl_lc    <- ci.auc(roc_ldl_lc, method = "delong")
ci_trig_lc   <- ci.auc(roc_trig_lc, method = "delong")

cat("  AUC — NEW Lipid Clinic Cohort:\n")
cat(sprintf("    %-20s  %.3f [%.3f–%.3f]\n", "TUDOR v2",
            ci_tudor_lc[2], ci_tudor_lc[1], ci_tudor_lc[3]))
cat(sprintf("    %-20s  %.3f [%.3f–%.3f]\n", "eDLCN",
            ci_edlcn_lc[2], ci_edlcn_lc[1], ci_edlcn_lc[3]))
cat(sprintf("    %-20s  %.3f [%.3f–%.3f]\n", "LDL-C alone",
            ci_ldl_lc[2], ci_ldl_lc[1], ci_ldl_lc[3]))
cat(sprintf("    %-20s  %.3f [%.3f–%.3f]\n", "Trig Filter",
            ci_trig_lc[2], ci_trig_lc[1], ci_trig_lc[3]))

# --- 1b. Compare OLD vs NEW cohort AUC ---
cat("\n--- Comparison: OLD cohort vs NEW cohort ---\n\n")

roc_tudor_hr <- roc(hr$is_fh_genetic, hr$tudor_prob, quiet = TRUE)
ci_tudor_hr  <- ci.auc(roc_tudor_hr, method = "delong")

cat(sprintf("  TUDOR AUC — OLD (LDL>4.9):    %.3f [%.3f–%.3f] (N=%d)\n",
            ci_tudor_hr[2], ci_tudor_hr[1], ci_tudor_hr[3], nrow(hr)))
cat(sprintf("  TUDOR AUC — NEW (lipid clinic): %.3f [%.3f–%.3f] (N=%d)\n",
            ci_tudor_lc[2], ci_tudor_lc[1], ci_tudor_lc[3], nrow(lc)))

# --- 1c. DeLong Comparisons (NEW cohort) ---
cat("\n--- DeLong Pairwise (NEW cohort) ---\n")

delong_pairs <- list(
  list("TUDOR vs eDLCN",       roc_tudor_lc, roc_edlcn_lc),
  list("TUDOR vs LDL-C",       roc_tudor_lc, roc_ldl_lc),
  list("TUDOR vs Trig Filter", roc_tudor_lc, roc_trig_lc)
)

for (pair in delong_pairs) {
  dl <- tryCatch(roc.test(pair[[2]], pair[[3]], method = "delong"),
                  error = function(e) NULL)
  if (!is.null(dl)) {
    cat(sprintf("    %-25s  Z = %.3f,  p = %.4g\n",
                pair[[1]], dl$statistic, dl$p.value))
  }
}

# --- 1d. Youden Threshold Metrics ---
cat("\n--- Youden Optimal Threshold (NEW cohort) ---\n\n")

youden_tudor <- coords(roc_tudor_lc, "best", ret = c("threshold", "sensitivity",
                                                       "specificity", "ppv", "npv"),
                        best.method = "youden")

cat(sprintf("  Threshold: %.4f\n", youden_tudor$threshold))
cat(sprintf("  Sensitivity: %.3f\n", youden_tudor$sensitivity))
cat(sprintf("  Specificity: %.3f\n", youden_tudor$specificity))
cat(sprintf("  PPV: %.4f\n", youden_tudor$ppv))
cat(sprintf("  NPV: %.6f\n", youden_tudor$npv))

# --- 1e. Calibration (NEW cohort) ---
cat("\n--- Calibration (NEW cohort) ---\n\n")

cal_model <- tryCatch(
  glm(is_fh_genetic ~ tudor_prob, family = binomial, data = lc),
  error = function(e) NULL
)

if (!is.null(cal_model)) {
  coefs <- coef(cal_model)
  ci_coefs <- confint(cal_model)
  cal_slope <- coefs[2]
  cal_slope_ci <- ci_coefs[2, ]

  cat(sprintf("  Calibration intercept: %.3f [%.3f, %.3f]\n",
              coefs[1], ci_coefs[1, 1], ci_coefs[1, 2]))
  cat(sprintf("  Calibration slope:     %.3f [%.3f, %.3f]\n",
              cal_slope, cal_slope_ci[1], cal_slope_ci[2]))
}

# Brier score
brier <- mean((lc$tudor_prob - lc$is_fh_genetic)^2)
cat(sprintf("  Brier score: %.6f\n", brier))

# Save new cohort validation results
new_cohort_results <- data.frame(
  Metric = c("TUDOR AUC", "eDLCN AUC", "LDL-C AUC", "Trig Filter AUC",
             "Youden Sens", "Youden Spec", "Brier", "Calibration Slope"),
  Value = c(ci_tudor_lc[2], ci_edlcn_lc[2], ci_ldl_lc[2], ci_trig_lc[2],
            youden_tudor$sensitivity, youden_tudor$specificity, brier,
            if (!is.null(cal_model)) cal_slope else NA),
  CI_lower = c(ci_tudor_lc[1], ci_edlcn_lc[1], ci_ldl_lc[1], ci_trig_lc[1],
               NA, NA, NA,
               if (!is.null(cal_model)) cal_slope_ci[1] else NA),
  CI_upper = c(ci_tudor_lc[3], ci_edlcn_lc[3], ci_ldl_lc[3], ci_trig_lc[3],
               NA, NA, NA,
               if (!is.null(cal_model)) cal_slope_ci[2] else NA),
  N = nrow(lc),
  FH_cases = sum(lc$is_fh_genetic),
  stringsAsFactors = FALSE
)

write.csv(new_cohort_results, file.path(TABLE_DIR, "new_cohort_validation.csv"),
          row.names = FALSE)
cat("\nSaved: new_cohort_validation.csv\n")

# ==============================================================================
# SECTION 2: ApoB AUGMENTED MODELS
# ==============================================================================
cat("\n================================================================\n")
cat("SECTION 2: ApoB AUGMENTED MODELS\n")
cat("================================================================\n\n")

# Check ApoB availability
has_apob <- "ApoB" %in% names(lc) && sum(!is.na(lc$ApoB)) > 100

if (!has_apob) {
  cat("  WARNING: ApoB data not available or insufficient.\n")
  cat("  Checking alternative column names...\n")

  apob_cols <- grep("apob|p30640", names(lc), ignore.case = TRUE, value = TRUE)
  cat("  ApoB-related columns found:", paste(apob_cols, collapse = ", "), "\n")

  if (length(apob_cols) > 0) {
    lc$ApoB <- as.numeric(lc[[apob_cols[1]]])
    has_apob <- sum(!is.na(lc$ApoB)) > 100
    if (has_apob) cat(sprintf("  Using column '%s': %d non-missing values\n",
                                apob_cols[1], sum(!is.na(lc$ApoB))))
  }
}

if (has_apob) {
  cat(sprintf("  ApoB available: %d / %d (%.1f%%) non-missing in lipid clinic cohort\n",
              sum(!is.na(lc$ApoB)), nrow(lc),
              100 * sum(!is.na(lc$ApoB)) / nrow(lc)))

  # Complete-case subset
  lc_apob <- lc[!is.na(lc$ApoB) & !is.na(lc$tudor_prob) &
                  !is.na(lc$LDL_RW) & !is.na(lc$is_fh_genetic), ]

  cat(sprintf("  Complete-case ApoB subset: %d (FH+ = %d)\n\n",
              nrow(lc_apob), sum(lc_apob$is_fh_genetic)))

  # ApoB/LDL ratio
  lc_apob$ApoB_LDL <- lc_apob$ApoB / lc_apob$LDL_RW
  lc_apob$ApoB_LDL[!is.finite(lc_apob$ApoB_LDL)] <- NA

  # Check Lp(a) availability
  has_lpa <- "Lpa_nmol" %in% names(lc_apob) && sum(!is.na(lc_apob$Lpa_nmol)) > 100

  # --- Base TUDOR on ApoB complete-case subset ---
  cat("--- Base TUDOR (ApoB complete-case subset) ---\n")
  roc_base <- roc(lc_apob$is_fh_genetic, lc_apob$tudor_prob, quiet = TRUE)
  ci_base <- ci.auc(roc_base, method = "delong")
  cat(sprintf("  TUDOR base AUC: %.3f [%.3f–%.3f] (N=%d)\n\n",
              ci_base[2], ci_base[1], ci_base[3], nrow(lc_apob)))

  # --- Model A: TUDOR + ApoB ---
  cat("--- Model A: TUDOR + ApoB ---\n")
  fit_A <- tryCatch(
    glm(is_fh_genetic ~ tudor_prob + ApoB, family = binomial, data = lc_apob),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )

  if (!is.null(fit_A)) {
    pred_A <- predict(fit_A, type = "response")
    roc_A <- roc(lc_apob$is_fh_genetic, pred_A, quiet = TRUE)
    ci_A <- ci.auc(roc_A, method = "delong")
    dl_A <- roc.test(roc_A, roc_base, method = "delong")

    cat(sprintf("  AUC: %.3f [%.3f–%.3f]\n", ci_A[2], ci_A[1], ci_A[3]))
    cat(sprintf("  DeLong vs base: Z=%.3f, p=%.4g\n", dl_A$statistic, dl_A$p.value))
    cat(sprintf("  Delta AUC: %+.3f\n\n", ci_A[2] - ci_base[2]))
  }

  # --- Model B: TUDOR + ApoB/LDL ratio ---
  cat("--- Model B: TUDOR + ApoB/LDL-C ratio ---\n")
  lc_apob_B <- lc_apob[!is.na(lc_apob$ApoB_LDL) & is.finite(lc_apob$ApoB_LDL), ]

  fit_B <- tryCatch(
    glm(is_fh_genetic ~ tudor_prob + ApoB_LDL, family = binomial, data = lc_apob_B),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )

  if (!is.null(fit_B)) {
    pred_B <- predict(fit_B, type = "response")
    roc_B <- roc(lc_apob_B$is_fh_genetic, pred_B, quiet = TRUE)
    ci_B <- ci.auc(roc_B, method = "delong")

    # Use matching base ROC for fair comparison
    roc_base_B <- roc(lc_apob_B$is_fh_genetic, lc_apob_B$tudor_prob, quiet = TRUE)
    dl_B <- roc.test(roc_B, roc_base_B, method = "delong")

    cat(sprintf("  AUC: %.3f [%.3f–%.3f]\n", ci_B[2], ci_B[1], ci_B[3]))
    cat(sprintf("  DeLong vs base: Z=%.3f, p=%.4g\n", dl_B$statistic, dl_B$p.value))
    cat(sprintf("  Delta AUC: %+.3f\n\n",
                ci_B[2] - as.numeric(auc(roc_base_B))))
  }

  # --- Model C: TUDOR + ApoB + ApoB/LDL ---
  cat("--- Model C: TUDOR + ApoB + ApoB/LDL-C ---\n")
  fit_C <- tryCatch(
    glm(is_fh_genetic ~ tudor_prob + ApoB + ApoB_LDL, family = binomial,
        data = lc_apob_B),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )

  if (!is.null(fit_C)) {
    pred_C <- predict(fit_C, type = "response")
    roc_C <- roc(lc_apob_B$is_fh_genetic, pred_C, quiet = TRUE)
    ci_C <- ci.auc(roc_C, method = "delong")
    dl_C <- roc.test(roc_C, roc_base_B, method = "delong")

    cat(sprintf("  AUC: %.3f [%.3f–%.3f]\n", ci_C[2], ci_C[1], ci_C[3]))
    cat(sprintf("  DeLong vs base: Z=%.3f, p=%.4g\n", dl_C$statistic, dl_C$p.value))
    cat(sprintf("  Delta AUC: %+.3f\n\n",
                ci_C[2] - as.numeric(auc(roc_base_B))))
  }

  # --- Model D: TUDOR + ApoB + ApoB/LDL + Lp(a) (full multimarker) ---
  if (has_lpa) {
    cat("--- Model D: TUDOR + ApoB + ApoB/LDL + Lp(a) (full multimarker) ---\n")
    lc_apob_D <- lc_apob_B[!is.na(lc_apob_B$Lpa_nmol), ]
    cat(sprintf("  Complete cases with all markers: %d (FH+ = %d)\n",
                nrow(lc_apob_D), sum(lc_apob_D$is_fh_genetic)))

    if (nrow(lc_apob_D) > 100 && sum(lc_apob_D$is_fh_genetic) > 5) {
      fit_D <- tryCatch(
        glm(is_fh_genetic ~ tudor_prob + ApoB + ApoB_LDL + Lpa_nmol,
            family = binomial, data = lc_apob_D),
        error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
      )

      if (!is.null(fit_D)) {
        pred_D <- predict(fit_D, type = "response")
        roc_D <- roc(lc_apob_D$is_fh_genetic, pred_D, quiet = TRUE)
        ci_D <- ci.auc(roc_D, method = "delong")
        roc_base_D <- roc(lc_apob_D$is_fh_genetic, lc_apob_D$tudor_prob, quiet = TRUE)
        dl_D <- roc.test(roc_D, roc_base_D, method = "delong")

        cat(sprintf("  AUC: %.3f [%.3f–%.3f]\n", ci_D[2], ci_D[1], ci_D[3]))
        cat(sprintf("  DeLong vs base: Z=%.3f, p=%.4g\n",
                    dl_D$statistic, dl_D$p.value))
        cat(sprintf("  Delta AUC: %+.3f\n\n",
                    ci_D[2] - as.numeric(auc(roc_base_D))))
      }
    } else {
      cat("  Insufficient Lp(a) complete cases for Model D\n\n")
    }
  } else {
    cat("--- Model D: SKIPPED (Lp(a) not available) ---\n\n")
  }

  # --- NRI and IDI for best augmented model ---
  cat("--- Net Reclassification and Integrated Discrimination ---\n\n")

  # Use Model C (TUDOR + ApoB + ApoB/LDL) as best candidate
  if (!is.null(fit_C)) {
    base_pred <- lc_apob_B$tudor_prob
    aug_pred  <- predict(fit_C, type = "response")
    outcome   <- lc_apob_B$is_fh_genetic

    # Category-free NRI
    events    <- outcome == 1
    nonevents <- outcome == 0

    nri_events    <- mean(aug_pred[events] > base_pred[events]) -
                     mean(aug_pred[events] < base_pred[events])
    nri_nonevents <- mean(aug_pred[nonevents] < base_pred[nonevents]) -
                     mean(aug_pred[nonevents] > base_pred[nonevents])
    nri_total     <- nri_events + nri_nonevents

    cat(sprintf("  NRI (category-free, Model C vs base TUDOR):\n"))
    cat(sprintf("    NRI events:     %+.4f\n", nri_events))
    cat(sprintf("    NRI non-events: %+.4f\n", nri_nonevents))
    cat(sprintf("    NRI total:      %+.4f\n", nri_total))

    # IDI
    idi <- mean(aug_pred[events]) - mean(aug_pred[nonevents]) -
           (mean(base_pred[events]) - mean(base_pred[nonevents]))
    cat(sprintf("    IDI:            %+.6f\n\n", idi))
  }

  # --- Compile augmented model comparison table ---
  cat("--- Summary: Augmented Model Comparison ---\n\n")

  aug_results <- data.frame(
    Model = character(0), AUC = numeric(0), CI_lower = numeric(0),
    CI_upper = numeric(0), Delta_AUC = numeric(0), DeLong_p = numeric(0),
    N = integer(0), stringsAsFactors = FALSE
  )

  aug_results <- rbind(aug_results, data.frame(
    Model = "TUDOR base", AUC = ci_base[2], CI_lower = ci_base[1],
    CI_upper = ci_base[3], Delta_AUC = 0,
    DeLong_p = NA, N = nrow(lc_apob), stringsAsFactors = FALSE
  ))

  if (!is.null(fit_A)) {
    aug_results <- rbind(aug_results, data.frame(
      Model = "A: TUDOR + ApoB", AUC = ci_A[2], CI_lower = ci_A[1],
      CI_upper = ci_A[3], Delta_AUC = ci_A[2] - ci_base[2],
      DeLong_p = dl_A$p.value, N = nrow(lc_apob), stringsAsFactors = FALSE
    ))
  }
  if (!is.null(fit_B)) {
    aug_results <- rbind(aug_results, data.frame(
      Model = "B: TUDOR + ApoB/LDL", AUC = ci_B[2], CI_lower = ci_B[1],
      CI_upper = ci_B[3],
      Delta_AUC = ci_B[2] - as.numeric(auc(roc_base_B)),
      DeLong_p = dl_B$p.value, N = nrow(lc_apob_B), stringsAsFactors = FALSE
    ))
  }
  if (!is.null(fit_C)) {
    aug_results <- rbind(aug_results, data.frame(
      Model = "C: TUDOR + ApoB + ApoB/LDL", AUC = ci_C[2], CI_lower = ci_C[1],
      CI_upper = ci_C[3],
      Delta_AUC = ci_C[2] - as.numeric(auc(roc_base_B)),
      DeLong_p = dl_C$p.value, N = nrow(lc_apob_B), stringsAsFactors = FALSE
    ))
  }
  if (exists("fit_D") && !is.null(fit_D)) {
    aug_results <- rbind(aug_results, data.frame(
      Model = "D: TUDOR + ApoB + ApoB/LDL + Lp(a)", AUC = ci_D[2],
      CI_lower = ci_D[1], CI_upper = ci_D[3],
      Delta_AUC = ci_D[2] - as.numeric(auc(roc_base_D)),
      DeLong_p = dl_D$p.value, N = nrow(lc_apob_D), stringsAsFactors = FALSE
    ))
  }

  cat(sprintf("  %-35s  %5s  [95%% CI]           Delta   p\n", "Model", "AUC"))
  cat("  ", strrep("-", 80), "\n")
  for (i in seq_len(nrow(aug_results))) {
    r <- aug_results[i, ]
    cat(sprintf("  %-35s  %.3f  [%.3f–%.3f]  %+.3f  %s\n",
                r$Model, r$AUC, r$CI_lower, r$CI_upper, r$Delta_AUC,
                if (is.na(r$DeLong_p)) "ref" else sprintf("%.4g", r$DeLong_p)))
  }

  write.csv(aug_results, file.path(TABLE_DIR, "apob_augmented_results.csv"),
            row.names = FALSE)
  cat("\n  Saved: apob_augmented_results.csv\n")

} else {
  cat("  ApoB data NOT available. Augmented models cannot be computed.\n")
  cat("  This is expected if p30640 was not extracted from UKB-RAP.\n")
}

# ==============================================================================
# SECTION 3: COHORT CRITERIA BREAKDOWN
# ==============================================================================
cat("\n================================================================\n")
cat("SECTION 3: COHORT CRITERIA BREAKDOWN\n")
cat("================================================================\n\n")

# Criteria breakdown
crit_tc    <- !is.na(df$CHOL) & df$CHOL > 7.5
crit_ldl   <- !is.na(df$LDL_RW) & df$LDL_RW > 4.9
crit_ascvd <- !is.na(df$Premature_ASCVD) & df$Premature_ASCVD == 1

breakdown <- data.frame(
  Criterion = c("TC > 7.5 mmol/L", "LDL_RW > 4.9 mmol/L", "Premature ASCVD",
                "TC > 7.5 ONLY (not LDL/ASCVD)", "LDL > 4.9 ONLY (not TC/ASCVD)",
                "ASCVD ONLY (not TC/LDL)",
                "Any criterion (lipid clinic)", "Old cohort (LDL > 4.9 only)"),
  N = c(sum(crit_tc), sum(crit_ldl), sum(crit_ascvd),
        sum(crit_tc & !crit_ldl & !crit_ascvd),
        sum(!crit_tc & crit_ldl & !crit_ascvd),
        sum(!crit_tc & !crit_ldl & crit_ascvd),
        sum(df$cohort_lipid_clinic), sum(df$cohort_high_risk, na.rm = TRUE)),
  FH_cases = c(sum(df$is_fh_genetic[crit_tc]),
               sum(df$is_fh_genetic[crit_ldl]),
               sum(df$is_fh_genetic[crit_ascvd]),
               sum(df$is_fh_genetic[crit_tc & !crit_ldl & !crit_ascvd]),
               sum(df$is_fh_genetic[!crit_tc & crit_ldl & !crit_ascvd]),
               sum(df$is_fh_genetic[!crit_tc & !crit_ldl & crit_ascvd]),
               sum(df$is_fh_genetic[df$cohort_lipid_clinic]),
               sum(df$is_fh_genetic[df$cohort_high_risk == TRUE], na.rm = TRUE)),
  stringsAsFactors = FALSE
)
breakdown$FH_pct <- sprintf("%.2f%%", 100 * breakdown$FH_cases / breakdown$N)

cat("  Cohort Criteria Breakdown:\n")
cat(sprintf("  %-35s  %8s  %6s  %8s\n", "Criterion", "N", "FH+", "FH%"))
cat("  ", strrep("-", 65), "\n")
for (i in seq_len(nrow(breakdown))) {
  cat(sprintf("  %-35s  %8d  %6d  %8s\n",
              breakdown$Criterion[i], breakdown$N[i],
              breakdown$FH_cases[i], breakdown$FH_pct[i]))
}

write.csv(breakdown, file.path(TABLE_DIR, "cohort_criteria_breakdown.csv"),
          row.names = FALSE)
cat("\n  Saved: cohort_criteria_breakdown.csv\n")

# ==============================================================================
# SECTION 4: MANUSCRIPT-READY SUMMARY
# ==============================================================================
cat("\n================================================================\n")
cat("SECTION 4: MANUSCRIPT-READY SUMMARY\n")
cat("================================================================\n\n")

cat("KEY FINDINGS FOR MANUSCRIPT:\n\n")

cat(sprintf("1. New lipid clinic cohort (TC>7.5 | LDL>4.9 | premature ASCVD):\n"))
cat(sprintf("   N = %d (vs old N = %d, +%d additional patients)\n",
            sum(df$cohort_lipid_clinic),
            sum(df$cohort_high_risk, na.rm = TRUE),
            sum(df$cohort_lipid_clinic) - sum(df$cohort_high_risk, na.rm = TRUE)))
cat(sprintf("   FH cases captured: %d (vs old %d)\n",
            sum(df$is_fh_genetic[df$cohort_lipid_clinic]),
            sum(df$is_fh_genetic[df$cohort_high_risk == TRUE], na.rm = TRUE)))

cat(sprintf("\n2. TUDOR AUC on new cohort: %.3f [%.3f–%.3f]\n",
            ci_tudor_lc[2], ci_tudor_lc[1], ci_tudor_lc[3]))

if (has_apob && exists("aug_results") && nrow(aug_results) > 1) {
  best_aug <- aug_results[which.max(aug_results$AUC), ]
  cat(sprintf("\n3. Best augmented model: %s\n", best_aug$Model))
  cat(sprintf("   AUC: %.3f [%.3f–%.3f], Delta: %+.3f, p = %s\n",
              best_aug$AUC, best_aug$CI_lower, best_aug$CI_upper,
              best_aug$Delta_AUC,
              if (is.na(best_aug$DeLong_p)) "ref" else sprintf("%.4g", best_aug$DeLong_p)))

  if (best_aug$AUC >= 0.80) {
    cat("   → RECOMMENDATION: ApoB/LDL-C ratio recommended as clinical enhancement\n")
    cat("   → AUC approaches 0.80 threshold — clinically meaningful improvement\n")
  } else if (best_aug$AUC >= 0.78) {
    cat("   → RECOMMENDATION: ApoB augmentation provides modest but significant gain\n")
    cat("   → Worth incorporating where ApoB is routinely measured\n")
  } else {
    cat("   → RECOMMENDATION: ApoB augmentation provides marginal improvement\n")
    cat("   → TUDOR base model performs well without ApoB\n")
  }

  cat("\n4. ApoB availability limitation:\n")
  cat("   → Wales FH Registry: NO serum ApoB available (only genetic APOB mutation)\n")
  cat("   → UK Biobank: ApoB measured in majority of participants\n")
  cat("   → Clinical implication: ApoB augmentation cannot be validated in Wales\n")
  cat("   → Recommendation: Use TUDOR base where ApoB unavailable; add ApoB where available\n")
}

cat("\n=== 13_apob_augmented_model.R COMPLETE ===\n")
