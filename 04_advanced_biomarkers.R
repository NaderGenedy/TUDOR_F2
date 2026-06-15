# ==============================================================================
# TUDOR PIPELINE: STEP 04 — ADVANCED BIOMARKER ANALYSIS
# ==============================================================================
# PURPOSE: Evaluate additional biomarkers (Lp(a), ApoB, ApoB/LDL ratio)
#          as potential adjuncts to the TUDOR model.
#
# REQUIRES: - tudor_analysis_ready.rds from 01_data_merge.R
#           - ukb_lpa.csv (Lp(a) from 00_extract_ukbrap.sh)
#
# OUTPUTS:  - Biomarker head-to-head comparison
#           - TUDOR + biomarker augmentation analysis
#           - Lp(a) correction impact
#           - 04_biomarker_results.rds
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

cat("=== TUDOR PIPELINE: 04_advanced_biomarkers.R ===\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
df <- readRDS(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))
if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
cat("Loaded:", nrow(df), "participants\n")

# --- Load Lp(a) if available ---
# BUGFIX: Script 01 already merges Lp(a) into the RDS. If we re-merge from
# CSV here, it creates Lpa_nmol.x / Lpa_nmol.y columns and 0 usable values.
# Solution: Check if Lpa_nmol already exists with data; only load from CSV if not.

if ("Lpa_nmol" %in% names(df) && sum(!is.na(df$Lpa_nmol)) > 0) {
  has_lpa <- TRUE
  cat("Lp(a) already present in RDS:", sum(!is.na(df$Lpa_nmol)), "non-missing values\n")
} else {
  lpa_file <- file.path(DATA_DIR, "ukb_lpa.csv")
  has_lpa <- file.exists(lpa_file)

  if (has_lpa) {
    lpa <- fread(lpa_file)
    if ("participant.eid" %in% names(lpa) && !"eid" %in% names(lpa)) {
      setnames(lpa, "participant.eid", "eid")
    }
    lpa_col <- grep("p30790", names(lpa), value = TRUE)[1]

    if (!is.null(lpa_col) && length(lpa_col) > 0) {
      lpa$Lpa_nmol <- as.numeric(lpa[[lpa_col]])
      df <- merge(df, lpa[, c("eid", "Lpa_nmol")], by = "eid", all.x = TRUE)
      cat("Lp(a) merged from CSV:", sum(!is.na(df$Lpa_nmol)), "non-missing values\n")
    } else {
      has_lpa <- FALSE
      cat("WARNING: Lp(a) column not found in file.\n")
    }
  } else {
    cat("WARNING: ukb_lpa.csv not found. Lp(a) analyses will be skipped.\n")
    cat("  Run 00_extract_ukbrap.sh on UKB-RAP to extract it.\n")
  }
}

# High-risk cohort
hr <- df[df$cohort_high_risk == TRUE, ]
cat("High-risk cohort:", nrow(hr), "participants\n")
cat("Genetic FH cases:", sum(hr$is_fh_genetic), "\n\n")

# ==============================================================================
# 2. BIOMARKER HEAD-TO-HEAD COMPARISON
# ==============================================================================
cat("==========================================================\n")
cat(" BIOMARKER HEAD-TO-HEAD: AUC COMPARISON (HIGH-RISK)\n")
cat("==========================================================\n")

# Helper function for AUC with bootstrap CI
get_auc_ci <- function(response, predictor, label, data_subset = NULL) {
  if (!is.null(data_subset)) {
    valid <- !is.na(predictor[data_subset]) & !is.na(response[data_subset])
    resp <- response[data_subset][valid]
    pred <- predictor[data_subset][valid]
  } else {
    valid <- !is.na(predictor) & !is.na(response)
    resp <- response[valid]
    pred <- predictor[valid]
  }

  if (sum(resp == 1) < 5 || sum(resp == 0) < 5) {
    return(list(label = label, auc = NA, ci_low = NA, ci_high = NA,
                n = length(resp), n_cases = sum(resp)))
  }

  r <- roc(resp, pred, quiet = TRUE)
  ci <- ci.auc(r, conf.level = 0.95, method = "bootstrap",
               boot.n = 2000, boot.stratified = TRUE)

  list(label = label, auc = as.numeric(ci[2]),
       ci_low = as.numeric(ci[1]), ci_high = as.numeric(ci[3]),
       n = length(resp), n_cases = sum(resp), roc_obj = r)
}

# Individual biomarkers
biomarkers <- list()

biomarkers[["TUDOR"]] <- get_auc_ci(hr$is_fh_genetic, hr$tudor_prob, "TUDOR v2")
biomarkers[["eDLCN"]] <- get_auc_ci(hr$is_fh_genetic, hr$edlcn_score, "eDLCN")
biomarkers[["LDL"]]   <- get_auc_ci(hr$is_fh_genetic, hr$LDL_RW, "LDL-C (RW)")

if ("ApoB" %in% names(hr)) {
  biomarkers[["ApoB"]] <- get_auc_ci(hr$is_fh_genetic, hr$ApoB, "ApoB")
}
if ("ApoB_LDL_Ratio" %in% names(hr)) {
  biomarkers[["ApoB_LDL"]] <- get_auc_ci(hr$is_fh_genetic, hr$ApoB_LDL_Ratio,
                                           "ApoB/LDL Ratio")
}
if ("TC_LDL_Ratio" %in% names(hr)) {
  biomarkers[["TC_LDL"]] <- get_auc_ci(hr$is_fh_genetic, hr$TC_LDL_Ratio,
                                         "TC/LDL Ratio")
}
if ("Trig_Filter_RW" %in% names(hr)) {
  biomarkers[["Trig_Filter"]] <- get_auc_ci(hr$is_fh_genetic, hr$Trig_Filter_RW,
                                              "Trig Filter (RW)")
}

# Print comparison table
cat(sprintf("%-20s | %6s | %9s | %s\n", "Biomarker", "N", "AUC", "95% CI"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (b in biomarkers) {
  if (is.na(b$auc)) {
    cat(sprintf("%-20s | %6d | %9s | %s\n", b$label, b$n, "N/A", "Insufficient data"))
  } else {
    cat(sprintf("%-20s | %6d | %9.3f | [%.3f - %.3f]\n",
                b$label, b$n, b$auc, b$ci_low, b$ci_high))
  }
}
cat(paste(rep("-", 60), collapse = ""), "\n\n")

# ==============================================================================
# 3. DELONG PAIRWISE COMPARISONS vs TUDOR
# ==============================================================================
cat("==========================================================\n")
cat(" DELONG TESTS: EACH BIOMARKER vs TUDOR\n")
cat("==========================================================\n")

if (!is.null(biomarkers[["TUDOR"]]$roc_obj)) {
  roc_tudor <- biomarkers[["TUDOR"]]$roc_obj

  for (name in names(biomarkers)) {
    if (name == "TUDOR") next
    b <- biomarkers[[name]]
    if (is.na(b$auc) || is.null(b$roc_obj)) next

    # DeLong test requires same cases — use complete cases for both
    test <- tryCatch(
      roc.test(roc_tudor, b$roc_obj, method = "delong"),
      error = function(e) NULL
    )

    if (!is.null(test)) {
      delta <- biomarkers[["TUDOR"]]$auc - b$auc
      cat(sprintf("TUDOR vs %-15s: delta = %+.3f, p = %.2e\n",
                  b$label, delta, test$p.value))
    }
  }
}
cat("\n")

# ==============================================================================
# 4. TUDOR + BIOMARKER AUGMENTATION
# ==============================================================================
cat("==========================================================\n")
cat(" TUDOR + BIOMARKER AUGMENTATION\n")
cat("==========================================================\n")
cat("(Logistic regression adding each biomarker to TUDOR score)\n\n")

augmentation_results <- list()

# Base model: TUDOR score alone (logistic regression to get calibrated probs)
base_formula <- is_fh_genetic ~ tudor_score

for (aug_name in c("ApoB", "ApoB_LDL_Ratio", "TC_LDL_Ratio")) {
  if (!(aug_name %in% names(hr))) next

  # Complete cases for this biomarker
  valid <- !is.na(hr[[aug_name]]) & !is.na(hr$tudor_score) & !is.na(hr$is_fh_genetic)
  sub <- hr[valid, ]

  if (nrow(sub) < 100 || sum(sub$is_fh_genetic) < 10) {
    cat(sprintf("%-20s: Insufficient data (N=%d, cases=%d)\n",
                aug_name, nrow(sub), sum(sub$is_fh_genetic)))
    next
  }

  # Fit base on this subset
  fit_base <- glm(is_fh_genetic ~ tudor_score, data = sub, family = binomial)
  pred_base <- predict(fit_base, type = "response")

  # Fit augmented
  aug_formula <- as.formula(paste("is_fh_genetic ~ tudor_score +", aug_name))
  fit_aug <- glm(aug_formula, data = sub, family = binomial)
  pred_aug <- predict(fit_aug, type = "response")

  # ROC comparison
  roc_base <- roc(sub$is_fh_genetic, pred_base, quiet = TRUE)
  roc_aug  <- roc(sub$is_fh_genetic, pred_aug, quiet = TRUE)

  ci_base <- ci.auc(roc_base, conf.level = 0.95, method = "bootstrap",
                     boot.n = 2000, boot.stratified = TRUE)
  ci_aug  <- ci.auc(roc_aug, conf.level = 0.95, method = "bootstrap",
                     boot.n = 2000, boot.stratified = TRUE)

  test <- tryCatch(roc.test(roc_base, roc_aug, method = "delong"),
                   error = function(e) list(p.value = NA))

  # Likelihood ratio test
  lr_test <- anova(fit_base, fit_aug, test = "Chisq")
  lr_p <- lr_test$`Pr(>Chi)`[2]

  cat(sprintf("%-20s: Base AUC = %.3f, Aug AUC = %.3f, delta = %+.3f\n",
              aug_name, auc(roc_base), auc(roc_aug),
              auc(roc_aug) - auc(roc_base)))
  cat(sprintf("                      DeLong p = %.2e, LR test p = %.2e\n",
              test$p.value, lr_p))
  cat(sprintf("                      N = %d, Cases = %d\n\n",
              nrow(sub), sum(sub$is_fh_genetic)))

  augmentation_results[[aug_name]] <- list(
    base_auc = auc(roc_base), aug_auc = auc(roc_aug),
    delta = auc(roc_aug) - auc(roc_base),
    delong_p = test$p.value, lr_p = lr_p,
    n = nrow(sub), n_cases = sum(sub$is_fh_genetic),
    coef_biomarker = coef(fit_aug)[aug_name]
  )
}

# ==============================================================================
# 5. Lp(a) ANALYSIS
# ==============================================================================
if (has_lpa && "Lpa_nmol" %in% names(hr)) {
  cat("==========================================================\n")
  cat(" Lp(a) ANALYSIS\n")
  cat("==========================================================\n")

  lpa_valid <- hr[!is.na(hr$Lpa_nmol), ]
  cat("High-risk participants with Lp(a):", nrow(lpa_valid), "\n")
  cat("FH cases with Lp(a):", sum(lpa_valid$is_fh_genetic), "\n\n")

  if (nrow(lpa_valid) >= 100 && sum(lpa_valid$is_fh_genetic) >= 10) {

    # --- 5a. Lp(a) alone ---
    roc_lpa <- roc(lpa_valid$is_fh_genetic, lpa_valid$Lpa_nmol, quiet = TRUE)
    ci_lpa <- ci.auc(roc_lpa, conf.level = 0.95, method = "bootstrap",
                      boot.n = 2000, boot.stratified = TRUE)
    cat(sprintf("Lp(a) alone AUC: %.3f (95%% CI: %.3f - %.3f)\n",
                ci_lpa[2], ci_lpa[1], ci_lpa[3]))

    # --- 5b. TUDOR + Lp(a) ---
    fit_tudor_lpa <- glm(is_fh_genetic ~ tudor_score + Lpa_nmol,
                         data = lpa_valid, family = binomial)
    pred_tudor_lpa <- predict(fit_tudor_lpa, type = "response")

    roc_tudor_lpa <- roc(lpa_valid$is_fh_genetic, pred_tudor_lpa, quiet = TRUE)
    ci_tudor_lpa <- ci.auc(roc_tudor_lpa, conf.level = 0.95, method = "bootstrap",
                            boot.n = 2000, boot.stratified = TRUE)

    roc_tudor_alone <- roc(lpa_valid$is_fh_genetic, lpa_valid$tudor_prob, quiet = TRUE)

    cat(sprintf("TUDOR alone AUC (Lp(a) subset): %.3f\n", auc(roc_tudor_alone)))
    cat(sprintf("TUDOR + Lp(a) AUC: %.3f (95%% CI: %.3f - %.3f)\n",
                ci_tudor_lpa[2], ci_tudor_lpa[1], ci_tudor_lpa[3]))
    cat(sprintf("Delta: %+.3f\n", auc(roc_tudor_lpa) - auc(roc_tudor_alone)))

    test_lpa <- tryCatch(
      roc.test(roc_tudor_alone, roc_tudor_lpa, method = "delong"),
      error = function(e) list(p.value = NA)
    )
    cat(sprintf("DeLong p: %.2e\n\n", test_lpa$p.value))

    # --- 5c. Lp(a)-corrected LDL ---
    # Holland et al. correction: LDL_corrected = LDL - Lp(a)_mass_mg_dL * 0.30
    # Lp(a) nmol/L to mg/dL: divide by ~2.4 (approximate)
    lpa_valid$Lpa_mg_dL <- lpa_valid$Lpa_nmol / 2.4
    lpa_valid$LDL_Lpa_corrected <- lpa_valid$LDL_RW - (lpa_valid$Lpa_mg_dL * 0.30 / 38.67)
    # Note: 38.67 converts mg/dL to mmol/L for LDL

    roc_ldl_corr <- roc(lpa_valid$is_fh_genetic, lpa_valid$LDL_Lpa_corrected, quiet = TRUE)
    roc_ldl_uncorr <- roc(lpa_valid$is_fh_genetic, lpa_valid$LDL_RW, quiet = TRUE)

    cat("--- Lp(a)-Corrected LDL Analysis ---\n")
    cat(sprintf("LDL (uncorrected) AUC: %.3f\n", auc(roc_ldl_uncorr)))
    cat(sprintf("LDL (Lp(a)-corrected) AUC: %.3f\n", auc(roc_ldl_corr)))
    cat(sprintf("Delta: %+.3f\n\n", auc(roc_ldl_corr) - auc(roc_ldl_uncorr)))

    # --- 5d. Lp(a) distribution in FH vs non-FH ---
    cat("--- Lp(a) Distribution ---\n")
    fh_lpa <- lpa_valid$Lpa_nmol[lpa_valid$is_fh_genetic == TRUE]
    nonfh_lpa <- lpa_valid$Lpa_nmol[lpa_valid$is_fh_genetic == FALSE]

    cat(sprintf("FH cases:     median %.1f nmol/L (IQR: %.1f - %.1f), N=%d\n",
                median(fh_lpa, na.rm = TRUE),
                quantile(fh_lpa, 0.25, na.rm = TRUE),
                quantile(fh_lpa, 0.75, na.rm = TRUE),
                length(fh_lpa)))
    cat(sprintf("Non-FH:       median %.1f nmol/L (IQR: %.1f - %.1f), N=%d\n",
                median(nonfh_lpa, na.rm = TRUE),
                quantile(nonfh_lpa, 0.25, na.rm = TRUE),
                quantile(nonfh_lpa, 0.75, na.rm = TRUE),
                length(nonfh_lpa)))

    wt <- wilcox.test(fh_lpa, nonfh_lpa)
    cat(sprintf("Wilcoxon p-value: %.2e\n\n", wt$p.value))

  } else {
    cat("Insufficient Lp(a) data for analysis.\n\n")
  }
}

# ==============================================================================
# 6. FULL MULTIMARKER MODEL (TUDOR + ApoB + Lp(a))
# ==============================================================================
cat("==========================================================\n")
cat(" FULL MULTIMARKER MODEL\n")
cat("==========================================================\n")

# Build the best possible model using all available biomarkers
available_markers <- c("tudor_score")
if ("ApoB" %in% names(hr)) available_markers <- c(available_markers, "ApoB")
if ("ApoB_LDL_Ratio" %in% names(hr)) available_markers <- c(available_markers, "ApoB_LDL_Ratio")
if (has_lpa && "Lpa_nmol" %in% names(hr)) available_markers <- c(available_markers, "Lpa_nmol")

cat("Available markers:", paste(available_markers, collapse = ", "), "\n")

# Complete cases across all markers
valid_all <- complete.cases(hr[, available_markers, with = FALSE])
hr_complete <- hr[valid_all, ]

cat("Complete cases:", nrow(hr_complete), "\n")
cat("FH cases:", sum(hr_complete$is_fh_genetic), "\n")

if (nrow(hr_complete) >= 100 && sum(hr_complete$is_fh_genetic) >= 10) {
  # Base: TUDOR alone
  fit_base <- glm(is_fh_genetic ~ tudor_score, data = hr_complete, family = binomial)
  roc_base <- roc(hr_complete$is_fh_genetic, predict(fit_base, type = "response"),
                  quiet = TRUE)

  # Full: All markers
  full_formula <- as.formula(paste("is_fh_genetic ~",
                                    paste(available_markers, collapse = " + ")))
  fit_full <- glm(full_formula, data = hr_complete, family = binomial)
  roc_full <- roc(hr_complete$is_fh_genetic, predict(fit_full, type = "response"),
                  quiet = TRUE)

  ci_full <- ci.auc(roc_full, conf.level = 0.95, method = "bootstrap",
                     boot.n = 2000, boot.stratified = TRUE)

  cat(sprintf("\nTUDOR alone AUC: %.3f\n", auc(roc_base)))
  cat(sprintf("Full model AUC:  %.3f (95%% CI: %.3f - %.3f)\n",
              ci_full[2], ci_full[1], ci_full[3]))
  cat(sprintf("Delta: %+.3f\n", auc(roc_full) - auc(roc_base)))

  test_full <- tryCatch(
    roc.test(roc_base, roc_full, method = "delong"),
    error = function(e) list(p.value = NA)
  )
  cat(sprintf("DeLong p: %.2e\n\n", test_full$p.value))

  cat("--- Full Model Coefficients ---\n")
  print(summary(fit_full)$coefficients)
  cat("\n")
} else {
  cat("Insufficient complete cases for multimarker model.\n\n")
}

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================
biomarker_results <- list(
  individual_aucs = lapply(biomarkers, function(b) {
    list(label = b$label, auc = b$auc, ci_low = b$ci_low, ci_high = b$ci_high,
         n = b$n, n_cases = b$n_cases)
  }),
  augmentation = augmentation_results,
  has_lpa = has_lpa
)

saveRDS(biomarker_results, file.path(OUTPUT_DIR, "04_biomarker_results.rds"))

cat("=== 04_advanced_biomarkers.R COMPLETE ===\n")
