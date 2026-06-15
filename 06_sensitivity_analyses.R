# ==============================================================================
# TUDOR PIPELINE: STEP 06 — SENSITIVITY ANALYSES
# ==============================================================================
# PURPOSE: 10 sensitivity analyses to test robustness of primary results.
#          Required for Lancet/BMJ supplementary materials.
#
# REQUIRES: tudor_analysis_ready.rds from 01_data_merge.R
#
# ANALYSES:
#   S1. LDL threshold sensitivity (4.0 - 6.0 mmol/L)
#   S2. Friedewald vs Direct LDL
#   S3. Statin correction factor sensitivity (+/-10%)
#   S4. Ethnicity stratification
#   S5. Outlier exclusion (Winsorisation)
#   S6. MEDPED comparison
#   S7. Prevalence-adjusted PPV/NPV
#   S8. Sex interaction analysis
#   S9. Age-sex interaction
#   S10. Statin-free subgroup (strongest validation)
#
# OUTPUTS: 06_sensitivity_results.rds
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

cat("=== TUDOR PIPELINE: 06_sensitivity_analyses.R ===\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================
df <- readRDS(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))
if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
cat("Loaded:", nrow(df), "participants\n\n")

# TUDOR Wales weights (fixed, for recalculation)
TUDOR_WEIGHTS <- list(
  intercept = 0.755722,
  beta_LDL  = 0.057911,
  beta_Trig = 0.492412,
  beta_HDL  = -1.128045,
  beta_Age  = -0.033393,
  beta_Sex  = -0.088550
)

# Real-world reduction factors
REDUCTION_FACTORS <- c(
  Atorvastatin = 0.38, Simvastatin = 0.35, Rosuvastatin = 0.34,
  Pravastatin = 0.25, Fluvastatin = 0.22
)

sensitivity_results <- list()

# Helper for AUC with bootstrap CI
auc_with_ci <- function(outcome, predictor, n_boot = 2000) {
  valid <- !is.na(outcome) & !is.na(predictor)
  if (sum(outcome[valid] == 1) < 5) return(list(auc = NA, ci = c(NA, NA), n = sum(valid)))
  r <- roc(outcome[valid], predictor[valid], quiet = TRUE)
  ci <- ci.auc(r, conf.level = 0.95, method = "bootstrap",
               boot.n = n_boot, boot.stratified = TRUE)
  list(auc = as.numeric(ci[2]), ci = as.numeric(ci[c(1, 3)]),
       n = sum(valid), n_cases = sum(outcome[valid] == 1), roc = r)
}

# ==============================================================================
# S1. LDL THRESHOLD SENSITIVITY
# ==============================================================================
cat("==========================================================\n")
cat(" S1. LDL THRESHOLD SENSITIVITY\n")
cat("==========================================================\n")

thresholds <- seq(4.0, 6.0, by = 0.5)
s1_results <- list()

cat(sprintf("%-12s | %7s | %6s | %9s | %s\n",
            "Threshold", "N", "Cases", "AUC", "95% CI"))
cat(paste(rep("-", 60), collapse = ""), "\n")

for (thresh in thresholds) {
  sub <- df[df$LDL_RW > thresh, ]
  res <- auc_with_ci(sub$is_fh_genetic, sub$tudor_prob)

  label <- sprintf(">%.1f", thresh)
  if (!is.na(res$auc)) {
    cat(sprintf("%-12s | %7d | %6d | %9.3f | [%.3f - %.3f]\n",
                label, res$n, res$n_cases, res$auc, res$ci[1], res$ci[2]))
  } else {
    cat(sprintf("%-12s | %7d | %6d | %9s | %s\n",
                label, res$n, res$n_cases, "N/A", "Insufficient cases"))
  }

  s1_results[[label]] <- list(threshold = thresh, n = res$n,
                               n_cases = res$n_cases, auc = res$auc, ci = res$ci)
}

cat(paste(rep("-", 60), collapse = ""), "\n\n")
sensitivity_results$s1_threshold <- s1_results

# ==============================================================================
# S2. FRIEDEWALD vs DIRECT LDL
# ==============================================================================
cat("==========================================================\n")
cat(" S2. FRIEDEWALD vs DIRECT LDL\n")
cat("==========================================================\n")

if ("LDL_Friedewald" %in% names(df)) {
  hr <- df[df$cohort_high_risk == TRUE, ]

  # Recalculate TUDOR score using Friedewald LDL
  hr$LDL_Friedewald_RW <- hr$LDL_Friedewald
  # Apply statin correction to Friedewald too
  correction <- ifelse(hr$statin_name != "None",
                       REDUCTION_FACTORS[hr$statin_name], 0)
  correction[is.na(correction)] <- 0
  hr$LDL_Friedewald_RW <- hr$LDL_Friedewald / (1 - correction)

  # Friedewald Trig_Filter
  # BUGFIX: Column names are TRG.1 and HDL.1 in the UKB dataset (not Trig/HDL)
  Trig_Filter_Fried <- hr$LDL_Friedewald_RW / (hr$TRG.1 + 0.1)

  tudor_fried <- TUDOR_WEIGHTS$intercept +
    TUDOR_WEIGHTS$beta_LDL * hr$LDL_Friedewald_RW +
    TUDOR_WEIGHTS$beta_Trig * Trig_Filter_Fried +
    TUDOR_WEIGHTS$beta_HDL * hr$HDL.1 +
    TUDOR_WEIGHTS$beta_Age * hr$Age_at_LDL1 +
    TUDOR_WEIGHTS$beta_Sex * hr$Gender_num

  tudor_prob_fried <- 1 / (1 + exp(-tudor_fried))

  res_direct <- auc_with_ci(hr$is_fh_genetic, hr$tudor_prob)
  res_fried  <- auc_with_ci(hr$is_fh_genetic, tudor_prob_fried)

  cat(sprintf("Direct LDL TUDOR AUC:    %.3f [%.3f - %.3f]\n",
              res_direct$auc, res_direct$ci[1], res_direct$ci[2]))
  cat(sprintf("Friedewald LDL TUDOR AUC: %.3f [%.3f - %.3f]\n",
              res_fried$auc, res_fried$ci[1], res_fried$ci[2]))

  if (!is.na(res_direct$auc) && !is.na(res_fried$auc)) {
    test <- tryCatch(roc.test(res_direct$roc, res_fried$roc, method = "delong"),
                     error = function(e) list(p.value = NA))
    cat(sprintf("DeLong p: %.2e\n", test$p.value))
  }

  sensitivity_results$s2_friedewald <- list(
    direct_auc = res_direct$auc, direct_ci = res_direct$ci,
    friedewald_auc = res_fried$auc, friedewald_ci = res_fried$ci
  )
} else {
  cat("LDL_Friedewald not available. Skipping.\n")
}
cat("\n")

# ==============================================================================
# S3. STATIN CORRECTION FACTOR SENSITIVITY (+/-10%)
# ==============================================================================
cat("==========================================================\n")
cat(" S3. STATIN CORRECTION FACTOR SENSITIVITY\n")
cat("==========================================================\n")

adjustments <- c(0.90, 0.95, 1.00, 1.05, 1.10)  # 90% to 110% of base factors
s3_results <- list()

cat(sprintf("%-12s | %9s | %s\n", "Adjustment", "AUC", "95% CI"))
cat(paste(rep("-", 45), collapse = ""), "\n")

for (adj in adjustments) {
  adj_factors <- REDUCTION_FACTORS * adj
  adj_factors <- pmin(adj_factors, 0.60)  # Cap at 60% reduction

  # Recalculate LDL_RW with adjusted factors
  correction <- ifelse(df$statin_name != "None",
                       adj_factors[df$statin_name], 0)
  correction[is.na(correction)] <- 0
  # BUGFIX: Column is LDL_treated (not LDL_direct) in the dataset
  ldl_rw_adj <- df$LDL_treated / (1 - correction)

  # Recalculate TUDOR
  # BUGFIX: Column names are TRG.1 and HDL.1 in UKB dataset (not Trig/HDL)
  trig_filter_adj <- ldl_rw_adj / (df$TRG.1 + 0.1)
  tudor_adj <- TUDOR_WEIGHTS$intercept +
    TUDOR_WEIGHTS$beta_LDL * ldl_rw_adj +
    TUDOR_WEIGHTS$beta_Trig * trig_filter_adj +
    TUDOR_WEIGHTS$beta_HDL * df$HDL.1 +
    TUDOR_WEIGHTS$beta_Age * df$Age_at_LDL1 +
    TUDOR_WEIGHTS$beta_Sex * df$Gender_num
  tudor_prob_adj <- 1 / (1 + exp(-tudor_adj))

  # High-risk with adjusted LDL
  hr_adj <- ldl_rw_adj > 4.9
  sub <- df[hr_adj, ]
  sub$tudor_prob_adj <- tudor_prob_adj[hr_adj]

  res <- auc_with_ci(sub$is_fh_genetic, sub$tudor_prob_adj)

  label <- sprintf("%.0f%%", adj * 100)
  if (!is.na(res$auc)) {
    cat(sprintf("%-12s | %9.3f | [%.3f - %.3f]  (N=%d)\n",
                label, res$auc, res$ci[1], res$ci[2], res$n))
  }

  s3_results[[label]] <- list(adjustment = adj, auc = res$auc,
                               ci = res$ci, n = res$n)
}

cat(paste(rep("-", 45), collapse = ""), "\n\n")
sensitivity_results$s3_statin_correction <- s3_results

# ==============================================================================
# S4. ETHNICITY STRATIFICATION
# ==============================================================================
cat("==========================================================\n")
cat(" S4. ETHNICITY STRATIFICATION\n")
cat("==========================================================\n")

hr <- df[df$cohort_high_risk == TRUE, ]

if ("Ethnicity" %in% names(hr)) {
  eth_groups <- unique(hr$Ethnicity[!is.na(hr$Ethnicity)])

  cat(sprintf("%-20s | %7s | %6s | %9s | %s\n",
              "Ethnicity", "N", "Cases", "AUC", "95% CI"))
  cat(paste(rep("-", 65), collapse = ""), "\n")

  s4_results <- list()
  for (eth in eth_groups) {
    sub <- hr[hr$Ethnicity == eth & !is.na(hr$Ethnicity), ]
    res <- auc_with_ci(sub$is_fh_genetic, sub$tudor_prob)

    if (!is.na(res$auc)) {
      cat(sprintf("%-20s | %7d | %6d | %9.3f | [%.3f - %.3f]\n",
                  eth, res$n, res$n_cases, res$auc, res$ci[1], res$ci[2]))
    } else {
      cat(sprintf("%-20s | %7d | %6d | %9s | %s\n",
                  eth, res$n, ifelse(is.null(res$n_cases), 0, res$n_cases),
                  "N/A", "Insufficient cases"))
    }

    s4_results[[eth]] <- list(n = res$n, n_cases = res$n_cases,
                               auc = res$auc, ci = res$ci)
  }

  cat(paste(rep("-", 65), collapse = ""), "\n")
  sensitivity_results$s4_ethnicity <- s4_results
} else {
  cat("Ethnicity variable not available.\n")
}
cat("\n")

# ==============================================================================
# S5. OUTLIER EXCLUSION (Winsorisation at 1st/99th percentile)
# ==============================================================================
cat("==========================================================\n")
cat(" S5. OUTLIER EXCLUSION (Winsorisation)\n")
cat("==========================================================\n")

hr <- df[df$cohort_high_risk == TRUE, ]

# Winsorise LDL at 1st and 99th percentile
ldl_q01 <- quantile(hr$LDL_RW, 0.01, na.rm = TRUE)
ldl_q99 <- quantile(hr$LDL_RW, 0.99, na.rm = TRUE)

hr$LDL_RW_wins <- pmax(pmin(hr$LDL_RW, ldl_q99), ldl_q01)

# Recalculate TUDOR with winsorised LDL
# BUGFIX: Column names are TRG.1 and HDL.1 in UKB dataset (not Trig/HDL)
trig_filter_wins <- hr$LDL_RW_wins / (hr$TRG.1 + 0.1)
tudor_wins <- TUDOR_WEIGHTS$intercept +
  TUDOR_WEIGHTS$beta_LDL * hr$LDL_RW_wins +
  TUDOR_WEIGHTS$beta_Trig * trig_filter_wins +
  TUDOR_WEIGHTS$beta_HDL * hr$HDL.1 +
  TUDOR_WEIGHTS$beta_Age * hr$Age_at_LDL1 +
  TUDOR_WEIGHTS$beta_Sex * hr$Gender_num
tudor_prob_wins <- 1 / (1 + exp(-tudor_wins))

res_orig <- auc_with_ci(hr$is_fh_genetic, hr$tudor_prob)
res_wins <- auc_with_ci(hr$is_fh_genetic, tudor_prob_wins)

cat(sprintf("Original AUC:     %.3f [%.3f - %.3f]\n",
            res_orig$auc, res_orig$ci[1], res_orig$ci[2]))
cat(sprintf("Winsorised AUC:   %.3f [%.3f - %.3f]\n",
            res_wins$auc, res_wins$ci[1], res_wins$ci[2]))
cat(sprintf("LDL Winsorisation range: [%.2f, %.2f] mmol/L\n",
            ldl_q01, ldl_q99))
cat(sprintf("Participants affected: %d (%.1f%%)\n\n",
            sum(hr$LDL_RW < ldl_q01 | hr$LDL_RW > ldl_q99),
            mean(hr$LDL_RW < ldl_q01 | hr$LDL_RW > ldl_q99) * 100))

sensitivity_results$s5_winsorisation <- list(
  original = list(auc = res_orig$auc, ci = res_orig$ci),
  winsorised = list(auc = res_wins$auc, ci = res_wins$ci),
  ldl_range = c(ldl_q01, ldl_q99)
)

# ==============================================================================
# S6. MEDPED COMPARISON
# ==============================================================================
cat("==========================================================\n")
cat(" S6. MEDPED COMPARISON\n")
cat("==========================================================\n")

hr <- df[df$cohort_high_risk == TRUE, ]

# MEDPED criteria (age-sex specific LDL thresholds, mmol/L)
# Age-dependent thresholds for general population (no known FH relative)
medped_threshold <- function(age) {
  ifelse(age < 20, 5.7,
  ifelse(age < 30, 6.2,
  ifelse(age < 40, 7.0,
                   7.5)))
}

hr$medped_threshold <- medped_threshold(hr$Age_at_LDL1)
hr$medped_positive <- hr$LDL_RW >= hr$medped_threshold

medped_score <- as.numeric(hr$medped_positive)

res_tudor <- auc_with_ci(hr$is_fh_genetic, hr$tudor_prob)
res_edlcn <- auc_with_ci(hr$is_fh_genetic, hr$edlcn_score)

# MEDPED is binary so we report sensitivity/specificity rather than AUC
medped_sens <- sum(hr$medped_positive & hr$is_fh_genetic) /
               sum(hr$is_fh_genetic)
medped_spec <- sum(!hr$medped_positive & !hr$is_fh_genetic) /
               sum(!hr$is_fh_genetic)
medped_ppv <- sum(hr$medped_positive & hr$is_fh_genetic) /
              max(sum(hr$medped_positive), 1)
medped_npv <- sum(!hr$medped_positive & !hr$is_fh_genetic) /
              max(sum(!hr$medped_positive), 1)

cat(sprintf("TUDOR AUC:   %.3f [%.3f - %.3f]\n",
            res_tudor$auc, res_tudor$ci[1], res_tudor$ci[2]))
cat(sprintf("eDLCN AUC:   %.3f [%.3f - %.3f]\n",
            res_edlcn$auc, res_edlcn$ci[1], res_edlcn$ci[2]))
cat(sprintf("\nMEDPED (binary):\n"))
cat(sprintf("  Sensitivity: %.1f%%\n", medped_sens * 100))
cat(sprintf("  Specificity: %.1f%%\n", medped_spec * 100))
cat(sprintf("  PPV:         %.1f%%\n", medped_ppv * 100))
cat(sprintf("  NPV:         %.1f%%\n", medped_npv * 100))
cat(sprintf("  Positive:    %d / %d\n\n", sum(hr$medped_positive), nrow(hr)))

sensitivity_results$s6_medped <- list(
  sensitivity = medped_sens, specificity = medped_spec,
  ppv = medped_ppv, npv = medped_npv,
  n_positive = sum(hr$medped_positive)
)

# ==============================================================================
# S7. PREVALENCE-ADJUSTED PPV/NPV
# ==============================================================================
cat("==========================================================\n")
cat(" S7. PREVALENCE-ADJUSTED PPV/NPV\n")
cat("==========================================================\n")

# Different assumed prevalences for FH in the general population
assumed_prevalences <- c(1/500, 1/250, 1/200, 1/100)

# Get TUDOR sensitivity and specificity at optimal threshold
roc_tudor <- roc(hr$is_fh_genetic, hr$tudor_prob, quiet = TRUE)
opt_coords <- coords(roc_tudor, "best", ret = c("threshold", "sensitivity", "specificity"))
sens <- opt_coords$sensitivity
spec <- opt_coords$specificity

cat(sprintf("TUDOR optimal threshold: %.4f\n", opt_coords$threshold))
cat(sprintf("Sensitivity: %.1f%%\n", sens * 100))
cat(sprintf("Specificity: %.1f%%\n\n", spec * 100))

cat(sprintf("%-15s | %8s | %8s\n", "Prevalence", "PPV", "NPV"))
cat(paste(rep("-", 40), collapse = ""), "\n")

s7_results <- list()
for (prev in assumed_prevalences) {
  ppv <- (sens * prev) / (sens * prev + (1 - spec) * (1 - prev))
  npv <- (spec * (1 - prev)) / (spec * (1 - prev) + (1 - sens) * prev)

  label <- sprintf("1/%d", round(1/prev))
  cat(sprintf("%-15s | %7.1f%% | %7.1f%%\n", label, ppv * 100, npv * 100))

  s7_results[[label]] <- list(prevalence = prev, ppv = ppv, npv = npv)
}
cat(paste(rep("-", 40), collapse = ""), "\n\n")

sensitivity_results$s7_prevalence <- s7_results

# ==============================================================================
# S8. SEX INTERACTION ANALYSIS
# ==============================================================================
cat("==========================================================\n")
cat(" S8. SEX INTERACTION ANALYSIS\n")
cat("==========================================================\n")

hr <- df[df$cohort_high_risk == TRUE, ]

# Test if TUDOR performance differs significantly between sexes
men <- hr[hr$Gender_num == 1, ]
women <- hr[hr$Gender_num == 0, ]

res_men <- auc_with_ci(men$is_fh_genetic, men$tudor_prob)
res_women <- auc_with_ci(women$is_fh_genetic, women$tudor_prob)

cat(sprintf("Men:   AUC = %.3f [%.3f - %.3f] (N=%d, Cases=%d)\n",
            res_men$auc, res_men$ci[1], res_men$ci[2],
            res_men$n, res_men$n_cases))
cat(sprintf("Women: AUC = %.3f [%.3f - %.3f] (N=%d, Cases=%d)\n",
            res_women$auc, res_women$ci[1], res_women$ci[2],
            res_women$n, res_women$n_cases))

# Interaction test: logistic regression with sex * tudor_score
fit_interaction <- glm(is_fh_genetic ~ tudor_score * Gender_num,
                       data = hr, family = binomial)
interaction_p <- summary(fit_interaction)$coefficients["tudor_score:Gender_num", "Pr(>|z|)"]

cat(sprintf("Interaction p-value (score x sex): %.3f\n", interaction_p))
cat(sprintf("(p > 0.05 = no significant interaction)\n\n"))

sensitivity_results$s8_sex_interaction <- list(
  men = list(auc = res_men$auc, ci = res_men$ci, n = res_men$n),
  women = list(auc = res_women$auc, ci = res_women$ci, n = res_women$n),
  interaction_p = interaction_p
)

# ==============================================================================
# S9. AGE-SEX INTERACTION
# ==============================================================================
cat("==========================================================\n")
cat(" S9. AGE x SEX INTERACTION\n")
cat("==========================================================\n")

fit_age_sex <- glm(is_fh_genetic ~ tudor_score * Age_at_LDL1 * Gender_num,
                   data = hr, family = binomial)

cat("Three-way interaction model:\n")
coefs <- summary(fit_age_sex)$coefficients
# Print only interaction terms
interaction_terms <- grep(":", rownames(coefs), value = TRUE)
for (term in interaction_terms) {
  cat(sprintf("  %-40s: coef = %+.4f, p = %.3f\n",
              term, coefs[term, "Estimate"], coefs[term, "Pr(>|z|)"]))
}

# LRT for interaction
fit_no_int <- glm(is_fh_genetic ~ tudor_score + Age_at_LDL1 + Gender_num,
                  data = hr, family = binomial)
lr_test <- anova(fit_no_int, fit_age_sex, test = "Chisq")
cat(sprintf("\nLR test (interaction terms): p = %.3f\n\n",
            lr_test$`Pr(>Chi)`[2]))

sensitivity_results$s9_age_sex <- list(
  interaction_coefs = coefs[interaction_terms, , drop = FALSE],
  lr_p = lr_test$`Pr(>Chi)`[2]
)

# ==============================================================================
# S10. STATIN-FREE SUBGROUP (Purest Validation)
# ==============================================================================
cat("==========================================================\n")
cat(" S10. STATIN-FREE SUBGROUP\n")
cat("==========================================================\n")

hr <- df[df$cohort_high_risk == TRUE, ]
statin_free <- hr[hr$statin_name == "None", ]
on_statin <- hr[hr$statin_name != "None", ]

res_free <- auc_with_ci(statin_free$is_fh_genetic, statin_free$tudor_prob)
res_statin <- auc_with_ci(on_statin$is_fh_genetic, on_statin$tudor_prob)

cat(sprintf("Statin-Free: AUC = %.3f [%.3f - %.3f] (N=%d, Cases=%d)\n",
            res_free$auc, res_free$ci[1], res_free$ci[2],
            res_free$n, res_free$n_cases))
cat(sprintf("On Statin:   AUC = %.3f [%.3f - %.3f] (N=%d, Cases=%d)\n",
            res_statin$auc, res_statin$ci[1], res_statin$ci[2],
            res_statin$n, res_statin$n_cases))

# Also compare eDLCN in statin-free
res_edlcn_free <- auc_with_ci(statin_free$is_fh_genetic, statin_free$edlcn_score)

if (!is.na(res_free$auc) && !is.na(res_edlcn_free$auc)) {
  test_free <- tryCatch(
    roc.test(res_free$roc, res_edlcn_free$roc, method = "delong"),
    error = function(e) list(p.value = NA)
  )
  cat(sprintf("\nStatin-free subgroup: TUDOR vs eDLCN DeLong p = %.2e\n",
              test_free$p.value))
}
cat("\n")

sensitivity_results$s10_statin_free <- list(
  statin_free = list(auc = res_free$auc, ci = res_free$ci,
                     n = res_free$n, n_cases = res_free$n_cases),
  on_statin = list(auc = res_statin$auc, ci = res_statin$ci,
                   n = res_statin$n, n_cases = res_statin$n_cases)
)

# ==============================================================================
# SAVE ALL RESULTS
# ==============================================================================
saveRDS(sensitivity_results, file.path(OUTPUT_DIR, "06_sensitivity_results.rds"))

cat("=== 06_sensitivity_analyses.R COMPLETE ===\n")
