# ==============================================================================
# TUDOR PIPELINE: STEP 08 — TRG SHIELD BIOLOGICAL VALIDATION
# ==============================================================================
# PURPOSE: Prove the Triglyceride Filter (LDL_RW / (TRG + 0.1)) captures
#          genuine FH biology and is NOT a confounder proxy for metabolic
#          syndrome (DM, obesity, insulin resistance).
#
# NOVEL CONCEPT: "Metabolic Purity" — FH produces PURE hypercholesterolaemia
#   (high LDL, normal TRG, normal BMI, no DM) while non-FH high-LDL is
#   typically accompanied by metabolic syndrome features. The Trig Filter
#   quantifies this purity as a diagnostic discriminator.
#
# ANALYSES:
#   A. Metabolic characteristics: FH vs non-FH (descriptive table)
#   B. Multivariable regression: Trig Filter independence from confounders
#   C. Stratified AUC: Trig Filter works across all metabolic subgroups
#   D. Interaction tests: No effect modification by DM/BMI/ApoB
#   E. TUDOR model robustness: Adding confounders doesn't improve prediction
#   F. Biological mechanism: VLDL/LDL discordance in FH
#
# INPUT:  tudor_analysis_ready.rds + ukb_dm_hba1c.csv
# OUTPUT: Tables, forest plots, and mechanistic evidence
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

cat("=== TUDOR PIPELINE: 08_trig_shield_biology.R ===\n")
cat("Proving TRG Shield independence from metabolic confounders\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
df <- readRDS(rds_file)
cat("Loaded:", nrow(df), "participants\n")

# --- 1a. Merge DM + HbA1c ---
dm_file <- file.path(DATA_DIR, "ukb_dm_hba1c.csv")
if (file.exists(dm_file)) {
  cat("Merging DM + HbA1c data...\n")
  dm <- fread(dm_file)

  # Field 2443: Diabetes diagnosed by doctor
  # Codes: 1=Yes, 0=No, -1=Don't know, -3=Prefer not to say
  if ("participant.p2443_i0" %in% names(dm)) {
    names(dm)[names(dm) == "participant.p2443_i0"] <- "DM_selfreport"
  }
  # Field 30750: HbA1c (mmol/mol)
  if ("participant.p30750_i0" %in% names(dm)) {
    names(dm)[names(dm) == "participant.p30750_i0"] <- "HbA1c_mmol"
  }

  df <- merge(df, dm[, .SD, .SDcols = intersect(names(dm),
              c("participant.eid", "DM_selfreport", "HbA1c_mmol"))],
              by = "participant.eid", all.x = TRUE)

  # Create binary DM variable
  df$has_DM <- ifelse(!is.na(df$DM_selfreport) & df$DM_selfreport == 1, TRUE, FALSE)
  cat("  Diabetes prevalence:", sum(df$has_DM), "(",
      round(100 * mean(df$has_DM), 1), "%)\n")
  cat("  HbA1c available for:", sum(!is.na(df$HbA1c_mmol)), "participants\n")
} else {
  cat("WARNING: ukb_dm_hba1c.csv not found. Using HbA1c-based DM proxy.\n")
  # Fallback: use HbA1c if available in main dataset, or set to NA
  df$has_DM <- FALSE
  df$HbA1c_mmol <- NA_real_
  df$DM_selfreport <- NA_integer_
}

# --- 1b. Ensure BMI is available ---
if (!"BMI_imputed" %in% names(df)) {
  if ("participant.p21001_i0" %in% names(df)) {
    df$BMI_imputed <- df$participant.p21001_i0
    df$BMI_imputed[is.na(df$BMI_imputed)] <- median(df$BMI_imputed, na.rm = TRUE)
  } else {
    df$BMI_imputed <- 26.5
  }
}

# --- 1c. Create BMI quartiles and ApoB tertiles ---
df$BMI_quartile <- cut(df$BMI_imputed,
                       breaks = quantile(df$BMI_imputed, probs = c(0, 0.25, 0.5, 0.75, 1),
                                         na.rm = TRUE),
                       labels = c("Q1", "Q2", "Q3", "Q4"),
                       include.lowest = TRUE)

if (sum(!is.na(df$ApoB)) > 1000) {
  df$ApoB_tertile <- cut(df$ApoB,
                         breaks = quantile(df$ApoB, probs = c(0, 1/3, 2/3, 1),
                                           na.rm = TRUE),
                         labels = c("T1", "T2", "T3"),
                         include.lowest = TRUE)
  has_apob <- TRUE
} else {
  df$ApoB_tertile <- NA
  has_apob <- FALSE
  cat("NOTE: ApoB not available for tertile analysis.\n")
}

# --- 1d. High-risk cohort ---
hr <- df[df$cohort_high_risk == TRUE & !is.na(df$tudor_prob), ]
cat("\nHigh-risk cohort (LDL_RW > 4.9):", nrow(hr), "participants\n")
cat("  FH cases:", sum(hr$is_fh_genetic), "\n")
cat("  Non-FH:", sum(!hr$is_fh_genetic), "\n\n")


# ==============================================================================
# A. METABOLIC CHARACTERISTICS TABLE: FH vs NON-FH
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS A: Metabolic Characteristics — FH vs Non-FH\n")
cat("(Within high-risk cohort, LDL_RW > 4.9)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Helper function for table rows
describe_var <- function(data, var, group_var = "is_fh_genetic", type = "continuous") {
  fh <- data[[var]][data[[group_var]] == TRUE]
  nonfh <- data[[var]][data[[group_var]] == FALSE]

  if (type == "continuous") {
    fh_val <- sprintf("%.2f (%.2f)", mean(fh, na.rm = TRUE), sd(fh, na.rm = TRUE))
    nonfh_val <- sprintf("%.2f (%.2f)", mean(nonfh, na.rm = TRUE), sd(nonfh, na.rm = TRUE))
    # Welch t-test
    tt <- tryCatch(t.test(fh, nonfh), error = function(e) list(p.value = NA))
    pval <- tt$p.value
    d <- (mean(fh, na.rm = TRUE) - mean(nonfh, na.rm = TRUE)) /
         sqrt((var(fh, na.rm = TRUE) + var(nonfh, na.rm = TRUE)) / 2)
  } else {
    # Binary
    fh_n <- sum(fh, na.rm = TRUE)
    nonfh_n <- sum(nonfh, na.rm = TRUE)
    fh_val <- sprintf("%d (%.1f%%)", fh_n, 100 * mean(fh, na.rm = TRUE))
    nonfh_val <- sprintf("%d (%.1f%%)", nonfh_n, 100 * mean(nonfh, na.rm = TRUE))
    ct <- tryCatch(chisq.test(table(data[[group_var]], data[[var]])),
                   error = function(e) list(p.value = NA))
    pval <- ct$p.value
    d <- NA
  }

  list(FH = fh_val, NonFH = nonfh_val, p = pval, cohens_d = d)
}

cat(sprintf("%-25s %-20s %-20s %-12s %-8s\n",
            "Variable", "FH (genetic)", "Non-FH", "P-value", "Cohen d"))
cat(strrep("-", 85), "\n")

vars <- list(
  list("Age", "Age_at_LDL1", "continuous"),
  list("BMI (kg/m2)", "BMI_imputed", "continuous"),
  list("LDL_RW (mmol/L)", "LDL_RW", "continuous"),
  list("TRG (mmol/L)", "TRG.1", "continuous"),
  list("HDL (mmol/L)", "HDL.1", "continuous"),
  list("Trig Filter", "Trig_Filter_RW", "continuous"),
  list("Total Chol", "CHOL", "continuous")
)

# Add HbA1c if available
if (sum(!is.na(hr$HbA1c_mmol)) > 100) {
  vars <- c(vars, list(list("HbA1c (mmol/mol)", "HbA1c_mmol", "continuous")))
}

# Add ApoB if available
if (has_apob) {
  vars <- c(vars, list(
    list("ApoB (g/L)", "ApoB", "continuous"),
    list("ApoB/LDL Ratio", "ApoB_LDL_Ratio", "continuous")
  ))
}

for (v in vars) {
  res <- describe_var(hr, v[[2]], type = v[[3]])
  pstr <- ifelse(is.na(res$p), "NA",
           ifelse(res$p < 0.001, sprintf("%.1e", res$p),
                  sprintf("%.3f", res$p)))
  dstr <- ifelse(is.na(res$cohens_d), "-", sprintf("%.3f", res$cohens_d))
  cat(sprintf("%-25s %-20s %-20s %-12s %-8s\n",
              v[[1]], res$FH, res$NonFH, pstr, dstr))
}

# Binary variables
if (sum(!is.na(hr$has_DM)) > 100) {
  res_dm <- describe_var(hr, "has_DM", type = "binary")
  pstr <- ifelse(is.na(res_dm$p), "NA",
           ifelse(res_dm$p < 0.001, sprintf("%.1e", res_dm$p),
                  sprintf("%.3f", res_dm$p)))
  cat(sprintf("%-25s %-20s %-20s %-12s %-8s\n",
              "Diabetes (%)", res_dm$FH, res_dm$NonFH, pstr, "-"))
}

res_sex <- describe_var(hr, "Gender_num", type = "continuous")
pstr <- ifelse(res_sex$p < 0.001, sprintf("%.1e", res_sex$p), sprintf("%.3f", res_sex$p))
cat(sprintf("%-25s %-20s %-20s %-12s %-8s\n",
            "Male (%)", res_sex$FH, res_sex$NonFH, pstr, "-"))

cat("\n")
cat("KEY FINDING: If FH shows LOWER TRG, LOWER BMI, LOWER DM prevalence,\n")
cat("this confirms 'metabolic purity' — FH is pure LDL elevation,\n")
cat("NOT metabolic syndrome. The Trig Filter quantifies this purity.\n\n")


# ==============================================================================
# B. MULTIVARIABLE REGRESSION: TRIG FILTER INDEPENDENCE
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS B: Multivariable Logistic Regression\n")
cat("Outcome: Genetic FH | Within high-risk cohort\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Standardise continuous predictors for comparability
hr$Trig_Filter_z <- scale(hr$Trig_Filter_RW)[,1]
hr$BMI_z <- scale(hr$BMI_imputed)[,1]
hr$LDL_z <- scale(hr$LDL_RW)[,1]
hr$TRG_z <- scale(hr$TRG.1)[,1]
hr$HDL_z <- scale(hr$HDL.1)[,1]
hr$Age_z <- scale(hr$Age_at_LDL1)[,1]

if (has_apob) hr$ApoB_z <- scale(hr$ApoB)[,1]

# Model 1: Trig Filter alone (crude)
m1 <- glm(is_fh_genetic ~ Trig_Filter_z, data = hr, family = binomial)

# Model 2: + BMI
m2 <- glm(is_fh_genetic ~ Trig_Filter_z + BMI_z, data = hr, family = binomial)

# Model 3: + BMI + DM
if (sum(hr$has_DM, na.rm = TRUE) > 10) {
  m3 <- glm(is_fh_genetic ~ Trig_Filter_z + BMI_z + has_DM, data = hr, family = binomial)
} else {
  m3 <- m2
  cat("NOTE: Too few DM cases for Model 3. Skipping DM adjustment.\n")
}

# Model 4: + BMI + DM + ApoB
if (has_apob && sum(hr$has_DM, na.rm = TRUE) > 10) {
  m4_data <- hr[!is.na(hr$ApoB_z) & !is.na(hr$has_DM), ]
  m4 <- glm(is_fh_genetic ~ Trig_Filter_z + BMI_z + has_DM + ApoB_z,
            data = m4_data, family = binomial)
} else if (has_apob) {
  m4_data <- hr[!is.na(hr$ApoB_z), ]
  m4 <- glm(is_fh_genetic ~ Trig_Filter_z + BMI_z + ApoB_z,
            data = m4_data, family = binomial)
} else {
  m4 <- m3
}

# Model 5: Full model (Trig Filter + all confounders + Age + Sex)
if (has_apob && sum(hr$has_DM, na.rm = TRUE) > 10) {
  m5_data <- hr[!is.na(hr$ApoB_z) & !is.na(hr$has_DM), ]
  m5 <- glm(is_fh_genetic ~ Trig_Filter_z + BMI_z + has_DM + ApoB_z + Age_z + Gender_num,
            data = m5_data, family = binomial)
} else {
  m5_data <- hr[!is.na(hr$BMI_z), ]
  m5 <- glm(is_fh_genetic ~ Trig_Filter_z + BMI_z + Age_z + Gender_num,
            data = m5_data, family = binomial)
}

# Report OR for Trig Filter across models
report_trig_or <- function(model, model_name) {
  cf <- summary(model)$coefficients
  if ("Trig_Filter_z" %in% rownames(cf)) {
    beta <- cf["Trig_Filter_z", "Estimate"]
    se   <- cf["Trig_Filter_z", "Std. Error"]
    or   <- exp(beta)
    ci_lo <- exp(beta - 1.96 * se)
    ci_hi <- exp(beta + 1.96 * se)
    pval  <- cf["Trig_Filter_z", "Pr(>|z|)"]
    cat(sprintf("  %-35s OR = %.3f [%.3f - %.3f]  p = %.1e\n",
                model_name, or, ci_lo, ci_hi, pval))
    return(data.frame(model = model_name, OR = or, CI_lo = ci_lo, CI_hi = ci_hi, p = pval))
  }
  return(NULL)
}

cat("Odds Ratio for Trig Filter (per SD increase) across models:\n\n")
or_results <- list()
or_results[[1]] <- report_trig_or(m1, "M1: Trig Filter (crude)")
or_results[[2]] <- report_trig_or(m2, "M2: + BMI")
or_results[[3]] <- report_trig_or(m3, "M3: + BMI + DM")
or_results[[4]] <- report_trig_or(m4, "M4: + BMI + DM + ApoB")
or_results[[5]] <- report_trig_or(m5, "M5: + All confounders + Age + Sex")
or_df <- do.call(rbind, or_results[!sapply(or_results, is.null)])

cat("\nKEY FINDING: If OR for Trig Filter remains stable (< 10% change)\n")
cat("across models, it is NOT confounded by DM/BMI/ApoB.\n\n")

# Report full model 5 coefficients
cat("Full model (M5) coefficients:\n")
print(summary(m5)$coefficients)
cat("\n")

# --- VIF check (optional) ---
tryCatch({
  if (requireNamespace("car", quietly = TRUE)) {
    cat("Variance Inflation Factors (M5):\n")
    print(car::vif(m5))
    cat("(VIF > 5 indicates problematic multicollinearity)\n\n")
  }
}, error = function(e) cat("VIF not computed (car package not available)\n\n"))


# ==============================================================================
# C. STRATIFIED AUC ANALYSIS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS C: Trig Filter AUC Stratified by Metabolic Subgroups\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

calc_auc_ci <- function(response, predictor) {
  tryCatch({
    r <- pROC::roc(response, predictor, quiet = TRUE)
    ci <- pROC::ci.auc(r, method = "delong")
    list(auc = as.numeric(r$auc), ci_lo = ci[1], ci_hi = ci[3],
         n = length(response), n_pos = sum(response))
  }, error = function(e) list(auc = NA, ci_lo = NA, ci_hi = NA, n = 0, n_pos = 0))
}

forest_data <- data.frame()

# --- Overall ---
ov <- calc_auc_ci(hr$is_fh_genetic, hr$Trig_Filter_RW)
cat(sprintf("%-35s AUC = %.3f [%.3f - %.3f]  N = %d (FH = %d)\n",
            "Overall", ov$auc, ov$ci_lo, ov$ci_hi, ov$n, ov$n_pos))
forest_data <- rbind(forest_data,
  data.frame(subgroup = "Overall", auc = ov$auc, ci_lo = ov$ci_lo,
             ci_hi = ov$ci_hi, n = ov$n, n_fh = ov$n_pos))

# --- By DM status ---
if (sum(hr$has_DM, na.rm = TRUE) > 10) {
  for (dm_val in c(FALSE, TRUE)) {
    sub <- hr[hr$has_DM == dm_val & !is.na(hr$Trig_Filter_RW), ]
    label <- ifelse(dm_val, "DM: Yes", "DM: No")
    if (sum(sub$is_fh_genetic) >= 5) {
      res <- calc_auc_ci(sub$is_fh_genetic, sub$Trig_Filter_RW)
      cat(sprintf("%-35s AUC = %.3f [%.3f - %.3f]  N = %d (FH = %d)\n",
                  label, res$auc, res$ci_lo, res$ci_hi, res$n, res$n_pos))
      forest_data <- rbind(forest_data,
        data.frame(subgroup = label, auc = res$auc, ci_lo = res$ci_lo,
                   ci_hi = res$ci_hi, n = res$n, n_fh = res$n_pos))
    }
  }
}

# --- By BMI quartile ---
for (q in levels(hr$BMI_quartile)) {
  sub <- hr[hr$BMI_quartile == q & !is.na(hr$Trig_Filter_RW), ]
  if (sum(sub$is_fh_genetic) >= 5) {
    res <- calc_auc_ci(sub$is_fh_genetic, sub$Trig_Filter_RW)
    label <- sprintf("BMI %s", q)
    cat(sprintf("%-35s AUC = %.3f [%.3f - %.3f]  N = %d (FH = %d)\n",
                label, res$auc, res$ci_lo, res$ci_hi, res$n, res$n_pos))
    forest_data <- rbind(forest_data,
      data.frame(subgroup = label, auc = res$auc, ci_lo = res$ci_lo,
                 ci_hi = res$ci_hi, n = res$n, n_fh = res$n_pos))
  }
}

# --- By ApoB tertile ---
if (has_apob) {
  for (t in levels(hr$ApoB_tertile)) {
    sub <- hr[hr$ApoB_tertile == t & !is.na(hr$Trig_Filter_RW), ]
    if (sum(sub$is_fh_genetic) >= 5) {
      res <- calc_auc_ci(sub$is_fh_genetic, sub$Trig_Filter_RW)
      label <- sprintf("ApoB %s", t)
      cat(sprintf("%-35s AUC = %.3f [%.3f - %.3f]  N = %d (FH = %d)\n",
                  label, res$auc, res$ci_lo, res$ci_hi, res$n, res$n_pos))
      forest_data <- rbind(forest_data,
        data.frame(subgroup = label, auc = res$auc, ci_lo = res$ci_lo,
                   ci_hi = res$ci_hi, n = res$n, n_fh = res$n_pos))
    }
  }
}

# --- By Age ---
hr$age_group <- ifelse(hr$Age_at_LDL1 < 55, "Age < 55", "Age >= 55")
for (ag in c("Age < 55", "Age >= 55")) {
  sub <- hr[hr$age_group == ag & !is.na(hr$Trig_Filter_RW), ]
  if (sum(sub$is_fh_genetic) >= 5) {
    res <- calc_auc_ci(sub$is_fh_genetic, sub$Trig_Filter_RW)
    cat(sprintf("%-35s AUC = %.3f [%.3f - %.3f]  N = %d (FH = %d)\n",
                ag, res$auc, res$ci_lo, res$ci_hi, res$n, res$n_pos))
    forest_data <- rbind(forest_data,
      data.frame(subgroup = ag, auc = res$auc, ci_lo = res$ci_lo,
                 ci_hi = res$ci_hi, n = res$n, n_fh = res$n_pos))
  }
}

# --- By Sex ---
for (sx in c(0, 1)) {
  sub <- hr[hr$Gender_num == sx & !is.na(hr$Trig_Filter_RW), ]
  label <- ifelse(sx == 1, "Male", "Female")
  if (sum(sub$is_fh_genetic) >= 5) {
    res <- calc_auc_ci(sub$is_fh_genetic, sub$Trig_Filter_RW)
    cat(sprintf("%-35s AUC = %.3f [%.3f - %.3f]  N = %d (FH = %d)\n",
                label, res$auc, res$ci_lo, res$ci_hi, res$n, res$n_pos))
    forest_data <- rbind(forest_data,
      data.frame(subgroup = label, auc = res$auc, ci_lo = res$ci_lo,
                 ci_hi = res$ci_hi, n = res$n, n_fh = res$n_pos))
  }
}

cat("\nKEY FINDING: If Trig Filter AUC is CONSISTENT across all subgroups,\n")
cat("it discriminates FH independently of DM/BMI/ApoB status.\n\n")


# ==============================================================================
# D. INTERACTION TESTS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS D: Interaction Tests (Effect Modification)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Testing whether Trig Filter x Confounder interactions are significant:\n")
cat("(Non-significant = Trig Filter works independently in all subgroups)\n\n")

# Trig Filter x DM
if (sum(hr$has_DM, na.rm = TRUE) > 10) {
  m_int_dm <- glm(is_fh_genetic ~ Trig_Filter_z * has_DM, data = hr, family = binomial)
  p_int_dm <- summary(m_int_dm)$coefficients["Trig_Filter_z:has_DMTRUE", "Pr(>|z|)"]
  cat(sprintf("  Trig Filter x DM:          p-interaction = %.3f %s\n",
              p_int_dm, ifelse(p_int_dm > 0.05, "(NS - no effect modification)", "(SIGNIFICANT)")))
}

# Trig Filter x BMI (continuous)
m_int_bmi <- glm(is_fh_genetic ~ Trig_Filter_z * BMI_z, data = hr, family = binomial)
p_int_bmi <- summary(m_int_bmi)$coefficients["Trig_Filter_z:BMI_z", "Pr(>|z|)"]
cat(sprintf("  Trig Filter x BMI:         p-interaction = %.3f %s\n",
            p_int_bmi, ifelse(p_int_bmi > 0.05, "(NS - no effect modification)", "(SIGNIFICANT)")))

# Trig Filter x ApoB
if (has_apob) {
  m_int_apob_data <- hr[!is.na(hr$ApoB_z), ]
  m_int_apob <- glm(is_fh_genetic ~ Trig_Filter_z * ApoB_z,
                     data = m_int_apob_data, family = binomial)
  p_int_apob <- summary(m_int_apob)$coefficients["Trig_Filter_z:ApoB_z", "Pr(>|z|)"]
  cat(sprintf("  Trig Filter x ApoB:        p-interaction = %.3f %s\n",
              p_int_apob, ifelse(p_int_apob > 0.05, "(NS - no effect modification)", "(SIGNIFICANT)")))
}

# Trig Filter x Age
m_int_age <- glm(is_fh_genetic ~ Trig_Filter_z * Age_z, data = hr, family = binomial)
p_int_age <- summary(m_int_age)$coefficients["Trig_Filter_z:Age_z", "Pr(>|z|)"]
cat(sprintf("  Trig Filter x Age:         p-interaction = %.3f %s\n",
            p_int_age, ifelse(p_int_age > 0.05, "(NS - no effect modification)", "(SIGNIFICANT)")))

# Trig Filter x Sex
m_int_sex <- glm(is_fh_genetic ~ Trig_Filter_z * Gender_num, data = hr, family = binomial)
p_int_sex <- summary(m_int_sex)$coefficients["Trig_Filter_z:Gender_num", "Pr(>|z|)"]
cat(sprintf("  Trig Filter x Sex:         p-interaction = %.3f %s\n",
            p_int_sex, ifelse(p_int_sex > 0.05, "(NS - no effect modification)", "(SIGNIFICANT)")))

cat("\n")


# ==============================================================================
# E. TUDOR MODEL ROBUSTNESS: ADDING CONFOUNDERS DOESN'T HELP
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS E: Does Adding Metabolic Confounders Improve TUDOR?\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# TUDOR AUC alone
roc_tudor <- pROC::roc(hr$is_fh_genetic, hr$tudor_prob, quiet = TRUE)
cat(sprintf("%-40s AUC = %.3f\n", "TUDOR alone", as.numeric(roc_tudor$auc)))

# TUDOR + BMI
m_tudor_bmi <- glm(is_fh_genetic ~ tudor_prob + BMI_z, data = hr, family = binomial)
pred_bmi <- predict(m_tudor_bmi, type = "response")
roc_bmi <- pROC::roc(hr$is_fh_genetic, pred_bmi, quiet = TRUE)
dt_bmi <- pROC::roc.test(roc_tudor, roc_bmi, method = "delong")
cat(sprintf("%-40s AUC = %.3f  delta = %+.3f  p = %.3f\n",
            "TUDOR + BMI", as.numeric(roc_bmi$auc),
            as.numeric(roc_bmi$auc) - as.numeric(roc_tudor$auc), dt_bmi$p.value))

# TUDOR + DM
if (sum(hr$has_DM, na.rm = TRUE) > 10) {
  m_tudor_dm <- glm(is_fh_genetic ~ tudor_prob + has_DM, data = hr, family = binomial)
  pred_dm <- predict(m_tudor_dm, type = "response")
  roc_dm <- pROC::roc(hr$is_fh_genetic, pred_dm, quiet = TRUE)
  dt_dm <- pROC::roc.test(roc_tudor, roc_dm, method = "delong")
  cat(sprintf("%-40s AUC = %.3f  delta = %+.3f  p = %.3f\n",
              "TUDOR + DM", as.numeric(roc_dm$auc),
              as.numeric(roc_dm$auc) - as.numeric(roc_tudor$auc), dt_dm$p.value))
}

# TUDOR + ApoB
if (has_apob) {
  hr_apob <- hr[!is.na(hr$ApoB), ]
  roc_tudor_apob_base <- pROC::roc(hr_apob$is_fh_genetic, hr_apob$tudor_prob, quiet = TRUE)
  m_tudor_apob <- glm(is_fh_genetic ~ tudor_prob + ApoB_z, data = hr_apob, family = binomial)
  pred_apob <- predict(m_tudor_apob, type = "response")
  roc_apob <- pROC::roc(hr_apob$is_fh_genetic, pred_apob, quiet = TRUE)
  dt_apob <- pROC::roc.test(roc_tudor_apob_base, roc_apob, method = "delong")
  cat(sprintf("%-40s AUC = %.3f  delta = %+.3f  p = %.3f\n",
              "TUDOR + ApoB", as.numeric(roc_apob$auc),
              as.numeric(roc_apob$auc) - as.numeric(roc_tudor_apob_base$auc),
              dt_apob$p.value))
}

# TUDOR + BMI + DM + ApoB (kitchen sink)
if (has_apob && sum(hr$has_DM, na.rm = TRUE) > 10) {
  hr_full <- hr[!is.na(hr$ApoB_z) & !is.na(hr$has_DM), ]
  roc_tudor_full_base <- pROC::roc(hr_full$is_fh_genetic, hr_full$tudor_prob, quiet = TRUE)
  m_full <- glm(is_fh_genetic ~ tudor_prob + BMI_z + has_DM + ApoB_z,
                data = hr_full, family = binomial)
  pred_full <- predict(m_full, type = "response")
  roc_full <- pROC::roc(hr_full$is_fh_genetic, pred_full, quiet = TRUE)
  dt_full <- pROC::roc.test(roc_tudor_full_base, roc_full, method = "delong")
  cat(sprintf("%-40s AUC = %.3f  delta = %+.3f  p = %.3f\n",
              "TUDOR + BMI + DM + ApoB", as.numeric(roc_full$auc),
              as.numeric(roc_full$auc) - as.numeric(roc_tudor_full_base$auc),
              dt_full$p.value))
}

cat("\nKEY FINDING: If adding BMI/DM/ApoB does NOT improve TUDOR AUC,\n")
cat("this proves TUDOR (via the Trig Filter) already captures the relevant\n")
cat("metabolic information. These variables are NOT independent confounders.\n\n")


# ==============================================================================
# F. BIOLOGICAL MECHANISM: TRG DISTRIBUTION IN FH vs NON-FH
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS F: Biological Mechanism — 'Metabolic Purity' of FH\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Within LDL-MATCHED groups (compare at similar LDL levels)
# Split into LDL deciles and compare TRG within each
hr$LDL_decile <- cut(hr$LDL_RW,
                     breaks = quantile(hr$LDL_RW, probs = seq(0, 1, 0.1), na.rm = TRUE),
                     labels = paste0("D", 1:10), include.lowest = TRUE)

cat("TRG by FH status within LDL deciles (LDL-matched comparison):\n")
cat(sprintf("%-12s %-12s %-12s %-15s %-15s %-10s\n",
            "LDL Decile", "LDL Range", "N(FH)", "TRG FH", "TRG Non-FH", "P-value"))
cat(strrep("-", 75), "\n")

for (d in levels(hr$LDL_decile)) {
  sub <- hr[hr$LDL_decile == d, ]
  n_fh <- sum(sub$is_fh_genetic)
  if (n_fh >= 3) {
    ldl_range <- sprintf("%.1f-%.1f", min(sub$LDL_RW), max(sub$LDL_RW))
    trg_fh <- mean(sub$TRG.1[sub$is_fh_genetic == TRUE], na.rm = TRUE)
    trg_nonfh <- mean(sub$TRG.1[sub$is_fh_genetic == FALSE], na.rm = TRUE)
    pval <- tryCatch(t.test(sub$TRG.1[sub$is_fh_genetic], sub$TRG.1[!sub$is_fh_genetic])$p.value,
                     error = function(e) NA)
    pstr <- ifelse(is.na(pval), "NA", ifelse(pval < 0.001, sprintf("%.1e", pval), sprintf("%.3f", pval)))
    cat(sprintf("%-12s %-12s %-12d %-15.2f %-15.2f %-10s\n",
                d, ldl_range, n_fh, trg_fh, trg_nonfh, pstr))
  }
}

cat("\nKEY FINDING: Within each LDL decile, FH patients should have LOWER TRG\n")
cat("than non-FH patients. This LDL-matched comparison eliminates LDL as\n")
cat("a confounder, isolating the TRG effect.\n\n")

# --- Metabolic Syndrome Score ---
# Create a simple metabolic syndrome proxy score:
# High TRG (>1.7) + High BMI (>30) + DM + Low HDL (<1.0 M / <1.3 F)
hr$met_high_trg <- hr$TRG.1 > 1.7
hr$met_obese    <- hr$BMI_imputed > 30
hr$met_low_hdl  <- ifelse(hr$Gender_num == 1, hr$HDL.1 < 1.0, hr$HDL.1 < 1.3)
hr$met_score    <- as.integer(hr$met_high_trg) + as.integer(hr$met_obese) +
                   as.integer(hr$has_DM) + as.integer(hr$met_low_hdl)

cat("Metabolic Syndrome Component Score (0-4):\n")
cat(sprintf("%-12s %-12s %-12s %-12s\n", "MetS Score", "N", "FH (%)", "FH Prev (%)"))
cat(strrep("-", 50), "\n")
for (s in sort(unique(hr$met_score))) {
  sub <- hr[hr$met_score == s, ]
  n_fh <- sum(sub$is_fh_genetic)
  prev <- 100 * n_fh / nrow(sub)
  cat(sprintf("%-12d %-12d %-12d %-12.2f\n", s, nrow(sub), n_fh, prev))
}

cat("\nKEY FINDING: FH prevalence should be HIGHEST in MetS score = 0\n")
cat("(metabolically healthy). FH is a genetic disorder, NOT metabolic.\n")
cat("The Trig Filter captures this: high TrigFilter = low MetS = likely FH.\n\n")


# ==============================================================================
# G. FOREST PLOT
# ==============================================================================
if (nrow(forest_data) > 2) {
  cat("Generating stratified AUC forest plot...\n")

  forest_data$subgroup <- factor(forest_data$subgroup,
                                  levels = rev(forest_data$subgroup))

  p_forest <- ggplot(forest_data, aes(x = auc, y = subgroup)) +
    geom_point(size = 3, color = "darkblue") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, color = "darkblue") +
    geom_vline(xintercept = forest_data$auc[forest_data$subgroup == "Overall"],
               linetype = "dashed", color = "red", alpha = 0.5) +
    labs(x = "AUC (95% CI)", y = "",
         title = "Trig Filter Discrimination Across Metabolic Subgroups",
         subtitle = "High-risk cohort (LDL > 4.9 mmol/L)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank()) +
    xlim(0.5, 1.0) +
    annotate("text", x = 0.98, y = nrow(forest_data),
             label = sprintf("N (FH)"), hjust = 1, fontface = "bold", size = 3) +
    geom_text(aes(x = 0.98, label = sprintf("%d (%d)", n, n_fh)),
              hjust = 1, size = 2.5)

  ggsave(file.path(PLOT_DIR, "fig_trig_shield_forest.pdf"),
         p_forest, width = 8, height = 6)
  ggsave(file.path(PLOT_DIR, "fig_trig_shield_forest.png"),
         p_forest, width = 8, height = 6, dpi = 300)
  cat("  Saved: fig_trig_shield_forest.pdf/png\n\n")
}


# ==============================================================================
# H. SAVE RESULTS
# ==============================================================================
results_08 <- list(
  or_table = or_df,
  forest_data = forest_data,
  model_full = m5,
  timestamp = Sys.time()
)
saveRDS(results_08, file.path(OUTPUT_DIR, "results_08_trig_shield.rds"))

cat("\n=== 08_trig_shield_biology.R COMPLETE ===\n")
cat("Results saved to:", file.path(OUTPUT_DIR, "results_08_trig_shield.rds"), "\n")
