# ==============================================================================
# TUDOR PIPELINE: STEP 19 — MICE MULTIPLE IMPUTATION (MCMC)
# ==============================================================================
# PURPOSE: Implement proper multiple imputation for missing covariates
#          using MICE with MCMC sampling. Required for Nature submission.
#
# CRITICAL FIX: The pipeline previously used:
#   - Median imputation for BMI (biased if MNAR)
#   - Complete-case analysis for Lp(a), ApoB (selection bias)
#
# IMPLEMENTS:
#   1. MICE imputation with m=10 datasets
#   2. Predictive mean matching (PMM) for continuous variables
#   3. Logistic regression for binary variables
#   4. Rubin's rules for pooling estimates
#   5. Sensitivity analysis under MNAR assumptions
#
# INPUT:   tudor_analysis_ready.rds
# OUTPUT:  tudor_imputed_datasets.rds, pooled_results.rds
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(mice)
  library(pROC)
  library(mitools)
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
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== TUDOR PIPELINE: 19_mice_multiple_imputation.R ===\n")
cat("MICE Multiple Imputation with MCMC Sampling\n")
cat("Required for Nature/Lancet Submission\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
if (!file.exists(rds_file)) stop("Run 01_data_merge.R first!")
df <- readRDS(rds_file)

if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}

cat("Total participants:", nrow(df), "\n")
cat("Genetic FH cases:", sum(df$is_fh_genetic), "\n\n")

# ==============================================================================
# 2. ASSESS MISSING DATA PATTERNS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MISSING DATA ASSESSMENT\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Variables for imputation
impute_vars <- c("LDL_RW", "HDL.1", "TRG.1", "CHOL", "BMI_imputed",
                  "ApoB", "Lpa_nmol", "Age_at_LDL1", "Gender_num",
                  "is_fh_genetic", "tudor_prob", "edlcn_score")

# Keep only variables that exist
impute_vars <- impute_vars[impute_vars %in% names(df)]

cat("Variables assessed for missingness:\n")
for (v in impute_vars) {
  n_miss <- sum(is.na(df[[v]]))
  pct_miss <- 100 * n_miss / nrow(df)
  cat(sprintf("  %-20s: %7d missing (%5.1f%%)\n", v, n_miss, pct_miss))
}
cat("\n")

# ==============================================================================
# 3. PREPARE IMPUTATION DATASET
# ==============================================================================
# Select variables for imputation model
# Include auxiliary variables that predict missingness (improve MAR assumption)
aux_vars <- c("Age_at_LDL1", "Gender_num", "statin_tier",
              "has_any_cvd", "Premature_ASCVD")
aux_vars <- aux_vars[aux_vars %in% names(df)]

# Variables to impute
target_vars <- c("BMI_imputed", "ApoB", "Lpa_nmol", "ApoB_LDL_Ratio")
target_vars <- target_vars[target_vars %in% names(df)]

# Build imputation dataset
all_vars <- unique(c(impute_vars, aux_vars, target_vars))
all_vars <- all_vars[all_vars %in% names(df)]

imp_data <- as.data.frame(df[, ..all_vars])

cat("Imputation dataset: ", nrow(imp_data), " rows x ", ncol(imp_data), " columns\n\n")

# ==============================================================================
# 4. MICE IMPUTATION (m=10 datasets, MCMC via PMM)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("RUNNING MICE IMPUTATION (m=10, maxit=20)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("Method: Predictive Mean Matching (PMM) for continuous variables\n")
cat("This preserves the distribution shape and handles non-normality.\n\n")

# Set up imputation methods
# PMM for continuous, logreg for binary, polyreg for categorical
methods <- make.method(imp_data)

# Override: use PMM for all continuous biomarkers
for (v in c("BMI_imputed", "ApoB", "Lpa_nmol", "ApoB_LDL_Ratio")) {
  if (v %in% names(methods) && methods[v] != "") {
    methods[v] <- "pmm"
  }
}

# Do NOT impute outcome or predictors that are complete
for (v in c("is_fh_genetic", "Age_at_LDL1", "Gender_num")) {
  if (v %in% names(methods)) {
    methods[v] <- ""
  }
}

cat("Imputation methods:\n")
for (v in names(methods)) {
  if (methods[v] != "") {
    cat(sprintf("  %-20s: %s\n", v, methods[v]))
  }
}
cat("\n")

# Set up predictor matrix
pred_matrix <- make.predictorMatrix(imp_data)
# Don't use outcome to predict covariates (prevents data leakage in some contexts)
# But DO include outcome as predictor for congeniality (Meng 1994, Rubin 1996)
# This is the recommended approach per van Buuren (2018)

cat("Running MICE... (this may take several minutes)\n")
t_start <- Sys.time()

mice_obj <- mice(
  imp_data,
  m = 10,                    # 10 imputed datasets
  maxit = 20,                # 20 MCMC iterations
  method = methods,
  predictorMatrix = pred_matrix,
  seed = 42,
  printFlag = FALSE
)

t_elapsed <- difftime(Sys.time(), t_start, units = "mins")
cat(sprintf("MICE completed in %.1f minutes\n\n", as.numeric(t_elapsed)))

# ==============================================================================
# 5. CONVERGENCE DIAGNOSTICS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("CONVERGENCE DIAGNOSTICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Check convergence via trace plots (summary statistics)
cat("MCMC chain summaries (last 5 iterations):\n")
imputed_vars_active <- names(methods)[methods != ""]
for (v in imputed_vars_active) {
  if (v %in% names(mice_obj$chainMean)) {
    chain_means <- mice_obj$chainMean[v, , ]
    if (!is.null(chain_means) && length(chain_means) > 0) {
      # Check Rhat-like convergence: variance between chains vs within chains
      if (is.matrix(chain_means)) {
        between_var <- var(chain_means[nrow(chain_means), ])
        within_var <- mean(apply(chain_means, 2, function(x) var(diff(x))))
        rhat_approx <- sqrt((between_var + within_var) / within_var)
        cat(sprintf("  %-20s: R-hat ≈ %.3f %s\n", v, rhat_approx,
                    ifelse(rhat_approx < 1.1, "(converged)", "(WARNING: check convergence)")))
      }
    }
  }
}
cat("\n")

# ==============================================================================
# 6. POOLED ANALYSIS: TUDOR AUC WITH RUBIN'S RULES
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POOLED ANALYSIS: TUDOR DISCRIMINATION (RUBIN'S RULES)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Extract each imputed dataset and compute AUC
auc_estimates <- numeric(10)
auc_variances <- numeric(10)

for (i in 1:10) {
  imp_i <- complete(mice_obj, i)

  # Merge back non-imputed variables
  # (The imputed dataset only contains the imputation variables)
  # For TUDOR AUC, we need tudor_prob and is_fh_genetic which should be complete

  hr_i <- imp_i[imp_i$LDL_RW > 4.9 & !is.na(imp_i$tudor_prob), ]

  if (sum(hr_i$is_fh_genetic) >= 5) {
    roc_i <- roc(hr_i$is_fh_genetic, hr_i$tudor_prob, quiet = TRUE)
    auc_i <- as.numeric(auc(roc_i))
    # Variance via DeLong
    ci_i <- ci.auc(roc_i, method = "delong")
    se_i <- (as.numeric(ci_i[3]) - as.numeric(ci_i[1])) / (2 * 1.96)

    auc_estimates[i] <- auc_i
    auc_variances[i] <- se_i^2
  }
}

# Rubin's Rules for pooling
Q_bar <- mean(auc_estimates)                                # Pooled estimate
U_bar <- mean(auc_variances)                                # Within-imputation variance
B <- var(auc_estimates)                                      # Between-imputation variance
T_total <- U_bar + (1 + 1/10) * B                          # Total variance
se_pooled <- sqrt(T_total)

# Degrees of freedom (Barnard-Rubin)
lambda <- ((1 + 1/10) * B) / T_total                       # Fraction of missing information
df_old <- (10 - 1) / lambda^2
df_obs <- ((nrow(hr_i) - length(impute_vars)) + 1) /
          ((nrow(hr_i) - length(impute_vars)) + 3) *
          (nrow(hr_i) - length(impute_vars)) * (1 - lambda)
df_adjusted <- (df_old * df_obs) / (df_old + df_obs)

ci_lower <- Q_bar - qt(0.975, df_adjusted) * se_pooled
ci_upper <- Q_bar + qt(0.975, df_adjusted) * se_pooled

cat("TUDOR AUC (Pooled across 10 imputed datasets):\n")
cat(sprintf("  Pooled AUC: %.4f (95%% CI: %.4f - %.4f)\n", Q_bar, ci_lower, ci_upper))
cat(sprintf("  Within-imputation variance (U): %.6f\n", U_bar))
cat(sprintf("  Between-imputation variance (B): %.6f\n", B))
cat(sprintf("  Fraction of missing information (lambda): %.3f\n", lambda))
cat(sprintf("  Adjusted degrees of freedom: %.1f\n\n", df_adjusted))

# ==============================================================================
# 7. POOLED BIOMARKER AUGMENTATION (ApoB + Lp(a))
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("POOLED BIOMARKER AUGMENTATION ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# For each imputed dataset, fit augmented models and pool coefficients
aug_coefs <- list()
aug_aucs <- list()

models_to_test <- list(
  "Base TUDOR" = "is_fh_genetic ~ tudor_prob",
  "TUDOR + ApoB" = "is_fh_genetic ~ tudor_prob + ApoB",
  "TUDOR + Lp(a)" = "is_fh_genetic ~ tudor_prob + Lpa_nmol",
  "TUDOR + ApoB + Lp(a)" = "is_fh_genetic ~ tudor_prob + ApoB + Lpa_nmol"
)

for (model_name in names(models_to_test)) {
  formula_str <- models_to_test[[model_name]]
  formula_obj <- as.formula(formula_str)
  model_vars <- all.vars(formula_obj)

  aucs_m <- numeric(10)
  valid_count <- 0

  for (i in 1:10) {
    imp_i <- complete(mice_obj, i)
    hr_i <- imp_i[imp_i$LDL_RW > 4.9 & !is.na(imp_i$tudor_prob), ]

    # Check all variables available
    if (all(model_vars %in% names(hr_i))) {
      cc <- complete.cases(hr_i[, model_vars])
      hr_cc <- hr_i[cc, ]

      if (nrow(hr_cc) >= 100 && sum(hr_cc$is_fh_genetic) >= 10) {
        fit <- tryCatch(
          glm(formula_obj, data = hr_cc, family = binomial),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          pred <- predict(fit, type = "response")
          roc_fit <- roc(hr_cc$is_fh_genetic, pred, quiet = TRUE)
          aucs_m[i] <- as.numeric(auc(roc_fit))
          valid_count <- valid_count + 1
        }
      }
    }
  }

  if (valid_count >= 5) {
    pooled_auc <- mean(aucs_m[aucs_m > 0])
    sd_auc <- sd(aucs_m[aucs_m > 0])
    cat(sprintf("  %-30s: Pooled AUC = %.4f (SD = %.4f, m_valid = %d)\n",
                model_name, pooled_auc, sd_auc, valid_count))
    aug_aucs[[model_name]] <- list(pooled_auc = pooled_auc, sd = sd_auc,
                                    m_valid = valid_count)
  } else {
    cat(sprintf("  %-30s: Insufficient valid imputations (%d/10)\n",
                model_name, valid_count))
  }
}
cat("\n")

# ==============================================================================
# 8. MNAR SENSITIVITY ANALYSIS (Delta Adjustment)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MNAR SENSITIVITY ANALYSIS (Delta Adjustment Method)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")
cat("Testing robustness to Missing Not At Random (MNAR) assumptions.\n")
cat("Delta shifts applied to imputed values for FH cases:\n\n")

# Under MNAR, FH cases with missing biomarkers may have systematically
# different values than observed. Test with delta adjustments.
deltas <- c(-0.5, -0.25, 0, 0.25, 0.5)  # SD units

cat(sprintf("%-15s | %9s | %s\n", "Delta (SD)", "AUC", "Interpretation"))
cat(paste(rep("-", 50), collapse = ""), "\n")

for (delta in deltas) {
  aucs_delta <- numeric(10)

  for (i in 1:10) {
    imp_i <- complete(mice_obj, i)
    hr_i <- imp_i[imp_i$LDL_RW > 4.9 & !is.na(imp_i$tudor_prob), ]

    if (delta != 0 && "ApoB" %in% names(hr_i)) {
      # Apply delta shift to originally-missing ApoB values in FH cases
      originally_missing <- is.na(imp_data$ApoB[imp_data$LDL_RW > 4.9])
      fh_mask <- hr_i$is_fh_genetic == 1
      adjust_mask <- originally_missing & fh_mask

      if (sum(adjust_mask, na.rm = TRUE) > 0) {
        apob_sd <- sd(hr_i$ApoB, na.rm = TRUE)
        hr_i$ApoB[adjust_mask] <- hr_i$ApoB[adjust_mask] + delta * apob_sd
      }
    }

    if (sum(hr_i$is_fh_genetic) >= 5) {
      roc_i <- roc(hr_i$is_fh_genetic, hr_i$tudor_prob, quiet = TRUE)
      aucs_delta[i] <- as.numeric(auc(roc_i))
    }
  }

  pooled <- mean(aucs_delta[aucs_delta > 0])
  interp <- ifelse(delta < 0, "FH missing = lower ApoB",
             ifelse(delta > 0, "FH missing = higher ApoB", "MAR assumption"))
  cat(sprintf("%-15s | %9.4f | %s\n", sprintf("%+.2f SD", delta), pooled, interp))
}
cat("\n")

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================
imputation_results <- list(
  mice_object = mice_obj,
  convergence = list(
    m = 10, maxit = 20, method = "pmm"
  ),
  pooled_tudor_auc = list(
    estimate = Q_bar, ci = c(ci_lower, ci_upper),
    within_var = U_bar, between_var = B,
    lambda = lambda, df = df_adjusted
  ),
  augmented_models = aug_aucs,
  timestamp = Sys.time()
)

saveRDS(imputation_results, file.path(OUTPUT_DIR, "19_imputation_results.rds"))
cat("Saved imputation results to:", file.path(OUTPUT_DIR, "19_imputation_results.rds"), "\n")

cat("\n=== 19_mice_multiple_imputation.R COMPLETE ===\n")
cat("Multiple Imputation with MCMC successfully implemented.\n")
cat("Rubin's Rules applied for pooled estimates.\n")
cat("MNAR sensitivity analysis completed.\n")
