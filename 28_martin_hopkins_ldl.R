# ==============================================================================
# TUDOR PIPELINE: STEP 28 — MARTIN-HOPKINS LDL & ADVANCED LDL METHODS
# ==============================================================================
# PURPOSE: Address Nature reviewer concern about LDL calculation methods.
#          Implement Martin-Hopkins equation (superior to Friedewald at
#          high triglycerides) and compare TUDOR performance across methods.
#
# BACKGROUND:
#   Friedewald equation: LDL = TC - HDL - (TG/2.2) [mmol/L]
#   - INVALID when TG > 4.5 mmol/L
#   - Systematically underestimates LDL at low LDL and high TG
#
#   Martin-Hopkins equation: LDL = TC - HDL - (TG / adjustable_factor)
#   - Uses strata-specific TG:VLDL-C ratio from 1.3M samples
#   - More accurate across all TG ranges
#   - Recommended by NLA 2020 guidelines
#
# SENSITIVITY ANALYSIS:
#   Compare TUDOR AUC using:
#   1. Direct-measured LDL (UKB enzymatic assay — gold standard)
#   2. Friedewald LDL
#   3. Martin-Hopkins LDL
#   4. Sampson equation (modified Friedewald for TG up to 9.0 mmol/L)
#
# AUTHORS: Tudor Pipeline Team
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
})

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

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("TUDOR PIPELINE: 28 — MARTIN-HOPKINS LDL SENSITIVITY ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. MARTIN-HOPKINS ADJUSTABLE FACTOR TABLE
# ==============================================================================
# Source: Martin et al., JAMA 2013; 310(19):2061-2068
# Table of TG:VLDL-C ratio by TG strata and non-HDL-C strata
# Values represent the denominator in: VLDL-C = TG / factor

# Simplified strata (full table has 180 cells; using representative values)
# TG in mg/dL strata: <100, 100-149, 150-199, 200-399, >=400
# Non-HDL-C in mg/dL strata: <100, 100-129, 130-159, 160-189, 190-219, >=220

martin_hopkins_factor <- function(tg_mgdl, nonhdl_mgdl) {
  # Simplified Martin-Hopkins adjustable factor
  # Full lookup table approximation
  factor <- ifelse(tg_mgdl < 100,
    ifelse(nonhdl_mgdl < 100, 3.5,
    ifelse(nonhdl_mgdl < 130, 4.0,
    ifelse(nonhdl_mgdl < 160, 4.3,
    ifelse(nonhdl_mgdl < 190, 4.5,
    ifelse(nonhdl_mgdl < 220, 4.7, 4.9))))),
  ifelse(tg_mgdl < 150,
    ifelse(nonhdl_mgdl < 100, 4.0,
    ifelse(nonhdl_mgdl < 130, 4.4,
    ifelse(nonhdl_mgdl < 160, 4.7,
    ifelse(nonhdl_mgdl < 190, 5.0,
    ifelse(nonhdl_mgdl < 220, 5.2, 5.4))))),
  ifelse(tg_mgdl < 200,
    ifelse(nonhdl_mgdl < 100, 4.5,
    ifelse(nonhdl_mgdl < 130, 4.8,
    ifelse(nonhdl_mgdl < 160, 5.1,
    ifelse(nonhdl_mgdl < 190, 5.4,
    ifelse(nonhdl_mgdl < 220, 5.6, 5.8))))),
  ifelse(tg_mgdl < 400,
    ifelse(nonhdl_mgdl < 100, 5.0,
    ifelse(nonhdl_mgdl < 130, 5.3,
    ifelse(nonhdl_mgdl < 160, 5.6,
    ifelse(nonhdl_mgdl < 190, 5.9,
    ifelse(nonhdl_mgdl < 220, 6.1, 6.4))))),
  # TG >= 400: Martin-Hopkins recommends direct measurement
  ifelse(nonhdl_mgdl < 100, 6.0,
    ifelse(nonhdl_mgdl < 130, 6.5,
    ifelse(nonhdl_mgdl < 160, 7.0,
    ifelse(nonhdl_mgdl < 190, 7.5,
    ifelse(nonhdl_mgdl < 220, 8.0, 8.5)))))))))

  return(factor)
}

# ==============================================================================
# 2. LDL CALCULATION FUNCTIONS
# ==============================================================================

# Convert mmol/L to mg/dL: multiply by 38.67 (cholesterol) or 88.57 (TG)
mmol_to_mgdl_chol <- function(x) x * 38.67
mmol_to_mgdl_tg <- function(x) x * 88.57
mgdl_to_mmol_chol <- function(x) x / 38.67

# Friedewald equation (mmol/L)
calc_ldl_friedewald <- function(tc, hdl, tg) {
  ldl <- tc - hdl - (tg / 2.2)
  # Mark as NA when TG > 4.5 mmol/L (unreliable)
  ldl[tg > 4.5] <- NA
  return(ldl)
}

# Martin-Hopkins equation (mmol/L)
calc_ldl_martin_hopkins <- function(tc, hdl, tg) {
  # Convert to mg/dL for lookup table
  tc_mg <- mmol_to_mgdl_chol(tc)
  hdl_mg <- mmol_to_mgdl_chol(hdl)
  tg_mg <- mmol_to_mgdl_tg(tg)
  nonhdl_mg <- tc_mg - hdl_mg

  # Get adjustable factor
  adj_factor <- martin_hopkins_factor(tg_mg, nonhdl_mg)

  # Calculate VLDL-C in mg/dL
  vldl_mg <- tg_mg / adj_factor

  # LDL = TC - HDL - VLDL (in mg/dL), then convert back
  ldl_mg <- tc_mg - hdl_mg - vldl_mg
  ldl_mmol <- mgdl_to_mmol_chol(ldl_mg)

  return(ldl_mmol)
}

# Sampson equation (2020) — extended Friedewald for TG up to 9.0 mmol/L
# Sampson et al., JAMA Cardiology 2020
calc_ldl_sampson <- function(tc, hdl, tg) {
  # Convert to mg/dL
  tc_mg <- mmol_to_mgdl_chol(tc)
  hdl_mg <- mmol_to_mgdl_chol(hdl)
  tg_mg <- mmol_to_mgdl_tg(tg)

  # Sampson equation (mg/dL)
  ldl_mg <- (tc_mg / 0.948) - (hdl_mg / 0.971) -
            ((tg_mg / 8.56) + (tg_mg * (tc_mg - hdl_mg) / 2140) - (tg_mg^2 / 16100)) - 9.44
  ldl_mmol <- mgdl_to_mmol_chol(ldl_mg)

  # Mark as NA when TG > 9.0 mmol/L
  ldl_mmol[tg > 9.0] <- NA

  return(ldl_mmol)
}

# ==============================================================================
# 3. APPLY TO DATA
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")

if (file.exists(rds_file)) {
  df <- readRDS(rds_file)
  setDT(df)
  if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
    setnames(df, "participant.eid", "eid")
  }
  cat("Loaded:", nrow(df), "participants\n\n")

  # Calculate LDL using all methods
  cat("Calculating LDL-C using multiple methods...\n")
  df[, LDL_Friedewald := calc_ldl_friedewald(CHOL, HDL.1, TRG.1)]
  df[, LDL_Martin_Hopkins := calc_ldl_martin_hopkins(CHOL, HDL.1, TRG.1)]
  df[, LDL_Sampson := calc_ldl_sampson(CHOL, HDL.1, TRG.1)]

  cat(sprintf("  Direct LDL:       N = %d, median = %.2f mmol/L\n",
              sum(!is.na(df$LDL_treated)), median(df$LDL_treated, na.rm = TRUE)))
  cat(sprintf("  Friedewald:       N = %d, median = %.2f mmol/L\n",
              sum(!is.na(df$LDL_Friedewald)), median(df$LDL_Friedewald, na.rm = TRUE)))
  cat(sprintf("  Martin-Hopkins:   N = %d, median = %.2f mmol/L\n",
              sum(!is.na(df$LDL_Martin_Hopkins)), median(df$LDL_Martin_Hopkins, na.rm = TRUE)))
  cat(sprintf("  Sampson:          N = %d, median = %.2f mmol/L\n\n",
              sum(!is.na(df$LDL_Sampson)), median(df$LDL_Sampson, na.rm = TRUE)))

  # Flag participants with TG > 4.5 (where Friedewald fails)
  df[, high_tg := TRG.1 > 4.5]
  cat(sprintf("Participants with TG > 4.5 mmol/L (Friedewald unreliable): %d (%.1f%%)\n",
              sum(df$high_tg, na.rm = TRUE), 100 * mean(df$high_tg, na.rm = TRUE)))
  cat(sprintf("Participants with TG > 9.0 mmol/L (Sampson unreliable): %d\n\n",
              sum(df$TRG.1 > 9.0, na.rm = TRUE)))

  # ==============================================================================
  # 4. TUDOR RECALCULATION WITH EACH LDL METHOD
  # ==============================================================================
  TUDOR_WEIGHTS <- list(
    intercept = 0.755722, beta_LDL = 0.057911, beta_Trig = 0.492412,
    beta_HDL = -1.128045, beta_Age = -0.033393, beta_Sex = -0.088550
  )

  # Apply statin correction to each LDL method
  REDUCTION_FACTORS <- c(
    Atorvastatin = 0.38, Simvastatin = 0.35, Rosuvastatin = 0.34,
    Pravastatin = 0.25, Fluvastatin = 0.22, None = 0
  )

  correction <- REDUCTION_FACTORS[df$statin_name]
  correction[is.na(correction)] <- 0

  calc_tudor <- function(ldl_treated, correction) {
    ldl_rw <- ldl_treated / (1 - correction)
    trig_filter <- ldl_rw / (df$TRG.1 + 0.1)
    score <- TUDOR_WEIGHTS$intercept +
      TUDOR_WEIGHTS$beta_LDL * ldl_rw +
      TUDOR_WEIGHTS$beta_Trig * trig_filter +
      TUDOR_WEIGHTS$beta_HDL * df$HDL.1 +
      TUDOR_WEIGHTS$beta_Age * df$Age_at_LDL1 +
      TUDOR_WEIGHTS$beta_Sex * df$Gender_num
    prob <- 1 / (1 + exp(-score))
    return(prob)
  }

  df[, tudor_prob_friedewald := calc_tudor(LDL_Friedewald, correction)]
  df[, tudor_prob_martin := calc_tudor(LDL_Martin_Hopkins, correction)]
  df[, tudor_prob_sampson := calc_tudor(LDL_Sampson, correction)]

  # ==============================================================================
  # 5. AUC COMPARISON ACROSS METHODS
  # ==============================================================================
  cat("=" |> rep(60) |> paste(collapse = ""), "\n")
  cat("AUC COMPARISON ACROSS LDL METHODS (HIGH-RISK COHORT)\n")
  cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

  hr <- df[cohort_high_risk == TRUE]

  methods <- list(
    list("Direct (UKB enzymatic)", "tudor_prob"),
    list("Friedewald", "tudor_prob_friedewald"),
    list("Martin-Hopkins", "tudor_prob_martin"),
    list("Sampson (2020)", "tudor_prob_sampson")
  )

  cat(sprintf("%-25s | %6s | %9s | %s\n", "LDL Method", "N", "AUC", "95% CI"))
  cat(strrep("-", 65), "\n")

  roc_objects <- list()
  for (m in methods) {
    label <- m[[1]]
    col <- m[[2]]
    valid <- !is.na(hr[[col]]) & !is.na(hr$is_fh_genetic)
    if (sum(valid) > 100 && sum(hr$is_fh_genetic[valid]) >= 10) {
      r <- roc(hr$is_fh_genetic[valid], hr[[col]][valid], quiet = TRUE)
      ci <- ci.auc(r, method = "delong")
      cat(sprintf("%-25s | %6d | %9.3f | [%.3f - %.3f]\n",
                  label, sum(valid), ci[2], ci[1], ci[3]))
      roc_objects[[label]] <- r
    } else {
      cat(sprintf("%-25s | %6d | %9s | %s\n",
                  label, sum(valid), "N/A", "Insufficient cases"))
    }
  }
  cat(strrep("-", 65), "\n\n")

  # DeLong comparisons vs Direct
  if (length(roc_objects) >= 2 && "Direct (UKB enzymatic)" %in% names(roc_objects)) {
    cat("DeLong tests vs Direct LDL:\n")
    ref_roc <- roc_objects[["Direct (UKB enzymatic)"]]
    for (name in names(roc_objects)) {
      if (name == "Direct (UKB enzymatic)") next
      test <- tryCatch(roc.test(ref_roc, roc_objects[[name]], method = "delong"),
                       error = function(e) list(p.value = NA))
      cat(sprintf("  Direct vs %-20s: delta = %+.4f, p = %.2e\n",
                  name,
                  auc(ref_roc) - auc(roc_objects[[name]]),
                  test$p.value))
    }
  }

  # ==============================================================================
  # 6. CONCORDANCE ANALYSIS
  # ==============================================================================
  cat("\n--- LDL Method Agreement ---\n\n")

  # Correlation and Bland-Altman metrics
  for (method_col in c("LDL_Friedewald", "LDL_Martin_Hopkins", "LDL_Sampson")) {
    valid <- !is.na(df$LDL_treated) & !is.na(df[[method_col]])
    if (sum(valid) > 0) {
      cor_val <- cor(df$LDL_treated[valid], df[[method_col]][valid])
      diff <- df$LDL_treated[valid] - df[[method_col]][valid]
      bias <- mean(diff)
      loa <- 1.96 * sd(diff)  # Limits of agreement

      cat(sprintf("  Direct vs %-20s: r = %.4f, bias = %+.3f, LOA = +/- %.3f mmol/L\n",
                  method_col, cor_val, bias, loa))
    }
  }

  # Save results
  ldl_results <- list(
    methods = c("Direct", "Friedewald", "Martin-Hopkins", "Sampson"),
    n_high_tg = sum(df$high_tg, na.rm = TRUE),
    timestamp = Sys.time()
  )
  saveRDS(ldl_results, file.path(OUTPUT_DIR, "28_ldl_methods_results.rds"))

} else {
  cat("tudor_analysis_ready.rds not found. Run 01_data_merge.R first.\n")
  cat("This script demonstrates the Martin-Hopkins implementation.\n")
}

cat("\n=== 28_martin_hopkins_ldl.R COMPLETE ===\n")
