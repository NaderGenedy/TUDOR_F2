# ==============================================================================
# TUDOR PIPELINE: STEP 03 — LONGITUDINAL STATIN VALIDATION
# ==============================================================================
# PURPOSE: Validate the statin reduction factors by comparing:
#          - Visit 0 LDL (off statin) → Visit 1 LDL (on statin)
#          - Observed reduction vs our assumed real-world reduction factors.
#
# REQUIRES: - tudor_analysis_ready.rds from 01_data_merge.R
#           - ukb_meds_v0_a.csv and ukb_meds_v0_b.csv (Visit 0 medications)
#           - A longitudinal CSV with Visit 0 + Visit 1 LDL + medications
#
# OUTPUTS:  - Observed vs expected LDL reductions for true new statin users
#           - 03_longitudinal_results.rds
# ==============================================================================

set.seed(42)

library(data.table)
library(dplyr)

DATA_DIR <- Sys.getenv("TUDOR_DATA_DIR", unset = "")
if (DATA_DIR == "") {
  if (file.exists(file.path(getwd(), "TUDOR_UKB_Features.csv"))) {
    DATA_DIR <- getwd()
  } else {
    DATA_DIR <- "C:/Users/nader/Downloads"
  }
}
OUTPUT_DIR <- file.path(DATA_DIR, "tudor_pipeline_output")

cat("=== TUDOR PIPELINE: 03_longitudinal_validation.R ===\n\n")

# Real-world reduction factors used in the pipeline
REDUCTION_FACTORS <- c(
  Atorvastatin = 0.38, Simvastatin = 0.35, Rosuvastatin = 0.34,
  Pravastatin = 0.25, Fluvastatin = 0.22
)

# Statin UKB medication codes
STATIN_CODES <- list(
  Atorvastatin  = 1141146234,
  Simvastatin   = 1140861958,
  Rosuvastatin  = 1141192410,
  Pravastatin   = 1141146138,
  Fluvastatin   = 1141192414
)

# ==============================================================================
# 1. LOAD MAIN DATA
# ==============================================================================
df <- readRDS(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))
if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
cat("Loaded main dataset:", nrow(df), "participants\n")

# ==============================================================================
# 2. LOAD VISIT 0 MEDICATIONS (to identify statin-free at Visit 0)
# ==============================================================================
meds0a_file <- file.path(DATA_DIR, "ukb_meds_v0_a.csv")
meds0b_file <- file.path(DATA_DIR, "ukb_meds_v0_b.csv")

if (!file.exists(meds0a_file) || !file.exists(meds0b_file)) {
  cat("WARNING: Visit 0 medication files not found.\n")
  cat("  Expected: ukb_meds_v0_a.csv and ukb_meds_v0_b.csv\n")
  cat("  Run 00_extract_ukbrap.sh on UKB-RAP to extract them.\n")
  cat("  Skipping longitudinal validation.\n")

  saveRDS(list(status = "skipped", reason = "medication files not found"),
          file.path(OUTPUT_DIR, "03_longitudinal_results.rds"))
  cat("\n=== 03_longitudinal_validation.R COMPLETE (skipped) ===\n")
  quit(save = "no")
}

meds0a <- fread(meds0a_file)
meds0b <- fread(meds0b_file)

# Normalise column names (UKB-RAP uses participant.eid)
if ("participant.eid" %in% names(meds0a) && !"eid" %in% names(meds0a)) {
  setnames(meds0a, "participant.eid", "eid")
}
if ("participant.eid" %in% names(meds0b) && !"eid" %in% names(meds0b)) {
  setnames(meds0b, "participant.eid", "eid")
}

cat("Loaded Visit 0 meds Part A:", nrow(meds0a), "rows,", ncol(meds0a), "cols\n")
cat("Loaded Visit 0 meds Part B:", nrow(meds0b), "rows,", ncol(meds0b), "cols\n")

# Merge A and B
meds0 <- merge(meds0a, meds0b, by = "eid", all = TRUE)
cat("Merged Visit 0 meds:", nrow(meds0), "participants\n")

# Identify any statin at Visit 0
# Medication columns: p20003_i0_a0, p20003_i0_a1, ... (may have participant. prefix)
med_cols_v0 <- grep("p20003_i0_a", names(meds0), value = TRUE)
cat("Visit 0 medication columns found:", length(med_cols_v0), "\n")

all_statin_codes <- unlist(STATIN_CODES)

detect_statin_v0 <- function(row, cols) {
  vals <- as.numeric(unlist(row[, ..cols]))
  vals <- vals[!is.na(vals)]
  any(vals %in% all_statin_codes)
}

meds0$on_statin_v0 <- apply(meds0[, med_cols_v0, with = FALSE], 1, function(x) {
  vals <- as.numeric(x)
  vals <- vals[!is.na(vals)]
  any(vals %in% all_statin_codes)
})

cat(sprintf("Statin-free at Visit 0: %d (%.1f%%)\n",
            sum(!meds0$on_statin_v0),
            mean(!meds0$on_statin_v0) * 100))
cat(sprintf("On statin at Visit 0:   %d (%.1f%%)\n\n",
            sum(meds0$on_statin_v0),
            mean(meds0$on_statin_v0) * 100))

# ==============================================================================
# 3. LOAD LONGITUDINAL DATA (Visit 0 + Visit 1 LDL + Visit 1 meds)
# ==============================================================================
# This file should contain Visit 1 LDL (p30780_i1) and Visit 1 medications
long_file <- file.path(DATA_DIR, "ukb_longitudinal.csv")

if (!file.exists(long_file)) {
  cat("WARNING: ukb_longitudinal.csv not found.\n")
  cat("  This file should contain Visit 1 LDL and medications.\n")
  cat("  Required fields: eid, p30780_i1, p20003_i1_a* columns.\n")
  cat("  Skipping longitudinal analysis.\n")

  # Still save Visit 0 statin status
  saveRDS(list(
    status = "partial",
    reason = "longitudinal file not found",
    visit0_statin_free = sum(!meds0$on_statin_v0),
    visit0_on_statin = sum(meds0$on_statin_v0)
  ), file.path(OUTPUT_DIR, "03_longitudinal_results.rds"))

  cat("\n=== 03_longitudinal_validation.R COMPLETE (partial) ===\n")
  quit(save = "no")
}

long <- fread(long_file)
if ("participant.eid" %in% names(long) && !"eid" %in% names(long)) {
  setnames(long, "participant.eid", "eid")
}
cat("Loaded longitudinal data:", nrow(long), "rows\n")

# Identify Visit 1 LDL column
ldl_v1_col <- grep("p30780_i1", names(long), value = TRUE)[1]
if (is.null(ldl_v1_col) || length(ldl_v1_col) == 0) {
  cat("ERROR: Visit 1 LDL column (p30780_i1) not found. Stopping.\n")
  saveRDS(list(status = "error", reason = "Visit 1 LDL column not found"),
          file.path(OUTPUT_DIR, "03_longitudinal_results.rds"))
  quit(save = "no")
}

long$LDL_v1 <- as.numeric(long[[ldl_v1_col]])
cat("Visit 1 LDL non-missing:", sum(!is.na(long$LDL_v1)), "\n")

# Detect statin at Visit 1
med_cols_v1 <- grep("p20003_i1_a", names(long), value = TRUE)
cat("Visit 1 medication columns found:", length(med_cols_v1), "\n")

if (length(med_cols_v1) > 0) {
  long$on_statin_v1 <- apply(long[, med_cols_v1, with = FALSE], 1, function(x) {
    vals <- as.numeric(x)
    vals <- vals[!is.na(vals)]
    any(vals %in% all_statin_codes)
  })

  # Detect WHICH statin at Visit 1
  long$statin_v1 <- apply(long[, med_cols_v1, with = FALSE], 1, function(x) {
    vals <- as.numeric(x)
    vals <- vals[!is.na(vals)]
    for (sname in names(STATIN_CODES)) {
      if (STATIN_CODES[[sname]] %in% vals) return(sname)
    }
    return("None")
  })

  cat(sprintf("On statin at Visit 1: %d\n", sum(long$on_statin_v1, na.rm = TRUE)))
} else {
  cat("WARNING: No Visit 1 medication columns found.\n")
  long$on_statin_v1 <- FALSE
  long$statin_v1 <- "None"
}

# ==============================================================================
# 4. IDENTIFY TRUE NEW STATIN USERS (statin-free at V0, on statin at V1)
# ==============================================================================
cat("\n--- Identifying True New Statin Users ---\n")

# Merge Visit 0 statin status and Visit 1 data
analysis <- merge(
  meds0[, c("eid", "on_statin_v0")],
  long[, c("eid", "LDL_v1", "on_statin_v1", "statin_v1")],
  by = "eid"
)

# Also need Visit 0 LDL from main dataset
analysis <- merge(analysis, df[, c("eid", "LDL_treated")], by = "eid")
analysis$LDL_v0 <- analysis$LDL_treated

# Filter to true new users
new_users <- analysis[
  analysis$on_statin_v0 == FALSE &
  analysis$on_statin_v1 == TRUE &
  !is.na(analysis$LDL_v0) &
  !is.na(analysis$LDL_v1),
]

cat(sprintf("Total merged participants: %d\n", nrow(analysis)))
cat(sprintf("Statin-free at V0 AND on statin at V1: %d\n", nrow(new_users)))

if (nrow(new_users) < 20) {
  cat("WARNING: Too few new statin users for meaningful analysis.\n")
  saveRDS(list(
    status = "insufficient_data",
    n_new_users = nrow(new_users),
    visit0_statin_free = sum(!meds0$on_statin_v0),
    visit0_on_statin = sum(meds0$on_statin_v0)
  ), file.path(OUTPUT_DIR, "03_longitudinal_results.rds"))

  cat("\n=== 03_longitudinal_validation.R COMPLETE (insufficient data) ===\n")
  quit(save = "no")
}

# ==============================================================================
# 5. CALCULATE OBSERVED vs EXPECTED REDUCTIONS
# ==============================================================================
cat("\n==========================================================\n")
cat(" OBSERVED vs EXPECTED LDL REDUCTIONS\n")
cat("==========================================================\n")

new_users$LDL_change <- new_users$LDL_v1 - new_users$LDL_v0
new_users$LDL_pct_reduction <- (new_users$LDL_v0 - new_users$LDL_v1) / new_users$LDL_v0

# Expected reduction based on our factors
new_users$expected_reduction <- REDUCTION_FACTORS[new_users$statin_v1]
new_users$expected_reduction[is.na(new_users$expected_reduction)] <- 0.30  # fallback

cat(sprintf("%-15s | %5s | %12s | %12s | %12s\n",
            "Statin", "N", "Obs Redn", "Exp Redn", "Ratio"))
cat(paste(rep("-", 70), collapse = ""), "\n")

statin_summary <- list()
for (sname in names(REDUCTION_FACTORS)) {
  sub <- new_users[new_users$statin_v1 == sname, ]
  if (nrow(sub) < 5) next

  obs_median <- median(sub$LDL_pct_reduction, na.rm = TRUE)
  obs_iqr <- quantile(sub$LDL_pct_reduction, c(0.25, 0.75), na.rm = TRUE)
  expected <- REDUCTION_FACTORS[sname]
  ratio <- obs_median / expected

  cat(sprintf("%-15s | %5d | %5.1f%% [%4.1f-%4.1f] | %11.1f%% | %11.2f\n",
              sname, nrow(sub),
              obs_median * 100, obs_iqr[1] * 100, obs_iqr[2] * 100,
              expected * 100, ratio))

  statin_summary[[sname]] <- list(
    n = nrow(sub),
    observed_median = obs_median,
    observed_iqr = as.numeric(obs_iqr),
    expected = expected,
    ratio = ratio
  )
}

cat(paste(rep("-", 70), collapse = ""), "\n")

# Overall
cat(sprintf("\nOverall (N=%d): Median observed reduction = %.1f%%\n",
            nrow(new_users),
            median(new_users$LDL_pct_reduction, na.rm = TRUE) * 100))
cat(sprintf("Overall expected (weighted): %.1f%%\n",
            mean(new_users$expected_reduction, na.rm = TRUE) * 100))

# Wilcoxon test: is the observed vs expected significantly different?
wt <- wilcox.test(new_users$LDL_pct_reduction, new_users$expected_reduction,
                  paired = TRUE)
cat(sprintf("Paired Wilcoxon p-value: %.2e\n", wt$p.value))
cat("(p > 0.05 suggests our factors are reasonable)\n\n")

# ==============================================================================
# 6. IMPACT ON TUDOR CLASSIFICATION
# ==============================================================================
cat("==========================================================\n")
cat(" IMPACT OF STATIN CORRECTION ON CLASSIFICATION\n")
cat("==========================================================\n")

# For new statin users: compare TUDOR using raw V1 LDL vs corrected V1 LDL
TUDOR_WEIGHTS <- list(
  intercept = 0.755722,
  beta_LDL  = 0.057911,
  beta_Trig = 0.492412,
  beta_HDL  = -1.128045,
  beta_Age  = -0.033393,
  beta_Sex  = -0.088550
)

# Merge needed variables from main df
new_users <- merge(new_users, df[, c("eid", "HDL.1", "TRG.1", "Age_at_LDL1",
                                       "Gender_num", "is_fh_genetic")],
                   by = "eid", all.x = TRUE)

# Raw (uncorrected) TUDOR using V1 LDL
trig_raw <- new_users$LDL_v1 / (new_users$TRG.1 + 0.1)
score_raw <- TUDOR_WEIGHTS$intercept +
  TUDOR_WEIGHTS$beta_LDL * new_users$LDL_v1 +
  TUDOR_WEIGHTS$beta_Trig * trig_raw +
  TUDOR_WEIGHTS$beta_HDL * new_users$HDL.1 +
  TUDOR_WEIGHTS$beta_Age * new_users$Age_at_LDL1 +
  TUDOR_WEIGHTS$beta_Sex * new_users$Gender_num
prob_raw <- 1 / (1 + exp(-score_raw))

# Corrected TUDOR using statin-corrected V1 LDL
correction <- REDUCTION_FACTORS[new_users$statin_v1]
correction[is.na(correction)] <- 0.30
ldl_corrected <- new_users$LDL_v1 / (1 - correction)

trig_corr <- ldl_corrected / (new_users$TRG.1 + 0.1)
score_corr <- TUDOR_WEIGHTS$intercept +
  TUDOR_WEIGHTS$beta_LDL * ldl_corrected +
  TUDOR_WEIGHTS$beta_Trig * trig_corr +
  TUDOR_WEIGHTS$beta_HDL * new_users$HDL.1 +
  TUDOR_WEIGHTS$beta_Age * new_users$Age_at_LDL1 +
  TUDOR_WEIGHTS$beta_Sex * new_users$Gender_num
prob_corr <- 1 / (1 + exp(-score_corr))

cat(sprintf("Mean TUDOR prob (raw V1 LDL):       %.4f\n", mean(prob_raw, na.rm = TRUE)))
cat(sprintf("Mean TUDOR prob (corrected V1 LDL): %.4f\n", mean(prob_corr, na.rm = TRUE)))
cat(sprintf("Mean absolute change:               %+.4f\n",
            mean(prob_corr - prob_raw, na.rm = TRUE)))

# Reclassification at 50th percentile threshold
thresh <- 0.5
cat(sprintf("\nAt threshold = %.2f:\n", thresh))
cat(sprintf("  Flagged (raw):       %d (%.1f%%)\n",
            sum(prob_raw > thresh, na.rm = TRUE),
            mean(prob_raw > thresh, na.rm = TRUE) * 100))
cat(sprintf("  Flagged (corrected): %d (%.1f%%)\n",
            sum(prob_corr > thresh, na.rm = TRUE),
            mean(prob_corr > thresh, na.rm = TRUE) * 100))

# ==============================================================================
# 7. SAVE RESULTS
# ==============================================================================
longitudinal_results <- list(
  status = "complete",
  n_new_users = nrow(new_users),
  statin_summary = statin_summary,
  overall_observed_median = median(new_users$LDL_pct_reduction, na.rm = TRUE),
  overall_expected_mean = mean(new_users$expected_reduction, na.rm = TRUE),
  wilcoxon_p = wt$p.value,
  classification_impact = list(
    mean_prob_raw = mean(prob_raw, na.rm = TRUE),
    mean_prob_corrected = mean(prob_corr, na.rm = TRUE),
    mean_absolute_change = mean(prob_corr - prob_raw, na.rm = TRUE)
  ),
  visit0_statin_free = sum(!meds0$on_statin_v0),
  visit0_on_statin = sum(meds0$on_statin_v0)
)

saveRDS(longitudinal_results, file.path(OUTPUT_DIR, "03_longitudinal_results.rds"))

cat("\n=== 03_longitudinal_validation.R COMPLETE ===\n")
