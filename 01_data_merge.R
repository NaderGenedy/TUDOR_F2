# ==============================================================================
# TUDOR PIPELINE: STEP 01 — DATA MERGE & FEATURE ENGINEERING
# ==============================================================================
# PURPOSE: Merge all UKB extractions, fix medication codes, calculate
#          TUDOR v2 scores and eDLCN comparator.
#
# INPUT:   TUDOR_UKB_Features.csv (main dataset from prior work)
#          + 5 optional CSVs from 00_extract_ukbrap.sh
#
# OUTPUT:  tudor_analysis_ready.rds (single analysis-ready file)
#
# CRITICAL FIX: Code 1141146234 = ATORVASTATIN (confirmed 18,028 reports),
#               NOT Rosuvastatin as previously labeled.
# ==============================================================================

set.seed(42)  # Reproducibility

library(data.table)
library(dplyr)

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

cat("=== TUDOR PIPELINE: 01_data_merge.R ===\n")
cat("Data directory:", DATA_DIR, "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# ==============================================================================
# 1. LOAD MAIN DATASET
# ==============================================================================
cat("Loading main features file...\n")
df <- fread(file.path(DATA_DIR, "TUDOR_UKB_Features.csv"))
cat("  Loaded:", nrow(df), "participants\n")
cat("  Columns:", ncol(df), "\n")
cat("  Genetic FH cases:", sum(df$is_fh_genetic, na.rm = TRUE), "\n\n")

# ==============================================================================
# 2. CORRECTED STATIN CODE MAPPING
# ==============================================================================
# CRITICAL: These codes are from UKB Data-Coding 4
# Verified by cross-checking expected UK prescribing frequencies:
#   Simvastatin ~49k (most prescribed), Atorvastatin ~13-16k, Rosuvastatin ~2.6k
#
# Code 1141146234 = Atorvastatin  (CONFIRMED: ~13-16k reports)
# Code 1140861958 = Simvastatin   (CONFIRMED: ~49k reports — UK's #1 statin)
# Code 1141146138 = Pravastatin   (expected: ~2k reports)
# Code 1141192414 = Fluvastatin   (expected: ~200 reports)
# Code 1141192410 = Rosuvastatin  (expected: ~2.6k reports)
#
# BUGFIX v2: Previous version had Simvastatin=1140888594 and Pravastatin=1140861958
#            which SWAPPED ~49,000 patients (Simvastatin→201, Pravastatin→48,901).
#            Now corrected: 1140861958=Simvastatin, 1141146138=Pravastatin.

STATIN_CODES <- list(
  "Atorvastatin"  = c(1141146234),
  "Simvastatin"   = c(1140861958, 1140888594),   # Primary + secondary code
  "Pravastatin"   = c(1141146138),
  "Fluvastatin"   = c(1141192414, 1140861970),    # Primary + secondary code
  "Rosuvastatin"  = c(1141192410)
)

ALL_STATIN_CODES <- unlist(STATIN_CODES)

# Canonical function: assign statin name from medication code columns
assign_statin_name <- function(dt, med_col_pattern = "p20003_i0") {
  med_cols <- grep(med_col_pattern, names(dt), value = TRUE)
  if (length(med_cols) == 0) {
    cat("  WARNING: No medication columns matching '", med_col_pattern, "' found\n")
    return(rep("None", nrow(dt)))
  }

  statin_names <- character(nrow(dt))
  for (i in seq_len(nrow(dt))) {
    patient_meds <- as.numeric(dt[i, ..med_cols])
    patient_meds <- patient_meds[!is.na(patient_meds)]

    found <- "None"
    for (drug_name in names(STATIN_CODES)) {
      if (any(STATIN_CODES[[drug_name]] %in% patient_meds)) {
        found <- drug_name
        break
      }
    }
    # Check for other statins not in our list
    if (found == "None" && any(patient_meds %in% ALL_STATIN_CODES)) {
      found <- "Other_Statin"
    }
    statin_names[i] <- found
  }
  return(statin_names)
}

# Re-assign statin names using corrected codes
# (Only if we have the raw medication columns; otherwise keep existing)
med_cols_i0 <- grep("p20003_i0", names(df), value = TRUE)
if (length(med_cols_i0) > 0) {
  cat("Re-assigning statin names with corrected codes...\n")
  df$statin_name_corrected <- assign_statin_name(df, "p20003_i0")

  # Report changes
  old_counts <- table(df$statin_name)
  new_counts <- table(df$statin_name_corrected)
  cat("  OLD statin counts:\n")
  print(old_counts)
  cat("\n  NEW statin counts (corrected codes):\n")
  print(new_counts)

  df$statin_name_original <- df$statin_name
  df$statin_name <- df$statin_name_corrected
  cat("\n")
} else {
  cat("NOTE: No raw medication columns found. Using existing statin_name column.\n")
  cat("  WARNING: Original statin_name may have incorrect Atorvastatin/Rosuvastatin mapping.\n\n")
}

# ==============================================================================
# 3. REAL-WORLD STATIN REDUCTION FACTORS
# ==============================================================================
# Derived from longitudinal validation (03_longitudinal_validation.R)
# These are UK population-level observed reductions, NOT clinical trial values.

REDUCTION_FACTORS <- list(
  "Atorvastatin"  = 0.38,  # Trial: 0.45-0.50; Observed ~38%
  "Simvastatin"   = 0.35,  # Trial: 0.37; Observed ~35%
  "Rosuvastatin"  = 0.34,  # Trial: 0.50-0.55; Observed ~34% (adherence gap)
  "Pravastatin"   = 0.25,  # Trial: 0.28; Observed ~25%
  "Fluvastatin"   = 0.22,  # Trial: 0.26; Observed ~22%
  "Lovastatin"    = 0.25,  # Approximate
  "Other_Statin"  = 0.30,  # Conservative default
  "None"          = 0.00
)

# Sanity check: all factors between 0 and 1
for (nm in names(REDUCTION_FACTORS)) {
  rf <- REDUCTION_FACTORS[[nm]]
  stopifnot(rf >= 0 && rf < 1)
}

# Calculate Real-World calibrated untreated LDL
calc_untreated_ldl <- function(ldl_treated, statin_name) {
  factor <- sapply(tolower(statin_name), function(sn) {
    # Case-insensitive matching
    matched <- FALSE
    for (drug in names(REDUCTION_FACTORS)) {
      if (tolower(drug) == sn) {
        return(REDUCTION_FACTORS[[drug]])
      }
    }
    return(0)  # Default: no correction
  })
  factor <- as.numeric(factor)
  # Avoid division by zero
  denom <- 1 - factor
  denom[denom == 0] <- 1
  return(ldl_treated / denom)
}

cat("Calculating Real-World calibrated LDL...\n")
df$LDL_RW <- calc_untreated_ldl(df$LDL_treated, df$statin_name)
df$Trig_Filter_RW <- df$LDL_RW / (df$TRG.1 + 0.1)

cat("  LDL_treated range:", round(range(df$LDL_treated, na.rm = TRUE), 2), "\n")
cat("  LDL_RW range:", round(range(df$LDL_RW, na.rm = TRUE), 2), "\n\n")

# ==============================================================================
# 4. MERGE ADDITIONAL EXTRACTIONS (if available)
# ==============================================================================

# --- 4a. Lp(a) ---
lpa_file <- file.path(DATA_DIR, "ukb_lpa.csv")
if (file.exists(lpa_file)) {
  cat("Merging Lp(a) data...\n")
  lpa <- fread(lpa_file)
  names(lpa)[names(lpa) == "participant.p30790_i0"] <- "Lpa_nmol"
  df <- merge(df, lpa[, .(participant.eid, Lpa_nmol)],
              by = "participant.eid", all.x = TRUE)
  cat("  Lp(a) available for:", sum(!is.na(df$Lpa_nmol)), "participants\n")
} else {
  cat("NOTE: ukb_lpa.csv not found. Skipping Lp(a).\n")
  df$Lpa_nmol <- NA_real_
}

# --- 4b. ASCVD History ---
# Field 6150: Vascular/heart problems diagnosed by doctor
# UKB-RAP format: bracket-encoded values e.g. "[-7]", "[1]", "[1,3]"
# Codes: 1=Heart attack, 2=Angina, 3=Stroke, 4=High BP, -7=None, -3=Prefer not to say
ascvd_file <- file.path(DATA_DIR, "ukb_ascvd.csv")
if (file.exists(ascvd_file)) {
  cat("Merging ASCVD history...\n")
  ascvd <- fread(ascvd_file)
  ascvd_cols <- grep("p6150", names(ascvd), value = TRUE)

  # Parser for bracket-encoded values: "[1,3]" -> c(1, 3)
  parse_bracket <- function(x) {
    if (is.na(x) || x == "" || x == "[]") return(integer(0))
    vals <- as.integer(strsplit(gsub("\\[|\\]", "", as.character(x)), ",")[[1]])
    vals[!is.na(vals)]
  }

  # Check each instance column for MI (1), Angina (2), Stroke (3)
  ascvd$has_mi <- FALSE
  ascvd$has_angina <- FALSE
  ascvd$has_stroke <- FALSE

  for (col in ascvd_cols) {
    parsed <- lapply(ascvd[[col]], parse_bracket)
    ascvd$has_mi     <- ascvd$has_mi     | sapply(parsed, function(v) 1L %in% v)
    ascvd$has_angina  <- ascvd$has_angina  | sapply(parsed, function(v) 2L %in% v)
    ascvd$has_stroke  <- ascvd$has_stroke  | sapply(parsed, function(v) 3L %in% v)
  }

  ascvd$has_any_cvd <- ascvd$has_mi | ascvd$has_angina | ascvd$has_stroke

  df <- merge(df, ascvd[, .(participant.eid, has_mi, has_angina, has_stroke, has_any_cvd)],
              by = "participant.eid", all.x = TRUE)

  # Fill NAs with FALSE (no event)
  for (col in c("has_mi", "has_angina", "has_stroke", "has_any_cvd")) {
    df[[col]][is.na(df[[col]])] <- FALSE
  }
  cat("  ASCVD events found:", sum(df$has_any_cvd), "\n")
  cat("    MI:", sum(df$has_mi), "| Angina:", sum(df$has_angina),
      "| Stroke:", sum(df$has_stroke), "\n")
} else {
  cat("NOTE: ukb_ascvd.csv not found. Setting ASCVD to FALSE.\n")
  df$has_mi <- FALSE
  df$has_angina <- FALSE
  df$has_stroke <- FALSE
  df$has_any_cvd <- FALSE
}

# --- 4c. MI/Angina Age ---
cvd_age_file <- file.path(DATA_DIR, "ukb_cvd_age.csv")
if (file.exists(cvd_age_file)) {
  cat("Merging CVD age data...\n")
  cvd_age <- fread(cvd_age_file)
  names(cvd_age)[names(cvd_age) == "participant.p3894_i0"] <- "MI_age"
  names(cvd_age)[names(cvd_age) == "participant.p3627_i0"] <- "Angina_age"
  df <- merge(df, cvd_age[, .(participant.eid, MI_age, Angina_age)],
              by = "participant.eid", all.x = TRUE)
} else {
  cat("NOTE: ukb_cvd_age.csv not found.\n")
  df$MI_age <- NA_real_
  df$Angina_age <- NA_real_
}

# Premature ASCVD: event < 55 (men) or < 60 (women)
mi_age_safe <- ifelse(is.na(df$MI_age), 999, df$MI_age)
angina_age_safe <- ifelse(is.na(df$Angina_age), 999, df$Angina_age)
min_cvd_age <- pmin(mi_age_safe, angina_age_safe)

df$Premature_ASCVD <- ifelse(
  df$has_any_cvd &
    ((df$Gender_num == 1 & min_cvd_age < 55) |
     (df$Gender_num == 0 & min_cvd_age < 60)),
  1, 0
)
cat("  Premature ASCVD cases:", sum(df$Premature_ASCVD), "\n\n")

# ==============================================================================
# 5. TUDOR v2 SCORE (External Validation — Wales Weights)
# ==============================================================================
# These weights are from the Wales training set (Elastic Net, lambda.min)
# They are FIXED and must NOT be re-estimated on UKB data.

TUDOR_WEIGHTS <- list(
  intercept = 0.755722,
  beta_LDL  = 0.057911,
  beta_Trig = 0.492412,
  beta_HDL  = -1.128045,
  beta_Age  = -0.033393,
  beta_Sex  = -0.088550
)

cat("Calculating TUDOR v2 score (Wales weights)...\n")
df$tudor_score <- with(TUDOR_WEIGHTS,
  intercept +
  (beta_LDL  * df$LDL_RW) +
  (beta_Trig * df$Trig_Filter_RW) +
  (beta_HDL  * df$HDL.1) +
  (beta_Age  * df$Age_at_LDL1) +
  (beta_Sex  * df$Gender_num)
)
df$tudor_prob <- 1 / (1 + exp(-df$tudor_score))

cat("  Score range:", round(range(df$tudor_score, na.rm = TRUE), 3), "\n")
cat("  Prob range:", round(range(df$tudor_prob, na.rm = TRUE), 4), "\n\n")

# ==============================================================================
# 6. eDLCN SCORE (Electronic Dutch Lipid Clinic Network)
# ==============================================================================
# Standard eDLCN scoring:
#   LDL >= 8.5 -> 8pts; >= 6.5 -> 5pts; >= 5.0 -> 3pts; >= 4.0 -> 1pt
#   Premature ASCVD -> 2pts
#   (Family history and physical signs not available electronically)

cat("Calculating eDLCN score...\n")
df$edlcn_ldl_pts <- ifelse(df$LDL_RW >= 8.5, 8,
                    ifelse(df$LDL_RW >= 6.5, 5,
                    ifelse(df$LDL_RW >= 5.0, 3,
                    ifelse(df$LDL_RW >= 4.0, 1, 0))))

df$edlcn_ascvd_pts <- df$Premature_ASCVD * 2

df$edlcn_score <- df$edlcn_ldl_pts + df$edlcn_ascvd_pts

cat("  eDLCN distribution:\n")
print(table(df$edlcn_score))
cat("\n")

# ==============================================================================
# 7. ADDITIONAL DERIVED VARIABLES
# ==============================================================================

# TC/LDL Ratio (LDL Purity)
df$TC_LDL_Ratio <- df$CHOL / df$LDL_RW

# LDL Purity (inverse: LDL/TC, higher = purer hypercholesterolaemia)
df$LDL_Purity <- df$LDL_RW / df$CHOL

# ApoB (if present)
if ("participant.p30640_i0" %in% names(df)) {
  df$ApoB <- df$participant.p30640_i0
} else if (!"ApoB" %in% names(df)) {
  df$ApoB <- NA_real_
}

# ApoB/LDL Ratio (particle density proxy)
df$ApoB_LDL_Ratio <- df$ApoB / df$LDL_RW

# Friedewald LDL (for sensitivity analysis comparison)
# Friedewald: LDL = TC - HDL - (TG/2.2) [mmol/L]
df$LDL_Friedewald <- df$CHOL - df$HDL.1 - (df$TRG.1 / 2.2)

# BMI (impute median if missing)
if (!"BMI_imputed" %in% names(df)) {
  if ("participant.p21001_i0" %in% names(df)) {
    df$BMI_imputed <- df$participant.p21001_i0
    df$BMI_imputed[is.na(df$BMI_imputed)] <- median(df$BMI_imputed, na.rm = TRUE)
  } else {
    df$BMI_imputed <- 26.5  # UK population median
  }
}

# Statin tier (intensity classification)
df$statin_tier <- ifelse(tolower(df$statin_name) == "rosuvastatin", 3,
                  ifelse(tolower(df$statin_name) %in% c("atorvastatin", "simvastatin"), 2,
                  ifelse(tolower(df$statin_name) == "none", 0, 1)))

# ==============================================================================
# 8. DEFINE ANALYSIS COHORTS
# ==============================================================================
df$cohort_high_risk <- df$LDL_RW > 4.9
df$cohort_moderate  <- df$LDL_RW >= 2.6 & df$LDL_RW <= 4.9
df$cohort_low_risk  <- df$LDL_RW < 2.6

cat("Cohort sizes:\n")
cat("  High Risk (LDL > 4.9):", sum(df$cohort_high_risk, na.rm = TRUE), "\n")
cat("  Moderate (2.6-4.9):", sum(df$cohort_moderate, na.rm = TRUE), "\n")
cat("  Low Risk (<2.6):", sum(df$cohort_low_risk, na.rm = TRUE), "\n\n")

# ==============================================================================
# 9. SAVE
# ==============================================================================
out_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
saveRDS(df, out_file)
cat("Saved analysis-ready dataset to:", out_file, "\n")
cat("  Rows:", nrow(df), "\n")
cat("  Key columns: tudor_score, tudor_prob, edlcn_score, LDL_RW,",
    "Trig_Filter_RW, Premature_ASCVD, Lpa_nmol, ApoB\n")

# Summary statistics
cat("\n=== SUMMARY ===\n")
cat("Total N:", nrow(df), "\n")
cat("Genetic FH cases:", sum(df$is_fh_genetic), "\n")
cat("On statin:", sum(df$statin_name != "None"), "\n")
cat("Statin distribution:\n")
print(table(df$statin_name))

cat("\n=== 01_data_merge.R COMPLETE ===\n")
