# ==============================================================================
# TUDOR PIPELINE: STEP 16 — SENSITIVITY ANALYSIS: EXCLUDING SECONDARY CAUSES
#                             OF HYPERCHOLESTEROLAEMIA
# ==============================================================================
# PURPOSE: Remove participants with secondary causes of elevated LDL-C to
#          ensure TUDOR's external validation reflects discrimination of
#          true monogenic FH vs polygenic/lifestyle hypercholesterolaemia,
#          not metabolic phenocopies.
#
# SECONDARY CAUSES EXCLUDED:
#   1. Hypothyroidism        — reduced LDL receptor expression
#   2. Type 2 Diabetes       — insulin resistance / dyslipidaemia
#   3. Nephrotic syndrome    — hepatic lipoprotein overproduction
#   4. Obstructive jaundice  — impaired bile acid excretion
#
# IDENTIFICATION:
#   - Field 20002 (self-reported non-cancer illness, UKB Data-Coding 6)
#   - Field 2443  (doctor-diagnosed diabetes, binary)
#   - Field 30750 (HbA1c, mmol/mol — WHO threshold >= 48)
#   - ICD-10 codes (Field 41270) if available
#
# INPUT:   tudor_analysis_ready.rds (from 01_data_merge.R)
#          ukb_selfreport_illness.csv (from 00d_extract_ukbrap_secondary.sh)
#          ukb_dm_hba1c.csv (from 00b, already used in pipeline)
#          ukb_icd10_secondary.csv (OPTIONAL — if ICD-10 dispensed)
#
# OUTPUT:  Comparative validation: original vs cleaned cohort
#          DeLong tests, exclusion counts, updated figures and tables
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
TABLE_DIR  <- file.path(OUTPUT_DIR, "tables")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("TUDOR PIPELINE: 16 — SECONDARY CAUSES EXCLUSION SENSITIVITY ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD MAIN DATASET
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
if (!file.exists(rds_file)) stop("Run 01_data_merge.R first!")
df <- readRDS(rds_file)

if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
eid_col <- ifelse("eid" %in% names(df), "eid", "participant.eid")

cat("Loaded:", nrow(df), "participants\n")
cat("Genetic FH cases:", sum(df$is_fh_genetic), "\n\n")

# ==============================================================================
# 2. UKB DATA-CODING 6: SELF-REPORTED NON-CANCER ILLNESS CODES
# ==============================================================================
# Reference: https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=6
#
# These codes identify secondary causes of hypercholesterolaemia.
# Each can independently elevate LDL-C and mimic an FH phenotype.

SECONDARY_CODES <- list(
  hypothyroidism = list(
    codes = c(1226),             # 1226 = Hypothyroidism/myxoedema
    label = "Hypothyroidism",
    rationale = "Reduced LDL receptor expression -> elevated LDL-C"
  ),
  t2dm = list(
    codes = c(1220, 1223),       # 1220 = Diabetes, 1223 = Type 2 diabetes
    label = "Type 2 Diabetes",
    rationale = "Insulin resistance -> dyslipidaemia, elevated LDL-C"
  ),
  nephrotic = list(
    codes = c(1519),             # 1519 = Nephrotic syndrome
    label = "Nephrotic Syndrome",
    rationale = "Hepatic lipoprotein overproduction -> elevated LDL-C"
  ),
  obstructive_jaundice = list(
    codes = c(1506, 1158, 1604), # 1506 = Primary biliary cirrhosis
    label = "Obstructive Jaundice / Cholestasis",  # 1158 = Cholangitis
    rationale = "Impaired bile acid excretion -> cholesterol accumulation"
    # 1604 = Liver failure/cirrhosis
  )
)

cat("Secondary cause codes (UKB Data-Coding 6):\n")
for (cond in names(SECONDARY_CODES)) {
  info <- SECONDARY_CODES[[cond]]
  cat(sprintf("  %-35s codes: %s\n",
              info$label, paste(info$codes, collapse = ", ")))
}
cat("\n")

# ==============================================================================
# 3. LOAD SELF-REPORTED ILLNESS DATA (Field 20002)
# ==============================================================================
# Strategy: Check if p20002 columns already exist in main df.
#           If not, load from separate extraction CSV.

p20002_cols <- grep("p20002", names(df), value = TRUE)

if (length(p20002_cols) >= 5) {
  cat("Field 20002 columns found in main dataset:", length(p20002_cols), "columns\n")
} else {
  # Load from separate extraction
  illness_file <- file.path(DATA_DIR, "ukb_selfreport_illness.csv")
  if (!file.exists(illness_file)) {
    stop(paste0(
      "Field 20002 not found in main data or as separate CSV.\n",
      "Run 00d_extract_ukbrap_secondary.sh on UKB-RAP first.\n",
      "Expected file: ", illness_file
    ))
  }

  cat("Loading self-reported illness data from:", illness_file, "\n")
  illness <- fread(illness_file)
  cat("  Loaded:", nrow(illness), "participants,", ncol(illness), "columns\n")

  # Merge into main dataset
  df <- merge(df, illness, by.x = eid_col, by.y = "participant.eid", all.x = TRUE)
  p20002_cols <- grep("p20002", names(df), value = TRUE)
  cat("  After merge:", length(p20002_cols), "illness columns available\n")
}
cat("\n")

# ==============================================================================
# 4. FLAG SECONDARY CAUSES FROM FIELD 20002
# ==============================================================================
cat("Scanning Field 20002 for secondary cause codes...\n\n")

# Helper: check if ANY of the p20002 columns contain ANY of the target codes
flag_illness <- function(dt, target_codes, col_pattern = "p20002") {
  illness_cols <- grep(col_pattern, names(dt), value = TRUE)
  if (length(illness_cols) == 0) return(rep(FALSE, nrow(dt)))

  flag <- rep(FALSE, nrow(dt))
  for (col in illness_cols) {
    vals <- as.numeric(dt[[col]])
    flag <- flag | (!is.na(vals) & vals %in% target_codes)
  }
  return(flag)
}

# Flag each condition
df$has_hypothyroidism <- flag_illness(df, SECONDARY_CODES$hypothyroidism$codes)
df$has_t2dm_selfreport <- flag_illness(df, SECONDARY_CODES$t2dm$codes)
df$has_nephrotic <- flag_illness(df, SECONDARY_CODES$nephrotic$codes)
df$has_obstructive_jaundice <- flag_illness(df, SECONDARY_CODES$obstructive_jaundice$codes)

# ==============================================================================
# 5. AUGMENT T2DM FLAG WITH FIELD 2443 AND HbA1c
# ==============================================================================
# Field 20002 alone may miss diabetes cases. Augment with:
#   - Field 2443 (doctor-diagnosed diabetes) = 1
#   - HbA1c >= 48 mmol/mol (WHO diagnostic threshold)

# Load DM data if not already merged
if (!"dm_diagnosed" %in% names(df)) {
  dm_file <- file.path(DATA_DIR, "ukb_dm_hba1c.csv")
  if (file.exists(dm_file)) {
    cat("Loading diabetes/HbA1c data...\n")
    dm <- fread(dm_file)
    names(dm)[names(dm) == "participant.p2443_i0"] <- "dm_diagnosed"
    names(dm)[names(dm) == "participant.p30750_i0"] <- "HbA1c"
    df <- merge(df, dm[, .(participant.eid, dm_diagnosed, HbA1c)],
                by.x = eid_col, by.y = "participant.eid", all.x = TRUE)
  }
}

# Also check for p2443 column names from prior merge
if (!"dm_diagnosed" %in% names(df)) {
  p2443_col <- grep("p2443", names(df), value = TRUE)
  if (length(p2443_col) > 0) {
    df$dm_diagnosed <- as.numeric(df[[p2443_col[1]]])
  } else {
    df$dm_diagnosed <- NA_real_
  }
}
if (!"HbA1c" %in% names(df)) {
  hba1c_col <- grep("p30750", names(df), value = TRUE)
  if (length(hba1c_col) > 0) {
    df$HbA1c <- as.numeric(df[[hba1c_col[1]]])
  } else {
    df$HbA1c <- NA_real_
  }
}

# Composite T2DM flag
df$has_t2dm <- df$has_t2dm_selfreport |
  (!is.na(df$dm_diagnosed) & df$dm_diagnosed == 1) |
  (!is.na(df$HbA1c) & df$HbA1c >= 48)

cat(sprintf("T2DM identification:\n"))
cat(sprintf("  Self-reported (Field 20002):    %d\n", sum(df$has_t2dm_selfreport, na.rm = TRUE)))
cat(sprintf("  Doctor-diagnosed (Field 2443):  %d\n", sum(!is.na(df$dm_diagnosed) & df$dm_diagnosed == 1)))
cat(sprintf("  HbA1c >= 48 mmol/mol:           %d\n", sum(!is.na(df$HbA1c) & df$HbA1c >= 48)))
cat(sprintf("  Composite T2DM flag:            %d\n\n", sum(df$has_t2dm, na.rm = TRUE)))

# ==============================================================================
# 6. LOAD ICD-10 DATA IF AVAILABLE (OPTIONAL AUGMENTATION)
# ==============================================================================
icd_file <- file.path(DATA_DIR, "ukb_icd10_secondary.csv")
if (file.exists(icd_file)) {
  cat("Loading ICD-10 hospital episode data...\n")
  icd <- fread(icd_file)
  icd_cols <- grep("p41270", names(icd), value = TRUE)

  if (length(icd_cols) > 0) {
    # ICD-10 codes for secondary causes:
    ICD_HYPOTHYROID <- c("E03", "E030", "E031", "E032", "E033", "E034", "E035",
                          "E038", "E039", "E890")
    ICD_T2DM        <- c("E11", "E110", "E111", "E112", "E113", "E114", "E115",
                          "E116", "E117", "E118", "E119", "E14")
    ICD_NEPHROTIC   <- c("N04", "N040", "N041", "N042", "N043", "N044", "N045",
                          "N046", "N047", "N048", "N049")
    ICD_JAUNDICE    <- c("K83", "K830", "K831", "K74", "K743", "K744", "K745",
                          "K710", "R17")  # Cholangitis, cirrhosis, jaundice

    # Helper: flag ICD codes across all columns
    flag_icd <- function(dt, icd_cols, target_prefixes) {
      flag <- rep(FALSE, nrow(dt))
      for (col in icd_cols) {
        vals <- as.character(dt[[col]])
        for (prefix in target_prefixes) {
          flag <- flag | (!is.na(vals) & startsWith(vals, prefix))
        }
      }
      return(flag)
    }

    # Merge ICD-10 into df
    df <- merge(df, icd, by.x = eid_col, by.y = "participant.eid", all.x = TRUE)
    icd_cols_in_df <- grep("p41270", names(df), value = TRUE)

    # Augment flags with ICD-10
    icd_hypo <- flag_icd(df, icd_cols_in_df, ICD_HYPOTHYROID)
    icd_dm   <- flag_icd(df, icd_cols_in_df, ICD_T2DM)
    icd_neph <- flag_icd(df, icd_cols_in_df, ICD_NEPHROTIC)
    icd_jaun <- flag_icd(df, icd_cols_in_df, ICD_JAUNDICE)

    cat(sprintf("  ICD-10 additional cases found:\n"))
    cat(sprintf("    Hypothyroidism: +%d\n", sum(icd_hypo & !df$has_hypothyroidism)))
    cat(sprintf("    T2DM:           +%d\n", sum(icd_dm & !df$has_t2dm)))
    cat(sprintf("    Nephrotic:      +%d\n", sum(icd_neph & !df$has_nephrotic)))
    cat(sprintf("    Jaundice:       +%d\n", sum(icd_jaun & !df$has_obstructive_jaundice)))

    df$has_hypothyroidism <- df$has_hypothyroidism | icd_hypo
    df$has_t2dm <- df$has_t2dm | icd_dm
    df$has_nephrotic <- df$has_nephrotic | icd_neph
    df$has_obstructive_jaundice <- df$has_obstructive_jaundice | icd_jaun
    cat("\n")
  }
} else {
  cat("NOTE: ICD-10 data not available. Using self-reported + DM fields only.\n\n")
}

# ==============================================================================
# 7. EXCLUSION SUMMARY
# ==============================================================================
df$has_any_secondary <- df$has_hypothyroidism | df$has_t2dm |
  df$has_nephrotic | df$has_obstructive_jaundice

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SECONDARY CAUSE EXCLUSION SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

excl_summary <- data.frame(
  Condition = c("Hypothyroidism", "Type 2 Diabetes (composite)",
                "Nephrotic Syndrome", "Obstructive Jaundice/Cholestasis",
                "ANY secondary cause"),
  Total_N = c(
    sum(df$has_hypothyroidism, na.rm = TRUE),
    sum(df$has_t2dm, na.rm = TRUE),
    sum(df$has_nephrotic, na.rm = TRUE),
    sum(df$has_obstructive_jaundice, na.rm = TRUE),
    sum(df$has_any_secondary, na.rm = TRUE)
  ),
  FH_positive = c(
    sum(df$has_hypothyroidism & df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_t2dm & df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_nephrotic & df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_obstructive_jaundice & df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_any_secondary & df$is_fh_genetic, na.rm = TRUE)
  ),
  FH_negative = c(
    sum(df$has_hypothyroidism & !df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_t2dm & !df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_nephrotic & !df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_obstructive_jaundice & !df$is_fh_genetic, na.rm = TRUE),
    sum(df$has_any_secondary & !df$is_fh_genetic, na.rm = TRUE)
  )
)
excl_summary$Pct_of_total <- sprintf("%.1f%%", 100 * excl_summary$Total_N / nrow(df))

print(excl_summary, row.names = FALSE)
cat("\n")

# Overlap matrix
cat("Overlap between secondary causes:\n")
overlap_mat <- matrix(0, nrow = 4, ncol = 4)
flags <- list(df$has_hypothyroidism, df$has_t2dm, df$has_nephrotic, df$has_obstructive_jaundice)
flag_names <- c("Hypothyroid", "T2DM", "Nephrotic", "Jaundice")
rownames(overlap_mat) <- flag_names
colnames(overlap_mat) <- flag_names
for (i in 1:4) {
  for (j in 1:4) {
    overlap_mat[i, j] <- sum(flags[[i]] & flags[[j]], na.rm = TRUE)
  }
}
print(overlap_mat)
cat("\n")

# ==============================================================================
# 8. CREATE CLEANED COHORT
# ==============================================================================
df_clean <- df[df$has_any_secondary == FALSE, ]

cat(sprintf("Original cohort:  N = %d  (FH = %d, prevalence = %.2f%%)\n",
            nrow(df), sum(df$is_fh_genetic),
            100 * mean(df$is_fh_genetic)))
cat(sprintf("Excluded:         N = %d  (FH = %d)\n",
            sum(df$has_any_secondary), sum(df$has_any_secondary & df$is_fh_genetic)))
cat(sprintf("Cleaned cohort:   N = %d  (FH = %d, prevalence = %.2f%%)\n\n",
            nrow(df_clean), sum(df_clean$is_fh_genetic),
            100 * mean(df_clean$is_fh_genetic)))

# ==============================================================================
# 9. RE-RUN TUDOR VALIDATION: ORIGINAL vs CLEANED
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPARATIVE VALIDATION: ORIGINAL vs CLEANED COHORT\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# --- Helper function for full validation on a cohort ---
run_validation <- function(data, label) {
  # Full population
  full <- data[!is.na(data$tudor_prob), ]
  # High-risk cohort (LDL > 4.9)
  hr <- full[full$LDL_RW > 4.9, ]

  # --- Lipid clinic cohort (matching script 11 criteria) ---
  lc <- full[full$LDL_RW > 4.9 |
               full$CHOL > 7.5 |
               full$Premature_ASCVD == 1, ]

  results <- list()
  for (cohort_name in c("Full population", "High-risk (LDL>4.9)", "Lipid clinic")) {
    sub <- switch(cohort_name,
                  "Full population" = full,
                  "High-risk (LDL>4.9)" = hr,
                  "Lipid clinic" = lc)

    if (sum(sub$is_fh_genetic) < 10) {
      cat(sprintf("  %-20s %-25s SKIPPED (FH < 10)\n", label, cohort_name))
      next
    }

    roc_t <- roc(sub$is_fh_genetic, sub$tudor_prob, quiet = TRUE)
    roc_e <- roc(sub$is_fh_genetic, sub$edlcn_score, quiet = TRUE)
    roc_l <- roc(sub$is_fh_genetic, sub$LDL_RW, quiet = TRUE)
    roc_f <- roc(sub$is_fh_genetic, sub$Trig_Filter_RW, quiet = TRUE)

    ci_t <- ci.auc(roc_t, method = "delong")
    ci_e <- ci.auc(roc_e, method = "delong")

    # Youden optimal
    youden_idx <- which.max(roc_t$sensitivities + roc_t$specificities - 1)
    youden_thresh <- roc_t$thresholds[youden_idx]
    youden_sens <- roc_t$sensitivities[youden_idx]
    youden_spec <- roc_t$specificities[youden_idx]

    # Brier score
    brier <- mean((sub$tudor_prob - sub$is_fh_genetic)^2)

    # NRI vs eDLCN
    tudor_correct <- (sub$tudor_prob >= youden_thresh) == sub$is_fh_genetic
    edlcn_correct <- (sub$edlcn_score >= 3) == sub$is_fh_genetic
    nri <- mean(tudor_correct) - mean(edlcn_correct)

    cat(sprintf("  %-20s %-25s TUDOR=%.3f [%.3f-%.3f]  eDLCN=%.3f  N=%d FH=%d\n",
                label, cohort_name, ci_t[2], ci_t[1], ci_t[3],
                ci_e[2], nrow(sub), sum(sub$is_fh_genetic)))

    results[[cohort_name]] <- list(
      cohort = cohort_name,
      n = nrow(sub),
      n_fh = sum(sub$is_fh_genetic),
      prevalence = mean(sub$is_fh_genetic),
      tudor_auc = as.numeric(ci_t[2]),
      tudor_ci_lo = as.numeric(ci_t[1]),
      tudor_ci_hi = as.numeric(ci_t[3]),
      edlcn_auc = as.numeric(ci_e[2]),
      edlcn_ci_lo = as.numeric(ci_e[1]),
      edlcn_ci_hi = as.numeric(ci_e[3]),
      ldl_auc = as.numeric(auc(roc_l)),
      trig_auc = as.numeric(auc(roc_f)),
      youden_thresh = youden_thresh,
      youden_sens = youden_sens,
      youden_spec = youden_spec,
      brier = brier,
      nri = nri,
      roc_tudor = roc_t,
      roc_edlcn = roc_e
    )
  }
  return(results)
}

# Run on both cohorts
cat("ORIGINAL COHORT:\n")
res_original <- run_validation(df, "Original")
cat("\nCLEANED COHORT (secondary causes excluded):\n")
res_cleaned <- run_validation(df_clean, "Cleaned")
cat("\n")

# ==============================================================================
# 10. DELONG COMPARISON: ORIGINAL vs CLEANED
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("STATISTICAL COMPARISON: CLEANED vs ORIGINAL\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

comparison_table <- data.frame()

for (cohort_name in names(res_original)) {
  if (!cohort_name %in% names(res_cleaned)) next

  orig <- res_original[[cohort_name]]
  clean <- res_cleaned[[cohort_name]]

  delta_auc <- clean$tudor_auc - orig$tudor_auc
  delta_edlcn <- clean$edlcn_auc - orig$edlcn_auc

  cat(sprintf("%-25s\n", cohort_name))
  cat(sprintf("  TUDOR AUC:  Original = %.3f -> Cleaned = %.3f  (delta = %+.3f)\n",
              orig$tudor_auc, clean$tudor_auc, delta_auc))
  cat(sprintf("  eDLCN AUC:  Original = %.3f -> Cleaned = %.3f  (delta = %+.3f)\n",
              orig$edlcn_auc, clean$edlcn_auc, delta_edlcn))
  cat(sprintf("  N:          Original = %d -> Cleaned = %d  (excluded = %d)\n",
              orig$n, clean$n, orig$n - clean$n))
  cat(sprintf("  FH:         Original = %d -> Cleaned = %d  (excluded = %d)\n",
              orig$n_fh, clean$n_fh, orig$n_fh - clean$n_fh))
  cat(sprintf("  Prevalence: Original = %.2f%% -> Cleaned = %.2f%%\n",
              100 * orig$prevalence, 100 * clean$prevalence))
  cat(sprintf("  Brier:      Original = %.4f -> Cleaned = %.4f\n",
              orig$brier, clean$brier))
  cat(sprintf("  Youden:     Sens %.1f%%/Spec %.1f%% -> Sens %.1f%%/Spec %.1f%%\n",
              100 * orig$youden_sens, 100 * orig$youden_spec,
              100 * clean$youden_sens, 100 * clean$youden_spec))
  cat("\n")

  comparison_table <- rbind(comparison_table, data.frame(
    Cohort = cohort_name,
    Original_N = orig$n,
    Cleaned_N = clean$n,
    Excluded_N = orig$n - clean$n,
    Original_FH = orig$n_fh,
    Cleaned_FH = clean$n_fh,
    TUDOR_AUC_Original = round(orig$tudor_auc, 3),
    TUDOR_AUC_Cleaned = round(clean$tudor_auc, 3),
    TUDOR_Delta = round(delta_auc, 3),
    eDLCN_AUC_Original = round(orig$edlcn_auc, 3),
    eDLCN_AUC_Cleaned = round(clean$edlcn_auc, 3),
    eDLCN_Delta = round(delta_edlcn, 3),
    Brier_Original = round(orig$brier, 4),
    Brier_Cleaned = round(clean$brier, 4)
  ))
}

# ==============================================================================
# 11. CONDITION-SPECIFIC EXCLUSION (ONE AT A TIME)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("INDIVIDUAL CONDITION EXCLUSION ANALYSIS\n")
cat("(Excluding one condition at a time to assess individual impact)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Use lipid clinic cohort as reference
lc_orig <- df[!is.na(df$tudor_prob) &
                (df$LDL_RW > 4.9 | df$CHOL > 7.5 | df$Premature_ASCVD == 1), ]

if (sum(lc_orig$is_fh_genetic) >= 10) {
  roc_orig <- roc(lc_orig$is_fh_genetic, lc_orig$tudor_prob, quiet = TRUE)

  conditions <- list(
    list(flag = "has_hypothyroidism", label = "Excl. Hypothyroidism"),
    list(flag = "has_t2dm", label = "Excl. T2DM"),
    list(flag = "has_nephrotic", label = "Excl. Nephrotic"),
    list(flag = "has_obstructive_jaundice", label = "Excl. Jaundice"),
    list(flag = "has_any_secondary", label = "Excl. ALL secondary")
  )

  cat(sprintf("%-30s %8s %6s %8s %8s %10s\n",
              "Analysis", "N", "FH", "AUC", "Delta", "DeLong p"))
  cat(strrep("-", 80), "\n")
  cat(sprintf("%-30s %8d %6d %8.3f %8s %10s\n",
              "Original (no exclusion)", nrow(lc_orig),
              sum(lc_orig$is_fh_genetic), auc(roc_orig), "-", "-"))

  individual_results <- data.frame()

  for (cond in conditions) {
    sub <- lc_orig[lc_orig[[cond$flag]] == FALSE, ]
    if (sum(sub$is_fh_genetic) < 10) {
      cat(sprintf("%-30s SKIPPED (FH < 10 after exclusion)\n", cond$label))
      next
    }

    roc_sub <- roc(sub$is_fh_genetic, sub$tudor_prob, quiet = TRUE)
    delta <- as.numeric(auc(roc_sub)) - as.numeric(auc(roc_orig))

    # DeLong test is not directly comparable (different sample sizes)
    # Report descriptively
    ci_sub <- ci.auc(roc_sub, method = "delong")

    cat(sprintf("%-30s %8d %6d %8.3f %+8.3f %10s\n",
                cond$label, nrow(sub), sum(sub$is_fh_genetic),
                auc(roc_sub), delta,
                sprintf("[%.3f-%.3f]", ci_sub[1], ci_sub[3])))

    individual_results <- rbind(individual_results, data.frame(
      Analysis = cond$label,
      N = nrow(sub),
      N_FH = sum(sub$is_fh_genetic),
      N_excluded = nrow(lc_orig) - nrow(sub),
      AUC = round(as.numeric(auc(roc_sub)), 3),
      CI_lower = round(as.numeric(ci_sub[1]), 3),
      CI_upper = round(as.numeric(ci_sub[3]), 3),
      Delta_AUC = round(delta, 3)
    ))
  }
  cat("\n")
}

# ==============================================================================
# 12. ROC COMPARISON FIGURE
# ==============================================================================
cat("Generating comparison ROC figure...\n")

# Use lipid clinic cohort for the figure
lc_clean <- df_clean[!is.na(df_clean$tudor_prob) &
                       (df_clean$LDL_RW > 4.9 | df_clean$CHOL > 7.5 |
                          df_clean$Premature_ASCVD == 1), ]

if (sum(lc_clean$is_fh_genetic) >= 10 && sum(lc_orig$is_fh_genetic) >= 10) {
  roc_clean <- roc(lc_clean$is_fh_genetic, lc_clean$tudor_prob, quiet = TRUE)
  roc_orig_plot <- roc(lc_orig$is_fh_genetic, lc_orig$tudor_prob, quiet = TRUE)

  roc_plot_data <- rbind(
    data.frame(
      sens = roc_orig_plot$sensitivities,
      fpr = 1 - roc_orig_plot$specificities,
      Model = sprintf("Original (AUC = %.3f, N = %s)",
                       auc(roc_orig_plot), format(nrow(lc_orig), big.mark = ","))
    ),
    data.frame(
      sens = roc_clean$sensitivities,
      fpr = 1 - roc_clean$specificities,
      Model = sprintf("Secondary excluded (AUC = %.3f, N = %s)",
                       auc(roc_clean), format(nrow(lc_clean), big.mark = ","))
    )
  )

  p_compare <- ggplot(roc_plot_data, aes(x = fpr, y = sens, color = Model)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("steelblue", "red3")) +
    labs(
      x = "1 - Specificity (False Positive Rate)",
      y = "Sensitivity (True Positive Rate)",
      title = "TUDOR Validation: Impact of Excluding Secondary Hypercholesterolaemia",
      subtitle = paste0("UK Biobank lipid clinic cohort | Excluded: hypothyroidism, T2DM, ",
                         "nephrotic syndrome, obstructive jaundice")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = c(0.65, 0.2),
      legend.background = element_rect(fill = "white", colour = "grey80")
    )

  ggsave(file.path(PLOT_DIR, "fig_roc_secondary_exclusion.pdf"),
         p_compare, width = 9, height = 7)
  ggsave(file.path(PLOT_DIR, "fig_roc_secondary_exclusion.png"),
         p_compare, width = 9, height = 7, dpi = 300)
  cat("  Saved: fig_roc_secondary_exclusion.pdf/png\n")
}

# ==============================================================================
# 13. SAVE TABLES
# ==============================================================================
cat("Saving tables...\n")

# Exclusion summary table
write.csv(excl_summary,
          file.path(TABLE_DIR, "secondary_exclusion_summary.csv"),
          row.names = FALSE)

# Comparison table
write.csv(comparison_table,
          file.path(TABLE_DIR, "secondary_exclusion_comparison.csv"),
          row.names = FALSE)

# Individual condition impact
if (exists("individual_results") && nrow(individual_results) > 0) {
  write.csv(individual_results,
            file.path(TABLE_DIR, "secondary_exclusion_individual.csv"),
            row.names = FALSE)
}

# Save cleaned dataset for downstream use
clean_rds <- file.path(OUTPUT_DIR, "tudor_analysis_clean_no_secondary.rds")
saveRDS(df_clean, clean_rds)
cat("  Saved cleaned dataset to:", clean_rds, "\n")

# ==============================================================================
# 14. MANUSCRIPT-READY SUMMARY
# ==============================================================================
cat("\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MANUSCRIPT-READY SUMMARY\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("SENSITIVITY ANALYSIS: EXCLUSION OF SECONDARY HYPERCHOLESTEROLAEMIA\n\n")
cat(sprintf("To assess whether secondary causes of hypercholesterolaemia\n"))
cat(sprintf("(hypothyroidism, type 2 diabetes, nephrotic syndrome, obstructive\n"))
cat(sprintf("jaundice) confounded the external validation, we repeated the\n"))
cat(sprintf("analysis after excluding %d participants (%.1f%%) with any of\n",
            sum(df$has_any_secondary),
            100 * sum(df$has_any_secondary) / nrow(df)))
cat(sprintf("these conditions.\n\n"))

cat(sprintf("Condition-specific exclusions:\n"))
cat(sprintf("  Hypothyroidism:              %6d  (%.1f%%)\n",
            sum(df$has_hypothyroidism), 100 * mean(df$has_hypothyroidism)))
cat(sprintf("  Type 2 Diabetes:             %6d  (%.1f%%)\n",
            sum(df$has_t2dm), 100 * mean(df$has_t2dm)))
cat(sprintf("  Nephrotic syndrome:          %6d  (%.1f%%)\n",
            sum(df$has_nephrotic), 100 * mean(df$has_nephrotic)))
cat(sprintf("  Obstructive jaundice:        %6d  (%.1f%%)\n",
            sum(df$has_obstructive_jaundice), 100 * mean(df$has_obstructive_jaundice)))

if (length(res_cleaned) > 0 && "Lipid clinic" %in% names(res_cleaned)) {
  orig_lc <- res_original[["Lipid clinic"]]
  clean_lc <- res_cleaned[["Lipid clinic"]]
  cat(sprintf("\nLipid clinic cohort performance:\n"))
  cat(sprintf("  Original:  TUDOR AUC = %.3f [%.3f-%.3f]  N = %d  FH = %d\n",
              orig_lc$tudor_auc, orig_lc$tudor_ci_lo, orig_lc$tudor_ci_hi,
              orig_lc$n, orig_lc$n_fh))
  cat(sprintf("  Cleaned:   TUDOR AUC = %.3f [%.3f-%.3f]  N = %d  FH = %d\n",
              clean_lc$tudor_auc, clean_lc$tudor_ci_lo, clean_lc$tudor_ci_hi,
              clean_lc$n, clean_lc$n_fh))
  cat(sprintf("  Delta AUC: %+.3f\n", clean_lc$tudor_auc - orig_lc$tudor_auc))
}

cat("\n=== 16_sensitivity_secondary_exclusion.R COMPLETE ===\n")
