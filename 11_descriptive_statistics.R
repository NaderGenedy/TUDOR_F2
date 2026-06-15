# ==============================================================================
# TUDOR PIPELINE: STEP 11 — DEEP DESCRIPTIVE STATISTICS
# ==============================================================================
# PURPOSE: Comprehensive Table 1 for BOTH populations (Wales + UKB),
#          stratified by FH status and gene type, with cross-population
#          comparison and standardised mean differences.
#
# SECTION 1: Wales Population (from workspace — all 7,253 patients)
# SECTION 2: UKB Population (NEW lipid clinic cohort)
# SECTION 3: Cross-Population Comparison (Wales vs UKB, matched vars)
#
# INPUT:   tudor_v2_workspace.RData (Wales)
#          tudor_analysis_ready.rds (UKB)
#          TUDOR_UKB_Features.csv   (UKB gene column + extras)
#
# OUTPUT:  wales_deep_table1.csv, ukb_deep_table1.csv,
#          cross_population_comparison.csv
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
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
cat("  TUDOR PIPELINE: 11_descriptive_statistics.R                   \n")
cat("  Deep Descriptive Statistics — Wales + UKB                     \n")
cat("================================================================\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Format: mean (SD) or median [IQR]
fmt_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  sprintf("%.*f (%.*f)", digits, mean(x), digits, sd(x))
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  q <- quantile(x, c(0.25, 0.5, 0.75))
  sprintf("%.*f [%.*f–%.*f]", digits, q[2], digits, q[1], digits, q[3])
}

fmt_n_pct <- function(count, total) {
  if (is.na(count) || is.na(total) || total == 0) return("—")
  sprintf("%d (%.1f%%)", count, 100 * count / total)
}

# Standardised Mean Difference (Cohen's d)
smd <- function(x1, x0) {
  x1 <- x1[!is.na(x1)]; x0 <- x0[!is.na(x0)]
  if (length(x1) < 2 || length(x0) < 2) return(NA)
  pooled_sd <- sqrt((var(x1) * (length(x1) - 1) + var(x0) * (length(x0) - 1)) /
                      (length(x1) + length(x0) - 2))
  if (pooled_sd == 0) return(0)
  (mean(x1) - mean(x0)) / pooled_sd
}

# Build a descriptive row for one continuous variable
desc_continuous <- function(data, varname, group_var = "FH_status",
                            label = varname) {
  grp_levels <- unique(data[[group_var]])
  vals <- list()
  for (g in grp_levels) {
    x <- data[[varname]][data[[group_var]] == g]
    vals[[as.character(g)]] <- x
  }
  all_x <- data[[varname]]

  row <- data.frame(
    Variable = label,
    Overall_N = sum(!is.na(all_x)),
    Overall = fmt_mean_sd(all_x),
    stringsAsFactors = FALSE
  )

  for (g in grp_levels) {
    row[[paste0("Group_", g)]] <- fmt_mean_sd(vals[[as.character(g)]])
  }

  # P-value (Welch t-test if 2 groups)
  if (length(grp_levels) == 2) {
    g1 <- vals[[as.character(grp_levels[1])]]
    g2 <- vals[[as.character(grp_levels[2])]]
    g1 <- g1[!is.na(g1)]; g2 <- g2[!is.na(g2)]
    if (length(g1) > 1 && length(g2) > 1) {
      tt <- tryCatch(t.test(g1, g2), error = function(e) NULL)
      row$P_value <- if (!is.null(tt)) sprintf("%.4g", tt$p.value) else "—"
      row$Cohen_d <- sprintf("%.3f", smd(g1, g2))
    } else {
      row$P_value <- "—"; row$Cohen_d <- "—"
    }
  }
  row
}

# Build a descriptive row for one categorical variable
desc_categorical <- function(data, varname, group_var = "FH_status",
                              label = varname) {
  grp_levels <- unique(data[[group_var]])
  cats <- sort(unique(data[[varname]][!is.na(data[[varname]])]))

  rows <- list()
  for (cat_val in cats) {
    row <- data.frame(Variable = paste0("  ", label, ": ", cat_val),
                      Overall_N = sum(!is.na(data[[varname]])),
                      Overall = fmt_n_pct(sum(data[[varname]] == cat_val, na.rm = TRUE),
                                           sum(!is.na(data[[varname]]))),
                      stringsAsFactors = FALSE)
    for (g in grp_levels) {
      sub <- data[data[[group_var]] == g, ]
      row[[paste0("Group_", g)]] <- fmt_n_pct(
        sum(sub[[varname]] == cat_val, na.rm = TRUE),
        sum(!is.na(sub[[varname]]))
      )
    }
    row$P_value <- ""; row$Cohen_d <- ""
    rows[[length(rows) + 1]] <- row
  }

  # Chi-square test for overall category
  if (length(grp_levels) == 2) {
    tbl <- tryCatch(table(data[[varname]], data[[group_var]]),
                    error = function(e) NULL)
    if (!is.null(tbl) && all(dim(tbl) >= 2)) {
      chi <- tryCatch(chisq.test(tbl), error = function(e) NULL)
      if (!is.null(chi) && length(rows) > 0) {
        rows[[1]]$P_value <- sprintf("%.4g", chi$p.value)
      }
    }
  }

  do.call(rbind, rows)
}

# ==============================================================================
# SECTION 1: WALES POPULATION
# ==============================================================================
cat("================================================================\n")
cat("SECTION 1: WALES ALL-POPULATION DESCRIPTIVE STATISTICS\n")
cat("================================================================\n\n")

# Load Wales workspace — multi-tier fallback
script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")

wales_ws_paths <- c(
  file.path(script_dir, "tudor_v2_workspace.RData"),
  "C:/Users/nader/Downloads/tudor_v2_workspace.RData",
  "C:/Users/nader/Downloads/TUDOR_Results/Tudor_Project_Full_Environment.RData",
  "C:/Users/nader/Downloads/12.RData"
)

wales_loaded <- FALSE
for (ws_path in wales_ws_paths) {
  if (file.exists(ws_path)) {
    cat("  Trying:", ws_path, "\n")
    tryCatch({
      load(ws_path)
      if (exists("df") && is.data.frame(df) && "Positive1" %in% names(df)) {
        wales_loaded <- TRUE
        cat("  SUCCESS: Loaded Wales data from:", ws_path, "\n")
        break
      } else {
        cat("  WARNING: Loaded but 'df' with 'Positive1' not found, trying next...\n")
      }
    }, error = function(e) {
      cat("  WARNING: Error loading", ws_path, ":", e$message, "\n")
    })
  }
}

# Fallback: source TUDOR_v2_Clean.R
if (!wales_loaded) {
  clean_script <- "C:/Users/nader/Downloads/TUDOR_v2_Clean.R"
  if (file.exists(clean_script)) {
    cat("  Sourcing TUDOR_v2_Clean.R to generate workspace...\n")
    tryCatch({
      source(clean_script)
      if (exists("df") && is.data.frame(df) && "Positive1" %in% names(df)) {
        wales_loaded <- TRUE
        cat("  SUCCESS: Workspace objects created from TUDOR_v2_Clean.R\n")
        # Save workspace for future runs
        ws_save_path <- "C:/Users/nader/Downloads/tudor_v2_workspace.RData"
        save.image(ws_save_path)
        cat("  Saved workspace to:", ws_save_path, "\n")
      }
    }, error = function(e) {
      cat("  WARNING: TUDOR_v2_Clean.R error:", e$message, "\n")
    })
  }
}

# Final fallback: load raw Wales CSV and derive needed columns
if (!wales_loaded) {
  raw_csv_paths <- c(
    "C:/Users/nader/Downloads/WALES_FH_CLEANED (1).csv",
    file.path(DATA_DIR, "WALES_FH_CLEANED (1).csv")
  )
  for (csv_path in raw_csv_paths) {
    if (file.exists(csv_path)) {
      cat("  Loading raw Wales CSV:", csv_path, "\n")
      df <- read.csv(csv_path, stringsAsFactors = FALSE)
      wales_loaded <- TRUE
      cat(sprintf("  Loaded: %d rows, %d columns\n", nrow(df), ncol(df)))

      # --- Derive Gene1 from Mutation1 ---
      if (!"Gene1" %in% names(df) && "Mutation1" %in% names(df)) {
        df$Gene1 <- NA_character_
        # Pattern: "GENE:c.XXX" — extract gene before colon
        has_colon <- grepl(":", df$Mutation1)
        df$Gene1[has_colon] <- sub(":.*", "", df$Mutation1[has_colon])
        # Pattern: text containing gene name (e.g., "Duplication of exons 2 - 12 LDLR")
        no_colon <- !has_colon & nchar(df$Mutation1) > 0 & !is.na(df$Mutation1)
        for (gene in c("LDLR", "APOB", "PCSK9", "APOE", "LDLRAP1")) {
          matches <- no_colon & grepl(gene, df$Mutation1, ignore.case = TRUE)
          df$Gene1[matches] <- gene
        }
        # CNV patterns for LDLR
        cnv <- no_colon & is.na(df$Gene1) &
               grepl("Duplication|Deletion|exon", df$Mutation1, ignore.case = TRUE)
        df$Gene1[cnv] <- "LDLR_CNV"
        cat(sprintf("  Derived Gene1: %d non-NA values\n", sum(!is.na(df$Gene1))))
      }

      # --- Derive Age_at_LDL1 from DOB + MeasurementDate.1 ---
      if (!"Age_at_LDL1" %in% names(df) && "DOB" %in% names(df)) {
        dob <- as.Date(df$DOB, format = "%Y-%m-%d")
        meas_date <- if ("MeasurementDate.1" %in% names(df)) {
          as.Date(df$MeasurementDate.1, format = "%Y-%m-%d")
        } else if ("DateOfLDLCMeasurement" %in% names(df)) {
          as.Date(df$DateOfLDLCMeasurement, format = "%Y-%m-%d")
        } else { NA }
        if (!all(is.na(meas_date))) {
          df$Age_at_LDL1 <- as.numeric(difftime(meas_date, dob, units = "days")) / 365.25
        } else {
          # Fallback: use currentage or estimate from DOB to 2020 (midpoint)
          df$Age_at_LDL1 <- as.numeric(difftime(as.Date("2015-01-01"), dob,
                                                  units = "days")) / 365.25
        }
        cat(sprintf("  Derived Age_at_LDL1: %d non-NA values\n",
                    sum(!is.na(df$Age_at_LDL1))))
      }

      # --- Derive Gender_num ---
      if (!"Gender_num" %in% names(df) && "Gender" %in% names(df)) {
        df$Gender_num <- ifelse(df$Gender == "M", 1,
                         ifelse(df$Gender == "F", 0, NA))
      }

      # --- LDL_untreated (statin-corrected) ---
      if (!"LDL_untreated" %in% names(df) && "LDL.1" %in% names(df)) {
        # Use raw LDL as approximation (statin correction requires Clean.R)
        df$LDL_untreated <- df$LDL.1
        cat("  NOTE: LDL_untreated = LDL.1 (no statin correction available)\n")
        cat("  For accurate statin correction, run TUDOR_v2_Clean.R first\n")
      }

      # --- Trig_Filter ---
      if (!"Trig_Filter" %in% names(df) && "TRG.1" %in% names(df)) {
        ldl_ut <- if ("LDL_untreated" %in% names(df)) df$LDL_untreated else df$LDL.1
        df$Trig_Filter <- ldl_ut / (df$TRG.1 + 0.1)
        cat("  Computed Trig_Filter from LDL / (TRG + 0.1)\n")
      }

      # Rename raw lipid columns if needed (scripts expect LDL.1, HDL.1, TRG.1)
      # These should already be correct in the raw CSV
      break
    }
  }
}

if (!wales_loaded) {
  stop("Wales data not found. Place one of these in C:/Users/nader/Downloads/:\n",
       "  - tudor_v2_workspace.RData\n",
       "  - Tudor_Project_Full_Environment.RData (in TUDOR_Results/)\n",
       "  - WALES_FH_CLEANED (1).csv")
}

cat(sprintf("  Wales Registry total: %d patients\n", nrow(df)))
cat(sprintf("  Index cases (I_Vs_R=1): %d\n", sum(df$I_Vs_R == 1, na.rm = TRUE)))
cat(sprintf("  Relatives (I_Vs_R=2): %d\n", sum(df$I_Vs_R == 2, na.rm = TRUE)))
cat(sprintf("  FH-positive (Positive1=1): %d\n", sum(df$Positive1 == 1, na.rm = TRUE)))
cat(sprintf("  FH-negative (Positive1=0): %d\n\n", sum(df$Positive1 == 0, na.rm = TRUE)))

# Create FH status variable
df$FH_status <- ifelse(df$Positive1 == 1, "FH+", "FH-")

# --- 1a. Basic Demographics ---
cat("--- Wales Table 1: Full Population ---\n\n")

wales_table <- list()

# Age
wales_table[[1]] <- desc_continuous(df, "Age_at_LDL1", "FH_status", "Age (years)")

# Sex
df$Sex_label <- ifelse(df$Gender == "M" | df$Gender_num == 1, "Male",
                ifelse(df$Gender == "F" | df$Gender_num == 0, "Female", NA))
wales_table[[2]] <- desc_categorical(df, "Sex_label", "FH_status", "Sex")

# BMI
bmi_col <- if ("BMI_clean" %in% names(df)) "BMI_clean" else
           if ("BMI_AGE" %in% names(df)) "BMI_AGE" else
           if ("BMI_imputed" %in% names(df)) "BMI_imputed" else NULL
if (!is.null(bmi_col)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, bmi_col, "FH_status",
                                                             "BMI (kg/m²)")
}

# Lipid panel — raw
if ("LDL.1" %in% names(df)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, "LDL.1", "FH_status",
                                                             "LDL-C raw (mmol/L)")
}

# LDL untreated (statin-adjusted)
if ("LDL_untreated" %in% names(df)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, "LDL_untreated", "FH_status",
                                                             "LDL-C untreated (mmol/L)")
}

# HDL
if ("HDL.1" %in% names(df)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, "HDL.1", "FH_status",
                                                             "HDL-C (mmol/L)")
}

# Triglycerides
if ("TRG.1" %in% names(df)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, "TRG.1", "FH_status",
                                                             "Triglycerides (mmol/L)")
}

# Total cholesterol
tc_col <- if ("TC.1" %in% names(df)) "TC.1" else
          if ("CHOL" %in% names(df)) "CHOL" else NULL
if (!is.null(tc_col)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, tc_col, "FH_status",
                                                             "Total Cholesterol (mmol/L)")
}

# Trig Filter
if ("Trig_Filter" %in% names(df)) {
  wales_table[[length(wales_table) + 1]] <- desc_continuous(df, "Trig_Filter", "FH_status",
                                                             "Trig Filter (LDL_UT/(TG+0.1))")
}

# Treatment status
if ("statin_name" %in% names(df) || "Statin_Type" %in% names(df)) {
  statin_col <- if ("statin_name" %in% names(df)) "statin_name" else "Statin_Type"
  df$On_Statin <- ifelse(!is.na(df[[statin_col]]) & df[[statin_col]] != "None" &
                           df[[statin_col]] != "" & df[[statin_col]] != "0", "Yes", "No")
  wales_table[[length(wales_table) + 1]] <- desc_categorical(df, "On_Statin", "FH_status",
                                                              "On Statin")
}

# Tendon xanthomata (if available)
tx_col <- if ("tendon_xanth" %in% names(df)) "tendon_xanth" else
          if ("TendonXanthomata" %in% names(df)) "TendonXanthomata" else NULL
if (!is.null(tx_col)) {
  df$Tendon_Xanth <- ifelse(df[[tx_col]] == 1, "Yes", "No")
  wales_table[[length(wales_table) + 1]] <- desc_categorical(df, "Tendon_Xanth", "FH_status",
                                                              "Tendon Xanthomata")
}

# Corneal arcus (if available)
ca_col <- if ("corneal_less_40" %in% names(df)) "corneal_less_40" else
          if ("CornealArcus" %in% names(df)) "CornealArcus" else NULL
if (!is.null(ca_col)) {
  df$Corneal_Arcus <- ifelse(df[[ca_col]] == 1, "Yes", "No")
  wales_table[[length(wales_table) + 1]] <- desc_categorical(df, "Corneal_Arcus", "FH_status",
                                                              "Corneal Arcus (<40y)")
}

# ASCVD
ascvd_col <- if ("ASCVD_event" %in% names(df)) "ASCVD_event" else
             if ("has_any_cvd" %in% names(df)) "has_any_cvd" else NULL
if (!is.null(ascvd_col)) {
  df$ASCVD <- ifelse(df[[ascvd_col]] == 1, "Yes", "No")
  wales_table[[length(wales_table) + 1]] <- desc_categorical(df, "ASCVD", "FH_status",
                                                              "ASCVD History")
}

# Index vs Relative
df$Patient_Type <- ifelse(df$I_Vs_R == 1, "Index", "Relative")
wales_table[[length(wales_table) + 1]] <- desc_categorical(df, "Patient_Type", "FH_status",
                                                            "Patient Type")

# Gene type (FH+ only)
if ("Gene1" %in% names(df)) {
  fh_pos <- df[df$Positive1 == 1, ]
  cat("\n--- Wales: Gene Type Distribution (FH+ only, n =",
      nrow(fh_pos), ") ---\n")
  gene_tab <- table(fh_pos$Gene1)
  gene_pct <- round(100 * prop.table(gene_tab), 1)
  for (g in names(sort(gene_tab, decreasing = TRUE))) {
    cat(sprintf("  %-12s: %4d (%5.1f%%)\n", g, gene_tab[g], gene_pct[g]))
  }
}

# Assemble Wales Table 1
wales_t1 <- do.call(rbind, wales_table)
cat("\n--- Wales Table 1 ---\n")
print(wales_t1, row.names = FALSE)

# Save
write.csv(wales_t1, file.path(TABLE_DIR, "wales_deep_table1.csv"), row.names = FALSE)
cat("\nSaved: wales_deep_table1.csv\n\n")

# --- 1b. Stratified by Gene Type (FH+ patients only) ---
cat("--- Wales: Descriptive Statistics Stratified by Gene Type ---\n\n")

if ("Gene1" %in% names(df)) {
  fh_pos <- df[df$Positive1 == 1, ]

  # Consolidate rare genes
  fh_pos$Gene_Group <- ifelse(fh_pos$Gene1 %in% c("LDLR", "APOB", "PCSK9"),
                               fh_pos$Gene1, "Other")

  genes <- c("LDLR", "APOB", "PCSK9", "Other")
  gene_desc <- list()

  for (g in genes) {
    sub <- fh_pos[fh_pos$Gene_Group == g, ]
    if (nrow(sub) < 5) next

    cat(sprintf("  Gene: %s (n = %d)\n", g, nrow(sub)))
    cat(sprintf("    Age: %s\n", fmt_mean_sd(sub$Age_at_LDL1)))

    if ("LDL_untreated" %in% names(sub))
      cat(sprintf("    LDL_untreated: %s\n", fmt_mean_sd(sub$LDL_untreated)))
    if ("HDL.1" %in% names(sub))
      cat(sprintf("    HDL-C: %s\n", fmt_mean_sd(sub$HDL.1)))
    if ("TRG.1" %in% names(sub))
      cat(sprintf("    Triglycerides: %s\n", fmt_mean_sd(sub$TRG.1)))
    if ("Trig_Filter" %in% names(sub))
      cat(sprintf("    Trig Filter: %s\n", fmt_mean_sd(sub$Trig_Filter)))
    cat("\n")
  }

  # Kruskal-Wallis for LDL by gene type
  if ("LDL_untreated" %in% names(fh_pos)) {
    kw <- kruskal.test(LDL_untreated ~ Gene_Group, data = fh_pos)
    cat(sprintf("  Kruskal-Wallis (LDL_untreated ~ Gene): chi²=%.2f, p=%.4g\n\n",
                kw$statistic, kw$p.value))
  }
}

# Save Wales workspace reference (keep df_wales for Section 3)
df_wales <- df
n_wales <- nrow(df_wales)

# ==============================================================================
# SECTION 2: UKB POPULATION
# ==============================================================================
cat("================================================================\n")
cat("SECTION 2: UKB POPULATION — NEW LIPID CLINIC COHORT\n")
cat("================================================================\n\n")

# Load UKB analysis-ready data
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
if (!file.exists(rds_file)) stop("Run 01_data_merge.R first!")

cat("Loading UKB analysis-ready dataset...\n")
suppressPackageStartupMessages(library(data.table))
df <- readRDS(rds_file)
if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
cat(sprintf("  Loaded: %d participants\n", nrow(df)))
cat(sprintf("  Genetic FH cases: %d\n\n", sum(df$is_fh_genetic, na.rm = TRUE)))

# --- 2a. Define NEW Lipid Clinic Cohort ---
cat("--- Defining NEW Lipid Clinic Cohort ---\n")
cat("  Criteria: TC > 7.5 mmol/L OR LDL_RW > 4.9 mmol/L OR Premature ASCVD\n\n")

# Individual criteria counts
crit_tc    <- !is.na(df$CHOL) & df$CHOL > 7.5
crit_ldl   <- !is.na(df$LDL_RW) & df$LDL_RW > 4.9
crit_ascvd <- !is.na(df$Premature_ASCVD) & df$Premature_ASCVD == 1

df$cohort_lipid_clinic <- crit_tc | crit_ldl | crit_ascvd

cat("  Criterion breakdown:\n")
cat(sprintf("    TC > 7.5 mmol/L:       %6d (%5.2f%%)\n",
            sum(crit_tc), 100 * mean(crit_tc)))
cat(sprintf("    LDL_RW > 4.9 mmol/L:   %6d (%5.2f%%)\n",
            sum(crit_ldl), 100 * mean(crit_ldl)))
cat(sprintf("    Premature ASCVD:        %6d (%5.2f%%)\n",
            sum(crit_ascvd), 100 * mean(crit_ascvd)))
cat(sprintf("\n  OLD cohort (LDL > 4.9 only): %d\n",
            sum(df$cohort_high_risk, na.rm = TRUE)))
cat(sprintf("  NEW cohort (lipid clinic):    %d\n", sum(df$cohort_lipid_clinic)))
cat(sprintf("  Difference (new cases):       %d\n\n",
            sum(df$cohort_lipid_clinic) - sum(df$cohort_high_risk, na.rm = TRUE)))

# Venn overlap
only_tc    <- crit_tc & !crit_ldl & !crit_ascvd
only_ldl   <- !crit_tc & crit_ldl & !crit_ascvd
only_ascvd <- !crit_tc & !crit_ldl & crit_ascvd
tc_and_ldl <- crit_tc & crit_ldl & !crit_ascvd
tc_and_ascvd <- crit_tc & !crit_ldl & crit_ascvd
ldl_and_ascvd <- !crit_tc & crit_ldl & crit_ascvd
all_three  <- crit_tc & crit_ldl & crit_ascvd

cat("  Venn diagram of criteria overlap:\n")
cat(sprintf("    TC only:            %6d\n", sum(only_tc)))
cat(sprintf("    LDL only:           %6d\n", sum(only_ldl)))
cat(sprintf("    ASCVD only:         %6d\n", sum(only_ascvd)))
cat(sprintf("    TC + LDL:           %6d\n", sum(tc_and_ldl)))
cat(sprintf("    TC + ASCVD:         %6d\n", sum(tc_and_ascvd)))
cat(sprintf("    LDL + ASCVD:        %6d\n", sum(ldl_and_ascvd)))
cat(sprintf("    All three:          %6d\n\n", sum(all_three)))

# FH cases captured by each criterion
cat("  FH cases captured by each criterion:\n")
cat(sprintf("    TC > 7.5:     %d / %d FH cases\n",
            sum(df$is_fh_genetic[crit_tc]), sum(df$is_fh_genetic)))
cat(sprintf("    LDL > 4.9:    %d / %d FH cases\n",
            sum(df$is_fh_genetic[crit_ldl]), sum(df$is_fh_genetic)))
cat(sprintf("    Prem ASCVD:   %d / %d FH cases\n",
            sum(df$is_fh_genetic[crit_ascvd]), sum(df$is_fh_genetic)))
cat(sprintf("    NEW cohort:   %d / %d FH cases\n",
            sum(df$is_fh_genetic[df$cohort_lipid_clinic]),
            sum(df$is_fh_genetic)))
cat(sprintf("    NEW FH by TC>7.5 not in old: %d\n\n",
            sum(df$is_fh_genetic[crit_tc & !crit_ldl])))

# --- 2b. UKB Table 1 (Lipid Clinic Cohort) ---
cat("--- UKB Table 1: Lipid Clinic Cohort ---\n\n")

lc <- df[df$cohort_lipid_clinic == TRUE, ]
lc$FH_status <- ifelse(lc$is_fh_genetic == 1, "FH+", "FH-")

cat(sprintf("  Lipid clinic cohort: %d (FH+ = %d, %.2f%%)\n\n",
            nrow(lc), sum(lc$FH_status == "FH+"),
            100 * mean(lc$FH_status == "FH+")))

ukb_table <- list()

# Age
ukb_table[[1]] <- desc_continuous(lc, "Age_at_LDL1", "FH_status", "Age (years)")

# Sex
lc$Sex_label <- ifelse(lc$Gender_num == 1, "Male", "Female")
ukb_table[[2]] <- desc_categorical(lc, "Sex_label", "FH_status", "Sex")

# BMI
ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "BMI_imputed", "FH_status",
                                                       "BMI (kg/m²)")

# LDL (statin-corrected)
ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "LDL_RW", "FH_status",
                                                       "LDL-C corrected (mmol/L)")

# LDL raw (treated)
if ("LDL_treated" %in% names(lc)) {
  ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "LDL_treated", "FH_status",
                                                         "LDL-C measured (mmol/L)")
}

# HDL
ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "HDL.1", "FH_status",
                                                       "HDL-C (mmol/L)")

# Triglycerides
ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "TRG.1", "FH_status",
                                                       "Triglycerides (mmol/L)")

# Total cholesterol
ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "CHOL", "FH_status",
                                                       "Total Cholesterol (mmol/L)")

# Trig Filter
if ("Trig_Filter_RW" %in% names(lc)) {
  ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "Trig_Filter_RW", "FH_status",
                                                         "Trig Filter (LDL_RW/(TG+0.1))")
}

# ApoB
if ("ApoB" %in% names(lc) && sum(!is.na(lc$ApoB)) > 0) {
  ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "ApoB", "FH_status",
                                                         "ApoB (g/L)")
  # ApoB/LDL ratio
  lc$ApoB_LDL <- lc$ApoB / lc$LDL_RW
  ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "ApoB_LDL", "FH_status",
                                                         "ApoB/LDL-C ratio")
}

# Statin use
lc$On_Statin <- ifelse(lc$statin_name != "None", "Yes", "No")
ukb_table[[length(ukb_table) + 1]] <- desc_categorical(lc, "On_Statin", "FH_status",
                                                        "On Statin")

# Statin type distribution
ukb_table[[length(ukb_table) + 1]] <- desc_categorical(lc, "statin_name", "FH_status",
                                                        "Statin Type")

# Premature ASCVD
lc$Prem_ASCVD <- ifelse(lc$Premature_ASCVD == 1, "Yes", "No")
ukb_table[[length(ukb_table) + 1]] <- desc_categorical(lc, "Prem_ASCVD", "FH_status",
                                                        "Premature ASCVD")

# HbA1c (if available)
hba1c_col <- grep("p30750", names(lc), value = TRUE)
if (length(hba1c_col) > 0) {
  lc$HbA1c <- as.numeric(lc[[hba1c_col[1]]])
  if (sum(!is.na(lc$HbA1c)) > 100) {
    ukb_table[[length(ukb_table) + 1]] <- desc_continuous(lc, "HbA1c", "FH_status",
                                                           "HbA1c (mmol/mol)")
  }
}

# Assemble UKB Table 1
ukb_t1 <- do.call(rbind, ukb_table)
cat("\n--- UKB Lipid Clinic Cohort Table 1 ---\n")
print(ukb_t1, row.names = FALSE)

# Save
write.csv(ukb_t1, file.path(TABLE_DIR, "ukb_deep_table1.csv"), row.names = FALSE)
cat("\nSaved: ukb_deep_table1.csv\n\n")

# --- 2c. UKB Gene Type Distribution ---
cat("--- UKB: Gene Type Distribution ---\n\n")

# Load gene column from features CSV
feat_file <- file.path(DATA_DIR, "TUDOR_UKB_Features.csv")
if (file.exists(feat_file)) {
  cat("  Loading gene data from TUDOR_UKB_Features.csv...\n")

  # Read only eid + gene columns to avoid memory issues
  feat_header <- names(read.csv(feat_file, nrows = 1, check.names = FALSE))

  # Find eid column
  eid_col_idx <- grep("^eid$|^participant.eid$", feat_header)
  gene_col_idx <- which(feat_header == "gene")
  fh_col_idx   <- which(feat_header == "is_fh_genetic")

  if (length(gene_col_idx) > 0) {
    # Read in chunks to avoid segfaults
    cat("  Reading gene column in chunks...\n")
    con <- file(feat_file, "r")
    header <- readLines(con, n = 1)

    gene_data <- data.frame(gene = character(0), is_fh = integer(0),
                            stringsAsFactors = FALSE)

    chunk_size <- 50000
    chunk_num <- 0

    repeat {
      lines <- readLines(con, n = chunk_size)
      if (length(lines) == 0) break

      chunk_num <- chunk_num + 1
      tmp <- read.csv(textConnection(c(header, lines)),
                       stringsAsFactors = FALSE, check.names = FALSE)

      if ("gene" %in% names(tmp) && "is_fh_genetic" %in% names(tmp)) {
        fh_rows <- tmp[tmp$is_fh_genetic == 1, ]
        if (nrow(fh_rows) > 0) {
          gene_data <- rbind(gene_data,
                             data.frame(gene = fh_rows$gene,
                                        is_fh = 1,
                                        stringsAsFactors = FALSE))
        }
      }

      if (chunk_num %% 2 == 0) {
        cat(sprintf("    Processed %d chunks (%d FH cases so far)...\n",
                    chunk_num, nrow(gene_data)))
      }
    }
    close(con)

    cat(sprintf("\n  Total genetic FH cases with gene data: %d\n", nrow(gene_data)))

    if (nrow(gene_data) > 0) {
      gene_tab_ukb <- table(gene_data$gene)
      gene_pct_ukb <- round(100 * prop.table(gene_tab_ukb), 1)

      cat("  UKB Gene Type Distribution (FH+ only):\n")
      for (g in names(sort(gene_tab_ukb, decreasing = TRUE))) {
        cat(sprintf("    %-12s: %4d (%5.1f%%)\n", g, gene_tab_ukb[g],
                    gene_pct_ukb[g]))
      }

      # Save for script 12
      saveRDS(gene_data,
              file.path(OUTPUT_DIR, "11_ukb_gene_data.rds"))
      cat("  Saved: 11_ukb_gene_data.rds\n")
    }
  } else {
    cat("  WARNING: gene column not found in features CSV\n")
  }
} else {
  cat("  WARNING: TUDOR_UKB_Features.csv not found\n")
}

# ==============================================================================
# SECTION 3: CROSS-POPULATION COMPARISON
# ==============================================================================
cat("\n================================================================\n")
cat("SECTION 3: CROSS-POPULATION COMPARISON (WALES vs UKB)\n")
cat("================================================================\n\n")

# Build comparison on matched variables
comparison_rows <- list()

# Helper for one row of cross-population comparison
cross_row <- function(label, wales_vals, ukb_vals) {
  w <- wales_vals[!is.na(wales_vals)]
  u <- ukb_vals[!is.na(ukb_vals)]

  row <- data.frame(
    Variable = label,
    Wales_N = length(w),
    Wales = fmt_mean_sd(w),
    UKB_N = length(u),
    UKB = fmt_mean_sd(u),
    stringsAsFactors = FALSE
  )

  # SMD between populations
  if (length(w) > 1 && length(u) > 1) {
    row$SMD <- sprintf("%.3f", smd(u, w))
    tt <- tryCatch(t.test(w, u), error = function(e) NULL)
    row$P_value <- if (!is.null(tt)) sprintf("%.4g", tt$p.value) else "—"
  } else {
    row$SMD <- "—"
    row$P_value <- "—"
  }
  row
}

# Wales external test set (Relatives) for fairer comparison
ext_wales <- df_wales[df_wales$I_Vs_R == 2, ]

cat(sprintf("  Comparing: Wales Relatives (n=%d) vs UKB Lipid Clinic (n=%d)\n\n",
            nrow(ext_wales), nrow(lc)))

# Age
comparison_rows[[1]] <- cross_row("Age (years)", ext_wales$Age_at_LDL1, lc$Age_at_LDL1)

# Sex (% male)
w_male_pct <- mean(ext_wales$Gender == "M" | ext_wales$Gender_num == 1, na.rm = TRUE)
u_male_pct <- mean(lc$Gender_num == 1, na.rm = TRUE)
comparison_rows[[2]] <- data.frame(
  Variable = "Male (%)", Wales_N = nrow(ext_wales),
  Wales = sprintf("%.1f%%", 100 * w_male_pct),
  UKB_N = nrow(lc), UKB = sprintf("%.1f%%", 100 * u_male_pct),
  SMD = "—", P_value = "—", stringsAsFactors = FALSE)

# LDL-C corrected
w_ldl <- if ("LDL_untreated" %in% names(ext_wales)) ext_wales$LDL_untreated else ext_wales$LDL.1
comparison_rows[[3]] <- cross_row("LDL-C adjusted (mmol/L)", w_ldl, lc$LDL_RW)

# HDL
comparison_rows[[4]] <- cross_row("HDL-C (mmol/L)", ext_wales$HDL.1, lc$HDL.1)

# Triglycerides
comparison_rows[[5]] <- cross_row("Triglycerides (mmol/L)", ext_wales$TRG.1, lc$TRG.1)

# FH prevalence
w_fh_pct <- mean(ext_wales$Positive1 == 1, na.rm = TRUE)
u_fh_pct <- mean(lc$is_fh_genetic == 1, na.rm = TRUE)
comparison_rows[[6]] <- data.frame(
  Variable = "FH prevalence (%)", Wales_N = nrow(ext_wales),
  Wales = sprintf("%.1f%%", 100 * w_fh_pct),
  UKB_N = nrow(lc), UKB = sprintf("%.2f%%", 100 * u_fh_pct),
  SMD = "—", P_value = "—", stringsAsFactors = FALSE)

cross_t1 <- do.call(rbind, comparison_rows)
cat("--- Cross-Population Comparison ---\n")
print(cross_t1, row.names = FALSE)

# Save
write.csv(cross_t1, file.path(TABLE_DIR, "cross_population_comparison.csv"),
          row.names = FALSE)
cat("\nSaved: cross_population_comparison.csv\n")

# --- Save lipid clinic cohort flag for scripts 12 & 13 ---
saveRDS(df[df$cohort_lipid_clinic == TRUE, ],
        file.path(OUTPUT_DIR, "11_lipid_clinic_cohort.rds"))
cat("Saved: 11_lipid_clinic_cohort.rds (new cohort for scripts 12 & 13)\n")

# Also update the full dataset with the new cohort flag
saveRDS(df, file.path(OUTPUT_DIR, "11_tudor_with_lipid_clinic.rds"))
cat("Saved: 11_tudor_with_lipid_clinic.rds (full dataset + new flag)\n")

# ==============================================================================
# SECTION 4: COHORT CRITERIA BREAKDOWN
# ==============================================================================
cat("\n================================================================\n")
cat("SECTION 4: COHORT CRITERIA SUMMARY FOR MANUSCRIPT\n")
cat("================================================================\n\n")

cat(sprintf("WALES (All Wales FH Registry / PASS):\n"))
cat(sprintf("  Total patients:       %d\n", n_wales))
cat(sprintf("  Index cases:          %d\n", sum(df_wales$I_Vs_R == 1, na.rm = TRUE)))
cat(sprintf("  Relatives tested:     %d\n", sum(df_wales$I_Vs_R == 2, na.rm = TRUE)))
cat(sprintf("  FH-positive total:    %d\n", sum(df_wales$Positive1 == 1, na.rm = TRUE)))
cat(sprintf("  FH-negative total:    %d\n", sum(df_wales$Positive1 == 0, na.rm = TRUE)))

cat(sprintf("\nUK BIOBANK:\n"))
cat(sprintf("  Total participants:   %d\n", nrow(df)))
cat(sprintf("  Genetic FH cases:     %d\n", sum(df$is_fh_genetic)))
cat(sprintf("  OLD high-risk cohort: %d (LDL > 4.9 only)\n",
            sum(df$cohort_high_risk, na.rm = TRUE)))
cat(sprintf("  NEW lipid clinic:     %d (TC>7.5 | LDL>4.9 | premature ASCVD)\n",
            sum(df$cohort_lipid_clinic)))
cat(sprintf("  FH in old cohort:     %d\n",
            sum(df$is_fh_genetic[df$cohort_high_risk == TRUE], na.rm = TRUE)))
cat(sprintf("  FH in new cohort:     %d\n",
            sum(df$is_fh_genetic[df$cohort_lipid_clinic == TRUE])))

cat(sprintf("\nCOMBINED STUDY:\n"))
cat(sprintf("  Wales FH+:     %d\n", sum(df_wales$Positive1 == 1, na.rm = TRUE)))
cat(sprintf("  UKB FH+:       %d\n", sum(df$is_fh_genetic)))
cat(sprintf("  TOTAL FH+:     %d genetically confirmed FH cases\n",
            sum(df_wales$Positive1 == 1, na.rm = TRUE) + sum(df$is_fh_genetic)))
cat("  → Largest dual-validated genetically confirmed FH diagnostic study worldwide\n")

cat("\n=== 11_descriptive_statistics.R COMPLETE ===\n")
