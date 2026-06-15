# ==============================================================================
# TUDOR PIPELINE: STEP 09 — NMR METABOLOMICS IN FH
# ==============================================================================
# PURPOSE: First large-scale NMR lipoprotein particle profiling of
#          genetically-confirmed FH in UK Biobank (N ~ 120-275k with NMR).
#
# NOVEL FINDINGS TARGETED:
#   1. FH "metabolic signature" — elevated LDL particles, NORMAL VLDL
#   2. Cholesterol-enriched LDL particles (larger size, lower ApoB/LDL)
#   3. Discordant ApoB-cholesterol relationship (fewer particles, more chol each)
#   4. VLDL normality as monogenic vs polygenic discriminator
#   5. NMR-level validation of the Trig Filter mechanism
#
# INPUT:  tudor_analysis_ready.rds + ukb_nmr_batch1-4.csv
# OUTPUT: NMR profile tables, correlation analyses, mechanistic evidence
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

cat("=== TUDOR PIPELINE: 09_nmr_metabolomics.R ===\n")
cat("NMR lipoprotein particle profiling in genetic FH\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
df <- readRDS(rds_file)
cat("Loaded:", nrow(df), "participants\n")

# --- 1a. Load and merge NMR batches ---
# v2: 5 small batches (10-11 fields each) to avoid DataTooLarge API errors
nmr_files <- list(
  file.path(DATA_DIR, "ukb_nmr_a.csv"),
  file.path(DATA_DIR, "ukb_nmr_b.csv"),
  file.path(DATA_DIR, "ukb_nmr_c.csv"),
  file.path(DATA_DIR, "ukb_nmr_d.csv"),
  file.path(DATA_DIR, "ukb_nmr_e.csv")
)

nmr_found <- FALSE
nmr_all <- NULL

for (nf in nmr_files) {
  if (file.exists(nf)) {
    cat("Loading NMR:", basename(nf), "...")
    tmp <- fread(nf)
    cat(" ", ncol(tmp) - 1, "metabolites,", nrow(tmp), "participants\n")
    if (is.null(nmr_all)) {
      nmr_all <- tmp
    } else {
      nmr_all <- merge(nmr_all, tmp, by = "participant.eid", all = TRUE)
    }
    nmr_found <- TRUE
  }
}

if (!nmr_found) {
  cat("\n")
  cat("============================================================\n")
  cat("WARNING: No NMR data files found.\n")
  cat("Run 00c_extract_nmr_dm.sh on UKB-RAP first, then download\n")
  cat("ukb_nmr_batch1-4.csv to:", DATA_DIR, "\n")
  cat("============================================================\n")
  cat("\n=== 09_nmr_metabolomics.R SKIPPED (no NMR data) ===\n")
  quit(save = "no")
}

# Merge NMR with main dataset
cat("\nMerging NMR with main dataset...\n")
df <- merge(df, nmr_all, by = "participant.eid", all.x = TRUE)

# Identify NMR columns (all participant.p23xxx_i0 columns)
nmr_cols <- grep("^participant\\.p23[4-6][0-9]{2}_i0$", names(df), value = TRUE)
cat("NMR columns available:", length(nmr_cols), "\n")

# Count participants with NMR data
nmr_coverage <- rowSums(!is.na(df[, ..nmr_cols])) > 0
cat("Participants with any NMR data:", sum(nmr_coverage), "\n")

# Filter to those with NMR data
df_nmr <- df[nmr_coverage, ]
cat("Analysis subset:", nrow(df_nmr), "participants\n")
cat("  FH cases with NMR:", sum(df_nmr$is_fh_genetic), "\n")
cat("  Non-FH with NMR:", sum(!df_nmr$is_fh_genetic), "\n\n")

# ==============================================================================
# 2. NMR FIELD MAPPING
# ==============================================================================
# Map UKB field IDs to human-readable metabolite names
# This is a best-effort mapping based on the Nightingale NMR platform
# The actual mapping depends on which fields exist in the dataset

# Create a lookup table (will be populated dynamically)
create_nmr_labels <- function(cols) {
  # Known Nightingale NMR field mappings (UKB Category 220)
  known_map <- list(
    # Lipoprotein subclass total lipids
    "p23400" = "XXL_VLDL_lipids",  "p23401" = "XL_VLDL_lipids",
    "p23402" = "L_VLDL_lipids",    "p23403" = "M_VLDL_lipids",
    "p23404" = "S_VLDL_lipids",    "p23405" = "XS_VLDL_lipids",
    "p23406" = "IDL_lipids",       "p23407" = "L_LDL_lipids",
    "p23408" = "M_LDL_lipids",     "p23409" = "S_LDL_lipids",
    "p23410" = "XL_HDL_lipids",    "p23411" = "L_HDL_lipids",
    "p23412" = "M_HDL_lipids",     "p23413" = "S_HDL_lipids",

    # Cholesterol in subclasses
    "p23414" = "XXL_VLDL_CE", "p23415" = "XL_VLDL_CE",
    "p23416" = "L_VLDL_CE",   "p23417" = "M_VLDL_CE",
    "p23418" = "S_VLDL_CE",   "p23419" = "XS_VLDL_CE",
    "p23420" = "IDL_CE",       "p23421" = "L_LDL_CE",
    "p23422" = "M_LDL_CE",    "p23423" = "S_LDL_CE",
    "p23424" = "XL_HDL_CE",   "p23425" = "L_HDL_CE",
    "p23426" = "M_HDL_CE",    "p23427" = "S_HDL_CE",

    # Free cholesterol
    "p23428" = "XXL_VLDL_FC", "p23429" = "XL_VLDL_FC",
    "p23430" = "L_VLDL_FC",   "p23431" = "M_VLDL_FC",
    "p23432" = "S_VLDL_FC",   "p23433" = "XS_VLDL_FC",
    "p23434" = "IDL_FC",       "p23435" = "L_LDL_FC",
    "p23436" = "M_LDL_FC",    "p23437" = "S_LDL_FC",
    "p23438" = "XL_HDL_FC",   "p23439" = "L_HDL_FC",
    "p23440" = "M_HDL_FC",    "p23441" = "S_HDL_FC",

    # Triglycerides in subclasses
    "p23442" = "XXL_VLDL_TG", "p23443" = "XL_VLDL_TG",
    "p23444" = "L_VLDL_TG",   "p23445" = "M_VLDL_TG",
    "p23446" = "S_VLDL_TG",   "p23447" = "XS_VLDL_TG",
    "p23448" = "IDL_TG",       "p23449" = "L_LDL_TG",
    "p23450" = "M_LDL_TG",    "p23451" = "S_LDL_TG",
    "p23452" = "XL_HDL_TG",   "p23453" = "L_HDL_TG",
    "p23454" = "M_HDL_TG",    "p23455" = "S_HDL_TG",

    # Phospholipids in subclasses
    "p23456" = "XXL_VLDL_PL", "p23457" = "XL_VLDL_PL",
    "p23458" = "L_VLDL_PL",   "p23459" = "M_VLDL_PL",
    "p23460" = "S_VLDL_PL",   "p23461" = "XS_VLDL_PL",
    "p23462" = "IDL_PL",       "p23463" = "L_LDL_PL",
    "p23464" = "M_LDL_PL",    "p23465" = "S_LDL_PL",
    "p23466" = "XL_HDL_PL",   "p23467" = "L_HDL_PL",
    "p23468" = "M_HDL_PL",    "p23469" = "S_HDL_PL",

    # Particle concentrations
    "p23470" = "XXL_VLDL_P", "p23471" = "XL_VLDL_P",
    "p23472" = "L_VLDL_P",   "p23473" = "M_VLDL_P",
    "p23474" = "S_VLDL_P",   "p23475" = "XS_VLDL_P",
    "p23476" = "IDL_P",       "p23477" = "L_LDL_P",
    "p23478" = "M_LDL_P",    "p23479" = "S_LDL_P",
    "p23480" = "XL_HDL_P",   "p23481" = "L_HDL_P",
    "p23482" = "M_HDL_P",    "p23483" = "S_HDL_P",

    # Particle sizes
    "p23484" = "VLDL_size",   "p23485" = "LDL_size",
    "p23486" = "HDL_size",

    # Clinical lipids (NMR-derived)
    "p23487" = "NMR_TotalChol",    "p23488" = "NMR_LDL_C",
    "p23489" = "NMR_HDL_C",        "p23490" = "NMR_TotalTG",
    "p23491" = "NMR_VLDL_C",       "p23492" = "Remnant_Chol",

    # Apolipoproteins (NMR)
    "p23493" = "NMR_ApoB",         "p23494" = "NMR_ApoA1",

    # GlycA (inflammation)
    "p23495" = "GlycA",

    # Amino acids
    "p23496" = "Alanine",      "p23497" = "Glutamine",
    "p23498" = "Glycine",      "p23499" = "Histidine",
    "p23500" = "Isoleucine",   "p23501" = "Leucine",
    "p23502" = "Valine",       "p23503" = "Phenylalanine",
    "p23504" = "Tyrosine",     "p23505" = "BCAA",

    # Fatty acids
    "p23506" = "TotalFA",      "p23507" = "Unsat_degree",
    "p23508" = "Omega3_FA",    "p23509" = "Omega6_FA",
    "p23510" = "PUFA",         "p23511" = "MUFA",
    "p23512" = "SFA",          "p23513" = "DHA",
    "p23514" = "LA",           "p23515" = "Omega3_pct",
    "p23516" = "Omega6_pct",   "p23517" = "PUFA_pct",
    "p23518" = "MUFA_pct",     "p23519" = "SFA_pct",
    "p23520" = "DHA_pct",      "p23521" = "LA_pct",

    # Glycolysis / ketone bodies
    "p23522" = "Glucose",      "p23523" = "Lactate",
    "p23524" = "Pyruvate",     "p23525" = "Citrate",
    "p23526" = "bOHbutyrate",  "p23527" = "Acetone",
    "p23528" = "Acetoacetate", "p23529" = "Acetate",
    "p23530" = "Creatinine",   "p23531" = "Albumin"
  )

  labels <- character(length(cols))
  for (i in seq_along(cols)) {
    # Extract field number: participant.p23400_i0 -> p23400
    field <- gsub("participant\\.(p[0-9]+)_i0", "\\1", cols[i])
    if (field %in% names(known_map)) {
      labels[i] <- known_map[[field]]
    } else {
      labels[i] <- field  # Use raw field ID if not mapped
    }
  }
  return(labels)
}

nmr_labels <- create_nmr_labels(nmr_cols)
names(nmr_labels) <- nmr_cols
cat("Mapped", sum(nmr_labels != gsub("participant\\.(p[0-9]+)_i0", "\\1", nmr_cols)),
    "of", length(nmr_cols), "NMR fields to metabolite names\n\n")


# ==============================================================================
# A. FH vs NON-FH NMR PROFILE COMPARISON
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS A: NMR Metabolite Profile — FH vs Non-FH\n")
cat("(Within high-risk cohort with NMR data)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

hr_nmr <- df_nmr[df_nmr$cohort_high_risk == TRUE & !is.na(df_nmr$tudor_prob), ]
cat("High-risk with NMR:", nrow(hr_nmr), "(FH:", sum(hr_nmr$is_fh_genetic), ")\n\n")

# Compare each NMR metabolite
nmr_results <- data.frame(
  field = character(),
  metabolite = character(),
  mean_fh = numeric(),
  sd_fh = numeric(),
  mean_nonfh = numeric(),
  sd_nonfh = numeric(),
  cohens_d = numeric(),
  pvalue = numeric(),
  n_available = integer(),
  stringsAsFactors = FALSE
)

for (col in nmr_cols) {
  vals <- hr_nmr[[col]]
  if (sum(!is.na(vals)) < 100) next  # Skip if too few values

  fh_vals <- vals[hr_nmr$is_fh_genetic == TRUE]
  nonfh_vals <- vals[hr_nmr$is_fh_genetic == FALSE]

  if (sum(!is.na(fh_vals)) < 5 || sum(!is.na(nonfh_vals)) < 50) next

  m_fh <- mean(fh_vals, na.rm = TRUE)
  s_fh <- sd(fh_vals, na.rm = TRUE)
  m_nonfh <- mean(nonfh_vals, na.rm = TRUE)
  s_nonfh <- sd(nonfh_vals, na.rm = TRUE)

  # Cohen's d
  pooled_sd <- sqrt((s_fh^2 + s_nonfh^2) / 2)
  d <- ifelse(pooled_sd > 0, (m_fh - m_nonfh) / pooled_sd, 0)

  # Welch t-test
  pval <- tryCatch(t.test(fh_vals, nonfh_vals)$p.value, error = function(e) NA)

  nmr_results <- rbind(nmr_results, data.frame(
    field = col,
    metabolite = nmr_labels[col],
    mean_fh = m_fh,
    sd_fh = s_fh,
    mean_nonfh = m_nonfh,
    sd_nonfh = s_nonfh,
    cohens_d = d,
    pvalue = pval,
    n_available = sum(!is.na(vals)),
    stringsAsFactors = FALSE
  ))
}

# Bonferroni correction
nmr_results$p_bonferroni <- p.adjust(nmr_results$pvalue, method = "bonferroni")
nmr_results$significant <- nmr_results$p_bonferroni < 0.05

# Sort by absolute Cohen's d
nmr_results <- nmr_results[order(-abs(nmr_results$cohens_d)), ]

cat("TOP 30 DISCRIMINATING NMR METABOLITES (by |Cohen's d|):\n\n")
cat(sprintf("%-25s %-12s %-12s %-10s %-12s %-8s\n",
            "Metabolite", "FH Mean", "Non-FH Mean", "Cohen d", "P-Bonf", "Sig"))
cat(strrep("-", 85), "\n")

top30 <- head(nmr_results, 30)
for (i in seq_len(nrow(top30))) {
  r <- top30[i, ]
  pstr <- ifelse(is.na(r$p_bonferroni), "NA",
           ifelse(r$p_bonferroni < 1e-10, sprintf("%.1e", r$p_bonferroni),
           ifelse(r$p_bonferroni < 0.001, sprintf("%.1e", r$p_bonferroni),
                  sprintf("%.4f", r$p_bonferroni))))
  sig <- ifelse(r$significant, "***", "")
  cat(sprintf("%-25s %-12.3f %-12.3f %-10.3f %-12s %-8s\n",
              r$metabolite, r$mean_fh, r$mean_nonfh, r$cohens_d, pstr, sig))
}

n_sig <- sum(nmr_results$significant, na.rm = TRUE)
cat(sprintf("\nSignificant after Bonferroni: %d of %d metabolites\n\n", n_sig, nrow(nmr_results)))


# ==============================================================================
# B. LIPOPROTEIN SUBCLASS ANALYSIS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS B: Lipoprotein Particle Subclass Profiles\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Identify VLDL, LDL, HDL subclass columns
vldl_cols <- nmr_cols[grepl("VLDL", nmr_labels[nmr_cols])]
ldl_cols  <- nmr_cols[grepl("^[LMS]_LDL|IDL", nmr_labels[nmr_cols])]
hdl_cols  <- nmr_cols[grepl("HDL", nmr_labels[nmr_cols])]

cat("VLDL subclass columns found:", length(vldl_cols), "\n")
cat("LDL/IDL subclass columns found:", length(ldl_cols), "\n")
cat("HDL subclass columns found:", length(hdl_cols), "\n\n")

# For each lipoprotein class, compare FH vs non-FH
compare_subclass <- function(cols, class_name) {
  if (length(cols) == 0) {
    cat("  No", class_name, "columns available.\n")
    return(NULL)
  }

  cat(class_name, "subclass comparison (FH vs Non-FH, high-risk cohort):\n")
  cat(sprintf("  %-25s %-12s %-12s %-10s %-12s\n",
              "Subclass", "FH Mean", "Non-FH Mean", "Ratio", "P-value"))
  cat("  ", strrep("-", 75), "\n")

  res <- data.frame()
  for (col in cols) {
    fh <- hr_nmr[[col]][hr_nmr$is_fh_genetic == TRUE]
    nonfh <- hr_nmr[[col]][hr_nmr$is_fh_genetic == FALSE]
    if (sum(!is.na(fh)) < 5 || sum(!is.na(nonfh)) < 50) next

    m_fh <- mean(fh, na.rm = TRUE)
    m_nonfh <- mean(nonfh, na.rm = TRUE)
    ratio <- ifelse(m_nonfh != 0, m_fh / m_nonfh, NA)
    pval <- tryCatch(t.test(fh, nonfh)$p.value, error = function(e) NA)
    pstr <- ifelse(is.na(pval), "NA",
             ifelse(pval < 0.001, sprintf("%.1e", pval), sprintf("%.3f", pval)))

    cat(sprintf("  %-25s %-12.4f %-12.4f %-10.2f %-12s\n",
                nmr_labels[col], m_fh, m_nonfh, ratio, pstr))
    res <- rbind(res, data.frame(
      subclass = nmr_labels[col], fh = m_fh, nonfh = m_nonfh, ratio = ratio, p = pval))
  }
  cat("\n")
  return(res)
}

vldl_res <- compare_subclass(vldl_cols, "VLDL")
ldl_res  <- compare_subclass(ldl_cols, "LDL/IDL")
hdl_res  <- compare_subclass(hdl_cols, "HDL")

cat("KEY HYPOTHESIS: FH should show ELEVATED LDL subclass measures but\n")
cat("NORMAL/LOW VLDL measures. This VLDL-LDL discordance is the particle-level\n")
cat("mechanism behind the Trig Filter's diagnostic power.\n\n")


# ==============================================================================
# C. PARTICLE SIZE ANALYSIS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS C: Lipoprotein Particle Sizes\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

size_cols <- nmr_cols[grepl("size|_D$", nmr_labels[nmr_cols], ignore.case = TRUE)]

# Also check for columns with "size" in the raw field mapping
if (length(size_cols) == 0) {
  # Try broader search
  size_candidates <- nmr_cols[grepl("p2348[4-6]", nmr_cols)]
  if (length(size_candidates) > 0) size_cols <- size_candidates
}

if (length(size_cols) > 0) {
  cat("Particle size comparison (FH vs Non-FH):\n")
  for (col in size_cols) {
    fh <- hr_nmr[[col]][hr_nmr$is_fh_genetic == TRUE]
    nonfh <- hr_nmr[[col]][hr_nmr$is_fh_genetic == FALSE]
    if (sum(!is.na(fh)) < 5) next

    tt <- tryCatch(t.test(fh, nonfh), error = function(e) NULL)
    if (!is.null(tt)) {
      d <- (mean(fh, na.rm = TRUE) - mean(nonfh, na.rm = TRUE)) /
           sqrt((var(fh, na.rm = TRUE) + var(nonfh, na.rm = TRUE)) / 2)
      cat(sprintf("  %-20s FH: %.2f +/- %.2f  |  Non-FH: %.2f +/- %.2f  |  d = %+.3f  p = %.1e\n",
                  nmr_labels[col],
                  mean(fh, na.rm = TRUE), sd(fh, na.rm = TRUE),
                  mean(nonfh, na.rm = TRUE), sd(nonfh, na.rm = TRUE),
                  d, tt$p.value))
    }
  }
  cat("\nKEY HYPOTHESIS: FH may show LARGER LDL particles (cholesterol-enriched\n")
  cat("due to prolonged circulation from impaired LDLR clearance).\n\n")
} else {
  cat("Particle size fields not found in this NMR extraction.\n\n")
}


# ==============================================================================
# D. ApoB CORRELATIONS: FH vs NON-FH
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS D: ApoB-Particle Correlations (Novel Finding)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("CONCEPT: In FH, LDL particles spend longer in circulation due to\n")
cat("impaired LDLR clearance. Each particle accumulates more cholesterol.\n")
cat("This means: FEWER particles (lower ApoB/LDL) but HIGHER cholesterol\n")
cat("per particle. This creates a DISCORDANT ApoB-LDL relationship unique to FH.\n\n")

# ApoB vs LDL correlation in FH vs non-FH (using biochemical ApoB)
if (sum(!is.na(hr_nmr$ApoB)) > 100) {
  fh_data <- hr_nmr[hr_nmr$is_fh_genetic == TRUE & !is.na(hr_nmr$ApoB) & !is.na(hr_nmr$LDL_RW), ]
  nonfh_data <- hr_nmr[hr_nmr$is_fh_genetic == FALSE & !is.na(hr_nmr$ApoB) & !is.na(hr_nmr$LDL_RW), ]

  cor_fh <- cor.test(fh_data$ApoB, fh_data$LDL_RW)
  cor_nonfh <- cor.test(nonfh_data$ApoB, nonfh_data$LDL_RW)

  cat("1. ApoB vs LDL-C Correlation:\n")
  cat(sprintf("   FH:     r = %.3f [%.3f, %.3f]  p = %.1e  N = %d\n",
              cor_fh$estimate, cor_fh$conf.int[1], cor_fh$conf.int[2],
              cor_fh$p.value, nrow(fh_data)))
  cat(sprintf("   Non-FH: r = %.3f [%.3f, %.3f]  p = %.1e  N = %d\n",
              cor_nonfh$estimate, cor_nonfh$conf.int[1], cor_nonfh$conf.int[2],
              cor_nonfh$p.value, nrow(nonfh_data)))

  # Fisher z-test to compare correlations
  z_fh <- atanh(cor_fh$estimate)
  z_nonfh <- atanh(cor_nonfh$estimate)
  se_diff <- sqrt(1/(nrow(fh_data) - 3) + 1/(nrow(nonfh_data) - 3))
  z_test <- (z_fh - z_nonfh) / se_diff
  p_diff <- 2 * pnorm(-abs(z_test))
  cat(sprintf("   Difference test (Fisher z): z = %.2f, p = %.3f\n", z_test, p_diff))
  cat("\n")

  # ApoB/LDL ratio comparison
  cat("2. ApoB/LDL Ratio (particle density proxy):\n")
  apob_ldl_fh <- fh_data$ApoB / fh_data$LDL_RW
  apob_ldl_nonfh <- nonfh_data$ApoB / nonfh_data$LDL_RW
  tt_ratio <- t.test(apob_ldl_fh, apob_ldl_nonfh)
  d_ratio <- (mean(apob_ldl_fh, na.rm = TRUE) - mean(apob_ldl_nonfh, na.rm = TRUE)) /
             sqrt((var(apob_ldl_fh, na.rm = TRUE) + var(apob_ldl_nonfh, na.rm = TRUE)) / 2)
  cat(sprintf("   FH:     %.4f +/- %.4f\n", mean(apob_ldl_fh, na.rm = TRUE), sd(apob_ldl_fh, na.rm = TRUE)))
  cat(sprintf("   Non-FH: %.4f +/- %.4f\n", mean(apob_ldl_nonfh, na.rm = TRUE), sd(apob_ldl_nonfh, na.rm = TRUE)))
  cat(sprintf("   Cohen's d = %.3f, p = %.1e\n", d_ratio, tt_ratio$p.value))
  cat("\n")

  cat("KEY FINDING: If FH has LOWER ApoB/LDL ratio, it means each LDL\n")
  cat("particle carries MORE cholesterol. This is the 'cholesterol-enriched\n")
  cat("particle' signature of impaired LDLR clearance.\n\n")
}

# --- ApoB correlations with NMR particle measures ---
# Find LDL particle concentration NMR columns
ldl_p_cols <- nmr_cols[grepl("LDL_P", nmr_labels[nmr_cols])]

if (length(ldl_p_cols) > 0 && sum(!is.na(hr_nmr$ApoB)) > 100) {
  cat("3. ApoB vs NMR LDL Particle Concentrations:\n")
  for (col in ldl_p_cols) {
    fh_sub <- hr_nmr[hr_nmr$is_fh_genetic == TRUE & !is.na(hr_nmr$ApoB) & !is.na(hr_nmr[[col]]), ]
    nonfh_sub <- hr_nmr[hr_nmr$is_fh_genetic == FALSE & !is.na(hr_nmr$ApoB) & !is.na(hr_nmr[[col]]), ]

    if (nrow(fh_sub) >= 10 && nrow(nonfh_sub) >= 50) {
      r_fh <- cor(fh_sub$ApoB, fh_sub[[col]], use = "complete.obs")
      r_nonfh <- cor(nonfh_sub$ApoB, nonfh_sub[[col]], use = "complete.obs")
      cat(sprintf("   %-20s  FH r = %.3f  |  Non-FH r = %.3f  |  diff = %+.3f\n",
                  nmr_labels[col], r_fh, r_nonfh, r_fh - r_nonfh))
    }
  }
  cat("\n")
}


# ==============================================================================
# E. REMNANT CHOLESTEROL ANALYSIS
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS E: Remnant Cholesterol in FH vs Non-FH\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Find remnant cholesterol column
rem_col <- nmr_cols[grepl("Remnant|remnant", nmr_labels[nmr_cols])]

if (length(rem_col) > 0) {
  for (rc in rem_col) {
    fh <- hr_nmr[[rc]][hr_nmr$is_fh_genetic == TRUE]
    nonfh <- hr_nmr[[rc]][hr_nmr$is_fh_genetic == FALSE]

    if (sum(!is.na(fh)) >= 5) {
      tt <- t.test(fh, nonfh)
      d <- (mean(fh, na.rm = TRUE) - mean(nonfh, na.rm = TRUE)) /
           sqrt((var(fh, na.rm = TRUE) + var(nonfh, na.rm = TRUE)) / 2)
      cat(sprintf("%-20s  FH: %.3f +/- %.3f  |  Non-FH: %.3f +/- %.3f  |  d = %+.3f  p = %.1e\n",
                  nmr_labels[rc],
                  mean(fh, na.rm = TRUE), sd(fh, na.rm = TRUE),
                  mean(nonfh, na.rm = TRUE), sd(nonfh, na.rm = TRUE),
                  d, tt$p.value))
    }
  }
  cat("\nKEY FINDING: FH should have LOWER/NORMAL remnant cholesterol.\n")
  cat("High remnant cholesterol = VLDL/IDL pathway = metabolic syndrome.\n")
  cat("FH is LDLR-driven, not VLDL-driven.\n\n")
} else {
  # Calculate remnant cholesterol from standard lipids if NMR not available
  cat("NMR remnant cholesterol not found. Calculating from standard lipids:\n")
  cat("Remnant-C = Total Chol - LDL - HDL\n")
  hr_nmr$Remnant_calc <- hr_nmr$CHOL - hr_nmr$LDL_RW - hr_nmr$HDL.1

  fh <- hr_nmr$Remnant_calc[hr_nmr$is_fh_genetic == TRUE]
  nonfh <- hr_nmr$Remnant_calc[hr_nmr$is_fh_genetic == FALSE]
  tt <- t.test(fh, nonfh)
  d <- (mean(fh, na.rm = TRUE) - mean(nonfh, na.rm = TRUE)) /
       sqrt((var(fh, na.rm = TRUE) + var(nonfh, na.rm = TRUE)) / 2)
  cat(sprintf("  FH: %.3f +/- %.3f  |  Non-FH: %.3f +/- %.3f  |  d = %+.3f  p = %.1e\n",
              mean(fh, na.rm = TRUE), sd(fh, na.rm = TRUE),
              mean(nonfh, na.rm = TRUE), sd(nonfh, na.rm = TRUE),
              d, tt$p.value))
  cat("\n")
}


# ==============================================================================
# F. NMR-BASED FH PREDICTION (EXPLORATORY)
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS F: Can NMR Metabolites Predict Genetic FH?\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Select NMR columns with sufficient coverage
nmr_for_model <- nmr_cols[sapply(nmr_cols, function(c) {
  sum(!is.na(hr_nmr[[c]])) > 0.5 * nrow(hr_nmr)
})]
cat("NMR features with >50% coverage:", length(nmr_for_model), "\n")

if (length(nmr_for_model) >= 5) {
  # Prepare complete-case matrix
  model_data <- hr_nmr[, c("is_fh_genetic", nmr_for_model), with = FALSE]
  model_data <- model_data[complete.cases(model_data), ]
  cat("Complete cases for NMR model:", nrow(model_data), "\n")
  cat("  FH:", sum(model_data$is_fh_genetic), "| Non-FH:", sum(!model_data$is_fh_genetic), "\n")

  if (sum(model_data$is_fh_genetic) >= 20) {
    # Simple logistic with top 10 NMR features (by univariate association)
    top_nmr <- head(nmr_results$field[nmr_results$field %in% nmr_for_model], 10)

    if (length(top_nmr) >= 3) {
      # Standardise
      for (col in top_nmr) {
        model_data[[paste0(col, "_z")]] <- scale(model_data[[col]])[,1]
      }
      z_cols <- paste0(top_nmr, "_z")

      formula_str <- paste("is_fh_genetic ~", paste(z_cols, collapse = " + "))
      nmr_glm <- glm(as.formula(formula_str), data = model_data, family = binomial)

      # AUC
      nmr_pred <- predict(nmr_glm, type = "response")
      roc_nmr <- pROC::roc(model_data$is_fh_genetic, nmr_pred, quiet = TRUE)
      ci_nmr <- pROC::ci.auc(roc_nmr)
      cat(sprintf("\nTop-%d NMR features AUC: %.3f [%.3f - %.3f]\n",
                  length(top_nmr), as.numeric(roc_nmr$auc), ci_nmr[1], ci_nmr[3]))

      # Compare to TUDOR
      tudor_in_model <- hr_nmr$tudor_prob[complete.cases(hr_nmr[, c("is_fh_genetic", nmr_for_model), with = FALSE])]
      roc_tudor_nmr <- pROC::roc(model_data$is_fh_genetic, tudor_in_model, quiet = TRUE)
      cat(sprintf("TUDOR AUC (same subset):  %.3f\n", as.numeric(roc_tudor_nmr$auc)))

      # TUDOR + NMR combined
      model_data$tudor_p <- tudor_in_model
      formula_combined <- paste("is_fh_genetic ~ tudor_p +", paste(z_cols, collapse = " + "))
      combined_glm <- glm(as.formula(formula_combined), data = model_data, family = binomial)
      combined_pred <- predict(combined_glm, type = "response")
      roc_combined <- pROC::roc(model_data$is_fh_genetic, combined_pred, quiet = TRUE)
      ci_combined <- pROC::ci.auc(roc_combined)
      cat(sprintf("TUDOR + NMR combined:     %.3f [%.3f - %.3f]\n",
                  as.numeric(roc_combined$auc), ci_combined[1], ci_combined[3]))

      dt_augment <- pROC::roc.test(roc_tudor_nmr, roc_combined, method = "delong")
      cat(sprintf("DeLong p (TUDOR vs TUDOR+NMR): %.4f\n", dt_augment$p.value))
      cat("\n")

      # Feature importance (absolute z-scores from logistic)
      cat("NMR feature coefficients in combined model:\n")
      cf <- summary(combined_glm)$coefficients
      nmr_coefs <- cf[z_cols, , drop = FALSE]
      nmr_coefs <- nmr_coefs[order(-abs(nmr_coefs[, "Estimate"])), ]
      for (i in seq_len(nrow(nmr_coefs))) {
        feat <- rownames(nmr_coefs)[i]
        orig_col <- gsub("_z$", "", feat)
        label <- nmr_labels[orig_col]
        cat(sprintf("  %-25s beta = %+.4f  OR = %.3f  p = %.1e\n",
                    label, nmr_coefs[i, "Estimate"],
                    exp(nmr_coefs[i, "Estimate"]), nmr_coefs[i, "Pr(>|z|)"]))
      }
    }
  } else {
    cat("Too few FH cases with complete NMR data for prediction model.\n")
  }
} else {
  cat("Insufficient NMR features for prediction analysis.\n")
}

cat("\n")


# ==============================================================================
# G. NMR VALIDATION OF TRIG FILTER MECHANISM
# ==============================================================================
cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS G: NMR-Level Validation of the Trig Filter\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("CONCEPT: The Trig Filter = LDL / (TRG + 0.1) captures the ratio of\n")
cat("LDL-pathway to VLDL-pathway lipids. At the NMR level, this should\n")
cat("correlate with LDL particle concentration / VLDL particle concentration.\n\n")

# Find total VLDL and LDL particle columns
total_vldl_cols <- nmr_cols[grepl("VLDL_P|VLDL_lipids", nmr_labels[nmr_cols])]
total_ldl_cols  <- nmr_cols[grepl("[LMS]_LDL_P|[LMS]_LDL_lipids|IDL_P|IDL_lipids",
                                   nmr_labels[nmr_cols])]

# Correlation between Trig Filter and NMR VLDL measures
if (length(total_vldl_cols) > 0) {
  cat("Trig Filter correlations with NMR VLDL measures:\n")
  for (col in total_vldl_cols) {
    r <- cor(hr_nmr$Trig_Filter_RW, hr_nmr[[col]], use = "complete.obs")
    if (!is.na(r)) {
      cat(sprintf("  %-25s r = %+.3f (expected: NEGATIVE)\n", nmr_labels[col], r))
    }
  }
  cat("\n")
}

if (length(total_ldl_cols) > 0) {
  cat("Trig Filter correlations with NMR LDL measures:\n")
  for (col in total_ldl_cols) {
    r <- cor(hr_nmr$Trig_Filter_RW, hr_nmr[[col]], use = "complete.obs")
    if (!is.na(r)) {
      cat(sprintf("  %-25s r = %+.3f (expected: POSITIVE)\n", nmr_labels[col], r))
    }
  }
  cat("\n")
}

cat("KEY FINDING: If Trig Filter correlates NEGATIVELY with VLDL particles\n")
cat("and POSITIVELY with LDL particles, this confirms it captures the\n")
cat("VLDL-LDL discordance at the particle level — the fundamental mechanism\n")
cat("distinguishing monogenic FH from polygenic/metabolic hyperlipidaemia.\n\n")


# ==============================================================================
# H. VOLCANO PLOT
# ==============================================================================
if (nrow(nmr_results) > 5) {
  cat("Generating NMR volcano plot...\n")

  plot_data <- nmr_results[!is.na(nmr_results$pvalue) & nmr_results$pvalue > 0, ]
  plot_data$neg_log10_p <- -log10(plot_data$pvalue)
  plot_data$label_show <- ifelse(abs(plot_data$cohens_d) > 0.1 | plot_data$significant,
                                  plot_data$metabolite, "")

  # Classify by lipoprotein type
  plot_data$lipo_class <- "Other"
  plot_data$lipo_class[grepl("VLDL", plot_data$metabolite)] <- "VLDL"
  plot_data$lipo_class[grepl("LDL|IDL", plot_data$metabolite)] <- "LDL/IDL"
  plot_data$lipo_class[grepl("HDL", plot_data$metabolite)] <- "HDL"

  p_volcano <- ggplot(plot_data, aes(x = cohens_d, y = neg_log10_p, color = lipo_class)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = -log10(0.05 / nrow(nmr_results)),
               linetype = "dashed", color = "red", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey50") +
    scale_color_manual(values = c("VLDL" = "#E69F00", "LDL/IDL" = "#D55E00",
                                   "HDL" = "#0072B2", "Other" = "grey50")) +
    labs(x = "Cohen's d (FH - Non-FH)",
         y = "-log10(P-value)",
         title = "NMR Metabolomics: FH vs Non-FH",
         subtitle = "Bonferroni threshold shown (dashed red line)",
         color = "Lipoprotein Class") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")

  ggsave(file.path(PLOT_DIR, "fig_nmr_volcano.pdf"), p_volcano, width = 10, height = 7)
  ggsave(file.path(PLOT_DIR, "fig_nmr_volcano.png"), p_volcano, width = 10, height = 7, dpi = 300)
  cat("  Saved: fig_nmr_volcano.pdf/png\n\n")
}


# ==============================================================================
# I. SAVE RESULTS
# ==============================================================================
results_09 <- list(
  nmr_results = nmr_results,
  vldl_comparison = vldl_res,
  ldl_comparison = ldl_res,
  hdl_comparison = hdl_res,
  nmr_labels = nmr_labels,
  n_with_nmr = sum(nmr_coverage),
  n_fh_with_nmr = sum(df_nmr$is_fh_genetic),
  timestamp = Sys.time()
)
saveRDS(results_09, file.path(OUTPUT_DIR, "results_09_nmr.rds"))

cat("\n=== 09_nmr_metabolomics.R COMPLETE ===\n")
cat("Results saved to:", file.path(OUTPUT_DIR, "results_09_nmr.rds"), "\n")
