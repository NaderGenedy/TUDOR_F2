# ==============================================================================
# TUDOR PIPELINE: STEP 25 — VUS RECLASSIFICATION & GENETIC FRAMEWORK
# ==============================================================================
# PURPOSE: Comprehensive Variant of Uncertain Significance (VUS) analysis
#          for Nature submission. Addresses how TUDOR probability can guide
#          VUS interpretation in FH genetic testing.
#
# SECTIONS:
#   1. VUS Classification Framework
#   2. TUDOR-Guided VUS Reclassification
#   3. Genotype-Phenotype Correlations by Variant Class
#   4. ACMG/AMP Criteria Integration
#   5. CNV Detection Limitations
#   6. Polygenic Score Comparison
#
# AUTHORS: Tudor Pipeline Team
# INSTITUTION: Cardiff and Vale University Health Board
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
FIG_DIR    <- file.path(OUTPUT_DIR, "figures")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("TUDOR PIPELINE: 25 — VUS RECLASSIFICATION & GENETIC FRAMEWORK\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# SECTION 1: VUS CLASSIFICATION FRAMEWORK
# ==============================================================================
cat("================================================================\n")
cat("SECTION 1: VUS CLASSIFICATION FRAMEWORK\n")
cat("================================================================\n\n")

# ClinVar classification levels (ordered by pathogenicity)
CLINVAR_CLASSES <- c(
  "Pathogenic",
  "Likely_Pathogenic",
  "VUS",                     # Variant of Uncertain Significance
  "Likely_Benign",
  "Benign",
  "Conflicting"              # Conflicting interpretations
)

cat("ClinVar Classification System for FH Variants:\n")
cat("  1. Pathogenic (P)        — Disease-causing, actionable\n")
cat("  2. Likely Pathogenic (LP) — >90% probability pathogenic, actionable\n")
cat("  3. VUS                    — Insufficient evidence, NOT actionable (current practice)\n")
cat("  4. Likely Benign (LB)    — >90% probability benign\n")
cat("  5. Benign (B)            — Not disease-causing\n")
cat("  6. Conflicting           — Discordant lab interpretations\n\n")

cat("PROBLEM: VUS in FH genes (LDLR, APOB, PCSK9) leave patients in\n")
cat("diagnostic limbo — no cascade screening, uncertain treatment intensity.\n")
cat("TUDOR can provide phenotypic evidence to support VUS reclassification.\n\n")

# ==============================================================================
# SECTION 2: TUDOR-GUIDED VUS RECLASSIFICATION MODEL
# ==============================================================================
cat("================================================================\n")
cat("SECTION 2: TUDOR-GUIDED VUS RECLASSIFICATION\n")
cat("================================================================\n\n")

# Load data
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
use_simulated <- FALSE

if (file.exists(rds_file)) {
  df <- readRDS(rds_file)
  setDT(df)
  if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
    setnames(df, "participant.eid", "eid")
  }
  cat("Loaded real data:", nrow(df), "participants\n")
} else {
  cat("Real data not available. Using simulated data based on known parameters.\n")
  use_simulated <- TRUE

  n_total <- 400000
  n_fh_pathogenic <- 800       # ~1/500 prevalence
  n_vus <- 2400                # ~3x pathogenic (typical VUS rate)
  n_negative <- n_total - n_fh_pathogenic - n_vus

  df <- data.table(
    eid = 1:n_total,
    variant_class = c(
      rep("Pathogenic", n_fh_pathogenic),
      rep("VUS", n_vus),
      rep("Negative", n_negative)
    ),
    is_fh_genetic = c(
      rep(1, n_fh_pathogenic),
      rbinom(n_vus, 1, 0.35),   # ~35% of VUS are truly pathogenic
      rep(0, n_negative)
    )
  )

  # Simulate TUDOR scores
  df[variant_class == "Pathogenic", tudor_prob := pmin(rbeta(.N, 3, 12), 1)]
  df[variant_class == "VUS" & is_fh_genetic == 1, tudor_prob := pmin(rbeta(.N, 2.5, 15), 1)]
  df[variant_class == "VUS" & is_fh_genetic == 0, tudor_prob := pmin(rbeta(.N, 1.5, 50), 1)]
  df[variant_class == "Negative", tudor_prob := pmin(rbeta(.N, 1, 100), 1)]

  # Simulate lipid values
  df[, LDL_RW := rnorm(.N, ifelse(is_fh_genetic == 1, 6.5, 4.2), 1.2)]
  df[, TRG.1 := rlnorm(.N, log(ifelse(is_fh_genetic == 1, 1.1, 1.8)), 0.5)]
  df[, HDL.1 := rnorm(.N, 1.4, 0.35)]
  df[, Age_at_LDL1 := rnorm(.N, 57, 8)]
  df[, Gender_num := rbinom(.N, 1, 0.46)]
  df[, Trig_Filter_RW := LDL_RW / (TRG.1 + 0.1)]
  df[, cohort_high_risk := LDL_RW > 4.9]
  df[, edlcn_score := ifelse(LDL_RW >= 8.5, 8,
                      ifelse(LDL_RW >= 6.5, 5,
                      ifelse(LDL_RW >= 5.0, 3,
                      ifelse(LDL_RW >= 4.0, 1, 0))))]

  # Simulate gene assignment for pathogenic/VUS
  df[variant_class == "Pathogenic", gene := sample(c("LDLR", "APOB", "PCSK9"),
    .N, replace = TRUE, prob = c(0.75, 0.20, 0.05))]
  df[variant_class == "VUS", gene := sample(c("LDLR", "APOB", "PCSK9"),
    .N, replace = TRUE, prob = c(0.60, 0.25, 0.15))]
  df[variant_class == "Negative", gene := NA_character_]

  # Simulate variant types
  df[variant_class %in% c("Pathogenic", "VUS"), variant_type := sample(
    c("Missense", "Frameshift", "Splice", "Nonsense", "In-frame_del"),
    .N, replace = TRUE, prob = c(0.55, 0.15, 0.12, 0.10, 0.08))]
}

# If real data but no variant_class column, create from existing genetic info
if (!use_simulated && !"variant_class" %in% names(df)) {
  cat("Deriving variant classification from existing genetic data...\n")
  df[, variant_class := ifelse(is_fh_genetic == 1, "Pathogenic", "Negative")]
  # Note: True VUS data requires ClinVar annotation pipeline
  cat("  WARNING: VUS not available in current dataset.\n")
  cat("  Analysis will use Pathogenic vs Negative comparison.\n")
  cat("  For full VUS analysis, merge ClinVar VUS annotations.\n\n")
}

# ==============================================================================
# 2a. TUDOR Score Distribution by Variant Class
# ==============================================================================
cat("--- TUDOR Score Distribution by Variant Class ---\n\n")

for (vc in unique(df$variant_class[!is.na(df$variant_class)])) {
  sub <- df[variant_class == vc]
  cat(sprintf("  %-15s: N = %6d | TUDOR prob: median = %.4f, IQR [%.4f - %.4f]\n",
              vc, nrow(sub),
              median(sub$tudor_prob, na.rm = TRUE),
              quantile(sub$tudor_prob, 0.25, na.rm = TRUE),
              quantile(sub$tudor_prob, 0.75, na.rm = TRUE)))
}

# Statistical test: Kruskal-Wallis across groups
if (length(unique(df$variant_class)) >= 3) {
  kw_test <- kruskal.test(tudor_prob ~ variant_class,
                           data = df[!is.na(variant_class) & !is.na(tudor_prob)])
  cat(sprintf("\nKruskal-Wallis test: chi2 = %.1f, df = %d, p = %.2e\n",
              kw_test$statistic, kw_test$parameter, kw_test$p.value))
}

# Pairwise Wilcoxon tests with Bonferroni correction
if ("VUS" %in% df$variant_class) {
  cat("\nPairwise Wilcoxon tests (Bonferroni-corrected):\n")

  pairs <- list(
    c("VUS", "Pathogenic"),
    c("VUS", "Negative"),
    c("Pathogenic", "Negative")
  )

  for (pair in pairs) {
    g1 <- df[variant_class == pair[1]]$tudor_prob
    g2 <- df[variant_class == pair[2]]$tudor_prob
    if (length(g1) > 0 && length(g2) > 0) {
      wt <- wilcox.test(g1, g2)
      p_adj <- min(wt$p.value * length(pairs), 1)
      cat(sprintf("  %-12s vs %-12s: p_raw = %.2e, p_Bonf = %.2e\n",
                  pair[1], pair[2], wt$p.value, p_adj))
    }
  }
}
cat("\n")

# ==============================================================================
# 2b. VUS Reclassification Using TUDOR Thresholds
# ==============================================================================
if ("VUS" %in% df$variant_class) {
  cat("--- VUS Reclassification Using TUDOR Probability ---\n\n")

  vus <- df[variant_class == "VUS"]
  n_vus <- nrow(vus)

  # Determine TUDOR-based reclassification thresholds
  # Use pathogenic group's 25th percentile as "likely pathogenic" threshold
  # Use negative group's 75th percentile as "likely benign" threshold
  patho_q25 <- quantile(df[variant_class == "Pathogenic"]$tudor_prob, 0.25, na.rm = TRUE)
  neg_q75 <- quantile(df[variant_class == "Negative"]$tudor_prob, 0.75, na.rm = TRUE)

  cat(sprintf("Reclassification thresholds:\n"))
  cat(sprintf("  Likely Pathogenic if TUDOR prob >= %.4f (P25 of pathogenic group)\n", patho_q25))
  cat(sprintf("  Likely Benign if TUDOR prob <= %.4f (P75 of negative group)\n", neg_q75))
  cat(sprintf("  Remains VUS if TUDOR prob between %.4f and %.4f\n\n", neg_q75, patho_q25))

  # Reclassify
  vus[, tudor_reclassification := ifelse(
    tudor_prob >= patho_q25, "Likely_Pathogenic",
    ifelse(tudor_prob <= neg_q75, "Likely_Benign", "Remains_VUS")
  )]

  reclass_table <- table(
    TUDOR_Class = vus$tudor_reclassification,
    True_FH = ifelse(vus$is_fh_genetic == 1, "True_FH", "Not_FH")
  )

  cat("VUS Reclassification Table:\n")
  print(reclass_table)

  # Performance metrics
  cat(sprintf("\nReclassification Summary:\n"))
  for (rc in c("Likely_Pathogenic", "Remains_VUS", "Likely_Benign")) {
    sub_rc <- vus[tudor_reclassification == rc]
    n_rc <- nrow(sub_rc)
    n_fh <- sum(sub_rc$is_fh_genetic)
    pct <- 100 * n_rc / n_vus
    fh_rate <- ifelse(n_rc > 0, 100 * n_fh / n_rc, 0)
    cat(sprintf("  %-20s: N = %4d (%.1f%%) | FH rate = %.1f%%\n",
                rc, n_rc, pct, fh_rate))
  }

  # Correctly reclassified
  correctly_patho <- sum(vus$tudor_reclassification == "Likely_Pathogenic" &
                          vus$is_fh_genetic == 1)
  correctly_benign <- sum(vus$tudor_reclassification == "Likely_Benign" &
                           vus$is_fh_genetic == 0)
  total_correct <- correctly_patho + correctly_benign
  total_reclassified <- sum(vus$tudor_reclassification != "Remains_VUS")

  cat(sprintf("\n  Total reclassified: %d / %d (%.1f%%)\n",
              total_reclassified, n_vus, 100 * total_reclassified / n_vus))
  cat(sprintf("  Correctly reclassified: %d / %d (%.1f%% accuracy)\n",
              total_correct, total_reclassified,
              ifelse(total_reclassified > 0, 100 * total_correct / total_reclassified, 0)))

  # Save VUS reclassification table
  vus_summary <- vus[, .(
    N = .N,
    N_FH = sum(is_fh_genetic),
    FH_rate = mean(is_fh_genetic),
    Median_TUDOR = median(tudor_prob, na.rm = TRUE),
    Median_LDL = median(LDL_RW, na.rm = TRUE),
    Median_TrigFilter = median(Trig_Filter_RW, na.rm = TRUE)
  ), by = tudor_reclassification]

  fwrite(vus_summary, file.path(TABLE_DIR, "vus_reclassification_summary.csv"))
  cat("\n  Saved: vus_reclassification_summary.csv\n")
}
cat("\n")

# ==============================================================================
# SECTION 3: GENOTYPE-PHENOTYPE CORRELATIONS BY VARIANT CLASS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 3: GENOTYPE-PHENOTYPE CORRELATIONS\n")
cat("================================================================\n\n")

if ("gene" %in% names(df)) {
  gene_data <- df[!is.na(gene) & gene %in% c("LDLR", "APOB", "PCSK9")]

  if (nrow(gene_data) > 0) {
    cat("--- Lipid Profile by Gene ---\n\n")
    cat(sprintf("%-8s | %5s | %8s | %8s | %8s | %8s | %10s\n",
                "Gene", "N", "LDL_RW", "TRG", "HDL", "TUDOR_p", "TrigFilter"))
    cat(strrep("-", 75), "\n")

    for (g in c("LDLR", "APOB", "PCSK9")) {
      sub <- gene_data[gene == g]
      if (nrow(sub) > 0) {
        cat(sprintf("%-8s | %5d | %8.1f | %8.2f | %8.2f | %8.4f | %10.1f\n",
                    g, nrow(sub),
                    median(sub$LDL_RW, na.rm = TRUE),
                    median(sub$TRG.1, na.rm = TRUE),
                    median(sub$HDL.1, na.rm = TRUE),
                    median(sub$tudor_prob, na.rm = TRUE),
                    median(sub$Trig_Filter_RW, na.rm = TRUE)))
      }
    }
    cat("\n")

    # Gene-specific TUDOR AUC (within variant carriers)
    cat("--- Gene-Specific TUDOR AUC ---\n")
    for (g in c("LDLR", "APOB", "PCSK9")) {
      sub <- df[(!is.na(gene) & gene == g) | variant_class == "Negative"]
      sub[, is_gene_fh := as.integer(!is.na(gene) & gene == g)]

      if (sum(sub$is_gene_fh) >= 10) {
        r <- roc(sub$is_gene_fh, sub$tudor_prob, quiet = TRUE)
        ci <- ci.auc(r, method = "delong")
        cat(sprintf("  %-8s: AUC = %.3f [%.3f - %.3f] (N_pos=%d, N_neg=%d)\n",
                    g, ci[2], ci[1], ci[3],
                    sum(sub$is_gene_fh), sum(!sub$is_gene_fh)))
      }
    }
    cat("\n")

    # Kruskal-Wallis for LDL differences across genes
    if (nrow(gene_data) > 30) {
      kw_ldl <- kruskal.test(LDL_RW ~ gene, data = gene_data)
      cat(sprintf("Kruskal-Wallis LDL across genes: p = %.2e\n", kw_ldl$p.value))

      kw_trig <- kruskal.test(TRG.1 ~ gene, data = gene_data)
      cat(sprintf("Kruskal-Wallis TRG across genes: p = %.2e\n\n", kw_trig$p.value))
    }
  }
}

# ==============================================================================
# SECTION 4: VARIANT TYPE ANALYSIS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 4: VARIANT TYPE ANALYSIS\n")
cat("================================================================\n\n")

if ("variant_type" %in% names(df)) {
  vt_data <- df[!is.na(variant_type)]

  cat("--- TUDOR Performance by Variant Type ---\n\n")
  cat(sprintf("%-15s | %5s | %8s | %8s | %8s\n",
              "Type", "N", "LDL_RW", "TUDOR_p", "TrigFilt"))
  cat(strrep("-", 55), "\n")

  for (vt in c("Nonsense", "Frameshift", "Splice", "Missense", "In-frame_del")) {
    sub <- vt_data[variant_type == vt]
    if (nrow(sub) > 0) {
      cat(sprintf("%-15s | %5d | %8.1f | %8.4f | %8.1f\n",
                  vt, nrow(sub),
                  median(sub$LDL_RW, na.rm = TRUE),
                  median(sub$tudor_prob, na.rm = TRUE),
                  median(sub$Trig_Filter_RW, na.rm = TRUE)))
    }
  }
  cat("\n")

  cat("Expected phenotype severity: Nonsense/Frameshift > Splice > Missense\n")
  cat("(Null variants → complete LDLR loss → higher LDL than missense)\n\n")

  # Null vs Missense comparison
  vt_data[, is_null := variant_type %in% c("Nonsense", "Frameshift")]
  if (sum(vt_data$is_null) >= 5 && sum(!vt_data$is_null) >= 5) {
    wt_ldl <- wilcox.test(LDL_RW ~ is_null, data = vt_data)
    wt_tudor <- wilcox.test(tudor_prob ~ is_null, data = vt_data)
    cat(sprintf("Null vs Missense LDL: p = %.2e\n", wt_ldl$p.value))
    cat(sprintf("Null vs Missense TUDOR: p = %.2e\n\n", wt_tudor$p.value))
  }
}

# ==============================================================================
# SECTION 5: CNV DETECTION LIMITATIONS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 5: CNV DETECTION LIMITATIONS\n")
cat("================================================================\n\n")

cat("CRITICAL LIMITATION: Copy Number Variant (CNV) Detection\n\n")
cat("UK Biobank uses the Affymetrix UK BiLEVE/UK Biobank Axiom arrays\n")
cat("for genotyping. These arrays have LIMITED ability to detect:\n\n")
cat("  1. LDLR large deletions/duplications (~10% of FH mutations)\n")
cat("     - Exon deletions (e.g., exons 2-6 deletion, Afrikaner founder)\n")
cat("     - Whole-gene deletions\n")
cat("     - Large duplications causing gene disruption\n\n")
cat("  2. Complex rearrangements\n")
cat("     - Alu-mediated recombination in LDLR (known hotspot)\n")
cat("     - Gene conversions between LDLR and pseudogenes\n\n")
cat("Estimated impact on TUDOR validation:\n")
cat("  - ~10% of true FH cases may be MISSED as cases (false negatives)\n")
cat("  - These appear as 'controls' but have FH → deflates specificity\n")
cat("  - Net effect: TUDOR's true AUC is likely HIGHER than reported\n\n")

cat("RECOMMENDATION: UKB whole-exome sequencing (WES) data covers ~50k\n")
cat("participants. Cross-reference WES results for CNV detection and\n")
cat("report sensitivity analysis with/without CNV-detected cases.\n\n")

# Estimate CNV impact
cnv_fraction <- 0.10  # 10% of FH mutations are CNVs
cat(sprintf("If %d%% of FH cases are CNVs missed by array genotyping:\n", cnv_fraction * 100))
cat(sprintf("  Estimated missed cases: ~%.0f\n",
            sum(df$is_fh_genetic, na.rm = TRUE) * cnv_fraction / (1 - cnv_fraction)))
cat(sprintf("  True FH count (corrected): ~%.0f\n",
            sum(df$is_fh_genetic, na.rm = TRUE) / (1 - cnv_fraction)))
cat(sprintf("  True prevalence (corrected): ~1/%.0f\n\n",
            nrow(df) / (sum(df$is_fh_genetic, na.rm = TRUE) / (1 - cnv_fraction))))

# ==============================================================================
# SECTION 6: ACMG/AMP CRITERIA INTEGRATION
# ==============================================================================
cat("================================================================\n")
cat("SECTION 6: ACMG/AMP EVIDENCE FRAMEWORK\n")
cat("================================================================\n\n")

cat("PROPOSED INTEGRATION: TUDOR as Supporting Evidence in ACMG/AMP Framework\n\n")

cat("ACMG/AMP Evidence Categories Relevant to TUDOR:\n")
cat("  PP4 (Supporting Pathogenic): Patient phenotype specific to disease\n")
cat("       → HIGH TUDOR probability (>75th percentile of pathogenic group)\n")
cat("         provides PP4-level evidence for VUS pathogenicity\n\n")
cat("  BS4 (Supporting Benign): Lack of segregation in affected family\n")
cat("       → LOW TUDOR probability (<25th percentile of pathogenic group)\n")
cat("         supports benign classification\n\n")
cat("  PP3 (Supporting Pathogenic): Computational evidence supports deleterious\n")
cat("       → TUDOR is a CLINICAL prediction model, complementary to\n")
cat("         sequence-based tools (REVEL, CADD, AlphaMissense)\n\n")

cat("PROPOSED TUDOR-ACMG INTEGRATION ALGORITHM:\n")
cat("  1. Patient has VUS in LDLR/APOB/PCSK9\n")
cat("  2. Calculate TUDOR probability from clinical data\n")
cat("  3. If TUDOR >= pathogenic P25: Apply PP4 (Supporting Pathogenic)\n")
cat("  4. If TUDOR <= negative P75:  Apply BS4-equivalent (Supporting Benign)\n")
cat("  5. If TUDOR intermediate:      No ACMG evidence, VUS remains\n")
cat("  6. Combine with other ACMG criteria for final classification\n\n")

cat("VALIDATION NEEDED:\n")
cat("  - Prospective study of VUS reclassification using TUDOR\n")
cat("  - Comparison with functional assays (saturation mutagenesis, LDLR uptake)\n")
cat("  - Expert panel review of TUDOR-reclassified variants\n\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================

vus_results <- list(
  framework = list(
    clinvar_classes = CLINVAR_CLASSES,
    cnv_fraction = cnv_fraction
  ),
  timestamp = Sys.time()
)

if ("VUS" %in% df$variant_class) {
  vus_results$reclassification <- list(
    patho_threshold = patho_q25,
    benign_threshold = neg_q75,
    table = reclass_table,
    summary = vus_summary,
    total_vus = n_vus,
    total_reclassified = total_reclassified,
    accuracy = ifelse(total_reclassified > 0, total_correct / total_reclassified, NA)
  )
}

saveRDS(vus_results, file.path(OUTPUT_DIR, "25_vus_results.rds"))
cat("Saved results to: 25_vus_results.rds\n")

cat("\n=== 25_vus_genetic_framework.R COMPLETE ===\n")
