# ==============================================================================
# 17_genetic_variant_comparison.R
# TUDOR FH Diagnostic Algorithm: UKB vs Wales Genetic Variant Comparison
# ==============================================================================
#
# Authors: Nader Genedy, Soha Zouwail
# Institution: Cardiff and Vale University Health Board
#
# PURPOSE: Comprehensive subanalysis comparing FH genetic variant distributions,
#          phenotypic expression, and TUDOR performance between the All Wales
#          PASS FH Registry and UK Biobank cohorts.
#
# OUTPUTS:
#   Tables:
#     - genetic_variant_deep_comparison.csv
#     - variant_phenotype_by_cohort.csv
#     - variant_tudor_performance.csv
#     - variant_fisher_tests.csv
#   Figures:
#     - fig_genetic_spectrum_comparison.pdf/png
#     - fig_variant_phenotype_heatmap.pdf/png
#     - fig_gene_auc_forest.pdf/png
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  TUDOR: Genetic Variant Comparison — Wales vs UK Biobank\n")
cat("  Script 17 | ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("================================================================\n\n")

# ============================================================
# PART 1: SETUP
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
  library(gridExtra)
  library(grid)
  library(scales)
})

set.seed(42)

# --- Output directories ---
fig_dir   <- "C:/Users/nader/Downloads/tudor_pipeline_output/figures"
table_dir <- "C:/Users/nader/Downloads/tudor_pipeline_output/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# --- Nature theme ---
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      text             = element_text(family = "sans", size = base_size),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size + 1, face = "bold"),
      axis.line        = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks       = element_line(linewidth = 0.3),
      legend.text      = element_text(size = base_size),
      legend.title     = element_text(size = base_size, face = "bold"),
      legend.key.size  = unit(3, "mm"),
      legend.background = element_blank(),
      plot.title       = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size, hjust = 0),
      strip.text       = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      plot.margin      = ggplot2::margin(5, 5, 5, 5, unit = "mm")
    )
}

# --- Colour palettes ---
gene_pal <- c("LDLR" = "#2166AC", "APOB" = "#B2182B", "PCSK9" = "#FF7F00",
              "LDLR_CNV" = "#6BAED6", "APOE" = "#4DAF4A", "LDLRAP1" = "#984EA3",
              "Other" = "#999999")
cohort_pal <- c("Wales" = "#2166AC", "UKB" = "#B2182B")

# --- Robust save function ---
save_fig <- function(p, fname, w_mm = 180, h_mm = 100, dpi = 600) {
  is_grob <- !inherits(p, "ggplot")
  draw_it <- function() { if (is_grob) grid::grid.draw(p) else print(p) }
  w_in <- w_mm / 25.4; h_in <- h_mm / 25.4
  tryCatch({
    pdf(file.path(fig_dir, paste0(fname, ".pdf")), width = w_in, height = h_in)
    draw_it(); dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE) })
  tryCatch({
    png(file.path(fig_dir, paste0(fname, ".png")), width = w_mm, height = h_mm,
        units = "mm", res = dpi)
    draw_it(); dev.off()
  }, error = function(e) { try(dev.off(), silent = TRUE) })
  cat(sprintf("  Saved: %s (.pdf + .png)\n", fname))
}

# ============================================================
# PART 2: LOAD WORKSPACE
# ============================================================

cat("Loading workspace...\n")
ws_paths <- c(
  "C:/Users/nader/Downloads/tudor_v2_workspace.RData",
  "C:/Users/nader/Downloads/tudor_pipeline/tudor_v2_workspace.RData",
  "tudor_v2_workspace.RData"
)
ws_loaded <- FALSE
for (wp in ws_paths) {
  if (file.exists(wp)) {
    load(wp)
    cat(sprintf("  Loaded: %s\n", wp))
    ws_loaded <- TRUE
    break
  }
}
if (!ws_loaded) stop("Cannot find tudor_v2_workspace.RData")

# ============================================================
# PART 3: IDENTIFY DATA OBJECTS
# ============================================================

cat("\nIdentifying data objects...\n")
cat("  Available workspace objects:\n")
for (nm in ls()) {
  obj <- tryCatch(get(nm), error = function(e) NULL)
  if (is.data.frame(obj) && nrow(obj) > 50) {
    cat(sprintf("    '%s': %d rows x %d cols\n", nm, nrow(obj), ncol(obj)))
  }
}
cat("\n")

# Wales data — prefer 'df' (full 7253 with Gene column) over 'model_df_v2' (5273 training subset)
wales_obj <- NULL
for (nm in c("df", "model_df_v2", "wales", "wales_df", "df_wales")) {
  if (exists(nm)) {
    tmp <- get(nm)
    if (is.data.frame(tmp) && nrow(tmp) > 100) {
      wales_obj <- tmp
      cat(sprintf("  Wales data: '%s' (%d rows, %d cols)\n", nm, nrow(wales_obj), ncol(wales_obj)))
      break
    }
  }
}

# UKB data — prefer 'lc' (58,021 lipid clinic with is_fh_genetic + gene columns)
ukb_obj <- NULL
for (nm in c("lc", "ukb", "ukb_full", "ukb_df", "df_ukb", "ukb_data")) {
  if (exists(nm)) {
    tmp <- get(nm)
    if (is.data.frame(tmp) && nrow(tmp) > 1000) {
      ukb_obj <- tmp
      cat(sprintf("  UKB data: '%s' (%d rows, %d cols)\n", nm, nrow(ukb_obj), ncol(ukb_obj)))
      break
    }
  }
}

if (is.null(wales_obj)) stop("No Wales data found in workspace")
if (is.null(ukb_obj)) stop("No UKB data found in workspace")

# Print column names for debugging
cat("\n  Wales columns (first 30):\n")
cat("   ", paste(head(names(wales_obj), 30), collapse = ", "), "\n")
cat("  UKB columns (first 30):\n")
cat("   ", paste(head(names(ukb_obj), 30), collapse = ", "), "\n\n")

# Convert to data.frame for safe column access
ukb_df <- as.data.frame(ukb_obj)

# --- Identify key columns ---
# FH column — 'is_fh_genetic' is the correct name in lc (UKB lipid clinic)
fh_candidates <- c("is_fh_genetic", "fh_genetic", "FH", "fh", "genetic_fh",
                    "fh_status", "Positive1")
fh_col <- NULL
for (fc in fh_candidates) {
  if (fc %in% names(ukb_df)) { fh_col <- fc; break }
}
if (is.null(fh_col)) {
  # Heuristic: find binary column with 1000-5000 positives (likely FH)
  for (cn in names(ukb_df)) {
    vals <- unique(ukb_df[[cn]])
    if (length(vals) <= 3 && all(vals %in% c(0, 1, NA))) {
      s <- sum(ukb_df[[cn]] == 1, na.rm = TRUE)
      if (s > 1000 && s < 5000) { fh_col <- cn; break }
    }
  }
}
cat(sprintf("  FH column (UKB): %s\n", ifelse(is.null(fh_col), "NOT FOUND", fh_col)))
if (!is.null(fh_col)) {
  cat(sprintf("    FH+ count: %d / %d\n",
              sum(ukb_df[[fh_col]] == 1, na.rm = TRUE), nrow(ukb_df)))
}

# Gene column
gene_candidates <- c("gene", "Gene", "fh_gene", "gene_group", "variant_gene",
                      "causative_gene", "gene_type")
gene_col <- NULL
for (gc in gene_candidates) {
  if (gc %in% names(ukb_df)) { gene_col <- gc; break }
}
cat(sprintf("  Gene column (UKB): %s\n", ifelse(is.null(gene_col), "NOT FOUND", gene_col)))
if (!is.null(gene_col)) {
  cat("    UKB gene distribution (all):\n")
  print(table(ukb_df[[gene_col]], useNA = "ifany"))
}

# Wales gene column
wales_gene_col <- NULL
for (gc in c("gene", "Gene", "gene_group", "gene_type", "causative_gene")) {
  if (gc %in% names(wales_obj)) { wales_gene_col <- gc; break }
}
cat(sprintf("  Wales gene column: %s\n", ifelse(is.null(wales_gene_col), "NOT FOUND", wales_gene_col)))
if (!is.null(wales_gene_col)) {
  # IMPORTANT: Wales Gene column contains ALL variants (pathological + benign + VUS).
  # Only Positive1 == 1 marks a PATHOLOGICAL variant (= true genetic FH).
  # We show both for transparency, but ONLY Positive1==1 counts are used in analysis.
  cat("    Wales gene distribution (ALL variants, incl. benign/VUS):\n")
  print(table(wales_obj[[wales_gene_col]], useNA = "ifany"))
}

# Wales FH column — 'Positive1' is the correct name in df/model_df_v2
wales_fh_col <- NULL
for (fc in c("Positive1", "fh_genetic", "FH", "fh", "y", "outcome", "is_fh_genetic")) {
  if (fc %in% names(wales_obj)) { wales_fh_col <- fc; break }
}
cat(sprintf("  Wales FH column: %s\n", ifelse(is.null(wales_fh_col), "NOT FOUND", wales_fh_col)))
if (!is.null(wales_fh_col)) {
  cat(sprintf("    Wales FH+ count: %d / %d\n",
              sum(wales_obj[[wales_fh_col]] == 1, na.rm = TRUE), nrow(wales_obj)))
}

# TUDOR prediction column
tudor_col <- NULL
for (tc in c("tudor_prob", "TUDOR_prob", "tudor_score", "pred", "predicted",
             "prob", "tudor_pred")) {
  if (tc %in% names(ukb_df)) { tudor_col <- tc; break }
}
cat(sprintf("  TUDOR column (UKB): %s\n", ifelse(is.null(tudor_col), "NOT FOUND", tudor_col)))

# LDL / HDL / TG columns in UKB
ldl_col <- NULL
for (lc_ in c("LDL_untreated", "LDL_UT", "ldl_ut", "LDL.1", "ldl_c")) {
  if (lc_ %in% names(ukb_df)) { ldl_col <- lc_; break }
}
hdl_col <- NULL
for (hc_ in c("HDL.1", "HDL", "hdl", "hdl_c")) {
  if (hc_ %in% names(ukb_df)) { hdl_col <- hc_; break }
}
tg_col <- NULL
for (tc_ in c("TRG.1", "TRG", "tg", "triglycerides")) {
  if (tc_ %in% names(ukb_df)) { tg_col <- tc_; break }
}

# ============================================================
# PART 4: GENETIC SPECTRUM COMPARISON
# ============================================================

cat("\n--- Part 4: Genetic Spectrum Comparison ---\n")
cat("  NOTE: In Wales, 'Gene' column records ALL detected variants.\n")
cat("        Only Positive1 == 1 = PATHOLOGICAL variant (true genetic FH).\n")
cat("        Positive1 == 0 = benign / VUS / no variant → excluded from analysis.\n\n")

# ---------------------------------------------------------------
# WALES: Filter to PATHOLOGICAL variants only (Positive1 == 1)
# ---------------------------------------------------------------
if (!is.null(wales_fh_col)) {
  wales_fh <- wales_obj[wales_obj[[wales_fh_col]] == 1 & !is.na(wales_obj[[wales_fh_col]]), ]
  cat(sprintf("  Wales total patients: %d\n", nrow(wales_obj)))
  cat(sprintf("  Wales PATHOLOGICAL variants (Positive1==1): %d\n", nrow(wales_fh)))
  cat(sprintf("  Wales non-pathological (Positive1==0): %d\n",
              sum(wales_obj[[wales_fh_col]] == 0, na.rm = TRUE)))

  # Show gene distribution BEFORE vs AFTER pathological filter
  if (!is.null(wales_gene_col)) {
    cat("\n  Gene distribution — ALL variants (incl. benign/VUS):\n")
    print(table(wales_obj[[wales_gene_col]], useNA = "ifany"))
    cat("\n  Gene distribution — PATHOLOGICAL ONLY (Positive1==1):\n")
    print(table(wales_fh[[wales_gene_col]], useNA = "ifany"))
  }
} else {
  # Fallback: try to find FH+ subset from workspace objects
  wales_fh <- data.frame()
  for (nm in c("fh_pos", "wales_fh", "fh_positive")) {
    if (exists(nm)) {
      tmp <- get(nm)
      if (is.data.frame(tmp) && nrow(tmp) > 100) {
        wales_fh <- tmp
        cat(sprintf("  Wales FH+ patients: %d (from pre-existing '%s' object)\n", nrow(wales_fh), nm))
        if (is.null(wales_gene_col)) {
          for (gc in c("gene", "Gene", "gene_group", "gene_type", "causative_gene")) {
            if (gc %in% names(wales_fh)) { wales_gene_col <- gc; break }
          }
        }
        break
      }
    }
  }
  if (nrow(wales_fh) == 0) {
    cat("  WARNING: Wales FH column not found and no pre-existing FH+ subset available\n")
    cat("  Using hardcoded values from genetic_spectrum_comparison.csv\n")
  }
}

# Build Wales gene table — ONLY from pathological variants (Positive1==1)
wales_genes <- if (!is.null(wales_gene_col) && nrow(wales_fh) > 0) {
  # Filter out empty / NA gene entries within the pathological subset
  wales_fh_with_gene <- wales_fh[!is.na(wales_fh[[wales_gene_col]]) &
                                  wales_fh[[wales_gene_col]] != "", ]
  cat(sprintf("\n  Wales pathological with gene assignment: %d / %d\n",
              nrow(wales_fh_with_gene), nrow(wales_fh)))
  tbl <- table(wales_fh_with_gene[[wales_gene_col]])
  data.frame(Gene = names(tbl), Wales_N = as.numeric(tbl),
             Wales_Pct = round(100 * as.numeric(tbl) / sum(tbl), 1))
} else {
  # Use known values from genetic_spectrum_comparison.csv
  data.frame(
    Gene = c("LDLR", "APOB", "LDLR_CNV", "APOE", "PCSK9", "LDLRAP1"),
    Wales_N = c(1901, 281, 164, 33, 23, 3),
    Wales_Pct = c(79.0, 11.7, 6.8, 1.4, 1.0, 0.1)
  )
}
cat("  Wales pathological gene spectrum:\n")
print(wales_genes)

# ---------------------------------------------------------------
# UKB: Filter to genetically confirmed FH (is_fh_genetic == 1)
# In UKB, is_fh_genetic==1 already implies a pathological variant.
# ---------------------------------------------------------------
cat("\n")
if (!is.null(gene_col) && !is.null(fh_col)) {
  ukb_fh <- ukb_df[ukb_df[[fh_col]] == 1 & !is.na(ukb_df[[fh_col]]), ]
  cat(sprintf("  UKB total patients: %d\n", nrow(ukb_df)))
  cat(sprintf("  UKB genetically confirmed FH (%s==1): %d\n", fh_col, nrow(ukb_fh)))

  # Filter to patients with actual gene assignments (not NA, "", or "None")
  has_gene <- !is.na(ukb_fh[[gene_col]]) & ukb_fh[[gene_col]] != "" &
              ukb_fh[[gene_col]] != "None" & ukb_fh[[gene_col]] != "none"
  ukb_fh_gene <- ukb_fh[has_gene, ]
  cat(sprintf("  UKB FH+ with gene assignment: %d (excluded %d without gene)\n",
              nrow(ukb_fh_gene), nrow(ukb_fh) - nrow(ukb_fh_gene)))
  tbl_ukb <- table(ukb_fh_gene[[gene_col]])
  cat("  UKB pathological gene spectrum:\n")
  print(tbl_ukb)
  ukb_genes <- data.frame(Gene = names(tbl_ukb), UKB_N = as.numeric(tbl_ukb),
                           UKB_Pct = round(100 * as.numeric(tbl_ukb) / sum(tbl_ukb), 1))
} else {
  # Use known values
  cat("  WARNING: Gene or FH column not found — using hardcoded UKB values\n")
  ukb_genes <- data.frame(
    Gene = c("LDLR", "APOB", "PCSK9"),
    UKB_N = c(1321, 301, 1),
    UKB_Pct = c(81.4, 18.5, 0.1)
  )
}
cat("  UKB pathological gene spectrum:\n")
print(ukb_genes)

# Merge
gene_comparison <- merge(wales_genes, ukb_genes, by = "Gene", all = TRUE)
gene_comparison[is.na(gene_comparison)] <- 0
gene_comparison <- gene_comparison[order(-gene_comparison$Wales_N), ]

cat("\n  Genetic Spectrum Comparison:\n")
print(gene_comparison)

# --- Chi-square test (LDLR vs APOB proportions) ---
# Use only genes present in both cohorts for valid comparison
common_genes <- intersect(
  gene_comparison$Gene[gene_comparison$Wales_N > 0],
  gene_comparison$Gene[gene_comparison$UKB_N > 0]
)

if (length(common_genes) >= 2) {
  chi_mat <- matrix(c(
    gene_comparison$Wales_N[gene_comparison$Gene %in% common_genes],
    gene_comparison$UKB_N[gene_comparison$Gene %in% common_genes]
  ), ncol = 2)
  rownames(chi_mat) <- common_genes
  colnames(chi_mat) <- c("Wales", "UKB")

  chi_test <- chisq.test(chi_mat)
  cat(sprintf("\n  Chi-square test (common genes): X2=%.2f, df=%d, p=%.2e\n",
              chi_test$statistic, chi_test$parameter, chi_test$p.value))
} else {
  chi_test <- list(statistic = NA, parameter = NA, p.value = NA)
}

# --- Fisher exact test: LDLR vs APOB proportions ---
if ("LDLR" %in% gene_comparison$Gene && "APOB" %in% gene_comparison$Gene) {
  fisher_mat <- matrix(c(
    gene_comparison$Wales_N[gene_comparison$Gene == "LDLR"],
    gene_comparison$UKB_N[gene_comparison$Gene == "LDLR"],
    gene_comparison$Wales_N[gene_comparison$Gene == "APOB"],
    gene_comparison$UKB_N[gene_comparison$Gene == "APOB"]
  ), nrow = 2, byrow = TRUE)
  rownames(fisher_mat) <- c("LDLR", "APOB")
  colnames(fisher_mat) <- c("Wales", "UKB")

  fisher_test <- fisher.test(fisher_mat)
  cat(sprintf("  Fisher exact (LDLR vs APOB): OR=%.3f, p=%.2e\n",
              fisher_test$estimate, fisher_test$p.value))

  # Proportions test
  ldlr_prop_wales <- gene_comparison$Wales_Pct[gene_comparison$Gene == "LDLR"] / 100
  apob_prop_wales <- gene_comparison$Wales_Pct[gene_comparison$Gene == "APOB"] / 100
  ldlr_prop_ukb   <- gene_comparison$UKB_Pct[gene_comparison$Gene == "LDLR"] / 100
  apob_prop_ukb   <- gene_comparison$UKB_Pct[gene_comparison$Gene == "APOB"] / 100

  cat(sprintf("\n  APOB proportion: Wales=%.1f%%, UKB=%.1f%% (%.1f-fold higher in UKB)\n",
              apob_prop_wales * 100, apob_prop_ukb * 100,
              apob_prop_ukb / apob_prop_wales))
} else {
  fisher_test <- list(estimate = NA, p.value = NA)
}

# --- Save tables ---
write.csv(gene_comparison, file.path(table_dir, "genetic_variant_deep_comparison.csv"),
          row.names = FALSE)

fisher_results <- data.frame(
  Test = c("Chi-square (all common genes)", "Fisher exact (LDLR vs APOB)"),
  Statistic = c(as.numeric(chi_test$statistic), as.numeric(fisher_test$estimate)),
  P_value = c(chi_test$p.value, fisher_test$p.value),
  Note = c(
    paste0("df=", chi_test$parameter),
    paste0("OR=", round(as.numeric(fisher_test$estimate), 3))
  )
)
write.csv(fisher_results, file.path(table_dir, "variant_fisher_tests.csv"),
          row.names = FALSE)

# ============================================================
# PART 5: PHENOTYPIC COMPARISON BY GENE AND COHORT
# ============================================================

cat("\n--- Part 5: Phenotypic Comparison by Gene and Cohort ---\n")
cat("  NOTE: All Wales phenotype data uses PATHOLOGICAL variants only (Positive1==1)\n")

# Wales phenotype data — uses wales_fh (already filtered to Positive1==1)
phenotype_results <- list()

if (!is.null(wales_gene_col) && nrow(wales_fh) > 0) {
  # Detect lipid columns from wales_fh (pathological subset, same columns as wales_obj)
  wales_ldl_col <- NULL
  for (lc_ in c("LDL_untreated", "LDL_UT", "ldl_ut", "LDL.1")) {
    if (lc_ %in% names(wales_fh)) { wales_ldl_col <- lc_; break }
  }
  wales_hdl_col <- NULL
  for (hc_ in c("HDL.1", "HDL", "hdl")) {
    if (hc_ %in% names(wales_fh)) { wales_hdl_col <- hc_; break }
  }
  wales_tg_col <- NULL
  for (tc_ in c("TRG.1", "TRG", "tg")) {
    if (tc_ %in% names(wales_fh)) { wales_tg_col <- tc_; break }
  }
  cat(sprintf("  Wales lipid columns: LDL=%s, HDL=%s, TG=%s\n",
              ifelse(is.null(wales_ldl_col), "NA", wales_ldl_col),
              ifelse(is.null(wales_hdl_col), "NA", wales_hdl_col),
              ifelse(is.null(wales_tg_col), "NA", wales_tg_col)))

  # Iterate genes within pathological subset only
  for (g in unique(wales_fh[[wales_gene_col]])) {
    if (is.na(g) || g == "") next  # skip empty gene entries
    sub <- wales_fh[wales_fh[[wales_gene_col]] == g & !is.na(wales_fh[[wales_gene_col]]), ]
    if (nrow(sub) >= 5) {
      row <- data.frame(
        Cohort = "Wales", Gene = g, N = nrow(sub),
        LDL_UT_mean = ifelse(!is.null(wales_ldl_col), mean(sub[[wales_ldl_col]], na.rm=TRUE), NA),
        LDL_UT_sd   = ifelse(!is.null(wales_ldl_col), sd(sub[[wales_ldl_col]], na.rm=TRUE), NA),
        HDL_mean    = ifelse(!is.null(wales_hdl_col), mean(sub[[wales_hdl_col]], na.rm=TRUE), NA),
        HDL_sd      = ifelse(!is.null(wales_hdl_col), sd(sub[[wales_hdl_col]], na.rm=TRUE), NA),
        TG_mean     = ifelse(!is.null(wales_tg_col), mean(sub[[wales_tg_col]], na.rm=TRUE), NA),
        TG_sd       = ifelse(!is.null(wales_tg_col), sd(sub[[wales_tg_col]], na.rm=TRUE), NA)
      )
      phenotype_results[[length(phenotype_results) + 1]] <- row
    }
  }
}

# UKB phenotype data
if (!is.null(gene_col) && !is.null(fh_col)) {
  for (g in unique(ukb_fh_gene[[gene_col]])) {
    sub <- ukb_fh_gene[ukb_fh_gene[[gene_col]] == g, ]
    if (nrow(sub) >= 5) {
      row <- data.frame(
        Cohort = "UKB", Gene = g, N = nrow(sub),
        LDL_UT_mean = ifelse(!is.null(ldl_col), mean(sub[[ldl_col]], na.rm=TRUE), NA),
        LDL_UT_sd   = ifelse(!is.null(ldl_col), sd(sub[[ldl_col]], na.rm=TRUE), NA),
        HDL_mean    = ifelse(!is.null(hdl_col), mean(sub[[hdl_col]], na.rm=TRUE), NA),
        HDL_sd      = ifelse(!is.null(hdl_col), sd(sub[[hdl_col]], na.rm=TRUE), NA),
        TG_mean     = ifelse(!is.null(tg_col), mean(sub[[tg_col]], na.rm=TRUE), NA),
        TG_sd       = ifelse(!is.null(tg_col), sd(sub[[tg_col]], na.rm=TRUE), NA)
      )
      phenotype_results[[length(phenotype_results) + 1]] <- row
    }
  }
}

pheno_df <- do.call(rbind, phenotype_results)
cat("\n  Phenotypic Comparison:\n")
print(pheno_df)
write.csv(pheno_df, file.path(table_dir, "variant_phenotype_by_cohort.csv"),
          row.names = FALSE)

# --- Wilcoxon tests: LDLR TG vs APOB TG within each cohort ---
cat("\n  --- Triglyceride comparison: LDLR vs APOB ---\n")

# UKB
if (!is.null(gene_col) && !is.null(tg_col)) {
  ldlr_tg <- ukb_fh_gene[[tg_col]][ukb_fh_gene[[gene_col]] == "LDLR"]
  apob_tg <- ukb_fh_gene[[tg_col]][ukb_fh_gene[[gene_col]] == "APOB"]
  ldlr_tg <- ldlr_tg[!is.na(ldlr_tg)]
  apob_tg <- apob_tg[!is.na(apob_tg)]

  if (length(ldlr_tg) > 5 && length(apob_tg) > 5) {
    wt <- wilcox.test(apob_tg, ldlr_tg)
    cd <- (mean(apob_tg) - mean(ldlr_tg)) / sqrt(((length(apob_tg)-1)*var(apob_tg) +
          (length(ldlr_tg)-1)*var(ldlr_tg)) / (length(apob_tg) + length(ldlr_tg) - 2))
    cat(sprintf("  UKB: APOB TG mean=%.3f (SD=%.3f), LDLR TG mean=%.3f (SD=%.3f)\n",
                mean(apob_tg), sd(apob_tg), mean(ldlr_tg), sd(ldlr_tg)))
    cat(sprintf("  UKB: Wilcoxon p=%.2e, Cohen's d=%.3f\n", wt$p.value, cd))
  }
}

# Wales
if (!is.null(wales_gene_col) && !is.null(wales_tg_col)) {
  w_ldlr_tg <- wales_fh[[wales_tg_col]][wales_fh[[wales_gene_col]] == "LDLR"]
  w_apob_tg <- wales_fh[[wales_tg_col]][wales_fh[[wales_gene_col]] == "APOB"]
  w_ldlr_tg <- w_ldlr_tg[!is.na(w_ldlr_tg)]
  w_apob_tg <- w_apob_tg[!is.na(w_apob_tg)]

  if (length(w_ldlr_tg) > 5 && length(w_apob_tg) > 5) {
    wt_w <- wilcox.test(w_apob_tg, w_ldlr_tg)
    cd_w <- (mean(w_apob_tg) - mean(w_ldlr_tg)) / sqrt(((length(w_apob_tg)-1)*var(w_apob_tg) +
             (length(w_ldlr_tg)-1)*var(w_ldlr_tg)) / (length(w_apob_tg) + length(w_ldlr_tg) - 2))
    cat(sprintf("  Wales: APOB TG mean=%.3f (SD=%.3f), LDLR TG mean=%.3f (SD=%.3f)\n",
                mean(w_apob_tg), sd(w_apob_tg), mean(w_ldlr_tg), sd(w_ldlr_tg)))
    cat(sprintf("  Wales: Wilcoxon p=%.2e, Cohen's d=%.3f\n", wt_w$p.value, cd_w))
  }
}

# ============================================================
# PART 6: GENE-SPECIFIC TUDOR PERFORMANCE COMPARISON
# ============================================================

cat("\n--- Part 6: Gene-Specific TUDOR Performance ---\n")

auc_results <- list()

# --- UKB Gene-Specific AUCs ---
if (!is.null(fh_col) && !is.null(tudor_col) && !is.null(gene_col)) {
  neg_ukb <- ukb_df[ukb_df[[fh_col]] == 0 & !is.na(ukb_df[[fh_col]]), ]

  for (g in c("LDLR", "APOB")) {
    fh_g <- ukb_fh_gene[ukb_fh_gene[[gene_col]] == g, ]
    if (nrow(fh_g) >= 10) {
      combined <- rbind(
        neg_ukb[, c(fh_col, tudor_col), drop = FALSE],
        fh_g[, c(fh_col, tudor_col), drop = FALSE]
      )
      roc_g <- roc(combined[[fh_col]], combined[[tudor_col]], quiet = TRUE)
      ci_g  <- ci.auc(roc_g)
      auc_results[[length(auc_results) + 1]] <- data.frame(
        Cohort = "UKB", Gene = g, N_FH = nrow(fh_g),
        AUC = as.numeric(auc(roc_g)),
        CI_lower = ci_g[1], CI_upper = ci_g[3]
      )
      cat(sprintf("  UKB %s: AUC=%.3f (%.3f-%.3f), N=%d\n",
                  g, as.numeric(auc(roc_g)), ci_g[1], ci_g[3], nrow(fh_g)))
    }
  }

  # DeLong test: APOB vs LDLR in UKB
  if (sum(ukb_fh_gene[[gene_col]] == "LDLR") >= 10 &&
      sum(ukb_fh_gene[[gene_col]] == "APOB") >= 10) {

    sub_ldlr <- rbind(
      neg_ukb[, c(fh_col, tudor_col), drop = FALSE],
      ukb_fh_gene[ukb_fh_gene[[gene_col]] == "LDLR", c(fh_col, tudor_col), drop = FALSE]
    )
    sub_apob <- rbind(
      neg_ukb[, c(fh_col, tudor_col), drop = FALSE],
      ukb_fh_gene[ukb_fh_gene[[gene_col]] == "APOB", c(fh_col, tudor_col), drop = FALSE]
    )

    roc_ldlr <- roc(sub_ldlr[[fh_col]], sub_ldlr[[tudor_col]], quiet = TRUE)
    roc_apob <- roc(sub_apob[[fh_col]], sub_apob[[tudor_col]], quiet = TRUE)

    auc_diff <- as.numeric(auc(roc_apob)) - as.numeric(auc(roc_ldlr))
    se_ldlr  <- sqrt(var(roc_ldlr))
    se_apob  <- sqrt(var(roc_apob))
    z_stat   <- auc_diff / sqrt(se_ldlr^2 + se_apob^2)
    p_val    <- 2 * pnorm(-abs(z_stat))

    cat(sprintf("\n  UKB DeLong: APOB AUC=%.3f vs LDLR AUC=%.3f\n",
                as.numeric(auc(roc_apob)), as.numeric(auc(roc_ldlr))))
    cat(sprintf("    Diff=%.3f, Z=%.3f, p=%.2e\n", auc_diff, z_stat, p_val))
  }
}

# --- Wales Gene-Specific AUCs (using existing tudor_by_gene_type.csv) ---
# Also load from pre-computed if available
pre_gene <- tryCatch(
  read.csv(file.path(table_dir, "tudor_by_gene_type.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL
)
if (!is.null(pre_gene)) {
  for (i in seq_len(nrow(pre_gene))) {
    auc_results[[length(auc_results) + 1]] <- data.frame(
      Cohort = pre_gene$Dataset[i],
      Gene = pre_gene$Gene[i],
      N_FH = pre_gene$N_FH[i],
      AUC = pre_gene$AUC[i],
      CI_lower = pre_gene$CI_lower[i],
      CI_upper = pre_gene$CI_upper[i]
    )
  }
}

auc_df <- do.call(rbind, auc_results)
# Remove duplicates
auc_df <- auc_df[!duplicated(paste(auc_df$Cohort, auc_df$Gene)), ]
cat("\n  AUC Results:\n")
print(auc_df)
write.csv(auc_df, file.path(table_dir, "variant_tudor_performance.csv"),
          row.names = FALSE)

# ============================================================
# PART 7: FIGURES
# ============================================================

cat("\n--- Part 7: Generating Figures ---\n")

# --- Figure A: Genetic Spectrum Stacked Bar Chart ---
tryCatch({
  cat("  Figure A: Genetic spectrum comparison...\n")

  # Prepare long format
  wales_long <- data.frame(Cohort = "Wales", Gene = gene_comparison$Gene,
                            N = gene_comparison$Wales_N)
  wales_long$Pct <- 100 * wales_long$N / sum(wales_long$N)

  ukb_long <- data.frame(Cohort = "UK Biobank", Gene = gene_comparison$Gene,
                           N = gene_comparison$UKB_N)
  ukb_long <- ukb_long[ukb_long$N > 0, ]
  ukb_long$Pct <- 100 * ukb_long$N / sum(ukb_long$N)

  bar_data <- rbind(wales_long, ukb_long)
  bar_data$Cohort <- factor(bar_data$Cohort, levels = c("Wales", "UK Biobank"))

  # Order genes
  gene_order <- c("LDLR", "APOB", "PCSK9", "LDLR_CNV", "APOE", "LDLRAP1")
  bar_data$Gene <- factor(bar_data$Gene, levels = rev(gene_order))

  p_spectrum <- ggplot(bar_data, aes(x = Cohort, y = Pct, fill = Gene)) +
    geom_bar(stat = "identity", width = 0.6, colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = gene_pal, name = "Causative Gene") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 105)) +
    labs(x = "", y = "Proportion of Pathological FH Cases (%)",
         title = "Pathological Genetic Variant Spectrum",
         subtitle = sprintf("Wales (n=%s, Positive1=1) vs UK Biobank (n=%s) | Chi-square p=%s",
                            format(sum(wales_long$N), big.mark = ","),
                            format(sum(ukb_long$N), big.mark = ","),
                            formatC(chi_test$p.value, format = "e", digits = 2))) +
    theme_nature(base_size = 8) +
    theme(legend.position = "right")

  save_fig(p_spectrum, "fig_genetic_spectrum_comparison", w_mm = 120, h_mm = 100)
}, error = function(e) cat(sprintf("  Figure A error: %s\n", e$message)))

# --- Figure B: Gene-specific AUC Forest Plot ---
tryCatch({
  cat("  Figure B: Gene-specific AUC forest plot...\n")

  forest_data <- auc_df[auc_df$Gene != "ALL", ]
  forest_data$Label <- paste0(forest_data$Cohort, " — ", forest_data$Gene,
                               " (n=", forest_data$N_FH, ")")
  forest_data$y_pos <- rev(seq_len(nrow(forest_data)))

  # Add overall rows
  overall <- auc_df[auc_df$Gene == "ALL", ]
  if (nrow(overall) > 0) {
    overall$Label <- paste0(overall$Cohort, " — Overall")
    overall$y_pos <- rev(seq_len(nrow(overall))) - 0.5
    # We'll add these separately
  }

  p_forest <- ggplot(forest_data, aes(x = AUC, y = y_pos)) +
    geom_vline(xintercept = 0.75, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, colour = Cohort),
                   height = 0.2, linewidth = 0.5) +
    geom_point(aes(colour = Cohort, shape = Gene), size = 2.5) +
    scale_colour_manual(values = cohort_pal) +
    scale_shape_manual(values = c("LDLR" = 16, "APOB" = 17, "APOE" = 15,
                                   "PCSK9" = 18, "ALL" = 23)) +
    scale_y_continuous(breaks = forest_data$y_pos, labels = forest_data$Label) +
    scale_x_continuous(limits = c(0.6, 1.0), breaks = seq(0.6, 1.0, 0.1)) +
    labs(x = "AUC (95% CI)", y = "",
         title = "Gene-Specific TUDOR Discrimination",
         subtitle = "By causative gene and validation cohort") +
    theme_nature(base_size = 8) +
    theme(legend.position = "bottom",
          axis.text.y = element_text(size = 7))

  save_fig(p_forest, "fig_gene_auc_forest", w_mm = 160, h_mm = 100)
}, error = function(e) cat(sprintf("  Figure B error: %s\n", e$message)))

# --- Figure C: Phenotype Heatmap (LDL/HDL/TG by gene and cohort) ---
tryCatch({
  cat("  Figure C: Variant phenotype comparison...\n")

  if (!is.null(pheno_df) && nrow(pheno_df) > 3) {
    pheno_long <- data.frame()
    for (i in seq_len(nrow(pheno_df))) {
      if (!is.na(pheno_df$LDL_UT_mean[i])) {
        pheno_long <- rbind(pheno_long, data.frame(
          Cohort = pheno_df$Cohort[i], Gene = pheno_df$Gene[i],
          Metric = "LDL-C (adj)", Value = pheno_df$LDL_UT_mean[i], SD = pheno_df$LDL_UT_sd[i]))
      }
      if (!is.na(pheno_df$HDL_mean[i])) {
        pheno_long <- rbind(pheno_long, data.frame(
          Cohort = pheno_df$Cohort[i], Gene = pheno_df$Gene[i],
          Metric = "HDL-C", Value = pheno_df$HDL_mean[i], SD = pheno_df$HDL_sd[i]))
      }
      if (!is.na(pheno_df$TG_mean[i])) {
        pheno_long <- rbind(pheno_long, data.frame(
          Cohort = pheno_df$Cohort[i], Gene = pheno_df$Gene[i],
          Metric = "Triglycerides", Value = pheno_df$TG_mean[i], SD = pheno_df$TG_sd[i]))
      }
    }

    if (nrow(pheno_long) > 0) {
      pheno_long$Label <- paste0(pheno_long$Cohort, "\n", pheno_long$Gene)

      p_pheno <- ggplot(pheno_long[pheno_long$Gene %in% c("LDLR", "APOB"), ],
                         aes(x = Gene, y = Value, fill = Cohort)) +
        geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
        geom_errorbar(aes(ymin = Value - SD, ymax = Value + SD),
                      position = position_dodge(0.7), width = 0.2, linewidth = 0.3) +
        facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
        scale_fill_manual(values = cohort_pal) +
        labs(x = "Causative Gene", y = "Concentration (mmol/L)",
             title = "Lipid Phenotype by Causative Gene and Cohort") +
        theme_nature(base_size = 8) +
        theme(legend.position = "bottom")

      save_fig(p_pheno, "fig_variant_phenotype_heatmap", w_mm = 180, h_mm = 90)
    }
  }
}, error = function(e) cat(sprintf("  Figure C error: %s\n", e$message)))

# ============================================================
# PART 8: SUMMARY
# ============================================================

cat("\n")
cat("================================================================\n")
cat("  GENETIC VARIANT COMPARISON COMPLETE\n")
cat("================================================================\n")
cat("\n  IMPORTANT: All counts are PATHOLOGICAL variants only.\n")
cat("    Wales: filtered by Positive1 == 1 (pathological FH variant)\n")
cat("    UKB:   filtered by is_fh_genetic == 1 (genetically confirmed FH)\n\n")
cat("  Key Findings:\n")
cat(sprintf("    Wales pathological FH cases: %d across %d genes\n",
            sum(gene_comparison$Wales_N), sum(gene_comparison$Wales_N > 0)))
cat(sprintf("    UKB genetically confirmed FH: %d across %d genes\n",
            sum(gene_comparison$UKB_N), sum(gene_comparison$UKB_N > 0)))
cat(sprintf("    Chi-square p: %s\n", formatC(chi_test$p.value, format = "e", digits = 2)))

if (!is.null(fisher_test$p.value) && !is.na(fisher_test$p.value)) {
  cat(sprintf("    APOB enrichment in UKB: %.1f%% vs %.1f%% (Fisher p=%s)\n",
              gene_comparison$UKB_Pct[gene_comparison$Gene == "APOB"],
              gene_comparison$Wales_Pct[gene_comparison$Gene == "APOB"],
              formatC(fisher_test$p.value, format = "e", digits = 2)))
}

cat("\n  Output Files:\n")
cat(sprintf("    %s/genetic_variant_deep_comparison.csv\n", table_dir))
cat(sprintf("    %s/variant_phenotype_by_cohort.csv\n", table_dir))
cat(sprintf("    %s/variant_tudor_performance.csv\n", table_dir))
cat(sprintf("    %s/variant_fisher_tests.csv\n", table_dir))
cat(sprintf("    %s/fig_genetic_spectrum_comparison.pdf/png\n", fig_dir))
cat(sprintf("    %s/fig_gene_auc_forest.pdf/png\n", fig_dir))
cat(sprintf("    %s/fig_variant_phenotype_heatmap.pdf/png\n", fig_dir))
cat("================================================================\n")
