# ==============================================================================
# 15_tables_to_word.R
# TUDOR FH Diagnostic Algorithm: Export All Tables to Word Documents
# ==============================================================================
#
# Authors: Nader Genedy, Soha Zouwail
# Institution: Cardiff and Vale University Health Board
#
# Creates publication-quality Word documents for all manuscript tables
# using the officer + flextable packages.
#
# INSTALL IF NEEDED:
#   install.packages(c("officer", "flextable", "ftExtra"))
# ==============================================================================

cat("\n================================================================\n")
cat("  TUDOR: Export Tables to Word Documents\n")
cat("  Script 15 | ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
cat("================================================================\n\n")

# --- Check / install packages ---
required_pkgs <- c("officer", "flextable")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(officer)
  library(flextable)
})

cat(sprintf("  officer version: %s\n", packageVersion("officer")))
cat(sprintf("  flextable version: %s\n", packageVersion("flextable")))

# --- Paths ---
table_dir <- "C:/Users/nader/Downloads/tudor_pipeline_output/tables"
out_dir   <- "C:/Users/nader/Downloads/tudor_pipeline_output/word_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Formatting helpers ---
set_flextable_defaults(
  font.family = "Arial",
  font.size   = 9,
  padding      = 3
)

style_table <- function(ft, title = NULL) {
  ft <- ft %>%
    theme_booktabs() %>%
    fontsize(size = 9, part = "all") %>%
    fontsize(size = 9.5, part = "header") %>%
    font(fontname = "Arial", part = "all") %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    align(j = 1, align = "left", part = "body") %>%
    border_inner_h(border = fp_border(color = "#D0D0D0", width = 0.5)) %>%
    border_outer(border = fp_border(color = "black", width = 1)) %>%
    autofit(add_w = 0.15, add_h = 0.05)

  if (!is.null(title)) {
    ft <- set_caption(ft, caption = title,
                      style = "Table Caption",
                      autonum = NULL)
  }
  ft
}

save_word <- function(ft, filename, title_text = NULL, footnote_text = NULL,
                      landscape = FALSE, w = 7.5) {
  doc <- read_docx()

  # Set landscape if needed
  if (landscape) {
    doc <- body_end_section_landscape(doc)
  }

  # Title
  if (!is.null(title_text)) {
    doc <- body_add_par(doc, title_text, style = "heading 2")
    doc <- body_add_par(doc, "")
  }

  # Table
  ft <- width(ft, width = w / ncol_keys(ft))
  doc <- body_add_flextable(doc, ft, align = "center")

  # Footnote
  if (!is.null(footnote_text)) {
    doc <- body_add_par(doc, "")
    doc <- body_add_par(doc, footnote_text, style = "Normal")
  }

  path <- file.path(out_dir, paste0(filename, ".docx"))
  print(doc, target = path)
  cat(sprintf("  Saved: %s\n", basename(path)))
  invisible(path)
}

cat("  Setup complete.\n\n")


# ============================================================
# TABLE 1a: Wales Patient Demographics
# ============================================================
cat("  Table 1a: Wales Patient Demographics\n")
tryCatch({
  wales_t1 <- read.csv(file.path(table_dir, "wales_deep_table1.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)

  # Clean column names for display
  names(wales_t1) <- c("Variable", "N", "Overall", "FH+", "FH-", "P-value", "Cohen's d")

  # Format p-values
  wales_t1$`P-value` <- ifelse(wales_t1$`P-value` == "", "",
                                ifelse(as.numeric(wales_t1$`P-value`) < 0.001,
                                       sprintf("%.2e", as.numeric(wales_t1$`P-value`)),
                                       sprintf("%.3f", as.numeric(wales_t1$`P-value`))))
  wales_t1$`P-value`[is.na(wales_t1$`P-value`)] <- ""

  ft1a <- flextable(wales_t1) %>%
    style_table(title = "Table 1a. Baseline Characteristics of the Wales Development Cohort") %>%
    italic(j = "P-value", part = "header") %>%
    italic(j = "Cohen's d", part = "header") %>%
    color(j = "P-value", color = "#333333") %>%
    width(j = 1, width = 2.5) %>%
    width(j = 2, width = 0.5) %>%
    width(j = 3:5, width = 1.2) %>%
    width(j = 6, width = 0.9) %>%
    width(j = 7, width = 0.7)

  save_word(ft1a, "Table_1a_Wales_Demographics",
            title_text = "Table 1a. Baseline Characteristics of the Wales Development Cohort (N = 7,253)",
            footnote_text = "Data presented as mean (SD) or n (%). P-values from Welch's t-test (continuous) or chi-squared test (categorical). Cohen's d: effect size for continuous variables. FH status determined by genetic sequencing.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 1b: UKB Patient Demographics
# ============================================================
cat("  Table 1b: UKB Patient Demographics\n")
tryCatch({
  ukb_t1 <- read.csv(file.path(table_dir, "ukb_deep_table1.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)

  names(ukb_t1) <- c("Variable", "N", "Overall", "FH-", "FH+", "P-value", "Cohen's d")
  # Reorder to FH+ before FH-
  ukb_t1 <- ukb_t1[, c("Variable", "N", "Overall", "FH+", "FH-", "P-value", "Cohen's d")]

  ukb_t1$`P-value` <- ifelse(ukb_t1$`P-value` == "", "",
                               ifelse(as.numeric(ukb_t1$`P-value`) < 0.001,
                                      sprintf("%.2e", as.numeric(ukb_t1$`P-value`)),
                                      sprintf("%.3f", as.numeric(ukb_t1$`P-value`))))
  ukb_t1$`P-value`[is.na(ukb_t1$`P-value`)] <- ""

  ft1b <- flextable(ukb_t1) %>%
    style_table(title = "Table 1b. Baseline Characteristics of the UK Biobank Lipid Clinic Cohort") %>%
    italic(j = "P-value", part = "header") %>%
    italic(j = "Cohen's d", part = "header") %>%
    color(j = "P-value", color = "#333333") %>%
    width(j = 1, width = 2.5) %>%
    width(j = 2, width = 0.5) %>%
    width(j = 3:5, width = 1.2) %>%
    width(j = 6, width = 0.9) %>%
    width(j = 7, width = 0.7)

  save_word(ft1b, "Table_1b_UKB_Demographics",
            title_text = "Table 1b. Baseline Characteristics of the UK Biobank Lipid Clinic Cohort (N = 58,021)",
            footnote_text = "Data presented as mean (SD) or n (%). P-values from Welch's t-test (continuous) or chi-squared test (categorical). Cohen's d: effect size for continuous variables. FH status from exome sequencing. LDL-C corrected: treatment-adjusted using statin-specific correction factors.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 2: TUDOR v2 Model Coefficients
# ============================================================
cat("  Table 2: Model Coefficients\n")
tryCatch({
  params <- read.csv(file.path(table_dir, "model_parameters.csv"),
                      stringsAsFactors = FALSE)

  # Clean feature names for display
  display_names <- c(
    "(Intercept)"    = "Intercept",
    "LDL_untreated"  = "Treatment-adjusted LDL-C",
    "Trig_Filter"    = "Triglyceride Filter (LDL_UT / [TG + 0.1])",
    "HDL.1"          = "HDL-C",
    "TRG.1"          = "Triglycerides",
    "Age_at_LDL1"    = "Age at first LDL measurement",
    "Gender_num"     = "Sex (Male = 1)",
    "I_Vs_R"         = "Index case x LDL_untreated",
    "corneal_less_40"= "Corneal arcus (< 40 years)",
    "ascvd_combine"  = "Premature ASCVD",
    "pers_residual"  = "Personal statin residual (LOOCV)"
  )

  params$Feature_Display <- ifelse(params$Feature %in% names(display_names),
                                    display_names[params$Feature], params$Feature)
  params$Coefficient <- sprintf("%.4f", params$Coefficient)
  params$Scale_Center <- ifelse(is.na(params$Scale_Center), "—",
                                 sprintf("%.4f", params$Scale_Center))
  params$Scale_SD <- ifelse(is.na(params$Scale_SD), "—",
                             sprintf("%.4f", params$Scale_SD))

  t2 <- params[, c("Feature_Display", "Coefficient", "Scale_Center", "Scale_SD")]
  names(t2) <- c("Feature", "Coefficient (beta)", "Training Mean", "Training SD")

  ft2 <- flextable(t2) %>%
    style_table(title = "Table 2. TUDOR v2 Elastic Net Model Coefficients and Standardisation Parameters") %>%
    width(j = 1, width = 3.0) %>%
    width(j = 2:4, width = 1.3) %>%
    align(j = 2:4, align = "center", part = "body") %>%
    bold(i = 1, part = "body")  # Bold the intercept row

  save_word(ft2, "Table_2_Model_Coefficients",
            title_text = "Table 2. TUDOR v2 Elastic Net Model Coefficients and Standardisation Parameters",
            footnote_text = "Elastic net (alpha = 0.5) trained on Wales development cohort (N = 7,253). Features standardised to zero mean, unit variance before coefficient estimation. The personal statin residual is set to 0 in the calculator (average statin responder assumption). P(FH) = 1 / (1 + exp(-(intercept + sum(beta_i * z_i)))).")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 3: External Validation Performance
# ============================================================
cat("  Table 3: External Validation Performance\n")
tryCatch({
  val <- read.csv(file.path(table_dir, "new_cohort_validation.csv"),
                   stringsAsFactors = FALSE)

  val$AUC_CI <- sprintf("%.3f [%.3f-%.3f]",
                          val$Value, val$CI_lower, val$CI_upper)
  val$AUC_CI <- ifelse(is.na(val$CI_lower),
                         sprintf("%.3f", val$Value), val$AUC_CI)

  t3 <- val[, c("Metric", "AUC_CI", "N", "FH_cases")]
  names(t3) <- c("Metric", "Value [95% CI]", "N", "FH+ Cases")
  t3$N <- formatC(t3$N, format = "d", big.mark = ",")

  ft3 <- flextable(t3) %>%
    style_table(title = "Table 3. TUDOR External Validation in UK Biobank Lipid Clinic Cohort") %>%
    width(j = 1, width = 2.0) %>%
    width(j = 2, width = 2.5) %>%
    width(j = 3:4, width = 1.0) %>%
    align(j = 2:4, align = "center", part = "body")

  save_word(ft3, "Table_3_External_Validation",
            title_text = "Table 3. TUDOR External Validation Performance in UK Biobank Lipid Clinic Cohort",
            footnote_text = "AUC: area under the receiver operating characteristic curve. eDLCN: electronic Dutch Lipid Clinic Network score. Calibration slope reflects prediction accuracy at 1.26% FH prevalence (vs 33% training prevalence). Youden's index used for optimal threshold selection.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 4: Gene-Specific Discrimination
# ============================================================
cat("  Table 4: Gene-Specific Discrimination\n")
tryCatch({
  gene <- read.csv(file.path(table_dir, "tudor_by_gene_type.csv"),
                    stringsAsFactors = FALSE)

  gene$AUC_CI <- sprintf("%.3f [%.3f-%.3f]", gene$AUC, gene$CI_lower, gene$CI_upper)

  t4 <- gene[, c("Dataset", "Gene", "N_FH", "AUC_CI")]
  names(t4) <- c("Dataset", "Gene", "N (FH+)", "AUROC [95% CI]")

  ft4 <- flextable(t4) %>%
    style_table(title = "Table 4. Gene-Specific TUDOR Discrimination") %>%
    merge_v(j = "Dataset") %>%
    width(j = 1, width = 1.0) %>%
    width(j = 2, width = 1.0) %>%
    width(j = 3, width = 0.8) %>%
    width(j = 4, width = 2.5) %>%
    align(j = 2:4, align = "center", part = "body") %>%
    bold(i = ~ Gene == "ALL", part = "body") %>%
    hline(i = 4, border = fp_border(color = "black", width = 0.5))

  save_word(ft4, "Table_4_Gene_Specific_Discrimination",
            title_text = "Table 4. Gene-Specific TUDOR Discrimination by Causative Gene",
            footnote_text = "AUROC: area under the receiver operating characteristic curve. Each gene-specific AUC computed by combining that gene's FH+ cases with all FH-negative controls. ALL: all genetically confirmed FH cases regardless of gene. APOB shows significantly higher discrimination than LDLR in UKB (DeLong test).")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 5: ApoB Augmentation Results
# ============================================================
cat("  Table 5: ApoB Augmentation Results\n")
tryCatch({
  apob <- read.csv(file.path(table_dir, "apob_augmented_results.csv"),
                    stringsAsFactors = FALSE)

  apob$AUC_CI <- sprintf("%.3f [%.3f-%.3f]", apob$AUC, apob$CI_lower, apob$CI_upper)
  apob$Delta  <- ifelse(apob$Delta_AUC == 0, "Ref",
                          sprintf("%+.3f", apob$Delta_AUC))
  apob$P_val  <- ifelse(is.na(apob$DeLong_p), "Ref",
                          ifelse(apob$DeLong_p < 0.001,
                                 sprintf("%.2e", apob$DeLong_p),
                                 sprintf("%.3f", apob$DeLong_p)))
  apob$N_fmt  <- formatC(apob$N, format = "d", big.mark = ",")

  t5 <- apob[, c("Model", "AUC_CI", "Delta", "P_val", "N_fmt")]
  names(t5) <- c("Model", "AUROC [95% CI]", "Delta AUC", "DeLong P", "N")

  ft5 <- flextable(t5) %>%
    style_table(title = "Table 5. ApoB Augmentation of TUDOR in UK Biobank Lipid Clinic Cohort") %>%
    width(j = 1, width = 2.8) %>%
    width(j = 2, width = 2.0) %>%
    width(j = 3, width = 0.8) %>%
    width(j = 4, width = 0.9) %>%
    width(j = 5, width = 0.8) %>%
    align(j = 2:5, align = "center", part = "body") %>%
    bold(i = 1, part = "body") %>%
    color(i = ~ `DeLong P` != "Ref" & as.numeric(gsub(".*e.*", "0", `DeLong P`)) < 0.05,
          j = "DeLong P", color = "#B2182B")

  save_word(ft5, "Table_5_ApoB_Augmentation",
            title_text = "Table 5. ApoB Augmentation of TUDOR in UK Biobank Lipid Clinic Cohort",
            footnote_text = "DeLong test comparing each augmented model to the base TUDOR model. Model A: adds ApoB alone. Model B: adds ApoB/LDL-C ratio. Model C: adds both ApoB and ApoB/LDL-C ratio. Model D: adds ApoB, ApoB/LDL-C ratio, and Lp(a). Ref: reference model. N differs for Model D due to Lp(a) missingness.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 6: Cross-Population Comparison
# ============================================================
cat("  Table 6: Cross-Population Comparison\n")
tryCatch({
  cross <- read.csv(file.path(table_dir, "cross_population_comparison.csv"),
                     stringsAsFactors = FALSE)

  names(cross) <- c("Variable", "Wales N", "Wales", "UKB N", "UKB", "SMD", "P-value")
  cross$`Wales N` <- formatC(cross$`Wales N`, format = "d", big.mark = ",")
  cross$`UKB N`   <- formatC(as.numeric(cross$`UKB N`), format = "d", big.mark = ",")

  ft6 <- flextable(cross) %>%
    style_table(title = "Table 6. Cross-Population Comparison: Wales Index Cases vs UK Biobank Lipid Clinic") %>%
    width(j = 1, width = 2.2) %>%
    width(j = 2:3, width = 1.0) %>%
    width(j = 4:5, width = 1.0) %>%
    width(j = 6:7, width = 0.7) %>%
    align(j = 2:7, align = "center", part = "body") %>%
    italic(j = "P-value", part = "header") %>%
    italic(j = "SMD", part = "header")

  save_word(ft6, "Table_6_Cross_Population",
            title_text = "Table 6. Cross-Population Comparison: Wales Index Cases vs UK Biobank Lipid Clinic",
            footnote_text = "Wales data limited to index cases (first referred family member). SMD: standardised mean difference. Comparisons highlight the distributional shift between specialist lipid clinic (Wales) and population-based lipid clinic (UKB) settings.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 7: Genotype-Phenotype by Gene
# ============================================================
cat("  Table 7: Genotype-Phenotype by Gene\n")
tryCatch({
  gp <- read.csv(file.path(table_dir, "genotype_phenotype.csv"),
                  stringsAsFactors = FALSE)

  gp$LDL_UT  <- sprintf("%.1f (%.1f)", gp$LDL_UT_mean, gp$LDL_UT_sd)
  gp$HDL     <- sprintf("%.2f (%.2f)", gp$HDL_mean, gp$HDL_sd)
  gp$TRG     <- sprintf("%.2f (%.2f)", gp$TRG_mean, gp$TRG_sd)

  t7 <- gp[, c("Dataset", "Gene", "N", "LDL_UT", "HDL", "TRG")]
  names(t7) <- c("Dataset", "Gene", "N", "LDL-C Adjusted, mean (SD)", "HDL-C, mean (SD)", "Triglycerides, mean (SD)")

  ft7 <- flextable(t7) %>%
    style_table(title = "Table 7. Genotype-Phenotype Associations by Causative Gene") %>%
    merge_v(j = "Dataset") %>%
    width(j = 1, width = 0.8) %>%
    width(j = 2, width = 0.8) %>%
    width(j = 3, width = 0.5) %>%
    width(j = 4:6, width = 1.6) %>%
    align(j = 2:6, align = "center", part = "body") %>%
    hline(i = which(gp$Dataset == "Wales" & c(gp$Dataset[-1], "") != "Wales"),
          border = fp_border(color = "black", width = 0.5))

  save_word(ft7, "Table_7_Genotype_Phenotype",
            title_text = "Table 7. Genotype-Phenotype Associations by Causative Gene",
            footnote_text = "LDL-C Adjusted: treatment-corrected using statin-specific factors. APOB carriers show lower triglycerides than LDLR carriers in both cohorts, reflecting preserved VLDL metabolism (the 'metabolic shield' effect). All values in mmol/L.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 8: Genetic Spectrum Comparison
# ============================================================
cat("  Table 8: Genetic Spectrum Comparison\n")
tryCatch({
  spec <- read.csv(file.path(table_dir, "genetic_spectrum_comparison.csv"),
                    stringsAsFactors = FALSE)

  names(spec) <- c("Gene", "Wales N", "Wales %", "UKB N", "UKB %")
  spec$`Wales %` <- sprintf("%.1f%%", spec$`Wales %`)
  spec$`UKB %`   <- sprintf("%.1f%%", spec$`UKB %`)

  ft8 <- flextable(spec) %>%
    style_table(title = "Table 8. Genetic Spectrum of FH Mutations: Wales vs UK Biobank") %>%
    width(j = 1, width = 1.5) %>%
    width(j = 2:5, width = 1.0) %>%
    align(j = 2:5, align = "center", part = "body") %>%
    bold(i = 1, part = "body")  # Bold LDLR as most common

  save_word(ft8, "Table_8_Genetic_Spectrum",
            title_text = "Table 8. Genetic Spectrum of FH Mutations: Wales vs UK Biobank",
            footnote_text = "Wales: full genetic panel (sequencing + MLPA for CNVs). UKB: exome sequencing only (no CNV detection, no APOE, limited PCSK9). LDLR_CNV: LDLR copy number variants detectable by MLPA but not exome sequencing. The higher APOB proportion in UKB (18.5% vs 11.7%) may reflect referral patterns and detection methodology.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 9: Cohort Criteria Breakdown
# ============================================================
cat("  Table 9: Cohort Criteria Breakdown\n")
tryCatch({
  crit <- read.csv(file.path(table_dir, "cohort_criteria_breakdown.csv"),
                    stringsAsFactors = FALSE)

  names(crit) <- c("Criterion", "N", "FH+ Cases", "FH Prevalence")
  crit$N <- formatC(crit$N, format = "d", big.mark = ",")

  ft9 <- flextable(crit) %>%
    style_table(title = "Table 9. UK Biobank Lipid Clinic Cohort Criteria Breakdown") %>%
    width(j = 1, width = 3.0) %>%
    width(j = 2, width = 1.0) %>%
    width(j = 3, width = 0.8) %>%
    width(j = 4, width = 1.0) %>%
    align(j = 2:4, align = "center", part = "body") %>%
    bold(i = nrow(crit) - 1, part = "body") %>%  # Bold the "Any criterion" row
    hline(i = nrow(crit) - 2, border = fp_border(color = "black", width = 0.5))

  save_word(ft9, "Table_9_Cohort_Criteria",
            title_text = "Table 9. UK Biobank Lipid Clinic Cohort Selection Criteria Breakdown",
            footnote_text = "Lipid clinic cohort defined by: TC > 7.5 mmol/L OR treatment-adjusted LDL-C > 4.9 mmol/L OR premature ASCVD (MI/stroke/PVD at < 55M/60F). 'ONLY' rows indicate patients meeting that single criterion and no others. FH prevalence: proportion of genetically confirmed FH cases within each criterion subgroup.")
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# TABLE 10: Statin-Specific LDL-C Correction Factors
# ============================================================
cat("  Table 10: Statin Correction Factors\n")
tryCatch({
  statin_df <- data.frame(
    Statin = c("Atorvastatin", "", "", "",
               "Rosuvastatin", "", "", "",
               "Simvastatin", "", "", "",
               "Pravastatin", "", "",
               "Fluvastatin", "", ""),
    Dose = c("10 mg", "20 mg", "40 mg", "80 mg",
             "5 mg", "10 mg", "20 mg", "40 mg",
             "10 mg", "20 mg", "40 mg", "80 mg",
             "10 mg", "20 mg", "40 mg",
             "20 mg", "40 mg", "80 mg"),
    `Expected Reduction` = c("37%", "43%", "49%", "55%",
                              "42%", "46%", "52%", "58%",
                              "27%", "32%", "37%", "42%",
                              "20%", "24%", "29%",
                              "17%", "23%", "33%"),
    Intensity = c("Low", "Medium", "High", "High",
                  "Low", "Medium", "High", "High",
                  "Low", "Medium", "Medium", "High",
                  "Low", "Low", "Medium",
                  "Low", "Low", "Medium"),
    stringsAsFactors = FALSE, check.names = FALSE
  )

  names(statin_df) <- c("Statin", "Dose", "Expected LDL-C Reduction (%)", "Intensity")

  ft10 <- flextable(statin_df) %>%
    style_table(title = "Table 10. Statin-Specific LDL-C Correction Factors") %>%
    merge_v(j = "Statin") %>%
    width(j = 1, width = 1.5) %>%
    width(j = 2, width = 0.8) %>%
    width(j = 3, width = 2.0) %>%
    width(j = 4, width = 1.0) %>%
    align(j = 2:4, align = "center", part = "body") %>%
    hline(i = c(4, 8, 12, 15), border = fp_border(color = "#999999", width = 0.5)) %>%
    bg(i = ~ Intensity == "High", bg = "#FFF3E0")

  # Add-ons as separate mini table
  addon_df <- data.frame(
    Therapy = c("Ezetimibe", "PCSK9 inhibitor", "Bempedoic acid"),
    `Additional Reduction` = c("20%", "60%", "18%"),
    Note = c("Combined multiplicatively with statin",
             "Combined multiplicatively with statin",
             "Combined multiplicatively with statin"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(addon_df) <- c("Add-on Therapy", "Additional Reduction (%)", "Note")

  ft10b <- flextable(addon_df) %>%
    style_table() %>%
    width(j = 1, width = 1.5) %>%
    width(j = 2, width = 1.8) %>%
    width(j = 3, width = 2.5) %>%
    align(j = 2, align = "center", part = "body")

  # Write both tables to one document
  doc <- read_docx()
  doc <- body_add_par(doc, "Table 10. Statin-Specific LDL-C Correction Factors Used in TUDOR Treatment Adjustment",
                      style = "heading 2")
  doc <- body_add_par(doc, "")
  doc <- body_add_flextable(doc, ft10, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Add-on Therapies", style = "heading 3")
  doc <- body_add_par(doc, "")
  doc <- body_add_flextable(doc, ft10b, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc,
    "Correction formula: LDL_untreated = LDL_measured / (1 - reduction * compliance). Reductions from statin and add-on therapies are combined multiplicatively, capped at 85% total reduction. Compliance: Poor = 50%, Moderate = 75%, Good = 100%.",
    style = "Normal")

  path <- file.path(out_dir, "Table_10_Statin_Corrections.docx")
  print(doc, target = path)
  cat(sprintf("  Saved: %s\n", basename(path)))
}, error = function(e) cat(sprintf("  Error: %s\n", e$message)))


# ============================================================
# COMBINED: All Tables in One Word Document
# ============================================================
cat("\n  Creating combined document with all tables...\n")
tryCatch({
  doc <- read_docx()

  doc <- body_add_par(doc, "TUDOR Manuscript Tables", style = "heading 1")
  doc <- body_add_par(doc, sprintf("Generated: %s", format(Sys.time(), "%d %B %Y")))
  doc <- body_add_par(doc, "Genedy N, Zouwail S. Cardiff and Vale University Health Board")
  doc <- body_add_par(doc, "")

  # Helper to add a table section
  add_table_section <- function(doc, csv_file, title, footnote, col_names = NULL,
                                 format_func = NULL, bold_rows = NULL) {
    doc <- body_add_par(doc, "")
    doc <- body_add_par(doc, title, style = "heading 2")
    doc <- body_add_par(doc, "")

    dd <- read.csv(file.path(table_dir, csv_file), stringsAsFactors = FALSE,
                    check.names = FALSE)
    if (!is.null(col_names)) names(dd) <- col_names
    if (!is.null(format_func)) dd <- format_func(dd)

    ft <- flextable(dd) %>%
      style_table() %>%
      autofit(add_w = 0.1) %>%
      width(width = 7.0 / ncol(dd))

    if (!is.null(bold_rows)) {
      ft <- bold(ft, i = bold_rows, part = "body")
    }

    doc <- body_add_flextable(doc, ft, align = "center")
    doc <- body_add_par(doc, "")
    doc <- body_add_par(doc, footnote, style = "Normal")
    doc <- body_add_break(doc, type = "page")
    doc
  }

  # --- Table 1a: Wales ---
  wales_t1 <- read.csv(file.path(table_dir, "wales_deep_table1.csv"),
                         stringsAsFactors = FALSE, check.names = FALSE)
  names(wales_t1) <- c("Variable", "N", "Overall", "FH+", "FH-", "P-value", "Cohen's d")
  doc <- body_add_par(doc, "Table 1a. Baseline Characteristics — Wales Development Cohort (N = 7,253)", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(wales_t1) %>% style_table() %>% autofit(add_w = 0.08)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Data as mean (SD) or n (%). P-values: Welch's t-test / chi-squared. Cohen's d for continuous variables.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 1b: UKB ---
  ukb_t1 <- read.csv(file.path(table_dir, "ukb_deep_table1.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  names(ukb_t1) <- c("Variable", "N", "Overall", "FH-", "FH+", "P-value", "Cohen's d")
  ukb_t1 <- ukb_t1[, c("Variable", "N", "Overall", "FH+", "FH-", "P-value", "Cohen's d")]
  doc <- body_add_par(doc, "Table 1b. Baseline Characteristics — UK Biobank Lipid Clinic (N = 58,021)", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(ukb_t1) %>% style_table() %>% autofit(add_w = 0.08)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Data as mean (SD) or n (%). FH status from exome sequencing. LDL-C corrected using statin-specific factors.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 2: Model Coefficients ---
  params <- read.csv(file.path(table_dir, "model_parameters.csv"), stringsAsFactors = FALSE)
  display_names <- c("(Intercept)"="Intercept","LDL_untreated"="Treatment-adjusted LDL-C",
    "Trig_Filter"="Triglyceride Filter","HDL.1"="HDL-C","TRG.1"="Triglycerides",
    "Age_at_LDL1"="Age at first LDL","Gender_num"="Sex (Male=1)",
    "I_Vs_R"="Index x LDL_untreated","corneal_less_40"="Corneal arcus (<40y)",
    "ascvd_combine"="Premature ASCVD","pers_residual"="Statin residual (LOOCV)")
  params$Feature <- ifelse(params$Feature %in% names(display_names),
                            display_names[params$Feature], params$Feature)
  params$Coefficient <- round(params$Coefficient, 4)
  params$Scale_Center <- round(params$Scale_Center, 4)
  params$Scale_SD <- round(params$Scale_SD, 4)
  names(params) <- c("Feature", "Coefficient", "Training Mean", "Training SD")
  doc <- body_add_par(doc, "Table 2. TUDOR v2 Model Coefficients", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(params) %>% style_table() %>% autofit(add_w = 0.08)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Elastic net (alpha=0.5). Features standardised before fitting. P(FH) = logistic(intercept + sum(beta * z)).", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 3: Validation ---
  val <- read.csv(file.path(table_dir, "new_cohort_validation.csv"), stringsAsFactors = FALSE)
  val$Value <- round(val$Value, 3)
  val$CI_lower <- round(val$CI_lower, 3)
  val$CI_upper <- round(val$CI_upper, 3)
  val$`Value [95% CI]` <- ifelse(is.na(val$CI_lower), sprintf("%.3f", val$Value),
                                  sprintf("%.3f [%.3f-%.3f]", val$Value, val$CI_lower, val$CI_upper))
  val$N <- formatC(val$N, format="d", big.mark=",")
  t3 <- val[, c("Metric", "Value [95% CI]", "N", "FH_cases")]
  names(t3)[4] <- "FH+ Cases"
  doc <- body_add_par(doc, "Table 3. External Validation — UK Biobank", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(t3) %>% style_table() %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Calibration slope 6.33 reflects prevalence shift (33% training vs 1.26% validation). Discrimination (AUC ranking) preserved.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 4: Gene-Specific ---
  gene <- read.csv(file.path(table_dir, "tudor_by_gene_type.csv"), stringsAsFactors = FALSE)
  gene$`AUROC [95% CI]` <- sprintf("%.3f [%.3f-%.3f]", gene$AUC, gene$CI_lower, gene$CI_upper)
  t4 <- gene[, c("Dataset", "Gene", "N_FH", "AUROC [95% CI]")]
  names(t4)[3] <- "N (FH+)"
  doc <- body_add_par(doc, "Table 4. Gene-Specific TUDOR Discrimination", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(t4) %>% style_table() %>% merge_v(j="Dataset") %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "APOB AUC 0.830 vs LDLR AUC 0.717 in UKB. The metabolic shield effect (preserved TG metabolism in APOB) explains this difference.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 5: ApoB Augmentation ---
  apob <- read.csv(file.path(table_dir, "apob_augmented_results.csv"), stringsAsFactors = FALSE)
  apob$`AUROC [95% CI]` <- sprintf("%.3f [%.3f-%.3f]", apob$AUC, apob$CI_lower, apob$CI_upper)
  apob$`Delta AUC` <- ifelse(apob$Delta_AUC == 0, "Ref", sprintf("%+.3f", apob$Delta_AUC))
  apob$`DeLong P` <- ifelse(is.na(apob$DeLong_p), "Ref",
                              ifelse(apob$DeLong_p < 0.001, sprintf("%.2e", apob$DeLong_p),
                                     sprintf("%.3f", apob$DeLong_p)))
  apob$N <- formatC(apob$N, format="d", big.mark=",")
  t5 <- apob[, c("Model", "AUROC [95% CI]", "Delta AUC", "DeLong P", "N")]
  doc <- body_add_par(doc, "Table 5. ApoB Augmentation of TUDOR", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(t5) %>% style_table() %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Combined ApoB + ApoB/LDL-C ratio (Model C) yields +0.019 AUC improvement (p < 0.001). Ref: Genedy & Zouwail 2025.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 6: Cross-Population ---
  cross <- read.csv(file.path(table_dir, "cross_population_comparison.csv"), stringsAsFactors = FALSE)
  names(cross) <- c("Variable", "Wales N", "Wales", "UKB N", "UKB", "SMD", "P-value")
  doc <- body_add_par(doc, "Table 6. Cross-Population Comparison", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(cross) %>% style_table() %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Wales: specialist registry index cases. UKB: population-based lipid clinic cohort. FH prevalence 52% vs 1.26%.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 7: Genotype-Phenotype ---
  gp <- read.csv(file.path(table_dir, "genotype_phenotype.csv"), stringsAsFactors = FALSE)
  gp$`LDL-C Adj` <- sprintf("%.1f (%.1f)", gp$LDL_UT_mean, gp$LDL_UT_sd)
  gp$`HDL-C` <- sprintf("%.2f (%.2f)", gp$HDL_mean, gp$HDL_sd)
  gp$TG <- sprintf("%.2f (%.2f)", gp$TRG_mean, gp$TRG_sd)
  t7 <- gp[, c("Dataset", "Gene", "N", "LDL-C Adj", "HDL-C", "TG")]
  doc <- body_add_par(doc, "Table 7. Genotype-Phenotype by Gene", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(t7) %>% style_table() %>% merge_v(j="Dataset") %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "APOB carriers: lower TG than LDLR (metabolic shield). Values mmol/L, mean (SD).", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 8: Genetic Spectrum ---
  spec <- read.csv(file.path(table_dir, "genetic_spectrum_comparison.csv"), stringsAsFactors = FALSE)
  names(spec) <- c("Gene", "Wales N", "Wales %", "UKB N", "UKB %")
  doc <- body_add_par(doc, "Table 8. Genetic Spectrum of FH Mutations", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(spec) %>% style_table() %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Wales: sequencing + MLPA. UKB: exome only (no CNVs). Higher APOB proportion in UKB likely reflects detection bias.", style = "Normal")
  doc <- body_add_break(doc, type = "page")

  # --- Table 9: Cohort Criteria ---
  crit <- read.csv(file.path(table_dir, "cohort_criteria_breakdown.csv"), stringsAsFactors = FALSE)
  names(crit) <- c("Criterion", "N", "FH+ Cases", "FH Prevalence")
  crit$N <- formatC(crit$N, format="d", big.mark=",")
  doc <- body_add_par(doc, "Table 9. Lipid Clinic Cohort Criteria Breakdown", style = "heading 2")
  doc <- body_add_par(doc, "")
  ft <- flextable(crit) %>% style_table() %>% autofit(add_w = 0.1)
  doc <- body_add_flextable(doc, ft, align = "center")
  doc <- body_add_par(doc, "")
  doc <- body_add_par(doc, "Lipid clinic: TC>7.5 OR LDL_RW>4.9 OR premature ASCVD. Overall FH prevalence 1.26%.", style = "Normal")

  # Save combined document
  path <- file.path(out_dir, "TUDOR_All_Tables_Combined.docx")
  print(doc, target = path)
  cat(sprintf("\n  Saved combined document: %s\n", basename(path)))

}, error = function(e) cat(sprintf("  Combined document error: %s\n", e$message)))


# ============================================================
# SUMMARY
# ============================================================
cat("\n\n================================================================\n")
cat("  SUMMARY\n")
cat("================================================================\n")
cat(sprintf("  Output directory: %s\n", out_dir))
files <- list.files(out_dir, pattern = "\\.docx$")
cat(sprintf("  Word documents generated: %d\n", length(files)))
for (f in files) cat(sprintf("    - %s\n", f))
cat("\n================================================================\n")
cat("  Script 15 complete.\n")
cat("================================================================\n\n")
