# ==============================================================================
# TUDOR PIPELINE: STEP 07 — PUBLICATION TABLES & FIGURES
# ==============================================================================
# PURPOSE: Generate Lancet/BMJ-quality tables and figures.
#
# REQUIRES: All RDS outputs from steps 02-06.
#
# OUTPUTS:  - Table 1: Baseline characteristics
#           - Table 2: Subgroup AUCs
#           - Table 3: Lancet statistics (NRI, IDI, Calibration, Brier)
#           - Table 4: Sensitivity analyses summary
#           - Figure 1: ROC curves (TUDOR vs eDLCN vs LDL)
#           - Figure 2: Decision Curve Analysis
#           - Figure 3: Calibration plot
#           - Supplementary: Forest plot, biomarker comparison
#           - TRIPOD checklist
# ==============================================================================

set.seed(42)

library(data.table)
library(dplyr)
library(pROC)
library(ggplot2)

DATA_DIR <- Sys.getenv("TUDOR_DATA_DIR", unset = "")
if (DATA_DIR == "") {
  if (file.exists(file.path(getwd(), "TUDOR_UKB_Features.csv"))) {
    DATA_DIR <- getwd()
  } else {
    DATA_DIR <- "C:/Users/nader/Downloads"
  }
}
OUTPUT_DIR <- file.path(DATA_DIR, "tudor_pipeline_output")
FIG_DIR   <- file.path(OUTPUT_DIR, "figures")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== TUDOR PIPELINE: 07_tables_figures.R ===\n\n")

# ==============================================================================
# LOAD ALL RESULTS
# ==============================================================================
df <- readRDS(file.path(OUTPUT_DIR, "tudor_analysis_ready.rds"))
if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
  setnames(df, "participant.eid", "eid")
}
val <- readRDS(file.path(OUTPUT_DIR, "02_validation_results.rds"))

lancet_file <- file.path(OUTPUT_DIR, "05_lancet_stats.rds")
sens_file   <- file.path(OUTPUT_DIR, "06_sensitivity_results.rds")
bio_file    <- file.path(OUTPUT_DIR, "04_biomarker_results.rds")

has_lancet <- file.exists(lancet_file)
has_sens   <- file.exists(sens_file)
has_bio    <- file.exists(bio_file)

if (has_lancet) lancet <- readRDS(lancet_file)
if (has_sens)   sens   <- readRDS(sens_file)
if (has_bio)    bio    <- readRDS(bio_file)

hr <- df[df$cohort_high_risk == TRUE, ]

# ==============================================================================
# TABLE 1: BASELINE CHARACTERISTICS
# ==============================================================================
cat("--- Generating Table 1: Baseline Characteristics ---\n")

# Helper: median [IQR] or mean (SD)
median_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%.1f [%.1f - %.1f]", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

mean_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%.1f (%.1f)", mean(x), sd(x))
}

n_pct <- function(x, total) {
  if (length(x) == 0 || all(is.na(x))) return("0 (0.0%)")
  sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), mean(x, na.rm = TRUE) * 100)
}

# Build table for each subgroup
make_table1_row <- function(sub, label) {
  data.frame(
    Characteristic = label,
    N = nrow(sub),
    Age = mean_sd(sub$Age_at_LDL1),
    Female_pct = sprintf("%.1f%%", mean(sub$Gender_num == 0) * 100),
    LDL_mmol = median_iqr(sub$LDL_treated),
    LDL_RW_mmol = median_iqr(sub$LDL_RW),
    HDL_mmol = median_iqr(sub$HDL.1),
    Trig_mmol = median_iqr(sub$TRG.1),
    TC_mmol = median_iqr(sub$CHOL),
    BMI = mean_sd(sub$BMI_imputed),
    On_Statin_pct = sprintf("%.1f%%", mean(sub$statin_name != "None") * 100),
    FH_Cases = n_pct(sub$is_fh_genetic, nrow(sub)),
    stringsAsFactors = FALSE
  )
}

table1 <- rbind(
  make_table1_row(df, "Full Cohort"),
  make_table1_row(hr, "High Risk (LDL>4.9)"),
  make_table1_row(df[df$cohort_moderate == TRUE, ], "Moderate (2.6-4.9)"),
  make_table1_row(df[df$cohort_low_risk == TRUE, ], "Low Risk (<2.6)"),
  make_table1_row(hr[hr$is_fh_genetic == 1, ], "FH Cases (High Risk)"),
  make_table1_row(hr[hr$is_fh_genetic == 0, ], "Non-FH (High Risk)")
)

write.csv(table1, file.path(TABLE_DIR, "table1_baseline.csv"), row.names = FALSE)
cat("  Saved: table1_baseline.csv\n")

# Print summary
cat("\nTable 1 Summary:\n")
print(table1[, c("Characteristic", "N", "Age", "LDL_mmol", "FH_Cases")], row.names = FALSE)
cat("\n")

# ==============================================================================
# TABLE 2: SUBGROUP AUCs (from validation results)
# ==============================================================================
cat("--- Generating Table 2: Subgroup AUCs ---\n")
table2 <- val$subgroups
if (!is.null(table2) && nrow(table2) > 0) {
  write.csv(table2, file.path(TABLE_DIR, "table2_subgroup_aucs.csv"), row.names = FALSE)
  cat("  Saved: table2_subgroup_aucs.csv\n\n")
} else {
  cat("  SKIPPED: No subgroup results found in validation output.\n\n")
  table2 <- data.frame()
}

# ==============================================================================
# TABLE 3: LANCET STATISTICS
# ==============================================================================
if (has_lancet) {
  cat("--- Generating Table 3: Lancet Statistics ---\n")

  table3 <- data.frame(
    Metric = c(
      "TUDOR AUC", "eDLCN AUC", "DeLong p (TUDOR vs eDLCN)",
      "Categorical NRI", "Continuous NRI",
      "IDI",
      "Calibration Intercept", "Calibration Slope",
      "Hosmer-Lemeshow p",
      "Brier Score (TUDOR)", "Brier Score (eDLCN)",
      "Scaled Brier (TUDOR)", "Scaled Brier (eDLCN)"
    ),
    Value = c(
      sprintf("%.3f [%.3f - %.3f]", val$primary$tudor_auc,
              val$primary$tudor_ci[1], val$primary$tudor_ci[2]),
      sprintf("%.3f [%.3f - %.3f]", val$primary$edlcn_auc,
              val$primary$edlcn_ci[1], val$primary$edlcn_ci[2]),
      ifelse(!is.null(val$primary$delong_p_tudor_edlcn),
             sprintf("%.2e", val$primary$delong_p_tudor_edlcn), "See script 02 output"),
      sprintf("%.3f [%.3f, %.3f]", lancet$nri$categorical$nri,
              lancet$nri$categorical_ci[1], lancet$nri$categorical_ci[2]),
      sprintf("%.3f", lancet$nri$continuous$total),
      sprintf("%.4f [%.4f, %.4f]", lancet$idi$idi,
              lancet$idi$ci[1], lancet$idi$ci[2]),
      sprintf("%.3f", lancet$calibration$intercept),
      sprintf("%.3f [%.3f, %.3f]", lancet$calibration$slope,
              lancet$calibration$slope_ci[1], lancet$calibration$slope_ci[2]),
      sprintf("%.3f", lancet$calibration$hl_p),
      sprintf("%.4f", lancet$brier$tudor),
      sprintf("%.4f", lancet$brier$edlcn),
      sprintf("%.3f", lancet$brier$scaled_tudor),
      sprintf("%.3f", lancet$brier$scaled_edlcn)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(table3, file.path(TABLE_DIR, "table3_lancet_stats.csv"), row.names = FALSE)
  cat("  Saved: table3_lancet_stats.csv\n\n")
}

# ==============================================================================
# TABLE 4: SENSITIVITY ANALYSES SUMMARY
# ==============================================================================
if (has_sens) {
  cat("--- Generating Table 4: Sensitivity Summary ---\n")

  s_rows <- list()

  # S1 threshold
  if (!is.null(sens$s1_threshold)) {
    for (name in names(sens$s1_threshold)) {
      s <- sens$s1_threshold[[name]]
      if (!is.na(s$auc)) {
        s_rows <- c(s_rows, list(data.frame(
          Analysis = paste0("S1: LDL ", name),
          AUC = sprintf("%.3f", s$auc),
          CI = sprintf("[%.3f - %.3f]", s$ci[1], s$ci[2]),
          N = s$n, stringsAsFactors = FALSE
        )))
      }
    }
  }

  # S2 Friedewald
  if (!is.null(sens$s2_friedewald)) {
    s <- sens$s2_friedewald
    if (!is.na(s$friedewald_auc)) {
      s_rows <- c(s_rows, list(data.frame(
        Analysis = "S2: Friedewald LDL",
        AUC = sprintf("%.3f", s$friedewald_auc),
        CI = sprintf("[%.3f - %.3f]", s$friedewald_ci[1], s$friedewald_ci[2]),
        N = NA, stringsAsFactors = FALSE
      )))
    }
  }

  # S5 Winsorisation
  if (!is.null(sens$s5_winsorisation)) {
    s <- sens$s5_winsorisation$winsorised
    if (!is.na(s$auc)) {
      s_rows <- c(s_rows, list(data.frame(
        Analysis = "S5: Winsorised LDL",
        AUC = sprintf("%.3f", s$auc),
        CI = sprintf("[%.3f - %.3f]", s$ci[1], s$ci[2]),
        N = NA, stringsAsFactors = FALSE
      )))
    }
  }

  # S10 Statin-free
  if (!is.null(sens$s10_statin_free)) {
    s <- sens$s10_statin_free$statin_free
    if (!is.na(s$auc)) {
      s_rows <- c(s_rows, list(data.frame(
        Analysis = "S10: Statin-Free Only",
        AUC = sprintf("%.3f", s$auc),
        CI = sprintf("[%.3f - %.3f]", s$ci[1], s$ci[2]),
        N = s$n, stringsAsFactors = FALSE
      )))
    }
  }

  if (length(s_rows) > 0) {
    table4 <- do.call(rbind, s_rows)
    write.csv(table4, file.path(TABLE_DIR, "table4_sensitivity.csv"), row.names = FALSE)
    cat("  Saved: table4_sensitivity.csv\n\n")
  }
}

# ==============================================================================
# FIGURE 1: ROC CURVES
# ==============================================================================
cat("--- Generating Figure 1: ROC Curves ---\n")

roc_tudor <- val$roc_tudor
roc_edlcn <- val$roc_edlcn
roc_ldl   <- val$roc_ldl

# Build plot data
roc_to_df <- function(roc_obj, label) {
  data.frame(
    Sensitivity = roc_obj$sensitivities,
    Specificity = 1 - roc_obj$specificities,
    Model = label,
    stringsAsFactors = FALSE
  )
}

plot_data <- rbind(
  roc_to_df(roc_tudor, sprintf("TUDOR v2 (AUC: %.3f)", auc(roc_tudor))),
  roc_to_df(roc_edlcn, sprintf("eDLCN (AUC: %.3f)", auc(roc_edlcn))),
  roc_to_df(roc_ldl, sprintf("LDL-C Only (AUC: %.3f)", auc(roc_ldl)))
)

fig1 <- ggplot(plot_data, aes(x = Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(
    title = "External Validation: TUDOR v2 vs eDLCN vs LDL-C Alone",
    subtitle = sprintf("High-Risk Cohort (LDL > 4.9 mmol/L), N = %d, FH Cases = %d",
                        nrow(hr), sum(hr$is_fh_genetic)),
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = c(0.65, 0.25),
    legend.background = element_rect(fill = "white", color = "grey80"),
    plot.title = element_text(face = "bold", size = 13),
    panel.grid.minor = element_blank()
  ) +
  coord_equal()

ggsave(file.path(FIG_DIR, "figure1_roc_curves.pdf"), fig1,
       width = 7, height = 7, dpi = 300)
ggsave(file.path(FIG_DIR, "figure1_roc_curves.png"), fig1,
       width = 7, height = 7, dpi = 300)
cat("  Saved: figure1_roc_curves.pdf/.png\n")

# ==============================================================================
# FIGURE 2: DECISION CURVE ANALYSIS
# ==============================================================================
if (has_lancet) {
  cat("--- Generating Figure 2: Decision Curve ---\n")

  dca <- lancet$dca

  dca_long <- data.frame(
    Threshold = rep(dca$threshold, 4),
    Net_Benefit = c(dca$nb_tudor, dca$nb_edlcn, dca$nb_treat_all, dca$nb_treat_none),
    Strategy = rep(c("TUDOR v2", "eDLCN", "Test All", "Test None"),
                   each = nrow(dca)),
    stringsAsFactors = FALSE
  )

  fig2 <- ggplot(dca_long, aes(x = Threshold * 100, y = Net_Benefit,
                                color = Strategy, linetype = Strategy)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c("#E41A1C", "#377EB8", "grey40", "grey60")) +
    scale_linetype_manual(values = c("solid", "solid", "dashed", "dotted")) +
    labs(
      title = "Decision Curve Analysis: TUDOR v2 vs eDLCN",
      x = "Threshold Probability (%)",
      y = "Net Benefit",
      color = NULL, linetype = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = c(0.75, 0.75),
      legend.background = element_rect(fill = "white", color = "grey80"),
      plot.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    ) +
    ylim(-0.01, NA)

  ggsave(file.path(FIG_DIR, "figure2_dca.pdf"), fig2,
         width = 8, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, "figure2_dca.png"), fig2,
         width = 8, height = 6, dpi = 300)
  cat("  Saved: figure2_dca.pdf/.png\n")
}

# ==============================================================================
# FIGURE 3: CALIBRATION PLOT
# ==============================================================================
if (has_lancet) {
  cat("--- Generating Figure 3: Calibration Plot ---\n")

  calib_data <- lancet$calibration$decile_table

  fig3 <- ggplot(calib_data, aes(x = exp_rate, y = obs_rate)) +
    geom_point(aes(size = n), color = "#377EB8") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    geom_smooth(method = "loess", se = TRUE, color = "#E41A1C",
                fill = "#E41A1C", alpha = 0.2) +
    scale_size_continuous(range = c(3, 8), name = "N per decile") +
    labs(
      title = "Calibration Plot: TUDOR v2 Predicted vs Observed FH Probability",
      subtitle = sprintf("Calibration slope = %.2f, H-L p = %.3f",
                          lancet$calibration$slope, lancet$calibration$hl_p),
      x = "Mean Predicted Probability",
      y = "Observed Proportion"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank()
    ) +
    coord_equal()

  ggsave(file.path(FIG_DIR, "figure3_calibration.pdf"), fig3,
         width = 7, height = 7, dpi = 300)
  ggsave(file.path(FIG_DIR, "figure3_calibration.png"), fig3,
         width = 7, height = 7, dpi = 300)
  cat("  Saved: figure3_calibration.pdf/.png\n")
}

# ==============================================================================
# SUPPLEMENTARY: SUBGROUP FOREST PLOT
# ==============================================================================
cat("--- Generating Supplementary: Forest Plot ---\n")

if (nrow(table2) > 0 && "tudor_auc" %in% names(table2) && any(!is.na(table2$tudor_auc))) {
  forest_data <- table2[!is.na(table2$tudor_auc), ]
  forest_data$subgroup <- factor(forest_data$subgroup,
                                  levels = rev(forest_data$subgroup))

  fig_forest <- ggplot(forest_data, aes(x = tudor_auc, y = subgroup)) +
    geom_point(size = 3, color = "#E41A1C") +
    geom_errorbarh(aes(xmin = tudor_ci_lo, xmax = tudor_ci_hi), height = 0.2,
                   color = "#E41A1C") +
    geom_vline(xintercept = val$primary$tudor_auc, linetype = "dashed",
               color = "grey50") +
    labs(
      title = "TUDOR v2 AUC by Subgroup",
      subtitle = sprintf("Dashed line = overall AUC (%.3f)", val$primary$tudor_auc),
      x = "Area Under the ROC Curve (AUC)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    ) +
    xlim(0.5, 1.0)

  ggsave(file.path(FIG_DIR, "supp_forest_plot.pdf"), fig_forest,
         width = 8, height = 5, dpi = 300)
  ggsave(file.path(FIG_DIR, "supp_forest_plot.png"), fig_forest,
         width = 8, height = 5, dpi = 300)
  cat("  Saved: supp_forest_plot.pdf/.png\n")
}

# ==============================================================================
# TRIPOD CHECKLIST
# ==============================================================================
cat("\n--- Generating TRIPOD Checklist ---\n")

tripod <- data.frame(
  Item = c(
    "1. Title", "2. Abstract", "3a. Background",
    "3b. Objectives", "4a. Study design", "4b. Setting",
    "5a. Participants", "5b. Eligibility", "6a. Outcome",
    "6b. Predictors", "7a. Sample size", "7b. Missing data",
    "8. Statistical methods", "9. Risk groups", "10a. Participants flow",
    "10b. Demographics", "11. Model development", "12. Model performance",
    "13a. Results", "13b. Sensitivity", "14. Discussion",
    "15. Limitations", "16. Implications"
  ),
  Description = c(
    "External validation of TUDOR FH diagnostic model in UK Biobank",
    "AUC with 95% CI, DeLong comparison vs eDLCN",
    "FH underdiagnosis; NICE guidelines; clinical scoring limitations",
    "Validate TUDOR v2 in an independent UK population-based cohort",
    "External validation of a pre-specified logistic regression model (TRIPOD Type 4)",
    "UK Biobank, recruited 2006-2010, assessment centres across UK",
    "N in scripts; genetic FH defined by pathogenic/likely-pathogenic variants",
    "All participants with lipid panel data; cohort stratification by LDL",
    "Genetic FH status (binary); defined by ClinVar pathogenic variants",
    "LDL-C, Triglycerides, HDL-C, Age, Sex (Wales-trained weights, FIXED)",
    "See 02_external_validation.R output for Ns",
    "BMI: median imputation documented; statin: coded as None if missing",
    "Bootstrap CIs (2000 reps), DeLong test, NRI, IDI, DCA, calibration",
    "High-risk (>4.9), Moderate (2.6-4.9), Low (<2.6) LDL strata",
    "See 01_data_merge.R filtering log",
    "Table 1 from this script",
    "Weights from Wales Lipid Clinic (NOT re-estimated on UKB)",
    "AUC, calibration slope/intercept, Brier score",
    "Table 2 subgroups, Table 3 Lancet stats",
    "10 sensitivity analyses in 06_sensitivity_analyses.R",
    "Strengths: large N, genetic gold standard. Compare to DLCN/MEDPED.",
    "Statin correction assumptions, no family history data in UKB, single ethnicity bias",
    "Pathway: TUDOR as pre-genetic screening tool to prioritize genetic testing"
  ),
  Status = rep("Script-generated", 23),
  stringsAsFactors = FALSE
)

write.csv(tripod, file.path(TABLE_DIR, "tripod_checklist.csv"), row.names = FALSE)
cat("  Saved: tripod_checklist.csv\n")

# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n==========================================================\n")
cat(" OUTPUT SUMMARY\n")
cat("==========================================================\n")

tables <- list.files(TABLE_DIR, full.names = FALSE)
figures <- list.files(FIG_DIR, full.names = FALSE)

cat("Tables:\n")
for (t in tables) cat("  ", t, "\n")
cat("\nFigures:\n")
for (f in figures) cat("  ", f, "\n")

cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))
cat("\n=== 07_tables_figures.R COMPLETE ===\n")
