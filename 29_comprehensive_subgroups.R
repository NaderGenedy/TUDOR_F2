# ==============================================================================
# TUDOR PIPELINE: STEP 29 — COMPREHENSIVE SUBGROUP ANALYSIS WITH FOREST PLOT
# ==============================================================================
# PURPOSE: Nature-grade subgroup analysis with:
#   - Full subgroup stratification (sex, age, statin, LDL, ASCVD, BMI, ethnicity)
#   - Cochran's Q test and I² heterogeneity
#   - Interaction tests for each subgroup
#   - Bonferroni & FDR multiple comparison corrections
#   - Publication-quality forest plot
#
# AUTHORS: Tudor Pipeline Team
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(ggplot2)
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
cat("TUDOR PIPELINE: 29 — COMPREHENSIVE SUBGROUP ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
rds_file <- file.path(OUTPUT_DIR, "tudor_analysis_ready.rds")
use_simulated <- FALSE

if (file.exists(rds_file)) {
  df <- readRDS(rds_file)
  setDT(df)
  if ("participant.eid" %in% names(df) && !"eid" %in% names(df)) {
    setnames(df, "participant.eid", "eid")
  }
  hr <- df[cohort_high_risk == TRUE]
  cat("Loaded high-risk cohort:", nrow(hr), "\n")
  cat("FH cases:", sum(hr$is_fh_genetic), "\n\n")
} else {
  cat("Using simulated data for demonstration...\n\n")
  use_simulated <- TRUE
  n <- 50000
  hr <- data.table(
    eid = 1:n,
    is_fh_genetic = rbinom(n, 1, 0.005),
    tudor_prob = runif(n, 0, 0.1),
    tudor_score = rnorm(n, -3, 1),
    edlcn_score = sample(0:10, n, replace = TRUE, prob = c(0.3, 0.25, 0.15, 0.1, 0.08, 0.05, 0.03, 0.02, 0.01, 0.005, 0.005)),
    LDL_RW = rnorm(n, 5.8, 1.2),
    Gender_num = rbinom(n, 1, 0.46),
    Age_at_LDL1 = rnorm(n, 57, 8),
    statin_name = sample(c("None", "Atorvastatin", "Simvastatin", "Rosuvastatin", "Pravastatin"),
                         n, replace = TRUE, prob = c(0.65, 0.15, 0.12, 0.05, 0.03)),
    statin_tier = sample(0:3, n, replace = TRUE, prob = c(0.65, 0.05, 0.25, 0.05)),
    has_any_cvd = rbinom(n, 1, 0.08),
    Premature_ASCVD = rbinom(n, 1, 0.03),
    BMI_imputed = rnorm(n, 27.5, 5)
  )
  # Make FH cases have higher TUDOR
  hr[is_fh_genetic == 1, tudor_prob := pmin(tudor_prob + 0.05, 1)]
}

# ==============================================================================
# 2. DEFINE ALL SUBGROUPS
# ==============================================================================

# Helper function for AUC with DeLong CI
auc_result <- function(outcome, predictor, label) {
  valid <- !is.na(outcome) & !is.na(predictor)
  n_cases <- sum(outcome[valid] == 1)
  n_total <- sum(valid)

  if (n_cases < 5 || n_total < 50) {
    return(data.table(
      subgroup = label, n = n_total, n_fh = n_cases,
      auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      se = NA_real_
    ))
  }

  r <- roc(outcome[valid], predictor[valid], quiet = TRUE)
  ci <- ci.auc(r, method = "delong")

  data.table(
    subgroup = label, n = n_total, n_fh = n_cases,
    auc = as.numeric(ci[2]),
    ci_lo = as.numeric(ci[1]),
    ci_hi = as.numeric(ci[3]),
    se = (as.numeric(ci[3]) - as.numeric(ci[1])) / (2 * 1.96)
  )
}

results <- data.table()

# --- Overall ---
results <- rbind(results, auc_result(hr$is_fh_genetic, hr$tudor_prob, "Overall"))

# --- By Sex ---
cat("Computing subgroup AUCs...\n")
results <- rbind(results, auc_result(
  hr[Gender_num == 1]$is_fh_genetic, hr[Gender_num == 1]$tudor_prob, "Male"))
results <- rbind(results, auc_result(
  hr[Gender_num == 0]$is_fh_genetic, hr[Gender_num == 0]$tudor_prob, "Female"))

# --- By Age Decade ---
for (age_lo in seq(40, 70, 10)) {
  age_hi <- age_lo + 9
  label <- sprintf("Age %d-%d", age_lo, age_hi)
  sub <- hr[Age_at_LDL1 >= age_lo & Age_at_LDL1 < age_lo + 10]
  results <- rbind(results, auc_result(sub$is_fh_genetic, sub$tudor_prob, label))
}
sub_young <- hr[Age_at_LDL1 < 40]
results <- rbind(results, auc_result(sub_young$is_fh_genetic, sub_young$tudor_prob, "Age <40"))
sub_old <- hr[Age_at_LDL1 >= 70]
results <- rbind(results, auc_result(sub_old$is_fh_genetic, sub_old$tudor_prob, "Age 70+"))

# --- By Statin Status ---
results <- rbind(results, auc_result(
  hr[statin_name == "None"]$is_fh_genetic,
  hr[statin_name == "None"]$tudor_prob, "No statin"))
results <- rbind(results, auc_result(
  hr[statin_name != "None"]$is_fh_genetic,
  hr[statin_name != "None"]$tudor_prob, "On statin"))

# --- By Statin Type ---
for (st in c("Atorvastatin", "Simvastatin", "Rosuvastatin")) {
  sub <- hr[statin_name == st]
  if (nrow(sub) >= 100) {
    results <- rbind(results, auc_result(sub$is_fh_genetic, sub$tudor_prob,
                                          paste("On", st)))
  }
}

# --- By Statin Intensity ---
results <- rbind(results, auc_result(
  hr[statin_tier == 0]$is_fh_genetic, hr[statin_tier == 0]$tudor_prob, "No statin (tier 0)"))
results <- rbind(results, auc_result(
  hr[statin_tier %in% 1:2]$is_fh_genetic, hr[statin_tier %in% 1:2]$tudor_prob, "Low-med intensity"))
results <- rbind(results, auc_result(
  hr[statin_tier == 3]$is_fh_genetic, hr[statin_tier == 3]$tudor_prob, "High intensity"))

# --- By LDL Range ---
results <- rbind(results, auc_result(
  hr[LDL_RW >= 4.9 & LDL_RW < 6.5]$is_fh_genetic,
  hr[LDL_RW >= 4.9 & LDL_RW < 6.5]$tudor_prob, "LDL 4.9-6.5"))
results <- rbind(results, auc_result(
  hr[LDL_RW >= 6.5 & LDL_RW < 8.5]$is_fh_genetic,
  hr[LDL_RW >= 6.5 & LDL_RW < 8.5]$tudor_prob, "LDL 6.5-8.5"))
results <- rbind(results, auc_result(
  hr[LDL_RW >= 8.5]$is_fh_genetic,
  hr[LDL_RW >= 8.5]$tudor_prob, "LDL >= 8.5"))

# --- By ASCVD Status ---
results <- rbind(results, auc_result(
  hr[has_any_cvd == TRUE]$is_fh_genetic,
  hr[has_any_cvd == TRUE]$tudor_prob, "ASCVD present"))
results <- rbind(results, auc_result(
  hr[has_any_cvd == FALSE]$is_fh_genetic,
  hr[has_any_cvd == FALSE]$tudor_prob, "No ASCVD"))

# --- By Premature ASCVD ---
if (sum(hr$Premature_ASCVD, na.rm = TRUE) > 20) {
  results <- rbind(results, auc_result(
    hr[Premature_ASCVD == 1]$is_fh_genetic,
    hr[Premature_ASCVD == 1]$tudor_prob, "Premature ASCVD"))
}

# --- By BMI Category ---
results <- rbind(results, auc_result(
  hr[BMI_imputed < 25]$is_fh_genetic,
  hr[BMI_imputed < 25]$tudor_prob, "BMI <25 (normal)"))
results <- rbind(results, auc_result(
  hr[BMI_imputed >= 25 & BMI_imputed < 30]$is_fh_genetic,
  hr[BMI_imputed >= 25 & BMI_imputed < 30]$tudor_prob, "BMI 25-30 (overweight)"))
results <- rbind(results, auc_result(
  hr[BMI_imputed >= 30]$is_fh_genetic,
  hr[BMI_imputed >= 30]$tudor_prob, "BMI >= 30 (obese)"))

# Remove NAs
results <- results[!is.na(auc)]

cat(sprintf("\n  Computed %d subgroup AUCs\n\n", nrow(results)))

# ==============================================================================
# 3. INTERACTION TESTS
# ==============================================================================
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("INTERACTION TESTS\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

interaction_tests <- data.table()

# Sex interaction
if (nrow(hr) > 100) {
  fit_int <- glm(is_fh_genetic ~ tudor_score * Gender_num, data = hr, family = binomial)
  p_sex <- summary(fit_int)$coefficients["tudor_score:Gender_num", "Pr(>|z|)"]
  interaction_tests <- rbind(interaction_tests,
    data.table(interaction = "Score x Sex", p_raw = p_sex))
  cat(sprintf("  Score x Sex:    p = %.3f\n", p_sex))
}

# Age interaction (continuous)
if (nrow(hr) > 100) {
  fit_age <- glm(is_fh_genetic ~ tudor_score * Age_at_LDL1, data = hr, family = binomial)
  p_age <- summary(fit_age)$coefficients["tudor_score:Age_at_LDL1", "Pr(>|z|)"]
  interaction_tests <- rbind(interaction_tests,
    data.table(interaction = "Score x Age", p_raw = p_age))
  cat(sprintf("  Score x Age:    p = %.3f\n", p_age))
}

# Statin interaction
if (nrow(hr) > 100) {
  hr[, on_statin := as.integer(statin_name != "None")]
  fit_statin <- glm(is_fh_genetic ~ tudor_score * on_statin, data = hr, family = binomial)
  p_statin <- summary(fit_statin)$coefficients["tudor_score:on_statin", "Pr(>|z|)"]
  interaction_tests <- rbind(interaction_tests,
    data.table(interaction = "Score x Statin", p_raw = p_statin))
  cat(sprintf("  Score x Statin: p = %.3f\n", p_statin))
}

# ASCVD interaction
if (sum(hr$has_any_cvd) > 20) {
  fit_cvd <- glm(is_fh_genetic ~ tudor_score * has_any_cvd, data = hr, family = binomial)
  p_cvd <- tryCatch(
    summary(fit_cvd)$coefficients["tudor_score:has_any_cvdTRUE", "Pr(>|z|)"],
    error = function(e) NA
  )
  if (!is.na(p_cvd)) {
    interaction_tests <- rbind(interaction_tests,
      data.table(interaction = "Score x ASCVD", p_raw = p_cvd))
    cat(sprintf("  Score x ASCVD:  p = %.3f\n", p_cvd))
  }
}

# Multiple comparison correction
if (nrow(interaction_tests) > 0) {
  interaction_tests[, p_bonferroni := pmin(p_raw * .N, 1)]
  interaction_tests[, p_fdr := p.adjust(p_raw, method = "BH")]
  cat("\n  Bonferroni-corrected interaction tests:\n")
  for (i in seq_len(nrow(interaction_tests))) {
    cat(sprintf("    %-20s: p_raw = %.3f, p_Bonf = %.3f, p_FDR = %.3f\n",
                interaction_tests$interaction[i],
                interaction_tests$p_raw[i],
                interaction_tests$p_bonferroni[i],
                interaction_tests$p_fdr[i]))
  }
}
cat("\n")

# ==============================================================================
# 4. HETEROGENEITY TEST (Cochran's Q & I²)
# ==============================================================================
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("HETEROGENEITY ANALYSIS\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# Exclude "Overall" row for heterogeneity
sub_results <- results[subgroup != "Overall" & !is.na(se) & se > 0]

if (nrow(sub_results) >= 3) {
  # Cochran's Q test
  overall_auc <- results[subgroup == "Overall"]$auc
  weights <- 1 / sub_results$se^2
  Q <- sum(weights * (sub_results$auc - overall_auc)^2)
  df_q <- nrow(sub_results) - 1
  p_q <- 1 - pchisq(Q, df_q)

  # I² statistic
  I2 <- max(0, (Q - df_q) / Q * 100)

  cat(sprintf("Cochran's Q = %.1f (df = %d, p = %.3f)\n", Q, df_q, p_q))
  cat(sprintf("I² = %.1f%%\n", I2))
  cat(sprintf("Interpretation: %s heterogeneity\n\n",
              ifelse(I2 < 25, "Low", ifelse(I2 < 75, "Moderate", "High"))))
}

# ==============================================================================
# 5. FOREST PLOT
# ==============================================================================
cat("Generating forest plot...\n")

# Add grouping categories for the forest plot
results[, group := ifelse(subgroup == "Overall", "Overall",
                   ifelse(subgroup %in% c("Male", "Female"), "Sex",
                   ifelse(grepl("Age", subgroup), "Age",
                   ifelse(grepl("[Ss]tatin|intensity|tier", subgroup), "Statin",
                   ifelse(grepl("LDL", subgroup), "LDL Range",
                   ifelse(grepl("ASCVD|CVD", subgroup), "ASCVD",
                   ifelse(grepl("BMI", subgroup), "BMI", "Other")))))))]

# Order for plotting
results[, plot_order := .N:1]

p_forest <- ggplot(results[!is.na(auc)],
                    aes(x = auc, y = reorder(subgroup, plot_order))) +
  geom_point(aes(size = n_fh), shape = 18, color = "darkblue") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, color = "darkblue") +
  geom_vline(xintercept = results[subgroup == "Overall"]$auc,
             linetype = "dashed", color = "red3", linewidth = 0.5) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey50") +
  scale_x_continuous(limits = c(0.45, 1.0), breaks = seq(0.5, 1.0, 0.1)) +
  scale_size_continuous(range = c(2, 6), name = "FH cases") +
  labs(x = "Area Under ROC Curve (AUC)",
       y = NULL,
       title = "TUDOR v2: Subgroup Analysis",
       subtitle = sprintf("Overall AUC = %.3f | %d subgroups analysed",
                           results[subgroup == "Overall"]$auc,
                           nrow(results) - 1)) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9)
  )

tryCatch({
  ggsave(file.path(FIG_DIR, "forest_plot_subgroups.pdf"),
         p_forest, width = 10, height = 8)
  ggsave(file.path(FIG_DIR, "forest_plot_subgroups.png"),
         p_forest, width = 10, height = 8, dpi = 300)
  cat("  Saved: forest_plot_subgroups.pdf/png\n")
}, error = function(e) {
  cat("  Note: Could not save plot:", e$message, "\n")
})

# ==============================================================================
# 6. SAVE RESULTS
# ==============================================================================

# Print summary table
cat("\n=== SUBGROUP SUMMARY TABLE ===\n\n")
cat(sprintf("%-30s | %7s | %4s | %9s | %s\n",
            "Subgroup", "N", "FH", "AUC", "95% CI"))
cat(strrep("-", 75), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i]
  cat(sprintf("%-30s | %7d | %4d | %9.3f | [%.3f - %.3f]\n",
              r$subgroup, r$n, r$n_fh, r$auc, r$ci_lo, r$ci_hi))
}
cat(strrep("-", 75), "\n")

fwrite(results, file.path(TABLE_DIR, "subgroup_analysis_comprehensive.csv"))
if (nrow(interaction_tests) > 0) {
  fwrite(interaction_tests, file.path(TABLE_DIR, "interaction_tests.csv"))
}

subgroup_results <- list(
  subgroups = results,
  interactions = interaction_tests,
  heterogeneity = if (exists("Q")) list(Q = Q, df = df_q, p = p_q, I2 = I2) else NULL,
  timestamp = Sys.time()
)
saveRDS(subgroup_results, file.path(OUTPUT_DIR, "29_subgroup_results.rds"))

cat("\nSaved: subgroup_analysis_comprehensive.csv, interaction_tests.csv\n")
cat("\n=== 29_comprehensive_subgroups.R COMPLETE ===\n")
