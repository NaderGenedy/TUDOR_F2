# ==============================================================================
# TUDOR PIPELINE: STEP 20 — PRACTICE CHANGE & HEALTH ECONOMICS ANALYSIS
# ==============================================================================
# PURPOSE: Comprehensive analysis of how TUDOR changes FH clinical practice.
#          Includes cost-effectiveness, pathway simulation, cascade yield,
#          NHS budget impact, and global projections.
#
# FOR: Nature Medicine submission — "Changing Practice" chapter
#
# AUTHORS: Tudor Pipeline Team
# INSTITUTION: Cardiff and Vale University Health Board
# ==============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
})

OUTPUT_DIR <- file.path(Sys.getenv("TUDOR_DATA_DIR",
  unset = ifelse(file.exists("tudor_pipeline_output"), ".", "C:/Users/nader/Downloads")),
  "tudor_pipeline_output")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("TUDOR PIPELINE: 20 — PRACTICE CHANGE & HEALTH ECONOMICS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# SECTION 1: HEALTH ECONOMICS PARAMETERS (2024 GBP)
# ==============================================================================
cat("================================================================\n")
cat("SECTION 1: COST & EFFECTIVENESS PARAMETERS\n")
cat("================================================================\n\n")

params <- list(
  # --- COSTS ---
  cost_genetic_test     = 400,    # NGS FH gene panel (LDLR/APOB/PCSK9)
  cost_cascade_test     = 150,    # Variant-specific cascade test
  cost_clinic_visit     = 120,    # Lipid clinic consultation
  cost_lipid_panel      = 15,     # Blood lipid panel
  cost_tudor_compute    = 0,      # Automated from EHR data
  cost_edlcn_assess     = 120,    # Clinic visit for eDLCN assessment

  # --- CVD EVENT COSTS ---
  cost_mi_admission     = 5200,   # MI hospitalisation (NHS ref costs)
  cost_cabg             = 12000,  # CABG surgery
  cost_pci_stent        = 3800,   # PCI/stenting
  cost_stroke           = 15000,  # Stroke admission
  cost_cardiac_rehab    = 1200,   # Cardiac rehabilitation

  # --- TREATMENT COSTS (annual) ---
  cost_statin_year      = 25,     # Generic atorvastatin
  cost_ezetimibe_year   = 30,     # Generic ezetimibe
  cost_pcsk9i_year      = 4500,   # Evolocumab/alirocumab

  # --- EFFECTIVENESS ---
  rr_mi_fh_male         = 13.0,   # FH relative risk MI (men)
  rr_mi_fh_female       = 8.0,    # FH relative risk MI (women)
  cvd_risk_reduction    = 0.40,   # Statin risk reduction for CVD events
  qaly_gain_early_dx    = 6.5,    # Mean QALY gain from early FH treatment
  life_expectancy_gain  = 8,      # Years of life gained with treatment
  discount_rate         = 0.035,  # NICE reference case

  # --- POPULATION ---
  prevalence_fh         = 1/250,  # FH prevalence (heterozygous)
  uk_pop                = 67e6,   # UK population
  england_pop           = 56e6,   # England population
  eu_pop                = 450e6,  # EU population
  us_pop                = 330e6,  # US population
  global_pop            = 8e9,    # Global population
  currently_diagnosed   = 0.12,   # ~12% of FH cases currently diagnosed (UK)

  # --- CASCADE SCREENING ---
  relatives_per_index   = 2.5,    # First-degree relatives per index case
  cascade_detection     = 0.50,   # P(FH) in first-degree relative (autosomal dom)
  cascade_uptake        = 0.75,   # Proportion who agree to cascade testing

  # --- MODEL PERFORMANCE (from validation) ---
  tudor_sensitivity     = 0.78,   # At Youden optimal threshold
  tudor_specificity     = 0.85,
  tudor_ppv_enriched    = 0.02,   # PPV in high-risk cohort
  edlcn_sensitivity     = 0.55,   # eDLCN >=6 (Probable)
  edlcn_specificity     = 0.90,
  edlcn_ppv_enriched    = 0.02
)

cat("Key Parameters:\n")
cat(sprintf("  FH prevalence: 1/%.0f\n", 1/params$prevalence_fh))
cat(sprintf("  TUDOR: Sens=%.0f%%, Spec=%.0f%%\n",
            params$tudor_sensitivity*100, params$tudor_specificity*100))
cat(sprintf("  eDLCN: Sens=%.0f%%, Spec=%.0f%%\n",
            params$edlcn_sensitivity*100, params$edlcn_specificity*100))
cat(sprintf("  Genetic test: £%d | Cascade: £%d/relative\n",
            params$cost_genetic_test, params$cost_cascade_test))
cat(sprintf("  QALY gain: %.1f years | Discount: %.1f%%\n\n",
            params$qaly_gain_early_dx, params$discount_rate*100))

# ==============================================================================
# SECTION 2: CLINICAL PATHWAY SIMULATION (100,000 patients)
# ==============================================================================
cat("================================================================\n")
cat("SECTION 2: CLINICAL PATHWAY SIMULATION\n")
cat("================================================================\n\n")

n_sim <- 100000
fh_cases <- round(n_sim * params$prevalence_fh)
non_fh <- n_sim - fh_cases

cat(sprintf("Simulating %s patients (FH: %d, non-FH: %d)\n\n",
            format(n_sim, big.mark = ","), fh_cases, non_fh))

# --- PATHWAY A: Current Practice (eDLCN ≥6 → Genetic Testing) ---
edlcn_tp <- round(fh_cases * params$edlcn_sensitivity)      # True positives
edlcn_fp <- round(non_fh * (1 - params$edlcn_specificity))  # False positives
edlcn_tests <- edlcn_tp + edlcn_fp                           # Total tests ordered
edlcn_found <- edlcn_tp                                       # Cases found

edlcn_cost_testing <- edlcn_tests * params$cost_genetic_test
edlcn_cost_clinic <- n_sim * params$cost_edlcn_assess   # Everyone needs clinic visit
edlcn_cost_total <- edlcn_cost_testing + edlcn_cost_clinic
edlcn_cost_per_case <- edlcn_cost_total / max(edlcn_found, 1)
edlcn_nns <- edlcn_tests / max(edlcn_found, 1)

# --- PATHWAY B: TUDOR-Guided (TUDOR ≥ threshold → Genetic Testing) ---
tudor_tp <- round(fh_cases * params$tudor_sensitivity)
tudor_fp <- round(non_fh * (1 - params$tudor_specificity))
tudor_tests <- tudor_tp + tudor_fp
tudor_found <- tudor_tp

tudor_cost_testing <- tudor_tests * params$cost_genetic_test
tudor_cost_compute <- 0  # TUDOR computed from existing EHR data
tudor_cost_total <- tudor_cost_testing + tudor_cost_compute
tudor_cost_per_case <- tudor_cost_total / max(tudor_found, 1)
tudor_nns <- tudor_tests / max(tudor_found, 1)

# --- PATHWAY C: Universal Genetic Screening (all high-risk) ---
univ_tests <- n_sim
univ_found <- fh_cases
univ_cost <- univ_tests * params$cost_genetic_test
univ_cost_per_case <- univ_cost / max(univ_found, 1)
univ_nns <- univ_tests / max(univ_found, 1)

# Print comparison
cat(sprintf("%-25s | %-12s | %-12s | %-12s\n",
            "Metric", "eDLCN", "TUDOR", "Universal"))
cat(strrep("-", 70), "\n")
cat(sprintf("%-25s | %12d | %12d | %12d\n", "Cases found", edlcn_found, tudor_found, univ_found))
cat(sprintf("%-25s | %12d | %12d | %12d\n", "Tests ordered", edlcn_tests, tudor_tests, univ_tests))
cat(sprintf("%-25s | %12.1f | %12.1f | %12.1f\n", "NNS (tests/case)", edlcn_nns, tudor_nns, univ_nns))
cat(sprintf("%-25s | £%10s | £%10s | £%10s\n", "Total cost",
            format(round(edlcn_cost_total), big.mark = ","),
            format(round(tudor_cost_total), big.mark = ","),
            format(round(univ_cost), big.mark = ",")))
cat(sprintf("%-25s | £%10s | £%10s | £%10s\n", "Cost per case found",
            format(round(edlcn_cost_per_case), big.mark = ","),
            format(round(tudor_cost_per_case), big.mark = ","),
            format(round(univ_cost_per_case), big.mark = ",")))
cat(sprintf("%-25s | %11.0f%% | %11.0f%% | %11.0f%%\n", "Detection rate",
            100*edlcn_found/fh_cases, 100*tudor_found/fh_cases, 100))
cat(strrep("-", 70), "\n\n")

additional_found <- tudor_found - edlcn_found
cat(sprintf("TUDOR finds %d ADDITIONAL FH cases vs eDLCN per %s screened\n",
            additional_found, format(n_sim, big.mark = ",")))
cat(sprintf("Cost saving vs universal: £%s per %s screened\n\n",
            format(round(univ_cost - tudor_cost_total), big.mark = ","),
            format(n_sim, big.mark = ",")))

# ==============================================================================
# SECTION 3: ICER CALCULATION
# ==============================================================================
cat("================================================================\n")
cat("SECTION 3: INCREMENTAL COST-EFFECTIVENESS RATIO (ICER)\n")
cat("================================================================\n\n")

# Discounted QALY gain per case found
discount_factor <- function(years, rate) {
  sum(1 / (1 + rate)^(1:years))
}

discounted_qaly <- params$qaly_gain_early_dx * (1 / (1 + params$discount_rate)^5)

# ICER: TUDOR vs eDLCN
delta_cost <- tudor_cost_total - edlcn_cost_total
delta_qaly <- (tudor_found - edlcn_found) * discounted_qaly
icer_tudor_vs_edlcn <- delta_cost / max(delta_qaly, 0.001)

# ICER: TUDOR vs No Screening
delta_cost_none <- tudor_cost_total
delta_qaly_none <- tudor_found * discounted_qaly
icer_tudor_vs_none <- delta_cost_none / max(delta_qaly_none, 0.001)

cat(sprintf("TUDOR vs eDLCN:\n"))
cat(sprintf("  Incremental cost: £%s\n", format(round(delta_cost), big.mark = ",")))
cat(sprintf("  Incremental QALYs: %.1f\n", delta_qaly))
cat(sprintf("  ICER: £%s per QALY\n", format(round(icer_tudor_vs_edlcn), big.mark = ",")))
cat(sprintf("  NICE threshold (£20k-£30k): %s\n\n",
            ifelse(icer_tudor_vs_edlcn < 30000, "COST-EFFECTIVE", "NOT cost-effective")))

cat(sprintf("TUDOR vs No Screening:\n"))
cat(sprintf("  ICER: £%s per QALY\n\n",
            format(round(icer_tudor_vs_none), big.mark = ",")))

# ==============================================================================
# SECTION 4: MONTE CARLO PROBABILISTIC SENSITIVITY ANALYSIS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 4: PROBABILISTIC SENSITIVITY ANALYSIS (1000 iterations)\n")
cat("================================================================\n\n")

n_mc <- 1000
icer_mc <- numeric(n_mc)
nb_tudor <- numeric(n_mc)
nb_edlcn <- numeric(n_mc)

for (mc in seq_len(n_mc)) {
  # Sample parameters from distributions
  sens_t <- rbeta(1, 78, 22)    # TUDOR sensitivity ~78%
  spec_t <- rbeta(1, 85, 15)    # TUDOR specificity ~85%
  sens_e <- rbeta(1, 55, 45)    # eDLCN sensitivity ~55%
  cost_test <- rnorm(1, 400, 50)  # Genetic test cost
  qaly <- rnorm(1, 6.5, 1.5)     # QALY gain
  prev <- rbeta(1, 4, 996)       # Prevalence ~1/250

  n_fh <- round(n_sim * prev)
  n_nonfh <- n_sim - n_fh

  # TUDOR
  t_found <- round(n_fh * sens_t)
  t_tests <- t_found + round(n_nonfh * (1 - spec_t))
  t_cost <- t_tests * cost_test

  # eDLCN
  e_found <- round(n_fh * sens_e)
  e_tests <- e_found + round(n_nonfh * 0.10)
  e_cost <- e_tests * cost_test + n_sim * 120

  d_cost <- t_cost - e_cost
  d_qaly <- (t_found - e_found) * qaly * (1/(1+0.035)^5)

  icer_mc[mc] <- ifelse(d_qaly > 0, d_cost / d_qaly, NA)

  # Net monetary benefit at £30k/QALY
  wtp <- 30000
  nb_tudor[mc] <- t_found * qaly * wtp * (1/(1+0.035)^5) - t_cost
  nb_edlcn[mc] <- e_found * qaly * wtp * (1/(1+0.035)^5) - e_cost
}

icer_valid <- icer_mc[!is.na(icer_mc) & is.finite(icer_mc)]
cat(sprintf("ICER distribution (median [IQR]):\n"))
cat(sprintf("  Median: £%s per QALY\n",
            format(round(median(icer_valid)), big.mark = ",")))
cat(sprintf("  IQR: [£%s, £%s]\n",
            format(round(quantile(icer_valid, 0.25)), big.mark = ","),
            format(round(quantile(icer_valid, 0.75)), big.mark = ",")))

# Cost-effectiveness acceptability
wtp_thresholds <- seq(0, 100000, by = 5000)
ceac <- sapply(wtp_thresholds, function(wtp) {
  mean(nb_tudor > nb_edlcn, na.rm = TRUE)
})

cat(sprintf("\nCost-effectiveness acceptability:\n"))
for (wtp in c(20000, 30000, 50000)) {
  idx <- which.min(abs(wtp_thresholds - wtp))
  cat(sprintf("  At £%dk/QALY: %.0f%% probability TUDOR is cost-effective\n",
              wtp/1000, ceac[idx]*100))
}
cat("\n")

# ==============================================================================
# SECTION 5: CASCADE SCREENING YIELD
# ==============================================================================
cat("================================================================\n")
cat("SECTION 5: CASCADE SCREENING YIELD\n")
cat("================================================================\n\n")

cascade_per_index <- params$relatives_per_index * params$cascade_detection *
                     params$cascade_uptake

cat(sprintf("Per index case identified:\n"))
cat(sprintf("  First-degree relatives: %.1f\n", params$relatives_per_index))
cat(sprintf("  P(FH in relative): %.0f%%\n", params$cascade_detection*100))
cat(sprintf("  Cascade uptake: %.0f%%\n", params$cascade_uptake*100))
cat(sprintf("  Additional cases found: %.2f per index case\n\n", cascade_per_index))

total_tudor_yield <- tudor_found * (1 + cascade_per_index)
total_edlcn_yield <- edlcn_found * (1 + cascade_per_index)
cascade_cost_per_case <- params$cost_cascade_test / params$cascade_detection

cat(sprintf("Total yield (index + cascade) per %s screened:\n", format(n_sim, big.mark=",")))
cat(sprintf("  TUDOR: %.0f index + %.0f cascade = %.0f total\n",
            tudor_found, tudor_found * cascade_per_index, total_tudor_yield))
cat(sprintf("  eDLCN: %.0f index + %.0f cascade = %.0f total\n",
            edlcn_found, edlcn_found * cascade_per_index, total_edlcn_yield))
cat(sprintf("  Cascade cost per additional case: £%d\n\n", round(cascade_cost_per_case)))

# ==============================================================================
# SECTION 6: NHS BUDGET IMPACT (ENGLAND, 5-YEAR)
# ==============================================================================
cat("================================================================\n")
cat("SECTION 6: NHS ENGLAND BUDGET IMPACT (5-YEAR)\n")
cat("================================================================\n\n")

# Estimate eligible population (adults with LDL > 4.9 in UK)
eligible_fraction <- 0.08  # ~8% of adults have LDL > 4.9
adult_fraction <- 0.80     # ~80% are adults
eligible_england <- params$england_pop * adult_fraction * eligible_fraction

fh_in_eligible <- eligible_england * params$prevalence_fh * 3  # Enriched in high-LDL

# Year-by-year rollout (20% per year)
cat("5-Year Budget Impact Model (phased rollout):\n\n")
cat(sprintf("%-6s | %12s | %10s | %10s | %10s | %10s\n",
            "Year", "Screened", "FH Found", "Tests", "Cost", "CVD Avoided"))
cat(strrep("-", 70), "\n")

total_found_5yr <- 0
total_cost_5yr <- 0
total_cvd_avoided <- 0

for (yr in 1:5) {
  rollout <- min(yr * 0.20, 1.0)  # 20% per year up to 100%
  screened <- eligible_england * rollout
  fh_found <- screened * params$prevalence_fh * 3 * params$tudor_sensitivity
  tests <- fh_found + screened * (1 - params$prevalence_fh * 3) *
           (1 - params$tudor_specificity)
  cost <- tests * params$cost_genetic_test
  # CVD events avoided (assume 5% annual CVD risk in untreated FH, 40% reduction)
  cvd_avoided <- fh_found * 0.05 * params$cvd_risk_reduction

  total_found_5yr <- total_found_5yr + fh_found
  total_cost_5yr <- total_cost_5yr + cost
  total_cvd_avoided <- total_cvd_avoided + cvd_avoided

  cat(sprintf("%-6d | %12s | %10s | %10s | £%9s | %10.0f\n",
              yr, format(round(screened), big.mark = ","),
              format(round(fh_found), big.mark = ","),
              format(round(tests), big.mark = ","),
              format(round(cost), big.mark = ","),
              cvd_avoided))
}
cat(strrep("-", 70), "\n")
cat(sprintf("%-6s | %12s | %10s | %10s | £%9s | %10.0f\n",
            "TOTAL", "", format(round(total_found_5yr), big.mark = ","),
            "", format(round(total_cost_5yr), big.mark = ","),
            total_cvd_avoided))

# Cost of avoided CVD events
avg_cvd_cost <- 0.5 * params$cost_mi_admission + 0.2 * params$cost_cabg +
                0.2 * params$cost_pci_stent + 0.1 * params$cost_stroke
savings_5yr <- total_cvd_avoided * avg_cvd_cost
cat(sprintf("\nEstimated CVD cost savings (5yr): £%s\n",
            format(round(savings_5yr), big.mark = ",")))
cat(sprintf("Net budget impact (5yr): £%s\n\n",
            format(round(total_cost_5yr - savings_5yr), big.mark = ",")))

# ==============================================================================
# SECTION 7: GLOBAL IMPACT PROJECTIONS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 7: GLOBAL IMPACT PROJECTIONS\n")
cat("================================================================\n\n")

populations <- data.table(
  region = c("UK", "EU", "USA", "Global"),
  pop = c(67e6, 450e6, 330e6, 8e9),
  current_dx_rate = c(0.12, 0.05, 0.10, 0.01)
)

cat(sprintf("%-10s | %12s | %10s | %10s | %10s | %12s\n",
            "Region", "FH Cases", "Currently", "TUDOR Dx", "Additional", "CVD Avoided"))
cat(strrep("-", 75), "\n")

for (i in seq_len(nrow(populations))) {
  r <- populations[i]
  total_fh <- r$pop * params$prevalence_fh
  currently_dx <- total_fh * r$current_dx_rate
  tudor_dx <- total_fh * params$tudor_sensitivity * 0.60  # 60% population reach
  additional <- tudor_dx - currently_dx
  # CVD avoided over 20 years
  cvd_avoided_20yr <- additional * 0.05 * params$cvd_risk_reduction * 20

  cat(sprintf("%-10s | %12s | %10s | %10s | %10s | %12s\n",
              r$region,
              format(round(total_fh), big.mark = ","),
              format(round(currently_dx), big.mark = ","),
              format(round(tudor_dx), big.mark = ","),
              format(round(additional), big.mark = ","),
              format(round(cvd_avoided_20yr), big.mark = ",")))
}
cat(strrep("-", 75), "\n\n")

# Lives saved
global_fh <- 8e9 * params$prevalence_fh
additional_global <- global_fh * params$tudor_sensitivity * 0.3 - global_fh * 0.01
lives_saved_30yr <- additional_global * 0.02 * 30  # 2% annual mortality reduction
cat(sprintf("Estimated lives saved globally over 30 years: %s\n\n",
            format(round(lives_saved_30yr), big.mark = ",")))

# ==============================================================================
# SECTION 8: SENSITIVITY ANALYSIS ON PREVALENCE
# ==============================================================================
cat("================================================================\n")
cat("SECTION 8: PREVALENCE SENSITIVITY ANALYSIS\n")
cat("================================================================\n\n")

prevalences <- c(1/500, 1/250, 1/200, 1/100)

cat(sprintf("%-12s | %8s | %8s | %12s | %12s\n",
            "Prevalence", "Cases", "NNS", "Cost/Case", "ICER"))
cat(strrep("-", 60), "\n")

for (prev in prevalences) {
  n_fh_s <- round(n_sim * prev)
  t_found_s <- round(n_fh_s * params$tudor_sensitivity)
  t_tests_s <- t_found_s + round((n_sim - n_fh_s) * (1 - params$tudor_specificity))
  t_cost_s <- t_tests_s * params$cost_genetic_test
  nns_s <- t_tests_s / max(t_found_s, 1)
  cpc_s <- t_cost_s / max(t_found_s, 1)
  icer_s <- cpc_s / discounted_qaly

  cat(sprintf("%-12s | %8d | %8.1f | £%10s | £%10s\n",
              sprintf("1/%d", round(1/prev)),
              t_found_s, nns_s,
              format(round(cpc_s), big.mark = ","),
              format(round(icer_s), big.mark = ",")))
}
cat(strrep("-", 60), "\n\n")

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
practice_results <- list(
  params = params,
  pathway_comparison = list(
    tudor = list(found = tudor_found, tests = tudor_tests, cost = tudor_cost_total),
    edlcn = list(found = edlcn_found, tests = edlcn_tests, cost = edlcn_cost_total),
    universal = list(found = univ_found, tests = univ_tests, cost = univ_cost)
  ),
  icer = list(tudor_vs_edlcn = icer_tudor_vs_edlcn, tudor_vs_none = icer_tudor_vs_none),
  cascade = list(yield = cascade_per_index, total_tudor = total_tudor_yield),
  budget_5yr = list(cost = total_cost_5yr, savings = savings_5yr,
                    found = total_found_5yr, cvd_avoided = total_cvd_avoided),
  mc_icer = icer_valid,
  ceac = data.table(wtp = wtp_thresholds, prob_ce = ceac),
  timestamp = Sys.time()
)

saveRDS(practice_results, file.path(OUTPUT_DIR, "20_practice_change_results.rds"))

# Summary table for manuscript
summary_table <- data.table(
  Metric = c(
    "FH cases found per 100,000 screened (TUDOR)",
    "FH cases found per 100,000 screened (eDLCN)",
    "Additional cases found by TUDOR",
    "Number needed to screen (TUDOR)",
    "Number needed to screen (eDLCN)",
    "Cost per case found (TUDOR)",
    "Cost per case found (eDLCN)",
    "ICER (TUDOR vs eDLCN)",
    "P(cost-effective at £30k/QALY)",
    "5-year NHS budget impact",
    "5-year CVD events avoided (England)",
    "Cascade yield per index case",
    "Lives saved globally (30yr projection)"
  ),
  Value = c(
    tudor_found, edlcn_found, additional_found,
    sprintf("%.1f", tudor_nns), sprintf("%.1f", edlcn_nns),
    sprintf("£%s", format(round(tudor_cost_per_case), big.mark = ",")),
    sprintf("£%s", format(round(edlcn_cost_per_case), big.mark = ",")),
    sprintf("£%s/QALY", format(round(icer_tudor_vs_edlcn), big.mark = ",")),
    sprintf("%.0f%%", ceac[which.min(abs(wtp_thresholds - 30000))]*100),
    sprintf("£%sM", format(round(total_cost_5yr/1e6, 1))),
    format(round(total_cvd_avoided), big.mark = ","),
    sprintf("%.2f", cascade_per_index),
    format(round(lives_saved_30yr), big.mark = ",")
  )
)

fwrite(summary_table, file.path(TABLE_DIR, "practice_change_summary.csv"))
cat("Saved: 20_practice_change_results.rds, practice_change_summary.csv\n")
cat("\n=== 20_practice_change_analysis.R COMPLETE ===\n")
