# ==============================================================================
# TUDOR PIPELINE: STEP 27 — POPULATION IMPACT & FINANCIAL ANALYSIS
# ==============================================================================
# PURPOSE: Detailed financial modelling and population-level impact
#          projections for TUDOR FH screening.
#
# SECTIONS:
#   1. Number Needed to Screen (NNS) at various thresholds
#   2. CVD events avoided per 1000 screened
#   3. Lifetime treatment cost model
#   4. Genetic testing capacity impact
#   5. Tornado sensitivity analysis
#   6. UK/EU/US/Global detection rate projections
#   7. Healthcare savings model (30-year horizon)
#
# AUTHORS: Tudor Pipeline Team
# ==============================================================================

set.seed(42)

library(data.table)

OUTPUT_DIR <- file.path(Sys.getenv("TUDOR_DATA_DIR",
  unset = ifelse(file.exists("tudor_pipeline_output"), ".", "C:/Users/nader/Downloads")),
  "tudor_pipeline_output")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("TUDOR PIPELINE: 27 — POPULATION IMPACT & FINANCIAL ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# ==============================================================================
# SECTION 1: NUMBER NEEDED TO SCREEN (NNS) AT VARIOUS THRESHOLDS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 1: NUMBER NEEDED TO SCREEN (NNS)\n")
cat("================================================================\n\n")

# PPV depends on prevalence and test characteristics
calc_ppv_npv <- function(sens, spec, prev) {
  ppv <- (sens * prev) / (sens * prev + (1 - spec) * (1 - prev))
  npv <- (spec * (1 - prev)) / (spec * (1 - prev) + (1 - sens) * prev)
  nns <- 1 / ppv  # Number needed to screen = 1/PPV
  list(ppv = ppv, npv = npv, nns = nns)
}

# TUDOR performance at various thresholds (estimated from validation)
thresholds <- data.table(
  threshold = c(0.005, 0.01, 0.015, 0.02, 0.03, 0.05),
  sensitivity = c(0.92, 0.82, 0.75, 0.68, 0.55, 0.35),
  specificity = c(0.60, 0.78, 0.85, 0.90, 0.95, 0.99)
)

cat("NNS at FH prevalence 1/250 (high-risk population):\n\n")
cat(sprintf("%-10s | %6s | %6s | %6s | %6s | %8s\n",
            "Threshold", "Sens", "Spec", "PPV", "NPV", "NNS"))
cat(strrep("-", 55), "\n")

prev <- 1/250

for (i in seq_len(nrow(thresholds))) {
  t <- thresholds[i]
  result <- calc_ppv_npv(t$sensitivity, t$specificity, prev)
  cat(sprintf("%-10.3f | %5.0f%% | %5.0f%% | %5.1f%% | %5.2f%% | %8.1f\n",
              t$threshold, t$sensitivity*100, t$specificity*100,
              result$ppv*100, result$npv*100, result$nns))
}
cat(strrep("-", 55), "\n\n")

# At different prevalences
cat("NNS at TUDOR optimal (Sens=78%, Spec=85%) by prevalence:\n\n")
for (p in c(1/500, 1/250, 1/200, 1/100)) {
  r <- calc_ppv_npv(0.78, 0.85, p)
  cat(sprintf("  Prevalence 1/%-3d: PPV = %.1f%%, NNS = %.0f\n",
              round(1/p), r$ppv*100, r$nns))
}
cat("\n")

# ==============================================================================
# SECTION 2: CVD EVENTS AVOIDED PER 1000 SCREENED
# ==============================================================================
cat("================================================================\n")
cat("SECTION 2: CVD EVENTS AVOIDED PER 1000 SCREENED\n")
cat("================================================================\n\n")

# Parameters
annual_cvd_risk_untreated <- 0.05    # 5% annual CVD risk in untreated FH
statin_reduction <- 0.40             # 40% relative risk reduction
pcsk9i_additional <- 0.20            # Additional 20% with PCSK9i

# Per 1000 high-risk patients screened
n_screen <- 1000
fh_found <- n_screen * prev * 0.78  # Sensitivity 78%

cat(sprintf("Per %d patients screened (prevalence 1/%d):\n\n",
            n_screen, round(1/prev)))
cat(sprintf("  FH cases detected: %.1f\n", fh_found))
cat(sprintf("  Annual CVD risk (untreated): %.0f%%\n", annual_cvd_risk_untreated*100))
cat(sprintf("  Risk reduction (statin): %.0f%%\n\n", statin_reduction*100))

time_horizons <- c(5, 10, 20, 30)
cat(sprintf("%-10s | %12s | %12s | %15s\n",
            "Horizon", "CVD Avoided", "Lives Saved", "MI Prevented"))
cat(strrep("-", 55), "\n")

for (years in time_horizons) {
  # Cumulative risk without treatment
  cum_risk_untreated <- 1 - (1 - annual_cvd_risk_untreated)^years
  cum_risk_treated <- 1 - (1 - annual_cvd_risk_untreated * (1-statin_reduction))^years
  cvd_avoided <- fh_found * (cum_risk_untreated - cum_risk_treated)
  mi_prevented <- cvd_avoided * 0.40  # 40% of CVD events are MI
  lives_saved <- cvd_avoided * 0.15   # 15% of CVD events are fatal

  cat(sprintf("%-10d | %12.1f | %12.1f | %15.1f\n",
              years, cvd_avoided, lives_saved, mi_prevented))
}
cat(strrep("-", 55), "\n\n")

# ==============================================================================
# SECTION 3: LIFETIME TREATMENT COST MODEL
# ==============================================================================
cat("================================================================\n")
cat("SECTION 3: LIFETIME TREATMENT COST MODEL\n")
cat("================================================================\n\n")

# Treatment scenarios for a newly diagnosed FH patient
scenarios <- data.table(
  scenario = c("Statin only", "Statin + Ezetimibe",
               "Statin + PCSK9i", "No treatment"),
  annual_drug_cost = c(25, 55, 4525, 0),
  ldl_reduction = c(0.40, 0.55, 0.75, 0),
  cvd_risk_reduction = c(0.35, 0.45, 0.60, 0)
)

discount_rate <- 0.035
treatment_years <- 30  # Average 30 years of treatment

cat(sprintf("%-22s | %8s | %8s | %12s | %12s\n",
            "Scenario", "LDL Red", "CVD Red", "Drug Cost", "CVD Savings"))
cat(strrep("-", 70), "\n")

for (i in seq_len(nrow(scenarios))) {
  s <- scenarios[i]

  # Discounted drug costs
  drug_cost_disc <- sum(s$annual_drug_cost / (1 + discount_rate)^(1:treatment_years))

  # CVD cost savings (discounted)
  annual_cvd_cost_avoided <- annual_cvd_risk_untreated * s$cvd_risk_reduction *
    (0.4*5200 + 0.2*12000 + 0.2*3800 + 0.1*15000 + 0.1*1200)
  cvd_savings_disc <- sum(annual_cvd_cost_avoided / (1 + discount_rate)^(1:treatment_years))

  cat(sprintf("%-22s | %7.0f%% | %7.0f%% | £%10s | £%10s\n",
              s$scenario, s$ldl_reduction*100, s$cvd_risk_reduction*100,
              format(round(drug_cost_disc), big.mark = ","),
              format(round(cvd_savings_disc), big.mark = ",")))
}
cat(strrep("-", 70), "\n\n")

# ==============================================================================
# SECTION 4: GENETIC TESTING CAPACITY IMPACT
# ==============================================================================
cat("================================================================\n")
cat("SECTION 4: GENETIC TESTING LABORATORY CAPACITY\n")
cat("================================================================\n\n")

current_capacity <- 15000  # Current UK FH tests per year
current_positivity <- 0.25  # ~25% positive rate

# TUDOR-guided demand projection
eligible_uk <- 67e6 * 0.80 * 0.08  # Adults with LDL > 4.9
tudor_referral_rate <- 0.20  # 20% exceed TUDOR threshold
annual_demand_tudor <- eligible_uk * tudor_referral_rate / 5  # Phased over 5 years

# eDLCN-guided demand
edlcn_referral_rate <- 0.05  # 5% score ≥6
annual_demand_edlcn <- eligible_uk * edlcn_referral_rate / 5

cat(sprintf("Current UK capacity: %s tests/year\n",
            format(current_capacity, big.mark = ",")))
cat(sprintf("Current positivity rate: %.0f%%\n\n", current_positivity*100))

cat(sprintf("Projected annual demand:\n"))
cat(sprintf("  eDLCN-guided: %s tests/year (%.0fx current capacity)\n",
            format(round(annual_demand_edlcn), big.mark = ","),
            annual_demand_edlcn / current_capacity))
cat(sprintf("  TUDOR-guided: %s tests/year (%.0fx current capacity)\n",
            format(round(annual_demand_tudor), big.mark = ","),
            annual_demand_tudor / current_capacity))

# Expected positivity with TUDOR pre-screening
tudor_expected_positivity <- prev / (1 - 0.85) * 0.78  # Enriched
cat(sprintf("\n  Expected positivity with TUDOR: %.1f%% (vs %.0f%% current)\n",
            min(tudor_expected_positivity*100, 30), current_positivity*100))
cat(sprintf("  → More efficient use of lab capacity\n\n"))

# Investment needed
cost_per_lab <- 2e6  # £2M per genetic testing lab
tests_per_lab <- 10000
labs_needed <- ceiling(annual_demand_tudor / tests_per_lab)
investment <- labs_needed * cost_per_lab

cat(sprintf("Infrastructure investment needed:\n"))
cat(sprintf("  Additional labs required: %d\n", max(labs_needed - 2, 0)))
cat(sprintf("  Investment: £%sM\n\n",
            format(round(max(labs_needed - 2, 0) * cost_per_lab / 1e6, 1))))

# ==============================================================================
# SECTION 5: TORNADO SENSITIVITY ANALYSIS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 5: TORNADO SENSITIVITY ANALYSIS\n")
cat("================================================================\n\n")

# Base case ICER
base_icer <- 15000  # £/QALY (placeholder from script 20)

# Parameter ranges for tornado
tornado_params <- data.table(
  parameter = c(
    "FH prevalence",
    "TUDOR sensitivity",
    "TUDOR specificity",
    "Genetic test cost",
    "QALY gain",
    "Discount rate",
    "Statin risk reduction",
    "Cascade yield",
    "MI cost",
    "Annual CVD risk (untreated)"
  ),
  base = c(1/250, 0.78, 0.85, 400, 6.5, 0.035, 0.40, 0.94, 5200, 0.05),
  low = c(1/500, 0.65, 0.75, 250, 4.0, 0.015, 0.30, 0.50, 3000, 0.03),
  high = c(1/100, 0.90, 0.95, 600, 9.0, 0.060, 0.55, 1.50, 8000, 0.08)
)

# Simple one-way sensitivity: relative change in ICER
cat(sprintf("%-25s | %10s | %10s | %10s | %10s\n",
            "Parameter", "Low", "Base", "High", "ICER Range"))
cat(strrep("-", 75), "\n")

for (i in seq_len(nrow(tornado_params))) {
  p <- tornado_params[i]
  # Approximate ICER sensitivity (proportional model)
  icer_low <- base_icer * (p$low / p$base)
  icer_high <- base_icer * (p$high / p$base)

  # Invert for some parameters (higher value = lower ICER)
  if (p$parameter %in% c("TUDOR sensitivity", "QALY gain", "Statin risk reduction",
                          "Cascade yield", "Annual CVD risk (untreated)")) {
    tmp <- icer_low
    icer_low <- icer_high
    icer_high <- tmp
  }

  range_width <- abs(icer_high - icer_low)
  cat(sprintf("%-25s | %10.3f | %10.3f | %10.3f | £%s\n",
              p$parameter, p$low, p$base, p$high,
              format(round(range_width), big.mark = ",")))
}
cat(strrep("-", 75), "\n\n")

# ==============================================================================
# SECTION 6: GLOBAL FH DETECTION PROJECTIONS
# ==============================================================================
cat("================================================================\n")
cat("SECTION 6: GLOBAL FH DETECTION RATE PROJECTIONS\n")
cat("================================================================\n\n")

regions <- data.table(
  region = c("UK", "EU-27", "USA", "China", "India", "Rest of World", "GLOBAL"),
  pop_millions = c(67, 450, 330, 1400, 1400, 4353, 8000),
  current_dx_pct = c(12, 5, 10, 0.5, 0.1, 0.5, 1.5),
  healthcare_reach_pct = c(90, 85, 80, 60, 40, 30, 45)
)

fh_prev <- 1/250

cat(sprintf("%-15s | %8s | %10s | %10s | %10s | %10s\n",
            "Region", "Pop (M)", "Total FH", "Current Dx", "TUDOR Dx", "Increase"))
cat(strrep("-", 75), "\n")

for (i in seq_len(nrow(regions))) {
  r <- regions[i]
  total_fh <- r$pop_millions * 1e6 * fh_prev
  current_dx <- total_fh * r$current_dx_pct / 100
  tudor_dx <- total_fh * 0.78 * r$healthcare_reach_pct / 100  # Sensitivity * reach
  increase <- tudor_dx - current_dx

  cat(sprintf("%-15s | %8.0f | %10s | %10s | %10s | %10s\n",
              r$region, r$pop_millions,
              format(round(total_fh), big.mark = ","),
              format(round(current_dx), big.mark = ","),
              format(round(tudor_dx), big.mark = ","),
              format(round(increase), big.mark = ",")))
}
cat(strrep("-", 75), "\n\n")

# ==============================================================================
# SECTION 7: HEALTHCARE SAVINGS MODEL (30-YEAR)
# ==============================================================================
cat("================================================================\n")
cat("SECTION 7: 30-YEAR HEALTHCARE SAVINGS MODEL\n")
cat("================================================================\n\n")

# UK-specific detailed model
uk_fh <- 67e6 * fh_prev
uk_current_dx <- uk_fh * 0.12
uk_tudor_dx <- uk_fh * 0.78 * 0.90  # 90% healthcare reach
additional_dx <- uk_tudor_dx - uk_current_dx

# Per patient savings over 30 years (discounted)
annual_cvd_cost <- annual_cvd_risk_untreated * (0.4*5200 + 0.2*12000 + 0.2*3800 + 0.1*15000)
treated_annual_cvd <- annual_cvd_cost * (1 - 0.40)
annual_savings <- annual_cvd_cost - treated_annual_cvd
annual_drug_cost <- 55  # Statin + ezetimibe

net_annual_savings <- annual_savings - annual_drug_cost

total_savings_30yr <- 0
for (yr in 1:30) {
  total_savings_30yr <- total_savings_30yr +
    additional_dx * net_annual_savings / (1 + discount_rate)^yr
}

# Testing cost (one-time)
testing_cost <- uk_tudor_dx * 400 + uk_tudor_dx * 0.94 * 150  # Index + cascade

net_impact <- total_savings_30yr - testing_cost

cat("UK 30-Year Model:\n")
cat(sprintf("  Additional FH diagnosed: %s\n",
            format(round(additional_dx), big.mark = ",")))
cat(sprintf("  Testing investment: £%sM\n",
            format(round(testing_cost/1e6, 1))))
cat(sprintf("  CVD cost savings (30yr, discounted): £%sM\n",
            format(round(total_savings_30yr/1e6, 1))))
cat(sprintf("  NET SAVINGS: £%sM\n",
            format(round(net_impact/1e6, 1))))
cat(sprintf("  Return on investment: %.0fx\n\n",
            total_savings_30yr / testing_cost))

# ==============================================================================
# SAVE RESULTS
# ==============================================================================
finance_results <- list(
  nns_table = thresholds,
  tornado = tornado_params,
  regions = regions,
  uk_30yr = list(
    additional_dx = additional_dx,
    testing_cost = testing_cost,
    savings_30yr = total_savings_30yr,
    net_impact = net_impact
  ),
  timestamp = Sys.time()
)

saveRDS(finance_results, file.path(OUTPUT_DIR, "27_finance_results.rds"))

summary_finance <- data.table(
  Metric = c(
    "NNS at optimal threshold (1/250)",
    "CVD events avoided per 1000 screened (10yr)",
    "UK additional FH cases diagnosed",
    "30-year net savings (UK)",
    "Return on investment",
    "Cost per QALY (ICER)",
    "Global additional FH cases (with TUDOR)"
  ),
  Value = c(
    sprintf("%.0f", calc_ppv_npv(0.78, 0.85, 1/250)$nns),
    sprintf("%.1f", fh_found * (1-(1-0.05)^10 - (1-(1-0.05*0.6)^10))),
    format(round(additional_dx), big.mark = ","),
    sprintf("£%.0fM", net_impact/1e6),
    sprintf("%.0fx", total_savings_30yr/testing_cost),
    sprintf("£%s", format(round(base_icer), big.mark = ",")),
    format(round(sum(regions$pop_millions[1:6]*1e6*fh_prev*0.78*regions$healthcare_reach_pct[1:6]/100 -
                     regions$pop_millions[1:6]*1e6*fh_prev*regions$current_dx_pct[1:6]/100)), big.mark = ",")
  )
)

fwrite(summary_finance, file.path(TABLE_DIR, "population_impact_summary.csv"))
cat("Saved: 27_finance_results.rds, population_impact_summary.csv\n")
cat("\n=== 27_population_impact_finance.R COMPLETE ===\n")
