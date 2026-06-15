# ==============================================================================
# 18_tudor_c_grey_zone.R
# TUDOR-C: ApoB/LDL-C + Lp(a) Augmented Model for Grey Zone Reclassification
# ==============================================================================
#
# Authors: Nader Genedy, Soha Zouwail
# Institution: Cardiff and Vale University Health Board
#
# PURPOSE: Develop and evaluate TUDOR-C, a grey-zone-specific augmented model
#          combining TUDOR base probability with ApoB/LDL-C ratio and Lp(a)
#          for reclassification of patients in the intermediate probability
#          zone (25-75%). Uses the CLEANED UKB cohort (secondary causes excluded).
#
# INPUTS:
#   1. tudor_v2_workspace.RData (base TUDOR predictions + lipid data)
#   2. tudor_analysis_clean_no_secondary.rds (cleaned cohort, if available)
#
# OUTPUTS:
#   Tables:
#     - tudor_c_model_results.csv
#     - tudor_c_grey_zone_reclassification.csv
#     - tudor_c_nri_idi.csv
#     - tudor_c_apob_ratio_grey_zone.csv
#   Figures:
#     - fig_tudor_c_grey_zone_roc.pdf/png
#     - fig_tudor_c_reclassification.pdf/png
#     - fig_tudor_c_apob_ratio_distribution.pdf/png
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  TUDOR-C: Grey Zone Augmentation Model\n")
cat("  Script 18 | ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
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

fig_dir   <- "C:/Users/nader/Downloads/tudor_pipeline_output/figures"
table_dir <- "C:/Users/nader/Downloads/tudor_pipeline_output/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

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
# PART 2: LOAD DATA
# ============================================================

cat("Loading workspace...\n")
ws_paths <- c(
  "C:/Users/nader/Downloads/tudor_v2_workspace.RData",
  "C:/Users/nader/Downloads/tudor_pipeline/tudor_v2_workspace.RData",
  "tudor_v2_workspace.RData"
)
ws_loaded <- FALSE
for (wp in ws_paths) {
  if (file.exists(wp)) { load(wp); cat(sprintf("  Loaded: %s\n", wp)); ws_loaded <- TRUE; break }
}
if (!ws_loaded) stop("Cannot find workspace")

# Try to load cleaned cohort
clean_paths <- c(
  "C:/Users/nader/Downloads/tudor_pipeline_output/tudor_analysis_clean_no_secondary.rds",
  "tudor_pipeline_output/tudor_analysis_clean_no_secondary.rds"
)
clean_loaded <- FALSE
for (cp in clean_paths) {
  if (file.exists(cp)) {
    clean_data <- readRDS(cp)
    cat(sprintf("  Loaded cleaned cohort: %s (%d rows)\n", cp, nrow(clean_data)))
    clean_loaded <- TRUE
    break
  }
}

# ============================================================
# PART 3: IDENTIFY DATA
# ============================================================

cat("\nIdentifying data objects...\n")

# Find UKB object
ukb_obj <- NULL
for (nm in c("ukb", "ukb_full", "ukb_df", "df_ukb")) {
  if (exists(nm)) { ukb_obj <- get(nm); break }
}
if (is.null(ukb_obj)) stop("No UKB data found")
ukb_df <- as.data.frame(ukb_obj)
cat(sprintf("  UKB data: %d rows x %d cols\n", nrow(ukb_df), ncol(ukb_df)))

# Find lipid clinic subset
lc_obj <- NULL
for (nm in c("lc", "lipid_clinic", "ukb_lc")) {
  if (exists(nm)) { lc_obj <- as.data.frame(get(nm)); break }
}

# Key columns
fh_col <- NULL
for (fc in c("fh_genetic", "FH", "fh", "genetic_fh")) {
  if (fc %in% names(ukb_df)) { fh_col <- fc; break }
}

tudor_col <- NULL
for (tc in c("tudor_prob", "TUDOR_prob", "tudor_score", "pred", "predicted", "prob")) {
  if (tc %in% names(ukb_df)) { tudor_col <- tc; break }
}

apob_col <- NULL
for (ac in c("ApoB", "apob", "apoB", "ApoB.1", "p30640_i0")) {
  if (ac %in% names(ukb_df)) { apob_col <- ac; break }
}

ldl_col <- NULL
for (lc_name in c("LDL_untreated", "LDL_UT", "ldl_ut", "LDL.1")) {
  if (lc_name %in% names(ukb_df)) { ldl_col <- lc_name; break }
}

lpa_col <- NULL
for (lpc in c("Lpa", "lpa", "Lp_a", "p30790_i0")) {
  if (lpc %in% names(ukb_df)) { lpa_col <- lpc; break }
}

gene_col <- NULL
for (gc in c("gene", "Gene", "fh_gene", "gene_group")) {
  if (gc %in% names(ukb_df)) { gene_col <- gc; break }
}

cat(sprintf("  FH col: %s | TUDOR col: %s\n", fh_col, tudor_col))
cat(sprintf("  ApoB col: %s | LDL col: %s | Lp(a) col: %s\n", apob_col, ldl_col, lpa_col))

# ============================================================
# PART 4: PREPARE ANALYSIS DATASET
# ============================================================

cat("\n--- Part 4: Preparing Analysis Dataset ---\n")

# Use lipid clinic subset if available, otherwise construct it
if (!is.null(lc_obj) && !is.null(fh_col) && fh_col %in% names(lc_obj)) {
  analysis_df <- lc_obj
  cat(sprintf("  Using pre-existing lipid clinic object: %d rows\n", nrow(analysis_df)))
} else {
  analysis_df <- ukb_df
  cat(sprintf("  Using full UKB: %d rows\n", nrow(analysis_df)))
}

# If cleaned cohort available, filter analysis dataset
if (clean_loaded && "eid" %in% names(clean_data)) {
  # Get eids to keep
  eid_col <- NULL
  for (ec in c("eid", "participant.eid", "f.eid")) {
    if (ec %in% names(analysis_df)) { eid_col <- ec; break }
  }
  if (!is.null(eid_col)) {
    clean_eids <- clean_data$eid
    n_before <- nrow(analysis_df)
    analysis_df <- analysis_df[analysis_df[[eid_col]] %in% clean_eids, ]
    cat(sprintf("  Applied secondary cause exclusion: %d -> %d (-%d)\n",
                n_before, nrow(analysis_df), n_before - nrow(analysis_df)))
  }
}

# Compute ApoB/LDL-C ratio
if (!is.null(apob_col) && !is.null(ldl_col)) {
  analysis_df$apob_ldl_ratio <- analysis_df[[apob_col]] / analysis_df[[ldl_col]]
  # Remove infinite/extreme values
  analysis_df$apob_ldl_ratio[!is.finite(analysis_df$apob_ldl_ratio)] <- NA
  analysis_df$apob_ldl_ratio[analysis_df$apob_ldl_ratio > 2 | analysis_df$apob_ldl_ratio < 0] <- NA
  cat(sprintf("  ApoB/LDL-C ratio: %d non-missing\n",
              sum(!is.na(analysis_df$apob_ldl_ratio))))
}

# ============================================================
# PART 5: DEFINE GREY ZONE
# ============================================================

cat("\n--- Part 5: Grey Zone Definition ---\n")

if (!is.null(tudor_col) && !is.null(fh_col)) {
  analysis_df$tudor_tier <- ifelse(analysis_df[[tudor_col]] < 0.25, "Low",
                             ifelse(analysis_df[[tudor_col]] <= 0.75, "Intermediate", "High"))

  tier_tbl <- table(analysis_df$tudor_tier, analysis_df[[fh_col]])
  cat("  TUDOR Tier Distribution:\n")
  print(tier_tbl)

  # Grey zone = intermediate (25-75%)
  grey_zone <- analysis_df[analysis_df$tudor_tier == "Intermediate" &
                             !is.na(analysis_df[[tudor_col]]) &
                             !is.na(analysis_df[[fh_col]]), ]
  cat(sprintf("\n  Grey zone patients: %d\n", nrow(grey_zone)))
  cat(sprintf("  Grey zone FH+: %d (%.1f%%)\n",
              sum(grey_zone[[fh_col]] == 1),
              100 * mean(grey_zone[[fh_col]] == 1)))
}

# ============================================================
# PART 6: TUDOR-C MODEL — GREY ZONE AUGMENTATION
# ============================================================

cat("\n--- Part 6: TUDOR-C Model Development ---\n")

model_results <- list()

# --- Model 0: Base TUDOR in grey zone ---
if (nrow(grey_zone) > 0 && sum(grey_zone[[fh_col]] == 1) >= 10) {
  roc_base <- roc(grey_zone[[fh_col]], grey_zone[[tudor_col]], quiet = TRUE)
  ci_base  <- ci.auc(roc_base)
  cat(sprintf("  Base TUDOR in grey zone: AUC=%.3f (%.3f-%.3f)\n",
              as.numeric(auc(roc_base)), ci_base[1], ci_base[3]))
  model_results[["Base TUDOR"]] <- data.frame(
    Model = "Base TUDOR", AUC = as.numeric(auc(roc_base)),
    CI_lower = ci_base[1], CI_upper = ci_base[3],
    N = nrow(grey_zone), N_FH = sum(grey_zone[[fh_col]] == 1),
    DeLong_p = NA
  )
}

# --- Model A: ApoB/LDL-C ratio alone in grey zone ---
if ("apob_ldl_ratio" %in% names(grey_zone)) {
  gz_apob <- grey_zone[!is.na(grey_zone$apob_ldl_ratio), ]
  if (nrow(gz_apob) > 50 && sum(gz_apob[[fh_col]] == 1) >= 10) {
    roc_ratio <- roc(gz_apob[[fh_col]], gz_apob$apob_ldl_ratio, quiet = TRUE,
                     direction = ">")  # Lower ratio = MORE likely FH
    ci_ratio  <- ci.auc(roc_ratio)
    cat(sprintf("  ApoB/LDL-C ratio alone in grey zone: AUC=%.3f (%.3f-%.3f)\n",
                as.numeric(auc(roc_ratio)), ci_ratio[1], ci_ratio[3]))
    model_results[["ApoB/LDL ratio"]] <- data.frame(
      Model = "ApoB/LDL-C ratio alone", AUC = as.numeric(auc(roc_ratio)),
      CI_lower = ci_ratio[1], CI_upper = ci_ratio[3],
      N = nrow(gz_apob), N_FH = sum(gz_apob[[fh_col]] == 1),
      DeLong_p = NA
    )
  }
}

# --- Model B: TUDOR + ApoB/LDL-C (logistic combination) ---
if ("apob_ldl_ratio" %in% names(grey_zone)) {
  gz_complete <- grey_zone[!is.na(grey_zone$apob_ldl_ratio) &
                             !is.na(grey_zone[[tudor_col]]), ]
  if (nrow(gz_complete) > 50 && sum(gz_complete[[fh_col]] == 1) >= 10) {
    fit_b <- glm(as.formula(paste(fh_col, "~ tudor_col_val + apob_ldl_ratio")),
                 data = data.frame(
                   fh = gz_complete[[fh_col]],
                   tudor_col_val = gz_complete[[tudor_col]],
                   apob_ldl_ratio = gz_complete$apob_ldl_ratio
                 ) |> setNames(c(fh_col, "tudor_col_val", "apob_ldl_ratio")),
                 family = binomial)

    gz_complete$tudor_c_b <- predict(fit_b, newdata = data.frame(
      tudor_col_val = gz_complete[[tudor_col]],
      apob_ldl_ratio = gz_complete$apob_ldl_ratio
    ), type = "response")

    roc_b <- roc(gz_complete[[fh_col]], gz_complete$tudor_c_b, quiet = TRUE)
    ci_b  <- ci.auc(roc_b)

    # DeLong vs base
    roc_base_sub <- roc(gz_complete[[fh_col]], gz_complete[[tudor_col]], quiet = TRUE)
    delong_b <- tryCatch({
      dt <- roc.test(roc_base_sub, roc_b, method = "delong")
      dt$p.value
    }, error = function(e) NA)

    cat(sprintf("  TUDOR + ApoB/LDL-C: AUC=%.3f (%.3f-%.3f), DeLong p=%s\n",
                as.numeric(auc(roc_b)), ci_b[1], ci_b[3],
                ifelse(is.na(delong_b), "NA", formatC(delong_b, format = "e", digits = 2))))

    model_results[["TUDOR + ApoB/LDL"]] <- data.frame(
      Model = "TUDOR-C (TUDOR + ApoB/LDL-C)", AUC = as.numeric(auc(roc_b)),
      CI_lower = ci_b[1], CI_upper = ci_b[3],
      N = nrow(gz_complete), N_FH = sum(gz_complete[[fh_col]] == 1),
      DeLong_p = delong_b
    )
  }
}

# --- Model C: TUDOR + ApoB/LDL-C + Lp(a) ---
if (!is.null(lpa_col) && lpa_col %in% names(grey_zone)) {
  gz_lpa <- grey_zone[!is.na(grey_zone$apob_ldl_ratio) &
                        !is.na(grey_zone[[lpa_col]]) &
                        !is.na(grey_zone[[tudor_col]]), ]

  if (nrow(gz_lpa) > 50 && sum(gz_lpa[[fh_col]] == 1) >= 10) {
    # Log-transform Lp(a) for better distribution
    gz_lpa$log_lpa <- log(gz_lpa[[lpa_col]] + 1)

    fit_c <- glm(as.formula(paste(fh_col, "~ tudor_val + apob_ratio + log_lpa_val")),
                 data = data.frame(
                   fh = gz_lpa[[fh_col]],
                   tudor_val = gz_lpa[[tudor_col]],
                   apob_ratio = gz_lpa$apob_ldl_ratio,
                   log_lpa_val = gz_lpa$log_lpa
                 ) |> setNames(c(fh_col, "tudor_val", "apob_ratio", "log_lpa_val")),
                 family = binomial)

    gz_lpa$tudor_c_full <- predict(fit_c, newdata = data.frame(
      tudor_val = gz_lpa[[tudor_col]],
      apob_ratio = gz_lpa$apob_ldl_ratio,
      log_lpa_val = gz_lpa$log_lpa
    ), type = "response")

    roc_c <- roc(gz_lpa[[fh_col]], gz_lpa$tudor_c_full, quiet = TRUE)
    ci_c  <- ci.auc(roc_c)

    roc_base_lpa <- roc(gz_lpa[[fh_col]], gz_lpa[[tudor_col]], quiet = TRUE)
    delong_c <- tryCatch({
      dt <- roc.test(roc_base_lpa, roc_c, method = "delong")
      dt$p.value
    }, error = function(e) NA)

    cat(sprintf("  TUDOR-C full (+ ApoB/LDL + Lp(a)): AUC=%.3f (%.3f-%.3f), DeLong p=%s\n",
                as.numeric(auc(roc_c)), ci_c[1], ci_c[3],
                ifelse(is.na(delong_c), "NA", formatC(delong_c, format = "e", digits = 2))))

    model_results[["TUDOR-C full"]] <- data.frame(
      Model = "TUDOR-C full (TUDOR + ApoB/LDL-C + Lp(a))", AUC = as.numeric(auc(roc_c)),
      CI_lower = ci_c[1], CI_upper = ci_c[3],
      N = nrow(gz_lpa), N_FH = sum(gz_lpa[[fh_col]] == 1),
      DeLong_p = delong_c
    )

    # Print model coefficients
    cat("\n  TUDOR-C full model coefficients:\n")
    print(summary(fit_c)$coefficients)
  } else {
    cat(sprintf("  Lp(a) available but insufficient grey zone data (n=%d, FH+=%d)\n",
                nrow(gz_lpa), sum(gz_lpa[[fh_col]] == 1)))
  }
} else {
  cat("  Lp(a) column not found — TUDOR-C full model skipped\n")
}

# Save model results
if (length(model_results) > 0) {
  mr_df <- do.call(rbind, model_results)
  write.csv(mr_df, file.path(table_dir, "tudor_c_model_results.csv"), row.names = FALSE)
  cat("\n  Model comparison:\n")
  print(mr_df)
}

# ============================================================
# PART 7: APOB/LDL-C RATIO ANALYSIS IN GREY ZONE
# ============================================================

cat("\n--- Part 7: ApoB/LDL-C Ratio Grey Zone Analysis ---\n")

if ("apob_ldl_ratio" %in% names(grey_zone) && !is.null(fh_col)) {
  gz_r <- grey_zone[!is.na(grey_zone$apob_ldl_ratio), ]

  fh_pos <- gz_r[gz_r[[fh_col]] == 1, ]
  fh_neg <- gz_r[gz_r[[fh_col]] == 0, ]

  cat(sprintf("  Grey zone with ApoB/LDL-C data: %d (FH+: %d, FH-: %d)\n",
              nrow(gz_r), nrow(fh_pos), nrow(fh_neg)))

  # Ratio statistics
  cat(sprintf("  FH+ ratio: mean=%.4f (SD=%.4f), median=%.4f\n",
              mean(fh_pos$apob_ldl_ratio), sd(fh_pos$apob_ldl_ratio),
              median(fh_pos$apob_ldl_ratio)))
  cat(sprintf("  FH- ratio: mean=%.4f (SD=%.4f), median=%.4f\n",
              mean(fh_neg$apob_ldl_ratio), sd(fh_neg$apob_ldl_ratio),
              median(fh_neg$apob_ldl_ratio)))

  # Wilcoxon test
  wt <- wilcox.test(fh_pos$apob_ldl_ratio, fh_neg$apob_ldl_ratio)
  cat(sprintf("  Wilcoxon p=%s\n", formatC(wt$p.value, format = "e", digits = 2)))

  # Threshold analysis: ratio >= 0.31 (from Genedy & Zouwail 2025)
  cat("\n  --- Threshold >=0.31 (Genedy & Zouwail, ref #39) ---\n")
  gz_r$ratio_high <- ifelse(gz_r$apob_ldl_ratio >= 0.31, 1, 0)

  tbl_thresh <- table(gz_r$ratio_high, gz_r[[fh_col]])
  if (ncol(tbl_thresh) == 2 && nrow(tbl_thresh) == 2) {
    # Note: In FH, ratio is LOWER (cholesterol-enriched particles)
    # So ratio < 0.31 might be more specific for FH
    # Let's compute both ways
    tp <- tbl_thresh["0", "1"]  # ratio < 0.31 and FH+ (correct for FH)
    fp <- tbl_thresh["0", "0"]
    fn <- tbl_thresh["1", "1"]
    tn <- tbl_thresh["1", "0"]

    sens_low <- tp / (tp + fn)
    spec_low <- tn / (tn + fp)
    ppv_low  <- tp / (tp + fp)
    npv_low  <- tn / (tn + fn)

    cat(sprintf("  Ratio < 0.31 for FH (cholesterol-enriched):\n"))
    cat(sprintf("    Sensitivity: %.1f%%\n", 100 * sens_low))
    cat(sprintf("    Specificity: %.1f%%\n", 100 * spec_low))
    cat(sprintf("    PPV: %.1f%%, NPV: %.1f%%\n", 100 * ppv_low, 100 * npv_low))

    # Also try ratio >= 0.31 as risk marker (per previous paper — higher ratio = more events)
    tp2 <- tbl_thresh["1", "1"]
    fp2 <- tbl_thresh["1", "0"]
    fn2 <- tbl_thresh["0", "1"]
    tn2 <- tbl_thresh["0", "0"]

    sens_high <- tp2 / (tp2 + fn2)
    spec_high <- tn2 / (tn2 + fp2)

    cat(sprintf("\n  Ratio >= 0.31 (per Genedy & Zouwail ASCVD marker):\n"))
    cat(sprintf("    Sensitivity: %.1f%%\n", 100 * sens_high))
    cat(sprintf("    Specificity: %.1f%%\n", 100 * spec_high))
  }

  # Save ratio analysis
  ratio_results <- data.frame(
    Metric = c("Grey_zone_N", "FH_positive", "FH_negative",
               "FH_pos_ratio_mean", "FH_neg_ratio_mean",
               "FH_pos_ratio_median", "FH_neg_ratio_median",
               "Wilcoxon_p", "Ratio_ge_031_N", "Ratio_lt_031_N"),
    Value = c(nrow(gz_r), nrow(fh_pos), nrow(fh_neg),
              round(mean(fh_pos$apob_ldl_ratio), 4),
              round(mean(fh_neg$apob_ldl_ratio), 4),
              round(median(fh_pos$apob_ldl_ratio), 4),
              round(median(fh_neg$apob_ldl_ratio), 4),
              wt$p.value,
              sum(gz_r$ratio_high == 1),
              sum(gz_r$ratio_high == 0))
  )
  write.csv(ratio_results, file.path(table_dir, "tudor_c_apob_ratio_grey_zone.csv"),
            row.names = FALSE)
}

# ============================================================
# PART 8: RECLASSIFICATION ANALYSIS
# ============================================================

cat("\n--- Part 8: Reclassification Analysis ---\n")

if (exists("gz_complete") && "tudor_c_b" %in% names(gz_complete)) {
  # Reclassify using TUDOR-C
  gz_complete$tudor_tier_base <- ifelse(gz_complete[[tudor_col]] < 0.25, "Low",
                                  ifelse(gz_complete[[tudor_col]] <= 0.75, "Intermediate", "High"))
  gz_complete$tudor_c_tier <- ifelse(gz_complete$tudor_c_b < 0.25, "Low",
                               ifelse(gz_complete$tudor_c_b <= 0.75, "Intermediate", "High"))

  reclass_tbl <- table(Base = gz_complete$tudor_tier_base,
                        TUDOR_C = gz_complete$tudor_c_tier,
                        FH = gz_complete[[fh_col]])
  cat("  Reclassification table:\n")
  print(reclass_tbl)

  # NRI calculation (categorical)
  # Events (FH+)
  events <- gz_complete[gz_complete[[fh_col]] == 1, ]
  events_up   <- sum(events$tudor_c_tier == "High" & events$tudor_tier_base == "Intermediate")
  events_down <- sum(events$tudor_c_tier == "Low" & events$tudor_tier_base == "Intermediate")
  nri_events  <- (events_up - events_down) / nrow(events)

  # Non-events (FH-)
  nonevents <- gz_complete[gz_complete[[fh_col]] == 0, ]
  ne_up   <- sum(nonevents$tudor_c_tier == "High" & nonevents$tudor_tier_base == "Intermediate")
  ne_down <- sum(nonevents$tudor_c_tier == "Low" & nonevents$tudor_tier_base == "Intermediate")
  nri_nonevents <- (ne_down - ne_up) / nrow(nonevents)

  nri_total <- nri_events + nri_nonevents

  cat(sprintf("\n  NRI (categorical):\n"))
  cat(sprintf("    Events correctly reclassified up: %d, down: %d (NRI_events=%.3f)\n",
              events_up, events_down, nri_events))
  cat(sprintf("    Non-events correctly reclassified down: %d, up: %d (NRI_nonevents=%.3f)\n",
              ne_down, ne_up, nri_nonevents))
  cat(sprintf("    Total NRI: %.3f\n", nri_total))

  # IDI
  mean_prob_events_base <- mean(events[[tudor_col]])
  mean_prob_events_c    <- mean(events$tudor_c_b)
  mean_prob_ne_base     <- mean(nonevents[[tudor_col]])
  mean_prob_ne_c        <- mean(nonevents$tudor_c_b)

  idi <- (mean_prob_events_c - mean_prob_events_base) -
         (mean_prob_ne_c - mean_prob_ne_base)
  cat(sprintf("    IDI: %.4f\n", idi))

  # Save
  nri_df <- data.frame(
    Metric = c("NRI_events", "NRI_nonevents", "NRI_total", "IDI",
               "Events_reclassified_up", "Events_reclassified_down",
               "NonEvents_reclassified_up", "NonEvents_reclassified_down"),
    Value = c(nri_events, nri_nonevents, nri_total, idi,
              events_up, events_down, ne_up, ne_down)
  )
  write.csv(nri_df, file.path(table_dir, "tudor_c_nri_idi.csv"), row.names = FALSE)

  # Reclassification summary
  reclass_summary <- data.frame(
    Category = c("Stayed_Intermediate", "Reclassified_High", "Reclassified_Low"),
    FH_positive = c(
      sum(events$tudor_c_tier == "Intermediate"),
      events_up, events_down
    ),
    FH_negative = c(
      sum(nonevents$tudor_c_tier == "Intermediate"),
      ne_up, ne_down
    )
  )
  write.csv(reclass_summary, file.path(table_dir, "tudor_c_grey_zone_reclassification.csv"),
            row.names = FALSE)
}

# ============================================================
# PART 9: FIGURES
# ============================================================

cat("\n--- Part 9: Generating Figures ---\n")

# --- Figure A: Grey Zone ROC Comparison ---
tryCatch({
  cat("  Figure A: Grey zone ROC comparison...\n")

  if (exists("gz_complete") && "tudor_c_b" %in% names(gz_complete)) {
    roc_base_gz <- roc(gz_complete[[fh_col]], gz_complete[[tudor_col]], quiet = TRUE)
    roc_c_gz    <- roc(gz_complete[[fh_col]], gz_complete$tudor_c_b, quiet = TRUE)

    roc_df <- rbind(
      data.frame(Model = sprintf("Base TUDOR (AUC=%.3f)", as.numeric(auc(roc_base_gz))),
                 sens = roc_base_gz$sensitivities, spec = 1 - roc_base_gz$specificities),
      data.frame(Model = sprintf("TUDOR-C (AUC=%.3f)", as.numeric(auc(roc_c_gz))),
                 sens = roc_c_gz$sensitivities, spec = 1 - roc_c_gz$specificities)
    )

    # Add Lp(a) model if available
    if (exists("gz_lpa") && "tudor_c_full" %in% names(gz_lpa)) {
      roc_full <- roc(gz_lpa[[fh_col]], gz_lpa$tudor_c_full, quiet = TRUE)
      roc_df <- rbind(roc_df,
        data.frame(Model = sprintf("TUDOR-C+Lp(a) (AUC=%.3f)", as.numeric(auc(roc_full))),
                   sens = roc_full$sensitivities, spec = 1 - roc_full$specificities))
    }

    p_roc <- ggplot(roc_df, aes(x = spec, y = sens, colour = Model)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
      geom_line(linewidth = 0.6) +
      scale_colour_manual(values = c("#2166AC", "#B2182B", "#4DAF4A")) +
      labs(x = "1 - Specificity", y = "Sensitivity",
           title = "TUDOR-C: Grey Zone Augmentation",
           subtitle = sprintf("Patients with intermediate TUDOR probability (25-75%%), n=%d",
                              nrow(gz_complete))) +
      theme_nature(base_size = 8) +
      theme(legend.position = c(0.7, 0.25),
            legend.title = element_blank())

    save_fig(p_roc, "fig_tudor_c_grey_zone_roc", w_mm = 120, h_mm = 110)
  }
}, error = function(e) cat(sprintf("  Figure A error: %s\n", e$message)))

# --- Figure B: ApoB/LDL-C Ratio Distribution ---
tryCatch({
  cat("  Figure B: ApoB/LDL-C ratio distribution in grey zone...\n")

  if (exists("gz_r") && nrow(gz_r) > 0) {
    gz_r$FH_status <- factor(ifelse(gz_r[[fh_col]] == 1, "FH+", "FH-"),
                              levels = c("FH-", "FH+"))

    p_dist <- ggplot(gz_r, aes(x = apob_ldl_ratio, fill = FH_status)) +
      geom_density(alpha = 0.5, colour = NA) +
      geom_vline(xintercept = 0.31, linetype = "dashed", colour = "red", linewidth = 0.4) +
      annotate("text", x = 0.315, y = Inf, label = "0.31 threshold",
               hjust = 0, vjust = 1.5, size = 2.5, colour = "red") +
      scale_fill_manual(values = c("FH-" = "#4393C3", "FH+" = "#D6604D"), name = "") +
      scale_x_continuous(limits = c(0, 0.8)) +
      labs(x = "ApoB / LDL-C Ratio (g/mmol)", y = "Density",
           title = "ApoB/LDL-C Ratio Distribution in TUDOR Grey Zone",
           subtitle = sprintf("FH+ (n=%d) vs FH- (n=%d) | Wilcoxon p=%s",
                              nrow(fh_pos), nrow(fh_neg),
                              formatC(wt$p.value, format = "e", digits = 2))) +
      theme_nature(base_size = 8) +
      theme(legend.position = c(0.85, 0.85))

    save_fig(p_dist, "fig_tudor_c_apob_ratio_distribution", w_mm = 140, h_mm = 100)
  }
}, error = function(e) cat(sprintf("  Figure B error: %s\n", e$message)))

# --- Figure C: Reclassification Waterfall ---
tryCatch({
  cat("  Figure C: Reclassification summary...\n")

  if (exists("reclass_summary") && nrow(reclass_summary) > 0) {
    rc_long <- data.frame(
      Category = rep(reclass_summary$Category, 2),
      FH_Status = rep(c("FH+", "FH-"), each = nrow(reclass_summary)),
      Count = c(reclass_summary$FH_positive, reclass_summary$FH_negative)
    )

    rc_long$Category <- factor(rc_long$Category,
                                levels = c("Reclassified_High", "Stayed_Intermediate", "Reclassified_Low"))

    p_reclass <- ggplot(rc_long, aes(x = Category, y = Count, fill = FH_Status)) +
      geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
      geom_text(aes(label = Count), position = position_dodge(0.7), vjust = -0.3, size = 2.5) +
      scale_fill_manual(values = c("FH+" = "#D6604D", "FH-" = "#4393C3"), name = "") +
      scale_x_discrete(labels = c("Reclassified_High" = "Reclassified\nHigh Risk",
                                   "Stayed_Intermediate" = "Remained\nIntermediate",
                                   "Reclassified_Low" = "Reclassified\nLow Risk")) +
      labs(x = "", y = "Number of Patients",
           title = "TUDOR-C Grey Zone Reclassification",
           subtitle = sprintf("NRI = %.3f | IDI = %.4f", nri_total, idi)) +
      theme_nature(base_size = 8) +
      theme(legend.position = "top")

    save_fig(p_reclass, "fig_tudor_c_reclassification", w_mm = 140, h_mm = 100)
  }
}, error = function(e) cat(sprintf("  Figure C error: %s\n", e$message)))

# ============================================================
# PART 10: SUMMARY
# ============================================================

cat("\n")
cat("================================================================\n")
cat("  TUDOR-C ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat("\n  Key Results:\n")
if (exists("mr_df")) {
  for (i in seq_len(nrow(mr_df))) {
    cat(sprintf("    %s: AUC=%.3f (%.3f-%.3f)\n",
                mr_df$Model[i], mr_df$AUC[i], mr_df$CI_lower[i], mr_df$CI_upper[i]))
  }
}
cat("\n  Output files saved to:\n")
cat(sprintf("    %s/tudor_c_*.csv\n", table_dir))
cat(sprintf("    %s/fig_tudor_c_*.pdf/png\n", fig_dir))
cat("================================================================\n")
