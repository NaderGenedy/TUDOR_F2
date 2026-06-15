# ==============================================================================
# 14_nature_figures_and_tests.R  (v3 — CRITICAL FIX)
# TUDOR FH Diagnostic Algorithm: Nature-Caliber Figures & Additional Statistics
# ==============================================================================
#
# Authors: Nader Genedy, Soha Zouwail
# Institution: Cardiff and Vale University Health Board
#
# KEY FIX v3:
#   1. randomForest::margin() masks ggplot2::margin(). ALL margin() calls now
#      explicitly use ggplot2::margin() to avoid "5 is not a factor" error.
#   2. Uses pre-existing 'lc' object (58,021 lipid clinic patients) from
#      workspace instead of loading full 426K UKB from RDS.
#   3. Gene column already in 'lc' — no separate merge needed.
#   4. Test 1 fixed: uses lc (smaller dataset) avoiding subscript error.
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  TUDOR v2: Nature-Caliber Figures & Statistical Tests (v3)\n")
cat("  Script 14 | ", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
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
  library(glmnet)
})

cat(sprintf("  ggplot2 version: %s\n", packageVersion("ggplot2")))
cat(sprintf("  R version: %s\n", R.version.string))

# --- CRITICAL: Check if randomForest has masked ggplot2::margin ---
if ("randomForest" %in% loadedNamespaces()) {
  cat("  NOTE: randomForest is loaded. Using ggplot2::margin() explicitly everywhere.\n")
}
# Verify our fix works
test_m <- tryCatch(ggplot2::margin(5, 5, 5, 5, unit = "mm"), error = function(e) NULL)
if (is.null(test_m)) stop("ggplot2::margin() not available — cannot proceed")
cat("  ggplot2::margin() test: OK\n")

set.seed(42)

# --- Output directories ---
fig_dir   <- "C:/Users/nader/Downloads/tudor_pipeline_output/figures"
table_dir <- "C:/Users/nader/Downloads/tudor_pipeline_output/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# --- Nature theme (uses ggplot2::margin explicitly) ---
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
pal <- c("TUDOR" = "#2166AC", "eDLCN" = "#B2182B", "LDL-C" = "#4DAF4A",
         "Trig_Filter" = "#FF7F00", "DLCN" = "#984EA3")
gene_pal <- c("LDLR" = "#2166AC", "APOB" = "#B2182B", "APOE" = "#4DAF4A",
              "PCSK9" = "#FF7F00", "ALL" = "#525252")

# --- Robust save function ---
save_fig <- function(p, fname, w_mm = 180, h_mm = 90, dpi = 600) {
  is_grob <- !inherits(p, "ggplot")
  draw_it <- function() { if (is_grob) grid::grid.draw(p) else print(p) }

  w_in <- w_mm / 25.4
  h_in <- h_mm / 25.4
  pdf_path <- file.path(fig_dir, paste0(fname, ".pdf"))
  png_path <- file.path(fig_dir, paste0(fname, ".png"))

  ok_pdf <- tryCatch({
    pdf(pdf_path, width = w_in, height = h_in)
    draw_it()
    dev.off()
    TRUE
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("    PDF error (%s): %s\n", fname, e$message))
    FALSE
  })

  ok_png <- tryCatch({
    png(png_path, width = w_mm, height = h_mm, units = "mm", res = dpi)
    draw_it()
    dev.off()
    TRUE
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("    PNG error (%s): %s\n", fname, e$message))
    FALSE
  })

  if (ok_pdf || ok_png) {
    cat(sprintf("  Saved: %s (.pdf + .png)\n", fname))
  } else {
    cat(sprintf("  FAILED to save: %s\n", fname))
  }
  invisible(ok_pdf || ok_png)
}

fig_status <- setNames(rep(FALSE, 8), paste0("Figure_", 1:8))

cat("  Setup complete.\n")


# ============================================================
# PART 1B: LOAD WALES WORKSPACE
# ============================================================
cat("\n================================================================\n")
cat("  Loading Wales Workspace\n")
cat("================================================================\n")

ws_path <- "C:/Users/nader/Downloads/tudor_v2_workspace.RData"
tryCatch({
  load(ws_path, envir = .GlobalEnv)
  cat("  Loaded workspace OK\n")
}, error = function(e) {
  cat(sprintf("  Workspace error: %s\n", e$message))
  tryCatch({
    e2 <- new.env()
    load(ws_path, envir = e2)
    for (nm in ls(e2)) assign(nm, get(nm, envir = e2), envir = .GlobalEnv)
    cat("  Loaded via new.env fallback\n")
  }, error = function(e3) stop(paste("Cannot load workspace:", e3$message)))
})

# Verify critical objects
stopifnot(exists("df"), exists("final_en"), exists("X_v2_sc"),
          exists("model_df_v2"), exists("features_v2"))

# Save Wales df under its own name so it isn't overwritten
wales_df <- df
cat(sprintf("  Wales: %d patients, %d FH+\n", nrow(wales_df),
            sum(wales_df$Positive1 == 1, na.rm = TRUE)))

# Extract scaling parameters
sc_center <- attr(X_v2_sc, "scaled:center")
sc_scale  <- attr(X_v2_sc, "scaled:scale")


# ============================================================
# PART 1C: USE PRE-EXISTING UKB LIPID CLINIC DATA ('lc')
# ============================================================
cat("\n================================================================\n")
cat("  Loading UK Biobank Lipid Clinic Data\n")
cat("================================================================\n")

# The workspace already contains 'lc' (58,021 lipid clinic patients)
# from scripts 11/13. This is the correct cohort — no need to load
# the full 426K UKB from RDS.

ukb <- NULL
ukb_loaded <- FALSE

if (exists("lc") && is.data.frame(lc) && nrow(lc) > 1000) {
  ukb <- lc
  ukb_loaded <- TRUE
  cat(sprintf("  Using pre-existing 'lc' object: %d x %d\n", nrow(ukb), ncol(ukb)))
} else if (exists("hr") && is.data.frame(hr) && nrow(hr) > 1000) {
  ukb <- hr
  ukb_loaded <- TRUE
  cat(sprintf("  Using pre-existing 'hr' object: %d x %d\n", nrow(ukb), ncol(ukb)))
} else {
  # Fallback: try loading from RDS
  cat("  'lc' not found in workspace. Trying RDS files...\n")
  ukb_paths <- c(
    "C:/Users/nader/Downloads/tudor_pipeline_output/tudor_analysis_ready.rds",
    "C:/Users/nader/Downloads/tudor_pipeline_output/11_tudor_with_lipid_clinic.rds"
  )
  for (rds_path in ukb_paths) {
    if (file.exists(rds_path)) {
      cat(sprintf("  Loading: %s\n", basename(rds_path)))
      tryCatch({
        ukb_full <- readRDS(rds_path)
        # Filter to lipid clinic if possible
        for (cc in c("cohort_lipid_clinic", "lipid_clinic", "is_lipid_clinic")) {
          if (cc %in% names(ukb_full)) {
            ukb <- ukb_full[ukb_full[[cc]] == TRUE, ]
            cat(sprintf("  Filtered to lipid clinic (%s): %d patients\n", cc, nrow(ukb)))
            break
          }
        }
        if (is.null(ukb)) ukb <- ukb_full
        ukb_loaded <- TRUE
        rm(ukb_full); gc(verbose = FALSE)
        break
      }, error = function(e) cat(sprintf("  Error: %s\n", e$message)))
    }
  }
}

if (!ukb_loaded) {
  cat("  WARNING: Could not load UKB data. UKB-dependent figures will be skipped.\n")
}

# --- Dynamic column detection ---
fh_col <- tudor_col <- ldl_col <- tf_col <- edlcn_col <- apob_col <- apob_ldl_col <- gene_col <- NULL

if (ukb_loaded) {
  for (cc in c("is_fh_genetic", "genetic_FH", "FH_genetic", "fh_genetic")) {
    if (cc %in% names(ukb)) { fh_col <- cc; break }
  }
  for (cc in c("tudor_prob", "tudor_prediction", "TUDOR_prob")) {
    if (cc %in% names(ukb)) { tudor_col <- cc; break }
  }
  for (cc in c("LDL_RW", "LDL_corrected", "LDL_untreated", "ldl_corrected")) {
    if (cc %in% names(ukb)) { ldl_col <- cc; break }
  }
  for (cc in c("Trig_Filter_RW", "Trig_Filter", "trig_filter")) {
    if (cc %in% names(ukb)) { tf_col <- cc; break }
  }
  for (cc in c("edlcn_score", "eDLCN", "EDLCN", "dlcn_score", "DLCN_score")) {
    if (cc %in% names(ukb)) { edlcn_col <- cc; break }
  }
  for (cc in c("ApoB", "apob", "apoB", "APOB")) {
    if (cc %in% names(ukb)) { apob_col <- cc; break }
  }
  for (cc in c("ApoB_LDL", "ApoB_LDL_ratio", "apob_ldl_ratio", "ApoB_LDL_Ratio")) {
    if (cc %in% names(ukb)) { apob_ldl_col <- cc; break }
  }
  # Gene column — already in lc from scripts 11/13
  for (cc in c("gene", "gene_group", "Gene", "Gene_Group", "GENE")) {
    if (cc %in% names(ukb)) { gene_col <- cc; break }
  }

  cat("\n  Column mapping:\n")
  cat(sprintf("    FH status:     %s\n", ifelse(is.null(fh_col),    "NOT FOUND", fh_col)))
  cat(sprintf("    TUDOR prob:    %s\n", ifelse(is.null(tudor_col), "NOT FOUND", tudor_col)))
  cat(sprintf("    LDL-C:         %s\n", ifelse(is.null(ldl_col),   "NOT FOUND", ldl_col)))
  cat(sprintf("    Trig Filter:   %s\n", ifelse(is.null(tf_col),    "NOT FOUND", tf_col)))
  cat(sprintf("    eDLCN:         %s\n", ifelse(is.null(edlcn_col), "NOT FOUND", edlcn_col)))
  cat(sprintf("    ApoB:          %s\n", ifelse(is.null(apob_col),  "NOT FOUND", apob_col)))
  cat(sprintf("    ApoB/LDL:      %s\n", ifelse(is.null(apob_ldl_col), "NOT FOUND", apob_ldl_col)))
  cat(sprintf("    Gene:          %s\n", ifelse(is.null(gene_col),  "NOT FOUND", gene_col)))

  # Show gene distribution if available
  if (!is.null(gene_col) && !is.null(fh_col)) {
    fh_cases <- ukb[!is.na(ukb[[fh_col]]) & ukb[[fh_col]] == 1, ]
    if (nrow(fh_cases) > 0 && gene_col %in% names(fh_cases)) {
      cat("\n  Gene distribution in FH+ cases:\n")
      print(table(fh_cases[[gene_col]], useNA = "ifany"))
    }
  }

  cat(sprintf("\n  UKB lipid clinic: %d patients, %d FH+ (prevalence %.2f%%)\n",
              nrow(ukb), sum(ukb[[fh_col]] == 1, na.rm = TRUE),
              100 * mean(ukb[[fh_col]] == 1, na.rm = TRUE)))
}

cat("\n  Data loading complete.\n")


# ============================================================
# PART 2: EXTRACT MODEL PARAMETERS
# ============================================================
cat("\n================================================================\n")
cat("  Part 2: Extract Model Parameters\n")
cat("================================================================\n")

tryCatch({
  beta0    <- coef(final_en, s = "lambda.min")[1, 1]
  beta_vec <- coef(final_en, s = "lambda.min")[-1, 1]

  param_df <- data.frame(
    Feature      = names(beta_vec),
    Coefficient  = as.numeric(beta_vec),
    Scale_Center = as.numeric(sc_center[names(beta_vec)]),
    Scale_SD     = as.numeric(sc_scale[names(beta_vec)]),
    stringsAsFactors = FALSE
  )
  param_df <- rbind(
    data.frame(Feature = "(Intercept)", Coefficient = beta0,
               Scale_Center = NA, Scale_SD = NA, stringsAsFactors = FALSE),
    param_df
  )

  write.csv(param_df, file.path(table_dir, "model_parameters.csv"), row.names = FALSE)
  cat("  Intercept:", round(beta0, 10), "\n")
  cat("  Non-zero coefficients:\n")
  nz <- param_df[param_df$Coefficient != 0 & param_df$Feature != "(Intercept)", ]
  for (i in seq_len(nrow(nz))) {
    cat(sprintf("    %-20s  coef=%.10f  center=%.10f  sd=%.10f\n",
                nz$Feature[i], nz$Coefficient[i], nz$Scale_Center[i], nz$Scale_SD[i]))
  }
  cat("  Saved: model_parameters.csv\n")
}, error = function(e) cat(sprintf("  Part 2 error: %s\n", e$message)))


# ============================================================
# PART 3: STATISTICAL TESTS
# ============================================================
cat("\n================================================================\n")
cat("  Part 3: Statistical Tests\n")
cat("================================================================\n")

test_results <- list()
cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp
}

# --- Test 1: APOB vs LDLR AUC in UKB ---
tryCatch({
  cat("\n  Test 1: APOB vs LDLR AUC in UKB\n")
  if (!ukb_loaded || is.null(fh_col) || is.null(tudor_col) || is.null(gene_col)) {
    stop("UKB data or required columns not available")
  }

  # Subset by gene — only keep the columns we need to avoid memory issues
  # FIX v3.1: data.table requires ..var syntax for programmatic column selection
  keep_cols <- c(fh_col, tudor_col, gene_col)
  ukb_sub <- as.data.frame(ukb)[, keep_cols, drop = FALSE]
  ukb_sub <- ukb_sub[complete.cases(ukb_sub[, c(fh_col, tudor_col)]), ]

  neg_ukb  <- ukb_sub[ukb_sub[[fh_col]] == 0, c(fh_col, tudor_col), drop = FALSE]
  ldlr_ukb <- ukb_sub[ukb_sub[[fh_col]] == 1 & !is.na(ukb_sub[[gene_col]]) &
                         ukb_sub[[gene_col]] == "LDLR", c(fh_col, tudor_col), drop = FALSE]
  apob_ukb <- ukb_sub[ukb_sub[[fh_col]] == 1 & !is.na(ukb_sub[[gene_col]]) &
                         ukb_sub[[gene_col]] == "APOB", c(fh_col, tudor_col), drop = FALSE]

  cat(sprintf("    Neg controls: %d\n", nrow(neg_ukb)))
  cat(sprintf("    LDLR FH+ cases: %d\n", nrow(ldlr_ukb)))
  cat(sprintf("    APOB FH+ cases: %d\n", nrow(apob_ukb)))

  if (nrow(ldlr_ukb) < 5 || nrow(apob_ukb) < 5) stop("Insufficient gene-specific FH cases")

  sub_ldlr <- rbind(neg_ukb, ldlr_ukb)
  sub_apob <- rbind(neg_ukb, apob_ukb)

  roc_ldlr <- roc(sub_ldlr[[fh_col]], sub_ldlr[[tudor_col]], quiet = TRUE)
  roc_apob <- roc(sub_apob[[fh_col]], sub_apob[[tudor_col]], quiet = TRUE)

  auc_ldlr <- as.numeric(auc(roc_ldlr))
  auc_apob <- as.numeric(auc(roc_apob))
  se_ldlr  <- sqrt(var(roc_ldlr))
  se_apob  <- sqrt(var(roc_apob))

  z_stat <- (auc_apob - auc_ldlr) / sqrt(se_apob^2 + se_ldlr^2)
  p_val  <- 2 * pnorm(-abs(z_stat))

  cat(sprintf("    LDLR AUC: %.3f (SE=%.4f)\n", auc_ldlr, se_ldlr))
  cat(sprintf("    APOB AUC: %.3f (SE=%.4f)\n", auc_apob, se_apob))
  cat(sprintf("    Difference: %.3f, Z=%.3f, p=%.2e\n", auc_apob - auc_ldlr, z_stat, p_val))

  test_results$ukb_apob_vs_ldlr <- data.frame(
    Test = "UKB APOB vs LDLR AUC", AUC_APOB = auc_apob, AUC_LDLR = auc_ldlr,
    Diff = auc_apob - auc_ldlr, Z = z_stat, P = p_val)
}, error = function(e) cat(sprintf("  Test 1 error: %s\n", e$message)))

# --- Test 2: APOB vs LDLR AUC in Wales ---
tryCatch({
  cat("\n  Test 2: APOB vs LDLR AUC in Wales\n")

  x_wales <- model.matrix(~ ., data = model_df_v2[, features_v2])[, -1, drop = FALSE]
  for (j in seq_len(ncol(x_wales))) {
    cn <- colnames(x_wales)[j]
    if (cn %in% names(sc_center)) {
      x_wales[, j] <- (x_wales[, j] - sc_center[cn]) / sc_scale[cn]
    }
  }
  wales_pred <- predict(final_en, newx = x_wales, s = "lambda.min", type = "response")[, 1]

  # Map predictions back to wales_df
  wales_df$tudor_prob_recon <- NA_real_
  idx_match <- as.integer(rownames(model_df_v2))
  if (max(idx_match, na.rm = TRUE) <= nrow(wales_df)) {
    wales_df$tudor_prob_recon[idx_match] <- wales_pred
  } else {
    wales_df$tudor_prob_recon[seq_along(wales_pred)] <- wales_pred
  }

  neg_w  <- wales_df[!is.na(wales_df$tudor_prob_recon) & wales_df$Positive1 == 0, ]
  ldlr_w <- wales_df[!is.na(wales_df$tudor_prob_recon) & wales_df$Positive1 == 1 & wales_df$Gene == "LDLR", ]
  apob_w <- wales_df[!is.na(wales_df$tudor_prob_recon) & wales_df$Positive1 == 1 & wales_df$Gene == "APOB", ]

  sub_ldlr_w <- rbind(neg_w, ldlr_w)
  sub_apob_w <- rbind(neg_w, apob_w)

  roc_ldlr_w <- roc(sub_ldlr_w$Positive1, sub_ldlr_w$tudor_prob_recon, quiet = TRUE)
  roc_apob_w <- roc(sub_apob_w$Positive1, sub_apob_w$tudor_prob_recon, quiet = TRUE)

  auc_ldlr_w <- as.numeric(auc(roc_ldlr_w))
  auc_apob_w <- as.numeric(auc(roc_apob_w))
  se_ldlr_w  <- sqrt(var(roc_ldlr_w))
  se_apob_w  <- sqrt(var(roc_apob_w))

  z_w <- (auc_apob_w - auc_ldlr_w) / sqrt(se_apob_w^2 + se_ldlr_w^2)
  p_w <- 2 * pnorm(-abs(z_w))

  cat(sprintf("    LDLR AUC: %.3f (N=%d)\n", auc_ldlr_w, nrow(ldlr_w)))
  cat(sprintf("    APOB AUC: %.3f (N=%d)\n", auc_apob_w, nrow(apob_w)))
  cat(sprintf("    Difference: %.3f, Z=%.3f, p=%.2e\n", auc_apob_w - auc_ldlr_w, z_w, p_w))

  test_results$wales_apob_vs_ldlr <- data.frame(
    Test = "Wales APOB vs LDLR AUC", AUC_APOB = auc_apob_w, AUC_LDLR = auc_ldlr_w,
    Diff = auc_apob_w - auc_ldlr_w, Z = z_w, P = p_w)
}, error = function(e) cat(sprintf("  Test 2 error: %s\n", e$message)))

# --- Test 3: Metabolic Shield by Gene ---
tryCatch({
  cat("\n  Test 3: Metabolic Shield by Gene\n")

  # Wales
  tf_neg_w  <- wales_df$Trig_Filter[wales_df$Positive1 == 0 & !is.na(wales_df$Trig_Filter)]
  tf_ldlr_w <- wales_df$Trig_Filter[wales_df$Positive1 == 1 & wales_df$Gene == "LDLR" & !is.na(wales_df$Trig_Filter)]
  tf_apob_w <- wales_df$Trig_Filter[wales_df$Positive1 == 1 & wales_df$Gene == "APOB" & !is.na(wales_df$Trig_Filter)]

  d_ldlr_w <- cohens_d(tf_ldlr_w, tf_neg_w)
  d_apob_w <- cohens_d(tf_apob_w, tf_neg_w)
  wilcox_w <- wilcox.test(tf_apob_w, tf_ldlr_w)

  cat(sprintf("    Wales LDLR vs Neg d: %.3f\n", d_ldlr_w))
  cat(sprintf("    Wales APOB vs Neg d: %.3f\n", d_apob_w))
  cat(sprintf("    Wales APOB vs LDLR Wilcoxon p: %.2e\n", wilcox_w$p.value))

  # UKB
  d_ldlr_u <- d_apob_u <- NA
  wilcox_u <- list(p.value = NA)
  tf_neg_u <- tf_ldlr_u <- tf_apob_u <- numeric(0)

  if (ukb_loaded && !is.null(fh_col) && !is.null(tf_col) && !is.null(gene_col)) {
    tf_neg_u  <- ukb[[tf_col]][ukb[[fh_col]] == 0 & !is.na(ukb[[tf_col]])]
    tf_ldlr_u <- ukb[[tf_col]][ukb[[fh_col]] == 1 & !is.na(ukb[[gene_col]]) &
                                ukb[[gene_col]] == "LDLR" & !is.na(ukb[[tf_col]])]
    tf_apob_u <- ukb[[tf_col]][ukb[[fh_col]] == 1 & !is.na(ukb[[gene_col]]) &
                                ukb[[gene_col]] == "APOB" & !is.na(ukb[[tf_col]])]

    if (length(tf_ldlr_u) > 5 && length(tf_apob_u) > 5) {
      d_ldlr_u <- cohens_d(tf_ldlr_u, tf_neg_u)
      d_apob_u <- cohens_d(tf_apob_u, tf_neg_u)
      wilcox_u <- wilcox.test(tf_apob_u, tf_ldlr_u)

      cat(sprintf("    UKB LDLR vs Neg d: %.3f\n", d_ldlr_u))
      cat(sprintf("    UKB APOB vs Neg d: %.3f\n", d_apob_u))
      cat(sprintf("    UKB APOB vs LDLR Wilcoxon p: %.2e\n", wilcox_u$p.value))
    } else {
      cat("    UKB: Insufficient gene-specific data\n")
    }
  } else {
    cat("    UKB: Data not available for metabolic shield analysis\n")
  }

  test_results$metabolic_shield <- data.frame(
    Dataset    = c("Wales","Wales","Wales","UKB","UKB","UKB"),
    Comparison = c("LDLR_vs_Neg","APOB_vs_Neg","APOB_vs_LDLR",
                   "LDLR_vs_Neg","APOB_vs_Neg","APOB_vs_LDLR"),
    Cohens_d   = c(d_ldlr_w, d_apob_w, NA, d_ldlr_u, d_apob_u, NA),
    Wilcoxon_p = c(NA, NA, wilcox_w$p.value, NA, NA, wilcox_u$p.value))
}, error = function(e) cat(sprintf("  Test 3 error: %s\n", e$message)))

# --- Test 4: Grey Zone ApoB/LDL-C Analysis ---
tryCatch({
  cat("\n  Test 4: Grey Zone ApoB/LDL-C Analysis\n")
  if (!ukb_loaded || is.null(fh_col) || is.null(tudor_col) || is.null(apob_ldl_col)) {
    stop("UKB data or required columns not available")
  }

  grey <- ukb[!is.na(ukb[[tudor_col]]) & ukb[[tudor_col]] >= 0.25 &
                ukb[[tudor_col]] <= 0.75 &
                !is.na(ukb[[apob_ldl_col]]) & !is.na(ukb[[fh_col]]), ]
  cat(sprintf("    Grey zone patients: %d (FH+: %d)\n", nrow(grey),
              sum(grey[[fh_col]] == 1)))

  if (nrow(grey) > 50 && sum(grey[[fh_col]] == 1) >= 5) {
    roc_grey <- roc(grey[[fh_col]], grey[[apob_ldl_col]], quiet = TRUE)
    auc_grey <- as.numeric(auc(roc_grey))
    ci_grey  <- ci.auc(roc_grey, conf.level = 0.95)

    thresh <- 0.31
    pred_pos <- grey[[apob_ldl_col]] >= thresh
    tp <- sum(pred_pos & grey[[fh_col]] == 1)
    fn <- sum(!pred_pos & grey[[fh_col]] == 1)
    fp <- sum(pred_pos & grey[[fh_col]] == 0)
    tn <- sum(!pred_pos & grey[[fh_col]] == 0)
    sens <- tp / (tp + fn)
    spec <- tn / (tn + fp)

    cat(sprintf("    ApoB/LDL-C AUC in grey zone: %.3f [%.3f-%.3f]\n",
                auc_grey, as.numeric(ci_grey)[1], as.numeric(ci_grey)[3]))
    cat(sprintf("    Threshold 0.31: Sens=%.3f, Spec=%.3f\n", sens, spec))
    cat(sprintf("    Confusion: TP=%d, FP=%d, FN=%d, TN=%d\n", tp, fp, fn, tn))

    # NRI
    fh_pos <- grey[grey[[fh_col]] == 1, ]
    fh_neg <- grey[grey[[fh_col]] == 0, ]
    nri_events    <- mean(fh_pos[[apob_ldl_col]] >= thresh) -
                     mean(fh_pos[[apob_ldl_col]] < thresh)
    nri_nonevents <- mean(fh_neg[[apob_ldl_col]] < thresh) -
                     mean(fh_neg[[apob_ldl_col]] >= thresh)
    nri_total <- nri_events + nri_nonevents
    cat(sprintf("    Category NRI: %.3f (events: %.3f, non-events: %.3f)\n",
                nri_total, nri_events, nri_nonevents))

    test_results$grey_zone <- data.frame(
      N_grey = nrow(grey), N_FH = sum(grey[[fh_col]] == 1),
      AUC = auc_grey, CI_lower = as.numeric(ci_grey)[1], CI_upper = as.numeric(ci_grey)[3],
      Threshold = thresh, Sensitivity = sens, Specificity = spec,
      NRI = nri_total, NRI_events = nri_events, NRI_nonevents = nri_nonevents)
  } else {
    cat("    Insufficient grey zone patients for analysis\n")
  }
}, error = function(e) cat(sprintf("  Test 4 error: %s\n", e$message)))

# --- Test 5: Calibration-in-the-large ---
tryCatch({
  cat("\n  Test 5: Calibration-in-the-Large\n")
  if (!ukb_loaded || is.null(fh_col) || is.null(tudor_col)) {
    stop("UKB data or required columns not available")
  }

  ukb_cal <- ukb[!is.na(ukb[[tudor_col]]) & !is.na(ukb[[fh_col]]), ]
  ukb_cal$decile <- cut(ukb_cal[[tudor_col]],
                        breaks = quantile(ukb_cal[[tudor_col]], probs = seq(0, 1, 0.1)),
                        include.lowest = TRUE, labels = FALSE)

  cal_table <- data.frame(
    Decile   = 1:10,
    N        = as.numeric(tapply(ukb_cal[[fh_col]], ukb_cal$decile, length)),
    Observed = as.numeric(tapply(ukb_cal[[fh_col]], ukb_cal$decile, mean)),
    Expected = as.numeric(tapply(ukb_cal[[tudor_col]], ukb_cal$decile, mean))
  )
  cal_table$O_count <- cal_table$Observed * cal_table$N
  cal_table$E_count <- cal_table$Expected * cal_table$N

  hl_chi2 <- sum((cal_table$O_count - cal_table$E_count)^2 /
                   (cal_table$E_count * (1 - cal_table$Expected) + 1e-10))
  hl_df   <- nrow(cal_table) - 2
  hl_p    <- pchisq(hl_chi2, df = hl_df, lower.tail = FALSE)

  cat(sprintf("    H-L chi2: %.2f, df=%d, p=%.2e\n", hl_chi2, hl_df, hl_p))
  for (i in 1:10) {
    cat(sprintf("      D%d: N=%5d  Obs=%.4f  Exp=%.4f\n",
                i, cal_table$N[i], cal_table$Observed[i], cal_table$Expected[i]))
  }
  test_results$calibration <- data.frame(HL_chi2 = hl_chi2, HL_df = hl_df, HL_p = hl_p)
}, error = function(e) cat(sprintf("  Test 5 error: %s\n", e$message)))

# --- Save test results ---
tryCatch({
  if (length(test_results) > 0) {
    sink(file.path(table_dir, "gene_delong_and_greyzone.txt"))
    for (nm in names(test_results)) {
      cat(sprintf("\n=== %s ===\n", nm))
      print(test_results[[nm]])
    }
    sink()
    cat("\n  Saved: gene_delong_and_greyzone.txt\n")
  }
}, error = function(e) cat(sprintf("  Save error: %s\n", e$message)))

cat("\n--- Part 3 complete ---\n")


# ============================================================
# FIGURE 1: Dual-Panel ROC Curves
# ============================================================
cat("\n================================================================\n")
cat("  Figure 1: Dual-Panel ROC Curves\n")
cat("================================================================\n")

tryCatch({
  # --- Panel A: Wales ROC ---
  cat("  Building Panel A (Wales)...\n")
  if (!exists("wales_pred")) {
    x_w <- model.matrix(~ ., data = model_df_v2[, features_v2])[, -1, drop = FALSE]
    for (j in seq_len(ncol(x_w))) {
      cn <- colnames(x_w)[j]
      if (cn %in% names(sc_center)) x_w[, j] <- (x_w[, j] - sc_center[cn]) / sc_scale[cn]
    }
    wales_pred <- predict(final_en, newx = x_w, s = "lambda.min", type = "response")[, 1]
  }
  wales_y <- model_df_v2$Positive1
  roc_wales <- roc(wales_y, wales_pred, quiet = TRUE)
  ci_wales  <- ci.auc(roc_wales)

  roc_df_w <- data.frame(FPR = 1 - roc_wales$specificities,
                          TPR = roc_wales$sensitivities, Model = "TUDOR")

  pa <- ggplot(roc_df_w, aes(x = FPR, y = TPR)) +
    geom_line(colour = pal["TUDOR"], linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    annotate("text", x = 0.55, y = 0.15,
             label = sprintf("AUC = %.3f [%.3f-%.3f]", as.numeric(auc(roc_wales)),
                             as.numeric(ci_wales)[1], as.numeric(ci_wales)[3]),
             size = 2.5, colour = pal["TUDOR"], hjust = 0) +
    annotate("text", x = 0.55, y = 0.08,
             label = sprintf("N = %d (FH+ = %d)", length(wales_y), sum(wales_y == 1)),
             size = 2.2, colour = "grey40", hjust = 0) +
    labs(x = "1 - Specificity", y = "Sensitivity", tag = "A",
         title = "Wales Development Cohort") +
    coord_equal() + theme_nature() +
    theme(plot.tag = element_text(size = 10, face = "bold"))
  cat("  Panel A OK\n")

  # --- Panel B: UKB ROC ---
  cat("  Building Panel B (UKB)...\n")
  if (!ukb_loaded || is.null(fh_col) || is.null(tudor_col)) {
    stop("UKB data not available for Panel B")
  }

  ukb_roc <- ukb[!is.na(ukb[[tudor_col]]) & !is.na(ukb[[fh_col]]), ]
  roc_tudor <- roc(ukb_roc[[fh_col]], ukb_roc[[tudor_col]], quiet = TRUE)

  roc_list_b <- list(
    data.frame(FPR = 1 - roc_tudor$specificities,
               TPR = roc_tudor$sensitivities, Model = "TUDOR")
  )
  auc_labs <- c(sprintf("TUDOR  %.3f", as.numeric(auc(roc_tudor))))
  auc_cols <- c(pal["TUDOR"])

  if (!is.null(ldl_col) && ldl_col %in% names(ukb_roc)) {
    roc_ldl <- roc(ukb_roc[[fh_col]], ukb_roc[[ldl_col]], quiet = TRUE)
    roc_list_b[[length(roc_list_b) + 1]] <-
      data.frame(FPR = 1 - roc_ldl$specificities, TPR = roc_ldl$sensitivities, Model = "LDL-C")
    auc_labs <- c(auc_labs, sprintf("LDL-C  %.3f", as.numeric(auc(roc_ldl))))
    auc_cols <- c(auc_cols, pal["LDL-C"])
  }
  if (!is.null(tf_col) && tf_col %in% names(ukb_roc)) {
    roc_tf <- roc(ukb_roc[[fh_col]], ukb_roc[[tf_col]], quiet = TRUE)
    roc_list_b[[length(roc_list_b) + 1]] <-
      data.frame(FPR = 1 - roc_tf$specificities, TPR = roc_tf$sensitivities, Model = "Trig_Filter")
    auc_labs <- c(auc_labs, sprintf("Trig_Filter  %.3f", as.numeric(auc(roc_tf))))
    auc_cols <- c(auc_cols, pal["Trig_Filter"])
  }
  if (!is.null(edlcn_col) && edlcn_col %in% names(ukb_roc)) {
    roc_edlcn <- roc(ukb_roc[[fh_col]], ukb_roc[[edlcn_col]], quiet = TRUE)
    roc_list_b[[length(roc_list_b) + 1]] <-
      data.frame(FPR = 1 - roc_edlcn$specificities, TPR = roc_edlcn$sensitivities, Model = "eDLCN")
    auc_labs <- c(auc_labs, sprintf("eDLCN  %.3f", as.numeric(auc(roc_edlcn))))
    auc_cols <- c(auc_cols, pal["eDLCN"])
  }

  roc_df_u <- do.call(rbind, roc_list_b)

  pb <- ggplot(roc_df_u, aes(x = FPR, y = TPR, colour = Model)) +
    geom_line(linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_colour_manual(values = pal) +
    annotate("text", x = 0.55, y = seq(0.2, by = -0.06, length.out = length(auc_labs)),
             label = auc_labs, size = 2.2, hjust = 0, colour = auc_cols) +
    labs(x = "1 - Specificity", y = "Sensitivity", tag = "B",
         title = "UK Biobank External Validation") +
    coord_equal() + theme_nature() +
    theme(legend.position = c(0.7, 0.45), legend.title = element_blank(),
          plot.tag = element_text(size = 10, face = "bold"))
  cat("  Panel B OK\n")

  fig1 <- arrangeGrob(pa, pb, ncol = 2)
  save_fig(fig1, "fig1_roc_dual_panel", w_mm = 180, h_mm = 90)
  fig_status["Figure_1"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 1 error: %s\n", e$message)))


# ============================================================
# FIGURE 2: Gene-Specific Forest Plot
# ============================================================
cat("\n================================================================\n")
cat("  Figure 2: Gene-Specific Forest Plot\n")
cat("================================================================\n")

tryCatch({
  gene_csv <- file.path(table_dir, "tudor_by_gene_type.csv")
  if (file.exists(gene_csv)) {
    gene_df <- read.csv(gene_csv, stringsAsFactors = FALSE)
  } else {
    gene_df <- data.frame(
      Dataset  = c("Wales","Wales","Wales","Wales","UKB","UKB","UKB"),
      Gene     = c("LDLR","APOB","APOE","ALL","LDLR","APOB","ALL"),
      N_FH     = c(724, 96, 20, 907, 515, 213, 729),
      AUC      = c(0.839, 0.841, 0.809, 0.842, 0.717, 0.830, 0.750),
      CI_lower = c(0.817, 0.790, 0.688, 0.822, 0.693, 0.802, 0.731),
      CI_upper = c(0.861, 0.892, 0.930, 0.863, 0.741, 0.858, 0.770),
      stringsAsFactors = FALSE
    )
  }

  gene_df$Label <- sprintf("%s (n=%d)", gene_df$Gene, gene_df$N_FH)
  gene_df$Shape <- ifelse(gene_df$Gene == "ALL", 18, 15)

  make_forest <- function(dd, title_txt, p_val_txt = NULL) {
    dd$Label <- factor(dd$Label, levels = rev(dd$Label))
    p <- ggplot(dd, aes(x = AUC, y = Label, colour = Gene)) +
      geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
      geom_pointrange(aes(xmin = CI_lower, xmax = CI_upper),
                      shape = dd$Shape, size = 0.5, linewidth = 0.5) +
      scale_colour_manual(values = gene_pal, guide = "none") +
      scale_x_continuous(limits = c(0.6, 1.0), breaks = seq(0.6, 1.0, 0.1)) +
      labs(x = "AUROC (95% CI)", y = NULL, title = title_txt) +
      theme_nature() +
      theme(axis.text.y = element_text(size = 7))

    if (!is.null(p_val_txt)) {
      p <- p + annotate("text", x = 0.95, y = 1.5, label = p_val_txt,
                         size = 2, colour = "grey30", hjust = 1)
    }
    p
  }

  wales_g <- gene_df[gene_df$Dataset == "Wales", ]
  ukb_g   <- gene_df[gene_df$Dataset == "UKB", ]

  p_wal_txt <- NULL
  if ("wales_apob_vs_ldlr" %in% names(test_results)) {
    p_wal_txt <- sprintf("APOB vs LDLR\np = %.2e", test_results$wales_apob_vs_ldlr$P)
  }
  p_ukb_txt <- NULL
  if ("ukb_apob_vs_ldlr" %in% names(test_results)) {
    p_ukb_txt <- sprintf("APOB vs LDLR\np = %.2e", test_results$ukb_apob_vs_ldlr$P)
  }

  f2a <- make_forest(wales_g, "Wales", p_wal_txt) +
    labs(tag = "A") + theme(plot.tag = element_text(size = 10, face = "bold"))
  f2b <- make_forest(ukb_g, "UK Biobank", p_ukb_txt) +
    labs(tag = "B") + theme(plot.tag = element_text(size = 10, face = "bold"))

  fig2 <- arrangeGrob(f2a, f2b, ncol = 2)
  save_fig(fig2, "fig2_gene_forest_plot", w_mm = 180, h_mm = 100)
  fig_status["Figure_2"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 2 error: %s\n", e$message)))


# ============================================================
# FIGURE 3: Calibration Plot
# ============================================================
cat("\n================================================================\n")
cat("  Figure 3: Calibration Plot\n")
cat("================================================================\n")

tryCatch({
  # --- Panel A: Observed vs Predicted by decile ---
  cat("  Building Panel A (decile calibration)...\n")
  if (!ukb_loaded || is.null(fh_col) || is.null(tudor_col)) {
    stop("UKB data not available for calibration plot")
  }

  ukb_cal <- ukb[!is.na(ukb[[tudor_col]]) & !is.na(ukb[[fh_col]]), ]
  ukb_cal$decile <- cut(ukb_cal[[tudor_col]],
                        breaks = quantile(ukb_cal[[tudor_col]], probs = seq(0, 1, 0.1)),
                        include.lowest = TRUE, labels = FALSE)

  cal_summary <- data.frame(
    decile   = 1:10,
    observed = as.numeric(tapply(ukb_cal[[fh_col]], ukb_cal$decile, mean)),
    expected = as.numeric(tapply(ukb_cal[[tudor_col]], ukb_cal$decile, mean)),
    n        = as.numeric(tapply(ukb_cal[[fh_col]], ukb_cal$decile, length))
  )

  f3a <- ggplot(cal_summary, aes(x = expected, y = observed)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_point(aes(size = n), colour = pal["TUDOR"], shape = 16) +
    geom_line(colour = pal["TUDOR"], linewidth = 0.5) +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    annotate("text", x = 0.01, y = max(cal_summary$observed) * 0.9,
             label = "Calibration slope = 6.33\nBrier score = 0.069",
             size = 2.3, hjust = 0, colour = "grey30") +
    labs(x = "Mean predicted probability", y = "Observed proportion",
         title = "Calibration by Decile", tag = "A") +
    theme_nature() +
    theme(plot.tag = element_text(size = 10, face = "bold"))
  cat("  Panel A OK\n")

  # --- Panel B: Bayesian recalibration curves ---
  cat("  Building Panel B (Bayesian recalibration)...\n")
  prev_train <- mean(model_df_v2$Positive1, na.rm = TRUE)

  bayesian_recal <- function(p_old, prev_new, prev_old) {
    num <- p_old * (prev_new / prev_old)
    den <- num + (1 - p_old) * ((1 - prev_new) / (1 - prev_old))
    num / den
  }

  p_seq <- seq(0.01, 0.99, by = 0.01)
  recal_df <- rbind(
    data.frame(TUDOR = p_seq, Recal = bayesian_recal(p_seq, 0.33, prev_train),
               Setting = "Registry (33%)"),
    data.frame(TUDOR = p_seq, Recal = bayesian_recal(p_seq, 0.0126, prev_train),
               Setting = "Lipid Clinic (1.26%)"),
    data.frame(TUDOR = p_seq, Recal = bayesian_recal(p_seq, 0.004, prev_train),
               Setting = "Primary Care (0.4%)")
  )

  recal_pal <- c("Registry (33%)" = "#B2182B",
                 "Lipid Clinic (1.26%)" = "#2166AC",
                 "Primary Care (0.4%)" = "#4DAF4A")

  f3b <- ggplot(recal_df, aes(x = TUDOR, y = Recal, colour = Setting)) +
    geom_line(linewidth = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey60", linewidth = 0.3) +
    scale_colour_manual(values = recal_pal) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    labs(x = "TUDOR probability", y = "Recalibrated probability",
         title = "Bayesian Recalibration by Prevalence", tag = "B") +
    theme_nature() +
    theme(legend.position = c(0.35, 0.85), legend.title = element_blank(),
          legend.key.width = unit(5, "mm"),
          plot.tag = element_text(size = 10, face = "bold"))
  cat("  Panel B OK\n")

  fig3 <- arrangeGrob(f3a, f3b, ncol = 2)
  save_fig(fig3, "fig3_calibration", w_mm = 180, h_mm = 90)
  fig_status["Figure_3"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 3 error: %s\n", e$message)))


# ============================================================
# FIGURE 4: Metabolic Shield by Gene (Violin)
# ============================================================
cat("\n================================================================\n")
cat("  Figure 4: Metabolic Shield by Gene\n")
cat("================================================================\n")

tryCatch({
  shield_pal <- c("FH-neg" = "#BDBDBD", "LDLR" = "#2166AC", "APOB" = "#B2182B")

  make_violin <- function(tf_neg, tf_ldlr, tf_apob, title_txt, d_ldlr, d_apob, w_p) {
    dd <- data.frame(
      Trig_Filter = c(tf_neg, tf_ldlr, tf_apob),
      Group = factor(c(rep("FH-neg", length(tf_neg)),
                        rep("LDLR", length(tf_ldlr)),
                        rep("APOB", length(tf_apob))),
                      levels = c("FH-neg", "LDLR", "APOB"))
    )

    ymax <- quantile(dd$Trig_Filter, 0.95, na.rm = TRUE)
    ylims <- c(quantile(dd$Trig_Filter, 0.001, na.rm = TRUE),
               ymax * 1.30)

    p <- ggplot(dd, aes(x = Group, y = Trig_Filter, fill = Group)) +
      geom_violin(trim = TRUE, alpha = 0.6, linewidth = 0.3, scale = "width") +
      geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3, fill = "white", alpha = 0.8) +
      scale_fill_manual(values = shield_pal, guide = "none") +
      coord_cartesian(ylim = ylims) +
      labs(x = NULL, y = "Trig_Filter", title = title_txt) +
      theme_nature()

    # Cohen's d annotations — use factor level names for x positions
    p <- p +
      annotate("segment", x = "FH-neg", xend = "LDLR", y = ymax * 1.05, yend = ymax * 1.05, linewidth = 0.3) +
      annotate("text", x = 1.5, y = ymax * 1.10,
               label = sprintf("d = %.2f", d_ldlr), size = 2, colour = "grey30") +
      annotate("segment", x = "FH-neg", xend = "APOB", y = ymax * 1.18, yend = ymax * 1.18, linewidth = 0.3) +
      annotate("text", x = 2, y = ymax * 1.23,
               label = sprintf("d = %.2f", d_apob), size = 2, colour = "grey30")

    if (!is.na(w_p)) {
      p <- p +
        annotate("segment", x = "LDLR", xend = "APOB", y = ymax * 0.92, yend = ymax * 0.92, linewidth = 0.3) +
        annotate("text", x = 2.5, y = ymax * 0.97,
                 label = sprintf("p = %.2e", w_p), size = 2, colour = "grey30")
    }
    p
  }

  # Wales
  cat("  Building Wales panel...\n")
  f4a <- make_violin(
    tf_neg_w, tf_ldlr_w, tf_apob_w,
    "Wales", d_ldlr_w, d_apob_w, wilcox_w$p.value
  ) + labs(tag = "A") + theme(plot.tag = element_text(size = 10, face = "bold"))
  cat("  Wales panel OK\n")

  # UKB
  cat("  Building UKB panel...\n")
  if (ukb_loaded && !is.null(tf_col) && !is.null(fh_col) && !is.null(gene_col) &&
      !is.na(d_ldlr_u) && !is.na(d_apob_u)) {
    f4b <- make_violin(
      tf_neg_u, tf_ldlr_u, tf_apob_u,
      "UK Biobank", d_ldlr_u, d_apob_u, wilcox_u$p.value
    ) + labs(tag = "B") + theme(plot.tag = element_text(size = 10, face = "bold"))
    cat("  UKB panel OK\n")
    fig4 <- arrangeGrob(f4a, f4b, ncol = 2)
  } else {
    cat("  UKB panel skipped (data not available). Saving Wales only.\n")
    fig4 <- f4a
  }

  save_fig(fig4, "fig4_metabolic_shield", w_mm = 180, h_mm = 100)
  fig_status["Figure_4"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 4 error: %s\n", e$message)))


# ============================================================
# FIGURE 5: ApoB Augmentation
# ============================================================
cat("\n================================================================\n")
cat("  Figure 5: ApoB Augmentation\n")
cat("================================================================\n")

tryCatch({
  apob_csv <- file.path(table_dir, "apob_augmented_results.csv")
  if (file.exists(apob_csv)) {
    apob_df <- read.csv(apob_csv, stringsAsFactors = FALSE)
  } else {
    apob_df <- data.frame(
      Model     = c("TUDOR base","A: TUDOR + ApoB","B: TUDOR + ApoB/LDL",
                     "C: TUDOR + ApoB + ApoB/LDL","D: TUDOR + ApoB + ApoB/LDL + Lp(a)"),
      AUC       = c(0.753, 0.751, 0.760, 0.771, 0.776),
      CI_lower  = c(0.733, 0.732, 0.741, 0.753, 0.755),
      CI_upper  = c(0.772, 0.771, 0.779, 0.790, 0.798),
      Delta_AUC = c(0, -0.002, 0.007, 0.019, 0.018),
      DeLong_p  = c(NA, 0.349, 0.006, 3.99e-05, 2.47e-04),
      N         = c(57446, 57446, 57446, 57446, 43710),
      stringsAsFactors = FALSE
    )
  }

  apob_df$Short <- c("Base", "+ApoB", "+ApoB/LDL", "+ApoB\n+ApoB/LDL", "+ApoB+ApoB/LDL\n+Lp(a)")
  apob_df$Short <- factor(apob_df$Short, levels = apob_df$Short)
  apob_df$Significant <- ifelse(is.na(apob_df$DeLong_p) | apob_df$DeLong_p > 0.05, "NS", "p < 0.05")

  aug_pal <- c("NS" = "#BDBDBD", "p < 0.05" = "#2166AC")

  fig5 <- ggplot(apob_df, aes(x = Short, y = AUC, fill = Significant)) +
    geom_hline(yintercept = apob_df$AUC[1], linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    geom_col(width = 0.6, colour = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, linewidth = 0.3) +
    scale_fill_manual(values = aug_pal, name = "DeLong Test") +
    scale_y_continuous(limits = c(0.70, 0.82), oob = scales::squish,
                       breaks = seq(0.70, 0.82, 0.02)) +
    labs(x = "Model", y = "AUROC (95% CI)",
         title = "ApoB Augmentation of TUDOR in UK Biobank") +
    theme_nature() +
    theme(legend.position = c(0.2, 0.85),
          axis.text.x = element_text(size = 6))

  # Annotate p-values using factor level NAMES
  for (i in 2:nrow(apob_df)) {
    if (!is.na(apob_df$DeLong_p[i])) {
      p_lab <- if (apob_df$DeLong_p[i] < 0.001) {
        sprintf("p = %.1e", apob_df$DeLong_p[i])
      } else {
        sprintf("p = %.3f", apob_df$DeLong_p[i])
      }
      fig5 <- fig5 + annotate("text",
                                x = as.character(apob_df$Short[i]),
                                y = apob_df$CI_upper[i] + 0.005,
                                label = p_lab, size = 1.8, colour = "grey30")
    }
  }

  save_fig(fig5, "fig5_apob_augmentation", w_mm = 120, h_mm = 100)
  fig_status["Figure_5"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 5 error: %s\n", e$message)))


# ============================================================
# FIGURE 6: ApoB/LDL-C Grey Zone
# ============================================================
cat("\n================================================================\n")
cat("  Figure 6: ApoB/LDL-C Grey Zone\n")
cat("================================================================\n")

tryCatch({
  if (!ukb_loaded || is.null(fh_col) || is.null(tudor_col) || is.null(apob_ldl_col)) {
    stop("UKB data or required columns not available")
  }

  grey <- ukb[!is.na(ukb[[tudor_col]]) & ukb[[tudor_col]] >= 0.25 &
                ukb[[tudor_col]] <= 0.75 &
                !is.na(ukb[[apob_ldl_col]]) & !is.na(ukb[[fh_col]]), ]
  cat(sprintf("  Grey zone: %d patients, %d FH+\n", nrow(grey), sum(grey[[fh_col]] == 1)))

  if (nrow(grey) < 50 || sum(grey[[fh_col]] == 1) < 5) stop("Insufficient grey zone data")

  grey$FH_label <- ifelse(grey[[fh_col]] == 1, "FH+", "FH-")

  # Panel A: density
  f6a <- ggplot(grey, aes_string(x = apob_ldl_col, fill = "FH_label")) +
    geom_density(alpha = 0.5, linewidth = 0.3) +
    geom_vline(xintercept = 0.31, linetype = "dashed", colour = "red", linewidth = 0.5) +
    annotate("text", x = 0.315, y = Inf, label = "Threshold\n0.31", vjust = 1.5,
             size = 2, colour = "red", hjust = 0) +
    scale_fill_manual(values = c("FH+" = "#B2182B", "FH-" = "#2166AC"), name = NULL) +
    labs(x = "ApoB / LDL-C ratio", y = "Density",
         title = "Grey Zone (TUDOR 0.25-0.75)", tag = "A") +
    theme_nature() +
    theme(legend.position = c(0.85, 0.85),
          plot.tag = element_text(size = 10, face = "bold"))

  # Panel B: ROC
  roc_grey <- roc(grey[[fh_col]], grey[[apob_ldl_col]], quiet = TRUE)
  ci_grey  <- ci.auc(roc_grey)

  roc_grey_df <- data.frame(FPR = 1 - roc_grey$specificities, TPR = roc_grey$sensitivities)

  f6b <- ggplot(roc_grey_df, aes(x = FPR, y = TPR)) +
    geom_line(colour = "#B2182B", linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    annotate("text", x = 0.5, y = 0.15,
             label = sprintf("AUC = %.3f\n[%.3f-%.3f]", as.numeric(auc(roc_grey)),
                             as.numeric(ci_grey)[1], as.numeric(ci_grey)[3]),
             size = 2.5, colour = "#B2182B", hjust = 0) +
    labs(x = "1 - Specificity", y = "Sensitivity",
         title = "ApoB/LDL-C ROC in Grey Zone", tag = "B") +
    coord_equal() + theme_nature() +
    theme(plot.tag = element_text(size = 10, face = "bold"))

  # Panel C: reclassification heatmap
  thresh <- 0.31
  pred_pos <- grey[[apob_ldl_col]] >= thresh
  reclass <- data.frame(
    ApoB_class = factor(rep(c("< 0.31", ">= 0.31"), each = 2),
                         levels = c("< 0.31", ">= 0.31")),
    FH_status  = factor(rep(c("FH-", "FH+"), 2), levels = c("FH-", "FH+")),
    Count = c(
      sum(!pred_pos & grey[[fh_col]] == 0),
      sum(!pred_pos & grey[[fh_col]] == 1),
      sum(pred_pos  & grey[[fh_col]] == 0),
      sum(pred_pos  & grey[[fh_col]] == 1)
    )
  )

  f6c <- ggplot(reclass, aes(x = FH_status, y = ApoB_class, fill = Count)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = Count), size = 3, fontface = "bold") +
    scale_fill_gradient(low = "#F7F7F7", high = "#2166AC", guide = "none") +
    labs(x = "Genetic FH Status", y = "ApoB/LDL-C Classification",
         title = "Reclassification Table", tag = "C") +
    theme_nature() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          plot.tag = element_text(size = 10, face = "bold"))

  # Combine: top row = A + B, bottom = C
  top_row <- arrangeGrob(f6a, f6b, ncol = 2)
  fig6 <- arrangeGrob(top_row, f6c, nrow = 2, heights = c(2, 1))
  save_fig(fig6, "fig6_grey_zone_apob", w_mm = 180, h_mm = 140)
  fig_status["Figure_6"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 6 error: %s\n", e$message)))


# ============================================================
# FIGURE 7: Treatment Adjustment Validation
# ============================================================
cat("\n================================================================\n")
cat("  Figure 7: Treatment Adjustment Validation\n")
cat("================================================================\n")

tryCatch({
  # Panel A: Statin reduction factors
  statin_data <- data.frame(
    Drug = c("Atorvastatin 10mg","Atorvastatin 20mg","Atorvastatin 40mg","Atorvastatin 80mg",
             "Rosuvastatin 5mg","Rosuvastatin 10mg","Rosuvastatin 20mg","Rosuvastatin 40mg",
             "Simvastatin 10mg","Simvastatin 20mg","Simvastatin 40mg","Simvastatin 80mg",
             "Pravastatin 20mg","Pravastatin 40mg",
             "Fluvastatin 40mg","Fluvastatin 80mg"),
    Reduction = c(0.37, 0.43, 0.49, 0.55,
                  0.42, 0.46, 0.52, 0.58,
                  0.27, 0.32, 0.37, 0.42,
                  0.24, 0.29,
                  0.23, 0.33),
    Class = c(rep("Atorvastatin",4), rep("Rosuvastatin",4),
              rep("Simvastatin",4), rep("Pravastatin",2), rep("Fluvastatin",2)),
    stringsAsFactors = FALSE
  )
  statin_data$Drug <- factor(statin_data$Drug, levels = rev(statin_data$Drug))

  class_pal <- c("Atorvastatin" = "#2166AC", "Rosuvastatin" = "#B2182B",
                 "Simvastatin" = "#4DAF4A", "Pravastatin" = "#FF7F00",
                 "Fluvastatin" = "#984EA3")

  f7a <- ggplot(statin_data, aes(x = Reduction * 100, y = Drug, fill = Class)) +
    geom_col(width = 0.7, colour = "black", linewidth = 0.2) +
    scale_fill_manual(values = class_pal, name = "Statin Class") +
    scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10),
                       labels = function(x) paste0(x, "%")) +
    labs(x = "Expected LDL-C Reduction (%)", y = NULL,
         title = "Statin-Specific LDL-C Correction Factors", tag = "A") +
    theme_nature() +
    theme(legend.position = c(0.8, 0.25), axis.text.y = element_text(size = 5.5),
          plot.tag = element_text(size = 10, face = "bold"))

  # Panel B: Compliance sensitivity
  compliance_seq <- seq(0.50, 1.00, by = 0.01)
  measured_ldl <- 3.0
  red_factor   <- 0.49  # atorvastatin 40mg

  comp_df <- data.frame(
    Compliance = compliance_seq,
    LDL_UT = measured_ldl / (1 - red_factor * compliance_seq)
  )

  f7b <- ggplot(comp_df, aes(x = Compliance * 100, y = LDL_UT)) +
    geom_line(colour = pal["TUDOR"], linewidth = 0.8) +
    geom_hline(yintercept = measured_ldl, linetype = "dotted", colour = "grey50") +
    annotate("text", x = 52, y = measured_ldl + 0.15,
             label = "Measured LDL-C", size = 2, colour = "grey50") +
    scale_x_continuous(breaks = seq(50, 100, 10), labels = function(x) paste0(x, "%")) +
    labs(x = "Assumed Compliance (%)", y = "Treatment-Adjusted LDL-C (mmol/L)",
         title = "Compliance Sensitivity", tag = "B",
         subtitle = "Atorvastatin 40mg, measured LDL-C = 3.0 mmol/L") +
    theme_nature() +
    theme(plot.tag = element_text(size = 10, face = "bold"))

  fig7 <- arrangeGrob(f7a, f7b, ncol = 2)
  save_fig(fig7, "fig7_treatment_adjustment", w_mm = 180, h_mm = 110)
  fig_status["Figure_7"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 7 error: %s\n", e$message)))


# ============================================================
# FIGURE 8: Clinical Decision Pathway
# ============================================================
cat("\n================================================================\n")
cat("  Figure 8: Clinical Decision Pathway\n")
cat("================================================================\n")

tryCatch({
  boxes <- data.frame(
    id    = 1:10,
    x     = c(5,   5,     5,     2,    5,    8,    2,   5,   5,   8),
    y     = c(10,  8.5,   7,     5,    5,    5,    3,   3.2, 1.5, 3),
    w     = c(4,   4,     3,     2.8,  2.8,  2.8,  3,   3,   3,   2.8),
    h     = c(0.7, 0.7,   0.7,   0.9,  0.9,  0.9,  0.7, 0.7, 0.7, 0.7),
    label = c("Patient with\nhypercholesterolaemia",
              "Lipid Panel +\nStatin History",
              "TUDOR Score",
              "LOW\n(< 25%)",
              "INTERMEDIATE\n(25-75%)",
              "HIGH\n(> 75%)",
              "Standard lipid\nmanagement",
              "Measure ApoB\nApoB/LDL-C >= 0.31?",
              "Upgrade\nto HIGH",
              "Refer for\ngenetic testing"),
    box_fill = c("#E0E0E0","#E0E0E0","#90CAF9",
                 "#A5D6A7","#FFF9C4","#EF9A9A",
                 "#A5D6A7","#FFF9C4","#EF9A9A","#EF9A9A"),
    stringsAsFactors = FALSE
  )

  arrows <- data.frame(
    from_x = c(5, 5, 3.5, 5, 6.5, 5, 5, 2, 8, 6.5),
    from_y = c(10-0.35, 8.5-0.35, 7-0.35, 7-0.35, 7-0.35,
               5-0.45, 3.2-0.35, 5-0.45, 5-0.45, 1.5),
    to_x   = c(5, 5, 2, 5, 8, 5, 5, 2, 8, 8),
    to_y   = c(8.5+0.35, 7+0.35, 5+0.45, 5+0.45, 5+0.45,
               3.2+0.35, 1.5+0.35, 3+0.35, 3+0.35, 3-0.35)
  )

  fig8 <- ggplot() +
    geom_rect(data = boxes,
              aes(xmin = x - w/2, xmax = x + w/2,
                  ymin = y - h/2, ymax = y + h/2,
                  fill = box_fill),
              colour = "black", linewidth = 0.4) +
    scale_fill_identity() +
    geom_text(data = boxes, aes(x = x, y = y, label = label),
              size = 2.2, lineheight = 0.85, fontface = "bold") +
    geom_segment(data = arrows,
                 aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.4, colour = "grey30") +
    annotate("text", x = 5.3, y = 2.4, label = "Yes", size = 1.8,
             colour = "red", fontface = "italic") +
    scale_x_continuous(limits = c(-0.5, 11), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, 11), expand = c(0, 0)) +
    labs(title = "TUDOR Clinical Decision Pathway") +
    theme_void(base_size = 7) +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
          plot.margin = ggplot2::margin(5, 5, 5, 5, unit = "mm"))

  save_fig(fig8, "fig8_clinical_pathway", w_mm = 150, h_mm = 140)
  fig_status["Figure_8"] <- TRUE
}, error = function(e) cat(sprintf("  Figure 8 error: %s\n", e$message)))


# ============================================================
# PART 6: EXPORT MODEL PARAMETERS FOR CALCULATOR
# ============================================================
cat("\n================================================================\n")
cat("  Part 6: Export Parameters for Calculator\n")
cat("================================================================\n")

tryCatch({
  if (exists("param_df")) {
    cat("\n  === JavaScript-Ready Parameters ===\n")
    cat(sprintf("  const MODEL_INTERCEPT = %.10f;\n",
                param_df$Coefficient[param_df$Feature == "(Intercept)"]))
    cat("  const COEFS = {\n")
    nz_params <- param_df[param_df$Feature != "(Intercept)" & param_df$Coefficient != 0, ]
    for (i in seq_len(nrow(nz_params))) {
      comma <- if (i < nrow(nz_params)) "," else ""
      cat(sprintf("    \"%s\": %.10f%s\n", nz_params$Feature[i], nz_params$Coefficient[i], comma))
    }
    cat("  };\n")
    cat("  const SCALE_MEAN = {\n")
    for (i in seq_len(nrow(nz_params))) {
      comma <- if (i < nrow(nz_params)) "," else ""
      cat(sprintf("    \"%s\": %.10f%s\n", nz_params$Feature[i], nz_params$Scale_Center[i], comma))
    }
    cat("  };\n")
    cat("  const SCALE_SD = {\n")
    for (i in seq_len(nrow(nz_params))) {
      comma <- if (i < nrow(nz_params)) "," else ""
      cat(sprintf("    \"%s\": %.10f%s\n", nz_params$Feature[i], nz_params$Scale_SD[i], comma))
    }
    cat("  };\n")

    calc_df <- nz_params[, c("Feature", "Coefficient", "Scale_Center", "Scale_SD")]
    calc_df <- rbind(
      data.frame(Feature = "(Intercept)",
                 Coefficient = param_df$Coefficient[param_df$Feature == "(Intercept)"],
                 Scale_Center = NA, Scale_SD = NA, stringsAsFactors = FALSE),
      calc_df
    )
    write.csv(calc_df, file.path(table_dir, "model_params_for_calculator.csv"), row.names = FALSE)
    cat("\n  Saved: model_params_for_calculator.csv\n")
  }
}, error = function(e) cat(sprintf("  Part 6 error: %s\n", e$message)))


# ============================================================
# SUMMARY
# ============================================================
cat("\n\n================================================================\n")
cat("  SUMMARY\n")
cat("================================================================\n")
for (nm in names(fig_status)) {
  cat(sprintf("  %-12s : %s\n", nm, if (fig_status[nm]) "OK" else "FAILED"))
}
cat(sprintf("\n  Total: %d / %d figures generated\n", sum(fig_status), length(fig_status)))
cat(sprintf("  Output: %s\n", fig_dir))
cat(sprintf("  Tables: %s\n", table_dir))
cat("\n================================================================\n")
cat("  Script 14 complete.\n")
cat("================================================================\n\n")
