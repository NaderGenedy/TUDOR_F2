# TUDOR Model Specification — Locked

**Source of truth:** `TUDOR_coefficients_locked.csv`
**Reproducer:** 66/68 PASS (2 failures to resolve — see O-7)

## Model Form
Elastic-net penalised logistic regression

## Hyperparameters (locked)
- Alpha: 0.5
- C: 1.0
- Solver: saga
- Seed: 20260518
- Intercept: -1.851

## Prediction Target
`P(LDLR pathogenic/likely-pathogenic carrier) = 1 / (1 + exp(-logit))`

`logit = -1.851 + Σ βⱼ·zⱼ` where zⱼ = standardised predictor

## 11 Locked Coefficients (per 1-SD, standardised)

| # | Variable | β | Conceptual Idea | Clinical Meaning |
|---|----------|---|-----------------|------------------|
| 1 | `ldl_ul` | **+1.609** | Lipid Age | Back-calculated untreated LDL (strongest predictor) |
| 2 | `on_statin` | +0.357 | Lipid Age | On lipid-lowering therapy at index (binary) |
| 3 | `tc_chem` | -0.743 | Triglyceride Shield | Total cholesterol (direct assay) |
| 4 | `tg_chem` | -0.604 | Triglyceride Shield | Triglycerides (direct assay) |
| 5 | `non_hdl` | -0.151 | Triglyceride Shield | Non-HDL cholesterol |
| 6 | `hdl_chem` | -0.137 | Triglyceride Shield | HDL cholesterol |
| 7 | `age_per_decade` | -0.357 | Demographic | Age, per decade |
| 8 | `sex_F` | -0.070 | Demographic | Female sex (binary) |
| 9 | `premature_mace` | +0.064 | Proband Effect | Premature MACE (age-/sex-defined) |
| 10 | `t2dm_int` | +0.047 | Proband Effect | Type 2 diabetes indicator/interaction |
| 11 | `fam_hist_cvd` | 0.000 | Proband Effect | Family history of CVD (L1-eliminated) |

## Three Ideas — The Conceptual Architecture

### 1. Lipid Age
"What would this patient's LDL be without treatment?"
- Back-calculates untreated LDL from drug × dose × duration × adherence
- Statin intensity bands: atorvastatin 25-48%, rosuvastatin 35-55%, simvastatin 20-42%, pravastatin 15-29%, fluvastatin 15-22%
- Add-ons: ezetimibe +20%, bempedoic acid +25%, PCSK9i +65% (capped 85% total)

### 2. Triglyceride Shield
"Is this isolated LDL elevation, or part of a mixed dyslipidaemia?"
- `Trig_Filter = ldl_ul / (tg_chem + 0.1)`
- Monogenic FH = isolated LDL ↑ with normal TG (high Trig_Filter)
- Metabolic dyslipidaemia = LDL ↑ AND TG ↑ (low Trig_Filter)
- Biology: Brown-Goldstein selective LDLR lesion (Nobel 1985) — breaks LDL clearance, leaves TG pathway intact

### 3. Proband Effect
"Was this patient found by clinical ascertainment or cascade screening?"
- `Index_Effect = (1 - Is_Relative) × ldl_ul`
- Probands present dramatically; cascade relatives carry the variant but attenuated phenotype
- Model adjusts for how the patient entered the system

## Comparators (applied as published rules, NO refit)
- DLCN (Dutch Lipid Clinic Network)
- Simon Broome
- MEDPED
- FAMCAT: not computable (coefficient supplement unavailable)

## Development
- Design A primary: n=3,099 (554 LDLR carriers)
- Cardiff and Vale catchment, South Wales
- Family-level deduplication; zero shared identifiers with validation
- MICE imputation; z-score standardisation
