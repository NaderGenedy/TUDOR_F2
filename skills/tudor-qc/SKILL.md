---
name: tudor-qc
description: Use whenever the TUDOR diagnostic-algorithm Python pipeline needs an end-to-end audit from raw UK Biobank extraction through the final Elastic Net model, OR when evaluating whether candidate features (T2DM, age-per-decade, LDLR genotype-tier, ApoB-continuous, or any user-specified additions) would improve TUDOR's discrimination, calibration, or net reclassification. Dispatches five independent verification agents (Raw Data Integrity, Statistical Methodology, Clinical Plausibility, Molecular / Genetic, Feature-Augmentation Experiment) in parallel, then produces a unified audit report plus an augmentation-experiment report. Trigger explicitly on phrases like "QC the TUDOR code from raw data", "any bugs in TUDOR", "audit the full TUDOR pipeline", "does adding T2DM improve TUDOR", "should we augment TUDOR with [feature]", "5 agents on TUDOR", or before any TUDOR submission / response-to-reviewer / R2 revision. Designed to catch the class of bug analogous to the prior WRONG_FIELDS endpoint-mapping issue (where UKB date fields p131286–p131294 silently coded I10–I20 hypertension rather than ASCVD), as well as treatment-adjustment, family-leakage, calibration-misweighting, and variant-calling errors specific to the FH diagnostic context.
---

# TUDOR-QC — five-agent end-to-end audit + feature-augmentation experiment

This skill is the submission-gate for the TUDOR diagnostic algorithm. It exists because TUDOR sits at the intersection of UK Biobank raw-field engineering, statin-treatment reconstruction, Elastic Net classification, and FH-specific molecular genetics — every one of those layers has had a substantive bug somewhere in the cardiology literature in the last five years, and any single bug invalidates the published AUCs. Five independent verification agents work in parallel because no single agent can hold all four expert frames (raw-data, statistics, clinical, molecular) at once with sufficient depth.

## When to invoke

Trigger on any of:

- "QC the TUDOR code" / "audit TUDOR from the start" / "any bugs in TUDOR?"
- "5 agents on TUDOR" / "parallel audit on the TUDOR pipeline"
- "Run the TUDOR augmentation experiment" / "does adding [T2DM / age / LDLR-tier / ApoB] improve TUDOR?"
- Before submitting any TUDOR revision, R2, or response-to-reviewer letter
- Whenever the underlying raw UKB extract is regenerated, refreshed, or re-mapped (e.g., after a RAP pull) — the skill catches drift between extraction logic and downstream model

Pair with `cox-analysis` if any survival-outcome extension is being added; pair with `manuscript-qc` once the audit passes and the manuscript text needs final claim-by-claim verification.

## Five-agent architecture — overview

Each agent runs **independently** with no shared state, in its own context window, in parallel. They report back to a central reconciliation step that assembles the unified report.

| Agent | Primary frame | Verifies | Output |
|---|---|---|---|
| 1. Raw Data Integrity | data-analysis + ukb-inventory | Raw UKB extraction logic; field-ID mappings; cohort filters; data-leakage barriers | `agent1_raw_data_audit.md` |
| 2. Statistical Methodology | biostatistics + reproducibility-tests | Elastic Net spec; cross-validation; AUC / calibration / NRI / Brier / DCA computation; TRIPOD compliance | `agent2_stats_methodology_audit.md` |
| 3. Clinical Plausibility | cardiology + metabolic-medicine + nature-reviewer | Treatment-adjusted LDL clinical sensibility; Trig_Filter direction; gene-specific pattern (APOB > PCSK9 > LDLR); cardiovascular-mortality gradient monotonicity | `agent3_clinical_plausibility_audit.md` |
| 4. Molecular / Genetic | molecular-genomics + genomics | Variant calling pipeline; LDLR severity tiers; APOB R3500/R3527 specificity; PCSK9 gain-of-function; family-level structure (Is_Relative) | `agent4_molecular_genetic_audit.md` |
| 5. Feature-Augmentation Experiment | biostatistics + data-analysis | Actually runs the augmentation experiment using the bundled script; produces ΔAUC, ΔCalibration, ΔNRI, ΔDCA for each candidate feature increment | `agent5_augmentation_report.md` |

A sixth coordination step (NOT an independent agent) reads all five outputs and produces `TUDOR_QC_REPORT.md` — the unified deliverable for the user.

## Data-leakage prevention — the load-bearing guarantee

Every run of this skill begins with a **mandatory leakage sweep** that runs before any model fitting. The sweep catches five classes of leakage that can silently inflate AUC, NRI, and DCA in a way that survives every other QC step:

| Class | What it catches | Threshold | Action |
|---|---|---|---|
| **L1. Outcome-feature circularity** | a candidate predictor is near-equivalent to the outcome label (e.g., adding LDLR-tier when FH-positive is defined by LDLR variant calling) | phi ≥ 0.95 (binary) or solo-AUC ≥ 0.98 (continuous) | **HALT** that augmentation variant (and any composite containing it) |
| **L2. Patient-level train/test overlap** | the same `eid` appears in both training and external-validation cohorts | any overlap | **HALT** the entire run; force cohort-construction fix |
| **L3. Family-level train/test overlap** | different `eid` but same `FamilyNumber` in both splits — the cascade-screening pattern that contaminates FH validation | any overlap | **HALT** the entire run |
| **L4. Component collinearity** | augmentation features are near-linear-combinations of TUDOR baseline (e.g., ApoB strongly co-determined by treatment-adjusted LDL) | VIF > 5 | **WARN** only — collinearity is a stability concern, not a leakage one |
| **L5. Cascade-relative within-cohort** | TUDOR's `Index_Effect = (1 − Is_Relative) × LDL_Untreated` design must be present in the model spec | `is_relative` column missing | **WARN** — relies on cohort-construction documentation |

The leakage sweep produces `leakage_sweep.csv` (one row per check, with severity / metric / threshold / detail) as the FIRST artefact of any TUDOR-QC run. Any HALT in L1 propagates the refusal to every composite variant that depends on the offending feature — so refusing LDLR-tier automatically refuses "TUDOR + all four" if it includes LDLR-tier dummies.

**Why this matters for TUDOR specifically.** FH-positive in the TUDOR cohort is *defined* by a pathogenic variant in one of six genes (LDLR / APOB / PCSK9 / APOE / LDLRAP1 / LDLR CNVs). Adding LDLR-tier as a predictor is therefore *circular by construction* — the model would be asked to predict what the variant call has already told it. The L1 gate refuses this augmentation automatically and reports it as **REFUSED-LEAKAGE** in the augmentation table rather than producing a misleadingly large ΔAUC. The LDLR-tier signal is the right *severity stratifier within FH-positives*, not a discrimination feature.

The smoke-test invariants this skill maintains (regression-tested on synthetic + injected-leakage data):
1. A synthetic cohort with perfect LDLR-tier ↔ FH overlap **refuses** both `TUDOR + LDLR-tier` AND `TUDOR + all four`.
2. An injected 50-`eid` train/test overlap **halts** the entire run with exit code 2.
3. Legitimate augmentations (T2DM, age, ApoB) continue to be evaluated and report ADOPT / DO-NOT-ADOPT per the four-criterion gate.

## Step-by-step protocol

### Step 0 — Locate the TUDOR pipeline root

Default candidate locations (check in order):

1. `C:\Users\nader\Downloads\calon_ukb_pipeline\` (the public GitHub repo `NaderGenedy/calon-ukb-pipeline`)
2. `D:\Projects\Lpa_Multilevel\TUDOR\` (if a project-specific working copy exists)
3. Ask the user explicitly if neither is present

Once located, glob for all `*.py` and `*.R` files relevant to TUDOR. The canonical TUDOR pipeline typically comprises:

- `00_*` — UK Biobank RAP extraction scripts
- `01_*` to `09_*` — data preparation, feature engineering, treatment adjustment, model fitting
- `TUDOR_*.R` / `TUDOR_*.py` — the final model fitting and validation
- `*_validation_*.py` — external validation in UKB-LC and Wales

Produce a file-by-file inventory before dispatching the agents — each agent reads only the files relevant to its frame to keep its context window focused.

### Step 1 — Dispatch the five agents in parallel

Use the Task tool (or equivalent) to spawn five independent subagents in a single message. Each agent receives:

- The pipeline root path
- Its specific verification remit (see "Agent specifications" below)
- A reporting template
- An explicit instruction to NOT modify any files (read-only audit)
- A 350–500 word output budget

Spawn all five in the same turn — do NOT spawn sequentially.

### Step 2 — Reconcile and produce the unified report

When all five agents return:

1. Read all five `agent*_*.md` reports
2. Cross-check for any conflicting findings (e.g., Agent 1 says cohort n = X, Agent 2 says n = Y — flag and resolve)
3. Aggregate any FAIL or HIGH-RISK findings into a single remediation list at the top
4. Aggregate the Agent 5 augmentation experiment into a dedicated table
5. Write `TUDOR_QC_REPORT.md` with: master verdict, per-agent verdict table, cornerstone numerical checks, remediation list, augmentation experiment results, three actionable issues to close before R2 / submission

Provide the user with: (a) the unified report path; (b) a 5-line chat summary; (c) the augmentation-table headline conclusions (does adding T2DM / age / LDLR-tier / ApoB clear the threshold for clinically meaningful improvement?).

## Agent specifications

### Agent 1 — Raw Data Integrity

**Remit.** Audit the raw UK Biobank extraction logic and cohort-construction code. Does the data that enters the Elastic Net model faithfully represent what the manuscript claims?

**Specific checks.**

- WRONG_FIELDS regression: confirm `p131286 / p131288 / p131290 / p131292 / p131294` are NOT used as ASCVD endpoint anchors (they code I10–I20 hypertension). For any outcome-based biological-validation sub-analysis, confirm correct fields (`p131296 / p131298 / p131306`) plus HES first-admission dates.
- Friedewald LDL: applied only when TG < 4.5 mmol/L; alternative method (Martin–Hopkins or NIH equation 2) for TG ≥ 4.5.
- Treatment-adjusted LDL formula: statin-intensity → expected LDL reduction table is documented, and the inverse transform reconstructs a clinically plausible pre-treatment LDL distribution (typical mean ~ 5.5 mmol/L for FH-positive).
- Cohort filters: exclusion criteria are explicit, ordered, and produce documented n-attrition figures.
- Data-leakage barriers: training and validation cohorts share no eids; family-level structure handled via `Is_Relative` flag and/or `FamilyNumber` deduplication.
- Variant-call set provenance: WES path, ClinVar / ACMG filter version, and pathogenicity threshold are reproducible.
- Statin prescription source: GP records (`p42039` linkage flag), with date and DDD per prescription.

**Output template.** PASS / WARN / FAIL per check above; cohort-construction n-attrition diagram in markdown table; one-paragraph verdict; up to 3 actionable issues.

### Agent 2 — Statistical Methodology

**Remit.** Audit the statistical machinery. Are AUC / calibration / NRI / Brier / DCA computed correctly?

**Specific checks.**

- Elastic Net hyperparameter selection: α-grid (typically {0.1, 0.5, 0.9, 1.0}) and λ-grid (100 values, log scale) selected via *k*-fold cross-validation; folds preserve outcome prevalence (stratified).
- LOCO-CV implementation: each Welsh Health Board acts as a fold; coefficients NOT refit when computing external validation in UKB (TRIPOD Type 4 = frozen Wales coefficients).
- AUC + DeLong 95 % CI: implemented correctly (paired test on the same patients when comparing TUDOR vs DLCN).
- Calibration intercept + slope: assessed on the LINEAR predictor scale (logit space), not on probability scale.
- Brier score: raw + scaled (against null-model Brier).
- NRI bilateral bootstrap: B ≥ 100; **both** baseline and augmented models refit per resample (not just one). Categorical and continuous NRI reported separately.
- IDI: computed as mean difference in predicted probabilities for cases minus mean for non-cases.
- DCA: net benefit reported at 5 %, 10 %, 20 % thresholds (and at the clinical operating point).
- DeLong *p* for ΔAUC: paired test, not unpaired.
- Multiple testing: BH-FDR if > 10 comparisons.
- E-value sensitivity (if any causal claim).
- Reproducibility: random seeds set at module top; pyarrow / lifelines / scikit-learn versions pinned.

**Output template.** PASS / WARN / FAIL per check; flagged any methodology gap relative to TRIPOD-AI 2024; one-paragraph verdict; up to 3 actionable issues.

### Agent 3 — Clinical Plausibility

**Remit.** Audit the clinical sensibility of every reconstructed variable and every model finding. Could a senior lipid clinician read this and find anything implausible?

**Specific checks.**

- Treatment-adjusted LDL distribution: median for FH-positive should be ~ 5.5 mmol/L; for FH-negative ~ 4.0 mmol/L. Implausible if either is < 3 or > 8.
- Trig_Filter direction: high triglycerides ↓ probability of monogenic FH (because high-TG shifts diagnostic separation from monogenic toward polygenic / metabolic-syndrome biology).
- Metabolic shield effect: Cohen's *d* of LDL separation drops as TG rises (~ 0.89 normal-TG → ~ 0.18 high-TG). If direction is reversed, FLAG.
- Gene-specific AUC hierarchy: APOB > PCSK9 > LDLR. This is mechanistic — APOB R3500/R3527 are molecularly homogeneous; LDLR has > 2,000 allelic variants spanning ligand-binding through cytoplasmic domains.
- Cardiovascular mortality gradient across treatment-adjusted LDL deciles: monotonic increase (e.g., 0.65 % → 1.76 % from D1 to D10).
- Index_Effect handling: cascade-screened relatives (`Is_Relative` = 1) have their LDL coefficient effectively zeroed (`Index_Effect = (1 − Is_Relative) × LDL_Untreated`), reflecting that LDL is already known to be high in cascade screening.
- Operating-point coherence: at the chosen Youden threshold, sensitivity / specificity / PPV / NPV are reported AND make clinical sense (e.g., screening threshold prioritises sensitivity).

**Output template.** PASS / WARN / FAIL; senior-clinician verdict ("would a lipid-clinic consultant accept these as plausible?"); up to 3 actionable issues.

### Agent 4 — Molecular / Genetic

**Remit.** Audit the variant-calling and FH-classification pipeline.

**Specific checks.**

- WES call set: variant caller version, joint-calling pipeline, QC filters (depth ≥ 10, GQ ≥ 20, MAF threshold).
- LDLR severity-tier classification: severe / moderate / mild / null definitions are explicit and consistent with published FH variant-effect-prediction literature (e.g., Tabet 2025 DMS uptake; Iacocca 2018 ClinVar consensus).
- APOB pathogenic variants: R3500Q (rs5742904) and R3527W (rs144467873) explicitly included.
- PCSK9 gain-of-function: D374Y (rs28942111) and other published GOF variants.
- LDLRAP1 (autosomal-recessive hypercholesterolaemia): biallelic carriers identified.
- LDLR copy-number variants: deletion / duplication detection from WES read-depth, ExomeDepth or equivalent.
- APOE rare pathogenic variants: included if claimed in cohort breakdown.
- ACMG pathogenicity filter: PVS1 / PS1 / PS3 / PM1 / PM2 / PP3 / BS1 / BP4 logic documented.
- AlphaMissense / REVEL / CADD thresholds: declared and applied consistently.
- Family-level structure: `FamilyNumber` field used to flag cascade-screened relatives; `Is_Relative` derived from it.
- Founder effects (especially Welsh / Cape Coloured / Lebanese / Christian Lebanese / French-Canadian) flagged if relevant.

**Output template.** PASS / WARN / FAIL; gene-by-gene n-positives table; one-paragraph verdict; up to 3 actionable issues.

### Agent 5 — Feature-Augmentation Experiment

**Remit.** Run the actual augmentation experiment using the bundled `scripts/feature_augmentation.py`. Does adding T2DM (binary), age (per decade), LDLR-tier (4-level), and/or ApoB (continuous, per SD) improve TUDOR's AUC, calibration, or NRI?

**Specific protocol.**

For each cohort (Wales registry development set, Wales external validation, UKB-LC validation):

1. Refit TUDOR baseline (Elastic Net, frozen feature set).
2. Refit augmented variants:
   - TUDOR + T2DM
   - TUDOR + age (per decade)
   - TUDOR + LDLR-tier (severe / moderate / mild / null; FH-negative as reference)
   - TUDOR + ApoB (per SD, continuous)
   - TUDOR + all four
3. For each variant, compute:
   - AUC + DeLong 95 % CI; ΔAUC vs baseline + DeLong *p*
   - Calibration intercept and slope (logit scale)
   - Brier scaled
   - NRI bilateral bootstrap (B = 100); refit baseline and augmented per resample
   - DCA net benefit at 5 / 10 / 20 % thresholds; ΔNB vs baseline

4. Pre-specified thresholds for "clinically meaningful improvement":
   - ΔAUC ≥ +0.01 absolute (DeLong *p* < 0.05) AND
   - Calibration slope shifts toward 1.0 by ≥ 0.05 AND
   - ΔNRI (bilateral) ≥ +5 percentage points AND
   - ΔNB at the operating threshold is positive

   Report which augmentation passes ALL four; any single failure means the augmentation should NOT be adopted.

5. Output table (rows = augmentation variants × cohorts; columns = AUC + ΔAUC + DeLong p, calibration slope + shift, NRI + 95% CI, ΔNB at three thresholds, verdict).

**Note on LDLR-tier.** This requires variant-call data on test patients. If the augmentation is being run on a development-style cohort, LDLR-tier is observable. If on a population-screening cohort (UKB-LC), LDLR-tier is observable only for variant carriers (~ 0.6 % of UKB); the rest are FH-negative by default. This means the LDLR-tier augmentation tests whether *carrier-status genotype* adds information to *clinical phenotype* — which is the conceptually correct experiment for an FH diagnostic algorithm.

**Output template.** Full augmentation table per cohort; ranked recommendation (which augmentation, if any, is worth incorporating into a TUDOR v2.0); one-paragraph verdict.

## Reconciliation step (NOT an agent — single-thread)

After all five agents return:

1. Open all five reports.
2. Cross-check headline numbers across agents (e.g., does Agent 1's n match Agent 2's n match Agent 5's n?).
3. Write `TUDOR_QC_REPORT.md` with this structure:

```
# TUDOR Full-Pipeline QC Report

**Date:** YYYY-MM-DD
**Pipeline root:** <path>
**Manuscript version under audit:** <e.g., JCLINLIPID-D-25-01142R1>

## Master Verdict
[PASS / PASS with minor remediation / FAIL with major issues]

## Per-Agent Verdict Table
| Agent | Verdict | Issues raised |
|---|---|---|
| 1. Raw Data Integrity | ... | ... |
...

## Cornerstone Numerical Checks
[Wales AUC; UKB-LC AUC; DLCN comparator AUC; gene-specific AUCs; etc.]

## Remediation List (prioritised)
[A. (CRITICAL): ...; B. (HIGH): ...; C. (MEDIUM): ...]

## Feature-Augmentation Experiment
[Full table from Agent 5; explicit recommendation]

## Reproducibility Command
[Path to run the audit again]
```

4. Print a 5-line chat summary for the user. Highlight any CRITICAL or HIGH issue. Lead with the augmentation recommendation if non-trivial.

## Output deliverables

For every TUDOR-QC pass:

1. `agent1_raw_data_audit.md` (350–500 words)
2. `agent2_stats_methodology_audit.md` (350–500 words)
3. `agent3_clinical_plausibility_audit.md` (350–500 words)
4. `agent4_molecular_genetic_audit.md` (350–500 words)
5. `agent5_augmentation_report.md` (350–500 words plus a wide table)
6. `TUDOR_QC_REPORT.md` (master report; 800–1,500 words)
7. `scripts/feature_augmentation.py` (the actual experiment script, reusable)
8. A 5–10 line chat summary

Default output location: `D:\Projects\Lpa_Multilevel\nejm-email\tudor_qc\<YYYY-MM-DD>\` (per-run subdirectory keyed by date).

## Common pitfalls (avoid)

- **Spawning agents sequentially instead of in parallel.** Defeats the independence assumption — Agent 2 may unconsciously align to Agent 1's framing. Always parallel.
- **Letting one agent dominate the narrative.** The unified report must give equal weight to each agent's domain. If Agent 5 (augmentation) returns a striking ΔAUC, do not let it eclipse Agent 1's raw-data-integrity issues — a + 0.02 AUC built on a buggy cohort is worse than the original.
- **Skipping the feature-augmentation experiment when the user only asks for "QC".** The user has explicitly requested augmentation as part of this skill's contract; always run Agent 5 unless the user explicitly opts out.
- **Modifying any TUDOR source file during the audit.** Read-only. If a bug is found, the remediation list documents it; the user decides whether to lock-rerun.
- **Letting Agent 5 declare an augmentation worthwhile on AUC alone.** The pre-specified four-criterion threshold (AUC + calibration + NRI + DCA) prevents AUC-only over-claiming, which is the dominant failure mode in prediction-modelling literature.

## Reference exemplar

When TUDOR moved from XGBoost (n = 1,051) to dual-validated Elastic Net (n = 65,274) for the R1 revision at the Journal of Clinical Lipidology, the pipeline expanded across UK Biobank Application 1002450 (n = 58,021; 3,223 FH-positive across 6 causative genes) and the All-Wales FH Registry (n = 7,253; 805 FH-positive). The R1 manuscript reports AUC 0.842 (Wales) and AUC 0.750 (UKB-LC), with DLCN comparator AUC 0.636 (DeLong *p* = 6.73 × 10⁻²⁴). The dual-validation architecture is TRIPOD Type 2b (Wales) plus TRIPOD Type 4 (UKB-LC, frozen Wales coefficients). This skill's first run on the R1 manuscript should reproduce these headline AUCs to four decimal places.

## Companion skills

- `ukb-inventory` — locate raw UK Biobank fields in the working tree before any audit step
- `ukb-preflight` — verify Python 3.12 + pyarrow + console encoding before running Agent 5
- `cox-analysis` — if any survival-outcome extension is being added to TUDOR (currently TUDOR is a classifier, not a Cox model)
- `manuscript-qc` — after this skill passes, run manuscript-qc to certify every prose claim against the locked source CSVs
- `reproduce-paper` — to build a frozen ground-truth test suite that locks the audited pipeline against future drift
