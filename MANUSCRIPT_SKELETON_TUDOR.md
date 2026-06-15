# MANUSCRIPT SKELETON — TUDOR Primary Report
## "A treatment-adjusted, ascertainment-aware prediction model for familial hypercholesterolaemia: development and bidirectional external validation in 506,506 adults"

> **Status:** skeleton / working scaffold. Prose to be drafted section-by-section.
> **House rules (all skills):** British English; numbers/criteria before adjectives; every number traces to a source row (verify against locked files, never memory); blunt honest limitations; leakage discipline; EPV respected; claims matched to design.
> **Authoritative artefacts:** poster (BioRender), Slide 10 (`TUDOR_coefficients_locked.csv`), the locked reproducer suite (66/68 PASS — see Open Item O-7).

---

## 0. SKILL MAP — which skill owns which section

| Section | Owning skill(s) |
|---|---|
| Title / framing options | structured-brainstorming |
| Estimand & target-trial | preventive-cardio-epidemiologist |
| Model, validation, stats | cardiometabolic-biostatistician |
| Provenance & reproducibility | reproducibility-engineer |
| Comparator benchmarks | cardiometabolic-evidence-synthesis |
| Reviewer pre-empt | lipid-cardiology-reviewer |
| Prose architecture & journal | academic-medical-writer |
| Figures/tables | scientific-figure-designer |
| Submission roadmap | research-executive-planner |
| Future work / ceiling | autoresearch |

---

## 1. TITLE & FRAMING OPTIONS  *(structured-brainstorming — diverge→converge)*

**Ranked title candidates** (pick #1; keep #2 as cascade-journal alt):

1. **#1 — "Treatment-adjusted, ascertainment-aware prediction of familial hypercholesterolaemia: development and bidirectional external validation in 506,506 adults."** *(generator: first-principles — names the two novel mechanisms + the scale; declarative.)*
   - **Why it wins:** states the contribution (treatment + ascertainment adjustment) and the validation scale in one line; survives a structure-conscious editor.
   - **Biggest risk:** "ascertainment-aware" may read as jargon to a general journal — define in first sentence of abstract.
2. **#2 —** "Why classic FH criteria fail the treated patient: a 506,506-adult validation of a treatment-adjusted prediction model." *(generator: provocation — leads with the gap.)*
3. **#3 —** "Recovering the untreated LDL phenotype to find familial hypercholesterolaemia at population scale." *(generator: mechanism-forward.)*

**Kill list:** "machine-learning FH detector" (over-claims ML; it's penalised logistic — defensibility ✗); "TUDOR score" alone (topic, not finding ✗); anything with "AI" (reviewer-bait ✗).

**One-sentence message (lock this, every section serves it):**
> *In treated, metabolically complex, often cascade-ascertained modern patients, a treatment-adjusted, ascertainment-aware penalised-logistic model predicts monogenic FH carrier status better than every established clinical criterion, in every prespecified subgroup, in both validation directions — and refines who should proceed to genetic testing.*

---

## 2. ESTIMAND & TARGET-TRIAL FRAME  *(preventive-cardio-epidemiologist)*

**The prediction estimand (one sentence — write it before anything else):**
> *In adults referred or screened for hyperlipidaemia, the probability of carrying a pathogenic/likely-pathogenic LDLR variant (gold standard), estimated from routinely available lipid, treatment, demographic and ascertainment variables measured at/around the index lipid panel.*

| Target-trial element | Specification | Note |
|---|---|---|
| Population | Welsh FH registry (specialist) + UK Biobank (population) | two ascertainment worlds — the design's whole point |
| Index / time-zero | first qualifying lipid panel | confirm no post-baseline leakage (Open Item O-3) |
| "Exposure"/predictors | 11 locked variables (§3) + engineered Trig_Filter, Index Effect | all measured at/before index |
| Comparators | DLCN, Simon Broome, MEDPED (published rules, **no refit**); FAMCAT not computable | §6 |
| Outcome (label) | LDLR pathogenic/LP carrier status (ClinVar) | NOT ASCVD — this is case-finding, not event risk |
| Estimand type | predictive (discrimination + calibration), not causal | **verb discipline:** "predicts/identifies", never "causes" |
| Validation | bidirectional external (registry↔population), TRIPOD Type 1b→3→4 | §7 |

**Causal-claim gate:** this is a *prediction* paper. No causal verbs anywhere. The Trig_Filter "biology" (Slide 8–9) is supporting mechanism, framed as *why the predictor works*, not a causal effect estimate.

---

## 3. THE LOCKED MODEL — full specification  *(biostatistician + reproducibility-engineer)*

**Model form:** elastic-net penalised logistic regression. **Hyperparameters (locked):** α = 0.5, C = 1.0, saga solver, seed 20260518. **Intercept:** −1.851. **Source of truth:** `TUDOR_coefficients_locked.csv`.

`logit(P[carrier]) = −1.851 + Σ βⱼ·zⱼ`  (zⱼ = standardised predictor; βⱼ per 1-SD below)

| # | Variable (locked name) | Definition to write in Methods | β (1-SD) | Conceptual idea |
|---|---|---|---|---|
| 1 | `ldl_ul` | back-calculated **untreated** LDL (drug×dose×duration×adherence reconstruction) | **+1.609** | Lipid Age |
| 2 | `on_statin` | on lipid-lowering therapy at index (binary) | +0.357 | Lipid Age |
| 3 | `tc_chem` | total cholesterol (direct assay) | −0.743 | Triglyceride Shield |
| 4 | `tg_chem` | triglycerides (direct assay) | −0.604 | Triglyceride Shield |
| 5 | `non_hdl` | non-HDL cholesterol | −0.151 | Triglyceride Shield |
| 6 | `hdl_chem` | HDL cholesterol | −0.137 | Triglyceride Shield |
| 7 | `age_per_decade` | age, per decade | −0.357 | Demographic |
| 8 | `sex_F` | female sex (binary) | −0.070 | Demographic |
| 9 | `premature_mace` | premature MACE (age-/sex-defined) | +0.064 | Proband Effect |
| 10 | `t2dm_int` | type 2 diabetes indicator/interaction | +0.047 | Proband Effect |
| 11 | `fam_hist_cvd` | family history of CVD | **0.000 (L1-eliminated)** | Proband Effect |

**Engineered constructs to define explicitly in Methods (don't bury these — they are the novelty):**
- **Lipid Age / untreated-LDL reconstruction:** `ldl_ul` from published intensity bands — atorvastatin 25–48%, rosuvastatin 35–55%, simvastatin 20–42%, pravastatin 15–29%, fluvastatin 15–22%; ezetimibe +20%; bempedoic acid +25%; PCSK9i +65%; **capped at 85%**.
- **Trig_Filter** = `ldl_ul / (tg_chem + 0.1)` — metabolic-purity signal separating isolated-LDL FH from mixed dyslipidaemia.
- **Index Effect** = `(1 − Is_Relative) × ldl_ul` — ascertainment adjustment (proband vs cascade relative).

> **⚠ Reconciliation flag (O-1):** Slide 10 lists 11 named coefficients; the poster Methods also name `Trig_Filter` and `Index Effect` as engineered features. Confirm against `TUDOR_coefficients_locked.csv` whether Trig_Filter/Index Effect are (a) separate fitted terms, (b) encoded within the 11, or (c) pre-processing transforms feeding `ldl_ul`. The Methods variable list must match the locked CSV exactly — do not state from the slide alone.

**EPV check (biostatistician):** development Design A train n=3,099 with 554 carriers → ~50 events/predictor across 11 terms; comfortably above the EPV floor. State pmsampsize logic in Methods.

---

## 4. ABSTRACT SKELETON  *(academic-medical-writer — write LAST, polish MOST)*

Structured to target-journal headings. Every number carries its CI; conclusion states only what design licenses.

- **Background:** classic FH criteria calibrated on untreated probands; modern patients treated/metabolically complex/cascade-ascertained → criteria degrade.
- **Methods:** elastic-net model (11 variables, treatment-adjusted untreated LDL + ascertainment terms); developed in Welsh registry; **bidirectional** external validation vs DLCN/Simon Broome/MEDPED; gold standard LDLR carrier status; TRIPOD-AI.
- **Results (fill from locked outputs):** Welsh validation AUC **0.770** (vs DLCN 0.670, Δ+0.099, p<0.001; SB 0.570; MEDPED 0.553); full UK Biobank AUC **0.631** (vs DLCN 0.538, SB 0.510, MEDPED 0.520; all p<0.001); **25/25** prespecified subgroups won vs every comparator (11 Welsh + 14 UKB); type-2-diabetes subgroup **0.718 vs 0.538**, n=45,715; matched-positivity sensitivity AND specificity both exceed all comparators; ladder 0.770/0.732, pooled-frozen 0.746, recalibrated 0.756.
- **Conclusion:** treatment-adjusted, ascertainment-aware prediction outperforms established criteria across settings and subgroups; **refines pre-test probability for genetic testing — does not replace it.**
- **[Add CIs to every AUC from source before submission — Open Item O-2.]**

---

## 5. INTRODUCTION SKELETON  *(academic-medical-writer — 3 moves, ~4 paragraphs)*

1. **What's known:** FH ~1/250; lethal-if-missed, manageable-if-found (Luirink NEJM 2019); cascade screening leverage (NNT≈2).
2. **The gap:** criteria (Simon Broome 1991, MEDPED 1993, DLCN/WHO 1999) calibrated on untreated probands; statins/ezetimibe/PCSK9i compress LDL; cascade ascertainment + mixed dyslipidaemia distort the classic signature → criteria fail exactly where modern patients live (forward-reference the 0.435 DLCN-in-T2DM result).
3. **The specific aim:** develop a treatment-adjusted, ascertainment-aware prediction model recovering the untreated-LDL phenotype, and test it head-to-head, pre-registered, bidirectionally, against every established criterion.
4. **Objective sentence (explicit):** *"We aimed to develop and externally validate, in 506,506 adults, a prediction model for monogenic FH that is robust to treatment and ascertainment, and to compare it head-to-head with DLCN, Simon Broome and MEDPED."*

---

## 6. METHODS SKELETON  *(biostatistician + evidence-synthesis + reproducibility-engineer)*

- **6.1 Design & reporting:** development + bidirectional external validation; **TRIPOD-AI**; PROBAST self-appraisal; UK Biobank Application 1002450.
- **6.2 Cohorts & ascertainment:**
  - Welsh FH registry (PASS) — specialist; **family-level deduplication** (`wales_family_deduplicated_*`), zero shared family identifiers between development and validation.
  - UK Biobank — population; LDLR coding-variant carriers by WES.
  - **⚠ O-4: state every N from the live file.** Reconcile: development train n=3,099 (554 carriers, Slide 10) vs Welsh validation Design A n=1,471 (333 carriers, poster) vs total 506,506. Make the CONSORT/flow numbers reconcile (Figure 1).
- **6.3 Gold standard:** ClinVar pathogenic/likely-pathogenic LDLR (state version; note APOB/PCSK9 handling; acknowledge CNV blind spot — O-5).
- **6.4 Predictors & engineered features:** §3 table verbatim from locked CSV; untreated-LDL reconstruction bands; Trig_Filter; Index Effect; missing-data handling (multiple imputation; ±inf→NaN sanitisation guard).
- **6.5 Model development:** elastic-net (α=0.5, C=1.0, saga, seed 20260518); standardisation; L1 elimination of `fam_hist_cvd`; coefficients **frozen** then transported (no refit in validation).
- **6.6 Comparators (evidence-synthesis — retrieve, don't recall):** DLCN, Simon Broome, MEDPED as **published rules, no refit**; FAMCAT not computable (coefficient supplement unavailable — state explicitly); note electronic-DLCN cannot use xanthomata/arcus/family-history in population data (the structural-fairness caveat — must be in Limitations too).
- **6.7 Validation design (the headline):** bidirectional — Welsh→UKB and UKB→Welsh; validation ladder (Design A primary, A reciprocal, B pooled-frozen, C recalibrated, C-direct frozen, D full-UKB TRIPOD-4).
- **6.8 Statistics:** discrimination (AUC + **95% CI**, DeLong paired); calibration (slope/intercept + plot); decision-curve net benefit at 5/10/20%; matched-positivity sensitivity/specificity; family-block bootstrap; pre-registered 25-subgroup head-to-head with locked sign convention (ΔAUC = TUDOR − comparator, positive favours TUDOR).
- **6.9 Reproducibility:** single entry point; `_data_paths` resolver (no hard-coded paths); provenance ledger (raw column→variable); QC gate 0-FAIL; seed pinned; every reported number traceable. Reproducer status 66/68 — **resolve the 2 failures before submission (O-7).**

---

## 7. RESULTS SKELETON  *(map 1:1 to figures/tables; state, don't interpret)*

- **7.1 Cohorts & case-mix** (Table 1; Figure 1 flow) — the "contemporary patient": T2DM 27%, BMI ~29, TG >2, statin 62–81%.
- **7.2 Primary discrimination — Welsh validation:** TUDOR 0.770 vs DLCN 0.670 (Δ+0.099, p<0.001), SB 0.570 (Δ+0.200), MEDPED 0.553 (Δ+0.216). *(Figure 2: ROC + comparator operating points.)*
- **7.3 Population stress test — full UK Biobank:** TUDOR 0.631 vs DLCN 0.538, SB 0.510, MEDPED 0.520 (all p<0.001). **Frame as strength:** every comparator at/near chance; advantage persists. *(Figure 3b.)*
- **7.4 Validation ladder:** bidirectional 0.770/0.732; pooled-frozen 0.746; recalibrated 0.756. *(Figure 3a forest.)*
- **7.5 Subgroups — the clincher:** 25/25 prespecified won vs every comparator; 11/11 Welsh, 14/14 UKB; win/tie/loss 25/0/0 ×3 comparators. *(Figure 3c heatmap — **verify cell values vs source, O-6**.)*
- **7.6 Type-2-diabetes signal:** TUDOR 0.718 vs eDLCN 0.538, n=45,715 — supports Trig_Filter biology.
- **7.7 Operating characteristics:** matched-positivity sensitivity AND specificity exceed all comparators simultaneously.
- **7.8 Calibration & net benefit:** slope/intercept; DCA at 5/10/20%. *(Figures 4–5.)*

---

## 8. DISCUSSION SKELETON  *(academic-medical-writer + lipid-cardiology-reviewer)*

1. **Key finding (one sentence):** treatment-adjusted, ascertainment-aware prediction beats every classic criterion in every subgroup, both directions.
2. **Mechanism (why it works):** Lipid Age recovers the suppressed phenotype; Trig_Filter exploits Brown–Goldstein selective-LDLR biology (isolated LDL, intact TG pathway) — Slides 8–9 evidence.
3. **Comparison with literature:** contrast with DLCN/SB/MEDPED degradation in treated/ascertained populations; position vs EHR/ML FH case-finding (retrieve prior art — evidence-synthesis, O-8).
4. **Clinical implication:** triage to genetic testing; the three use-cases (treated patient, cascade relative, mixed dyslipidaemia); 3-column framing (lipid clinic / cardiologist-diabetics / primary-care no-xanthoma-needed).
5. **Limitations — OWN IT (do not bury):**
   - **O**bservational/retrospective — no prospective impact trial yet.
   - **W**eights frozen — population absolute risk needs local recalibration.
   - **N**ot all populations — UKB healthy-volunteer bias; non-UK/diverse validation needed.
   - **I**ncomplete eDLCN — structurally generous to TUDOR (no xanthomata/FH-history in population data); say so plainly.
   - **T**reatment reconstruction approximate — ±10% factor sensitivity tested; equity caveat (undertreatment ≠ lower risk); CNV gold-standard blind spot.
6. **Conclusion:** refines pre-test probability for genetic testing; does not replace it; prospective evaluation is the next step.

---

## 9. FIGURE & TABLE PLAN  *(scientific-figure-designer — one figure, one claim; regenerate from committed outputs)*

| Item | Claim it makes | Chart | Source (to confirm) | Honesty watch |
|---|---|---|---|---|
| Fig 1 | cohorts & exclusions reconcile | CONSORT flow | flow CSV | numbers = analytic N |
| Fig 2 | TUDOR > comparators (Welsh) | ROC + operating points, AUC+CI | `roc_welsh.csv` | AUC alone ≠ utility → pair w/ Fig 4 |
| Fig 3 | advantage persists & is universal | (a) ladder forest (b) UKB bars (c) subgroup heatmap | `ladder.csv`, `subgroups.csv` | colour-blind safe; annotate N |
| Fig 4 | well-calibrated | calibration plot + slope/intercept | `calibration.csv` | show line, not just H-L p |
| Fig 5 | clinically useful | decision curve, net benefit | `dca.csv` | treat-all/none reference lines |
| Table 1 | contemporary case-mix | baseline characteristics | master file | report N per cell |
| Table 2 | locked model | the 11 β + engineered features | `TUDOR_coefficients_locked.csv` | verbatim from locked CSV |

Defaults: vector (PDF/SVG); Okabe–Ito/viridis; ≥7 pt sans; single ≈90 mm / double ≈180 mm; one script regenerates all from committed outputs.

---

## 10. COMPARATOR / EVIDENCE TABLE TO FILL  *(cardiometabolic-evidence-synthesis — retrieve, never recall)*

Build claim→source map before drafting Intro/Discussion. Rows to populate from full text (scite/PubMed), not memory:

| Claim needing a source | Comparator / fact | Status |
|---|---|---|
| FH prevalence ~1/250 | Beheshti Circulation 2020 | ✅ have |
| Childhood-statin outcomes | Luirink NEJM 2019 | ✅ have |
| DLCN/Simon Broome/MEDPED published rules + thresholds | original criteria papers | ⬜ extract exact thresholds |
| FAMCAT variable set (why not computable) | FAMCAT papers | ⬜ document coefficient gap |
| Prior EHR/ML FH case-finding benchmarks | Akioyamen, Banda et al., etc. | ⬜ prior-art + benchmark AUCs |
| Brown–Goldstein selective-LDLR biology | 1985 Nobel work | ✅ cite |
| **"Could not be source"** list | — | ⬜ compile honestly |

---

## 11. ANTICIPATED REVIEWER ATTACKS + PRE-EMPTS  *(lipid-cardiology-reviewer)*

| 🔴/🟠 | Likely attack | Pre-empt built into the paper |
|---|---|---|
| 🟠 | "UKB AUC 0.631 is modest" | frame vs comparators at chance; report it as a population stress test, not the headline registry result |
| 🟠 | "eDLCN comparison is unfair (truncated)" | state explicitly in Methods §6.6 + Limitations; argue it's the realistic deployment setting |
| 🟠 | "It's logistic regression, not ML" | call it penalised regression; interpretability = safety feature; don't over-claim |
| 🟠 | "Untreated-LDL reconstruction is approximate" | ±10% factor sensitivity analysis; report robustness |
| 🟠 | "Polygenic hypercholesterolaemia in controls" | makes discrimination conservative; state it |
| 🟠 | "Gold standard misses CNVs" | acknowledge; biases toward null |
| 🟡 | "Family clustering inflates significance" | family-level dedup + family-block bootstrap |
| 🟡 | "Multiplicity across 25 subgroups" | pre-registration + locked sign convention; primary estimand designated |

---

## 12. REPRODUCIBILITY / PROVENANCE CHECKLIST  *(reproducibility-engineer — 0-FAIL before submission)*

- [ ] One entry point regenerates all outputs.
- [ ] `_data_paths` resolver; no hard-coded paths.
- [ ] Every reported N recomputed from live file (not stale constants).
- [ ] Provenance ledger: each of the 11 variables → raw column → transform.
- [ ] Leakage scan: all predictors pre-index; outcome not among predictors.
- [ ] UKB outcome provenance correct (if any ASCVD term used downstream): I20/I21/I25, never hypertension fields p131286–296.
- [ ] Seed pinned (20260518); env pinned; `run_meta.json` emitted.
- [ ] **Reproducer 66/68 → resolve the 2 failures (O-7).**
- [ ] Every manuscript/abstract/figure number reconciles to source CSV.

---

## 13. SUBMISSION ROADMAP  *(research-executive-planner — critical path)*

**Target + cascade:** Primary — *Circulation* or *EHJ* (strong CV science, population scale + mechanism). Cascade — *J Clin Lipidol* → *Heart/EJPC*. (Confirm fit; *Lancet*-family only if the prospective-impact angle is added.)

**Critical path (work backwards from submission date D):**
1. Resolve O-1 (model spec) + O-7 (reproducers) — **blocking, do first.**
2. Lock all numbers + CIs from source (O-2) → tables/figures regenerate.
3. Reconcile Ns + CONSORT (O-4).
4. Evidence/comparator table (§10) → draft Intro/Discussion.
5. Internal review (lipid-cardiology-reviewer pass) → revise.
6. Cover letter + declarations + TRIPOD-AI checklist → submit.

**WIP discipline:** finish this paper to submission-ready before opening the prospective-validation follow-up.

---

## 14. FUTURE WORK / HONEST CEILING  *(autoresearch)*

- **Front-runner next experiments (ranked by expected gain ÷ cost):** prospective lipid-clinic deployment (impact-defining); non-UK/diverse external validation; ApoB-gating of the Trig_Filter in borderline cases (the dyslipidaemia-tree integration); CNV-inclusive gold standard.
- **Logged as low-yield (don't re-spend):** adding more lipid ratios beyond the locked set (L1 already eliminated `fam_hist_cvd`); pushing UKB AUC by relaxing the honest validation = fabrication line, forbidden.
- **Honest ceiling statement:** registry discrimination ~0.77 is near the achievable frontier for routine-variable case-finding; further gain likely needs new measurements (ApoB/PRS/imaging), not more modelling of the same table.

---

## OPEN VERIFICATION ITEMS (resolve before any number is "final")

| ID | Item | Why it blocks |
|---|---|---|
| O-1 | Confirm model spec vs `TUDOR_coefficients_locked.csv` (are Trig_Filter/Index Effect separate fitted terms?) | Methods variable list must be exact |
| O-2 | Pull 95% CIs for every AUC from source | abstract/results incomplete without |
| O-3 | Leakage/temporal check on all 11 predictors | reviewer-fatal if any post-index |
| O-4 | Reconcile Ns: 3,099 train / 1,471 Design-A / 506,506 total / 333 vs 554 carriers | CONSORT must reconcile |
| O-5 | ClinVar version + APOB/PCSK9 + CNV handling | gold-standard definition |
| O-6 | Verify Fig 3c subgroup heatmap cell values vs source | per-cell numbers unconfirmed |
| O-7 | Resolve 2/68 failing reproducers | 0-FAIL gate before submission |
| O-8 | Prior-art / benchmark retrieval for Discussion | novelty claim + comparison |

---

*Skeleton built applying all 10 portfolio skills. Every headline number above is from the poster + Slide 10 (`TUDOR_coefficients_locked.csv`); confidence intervals and per-cell subgroup values must be drawn from the locked source files before drafting prose.*
