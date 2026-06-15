---
name: cox-analysis
description: Use when running any Cox proportional-hazards analysis on UKB / FH / cardiovascular cohort data. Triggers on "Cox PH", "survival analysis", "hazard ratio", "incident ASCVD", "risk model", "leak-free validation", "bi-external validation", or before any FH risk-prediction work. Enforces methodology TDD checklist: leak-free splits, family-level deduplication, NoAgeLDL sensitivity, OPCS-4/ICD-10 completeness, TRIPOD-compliant reporting.
---

# Cox-analysis — TDD-style methodology checklist

## When to invoke

Use this skill before commissioning ANY Cox proportional-hazards analysis or risk-prediction model on cardiovascular / FH / UKB data. Specifically when the user asks:
- "run Cox PH" / "survival analysis" / "hazard ratio"
- "incident ASCVD" / "composite outcome"
- "build a risk model" / "predict cardiovascular events"
- "bi-external validation" / "leave-one-cohort-out"
- "FH risk prediction" / "CALON" / "SAFEHEART comparison"
- "validate against [comparator]"

## TDD-style methodology checklist

Before running the analysis, write `methodology_checks.py` with executable assertions for each item below. Run the analysis, then run `methodology_checks.py`. If any assertion fails, autonomously diagnose, fix, and re-run. Loop until all checks pass green.

### 1 — Leak-free cohort splits

```python
def assert_no_train_test_leakage(train_ids, test_ids):
    overlap = set(train_ids) & set(test_ids)
    assert len(overlap) == 0, f"LEAKAGE: {len(overlap)} ids in both train and test"
```

For UK Biobank: split by `eid`. For multi-cohort (e.g. SW Dragon-3 ↔ Wales PASS): split by `DatabaseNumber` AND `FamilyNumber` (see check 2).

### 2 — Family-level deduplication

```python
def assert_family_deduped(train_df, test_df, family_col='FamilyNumber'):
    train_fams = set(train_df[family_col].dropna()) - {''}
    test_fams = set(test_df[family_col].dropna()) - {''}
    overlap = train_fams & test_fams
    assert len(overlap) == 0, (
        f"FAMILY LEAKAGE: {len(overlap)} families in both train and test. "
        "South Wales (Dragon-3) IS a subset of All-Wales PASS by FamilyNumber. "
        "Drop overlapping families from test before validation."
    )
```

Default behaviour: when both Dragon-3 and Wales PASS are in the same analysis, run `dedup vs SW: removed N rows by patient ID + M more by family ID; total removed: N+M`.

### 3 — NoAgeLDL sensitivity variant

```python
def assert_no_age_ldl_variant_present(model_variants):
    has_no_age_ldl = any('NoAgeLDL' in v for v in model_variants)
    assert has_no_age_ldl, (
        "MISSING SENSITIVITY: NoAgeLDL variant required. "
        "Tests whether novel predictors carry signal beyond the dominant "
        "traditional age + LDL backbone."
    )
```

Always include a model variant that drops age and LDL — this is the "what does CALON add beyond age+LDL?" question.

### 4 — Endpoint construction provenance (UKB-specific)

```python
WRONG_FIELDS = {'p131286', 'p131288', 'p131290', 'p131292', 'p131294'}  # I10-I20 (HTN)
RIGHT_FIELDS = {'p131296', 'p131298', 'p131306'}                          # I20, I21, I25

def assert_endpoint_uses_correct_fields(endpoint_construction_code):
    used_wrong = WRONG_FIELDS & set(endpoint_construction_code)
    assert not used_wrong, (
        f"ENDPOINT BUG: fields {used_wrong} code I10-I20 (hypertension), NOT ASCVD. "
        f"Use {RIGHT_FIELDS} plus HES first dates instead."
    )
```

### 5 — OPCS-4 / ICD-10 code list completeness

```python
EXPECTED_OPCS_FOR_SEVERE_AS = {'K611', 'K261', 'K262', 'K263', 'K264'}  # K611 = balloon valvuloplasty
EXPECTED_OPCS_FOR_PCI = {'K49', 'K50', 'K75'}
EXPECTED_OPCS_FOR_CABG = {'K40', 'K41', 'K42', 'K43', 'K44', 'K45', 'K46'}
EXPECTED_ICD_FOR_ASCVD_COMPOSITE = {'I20', 'I21', 'I22', 'I25', 'I63', 'G45', 'I70', 'I73', 'I74'}

def assert_code_list_complete(used_codes, expected_codes, context):
    missing = expected_codes - set(used_codes)
    assert not missing, f"MISSING CODES for {context}: {missing}"
```

### 6 — Cox proportional-hazards assumption check

```python
def assert_PH_assumption_holds(cox_fit):
    schoenfeld = cox_fit.check_assumptions(...)
    failing = schoenfeld[schoenfeld['p'] < 0.05]
    if len(failing) > 0:
        raise AssertionError(
            f"PH violated for: {list(failing.index)}. "
            "Use stratified Cox or time-varying coefficient or report HR with caution."
        )
```

### 7 — Competing risks (where applicable)

For elderly / FH cohorts where non-CV death is non-trivial: report Fine-Gray sub-distribution hazard alongside the standard Cox HR.

### 8 — TRIPOD-compliant reporting

For any prediction model:

```python
def assert_tripod_compliant_metrics_reported(results):
    required = {
        'C_statistic', 'C_statistic_CI',
        'cal_intercept', 'cal_slope',
        'brier', 'brier_scaled',
        'NRI', 'NRI_CI',
        'IDI',
        'DCA_NB_5pct', 'DCA_NB_10pct', 'DCA_NB_20pct',
        'DeLong_p_vs_comparator',
    }
    missing = required - set(results.keys())
    assert not missing, f"TRIPOD INCOMPLETE: missing {missing}"
```

NRI: bilateral bootstrap, B ≥ 100 (refit BOTH baseline and augmented Cox per resample). Comparator: frozen-coefficient external (TRIPOD Type 4) where possible.

### 9 — Multiple-testing correction

For any panel of > 10 simultaneous tests (NMR fields, PRS panel, subgroup forest):

```python
from statsmodels.stats.multitest import multipletests
fdr_qs = multipletests(p_values, alpha=0.05, method='fdr_bh')[1]
```

### 10 — E-value sensitivity

For every headline OR / HR / RR:

```python
def evalue_RR(rr, lower_CI):
    """VanderWeele & Ding 2017 E-value formula."""
    e_point = rr + (rr * (rr - 1))**0.5 if rr >= 1 else 1/rr + (1/rr * (1/rr - 1))**0.5
    e_lower = lower_CI + (lower_CI * (lower_CI - 1))**0.5 if lower_CI >= 1 else 1/lower_CI + (1/lower_CI * (1/lower_CI - 1))**0.5
    return {'point': e_point, 'lower_CI': e_lower}
```

Threshold: E-value > 2.0 for "robust to plausible unmeasured confounding".

## Pre-flight environment check

Before running the analysis pipeline:

```bash
python --version | grep -q "3.12" || echo "WARN: Python 3.12 required for pyarrow"
python -c "import pyarrow; import lifelines; import statsmodels; print('OK')" || pip install pyarrow lifelines statsmodels
```

Windows console: use ASCII-only output (no Unicode arrows / emoji). Set `PYTHONIOENCODING=utf-8` for any subprocess call.

## Output discipline

For every Cox analysis, save:
1. `*_cox_results.csv` — one row per Cox specification, with HR, lo, hi, p, n, events
2. `*_calibration.csv` — calibration intercept + slope per direction
3. `*_NRI.csv` — bilateral bootstrap median + CI
4. `*_subgroups.csv` — pre-specified subgroup forest
5. `*_evalues.csv` — E-value table
6. `*_provenance.csv` — file:row:n:events:date for every headline number
7. `methodology_checks.py` — executable assertion list
8. `*_TRACEABILITY_REPRODUCER.py` — auto-runnable audit (see `manuscript-qc` skill)

## Reference exemplars

- **CALON-9 leak-free 4-cohort validation** (May 2026): `C:/Users/nader/Downloads/calon_ukb_pipeline/CALON_NEW_v2.R` — 11 model variants × 9 directions × 2 outcomes; family-level dedup verified; NoAgeLDL sensitivity present; achieved 110/110 PASS in per-paper reproducer audit.
- **LOCK_RERUN exemplar** for endpoint-bug correction: `D:/Projects/Lpa_Multilevel/nejm-email/LOCK_RERUN_composite_ASCVD.py` — corrected p131286 → p131298 endpoint mapping, refreshed Cox + NRI(B=100) end-to-end, total runtime 2.5h.
- **CALON_NEWv2_side_by_side.csv** schema: `test_cohort, outcome_label, model, n_predictors, auc, ci_low, ci_high, delta_vs_SRE, cal_slope, n_test`.

## Common pitfalls (catch automatically)

- "FH paradox" — if FH stratum HR for a biomarker is null while non-FH HR is highly significant, suspect a UKB field-mapping bug (the original `v21_*` chain mapped p131286→I21 when those fields code hypertension).
- "Treatment paradox" — high observed re_LDL after a previous event reflects post-event drug intensification, not baseline severity. For incident-event Cox: include `prior_ascvd` as a covariate (matches SAFEHEART original design).
- "Quasi-separation" — when prior_ascvd predicts outcome perfectly in a subset, glmnet returns SE=NaN. Switch to Firth's penalised likelihood (`logistf` package).
- "Lipid-clinic ceiling" — UKB AUC plateaus around 0.72-0.74 for ASCVD prediction; this is biology (population screening), not model failure. Confirm by restricting to severe-phenotype subset (re_LDL ≥ 4.9 mmol/L) and noting if AUC does not rise.
