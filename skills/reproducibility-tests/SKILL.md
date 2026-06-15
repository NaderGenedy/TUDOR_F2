---
name: reproducibility-tests
description: Use when encoding statistical/methodological invariants as pytest tests so they run continuously rather than being re-checked manually each session. Triggers on "build a test harness", "continuous QC", "lock the methodology", "encode my invariants", "regression tests for the analysis", or when methodological gaps keep slipping through (family-level dedup, NoAgeLDL variant, OPCS K611, leak-free splits, TRIPOD completeness). Converts the 110/110 reproducer pattern into a pytest suite that fails CI when an invariant breaks.
---

# Reproducibility-tests — methodological invariants as executable assertions

## When to invoke

Use this skill whenever the user says any of:
- "build a test harness" / "encode the invariants"
- "lock the methodology" / "continuous QC"
- "regression tests for the analysis"
- "make the reproducer enforceable"
- "I keep having to flag the same methodological gaps"
- "TDD for my Cox analysis"

This skill converts the **one-shot `*_TRACEABILITY_REPRODUCER.py` pattern** from `manuscript-qc` into a **pytest suite** that runs every time the analysis re-runs, blocking regressions automatically. Manuscript-qc gives you 110/110 once; this skill gives you 110/110 forever.

## When NOT to invoke

- Just to verify a single number once → use `manuscript-qc` instead.
- For figure regeneration → use `reproduce-paper` instead.
- For pre-analysis methodology checks (single-pass) → use `cox-analysis` instead.

## Layout

```
<paper>/tests/
├── conftest.py                          # shared fixtures (cohort loaders, expected Ns)
├── test_01_cohort_invariants.py         # cohort sizes, family dedup, leak-free splits
├── test_02_endpoint_provenance.py       # OPCS/ICD code completeness, UKB field mapping
├── test_03_cox_model_invariants.py      # HR ranges, PH assumption, NoAgeLDL present
├── test_04_calibration_invariants.py    # TRIPOD-compliant metrics present + within range
├── test_05_subgroup_invariants.py       # subgroup forest n's match expected
├── test_06_manuscript_claims.py         # each claim in PROVENANCE_TABLE.csv matches CSV
├── test_07_figure_regeneration.py       # figures regenerate to byte-identical (allow tolerance for non-deterministic axes)
└── pytest.ini                            # config — fail fast, verbose, no warnings
```

## Five-step protocol

### Step 1 — Audit the current analysis for invariants

Scan the analysis scripts in the repo and propose 50–80 candidate tests. Categorise as:
- **Cohort invariants** (Ns, dedup, exclusions, leak-free)
- **Code-list invariants** (OPCS/ICD completeness)
- **Model invariants** (NoAgeLDL present, PH assumption, HR plausibility)
- **Calibration invariants** (TRIPOD-required metrics present + within plausible range)
- **Manuscript-claim invariants** (every number in PROVENANCE_TABLE.csv equals its CSV source)
- **Figure invariants** (regenerate byte-identical or within tolerance)

Save the proposal as `tests/INVARIANT_AUDIT.md` before writing any test code, and ask the user to confirm or amend the list.

### Step 2 — Build `conftest.py` with shared fixtures

```python
# tests/conftest.py
import pytest
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
RESULTS = ROOT / 'results'

@pytest.fixture(scope='session')
def cohort():
    return pd.read_csv(DATA / 'analytical_cohort.csv')

@pytest.fixture(scope='session')
def cox_results():
    return pd.read_csv(RESULTS / 'cox_results.csv')

@pytest.fixture(scope='session')
def provenance_table():
    return pd.read_csv(RESULTS / 'PROVENANCE_TABLE.csv')
```

### Step 3 — Encode the canonical invariants

```python
# tests/test_01_cohort_invariants.py
EXPECTED_N_PRE = 63939
EXPECTED_N_POST = 165276
EXPECTED_N_TOTAL = 273036
TOLERANCE = 10  # row-count tolerance for QC merges

def test_premenopausal_n(cohort):
    n = (cohort['menopausal_status'] == 'Premenopausal').sum()
    assert abs(n - EXPECTED_N_PRE) <= TOLERANCE, f"n_pre={n}, expected {EXPECTED_N_PRE}"

def test_no_family_leakage(train_cohort, test_cohort):
    train_fams = set(train_cohort['FamilyNumber'].dropna()) - {''}
    test_fams = set(test_cohort['FamilyNumber'].dropna()) - {''}
    overlap = train_fams & test_fams
    assert not overlap, (
        f"FAMILY LEAKAGE: {len(overlap)} families in both train and test. "
        "South Wales (Dragon-3) is a subset of All-Wales PASS by FamilyNumber."
    )

def test_no_eid_leakage(train_cohort, test_cohort):
    overlap = set(train_cohort['eid']) & set(test_cohort['eid'])
    assert not overlap, f"EID LEAKAGE: {len(overlap)} eids in both splits"
```

```python
# tests/test_02_endpoint_provenance.py
WRONG_UKB_FIELDS = {'p131286', 'p131288', 'p131290', 'p131292', 'p131294'}  # code I10-I20 (HTN)
RIGHT_UKB_FIELDS = {'p131296', 'p131298', 'p131306'}                          # I20, I21, I25

def test_ascvd_endpoint_uses_correct_ukb_fields(endpoint_script_text):
    used_wrong = WRONG_UKB_FIELDS & set(endpoint_script_text.split())
    assert not used_wrong, (
        f"ENDPOINT BUG: {used_wrong} code I10-I20 (hypertension), NOT ASCVD. "
        "Use p131296 (I20), p131298 (I21), p131306 (I25) plus HES first-admission dates."
    )

EXPECTED_OPCS_SEVERE_AS = {'K611', 'K261', 'K262', 'K263', 'K264'}  # K611 = balloon valvuloplasty
EXPECTED_ICD_ASCVD = {'I20', 'I21', 'I22', 'I25', 'I63', 'G45', 'I70', 'I73', 'I74'}

def test_opcs_completeness_severe_as(used_opcs_codes):
    missing = EXPECTED_OPCS_SEVERE_AS - set(used_opcs_codes)
    assert not missing, f"MISSING OPCS for severe-AS: {missing} (K611 covers 3,088 balloon valvuloplasty cases)"
```

```python
# tests/test_03_cox_model_invariants.py
def test_noageldl_variant_present(cox_results):
    variants = set(cox_results['model'].unique())
    has_noageldl = any('NoAgeLDL' in v or 'no_age_ldl' in v.lower() for v in variants)
    assert has_noageldl, "MISSING SENSITIVITY: NoAgeLDL variant required (does the predictor carry signal beyond age+LDL?)"

def test_ph_assumption_per_predictor(cox_ph_check):
    failing = cox_ph_check[cox_ph_check['p_schoenfeld'] < 0.05]
    assert len(failing) == 0, f"PH violated for: {list(failing['predictor'])}. Use stratified Cox or time-varying coef."

def test_hr_plausibility(cox_results):
    out_of_range = cox_results[(cox_results['HR'] < 0.1) | (cox_results['HR'] > 10)]
    assert len(out_of_range) == 0, f"Implausible HRs (likely quasi-separation): {out_of_range[['model','HR']].to_dict()}"
```

```python
# tests/test_04_calibration_invariants.py
TRIPOD_REQUIRED = {
    'C_statistic', 'C_statistic_CI_low', 'C_statistic_CI_high',
    'cal_intercept', 'cal_slope',
    'brier', 'brier_scaled',
    'NRI', 'NRI_CI_low', 'NRI_CI_high',
    'IDI',
    'DCA_NB_5pct', 'DCA_NB_10pct', 'DCA_NB_20pct',
    'DeLong_p_vs_comparator',
}

def test_tripod_metrics_present(calibration_results):
    missing = TRIPOD_REQUIRED - set(calibration_results.columns)
    assert not missing, f"TRIPOD INCOMPLETE: missing {missing}"

def test_calibration_slope_plausible(calibration_results):
    for _, row in calibration_results.iterrows():
        assert 0.5 <= row['cal_slope'] <= 1.5, (
            f"Calibration slope {row['cal_slope']:.2f} for model {row['model']} "
            "suggests miscalibration (recommend isotonic recalibration before reporting)."
        )

def test_nri_bootstrap_size(calibration_results):
    for _, row in calibration_results.iterrows():
        assert row['NRI_B'] >= 100, f"NRI bootstrap B={row['NRI_B']} < 100; refit BOTH baseline and augmented Cox per resample"
```

```python
# tests/test_06_manuscript_claims.py
def test_every_manuscript_claim_traceable(provenance_table):
    untraceable = provenance_table[provenance_table['source_csv'].isna()]
    assert len(untraceable) == 0, (
        f"{len(untraceable)} manuscript claims have no source CSV. "
        f"First: {untraceable.iloc[0]['claim_id']}"
    )

def test_every_claim_matches_locked_value(provenance_table, tol=1e-3):
    failures = []
    for _, row in provenance_table.iterrows():
        ms_val = float(row['manuscript_value'])
        locked_csv = pd.read_csv(row['source_csv'])
        csv_val = float(locked_csv.query(row['row_filter']).iloc[0][row['value_col']])
        if abs(ms_val - csv_val) > tol:
            failures.append((row['claim_id'], ms_val, csv_val))
    assert not failures, f"Manuscript drift: {failures}"
```

### Step 4 — Configure pytest for autonomous iteration

```ini
# pytest.ini
[pytest]
addopts = -v --tb=short -x --strict-markers
testpaths = tests
log_cli = true
log_cli_level = INFO
```

`-x` makes pytest exit on first failure so the autonomous loop can fix one bug at a time.

### Step 5 — Autonomous fix-loop pattern

When commissioned to "iterate until all tests pass":

```python
# scripts/iterate_until_green.py
import subprocess, sys

while True:
    result = subprocess.run(
        ['pytest', '-x'],
        encoding='utf-8', errors='replace',
        env={**os.environ, 'PYTHONIOENCODING': 'utf-8'},
        capture_output=True,
    )
    if result.returncode == 0:
        print("ALL TESTS PASS - reproducibility certified")
        sys.exit(0)
    # parse the failure, propose a fix, apply, loop
    failure_line = next(l for l in result.stdout.splitlines() if 'FAILED' in l)
    # ... (Claude diagnoses, edits, commits, re-runs)
```

The autonomous loop is the killer feature: encode each methodological invariant once and Claude enforces it forever without you having to flag the gap.

## Non-negotiables (encode by default for every cardiology paper)

| Invariant | Test file | Why |
|---|---|---|
| Family-level deduplication | `test_01_cohort_invariants.py` | Dragon-3 ⊂ Wales PASS; eid-only dedup misses it |
| Leak-free train/test splits | `test_01_cohort_invariants.py` | Catches accidental data reuse |
| NoAgeLDL sensitivity variant | `test_03_cox_model_invariants.py` | Tests novelty beyond age+LDL backbone |
| OPCS K611 in severe-AS code list | `test_02_endpoint_provenance.py` | 3,088 balloon valvuloplasty cases otherwise lost |
| Correct UKB ASCVD fields (p131296+) | `test_02_endpoint_provenance.py` | p131286-p131294 code hypertension, not ASCVD |
| Cox PH assumption | `test_03_cox_model_invariants.py` | If violated, stratify or use time-varying |
| TRIPOD metric completeness | `test_04_calibration_invariants.py` | C-stat+CI, cal_int+slope, Brier, NRI(B≥100), IDI, DCA, DeLong |
| BH-FDR for >10 simultaneous tests | `test_03_cox_model_invariants.py` | Multiple-testing correction |
| E-value for every headline OR/HR | `test_04_calibration_invariants.py` | Threshold > 2.0 for "robust to plausible unmeasured confounding" |
| Every manuscript claim has source CSV | `test_06_manuscript_claims.py` | No untraceable numbers |

## Output

For every reproducibility-tests pass, deliver:

1. `tests/INVARIANT_AUDIT.md` — full list of proposed invariants (user reviews)
2. `tests/conftest.py` — shared fixtures
3. `tests/test_0[1-7]_*.py` — one file per invariant class
4. `pytest.ini` — config
5. `scripts/iterate_until_green.py` — autonomous fix-loop driver
6. A chat summary: *"K/K tests passing, M skipped, 0 failing — reproducibility certified."*

## Pre-flight environment check

```bash
python --version | grep -q "3.12" || echo "WARN: Python 3.12 required for pyarrow"
python -c "import pytest; import pandas; print('OK')" || pip install pytest pandas
```

Windows: set `PYTHONIOENCODING=utf-8` for any subprocess call to pytest.

## Reference exemplar

The Cardiff MD-by-Published-Works dissertation (May 2026) shipped with:
- `tests/test_*.py` — 52 invariants across 5 papers
- `THESIS_TRACEABILITY_REPRODUCER.py` — 56/56 PASS one-shot
- `per_paper_reproducers/*` — 110/110 PASS thesis-wide

This skill is the "convert that to pytest" upgrade so the 110/110 becomes a CI gate, not a one-shot certification.

## Common pitfalls

- **Tests pass locally but fail on a clean run** — your conftest fixtures load from results CSVs that don't exist on a fresh clone. Add `scripts/regenerate_all_results.py` as a precondition.
- **Test passes for the wrong reason** — assertion checks the presence of a column but not its value. Always assert both presence AND a plausibility range.
- **Test runtime > 5 minutes** — split into `pytest -m fast` vs `pytest -m slow` markers; only `fast` runs on every save.
- **`-x` exits too aggressively in development** — flip to `--maxfail=5` while iterating, back to `-x` for CI.
