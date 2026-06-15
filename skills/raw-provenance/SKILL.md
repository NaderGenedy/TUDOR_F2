---
name: raw-provenance
description: Use when verifying that manuscript numbers trace ALL THE WAY BACK to the original raw source extracts (raw UKB batch CSVs, raw PASS/Dragon-3/All-Wales registry files), not merely to the derived analysis CSVs. Triggers on "trace back to raw", "are these from the original UKB CSVs", "full provenance", "two-level provenance", "recompute from raw", "level-2 verification", or whenever a QC pass has certified manuscript-vs-derived-CSV but not derived-CSV-vs-raw. Closes the gap where a bug in the cohort-build or analysis script would pass a Level-1 QC undetected.
---

# Raw-provenance — trace every number to the original source, not just the derived CSV

## The core idea: provenance has TWO links, most QC checks only ONE

```
RAW SOURCE  ──[build + analysis scripts]──>  DERIVED CSV  ──[hand-transfer]──>  MANUSCRIPT
  (UKB batch1/2/3,                              (part*.csv,                       (prose, tables)
   PASS / Dragon-3,                              cox_results.csv,
   All-Wales registry)                          mediation.csv)

            └────────── LINK 2 ──────────┘    └────────── LINK 1 ──────────┘
            (raw-provenance, THIS skill)       (manuscript-qc / final QC)
```

**`manuscript-qc` and most "32/32 certified" passes verify LINK 1 only** — manuscript prose equals the derived CSV value. They do NOT re-derive the CSV from raw. So a bug in `01_build_cohort.py` (wrong filter, wrong field mapping, wrong menopause classification) or in `03_analysis.py` (wrong groupby, wrong median) produces a wrong derived CSV that the manuscript faithfully copies — and LINK-1 QC passes green while the science is wrong.

**This skill verifies LINK 2**: recompute the derived value directly from the raw extract and assert they match. Only when BOTH links pass is a number genuinely traceable to source.

## When to invoke

- "are all the numbers traced back to the original UKB / Wales CSVs?"
- "full provenance" / "two-level provenance" / "level-2 verification"
- "recompute from raw" / "rebuild the cohort from the batches and check"
- After a `manuscript-qc` pass, when the user asks whether the derived CSVs themselves are trustworthy
- Before any high-stakes submission (NEJM/Lancet/Circulation) where a reviewer may request the raw-to-result audit trail

## Five-step protocol

### Step 1 — Map the provenance chain for each headline number

For every headline claim, identify BOTH links:

| Claim | Raw source | Build/analysis script | Derived CSV | Manuscript |
|---|---|---|---|---|
| Cohort n | `batch1_menopause.csv` | `01_build_cohort.py` | `cohort.csv` | Methods |
| Lipid median | `batch2_lipids.csv` | `03_analysis.py` | `part2_*.csv` | Results |
| Cox HR | cohort + ICD/OPCS | `17_cox.py` | `cox_results.csv` | Results |
| External validation | `matched_fsh_lipid_profile.csv` (PASS) | `02_pass.py` | `pass_*.csv` | Results |

Write this as `provenance_chain.csv` — one row per claim, columns: `claim, raw_file, build_script, derived_csv, manuscript_section`.

### Step 2 — Re-read the build script to extract the EXACT filter + transform logic

Do not guess the cohort filter. Open `01_build_cohort.py` (or equivalent) and transcribe:
- The join logic (inner vs left, on which key)
- The inclusion filter (e.g. `df[df['p31'] == 0]` for females)
- The classification logic (e.g. the `p2724_i0 -> i1 -> i2` fallback for menopause status)
- Any derived-variable formulas (non-HDL = TC - HDL, eGFR CKD-EPI, etc.)

Replicate this logic EXACTLY in the verification script. A verification that uses a different filter than the build is not a verification.

### Step 3 — Recompute from raw, memory-safe

Raw UKB batches are ~500k rows. Always use `usecols` to load only the needed columns and strip the `participant.` prefix:

```python
import pandas as pd

# Only the columns we need (sex + menopause + lipids)
b1 = pd.read_csv('batch1_menopause.csv',
                 usecols=['participant.eid', 'participant.p31',
                          'participant.p2724_i0', 'participant.p2724_i1', 'participant.p2724_i2'],
                 low_memory=False)
b1.columns = [c.replace('participant.', '') for c in b1.columns]

b2 = pd.read_csv('batch2_lipids_biochemistry.csv',
                 usecols=['participant.eid', 'participant.p30690_i0', 'participant.p30780_i0'],
                 low_memory=False)
b2.columns = ['eid', 'total_cholesterol', 'ldl_c']

# Replicate the build's INNER join + female filter EXACTLY
raw = b1.merge(b2, on='eid', how='inner')
raw = raw[raw['p31'] == 0].copy()

# Replicate the menopause classification fallback chain EXACTLY
raw['meno_raw'] = raw['p2724_i0']
for instance in ['p2724_i1', 'p2724_i2']:
    mask = raw['meno_raw'].isna() | (raw['meno_raw'] == -3)
    raw.loc[mask, 'meno_raw'] = raw.loc[mask, instance]
```

For files too large even with usecols, chunk and aggregate incrementally.

### Step 4 — Assert recomputed == claimed, with source attribution

```python
def add(test, expected, actual, tol, raw_source):
    ok = abs(float(expected) - float(actual)) <= tol
    checks.append((test, expected, actual, raw_source, 'PASS' if ok else 'FAIL'))

pre = raw[raw['status'] == 'Premenopausal']
add("RAW->TC premenopausal median", 5.32, pre['total_cholesterol'].median(),
    tol=0.05, raw_source="batch2 lipid assay")
```

Use realistic tolerances: medians rounded to 2 dp tolerate 0.05; percentages tolerate 1.0; counts tolerate a small N for QC-merge differences (but cohort sizes should match within ~50–200).

### Step 5 — Produce the combined two-level certificate

Merge the LINK-1 (manuscript-qc) and LINK-2 (this skill) results into one certificate:

```
COMBINED PROVENANCE CERTIFICATE
Level 1 (manuscript <-> derived CSV):  K1/K1 PASS   [final_qc_verification.py]
Level 2 (derived <-> recompute-raw):   K2/K2 PASS   [raw_provenance_verification.py]
---
TOTAL: (K1+K2)/(K1+K2) traced to original source
```

Only claim "fully traceable to raw source" when BOTH levels pass.

## Output deliverables

1. `provenance_chain.csv` — the two-link map per claim
2. `raw_provenance_verification.py` — re-runnable Level-2 script
3. `RAW_PROVENANCE_REPORT.md` — Level-2 PASS/FAIL table with raw-source attribution
4. `COMBINED_PROVENANCE_CERTIFICATE.md` — Level-1 + Level-2 unified
5. A chat summary: *"Level 1: K1/K1, Level 2: K2/K2 — every headline number recomputes from the original raw extracts."*

## Cohort-specific cross-checks (apply when relevant)

| Cohort | Raw file | Filter to replicate exactly | Gotcha |
|---|---|---|---|
| UKB females | `batch1_menopause.csv` | `p31 == 0` AND inner-join with lipid batch | menopause status uses `p2724_i0 -> i1 -> i2` fallback, treat `-3` as missing |
| PASS perimenopause | `matched_fsh_lipid_profile.csv` | FSH-confirmed status column | distinct from the FH-project PASS — no family structure here |
| FH Dragon-3 / All-Wales | registry CSVs | **family-level dedup by FamilyNumber** | Dragon-3 (SW Wales) IS a subset of All-Wales PASS — dedup by `FamilyNumber`, not just `DatabaseNumber` (see `cox-analysis` skill) |

**Do not conflate the two "PASS" cohorts:** the perimenopause PASS (330 FSH-confirmed women) has no FH family structure and needs no FamilyNumber dedup; the FH PASS (All-Wales registry) does.

## Reference exemplar

**Perimenopause manuscript (May 2026):** the original 32/32 certification verified LINK 1 only (manuscript == `data/part*.csv`). When the PI asked "are all traced back to the original UKB and Wales CSVs?", a Level-2 pass was built (`submission_v2/raw_provenance_verification.py`):

- Reconstructed the cohort straight from `batch1_menopause.csv` (501,936 rows) → filter `p31==0` → classify `p2724` → **exactly 273,036 / 63,939 / 165,276** (matched the build to the row).
- Recomputed every primary lipid median from raw `batch2` assay columns → all matched to 3 dp (TC 5.32→6.02 +13.1%, LDL 3.23→3.71 +14.9%, ApoB, TG, Lp(a)).
- Recomputed PASS TC +9.8% from raw `matched_fsh_lipid_profile.csv`.
- Result: **Level 2 = 19/19 PASS**, combined with Level 1 = **51/51 fully traceable to raw source.**

## Pre-flight environment check

```bash
python --version | grep -q "3.12" || echo "WARN: Python 3.12 required (pyarrow)"
PYTHONIOENCODING=utf-8 python -c "import pandas; print('OK')"
```

Large raw extracts: always `usecols` + `low_memory=False`; chunk if a single needed column set still OOMs. Never read a 500k x 100 raw extract with bare `pd.read_csv`.

## Common pitfalls

- **Verification uses a different filter than the build** — then a green pass is meaningless. Always transcribe the build's exact filter, never re-invent it.
- **Tolerance too tight on counts** — QC merges can drop a handful of rows; cohort sizes should match within ~50–200, not exactly, unless you replicate every join.
- **Forgetting the `participant.` prefix** — raw UKB RAP extracts prefix every column; strip it before merging.
- **Treating derived CSVs as ground truth** — they are the thing under test in Level 2, not the reference. The reference is the RAW extract.
- **Conflating same-named cohorts across projects** — the perimenopause PASS and the FH All-Wales PASS are different datasets with different dedup requirements.
