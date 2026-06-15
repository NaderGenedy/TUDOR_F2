---
name: tudor-reproduce
description: >
  Use to download, restart, or fully reproduce the TUDOR FH diagnostic project end-to-end —
  from a clean machine or Claude Code on the web. Triggers on "restart TUDOR", "reproduce
  TUDOR from scratch", "set up the TUDOR repo", "run the whole TUDOR pipeline", "rebuild
  TUDOR on the web", or "get TUDOR running on a new machine". Sets up Python 3.12, generates
  web-safe synthetic data, runs the full method (stages 0-8), prints the locked headline
  numbers, and enforces the data-governance firewall (UK Biobank / All-Wales data must never
  be uploaded). For the REAL numbers it points to the governed-data path; everywhere else it
  reproduces the method on synthetic data.
metadata:
  type: reference
  project: TUDOR_F
  owner: Dr Nader Genedy (UKB App 1002450)
---

# tudor-reproduce — restart / reproduce the TUDOR project

## When to use
- Cloning TUDOR_F onto a new machine, or opening it in Claude Code on the web.
- Re-running the pipeline after an environment reset.
- Verifying the headline numbers (head-to-head, NRI, 25/25 subgroups) reproduce.

## Hard rule (data governance — non-negotiable)
UK Biobank (App 1002450) and All-Wales registry data are controlled-access. **Never** upload,
paste, or commit individual-level data anywhere outside the secure environment — not even to a
private GitHub repo or a web sandbox (breaches the UKB MTA / NHS IG). The repo's `.gitignore`
is default-deny on all data formats; keep it that way. Real numbers regenerate ONLY where the
governed data physically lives (local workstation / UKB RAP).

## Steps

### 1. Environment (Python 3.12 — not 3.14)
```bash
python -m pip install -r code/requirements.txt
```

### 2. Reproduce the METHOD (anywhere, incl. web)
```bash
python code/make_synthetic_tudor_data.py        # synthetic Wales + UKB (NO real patients)
python code/TUDOR_MASTER.py --check             # inventory; governed data reads [ABSENT] off-site
python code/TUDOR_MASTER.py --synthetic --all   # stages 0-8 end-to-end
python code/TUDOR_MASTER.py --reproduce         # print locked headline targets
```
Expect: TUDOR > eDLCN > MEDPED > Simon Broome, TUDOR winning every subgroup (synthetic AUCs
run high by design — method demonstration, not real effect sizes).

### 3. Reproduce the REAL numbers (local / RAP only)
```bash
python code/TUDOR_MASTER.py --check             # all three governed files must read [PRESENT]
python code/TUDOR_MASTER.py --real --all
# then the full published analyses:
python code/CALON_rich_wales.py
python code/analysis/head_to_head_calon_vs_published_safeheart.py
python code/analysis/comprehensive_subgroup_pooled.py
python code/analysis/step32_cox_ascvd_nri.py
```
Cross-check every value against `results/verified_numbers_locked.json` and
`results/HEAD_TO_HEAD_RESULTS.md`.

### 4. Methodology invariants to confirm (baked into the code)
- **Family-level deduplication** (one event per FamilyNumber) — Stage 1.
- **NoAgeLDL sensitivity** variant alongside the primary model — Stage 4.
- **Frozen-coefficient** external comparators (TRIPOD Type 4) — head-to-head scripts.
- **Bidirectional** Wales geographic validation (train N→test S, then reverse) — Stage 5/6.

### 5. QC before any claim of "reproduced"
Invoke companion skills: `tudor-qc` (5-agent pipeline audit), `manuscript-qc` (trace every
number to source), `comparator-model-validation` (head-to-head discipline), `raw-provenance`
(numbers → raw extracts). Goal: 100% PASS, 0 FAIL.

## Outputs
`output/` — cohort, coefficients, head-to-head scored table, subgroup deltas.
Truth files: `results/verified_numbers_locked.json`, `results/HEAD_TO_HEAD_RESULTS.md`.
