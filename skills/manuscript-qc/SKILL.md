---
name: manuscript-qc
description: Use when validating every numerical claim in a cardiology / FH / UKB manuscript against source data CSVs. Triggers on phrases like "QC the manuscript", "audit numbers", "are all results traceable", "check every number", "submission-ready audit", or before any manuscript submission. Produces an automated reproducer script and a per-paper provenance table; goal is a 100% PASS / 0 FAIL audit.
---

# Manuscript-QC — every number maps to source

## When to invoke

Use this skill before any manuscript submission, response-letter assembly, or thesis lock. Specifically when the user says any of:
- "QC the manuscript" / "audit the numbers"
- "are all results traceable to raw [data source]"
- "check every number"
- "submission-ready audit"
- "build a reproducer"
- "lock the analysis"
- "we had a bug — re-verify everything"

## Five-step protocol

### Step 1 — Parse the manuscript and extract every numerical claim

Read the target manuscript (.md, .docx, or .txt) and build a structured ledger of every numerical claim. For each claim capture:
- `claim_id` — short identifier (e.g. `COX_C1_HR`, `META_FH_pooled`)
- `manuscript_section` — section number / chapter
- `headline_text` — the prose around the number
- `manuscript_value` — the value as it appears in the manuscript
- `category` — cohort N, hazard ratio, odds ratio, percentage, p-value, calibration slope, etc.

For a 7,000-word cardiology manuscript expect 50–100 distinct numerical claims. Save as `PROVENANCE_TABLE.csv`.

### Step 2 — Map each claim to its source CSV row

For every claim, locate the source CSV in `/results/`, `/nejm-email/result_csvs/`, or the analysis output directory. Add columns to the provenance table:
- `source_csv` — filename
- `row_filter` — e.g. `analysis=='C1_discord_only_ASCVD'`
- `n` — analytic n
- `events` — event count
- `date_locked` — when the CSV was last regenerated

If a claim has NO matching CSV, flag it as a critical issue — no number should appear in a manuscript without a source row.

### Step 3 — Build an automated reproducer

Write a Python script `<paper>_TRACEABILITY_REPRODUCER.py` that:
1. Loads each source CSV
2. Extracts the locked value for every claim
3. Asserts `manuscript_value == locked_value` within numerical tolerance (default 1e-3)
4. Prints PASS / FAIL / SKIP per claim
5. Writes a markdown audit report `<paper>_audit.md`
6. Returns exit code 0 if all PASS, non-zero otherwise

Use the helper pattern:

```python
def check(claim, ms, csv, source, note='', tol=1e-3):
    if pd.isna(ms) or pd.isna(csv):
        status = 'SKIP'
    else:
        diff = abs(float(ms) - float(csv))
        status = 'PASS' if diff <= tol else f'FAIL({diff:.4f})'
    AUDIT.append({'claim':claim,'ms':ms,'csv':csv,'source':source,'status':status})
    print(f'  [{status:10s}] {claim:24s}  ms={ms!s:14s}  csv={csv!s:18s}')
```

Estimated reproducer length: ~600 lines for a single paper, ~1200 lines for a 5-paper thesis.

### Step 4 — Run the reproducer and report

Execute the reproducer with `python <paper>_TRACEABILITY_REPRODUCER.py` and confirm:
- 100% PASS
- 0 FAIL
- 0 SKIP (any SKIP indicates missing source — must fix)

If any FAIL: present the failing rows to the user with both manuscript and CSV values, and ask whether to (a) update the manuscript to match the locked value, (b) re-run the analysis to refresh the CSV, or (c) investigate a deeper bug.

### Step 5 — Add reproducer entry-point to CLAUDE.md / docs

Reference the reproducer in the manuscript's reproducibility appendix:

```
## Reproducibility command
cd <paper-directory>
python <paper>_TRACEABILITY_REPRODUCER.py
Expected: N/N PASS / 0 FAIL. Runtime ~60 seconds.
```

## Non-negotiables (apply by default)

- **Family-level deduplication** for any FH cohort sharing registries (Dragon-3 ⊂ Wales PASS by FamilyNumber).
- **NoAgeLDL sensitivity** variant for any discordance / FH risk model.
- **OPCS-4 / ICD-10 completeness check** — include K611 (balloon valvuloplasty) for severe-AS; full I20-I25/I63/G45/I70/I73 for ASCVD composite.
- **Endpoint provenance** — UKB p131286-p131296 codes I10-I20 (hypertension), NOT ASCVD. Use p131296 (I20), p131298 (I21), p131306 (I25) plus HES first dates.
- **TRIPOD-compliant reporting** — C-statistic + 95% CI, calibration intercept + slope, Brier, NRI bilateral bootstrap (B≥100), IDI, DCA at 5/10/20%, DeLong p.
- **Multiple-testing correction** — Benjamini–Hochberg FDR for any panel of >10 simultaneous tests.
- **E-value sensitivity** — VanderWeele & Ding 2017; threshold > 2.0 for "robust to plausible unmeasured confounding".

## Output deliverables

For every QC pass, produce:
1. `PROVENANCE_TABLE.csv` (one row per claim → file:row:n:events:date)
2. `<paper>_TRACEABILITY_REPRODUCER.py` (auto-runnable audit)
3. `<paper>_audit.md` (auto-generated PASS/FAIL table)
4. A short chat summary: "X / X PASS, 0 FAIL — submission-ready"

## Reference exemplar

The Cardiff MD-by-Published-Works dissertation (May 2026) achieved:
- Single-paper: `PAPER1_TRACEABILITY_REPRODUCER.py` — 66 / 66 PASS
- Five-paper integrated: `THESIS_TRACEABILITY_REPRODUCER.py` — 56 / 56 PASS
- Per-paper suite: `run_all_reproducers.py` — 110 / 110 PASS thesis-wide

Files at: `D:/Projects/Lpa_Multilevel/nejm-email/`.

## Common failures and fixes

| Symptom | Likely cause | Fix |
|---|---|---|
| Reproducer reports CSV value differs from manuscript | Manuscript was written from a buggy / stale analysis | Lock the corrected analysis with an end-to-end rerun script (see `LOCK_RERUN_composite_ASCVD.py` exemplar); rewrite manuscript text. |
| `SKIP` for any claim | No matching source CSV found | Either find the missing CSV or remove the claim from the manuscript. |
| `python-docx` reports 0 substitutions when find-replacing in .docx | Word file is locked OR text is fragmented across runs | Close Word; iterate per-run; verify substitution count > 0. |
| Subprocess call to reproducer crashes with cp1252 decode error | Reproducer prints Unicode | Pass `encoding='utf-8', errors='replace', env={'PYTHONIOENCODING':'utf-8'}` to subprocess. |
