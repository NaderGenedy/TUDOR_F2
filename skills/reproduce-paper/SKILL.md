---
name: reproduce-paper
description: Use to build OR run a self-contained reproducibility test for a specific paper. Triggers on phrases like "build a reproducer for [paper]", "reproduce the manuscript", "lock the analysis", "make the X/X reproducer", "regression-test the paper", "block submission until reproduced". Distinct from manuscript-qc — that audits a finished manuscript against existing results; this builds a frozen ground-truth test suite that re-runs the analysis end-to-end and asserts every published number against fresh computation.
---

# Reproduce-paper — frozen ground-truth test suite per paper

## When to invoke

Use the moment the user says any of:
- "build a reproducer for [paper / manuscript]"
- "reproduce the manuscript"
- "lock the analysis"
- "make the X/X reproducer"
- "regression-test the paper"
- "block submission until 100% reproduced"
- "every number must come from a fresh computation"

This skill is **distinct from `manuscript-qc`**:
- `manuscript-qc` checks the manuscript text against existing result CSVs (assumes the CSVs are correct)
- `reproduce-paper` rebuilds the result CSVs from the slim cohort + locked engine, then asserts manuscript values against the fresh computation (catches stale CSVs)

Both run together for full submission-grade rigour.

## Six-step protocol

### Step 1 — Establish the paper's analytical contract

For the target paper, list:
- **Inputs**: every CSV the analysis depends on (slim, locked engine, raw extracts) with absolute path
- **Code**: every script that produces a published number, in dependency order
- **Outputs**: every result CSV that feeds a manuscript table or figure
- **Frozen claims**: the headline numbers (HR, AUC, NRI, IDI, n, events, p) the manuscript prints

Write the contract as `<paper>_CONTRACT.yml`. Every claim is named (e.g. `MVMR_AS_Lpa_OR`, `LDL_correction_pct_125`, `bed_Q_p_het`) with expected value + tolerance.

### Step 2 — Wrap each analysis pipeline as a pytest-style harness

Create `<paper>_REPRODUCER.py` with one function per claim:

```python
def test_MVMR_AS_Lpa_OR(slim, eng):
    """Real-cross-effects MVMR: Lp(a) on AS, conditional on LDL.  Expect 1.559 (1.480-1.642)."""
    # ... rerun fit_mvmr with frozen instrument matrix ...
    return {'OR': 1.559, 'lo': 1.480, 'hi': 1.642, 'p': 0.0}

def assert_claim(name, expected, fresh, tol):
    diff = abs(expected - fresh)
    status = 'PASS' if diff <= tol else f'FAIL (diff={diff:.4f})'
    AUDIT.append({'claim': name, 'expected': expected, 'fresh': fresh, 'status': status})
```

For each claim assert `manuscript_value == fresh_value` within tolerance. Default tolerances:
- HR / OR: ±0.005
- AUC: ±0.001
- NRI: ±0.0005
- IDI: ±0.00005
- n / events: exact match (0 tolerance)
- p-value: ratio within 10× (log-scale tolerance)

### Step 3 — Add proactive checks for known failure modes

Every reproducer in this stack must check for the historical bugs Dr Genedy has paid for:
- Family-level deduplication applied (FH cohorts share registries)
- NoAgeLDL sensitivity variant present
- OPCS-4 includes K611 (3,088 cases) when AS endpoint involved
- ICD-10 ASCVD composite uses p131296 / p131298 / p131306 — not p131286 (hypertension)
- Slim columns `lvm`, `mean_myo_t1` were renamed (they are LV stroke volume + Cardiac index, not LV mass + T1)
- Locked endpoint engine `t_<E>_years` per-endpoint, not coronary-FOC composite re-used

Hard-code these as boolean assertions at the top of every reproducer.

### Step 4 — Run the reproducer in clean state

Execute as a stand-alone script:
```bash
python <paper>_REPRODUCER.py
```

Expected return:
- exit code 0 if all PASS
- exit code 1 if any FAIL

Stdout produces a per-claim PASS/FAIL/SKIP table and a final certification line:
```
=== REPRODUCER REPORT ===
<paper>: 67 / 67 PASS, 0 FAIL, 0 SKIP — submission-ready
```

### Step 5 — Wire the reproducer into the project's run-all orchestrator

Append a phase to `run_all_<project>.py`:
```python
PIPELINE.append((9, 'py', f'{paper}_REPRODUCER.py', f'Reproducer: {paper} ({n}/{n} PASS expected)'))
```

So the reproducer is a real phase that runs after all the science scripts. If anything in the upstream pipeline drifts, the reproducer catches it.

### Step 6 — Add a self-healing diagnostic agent (optional — Claude Code subagent)

For papers where re-runs are expensive (>10 min), wrap the reproducer in a "doctor" loop:

```
On any FAIL:
  spawn Agent(subagent_type='general-purpose',
              prompt=f"Reproducer claim '{name}' fails: expected {exp}, fresh {got}. "
                     f"Diagnose root cause: env drift, data update, code bug, methodology gap. "
                     f"Propose ONE fix and the patch. Don't apply — report only.")
```

User reviews the proposed fix before applying. This keeps the loop human-in-the-loop while saving the user the diagnostic time.

## Reference exemplars

The Cardiff MD-by-Published-Works thesis demonstrated this pattern at scale:
- `PAPER1_TRACEABILITY_REPRODUCER.py` — 66 / 66 PASS
- `THESIS_TRACEABILITY_REPRODUCER.py` — 56 / 56 PASS (cross-paper integrated)
- `run_all_reproducers.py` — 110 / 110 PASS thesis-wide

Files at: `D:/Projects/Lpa_Multilevel/nejm-email/`.

For the v21 Lp(a) submission specifically, the contract should cover:
- bed-selectivity Cochran Q (54.6, p=1.77e-9)
- MVMR Lp(a) → AS conditional on LDL (OR 1.559, real cross-effects)
- Real Fine-Gray AS (sub-HR 1.138, agrees with cause-specific Cox 1.140)
- Lp(a) within-person stability (Pearson 0.952)
- LDL-correction reclassification (29.2% at Lp(a) ≥ 175)
- ASCVD risk score refinements (R1 NRI +0.0184, R2 apoB ≈ LDL, R3 Black HR 0.79)

## Output deliverables

Per paper:
1. `<paper>_CONTRACT.yml` — frozen ground-truth manifest
2. `<paper>_REPRODUCER.py` — runnable test suite
3. `<paper>_audit.md` — last-run PASS/FAIL table
4. Entry in `run_all_<project>.py`
5. Optional: a "doctor" agent dispatcher for self-healing

## Common failures and what they mean

| Symptom | Diagnosis |
|---|---|
| Expected n drops by ≥1% between runs | Slim cohort regenerated with new exclusion criteria — review |
| HR drift from 1.559 → 1.601 | Instrument matrix likely re-pulled (e.g. cross-effects updated); confirm intentional |
| Reproducer reports SKIP for a claim | Source CSV missing; either restore or remove the claim |
| Reproducer crashes with `UnicodeEncodeError` | Add `sys.stdout.reconfigure(encoding='utf-8')` at top |
| ALL claims FAIL after a re-run | Locked engine date may have advanced; check censor date and rerun upstream |
| 100% PASS but reviewer queries a number | Tolerance too loose; tighten and re-test |
