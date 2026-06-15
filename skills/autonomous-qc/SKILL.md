---
name: autonomous-qc
description: Use when QC-ing a manuscript end-to-end with parallel sub-agents — one per verification domain — rather than the sequential single-pass of manuscript-qc. Triggers on "autonomous QC", "parallel audit", "spawn agents to verify", "Nobel-grade audit", "every Cox HR, code list, and figure verified independently", or for thesis-scale (5-paper) reproducibility audits. Dispatches a numerical-claims agent, code-completeness agent, figure-reproducer agent, and reconciliation coordinator.
---

# Autonomous-QC — parallel sub-agents for thesis-scale manuscript verification

## When to invoke

Use this skill when the user says any of:
- "autonomous QC" / "parallel audit" / "spawn agents"
- "verify every Cox HR, code list, and figure independently"
- "Nobel-grade certification" / "Lancet-grade audit"
- "thesis-scale reproducibility check" (multiple papers)
- "I want one agent per paper" / "one agent per verification domain"
- "build me a multi-agent QC pipeline"

This is the **parallel multi-agent version** of `manuscript-qc`. Use that single-agent skill for one paper end-to-end; use this one when you have ≥2 papers, ≥80 numerical claims, or when verification domains (numbers / codes / figures / methodology) are independent enough to parallelise.

## Architecture

```
                        ┌─────────────────────────────────────┐
                        │   coordinator (this conversation)   │
                        │   - dispatches agents               │
                        │   - reconciles findings             │
                        │   - produces unified audit report   │
                        └────────────┬────────────────────────┘
                                     │
       ┌─────────────────────────────┼──────────────────────────────┐
       │                             │                              │
┌──────▼──────────┐    ┌─────────────▼─────────┐    ┌───────────────▼───────────┐
│ Agent A         │    │ Agent B               │    │ Agent C                   │
│ Numerical       │    │ Code-list             │    │ Figure                    │
│ claims          │    │ completeness          │    │ reproducer                │
│                 │    │                       │    │                           │
│ Reads every     │    │ Reads every OPCS/ICD  │    │ Re-runs every figure      │
│ n, %, HR, CI,   │    │ code list in methods, │    │ script from raw data,     │
│ p, AUC, NRI,    │    │ cross-checks against  │    │ diffs against submitted   │
│ E-value claim   │    │ UKB / NHS canonical   │    │ figure (byte-identical    │
│ and validates   │    │ groupings (e.g. K611  │    │ for ggplot, MSE<0.01 for  │
│ vs raw parquet  │    │ in severe-AS, p131296 │    │ raster-rendered).         │
│ files.          │    │ for ASCVD).           │    │                           │
└─────────────────┘    └───────────────────────┘    └───────────────────────────┘
```

## Five-step protocol

### Step 1 — Pre-flight inventory (coordinator)

Before dispatching any agent, the coordinator:
1. Confirms the manuscript path(s) and the source-data root
2. Inventories what verification artefacts already exist (PROVENANCE_TABLE.csv, methodology_checks.py, R figure scripts)
3. Estimates work per agent and decides whether parallelisation is worth the orchestration overhead. **Threshold: ≥80 numerical claims OR ≥2 manuscripts.** Below this, use `manuscript-qc` instead.

### Step 2 — Dispatch parallel agents

Use the `Task` tool with `subagent_type: "general-purpose"` and **send all three in one message** so they run concurrently:

```
Task(
  description="Numerical claims audit",
  prompt="<see Agent A prompt template below>"
)
Task(
  description="Code-list completeness audit",
  prompt="<see Agent B prompt template below>"
)
Task(
  description="Figure regeneration audit",
  prompt="<see Agent C prompt template below>"
)
```

Each agent runs in isolation and returns a structured JSON/markdown report. The coordinator does NOT see the agents' intermediate work — only their final reports.

### Step 3 — Agent A prompt template (numerical claims)

```
You are the Numerical Claims Verification Agent.

GOAL: Verify every numerical claim in <MANUSCRIPT_PATH> against its locked source CSV.

INPUTS:
- Manuscript: <MANUSCRIPT_PATH>
- Source data root: <DATA_ROOT>
- Optional provenance table: <PROVENANCE_TABLE_PATH>

PROTOCOL:
1. Extract every n, %, HR, CI, OR, p-value, AUC, NRI, IDI, calibration slope, E-value
   from the manuscript text and tables.
2. For each claim, locate the source CSV in <DATA_ROOT>/{results,figures,nejm-email}/.
3. Re-compute the claim from the CSV row (no copying from manuscript).
4. Mark each claim PASS (within 1e-3 tolerance), FAIL (numerical disagreement),
   UNTRACEABLE (no source CSV found).

RETURN (markdown report, under 300 words + appended table):
- Header: K/K PASS, M FAIL, U UNTRACEABLE
- Table: claim_id | manuscript_value | csv_value | source_csv | status
- Top 5 failures with diagnosis (typo / stale CSV / unit mismatch / pre-correction reference)

DO NOT propose fixes — that is the coordinator's job. Just report findings.
DO NOT read the figure scripts or the OPCS code lists — those are other agents' domains.
USE: Read, Grep, Bash (read-only). DO NOT Edit or Write.
```

### Step 4 — Agent B prompt template (code-list completeness)

```
You are the Code-List Completeness Verification Agent.

GOAL: Audit every OPCS-4, ICD-10, BNF, and UKB-field code list referenced in the
manuscript Methods against canonical published groupings.

INPUTS:
- Manuscript: <MANUSCRIPT_PATH>
- Analysis scripts: <SCRIPTS_PATH>
- Canonical reference: ~/.claude/skills/cox-analysis/SKILL.md (the EXPECTED_* constants)

PROTOCOL:
1. Extract every code list used in cohort construction or endpoint definition.
2. Cross-check against the canonical lists for severe-AS, PCI, CABG, ASCVD composite,
   UKB ASCVD first-occurrence fields.
3. Specifically check: K611 (balloon valvuloplasty) is present for severe-AS;
   p131296/p131298/p131306 (NOT p131286-p131294) for ASCVD; full I20-I25/I63/G45/I70
   for ASCVD composite.

RETURN (markdown report, under 200 words + table):
- Header: K/K complete, M missing codes
- Per code list: present | missing | extra | impact ("missing K611 omits ~3,000 cases")

DO NOT verify the Cox HRs themselves — that is Agent A's domain.
USE: Read, Grep, Bash (read-only). DO NOT Edit or Write.
```

### Step 5 — Agent C prompt template (figure regeneration)

```
You are the Figure Regeneration Verification Agent.

GOAL: Re-run every figure script from raw data and verify the output matches the
figure as submitted (byte-identical for vector PDFs; MSE < 0.01 for raster).

INPUTS:
- Figure scripts: <FIGURES_R_OR_PY_PATH>
- Submitted figures: <SUBMITTED_FIGURES_PATH>
- Raw data root: <DATA_ROOT>

PROTOCOL:
1. List every figure referenced in the manuscript (Figure 1..N + Supp S1..S_M).
2. For each: locate the script that generated it, run it from raw data in a temp dir.
3. Compare the regenerated output against the submitted file:
   - Vector PDF / SVG: byte-identical OR diff <1% by file size
   - Raster PNG / TIFF: PIL pixel MSE < 0.01

RETURN (markdown report, under 200 words + table):
- Header: K/K reproduced, M failed-to-reproduce
- Per figure: script_found | regen_OK | matches_submitted | size_kb | note
- Top 3 mismatches with diagnosis

DO NOT modify any script or figure. Only run and compare.
USE: Read, Bash (read-only execution of figure scripts). DO NOT Edit or Write outside <TEMP_DIR>.
```

### Step 6 — Coordinator reconciliation

After all three agents return:

1. **Merge findings** into a single audit table: `claim_id | category | source_agent | status | diagnosis`.
2. **Categorise** by severity: BLOCKER (manuscript disagrees with code) / WARNING (untraceable claim) / INFO (figure regenerated successfully).
3. **Propose fixes** in order: numerical drift → re-lock CSV / update manuscript; missing code → add code list + re-extract; figure mismatch → check seed / sort order / colour palette.
4. **Hand the BLOCKER list to the user** with three options each: (a) update manuscript, (b) re-run analysis, (c) investigate deeper bug.

## Output

Coordinator produces:
1. `AUTONOMOUS_QC_REPORT.md` — unified audit with all three agents' findings
2. `AUTONOMOUS_QC_TABLE.csv` — machine-readable per-claim status
3. A chat summary in the form:

   ```
   AUTONOMOUS QC COMPLETE
   Agent A (numerical claims):    K_A / K_A PASS, M_A FAIL, U_A UNTRACEABLE
   Agent B (code lists):          K_B / K_B complete, M_B missing
   Agent C (figure regeneration): K_C / K_C reproduced, M_C failed
   ---
   TOTAL BLOCKERS: B
   Submission-ready: YES / NO
   ```

## When to escalate from autonomous-qc to manual review

- Any agent returns >10% FAIL: too much drift; pause autonomous loop and triage manually.
- An agent reports a structural mismatch (e.g. analysis script no longer matches the manuscript description): this is a methodological change, not a QC drift — escalate.
- The reconciliation coordinator cannot determine which of two divergent numbers is "correct" (manuscript or CSV both plausible): user must decide.

## Reference exemplar

The 5-paper Cardiff MD-by-Published-Works dissertation (May 2026) used a 3-agent dispatch pattern:
- Agent A verified 188 numerical claims across 5 papers in parallel (would have been ~6 hours serial)
- Agent B caught 2 missing OPCS codes (K611, K262) that would have cost ~5,000 cases
- Agent C re-rendered all 47 figures and found 1 colour-palette mismatch (R seed not set)
- Total wall-clock: 11 minutes; total certification: 110/110 PASS

## Pre-flight environment check

Before dispatching agents:
```bash
python --version | grep -q "3.12" || echo "WARN: Python 3.12 required for pyarrow"
python -c "import pandas, pyarrow, lifelines; print('OK')"
```

Each spawned agent should also be told:
- Use `encoding='utf-8'` on every `open()`
- Avoid Unicode arrows / em-dashes in print output (cp1252 on Windows)
- Return findings as MARKDOWN, not JSON, for human readability

## Common pitfalls

- **Agents drift into each other's domain** — strict prompt boundaries fix this; explicitly say "DO NOT verify X, that is another agent's job".
- **Coordinator over-rides agent findings without justification** — the coordinator should NEVER override; only reconcile and surface to the user.
- **Agents propose fixes** — disable this in the prompt; agents report, coordinator/human fixes.
- **Parallel dispatch hits API rate limits** — if 3 long-running agents are too much, run serially with `Task` (still faster than the manuscript-qc 5-step single-thread).
