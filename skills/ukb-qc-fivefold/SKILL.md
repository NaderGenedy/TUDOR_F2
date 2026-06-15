---
name: ukb-qc-fivefold
description: Use whenever a UK Biobank-derived manuscript, statistical analysis, or risk-prediction model is being prepared for submission, when validating that numerical claims trace to raw data, when investigating suspected pipeline bugs or unreproducible results, when the user mentions "QC", "double-check", "validate", "fully traceable", "reproducer", "every number", or any FH / TUDOR / CALON / cardiology-cohort pipeline. Five independent specialised agents run in parallel — Raw Data Integrity, Cohort Validator, Feature Engineering Auditor, Statistical Reproducer, Provenance Detector — each producing a structured PASS/DRIFT/FAIL report. Catches manuscript-vs-live drift, hard-coded literals masquerading as live values, cohort labelling errors, missing ICD-10/OPCS codes, family-level deduplication failures, treatment-adjustment bugs, calibration miscomputation, and NRI/IDI methodology errors. Designed for the specific failure mode in which a number is correctly computed at one point in pipeline history, hand-transferred into prose as a literal, and then the underlying pipeline evolves without updating the prose. Trigger this BEFORE submitting any UKB analysis externally, not after.
---

# UKB QC Fivefold — Five Independent Audit Agents

## Why this skill exists

UK Biobank analyses fail in characteristic, recurring ways:

1. **Cohort labelling drift** — "Wales" in the manuscript means All-Wales PASS, but the cached value was computed on SouthWales / CAVUHB. They are different cohorts.
2. **Hand-transferred literals** — an interactive R session computes NRI=0.358, the analyst types that into the manuscript prose, the pipeline is later re-run with a different LOCO setup producing NRI=−0.097, and the prose never gets updated. The manuscript ships with a stale value.
3. **Methodology mismatch** — calibration slope reported is 6.33 but the actual TRIPOD-standard `glm(y ~ logit_pred)` produces 1.23. The 6.33 came from some other regression that was never specified.
4. **Family-level dedup forgotten** — South Wales (Dragon-3) is a subset of All-Wales PASS by FamilyNumber. Removing only DatabaseNumber overlap leaves siblings double-counted.
5. **OPCS-4/ICD-10 code omissions** — K611 (balloon valvuloplasty) absent from a severe-AS code list loses ~3,000 cases. I35 (aortic stenosis) wrongly included as ASCVD inflates events.
6. **Treatment-adjusted LDL miscomputed** — uniform 1.43 correction factor applied instead of drug-specific reductions (15–55% range). Discriminative performance reported is wrong.
7. **NMR field cross-talk** — `p131286–p131296` code I10-I20 (hypertension), NOT ASCVD. Using these as ASCVD endpoints invalidates the entire outcome definition.

Any single one of these can break a manuscript at peer review or, worse, after publication. **One careful pair of eyes is not enough — five independent specialists, each looking only at their own domain, catch what humans and single-agent reviewers miss.**

## How the skill runs

Five agents launch in parallel via the `Agent` tool, each with a dedicated instruction file in `agents/`. Each agent receives:

- Path to the project root (where raw UKB CSVs, scripts, and manuscript live)
- Path to a "claims file" (optional)
- Tolerance settings (default: 1×10⁻³ for AUCs, 5×10⁻³ for proportions, 5% for counts)

Each agent writes its report to `qc_output/agent_<N>_<name>_report.json` plus a markdown summary. The orchestrator script `scripts/run_fivefold_qc.py` collects all five and produces `qc_output/MASTER_QC_REPORT.md`.

Agents do **not** share state. If two agents agree, the finding is robust; if they disagree, that itself is informative.

## When to launch (workflow)

```
1. User mentions QC / validate / "double-check" / "every number" / submission-ready / "any bug"
2. Read the user's project layout (where is the manuscript? raw CSVs? cached outputs?)
3. Create qc_output/ alongside the project
4. Launch all five agents in PARALLEL (one Agent tool call per agent in the same message)
5. When all five complete, run `python scripts/run_fivefold_qc.py qc_output/` to aggregate
6. Present MASTER_QC_REPORT.md to the user
7. For each FAIL or DRIFT, propose a concrete fix and wait for user direction
```

**Critical rule — never edit manuscript text in response to a DRIFT without explicit user confirmation.** Drift is a flag, not a fix.

## Agent dispatch

| # | Agent | Domain | Reads | Writes |
|---|---|---|---|---|
| 1 | Raw Data Integrity Auditor | UKB CSV files | raw `ukb_*.csv` | `agent_1_raw_data.{json,md}` |
| 2 | Cohort Validator | inclusion/exclusion logic | raw CSVs + cohort def | `agent_2_cohort.{json,md}` |
| 3 | Feature Engineering Auditor | engineered features | raw lipids/treatment + cached features | `agent_3_features.{json,md}` |
| 4 | Statistical Reproducer | every numerical claim | predictions CSV + manuscript | `agent_4_stats.{json,md}` |
| 5 | Provenance Detector | hard-coded vs traced numbers | manuscript-generator script + all CSVs + .docx | `agent_5_provenance.{json,md}` |

Read each agent's spec from `agents/agent_<N>_*.md` before spawning that agent.

## Tolerance defaults

```python
TOLERANCE = {
    'AUC': 1e-3,
    'proportion': 5e-3,
    'count': 0.05,
    'p_value_log': 0.5,
    'OR_HR': 0.05,
    'CI_width': 0.1,
}
```

## Spawning the five agents — example invocation

```python
# In one message, five parallel Agent tool calls:
Agent(description="Agent 1: Raw UKB data integrity", subagent_type="general-purpose",
      prompt=read_file('agents/agent_1_raw_data_integrity.md') + project_context)
Agent(description="Agent 2: Cohort validation", subagent_type="general-purpose",
      prompt=read_file('agents/agent_2_cohort_validation.md') + project_context)
Agent(description="Agent 3: Feature engineering", subagent_type="general-purpose",
      prompt=read_file('agents/agent_3_feature_engineering.md') + project_context)
Agent(description="Agent 4: Statistical reproducer", subagent_type="general-purpose",
      prompt=read_file('agents/agent_4_statistical_reproducer.md') + project_context)
Agent(description="Agent 5: Provenance detector", subagent_type="general-purpose",
      prompt=read_file('agents/agent_5_provenance_detector.md') + project_context)
```

## Aggregation and presentation

```bash
python scripts/run_fivefold_qc.py <project_root> --out qc_output/
```

Produces `MASTER_QC_REPORT.md` with traffic-light summary per agent + critical findings + drifts requiring adjudication.

## What this skill does NOT do

- It does not edit the manuscript. It reports.
- It does not silently "fix" numbers.
- It does not replace TRIPOD or STROBE — those are additional checklists.

## References (bundled)

- `references/ukb_field_reference.md` — UKB field map
- `references/icd10_opcs4_codes.md` — cardiology code lists, with K611/I35/p131286 traps flagged
- `references/tripod_checklist.md` — TRIPOD reporting standard
- `references/nri_idi_methodology.md` — exact R-equivalent NRI/IDI Python implementation

## Helper scripts (bundled)

- `scripts/run_fivefold_qc.py` — orchestrator
- `scripts/helpers/ukb_field_audit.py` — Agent 1 worker
- `scripts/helpers/icd10_completeness.py` — code-list completeness
- `scripts/helpers/nri_idi_reference.py` — NRI/IDI Python verified bit-identical to R
- `scripts/helpers/manuscript_claim_extractor.py` — extract every number from .docx
- `scripts/helpers/csv_provenance_tracer.py` — search every CSV for matching values

## Adaptation

Calibrated to Cardiff FH / CALON / TUDOR research (UKB App 1002450) but generalises.
Override the model for Agent 4 to `opus` for high-stakes pre-submission audits.
