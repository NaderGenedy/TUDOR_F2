---
name: ukb-data-audit
description: Use this skill when the user wants to audit UK Biobank raw data extracts, hunt for bugs in analytic code lists, trace numerical claims in a manuscript back to source CSV rows, verify cohort filters, or check that nothing is "missing" before raising a UKB portal amendment. Triggers explicitly on requests like "run UKB audit", "audit my UKB data", "any bugs?", "trace this manuscript", "check what's actually on disk", "verify the cohort", "before I file the amendment, find everything". Five independent agents run in parallel: inventory, field verification, code-list completeness, manuscript provenance, and cohort consistency. Catches the class of bug where the ledger says "missing" but the file is in another project folder; where p24100 is described as cardiac-CT but is actually CMR; where K611 (BAV, 3,000+ events) is omitted from an OPCS-4 list; and where a manuscript HR doesn't match any CSV row to within tolerance.
---

# UKB Data Audit — 5-agent parallel audit

This skill exists because UKB analytic projects accumulate bugs that cannot be caught by any single linter or single read of any single script. Files migrate between project folders; UKB Showcase field IDs get mis-labelled in extraction scripts; code lists drift from the underlying HES distribution; manuscript values drift from their source CSVs; cohort filters silently differ between scripts that claim to use the same cohort.

The audit dispatches five independent agents that each look at a different facet, then a sixth consolidator surfaces only the contradictions and critical findings. Independence matters — running the checks in series in one context window means later checks anchor on earlier ones and miss the contradictions.

## When to use

The user says one of:
- "run UKB audit" / "audit my UKB data"
- "any bugs?" / "what's still missing?"
- "trace this manuscript" / "verify every number"
- "before I file the amendment, find everything that's already on disk"
- "is the cohort consistent across the v22 scripts?"

Or the context shows the user about to commit a major UKB deliverable (manuscript submission, amendment filing, cohort lock) and they would benefit from a final audit pass — in that case, propose running this skill and let the user opt in.

## What the five agents do

Each agent has detailed instructions in `agents/`. Read the agent prompt before dispatching. The agents are independent — they share no state, only the canonical inventory from A1 (when downstream agents need it).

| # | Agent | Inputs | Output |
|---|---|---|---|
| **A1** | Inventory — walk D:, E:, C: drives + backups; produce canonical CSV of every extract with `path, size_MB, n_pfields, pfields`. Uses `D:/Projects/Lpa_Multilevel/data/UKB_amendment_audit/01_inventory_all_extracts.csv` as cache if mtime < 7 days old. | none | `iteration-N/A1_inventory/inventory.csv` + `report.md` |
| **A2** | Field verification — for every unique p-field, cross-check against the UKB Showcase ground-truth in `references/ukb_field_dictionary.md`; flag mis-labelled fields like p24100 (CMR LV, NOT cardiac CT) | A1 output | `iteration-N/A2_field_verify/verified.csv` + `report.md` |
| **A3** | Code-list completeness — frequency-rank ICD-10 / OPCS-4 / ATC codes in raw HES; flag high-frequency codes missing from analytic scripts (catches the K611 / BAV omission class) | A1 output + path to analytic scripts dir | `iteration-N/A3_code_list/missing_codes.csv` + `report.md` |
| **A4** | Manuscript provenance — extract every numerical claim from a target manuscript (HR/CI/p/n); for each, find the source CSV row that produced it within tolerance | path to manuscript | `iteration-N/A4_provenance/per_claim_audit.csv` + `report.md` |
| **A5** | Cohort consistency — verify v22 cohort filters (`analytic_nonFH==1`, `prevalent_X==0`, FamilyNumber dedup); replicate v22 A1 baseline HR 1.121 as sanity gate; check train/val FamilyNumber overlap | A1 output + path to cohort builder | `iteration-N/A5_cohort/consistency.csv` + `report.md` |
| **C** | Consolidator — read all five reports; surface only contradictions and critical findings. NEVER repeat what individual agents already said cleanly — only what NEEDS user attention. | A1-A5 outputs | `iteration-N/summary.md` |

## How to dispatch

This skill assumes parallel subagents are available (Task tool with `subagent_type=general-purpose` or specialised types where appropriate). Dispatch all five A-agents **in a single message with five Task tool calls** — independence requires no shared context.

The consolidator runs AFTER all five complete. Do not dispatch the consolidator in parallel; it needs all five reports.

Before dispatch:
1. Create `workspace/iteration-N/` (N = next free integer)
2. Read the active manuscript path from the user, or default to the most recent `.docx` under `D:/Projects/Lpa_Multilevel/manuscript_NEJM/`
3. Read the active scripts dir, or default to `D:/Projects/Lpa_Multilevel/scripts/`

After dispatch:
1. Wait for all five to complete
2. Dispatch consolidator with all five `report.md` paths
3. Print consolidator summary to the user
4. Offer next action (filing amendment, locking manuscript, etc.)

## Output discipline

Reports are written to disk, not printed in chat. Chat output is a 5-line top-level summary plus a pointer to `workspace/iteration-N/summary.md`. The user reads the detailed reports themselves if they want.

Per the user's CLAUDE.md, manuscript-touching deliverables produce md + docx side-by-side. The consolidator does this when the user explicitly asks to "lock" the audit.

## Bug class this skill exists to catch

1. **Cross-project hiding** — file claimed missing is actually held in a different project. (Real example: `batch7_opcs4.csv` in `perimeopause/` while v9 ledger said OPCS-4 missing.)
2. **Field mis-labelling** — extraction script names a p-field for what the user wants, but UKB Showcase says it's something else. (Real example: p24100 labelled "Aortic valve calcium" in script; actually LVEDV CMR.)
3. **Code-list incompleteness** — frequency-distribution of HES codes shows high-frequency codes absent from analytic lists. (Real example: K611 balloon valvuloplasty, 3,088 events, missing from initial OPCS-4 list.)
4. **Manuscript drift** — value in prose ≠ value in source CSV after a re-run. (Per CLAUDE.md non-negotiable #5.)
5. **Cohort drift** — `analytic_nonFH==1` applied in script A but not in script B; FamilyNumber overlap unintentionally retained.

## Limitations and known issues

- Agent A4 (provenance) requires the manuscript to be in `.docx` or `.md`. PDF parsing is out of scope; convert first.
- Agent A2 ground-truth is the offline dictionary in `references/ukb_field_dictionary.md`. Live UKB Showcase fetches are not in scope (rate-limited, often hangs). Keep the dictionary up to date.
- Agent A5 sanity-replicates v22 A1 only. Adding other replication targets is a per-project extension.

## Iteration discipline

Per CLAUDE.md "Diminishing returns" red flag: after iteration 3 of audit fixes, stop. If something is still surfacing as a contradiction it is either a real ongoing issue worth manuscript discussion or a parametric blind spot the skill cannot resolve. Add to the manuscript's pre-stated limitations rather than running iteration 4.
