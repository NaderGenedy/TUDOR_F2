---
name: tudor-submission-coordinator
description: End-to-end TUDOR manuscript submission coordinator. Use whenever the user needs to COORDINATE the submission process, check what's ready vs what's blocking, manage the TRIPOD-AI checklist, prepare cover letters, format for a target journal, track open verification items, or get a status of the whole submission pipeline. Knows the 8 open verification items and the critical path to submission.
---

# TUDOR Submission Coordinator

You coordinate the entire path from manuscript skeleton to submitted paper.

## Current Status
- Manuscript skeleton: COMPLETE (MANUSCRIPT_SKELETON_TUDOR.md)
- R pipeline: 30 scripts COMPLETE
- Skills: 13/13 COMPLETE
- Model specification: LOCKED (11 variables, TUDOR_coefficients_locked.csv)
- Presentation materials: COMPLETE (15-min + 30-min + Word docs)

## 8 Open Verification Items (resolve before submission)

| ID | Item | Severity | Status |
|---|---|---|---|
| O-1 | Confirm model spec vs locked CSV (are Trig_Filter/Index Effect separate fitted terms?) | BLOCKING | OPEN |
| O-2 | Pull 95% CIs for every AUC from source | BLOCKING | OPEN |
| O-3 | Leakage/temporal check on all 11 predictors | BLOCKING | OPEN |
| O-4 | Reconcile Ns: 3,099 train / 1,471 Design-A / 506,506 total | BLOCKING | OPEN |
| O-5 | ClinVar version + APOB/PCSK9 + CNV handling | IMPORTANT | OPEN |
| O-6 | Verify Fig 3c subgroup heatmap cell values vs source | IMPORTANT | OPEN |
| O-7 | Resolve 2/68 failing reproducers | BLOCKING | OPEN |
| O-8 | Prior-art / benchmark retrieval for Discussion | IMPORTANT | OPEN |

## Critical Path to Submission
1. Resolve O-1 + O-7 (model spec + reproducers) — BLOCKING
2. Lock all numbers + CIs from source (O-2) → tables/figures regenerate
3. Reconcile Ns + CONSORT (O-4)
4. Evidence/comparator table (O-8) → draft Intro/Discussion
5. Internal review (lipid-cardiology-reviewer pass) → revise
6. Cover letter + declarations + TRIPOD-AI checklist → submit

## Target Journals (primary + cascade)
1. **Primary:** Circulation or EHJ (strong CV science + population scale + mechanism)
2. **Cascade 1:** J Clin Lipidol (lipid-focused, strong fit)
3. **Cascade 2:** Heart / EJPC (pragmatic, good for solid work)

## Pre-Submission Gate (all must pass)
- [ ] Every manuscript number traces to a source row
- [ ] Reproducer suite: 0 FAIL
- [ ] TRIPOD-AI checklist completed
- [ ] Cover letter + declarations ready
- [ ] Target-journal formatting (word/figure limits, reference style)
- [ ] All co-authors have reviewed and approved

## Reporting Standards
- TRIPOD+AI for prediction model
- PROBAST self-appraisal for risk of bias
- RECORD for EHR/biobank data
- Pre-registered head-to-head with locked sign convention
