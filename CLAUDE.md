# TUDOR — Project Instructions

## What this project is
TUDOR is a treatment-adjusted, ascertainment-aware FH prediction model for familial hypercholesterolaemia. It predicts the probability of genetically confirmed FH / LDLR carrier status and helps prioritise genetic testing. It does NOT replace genetic testing and is NOT an ASCVD event-risk model.

## Key results
- Bidirectional external validation in 506,506 adults (Welsh registry + UK Biobank)
- 25/25 prespecified subgroups won vs DLCN, MEDPED, Simon Broome
- Welsh AUC 0.770 vs DLCN 0.670; UKB AUC 0.631 vs DLCN 0.538

## House rules
- British English throughout
- Numbers/criteria before adjectives
- NEVER assert cohort sizes, coefficients, or AUCs from memory — verify against locked files
- Every claim must match the design (prediction, not causal)
- Verb discipline: "predicts/identifies", never "causes/drives"
- Concede limitations honestly; own the eDLCN unfairness

## Data governance
- UK Biobank Application 1002450 — participant data NEVER committed to git
- NHS Wales — same governance applies
- Set `TUDOR_DATA_DIR` env var to point at your data folder

## The locked model (from TUDOR_coefficients_locked.csv)
11 variables, elastic-net logistic regression (alpha=0.5, C=1.0, saga, seed 20260518, intercept -1.851). Coefficients are FROZEN — never re-estimated on validation data.

## Skills available (13)
Skills in `.claude/skills/` auto-load. Use `/skill-name` or describe what you need.

## Pipeline
30 R scripts (01-30) run sequentially. Each depends on prior outputs. Set `TUDOR_DATA_DIR` before running.
