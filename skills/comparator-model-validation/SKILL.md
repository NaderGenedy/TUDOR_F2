---
name: comparator-model-validation
description: >-
  Use whenever a risk model, score, or prediction equation is being compared against a published
  comparator and any claim of the form "our model beats / outperforms / is better than X" is being
  made, computed, or written into a manuscript or abstract. Triggers on "outperforms", "beats
  SAFEHEART", "beats Montreal-FH-SCORE", "better than DLCN/FAMCAT/PCE/SCORE2", "AUC advantage",
  "+0.1 AUC", "improved NRI", "head-to-head", "comparator model", "external validation against",
  "superior in all subgroups", or before any "our model is better" sentence ships. Enforces frozen
  published comparator coefficients, like-for-like out-of-sample comparison, an explicit pre-registered
  sign convention, full (uncherry-picked) subgroup tables, real NRI/IDI computation, and traceability
  of every reported metric to a live results file. Invoke this BEFORE writing the comparison, not after.
---

# Comparator Model Validation

## Why this skill exists

"Our model outperforms the established score" is one of the most damaging claims a prediction
paper can get wrong. It is also one of the easiest to get wrong *by accident* — not through
fraud, but through a chain of small, individually plausible mistakes that all push the result
the same way: in your model's favour.

This skill was distilled from a real five-agent QC of an FH MACE-prediction model. The manuscript
claimed the model "outperformed SAFEHEART and Montreal-FH-SCORE in all subgroups with an AUC
advantage of +0.1 and impressive NRI." Independent recomputation found:

- the model **lost** to SAFEHEART in 31 of 34 subgroups (mean delta -0.086, not +0.1);
- Montreal-FH-SCORE had **never been implemented** as a head-to-head comparator;
- the "impressive NRI" had **never been computed** (the manuscript itself said so elsewhere);
- the headline C-statistics were **hard-coded literals** that did not match the live results file;
- the "published" model's coefficient table **traced to no script and no CSV**.

None of that required bad intent. It required nobody checking. This skill is the check.

## The seven failure modes — and the rule that kills each

| # | Failure mode | What it looks like | The rule |
|---|---|---|---|
| 1 | **In-sample comparator re-fitting** | The comparator is re-fitted on your own data ("SAFEHEART-like model, refit") so it is handicapped against your model, which was tuned on the same data. | The comparator MUST use the **frozen, published coefficients** from its original paper. Never re-fit a comparator. If you re-fit it, it is no longer that comparator. |
| 2 | **Sign-flipped deltas** | The manuscript reports `delta = +0.02` favouring your model; the results CSV holds `-0.02`. Somewhere a subtraction ran the other way. | Pre-register one sign convention: `delta = metric(your_model) - metric(comparator)`, positive = your model wins. Write it down before computing. Verify every reported delta against the CSV's own delta column. |
| 3 | **Orphan published model** | The manuscript's Table of coefficients (the central deliverable) is reproduced by no script and stored in no file. A predictor appears in the table that exists in no model artefact. | Every published coefficient MUST come from a single frozen coefficient file, produced by one named script, re-runnable end to end. If a coefficient cannot be traced, the model is not real yet. |
| 4 | **Asserted-but-uncomputed NRI/IDI** | "Impressive NRI" / "significant reclassification" in the abstract, but no NRI computation exists anywhere in the repo. | If you claim NRI or IDI, the computation must exist, write its output to a file, and that file is the only source for the number. No computation, no claim. |
| 5 | **"Beats X in all subgroups" that reverses** | A sweeping superiority claim that, on recomputation, is the opposite of the data — often because only the favourable strata were inspected. | "All subgroups" means every cell of the full subgroup table favours your model. Show the whole table. One losing cell kills the word "all". |
| 6 | **Crippled comparator** | The comparator is implemented with a reduced predictor set (e.g. SAFEHEART without its Lp(a) term), so it underperforms its true published self. | Implement the comparator's **complete published predictor set**. A comparator missing predictors is a strawman, and reviewers who know the score will catch it. |
| 7 | **Hard-coded headline metrics** | The C-statistic in the abstract is a string literal typed into a figure script; the live results file says something else. | Every headline metric in the manuscript MUST be read at build time from a live results file. If a number is typed by hand, it is wrong the moment the pipeline changes. |

## The workflow

Run these steps in order. Do not skip step 0 — pre-registration is what makes the rest honest.

### Step 0 — Pre-register the comparison (before computing anything)

Write a short pre-registration file (`comparison_preregistration.md`) stating, *before any metric is
computed*:

- **Comparators**: exact named published scores, with citation (e.g. "SAFEHEART-RE, Perez de Isla
  2017 Circulation"; "Montreal-FH-SCORE, Paquette 2017 J Clin Lipidol").
- **Primary metric** and how superiority is defined (e.g. "delta C-statistic, out-of-sample,
  positive favours our model; superiority = delta > 0 with 95% CI excluding 0").
- **Subgroups**: the exact list, fixed now, so they cannot be chosen after seeing results.
- **Cohorts and the validation scheme**: which cohort develops, which validates, and confirmation
  that headline comparisons are **out-of-sample**.
- **Sign convention**: one sentence, e.g. `delta = AUC(ours) - AUC(comparator)`.

This file is the contract. If the result later contradicts the abstract's framing, the
pre-registration is what tells you the abstract is wrong, not the data.

### Step 1 — Freeze the comparator(s)

- Obtain the **exact published coefficients / points** from the comparator's primary paper
  (usually the supplement). See `references/published_fh_risk_equations.md` for predictor sets
  and citations for the common FH scores. **Do not reconstruct coefficients from memory** —
  extract them from the source and cite the table/equation number.
- Implement the **full** predictor set. If a predictor is unavailable in your data, that is a
  limitation to report — not a licence to drop it silently.
- Write the comparator as a pure function: data in, risk out, no fitting. Save the frozen
  coefficients to a committed file (`comparator_<name>_coefficients.csv`).

### Step 2 — Freeze your own model

- Fit your model with a named, re-runnable script.
- Write its coefficients to one committed file (`<model>_coefficients.csv`).
- The manuscript's coefficient table must be generated *from that file*, never typed.
- Confirm events-per-variable is adequate (>=10 is the usual floor); an overfit model on a small
  cohort inflates apparent discrimination and is the commonest reason a "win" evaporates on
  external data.

### Step 3 — Like-for-like evaluation

For every reported comparison, both models must be scored on:

- the **same patients** (identical inclusion, identical complete-case / imputation handling),
- the **same outcome** (identical event definition and follow-up window),
- the **same n** (report it; if n differs between models, the comparison is invalid),
- **out-of-sample** data for any headline claim (leave-one-cohort-out or external validation).
  An in-sample AUC for your model versus anything is not a headline result.

Check family-level structure: if cohorts share families (e.g. cascade-screened relatives, or a
sub-cohort nested inside a registry), deduplicate by family before pooling, and ensure no family
appears in both training and validation folds — a leak inflates the apparent win.

### Step 4 — Compute deltas with the pre-registered sign

- `delta = metric(ours) - metric(comparator)` for every cell.
- Bootstrap the 95% CI of the delta (paired, B >= 2000) — a delta without a CI cannot support
  "superior".
- Verify each delta's sign against the source CSV's own delta column. A disagreement here is
  failure mode 2.

### Step 5 — Build the full subgroup table

- One row per subgroup, every subgroup from the pre-registration, no omissions.
- Columns: subgroup, n, events, metric(ours), metric(comparator), delta, 95% CI.
- Count the cells where your model wins, ties, loses. State all three counts.
- The phrase "in all subgroups" is permitted **only** if the loss count is zero. Otherwise report
  honestly: "in X of Y subgroups", and describe where and why it loses.

### Step 6 — NRI / IDI (only if you will claim them)

- If the abstract or results will mention reclassification, NRI, or IDI, the computation must
  exist as code, run, and write to a file.
- Use a categorical-NRI implementation cross-checked against a reference (R `nricens` /
  `PredictABEL`, or a verified Python equivalent). Watch the two classic errors: probabilities
  not on the same scale before cutting into risk categories, and event/non-event components
  added with the wrong sign.
- If NRI was not computed, the word NRI does not appear in the manuscript.

### Step 7 — Traceability sweep

- Every metric in the manuscript (every AUC, C-statistic, delta, NRI, calibration slope, n,
  event count) must be readable from a live results file.
- Grep the figure / table / manuscript-generator scripts for numeric string literals. Any
  headline metric found as a literal is a hard-coded value (failure mode 7) — replace it with a
  read from the results file.
- Calibration slope is the coefficient of `glm(outcome ~ linear_predictor)`; a slope far from 1
  reported as evidence of good calibration is itself a finding.

### Step 8 — Honest reporting

- If the model wins cleanly, like-for-like, out-of-sample, with CIs excluding zero — say so
  plainly; it is a strong result.
- If it wins overall but loses in some strata — report exactly that; partial superiority is still
  publishable and is far stronger than an "all subgroups" claim that collapses at review.
- If the comparator wins — that is a legitimate, publishable finding too. A well-validated score
  beating a new one in a specific cohort is real evidence. Reframing a loss as a win is the one
  outcome that ends careers; reporting a loss honestly does not.

## Output

Produce `comparator_validation_report.md` with:

- the pre-registration (verbatim, from step 0),
- the full subgroup table (step 5) with win/tie/loss counts,
- a traceability ledger: every headline metric -> source file -> value,
- a one-line verdict per claim: SUPPORTED / PARTIALLY SUPPORTED / NOT SUPPORTED / NOT COMPUTABLE,
- and, if any claim is NOT SUPPORTED, the corrected wording the manuscript should use instead.

## What this skill does not do

- It does not decide that your model is good or bad — the data does.
- It does not edit the manuscript without explicit user direction — it reports and proposes wording.
- It does not replace TRIPOD; it is the comparison-specific layer on top of TRIPOD reporting.

## Reference

- `references/published_fh_risk_equations.md` — predictor sets, citations, and coefficient-extraction
  guidance for SAFEHEART-RE, Montreal-FH-SCORE, DLCN, FAMCAT and related FH risk tools.
