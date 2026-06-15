---
name: ukb-inventory
description: Use whenever a UK Biobank field, cohort size, outcome variable, NMR / PRS / imaging measurement, or result CSV needs to be located in Dr Genedy's working tree. Triggers on "where is [field]", "do I have data for", "load the master", "which tier has", "cohort n for", "is X in the master", and at the start of any new UK Biobank analysis. Single source of truth for the three-tier storage layout, the 501,936-row master v2 file (247 columns), the 32 result CSVs, headline cohort numbers, and the two currently-outstanding RAP extractions (LPA locus genotypes; UKB kinship / relatedness — fields 22011/22012/22018/22020/22021 + KING pairs file). Always cross-checks the candidate field against `D:\Projects\CALON_AlphaFold_Rebuild\New folder\fieldsum.tsv` (the UKB field catalogue) when in doubt about availability. Complements `ukb-preflight` (which handles environment + on-disk verification before running a pipeline).
---

# UKB-Inventory — single source of truth for the working tree

## When to invoke

Use this skill **before** any UK Biobank analysis step where you would otherwise ask the user "where is X". Triggers include:

- "Where is the [field / variable / outcome]?"
- "Do I have data for [exposure / outcome / mediator]?"
- "Load the master" / "load the UKB CSV" / "read the inventory"
- "What's the cohort n for [outcome]?"
- "Which tier has the [NMR / PRS / imaging] data?"
- "Is X already on disk?" (before scheduling any RAP extraction)
- Start of any new UK Biobank analysis or risk-prediction model

Pair with `ukb-preflight` when the next step is actually running a pipeline.

## Three storage tiers

| Tier | Path | What's there | Use it when |
|---|---|---|---|
| **Tier 1 — Working** | `D:\Projects\Lpa\` | Master v2 + project results | Daily analysis |
| **Tier 2 — Raw extracts (rich)** | `D:\Projects\CALON_AlphaFold_Rebuild\New folder\` | "Treasure trove" — 25 NMR chunks, 40 PCs, all first-occurrence dates, medications | When a UKB field is not in the master |
| **Tier 3 — Backup** | `D:\CALON_FH_BACKUP_FULL\ukb_reviewer\` | Older but complete reviewer-grade extracts | Sensitivity / sanity checks |

**Staging area** for transfer or RAP: `D:\Projects\Lpa\rap_staging\` (31 files / 5.55 GB; manifest in `MANIFEST.csv`).

## The single source of truth

```
D:\Projects\Lpa\ukb_FULL_MASTER_v2.csv     501,936 rows × 247 columns
```

This is what loads for analysis. Default to it unless a specific field is missing — then fall back to Tier 2.

### Column families in the master v2

| Family | Example columns | Source |
|---|---|---|
| Demographics | `sex_F`, `age_at_recruit`, `townsend`, `ethnicity_code` | UKB phenotype |
| Lipids (chemistry) | `tc_chem`, `ldl_chem`, `hdl_chem`, `tg_chem`, `lpa_chem`, `apob_chem` | Roche assays |
| NMR core (24 cols) | `nmr_tc`, `nmr_ldl`, `nmr_hdl`, `nmr_tg`, `nmr_apob`, `nmr_omega3` | Nightingale |
| NMR particles (16 cols) | `nmr_l_ldl_p`, `nmr_s_ldl_p`, `nmr_idl_p` and the rest | Nightingale (added) |
| PRS | `ldl_prs`, `cad_prs`, `bp_prs`, `hba1c_prs` | UKB standard |
| Principal components (40) | `pc_1` ... `pc_40` | UKB imputation (added) |
| Composite ASCVD outcomes | `first_ascvd`, `incident_ascvd`, `prevalent_ascvd`, `t_event_years` | HES + first-occurrence |
| Aortic stenosis outcomes | `as_first_date`, `as_ever`, `as_incident`, `t_as_years` | I35 first occurrence (added) |
| CV death | `cv_death` (T/F from ICD-10 cause-of-death array) | Death cause array (added) |
| Comorbidities | `t1dm`, `t2dm`, `has_af`, `has_hf`, `has_ckd`, `has_cancer` | HES |
| Imaging | `lvef`, `lvedv`, `lvm`, `cimt_120`, `liver_pdff_i2`, `vat_vol` | UKB MRI / ultrasound |
| Carrier-specific | `ldlr_carrier`, `sss_v3`, `ldlr_tier4`, `variant_id`, `domain_clean`, `f_ddG`, `f_dms_uptake`, `f_am` | CALON / AlphaFold |
| Statin (rich) | `statin_first_date_v2`, `statin_max_dose_mg`, `statin_intensity` | GP prescriptions (added) |
| Harmonised Lp(a) | `lpa_harm`, `lpa_harm_z`, `lpa_source` | Combined chem + NMR (added) |

If a field is in the column families above, **load it from `ukb_FULL_MASTER_v2.csv` directly**. If not, escalate to Tier 2 raw extracts.

### Canonical R loader

```r
library(data.table)
m <- fread("D:/Projects/Lpa/ukb_FULL_MASTER_v2.csv")
stopifnot(nrow(m) == 501936, ncol(m) >= 247)
```

### Canonical Python loader

```python
import pandas as pd
m = pd.read_csv(r'D:\Projects\Lpa\ukb_FULL_MASTER_v2.csv')
assert len(m) == 501_936 and m.shape[1] >= 247, f"unexpected shape {m.shape}"
```

For Tier 2 raw extracts (parquet), pin Python 3.12 with pyarrow available (see `ukb-preflight` skill).

## Headline cohort numbers (use verbatim in any methods section)

| Cohort | n |
|---|---:|
| Full UK Biobank | **501,936** |
| LDLR carriers (total) | **3,540** |
| ↳ severe-tier | 819 |
| ↳ moderate | 1,002 |
| ↳ mild | 615 |
| ↳ null | 458 |
| ↳ unclassified | 646 |
| Lp(a) measured | **375,200** |
| Incident aortic stenosis (any I35) | **41,796** |
| Incident composite ASCVD | 21,008 |
| All-cause death | 56,961 |
| Cardiovascular death | 21,889 |
| Median follow-up | 14.7 years |
| Total person-years | 7,067,921 |

Numbers above are authoritative — do not recompute on the fly without flagging that the canonical figure is being challenged.

## Pre-built results

`D:\Projects\Lpa\results\` contains **32 result CSVs**. Headline files:

| File | Content |
|---|---|
| `lpa_consequence_atlas_ALL.csv` | 504 rows: Lp(a) × outcome × stratum (5 families, full UKB) |
| `lpa_carrier_interaction.csv` | 72 rows: Lp(a) × LDLR-carrier interaction tests |
| `tissue_nmr_fullrisk_atlas.csv` | Carrier-focused atlas (v14 manuscript source) |
| `functional_retier_AS_interaction.csv` | 6-axis Fine-Gray (v14 Table 2) |
| `wales_dragon3_per_tertile_OR.csv` | Wales external replication |
| `G4_age_timescale_cox.csv` | Age-timescale Cox (v14 §3.8) |
| `NMR_imaging_lpa_partialcor.csv` | NMR × Lp(a) partial-correlations |

If a manuscript claim cites a number, the supporting row is in one of these CSVs (or in the per-paper reproducer outputs at `D:\Projects\Lpa_Multilevel\nejm-email\result_csvs\`).

## What is NOT yet on disk — outstanding RAP extractions

Two extractions remain after the May 2026 second-pass audit (the first audit missed the kinship fields).

### (1) LPA locus genotypes — for one-sample MR / colocalisation
- **`.bgen` chromosome 6, region 160,500,000 – 161,600,000** (covers the *LPA* locus)
- Minimum 8 SNPs: `rs10455872`, `rs3798220`, `rs140570886`, `rs41272110`, `rs143431368`, `rs41259144`, `rs55730499`, `rs186696265`
- Needed for: one-sample Mendelian randomisation, colocalisation, allele-specific sensitivity
- RAP cost: ~£1 compute, ~10 minutes
- Output target: `D:\Projects\Lpa\rap_staging\10_lpa_extras\ukb_lpa_genotypes.csv`

### (2) UKB kinship / relatedness data — for FH-relative identification in UKB
The fieldsum.tsv catalogue confirms all five fields are extractable but **none are yet on disk** in master v2 or Tier 2:

| Field ID | Title | Why needed |
|---:|---|---|
| **22011** | Genetic relatedness pairing | Pairs of related EIDs — enables identification of UKB FH carriers' relatives |
| **22012** | Genetic relatedness factor | KING kinship coefficient per pair |
| **22018** | Genetic relatedness exclusions | Standard exclusion flag (one-from-each-kin-pair) |
| **22020** | Used in genetic principal components | Indicates relatedness-eligible status |
| **22021** | Genetic kinship to other participants | Binary "has kin in UKB" flag |

Plus the canonical KING-derived **pairs file** at `/Bulk/Genotype Results/Genotype calls/ukbXXXXX_rel_sNNNNN.dat` on RAP.

- Needed for: TUDOR augmentation experiments (UKB Is_Relative flag); any analysis at risk of effective-N inflation from related individuals
- RAP cost: ~£2 compute, ~5–10 minutes
- Output target: `D:\Projects\Lpa\rap_staging\11_ukb_kinship\ukb_kinship_pairs.csv` + `ukb_kinship_fields.csv`
- Extraction script: `D:\Projects\Lpa_Multilevel\nejm-email\tudor_qc\2026-05-14-pass-ukb-relatives\extract_ukb_kinship.sh`

**Do not schedule any other RAP extraction** without first confirming the field is genuinely absent from both the master v2 *and* Tier 2 raw extracts. The previous "treasure-trove" audit revealed that 4 of 5 originally-planned extractions were already on disk — but it also missed kinship, demonstrating that the audit must check fieldsum.tsv (the UKB field catalogue at `D:\Projects\CALON_AlphaFold_Rebuild\New folder\fieldsum.tsv`) rather than just local-disk filename patterns.

## Companion skills (which to chain with this one)

| When | Then call |
|---|---|
| Before running any 100k+ row pipeline | `ukb-preflight` (environment + on-disk inventory verification) |
| To trace a manuscript number back to source | `ukb-data-audit` |
| At submission gate, certify every claim | `manuscript-qc` or `verify-manuscript` |
| To bundle files for transfer | `stage-files` |
| For any Cox / risk-prediction work | `cox-analysis` |

## Standard workflow lookup

| Scenario | Action |
|---|---|
| New analysis idea | Read the relevant outcome / exposure column-family block above; write the model against the master v2 |
| Field not in master | Check `D:\Projects\Lpa\results\UKB_INVENTORY.md` (443-line catalogue) or escalate to Tier 2 |
| QC a manuscript | `manuscript-qc` (build / run reproducer against locked CSVs) |
| Pre-flight a pipeline run | `ukb-preflight` |
| Stage files for transfer | `stage-files`, then Robocopy from `rap_staging\` |
| Pull a new RAP field | Add a job to `D:\Projects\Lpa\scripts\` modelled on `rap_extract_*.sh`; output to `rap_staging\10_lpa_extras\` |

## One-page mental map

```
ukb_FULL_MASTER_v2.csv  (501,936 x 247)
   |
   |-- lpa_chem / lpa_harm / lpa_harm_z       <-- exposure
   |-- ldlr_carrier (0/1) + ldlr_tier4        <-- modifier
   |-- sss_v3 / sss_v3_z + 6 severity axes    <-- carrier severity
   |-- as_incident + t_as_years               <-- primary outcome (Cox)
   |-- incident_ascvd + t_event_years         <-- composite outcome
   |-- cv_death / death_date                  <-- competing event
   |-- 16 NMR particle counts + 24 core       <-- mechanism
   |-- 40 PCs + genetic_sex                   <-- genetic adjustment
   |-- statin_first_date_v2 + intensity       <-- treatment
   `-- all imaging + chemistry                <-- biomarker atlas inputs
                              |
                              v
        D:\Projects\Lpa\results\*.csv   (32 files; atlas + interactions + per-tier + sensitivity)
                              |
                              v
        D:\Projects\Lpa\manuscript\v14_FINAL.md / .docx
                              |
                              v
        D:\Projects\Lpa\qc\<timestamp>\audit_report.md   (via manuscript-qc / verify-manuscript)
```

## Common pitfalls (avoid)

- **Reaching for raw extracts before checking the master.** Default to `ukb_FULL_MASTER_v2.csv`; escalate only when a column genuinely isn't there.
- **Scheduling a RAP extraction before the local-disk audit.** The previous round saved days of compute and storage by discovering most planned pulls were already local.
- **Recomputing cohort headline numbers ad hoc.** The numbers in the table above are the canonical figures — if a fresh computation disagrees, flag the discrepancy explicitly.
- **Mixing tier 1 and tier 3 in the same analysis.** Tier 3 is for sensitivity / sanity checks only; the primary analysis must use tier 1 unless explicitly stated.

## Reference exemplar

The Lp(a) consequence atlas (v14, April 2026) was assembled entirely from the master v2 plus the 32 result CSVs above, with one downstream certification pass via `manuscript-qc`. No additional RAP pull was required to deliver the manuscript.
