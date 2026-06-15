# TUDOR_F2 — Complete Project Setup Guide

## How to Download and Restart the Entire Project

### Quick Start (one command)
```bash
git clone https://github.com/NaderGenedy/TUDOR_F2.git
cd TUDOR_F2
```

### What's in this repository

```
TUDOR_F2/
├── .claude/skills/              # 13 research skills (Claude Code auto-loads these)
│   ├── structured-brainstorming.md
│   ├── preventive-cardio-epidemiologist.md
│   ├── cardiometabolic-biostatistician.md
│   ├── reproducibility-engineer.md
│   ├── cardiometabolic-evidence-synthesis.md
│   ├── lipid-cardiology-reviewer.md
│   ├── academic-medical-writer.md
│   ├── scientific-figure-designer.md
│   ├── research-executive-planner.md
│   ├── autoresearch.md
│   ├── tudor-data-engineer.md          # NEW: data pipeline specialist
│   ├── tudor-head-to-head-analyst.md   # NEW: subgroup/comparator analysis
│   └── tudor-submission-coordinator.md # NEW: end-to-end submission workflow
│
├── R_pipeline/                  # 30 R analysis scripts (numbered 01-30)
│   ├── 00_extract_ukbrap.sh     # UKB-RAP data extraction
│   ├── 00b_extract_ukbrap_extras.sh
│   ├── 00c_extract_nmr_dm.sh
│   ├── 00d_extract_ukbrap_secondary.sh
│   ├── 01_data_merge.R          # Data merge & feature engineering
│   ├── 02_external_validation.R # TRIPOD Type 4 validation
│   ├── ...                      # (full pipeline 03-30)
│   └── 30_reclassification_full.R
│
├── manuscript/
│   ├── MANUSCRIPT_SKELETON_TUDOR.md   # Full skeleton with all 10 skills applied
│   └── OPEN_VERIFICATION_ITEMS.md     # 8 items to resolve before submission
│
├── presentation/
│   ├── TUDOR_15min_Prize_Talk_Speaker_Notes.md
│   ├── TUDOR_NLA_30min_Talk_Training_Guide.md
│   ├── TUDOR_Speaker_Script_Fable_Cut.docx
│   └── TUDOR_Presentation_Notes.docx
│
├── model/
│   └── TUDOR_MODEL_SPECIFICATION.md   # Locked 11-variable model spec
│
├── results/
│   └── HEAD_TO_HEAD_SUMMARY.md        # All AUCs, deltas, subgroup wins
│
├── SETUP_GUIDE.md                     # THIS FILE
├── CLAUDE.md                          # Project instructions for Claude Code
└── .gitignore
```

### Required Data (NOT in the repo — you supply these)

These files contain participant-level data under UK Biobank DUA / NHS Wales governance and must never be committed to git.

```
# Set this env var to point to your data folder:
export TUDOR_DATA_DIR="D:\TUDOR_SUBMISSION"   # Windows
export TUDOR_DATA_DIR="/path/to/data"          # Mac/Linux

# Required files in that folder:
TUDOR_UKB_Features.csv          # Main UKB dataset (~400k participants)
tudor_v2_workspace.RData        # Wales FH Registry workspace
TUDOR_coefficients_locked.csv   # Locked model coefficients (11 variables)

# Optional (for extended analyses):
ukb_lpa.csv                     # Lipoprotein(a) data
ukb_ascvd.csv                   # ASCVD history
ukb_cvd_age.csv                 # MI/Angina age
ukb_meds_v0_a.csv               # Visit 0 medications Part A
ukb_meds_v0_b.csv               # Visit 0 medications Part B
```

### Running the Pipeline

#### R Pipeline (scripts 01-30)
```bash
# Install required R packages
Rscript -e 'install.packages(c("data.table","dplyr","pROC","ggplot2","mice","officer","flextable","randomForest","xgboost","e1071"))'

# Set data directory
export TUDOR_DATA_DIR="/path/to/your/data"

# Run sequentially (each script depends on prior outputs)
cd R_pipeline
Rscript 01_data_merge.R
Rscript 02_external_validation.R
# ... continue through 30
```

#### Using the Research Skills (Claude Code)
```bash
# Open the project in Claude Code
cd TUDOR_F2
claude

# Skills auto-load from .claude/skills/
# Invoke directly:
#   /lipid-cardiology-reviewer    — review a draft
#   /structured-brainstorming     — generate research ideas
#   /autoresearch                 — self-improving model loop

# Or describe what you need — skills activate by trigger words:
#   "review my methods section"   → lipid-cardiology-reviewer
#   "which test should I use"     → cardiometabolic-biostatistician
#   "plan my submissions"         → research-executive-planner
```

### The Winning Model — 11 Locked Coefficients

```
Intercept: -1.851 | Alpha: 0.5 | C: 1.0 | Solver: saga | Seed: 20260518

Variable              β (per 1-SD)   Idea
──────────────────────────────────────────────────
ldl_ul               +1.609          Lipid Age
on_statin            +0.357          Lipid Age
tc_chem              -0.743          Triglyceride Shield
tg_chem              -0.604          Triglyceride Shield
non_hdl              -0.151          Triglyceride Shield
hdl_chem             -0.137          Triglyceride Shield
age_per_decade       -0.357          Demographic
sex_F                -0.070          Demographic
premature_mace       +0.064          Proband Effect
t2dm_int             +0.047          Proband Effect
fam_hist_cvd          0.000          Proband Effect (L1-eliminated)

Engineered: Trig_Filter = ldl_ul / (tg_chem + 0.1)
            Index_Effect = (1 - Is_Relative) × ldl_ul
```

### Key Results Summary

```
Welsh validation:    TUDOR 0.770 vs DLCN 0.670 (Δ+0.099, p<0.001)
                     vs Simon Broome 0.570 (Δ+0.200)
                     vs MEDPED 0.553 (Δ+0.216)

UK Biobank (501,936): TUDOR 0.631 vs DLCN 0.538 vs SB 0.510 vs MEDPED 0.520

Subgroups:           25/25 prespecified — ALL positive
                     75/75 head-to-head comparisons — ZERO losses

T2DM subgroup:       TUDOR 0.718 vs DLCN 0.538 (n=45,715)
```

### Reproducing from Scratch

```bash
# 1. Clone
git clone https://github.com/NaderGenedy/TUDOR_F2.git && cd TUDOR_F2

# 2. Supply data
export TUDOR_DATA_DIR="/your/data/path"

# 3. Run pipeline
cd R_pipeline && for f in $(ls *.R | sort); do echo "Running $f..."; Rscript "$f"; done

# 4. Open in Claude Code for skills
cd .. && claude
```

### Contact
Nader Genedy MBBCh MRCP MRCPI SCE-AIM MRCGP PgDip
Department of Metabolic Medicine, Cardiff and Vale University Health Board
genedyn1@cardiff.ac.uk
