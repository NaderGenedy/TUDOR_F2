#!/bin/bash
# ==============================================================================
# TUDOR PIPELINE: STEP 00b — ADDITIONAL UKB-RAP EXTRACTIONS
# ==============================================================================
# PURPOSE: Extract fields needed by scripts 03 and 08 that were not included
#          in the original 00_extract_ukbrap.sh extraction.
#          Run this ON THE UKB-RAP PLATFORM (JupyterLab Bash terminal).
#
# OUTPUTS: 2 CSV files uploaded to your project root for download.
#
#   Script 03 (Longitudinal):  ukb_longitudinal.csv
#   Script 08 (TRG-Shield):    ukb_dm_hba1c.csv
#
# NOTE ON SCRIPT 09 (NMR Metabolomics):
#   Fields p23400–p23531 (Nightingale NMR, Category 220) are NOT dispensed
#   to this project. You must request NMR access via the UKB Access
#   Management System before script 09 can run. Once dispensed, re-run
#   the NMR extraction section at the bottom of this script.
#
# AFTER RUNNING: Download files from the Manage tab to C:/Users/nader/Downloads/
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

# ============================================================================
#  SCRIPT 03 — Longitudinal Validation
#  Needs: Visit 1 LDL, Visit 1 Cholesterol, Visit 1 Medications (40 arrays)
#
#  Split into 3 extractions to avoid "DataTooLarge" error, then merged
#  into a single ukb_longitudinal.csv via Python.
# ============================================================================

# --------------------------------------------------------------------------
# 1a. Visit 1 LDL + Cholesterol (2 fields — small, always succeeds)
# --------------------------------------------------------------------------
echo "=== [Script 03] Part 1/3: Visit 1 LDL + Cholesterol ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,participant.p30780_i1,participant.p30690_i1 \
  --output /tmp/ukb_v1_ldl.csv --delimiter ","
echo "Done: /tmp/ukb_v1_ldl.csv"

# --------------------------------------------------------------------------
# 1b. Visit 1 Medications arrays 0–19 (20 fields)
# --------------------------------------------------------------------------
echo "=== [Script 03] Part 2/3: Visit 1 Meds a0-a19 ==="
python3 -c "
fields = ['participant.eid']
for a in range(20):
    fields.append(f'participant.p20003_i1_a{a}')
print('\n'.join(fields))
" > /tmp/fields_meds_v1a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_meds_v1a.txt \
  --output /tmp/ukb_meds_v1_a.csv --delimiter ","
echo "Done: /tmp/ukb_meds_v1_a.csv"

# --------------------------------------------------------------------------
# 1c. Visit 1 Medications arrays 20–39 (20 fields)
# --------------------------------------------------------------------------
echo "=== [Script 03] Part 3/3: Visit 1 Meds a20-a39 ==="
python3 -c "
fields = ['participant.eid']
for a in range(20, 40):
    fields.append(f'participant.p20003_i1_a{a}')
print('\n'.join(fields))
" > /tmp/fields_meds_v1b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_meds_v1b.txt \
  --output /tmp/ukb_meds_v1_b.csv --delimiter ","
echo "Done: /tmp/ukb_meds_v1_b.csv"

# --------------------------------------------------------------------------
# 1d. Merge the 3 parts into ukb_longitudinal.csv on participant.eid
# --------------------------------------------------------------------------
echo "=== [Script 03] Merging 3 parts into ukb_longitudinal.csv ==="
python3 -c "
import pandas as pd

ldl  = pd.read_csv('/tmp/ukb_v1_ldl.csv')
meda = pd.read_csv('/tmp/ukb_meds_v1_a.csv')
medb = pd.read_csv('/tmp/ukb_meds_v1_b.csv')

merged = ldl.merge(meda, on='participant.eid', how='outer') \
            .merge(medb, on='participant.eid', how='outer')

merged.to_csv('ukb_longitudinal.csv', index=False)
print(f'Merged: {merged.shape[0]} rows x {merged.shape[1]} cols')
"
echo "Done: ukb_longitudinal.csv"

# ============================================================================
#  SCRIPT 08 — TRG-Shield Biology
#  Needs: Diabetes diagnosis, HbA1c
# ============================================================================

# --------------------------------------------------------------------------
# 2. ukb_dm_hba1c.csv  (already succeeded — re-run is safe/idempotent)
#    Fields: p2443_i0   (Diabetes diagnosed by doctor, Visit 0)
#            p30750_i0  (HbA1c, mmol/mol, Visit 0)
# --------------------------------------------------------------------------
echo "=== [Script 08] Extracting Diabetes + HbA1c ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,participant.p2443_i0,participant.p30750_i0 \
  --output ukb_dm_hba1c.csv --delimiter ","
echo "Done: ukb_dm_hba1c.csv"

# ============================================================================
#  UPLOAD TO PROJECT ROOT
# ============================================================================
echo ""
echo "=== Uploading files to project ==="
for f in ukb_longitudinal.csv ukb_dm_hba1c.csv; do
  dx upload "$f" --destination / --brief
done

echo ""
echo "============================================"
echo " EXTRACTIONS COMPLETE"
echo " Download these 2 files from the Manage tab:"
echo ""
echo "   Script 03:  ukb_longitudinal.csv"
echo "   Script 08:  ukb_dm_hba1c.csv"
echo ""
echo " NOTE: Script 09 NMR fields (p23400-p23531)"
echo " are NOT dispensed to this project."
echo " Request Category 220 access from UKB, then"
echo " re-run the NMR section below."
echo "============================================"

exit 0

# ============================================================================
#  SCRIPT 09 — NMR Metabolomics (DISABLED — fields not dispensed)
#  Uncomment and run AFTER Category 220 is dispensed to your project.
#  Fields p23400–p23531 (132 Nightingale NMR fields), 5 batches.
# ============================================================================

# echo "=== [Script 09] Extracting NMR Batch A (p23400-p23425) ==="
# python3 -c "
# fields = ['participant.eid']
# for f in range(23400, 23426):
#     fields.append(f'participant.p{f}')
# print('\n'.join(fields))
# " > /tmp/fields_nmr_a.txt
# dx extract_dataset "$RECORD" \
#   --fields-file /tmp/fields_nmr_a.txt \
#   --output ukb_nmr_a.csv --delimiter ","
# echo "Done: ukb_nmr_a.csv"

# echo "=== [Script 09] Extracting NMR Batch B (p23426-p23451) ==="
# python3 -c "
# fields = ['participant.eid']
# for f in range(23426, 23452):
#     fields.append(f'participant.p{f}')
# print('\n'.join(fields))
# " > /tmp/fields_nmr_b.txt
# dx extract_dataset "$RECORD" \
#   --fields-file /tmp/fields_nmr_b.txt \
#   --output ukb_nmr_b.csv --delimiter ","
# echo "Done: ukb_nmr_b.csv"

# echo "=== [Script 09] Extracting NMR Batch C (p23452-p23477) ==="
# python3 -c "
# fields = ['participant.eid']
# for f in range(23452, 23478):
#     fields.append(f'participant.p{f}')
# print('\n'.join(fields))
# " > /tmp/fields_nmr_c.txt
# dx extract_dataset "$RECORD" \
#   --fields-file /tmp/fields_nmr_c.txt \
#   --output ukb_nmr_c.csv --delimiter ","
# echo "Done: ukb_nmr_c.csv"

# echo "=== [Script 09] Extracting NMR Batch D (p23478-p23505) ==="
# python3 -c "
# fields = ['participant.eid']
# for f in range(23478, 23506):
#     fields.append(f'participant.p{f}')
# print('\n'.join(fields))
# " > /tmp/fields_nmr_d.txt
# dx extract_dataset "$RECORD" \
#   --fields-file /tmp/fields_nmr_d.txt \
#   --output ukb_nmr_d.csv --delimiter ","
# echo "Done: ukb_nmr_d.csv"

# echo "=== [Script 09] Extracting NMR Batch E (p23506-p23531) ==="
# python3 -c "
# fields = ['participant.eid']
# for f in range(23506, 23532):
#     fields.append(f'participant.p{f}')
# print('\n'.join(fields))
# " > /tmp/fields_nmr_e.txt
# dx extract_dataset "$RECORD" \
#   --fields-file /tmp/fields_nmr_e.txt \
#   --output ukb_nmr_e.csv --delimiter ","
# echo "Done: ukb_nmr_e.csv"

# echo "=== Uploading NMR files ==="
# for f in ukb_nmr_a.csv ukb_nmr_b.csv ukb_nmr_c.csv ukb_nmr_d.csv ukb_nmr_e.csv; do
#   dx upload "$f" --destination / --brief
# done
# echo "NMR upload complete."
