#!/bin/bash
# ==============================================================================
# TUDOR PIPELINE: STEP 00 — UKB-RAP DATA EXTRACTION
# ==============================================================================
# PURPOSE: Extract additional fields from UK Biobank RAP that were missing
#          from the original analysis. Run this ON THE UKB-RAP PLATFORM
#          (JupyterLab Bash terminal or %%bash cell), NOT locally.
#
# OUTPUTS: 5 CSV files uploaded to your project root for download.
#
# AFTER RUNNING: Download all 5 files to C:/Users/nader/Downloads/
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

# --------------------------------------------------------------------------
# 1. Lp(a) — Field 30790 (Lipoprotein(a), nmol/L)
# --------------------------------------------------------------------------
echo "=== Extracting Lp(a) ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,participant.p30790_i0 \
  --output ukb_lpa.csv --delimiter ","
echo "Done: ukb_lpa.csv"

# --------------------------------------------------------------------------
# 2. ASCVD History — Field 6150 (Vascular/heart problems, arrays 0-3)
# --------------------------------------------------------------------------
echo "=== Extracting ASCVD History ==="
python3 -c "
fields = ['participant.eid']
for a in range(4):
    fields.append(f'participant.p6150_i0_a{a}')
print('\n'.join(fields))
" > /tmp/fields_ascvd.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_ascvd.txt \
  --output ukb_ascvd.csv --delimiter ","
echo "Done: ukb_ascvd.csv"

# --------------------------------------------------------------------------
# 3. MI/Angina Age — Fields 3894 (MI age) and 3627 (Angina age)
# --------------------------------------------------------------------------
echo "=== Extracting MI/Angina Age ==="
dx extract_dataset "$RECORD" \
  --fields participant.eid,participant.p3894_i0,participant.p3627_i0 \
  --output ukb_cvd_age.csv --delimiter ","
echo "Done: ukb_cvd_age.csv"

# --------------------------------------------------------------------------
# 4. Visit 0 Medications Part A (arrays 0-19)
#    Needed to identify who was statin-FREE at baseline
# --------------------------------------------------------------------------
echo "=== Extracting Visit 0 Meds (Part A) ==="
python3 -c "
fields = ['participant.eid']
for a in range(20):
    fields.append(f'participant.p20003_i0_a{a}')
print('\n'.join(fields))
" > /tmp/fields_meds0a.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_meds0a.txt \
  --output ukb_meds_v0_a.csv --delimiter ","
echo "Done: ukb_meds_v0_a.csv"

# --------------------------------------------------------------------------
# 5. Visit 0 Medications Part B (arrays 20-39)
# --------------------------------------------------------------------------
echo "=== Extracting Visit 0 Meds (Part B) ==="
python3 -c "
fields = ['participant.eid']
for a in range(20, 40):
    fields.append(f'participant.p20003_i0_a{a}')
print('\n'.join(fields))
" > /tmp/fields_meds0b.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_meds0b.txt \
  --output ukb_meds_v0_b.csv --delimiter ","
echo "Done: ukb_meds_v0_b.csv"

# --------------------------------------------------------------------------
# 6. Upload all to project root
# --------------------------------------------------------------------------
echo "=== Uploading to project ==="
for f in ukb_lpa.csv ukb_ascvd.csv ukb_cvd_age.csv ukb_meds_v0_a.csv ukb_meds_v0_b.csv; do
  dx upload "$f" --destination / --brief
done

echo ""
echo "============================================"
echo " ALL EXTRACTIONS COMPLETE"
echo " Download these 5 files from the Manage tab:"
echo "   1. ukb_lpa.csv"
echo "   2. ukb_ascvd.csv"
echo "   3. ukb_cvd_age.csv"
echo "   4. ukb_meds_v0_a.csv"
echo "   5. ukb_meds_v0_b.csv"
echo "============================================"
