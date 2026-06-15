#!/bin/bash
# ==============================================================================
# TUDOR PIPELINE: STEP 00d — EXTRACT SECONDARY CAUSES OF HYPERCHOLESTEROLAEMIA
# ==============================================================================
# PURPOSE: Extract Field 20002 (self-reported non-cancer illness) to identify
#          secondary causes of hypercholesterolaemia for exclusion:
#            1. Hypothyroidism     (code 1226)
#            2. Type 2 Diabetes    (code 1223, also 1220 unspecified)
#            3. Nephrotic syndrome (code 1519)
#            4. Obstructive jaundice / cholestasis / biliary disease
#
#          Also extracts ICD-10 hospital data (Field 41270) if dispensed,
#          for more robust identification.
#
# RUN ON: UKB-RAP platform (JupyterLab Bash terminal)
#
# OUTPUTS:
#   1. ukb_selfreport_illness.csv — Field 20002, Instance 0, arrays 0-35
#   2. ukb_icd10_secondary.csv   — Field 41270 (HES ICD-10) if available
#
# AFTER RUNNING: Download to C:/Users/nader/Downloads/
#                Then run: 16_sensitivity_secondary_exclusion.R
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

echo "============================================"
echo " TUDOR: Secondary Causes Extraction"
echo " Field 20002 (Self-reported non-cancer illness)"
echo "============================================"
echo ""

# ==========================================================================
#  PART 1: Self-reported non-cancer illness (Field 20002)
#  Instance 0 (baseline visit), 36 arrays (a0-a35)
#  Split into 2 batches to avoid DataTooLarge errors
# ==========================================================================

# --------------------------------------------------------------------------
# 1a. Field 20002, arrays 0-17 (18 fields)
# --------------------------------------------------------------------------
echo "=== Part 1a: Field 20002, arrays 0-17 ==="
python3 -c "
fields = ['participant.eid']
for a in range(18):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/fields_illness_a.txt
echo "  Fields: $(wc -l < /tmp/fields_illness_a.txt) (including eid)"

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_illness_a.txt \
  --output /tmp/ukb_illness_a.csv --delimiter ","
echo "Done: /tmp/ukb_illness_a.csv ($(wc -l < /tmp/ukb_illness_a.csv) lines)"

# --------------------------------------------------------------------------
# 1b. Field 20002, arrays 18-35 (18 fields)
# --------------------------------------------------------------------------
echo "=== Part 1b: Field 20002, arrays 18-35 ==="
python3 -c "
fields = ['participant.eid']
for a in range(18, 36):
    fields.append(f'participant.p20002_i0_a{a}')
print('\n'.join(fields))
" > /tmp/fields_illness_b.txt
echo "  Fields: $(wc -l < /tmp/fields_illness_b.txt) (including eid)"

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_illness_b.txt \
  --output /tmp/ukb_illness_b.csv --delimiter ","
echo "Done: /tmp/ukb_illness_b.csv ($(wc -l < /tmp/ukb_illness_b.csv) lines)"

# --------------------------------------------------------------------------
# 1c. Merge parts A + B
# --------------------------------------------------------------------------
echo "=== Merging Field 20002 parts ==="
python3 -c "
import pandas as pd

a = pd.read_csv('/tmp/ukb_illness_a.csv')
b = pd.read_csv('/tmp/ukb_illness_b.csv')
merged = a.merge(b, on='participant.eid', how='outer')
merged.to_csv('ukb_selfreport_illness.csv', index=False)
print(f'Merged: {merged.shape[0]} rows x {merged.shape[1]} cols')
"
echo "Done: ukb_selfreport_illness.csv"

# ==========================================================================
#  PART 2: ICD-10 Hospital Episode Statistics (Field 41270)
#  This field may or may not be dispensed to your project.
#  If extraction fails, the R script will fall back to self-reported only.
# ==========================================================================

echo ""
echo "=== Part 2: ICD-10 codes (Field 41270) — may fail if not dispensed ==="

# Field 41270 has many arrays (up to ~250). Extract first 100 arrays.
python3 -c "
fields = ['participant.eid']
for a in range(100):
    fields.append(f'participant.p41270_a{a}')
print('\n'.join(fields))
" > /tmp/fields_icd10.txt

dx extract_dataset "$RECORD" \
  --fields-file /tmp/fields_icd10.txt \
  --output ukb_icd10_secondary.csv --delimiter "," 2>/dev/null

if [ -f "ukb_icd10_secondary.csv" ] && [ "$(wc -l < ukb_icd10_secondary.csv)" -gt 1 ]; then
  echo "  OK: ukb_icd10_secondary.csv ($(wc -l < ukb_icd10_secondary.csv) lines)"
else
  echo "  NOTE: ICD-10 field not dispensed or extraction failed."
  echo "        R script will use self-reported data only (Field 20002)."
fi

# ==========================================================================
#  UPLOAD
# ==========================================================================
echo ""
echo "=== Uploading to project ==="
for f in ukb_selfreport_illness.csv ukb_icd10_secondary.csv; do
  if [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]; then
    dx upload "$f" --destination / --brief
    echo "  Uploaded: $f"
  else
    echo "  Skipped: $f (not available)"
  fi
done

echo ""
echo "============================================"
echo " EXTRACTION COMPLETE"
echo " Download from the Manage tab:"
echo "   1. ukb_selfreport_illness.csv  (REQUIRED)"
echo "   2. ukb_icd10_secondary.csv     (OPTIONAL)"
echo ""
echo " Place in C:/Users/nader/Downloads/"
echo " Then run: 16_sensitivity_secondary_exclusion.R"
echo "============================================"
