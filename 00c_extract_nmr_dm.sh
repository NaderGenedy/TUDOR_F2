#!/bin/bash
# ==============================================================================
# TUDOR PIPELINE: STEP 00c — NMR METABOLOMICS + DIABETES EXTRACTION (v2)
# ==============================================================================
# PURPOSE: Extract NMR lipoprotein subclass data + diabetes/HbA1c fields
#          for TRG Shield biology validation and NMR metabolomics analysis.
#
# FIX (v2): Original 4 batches of 50 fields each hit DataTooLarge API errors.
#           Now uses 5 small batches of 10-11 ESSENTIAL fields each.
#           Only extracts the ~53 NMR fields actually needed by script 09.
#
# RUN ON: UKB-RAP platform (JupyterLab Bash terminal)
#
# OUTPUTS:
#   1. ukb_dm_hba1c.csv  — Diabetes diagnosis + HbA1c (ALREADY DONE)
#   2. ukb_nmr_a.csv     — Subclass total lipids: XXL-VLDL to S-LDL (10 fields)
#   3. ukb_nmr_b.csv     — Subclass lipids: HDL + Subclass TG VLDL (11 fields)
#   4. ukb_nmr_c.csv     — Subclass TG LDL/HDL + Particle conc VLDL (10 fields)
#   5. ukb_nmr_d.csv     — Particle concentrations: M-VLDL to S-HDL (11 fields)
#   6. ukb_nmr_e.csv     — Particle sizes + Clinical NMR lipids (11 fields)
# ==============================================================================

RECORD="project-J6K175jJZ01XppV5477pkYvJ:record-J6K32f8JgZ4JX4gYF39zjBQz"

echo "============================================"
echo " NMR EXTRACTION v2: Small batches (10 fields each)"
echo " Fixes DataTooLarge errors from v1"
echo "============================================"
echo ""

# ==========================================================================
# SKIP DM — already extracted successfully
# ==========================================================================
if [ -f "ukb_dm_hba1c.csv" ]; then
  echo "ukb_dm_hba1c.csv already exists ($(wc -l < ukb_dm_hba1c.csv) lines). Skipping."
else
  echo "=== Extracting Diabetes + HbA1c ==="
  python3 -c "
fields = ['participant.eid', 'participant.p2443_i0', 'participant.p30750_i0']
print('\n'.join(fields))
" > /tmp/dm_fields.txt
  dx extract_dataset "$RECORD" --fields-file /tmp/dm_fields.txt --output ukb_dm_hba1c.csv --delimiter ","
  echo "Done: ukb_dm_hba1c.csv ($(wc -l < ukb_dm_hba1c.csv) lines)"
fi
echo ""

# ==========================================================================
# BATCH A: Subclass total lipids (XXL-VLDL to S-LDL) — 10 fields
# p23400-p23409
# ==========================================================================
echo "=== Batch A: Subclass total lipids (p23400-p23409) ==="
python3 -c "
fields = ['participant.eid']
for f in range(23400, 23410):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_a.txt
echo "  Fields: $(wc -l < /tmp/nmr_a.txt) (including eid)"

dx extract_dataset "$RECORD" --fields-file /tmp/nmr_a.txt --output ukb_nmr_a.csv --delimiter ","

if [ -f "ukb_nmr_a.csv" ] && [ "$(wc -l < ukb_nmr_a.csv)" -gt 1 ]; then
  echo "  OK: ukb_nmr_a.csv ($(wc -l < ukb_nmr_a.csv) lines)"
else
  echo "  FAILED: ukb_nmr_a.csv not created or empty"
fi
echo ""

# ==========================================================================
# BATCH B: Subclass total lipids HDL + Subclass TG VLDL — 11 fields
# p23410-p23413 (HDL subclass lipids) + p23442-p23448 (VLDL TG)
# ==========================================================================
echo "=== Batch B: HDL lipids + VLDL triglycerides (p23410-13, p23442-48) ==="
python3 -c "
fields = ['participant.eid']
for f in range(23410, 23414):
    fields.append(f'participant.p{f}_i0')
for f in range(23442, 23449):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_b.txt
echo "  Fields: $(wc -l < /tmp/nmr_b.txt) (including eid)"

dx extract_dataset "$RECORD" --fields-file /tmp/nmr_b.txt --output ukb_nmr_b.csv --delimiter ","

if [ -f "ukb_nmr_b.csv" ] && [ "$(wc -l < ukb_nmr_b.csv)" -gt 1 ]; then
  echo "  OK: ukb_nmr_b.csv ($(wc -l < ukb_nmr_b.csv) lines)"
else
  echo "  FAILED: ukb_nmr_b.csv not created or empty"
fi
echo ""

# ==========================================================================
# BATCH C: Subclass TG LDL/HDL + Particle conc large VLDL — 10 fields
# p23449-p23455 (LDL/HDL TG) + p23470-p23472 (XXL/XL/L VLDL particles)
# ==========================================================================
echo "=== Batch C: LDL/HDL TG + large VLDL particles (p23449-55, p23470-72) ==="
python3 -c "
fields = ['participant.eid']
for f in range(23449, 23456):
    fields.append(f'participant.p{f}_i0')
for f in range(23470, 23473):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_c.txt
echo "  Fields: $(wc -l < /tmp/nmr_c.txt) (including eid)"

dx extract_dataset "$RECORD" --fields-file /tmp/nmr_c.txt --output ukb_nmr_c.csv --delimiter ","

if [ -f "ukb_nmr_c.csv" ] && [ "$(wc -l < ukb_nmr_c.csv)" -gt 1 ]; then
  echo "  OK: ukb_nmr_c.csv ($(wc -l < ukb_nmr_c.csv) lines)"
else
  echo "  FAILED: ukb_nmr_c.csv not created or empty"
fi
echo ""

# ==========================================================================
# BATCH D: Particle concentrations M-VLDL to S-HDL — 11 fields
# p23473-p23483
# ==========================================================================
echo "=== Batch D: Particle concentrations (p23473-p23483) ==="
python3 -c "
fields = ['participant.eid']
for f in range(23473, 23484):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_d.txt
echo "  Fields: $(wc -l < /tmp/nmr_d.txt) (including eid)"

dx extract_dataset "$RECORD" --fields-file /tmp/nmr_d.txt --output ukb_nmr_d.csv --delimiter ","

if [ -f "ukb_nmr_d.csv" ] && [ "$(wc -l < ukb_nmr_d.csv)" -gt 1 ]; then
  echo "  OK: ukb_nmr_d.csv ($(wc -l < ukb_nmr_d.csv) lines)"
else
  echo "  FAILED: ukb_nmr_d.csv not created or empty"
fi
echo ""

# ==========================================================================
# BATCH E: Particle sizes + Clinical NMR lipids — 11 fields
# p23484 (VLDL size), p23485 (LDL size), p23486 (HDL size)
# p23487 (Total Chol), p23488 (LDL-C), p23489 (HDL-C), p23490 (Total TG)
# p23491 (VLDL-C), p23492 (Remnant Chol), p23493 (ApoB), p23494 (ApoA1)
# ==========================================================================
echo "=== Batch E: Particle sizes + Clinical NMR (p23484-p23494) ==="
python3 -c "
fields = ['participant.eid']
for f in range(23484, 23495):
    fields.append(f'participant.p{f}_i0')
print('\n'.join(fields))
" > /tmp/nmr_e.txt
echo "  Fields: $(wc -l < /tmp/nmr_e.txt) (including eid)"

dx extract_dataset "$RECORD" --fields-file /tmp/nmr_e.txt --output ukb_nmr_e.csv --delimiter ","

if [ -f "ukb_nmr_e.csv" ] && [ "$(wc -l < ukb_nmr_e.csv)" -gt 1 ]; then
  echo "  OK: ukb_nmr_e.csv ($(wc -l < ukb_nmr_e.csv) lines)"
else
  echo "  FAILED: ukb_nmr_e.csv not created or empty"
fi
echo ""

# ==========================================================================
# UPLOAD ALL FILES
# ==========================================================================
echo "=== Uploading to project ==="
FAIL_COUNT=0
for f in ukb_dm_hba1c.csv ukb_nmr_a.csv ukb_nmr_b.csv ukb_nmr_c.csv ukb_nmr_d.csv ukb_nmr_e.csv; do
  if [ -f "$f" ] && [ "$(wc -l < "$f")" -gt 1 ]; then
    dx upload "$f" --destination / --brief
    echo "  Uploaded: $f ($(wc -l < "$f") lines)"
  else
    echo "  FAILED: $f missing or empty"
    FAIL_COUNT=$((FAIL_COUNT + 1))
  fi
done

echo ""
echo "============================================"
if [ "$FAIL_COUNT" -eq 0 ]; then
  echo " ALL NMR EXTRACTIONS SUCCEEDED"
else
  echo " WARNING: $FAIL_COUNT file(s) failed"
fi
echo ""
echo " Download from the Manage tab:"
echo "   1. ukb_dm_hba1c.csv  (DM + HbA1c)"
echo "   2. ukb_nmr_a.csv     (Subclass lipids VLDL-LDL)"
echo "   3. ukb_nmr_b.csv     (HDL lipids + VLDL TG)"
echo "   4. ukb_nmr_c.csv     (LDL/HDL TG + VLDL particles)"
echo "   5. ukb_nmr_d.csv     (Particle concentrations)"
echo "   6. ukb_nmr_e.csv     (Particle sizes + Clinical NMR)"
echo ""
echo " Place all files in C:/Users/nader/Downloads/"
echo " Then run: 09_nmr_metabolomics.R"
echo "============================================"
