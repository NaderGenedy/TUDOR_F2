---
name: ukb-preflight
description: Use BEFORE launching any UKB / FH / Lp(a) analysis pipeline that loads CSVs, fits Cox models, or writes results. Triggers on phrases like "run the pipeline", "kick off the analysis", "rerun all", "regenerate results", "execute v21", "build the slim", or any plan that touches a 100k+ row UKB extract. Verifies environment (Python 3.12 + required packages + console encoding) AND inventories what data is already on disk to avoid duplicate RAP extractions. Saves you the retry cycle that happens when pyarrow is missing or the master CSV has not been re-extracted.
---

# UKB Preflight — verify environment and disk before any pipeline run

## When to invoke

Use this skill BEFORE any UKB analysis pipeline run, including:
- "run the v21 pipeline" / "rerun all"
- "regenerate the slim" / "rebuild the cohort"
- "kick off the Cox / MVMR / TRIPOD"
- "execute the LOCK_RERUN script"
- any task that will load an existing UKB extract or write a new analytic CSV

**Do NOT** assume the environment is correct. Dr Genedy's setup has hit the same three failures repeatedly:
1. Python 3.14 active where 3.12 is needed for pyarrow / fastparquet
2. Console cp1252 encoding crashing on Unicode arrows / Greek letters
3. Re-extracting from RAP what is already sitting on D:/

The whole skill is ~30 seconds of checks that prevents 10–60 minutes of pipeline retry.

## Five-block preflight protocol

### Block 1 — Python environment

Run this Python one-liner FIRST and report the result back:

```bash
python -c "
import sys, importlib, os
print('python:', sys.version.split()[0], 'OK' if sys.version_info[:2]==(3,12) else 'WARN: not 3.12')
for pkg in ['pandas','numpy','statsmodels','lifelines','sklearn','scipy','pyarrow','docx']:
    try:
        m = importlib.import_module(pkg if pkg != 'docx' else 'docx')
        print(f'  {pkg}: OK ({getattr(m,\"__version__\",\"?\")})')
    except Exception as e:
        print(f'  {pkg}: MISSING — pip install {pkg}')
print('PYTHONIOENCODING:', os.environ.get('PYTHONIOENCODING','(unset — set utf-8 for Windows)'))
"
```

If Python is not 3.12, route to:
```
/c/Users/nader/AppData/Local/Programs/Python/Python312/python.exe
```

### Block 2 — Console encoding sanity

Confirm UTF-8 is active for the session:

```bash
python -c "import sys; sys.stdout.reconfigure(encoding='utf-8'); print('β ρ χ² → ✓')"
```

If this raises `UnicodeEncodeError`, the script that follows MUST add at top:
```python
import sys, os
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
os.environ['PYTHONIOENCODING'] = 'utf-8'
```

### Block 3 — R / Rscript availability (only if pipeline uses R figures)

```bash
"/c/Program Files/R/R-4.5.2/bin/Rscript.exe" --version
```

R is rarely on PATH — call by absolute path. If missing, install or skip R-figure steps with a flag.

### Block 4 — Disk inventory before RAP extraction

For any UKB field the user wants, scan local disk FIRST:

```bash
# Recursive scan for any extract that already contains the target field
grep -lE 'p30790|p131296|p131324' \
  D:/Projects/CALON_AlphaFold_Rebuild/data/*.csv \
  D:/Projects/Lpa_Multilevel/data/*.csv \
  D:/Projects/Lpa/*.csv 2>/dev/null | head
```

Also check the master inventory if it exists:
```
D:/Projects/Lpa_Multilevel/data/UKB_amendment_audit/01_inventory_all_extracts.csv
D:/Projects/Lpa_Multilevel/data/UKB_READY_DATA/INVENTORY.md
```

**Rule**: do NOT request RAP extraction without first confirming the field is not already on disk. Most "missing" UKB fields turn out to be already extracted, just unmerged into the slim.

### Block 5 — Lock-state of any output Word file

If the pipeline writes a `.docx`, check it is not currently open in Word:

```bash
python -c "
import os
target = 'D:/Projects/.../manuscript.docx'
try:
    f = open(target, 'r+b'); f.close(); print('OK — not locked')
except PermissionError:
    print('LOCKED — close Word before running pipeline')
"
```

`python-docx` find-and-replace silently writes 0 substitutions when the file is locked.

## Output deliverable

Produce a single one-screen preflight report and present it BEFORE any pipeline command runs:

```
=== UKB preflight ===
python      : 3.12.5  OK
pyarrow     : 17.0.0  OK
console     : UTF-8 OK
Rscript     : 4.5.2 at  C:/Program Files/R/...
slim cohort : D:/Projects/Lpa_Multilevel/data/v21_FIXED_analytic_slim.csv.gz  (612 MB, 502,505 rows)
locked engine: present (37,478 composite ASCVD events)
docx targets: not locked
RAP needed  : no — all 9 fields already on disk
=== ready to run ===
```

If any line is not green, fix it before launching. Never start a 60-minute pipeline on a half-broken environment.

## Companion files

- `D:/Projects/Lpa_Multilevel/data/UKB_READY_DATA/INVENTORY.md` — canonical Showcase-verified field map
- `D:/Projects/Lpa_Multilevel/data/UKB_amendment_audit/01_inventory_all_extracts.csv` — every extract × every p-field
- `~/.claude/hooks/python-env-check.sh` — already auto-fires on any Bash that runs Python; this skill is the manual / explicit version

## Common failures and how to recover

| Symptom | Likely cause | Fix |
|---|---|---|
| `ModuleNotFoundError: pyarrow` | Python 3.14 active | Switch to 3.12: `/c/Users/nader/AppData/Local/Programs/Python/Python312/python.exe` |
| `UnicodeEncodeError: 'charmap'` mid-run | Windows cp1252 console + Greek/arrow chars | Add `sys.stdout.reconfigure(encoding='utf-8')` at top of script |
| `python-docx` find-replace returns 0 subs | Word file is open | Close Word; or fall back to direct `<w:t>` XML edit |
| RAP extraction takes 4 hours | Already on disk | Check `INVENTORY.md` first |
| Subprocess crash from R script | encoding mismatch | Pass `encoding='utf-8', errors='replace', env={'PYTHONIOENCODING':'utf-8'}` to subprocess |
