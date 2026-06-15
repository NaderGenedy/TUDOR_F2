---
name: stage-files
description: Use when the user asks to "collect", "stage", "gather", "bundle", or "prepare" files for transfer / submission / analysis. Triggers BEFORE building any inventory — collect-then-inventory is the canonical user pattern (per the all-caps redirect "NO I WANT U TO COLLECT THE AVAILABLE FILES SITES READY TO BE COPIED"). Also handles locked-file detection for .docx, hard-link/copy fallback chains, and dedup against existing destinations.
---

# Stage-files — collect first, inventory second

## When to invoke

Use this skill the moment the user says any of:
- "collect [files]" / "stage [files]" / "gather [files]"
- "bundle the submission" / "prepare for transfer"
- "copy these to [dest]"
- "get the data ready"
- "package up the [analysis / manuscript / cohort]"

**Do NOT** start with `ls` / `find` / `find . -name '*' | wc -l` / "let me first inventory what's there". The canonical user pattern is **stage first, inventory second**. Inventory is a side-effect of staging, not a precursor.

## Five-step protocol

### Step 1 — Identify the staging targets

Before any disk operation, confirm three things:
- **Source**: which directories / drives to scan (default if unstated: D: for projects, C: for tools, E: for backups)
- **Destination**: where files should land (ask if unstated; do not default)
- **Scope**: which subset (extension, manifest, modality)

If the user has not specified destination, **ask once** — do not invent one. Example:

> Where should I stage these to — `D:/staging/`, a new `out/` folder, or the existing `<project>/submission/`?

Once confirmed, proceed to step 2 immediately.

### Step 2 — Stage with hard-link → copy fallback chain

Try in order:
1. **Hard-link** (`os.link`) — instantaneous, no extra disk space. Fails across drives or filesystems.
2. **Reflink / cp --reflink** — copy-on-write where supported.
3. **Copy** (`shutil.copy2`) — actual byte copy, preserves metadata.

```python
import os, shutil
from pathlib import Path

def stage_one(src: Path, dst: Path) -> str:
    """Return the method that succeeded: 'link' | 'copy'."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return 'skip-exists'
    try:
        os.link(src, dst)
        return 'link'
    except OSError:
        shutil.copy2(src, dst)
        return 'copy'
```

Track method per file and report at the end (e.g. *"linked 47, copied 12, skipped 3 (already present)"*).

### Step 3 — Detect locked / open files (Word, Excel, etc.)

`.docx` and `.xlsx` files locked by Word / Excel cannot be modified or sometimes copied silently. Check before staging:

```python
def is_locked(path: Path) -> bool:
    """Detect Word/Excel lock files in the same directory."""
    name = path.name
    lockfile = path.parent / f'~${name}'
    return lockfile.exists()
```

If a target file is locked, prompt the user *"`<file>` appears open in Word / Excel — close it before staging?"* — do NOT silently fail.

### Step 4 — Dedup against existing destination

If files already exist at the destination, compare by `sha256` of the first 1 MB (fast) before overwriting. Skip identical files; flag size / hash mismatches as `OVERWRITE_PROPOSED` and ask for confirmation.

```python
import hashlib
def quick_hash(path: Path, n_bytes: int = 1_048_576) -> str:
    return hashlib.sha256(path.read_bytes()[:n_bytes]).hexdigest()
```

### Step 5 — Inventory AS A SIDE-EFFECT of staging

Build the inventory only AFTER files have been staged. Output:
- `staging_manifest.csv` — `src, dst, method, size_bytes, sha256_first_1MB, status`
- A 5–10 line chat summary: total files, total bytes, link count, copy count, skip count, any lock-file warnings.

## Pre-flight discipline (UKB / Lp(a) atlas pattern)

Before commissioning a UK Biobank Research Analysis Platform extraction, **always** check existing local holdings:

```python
def preflight_disk_holdings(target_manifest: list[str]) -> dict:
    """For each requested file, check if it exists on local SSD / external / cloud mount."""
    holdings = {}
    for target in target_manifest:
        for search_dir in [Path('D:/'), Path('E:/'), Path('C:/Users/nader/Downloads/')]:
            matches = list(search_dir.rglob(target))
            if matches:
                holdings[target] = matches[0]
                break
        else:
            holdings[target] = None
    return holdings
```

The Lp(a) consequence atlas exemplar (April 2026) discovered that most "missing" RAP extractions were already on disk, saving days of compute and storage. Check before you fetch.

## Output discipline

For every staging operation, produce:
1. `staging_manifest.csv` — full per-file lineage
2. A chat summary: *"Staged N files (X bytes) to `<dst>`. Linked: A, copied: B, skipped: C. Any locked files: D."*
3. **Do NOT** dump the full manifest to chat — write to file and reference its path.

## Common pitfalls (avoid)

- **Inventorying before staging** — stop, the user wants files moved.
- **Defaulting to a destination** — ask once if not specified.
- **Silent overwrite** — confirm before clobbering existing destination files.
- **Ignoring `~$<filename>` Word lock files** — Word holds an exclusive lock that may prevent copy on some Windows configurations.
- **Hard-link across drives** — fails silently; always wrap in try/except and fall back to copy.
- **`shutil.copytree` with existing dest** — use `dirs_exist_ok=True` (Python 3.8+) or pre-merge.

## Reference exemplar

The 5.55 GB Lp(a) consequence atlas staging (April 2026) used this exact protocol: hard-links attempted first, fell back to copies after two failures (cross-drive), produced a clean manifest, and saved a full RAP re-extraction by detecting that most files already existed locally.
