---
name: md2docx
description: Use whenever the user asks for a Word (.docx) version of any markdown manuscript, reply, feedback file, cover letter, or response document. Triggers on "make a docx", "convert to Word", "give me both md and docx", "produce a Word version", "submission-ready docx", or any deliverable destined for journal portals / Cardiff University / clinical colleagues. Enforces dual-format default (md + docx side-by-side), Cambria/Arial body fonts, robust find-and-replace handling of fragmented runs, and embedded image preservation.
---

# md2docx — markdown to Word, every time, identical content

## When to invoke

Use this skill the moment the user says any of:
- "convert to docx" / "give me a Word version" / "produce both md and docx"
- "submission-ready Word file" / "for the portal"
- "send the cover letter as docx"
- "feedback file in Word"
- "the .md is fine, also export Word"

This is the canonical pattern: **manuscripts, replies, feedback, and cover letters always ship in BOTH .md AND .docx**. The .md is the source of truth; the .docx is what humans open.

## Toolchain

Use `python-docx` (already installed in 3.12). Never use Pandoc — it mis-handles Cambria substitution and silently drops Greek letters (β, ρ, χ²) on Windows.

```python
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
```

## Five-step protocol

1. **Read source.** Open the .md file. Note YAML frontmatter, headings, tables, code blocks, inline `code` spans, **bold**, *italic*, and Greek/Unicode characters.
2. **Set page geometry.** A4, 2.0 cm margins all sides, single-column unless the user specified two-column.
3. **Set defaults.** Body Cambria 11 pt; headings Cambria Bold (H1 16, H2 13, H3 11); tables Arial 9 pt; code blocks Consolas 9 pt with grey background.
4. **Render section by section.** Map heading levels, render lists, render markdown tables to `Document.add_table(...)` with autofit, preserve inline emphasis.
5. **Save BOTH files side-by-side.** If source is `paper1_v6.md`, output is `paper1_v6.docx` in the same directory. Never overwrite the .md.

## The fragmented-run problem (always check this)

After a `python-docx` find-and-replace returns 0 substitutions despite the string being visibly present, the cause is **always** that Word has split the run across multiple `<w:r>` XML elements (typically because of inline italic or font changes). Two fixes, in order of preference:

```python
# Fix A — merge runs in each paragraph before replacing
def merge_runs(paragraph):
    if not paragraph.runs:
        return
    base = paragraph.runs[0]
    for r in paragraph.runs[1:]:
        base.text += r.text
        r.text = ""

# Fix B — fall back to direct XML edit
from docx.oxml.ns import qn
for p in doc.paragraphs:
    for t in p._element.iter(qn('w:t')):
        if 'OLDTEXT' in t.text:
            t.text = t.text.replace('OLDTEXT', 'NEWTEXT')
```

## Greek and Unicode preservation

`python-docx` preserves UTF-8 by default, but **only if the source string is read as UTF-8**. Always:

```python
with open(src_md, 'r', encoding='utf-8') as f:
    md = f.read()
```

Never rely on the system default codec on Windows — it will silently mangle β, ρ, ², ⁻¹.

## Tables — always autofit + repeat header row

```python
table = doc.add_table(rows=1, cols=len(headers))
table.style = 'Light Grid Accent 1'
table.autofit = True
hdr = table.rows[0].cells
for i, h in enumerate(headers):
    hdr[i].text = h
    for run in hdr[i].paragraphs[0].runs:
        run.bold = True
table.rows[0].is_header = True   # repeat on page break
```

## Final QC (always, before declaring done)

- File opens in Word 2016+ without "recovery" prompt.
- Page count is plausible (≈ md word count / 350).
- A spot-check of three Greek/Unicode characters survives round-trip.
- No "Style 'Heading 1' is not defined" warnings on open.

## Output

Print the path of BOTH files written, with sizes:

```
WROTE: D:\...\paper1_v6.md         (24 KB)
WROTE: D:\...\paper1_v6.docx       (78 KB)
```

Never declare done if either file is missing or under 1 KB.
