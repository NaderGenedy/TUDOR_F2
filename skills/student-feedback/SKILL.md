---
name: student-feedback
description: Use when generating academic feedback for medical students at Cardiff University — forum posts, essays, reflective pieces, OSCE self-evaluations, end-of-module summaries. Triggers on "give feedback to students", "mark these submissions", "feedback for [student name]", "module summary", "moderate the cohort", "3-paragraph academic feedback". Enforces calibrated tone (Cardiff academic register), 3-paragraph structure (strength → development point → action), British English, and dual md+docx output for the entire cohort.
---

# Student-feedback — calibrated 3-paragraph academic register at scale

## When to invoke

Use this skill the moment the user asks for:
- "feedback on these student submissions"
- "mark these forum posts" / "feedback on the essays"
- "feedback for [student name]" / "feedback for the cohort"
- "module summary" / "end-of-module report"
- "give 3-paragraph feedback"
- "moderate the marks"

This skill replaces ad-hoc per-student responses with the canonical Cardiff three-paragraph structure.

## The three-paragraph structure (every student, every time)

**Paragraph 1 — Specific strength.** Open with one concrete strength from the submission, quoting or paraphrasing one phrase. Do not start with "Well done" or "Great work" — start with the substance.

**Paragraph 2 — Development point with reasoning.** Identify one specific area for development. Explain WHY it matters in clinical or academic terms. Reference the relevant module learning outcome where possible.

**Paragraph 3 — Action for next time.** Give one concrete next-step action the student can take in the next submission. Frame as an opportunity, not a deficit. Close with a short forward-looking sentence (one line, no exclamation marks).

**Length:** 180–260 words per student. Anything shorter feels dismissive, anything longer is unread.

## Tone calibration (Cardiff register)

- British English throughout (analyse, behaviour, organisation, colour).
- Second person ("you demonstrate", "your reflection") — never third person.
- One hedge word maximum per paragraph ("perhaps", "might", "could").
- No emojis, no exclamation marks, no "Great job!" / "Awesome!" / "Keep it up!".
- No fawning. Treat the student as a future colleague.
- One clinical or evidence-based reference per feedback if natural — never forced.

## Calibrated grade phrasing (when grades are present)

| Mark band | Opening phrase |
|---|---|
| 80–100 (1st, distinction) | "This is a strong piece of work that demonstrates..." |
| 70–79 (1st) | "You have produced a thoughtful and well-structured response that..." |
| 60–69 (2:1) | "Your submission demonstrates good understanding of..." |
| 50–59 (2:2) | "Your response shows engagement with..." |
| 40–49 (3rd / pass) | "Your submission engages with the question, and..." |
| <40 (fail / referred) | "Your submission attempts the question; however,..." |

## Batch protocol

1. **Read all submissions first.** Get a sense of cohort range before writing any feedback — this calibrates the tone of "above average" vs "typical" comments.
2. **Use TodoWrite.** One todo per student so progress is visible.
3. **Write per-student feedback.** Save each as `feedback_<student_id>.md` AND `feedback_<student_id>.docx` (invoke the `md2docx` skill).
4. **Generate cohort summary.** A single `cohort_summary.md` + `.docx` at the end with: n submitted, mark distribution table, three commonest strengths, three commonest development points, two recommended cohort-level interventions for next module.
5. **Pivot to end-of-module if requested.** If the user says "make this end-of-module" mid-task, switch from per-week files to a single comprehensive document per student covering all weeks — do not re-explain the change, just do it.

## Output structure (per student)

```
D:\...\feedback\
├── feedback_<id_001>.md
├── feedback_<id_001>.docx
├── feedback_<id_002>.md
├── feedback_<id_002>.docx
├── ...
├── cohort_summary.md
└── cohort_summary.docx
```

## Anti-patterns (stop yourself)

- "Well done on submitting on time." — vacuous, delete.
- "Great use of references!" — no exclamation, and explain WHY the references worked.
- Identical opening sentence across two students — every opening must reflect that student's actual submission.
- Marking with the grade in the first sentence — grade goes at the END or in a separate field.
- Using the same development point for >30 % of the cohort — if you find yourself writing it that often, surface it as a cohort-level intervention in `cohort_summary.md` instead.

## Final QC (before declaring done)

- Every student in the input list has both `.md` and `.docx`.
- No two feedback files share an opening sentence.
- `cohort_summary.docx` opens in Word without warnings.
- Total word count ≈ n_students × 220 ± 20%.
