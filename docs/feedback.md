# Feedback and Corrections

This site is a living document. If you spot an error, a typo, a broken link,
or an explanation that isn't clear, please report it — the course improves
when students tell me what doesn't work.

---

## What to report

| Category | Examples |
|---|---|
| **Typo** | Misspelled word, wrong punctuation, inconsistent symbol ($\lambda$ vs. $\Lambda$) |
| **Mathematical error** | Wrong sign, missing factor, incorrect formula, incorrect numeric answer |
| **Broken link** | A chapter link that returns 404, a reference to a textbook chapter that doesn't exist |
| **Rendering problem** | Equation not displaying, table misaligned, code block unreadable on mobile |
| **Clarity issue** | A worked example whose steps don't follow; a definition that is ambiguous; a solution step that skips justification |
| **Tutorial problem** | A problem that is under-specified, has no solution, or contradicts the chapter |
| **Attribution** | A borrowed figure or example that isn't credited |
| **Accessibility** | Missing alt text, poor colour contrast, keyboard-inaccessible interaction |

## What *not* to report here

- **Grade queries** — email the instructor directly (see [Assessment](assessment.md)).
- **Extension requests** — email the instructor with documentation.
- **General questions** — post on the course forum.
- **Requests to add material** — email the instructor; scope changes are
  considered at the start of each semester.

---

## How to report

### Option A — open a GitHub issue (preferred)

1. Visit the repository:
   [github.com/tesfamichael-g/ECEG-DSP-Courses](https://github.com/tesfamichael-g/ECEG-6201-Analytical-optimization-and-Stochastic-Process)
2. Click **Issues → New issue**.
3. Choose the **"Content correction"** template.
4. Fill in:

   - **Page** (e.g. `chapters/03-constrained-optimization` or `solutions/05-spectral-representation`)
   - **Section** (e.g. "Example 3.4" or "Problem 5.3")
   - **What is wrong** (quote the text or paste a screenshot)
   - **What it should be** (your correction or reasoning)
   - **Severity** (typo / minor / major / blocking)

### Option B — email the instructor

For sensitive reports (e.g. a solution that reveals too much detail about an
upcoming assessment), email **with subject line `[ECEG-DSP-Courses feedback]`**
so it doesn't get lost:


---


## Response times

| Severity | Response target | Fix target |
|---|---|---|
| **Blocking** (site broken, wrong formula in a major example) | 24 hours | 48 hours |
| **Major** (wrong numeric answer, unreadable example) | 3 days | 1 week |
| **Minor** (unclear explanation, small typo) | 1 week | next content pass |
| **Enhancement** (new example, expanded discussion) | acknowledgment within 1 week; implemented at the start of the next semester |

Fixes are deployed automatically via GitHub Actions when pushed to the `main`
branch. The **"last updated"** date in the footer of each page reflects the
most recent change.

---

## Acknowledged contributions

Contributors who report errors that lead to a substantive fix are **acknowledged
by name** on the affected page (with permission). This is a small way of
recognising the effort of the students and colleagues who help improve the
course.

If you would rather stay anonymous, just say so when you report — the
acknowledgment is entirely opt-in.

---

## How this page is maintained

The site is built with **MkDocs** + **Material theme**, deployed via **GitHub
Pages**. All content lives in the `docs/` directory of the repository. Corrections
are made by editing the corresponding `.md` file and pushing to `main` — the
GitHub Actions workflow rebuilds and redeploys in under a minute.

The `# file name:` header at the top of each chat exchange is **not** part of
the published site; it's an artefact of how the course materials are drafted
and reviewed.

---

## Thank you

A course is never finished. Every correction you send makes the next cohort's
experience a little better. Thank you for taking the time to write one.
