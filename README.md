# MSc DSP Courses — ECEG-6303 & ECEG-6423

Source repository for the **MSc Digital Signal Processing course site**:

> **📖 Read the rendered site here →** <https://tesfamichael-g.github.io/ECEG-DSP-Courses/>

The rendered site is the canonical course material — it includes properly
typeset math, worked examples, code, quizzes, solutions, and exam banks. The
`.md` files in `docs/` are the source those pages are built from. They do not
render cleanly on github.com because the MkDocs extensions (MathJax,
admonitions, Material icons, grid cards) are not understood by GitHub's
Markdown renderer.

---

## Courses covered

| Code | Title | Semester | Status |
|---|---|---|---|
| **ECEG-6303** | Advanced Digital Signal Processing | Year 1, Sem. I | Compulsory Core |
| **ECEG-6423** | Multimedia Signal Processing | Year 1, Sem. II | Elective II |

Both are part of the **MSc in Communication Engineering** programme.

Each course has, in the rendered site:

- Lecture notes with worked examples
- Computational implementations in **Python / Java / C#**
- Solutions to the tutorial problems
- A three-part quiz with answer key (MCQ, fill-in, discuss)
- **Exam bank** — 100 marks of exam-style questions with rubric
- Syllabus, schedule, assessment, resources

Shared across both courses:

- **DSP Primer** — notation, transforms, formulas
- **Software setup** — installation and tooling
- **References** — common DSP textbook shelf

---

## Repository structure

```
.
├── mkdocs.yml                    # MkDocs + Material theme configuration
├── requirements.txt              # Python dependencies for local build
├── README.md                     # this file
├── .github/
│   └── workflows/
│       └── deploy.yml            # GitHub Actions: build + deploy to Pages
└── docs/
    ├── index.md                  # landing page (DSP track)
    ├── feedback.md               # correction/feedback channel
    ├── shared/
    │   ├── dsp-primer.md         # notation, transforms, formulas
    │   ├── software.md           # installation, tools, IDEs
    │   └── references.md         # shared textbook shelf
    ├── eceg-6303/
    │   ├── index.md
    │   ├── syllabus.md
    │   ├── schedule.md
    │   ├── assessment.md
    │   ├── resources.md
    │   ├── exam-bank.md
    │   ├── chapters/
    │   │   ├── 01-discrete-random-signal-processing.md
    │   │   ├── 02-spectrum-estimation.md
    │   │   ├── 03-linear-estimation-prediction.md
    │   │   ├── 04-adaptive-filters.md
    │   │   └── 05-multirate-dsp.md
    │   └── solutions/
    │       ├── index.md
    │       └── 01-… 05-….md
    └── eceg-6423/
        ├── index.md
        ├── syllabus.md
        ├── schedule.md
        ├── assessment.md
        ├── resources.md
        ├── exam-bank.md
        ├── chapters/
        │   ├── 01-continuous-discrete-images.md
        │   ├── 02-image-processing-techniques.md
        │   ├── 03-speech-processing-intro.md
        │   ├── 04-speech-processing-methods.md
        │   ├── 05-speech-applications.md
        │   └── 06-video-processing.md
        └── solutions/
            ├── index.md
            └── 01-… 06-….md
```

---

## Local preview

Requires Python 3.10+ and `pip`.

```bash
# Install dependencies
pip install -r requirements.txt

# Build and serve locally
mkdocs serve
```

Open <http://127.0.0.1:8000/> in a browser.

`requirements.txt` contains:

```
mkdocs-material>=9.5
pymdown-extensions>=10.7
```

To build a static version without serving:

```bash
mkdocs build --clean
# Output goes to site/
```

---

## Deployment

Pushing to the `main` branch triggers the GitHub Actions workflow
(`.github/workflows/deploy.yml`), which builds the site and publishes it to
the `gh-pages` branch. **GitHub Pages must be set to serve from the `gh-pages`
branch** (Settings → Pages → Source = "Deploy from a branch: `gh-pages` / root").

If you prefer the newer "GitHub Actions" Pages source instead of the
`gh-pages` branch, update the workflow to use `actions/upload-pages-artifact`
and `actions/deploy-pages`, and set Pages → Source = "GitHub Actions". The
two approaches are exclusive — pick one.

---

## Content conventions

These conventions keep the site consistent and make the source readable when
browsed on github.com.

- **Math**: uses `$...$` for inline and `$$...$$` for display, rendered by MathJax.
  Block math must be on its own lines with blank lines above and below.
- **Admonitions**: `!!! note`, `!!! tip`, `!!! warning` — rendered as coloured
  callouts on the site, as literal `!!!` blocks on github.com.
- **Icons**: `:material-*:` and `:octicons-*:` shortcodes, rendered as SVG icons
  on the site, as literal `:material-*:` text on github.com.
- **Grid cards**: `<div class="grid cards" markdown>...</div>` — rendered as
  Material grid on the site, as raw HTML on github.com.
- **Answer keys**: `??? success "…"` — collapsible blocks in the rendered site,
  always visible on github.com.

**Expect to see rough edges on github.com.** Read the site instead.

---

## Contributing

Corrections, additions, and improvements are welcome. Two channels:

- **Anonymous / any report**: use the feedback form linked from the
  [Feedback](https://tesfamichael-g.github.io/ECEG-DSP-Courses/feedback/) page.
- **GitHub users**: open a
  [content-correction issue](https://github.com/tesfamichael-g/ECEG-DSP-Courses/issues/new?template=content-correction.yml)
  or a pull request.

When contributing, please respect the content conventions above (blank lines
around block math, no spaces inside `$...$`, Material icon names, etc.). A
`mkdocs build --clean` should succeed without warnings before submitting a PR.

---

## Related courses

The MSc in Communication Engineering includes several other courses, each
hosted in its own repository. A **centralized programme hub** linking all
MSc course sites is planned. Until it exists, this README deliberately
links only to the two DSP course sites served from this repository.

---

## Licence

Course content is licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code snippets
are licensed under [MIT](https://opensource.org/licenses/MIT). Individual
third-party figures and datasets retain their own licenses — see
[Resources](https://tesfamichael-g.github.io/ECEG-DSP-Courses/shared/references/)
for attributions.

---

## Contact

**Instructor:** Dr. Tesfamichael Gebregziabher
**Course site:** <https://tesfamichael-g.github.io/ECEG-DSP-Courses/>
