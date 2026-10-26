# Syllabus — ECEG-6303 Advanced Digital Signal Processing

**MSc in Communication Engineering · Year 1, Semester I · Compulsory Core**

---

## Course identification

| Field | Value |
|---|---|
| **Course code** | ECEG-6303 |
| **Course title** | Advanced Digital Signal Processing |
| **Degree program** | Master of Science in Communication Engineering |
| **Credit hours (ECTS)** | 3 (5) |
| **Semester** | Year 1, Semester I |
| **Status** | Compulsory Core |
| **Pre-requisites** | None (background in signals and systems, probability, and linear algebra strongly recommended) |
| **Instructor** | Dr. Tesfamichael Gebregziabher |
| **Contact** | *(office, email, office hours — see course page)* |

---

## Student workload

| Lecture | Tutorial | Practice / Laboratory | Home study |
|:---:|:---:|:---:|:---:|
| 2 h | – | 3 h | 5 h |

**Total per week:** ~10 hours. **Over 16 teaching weeks:** ~160 hours, consistent with 5 ECTS.

**Note on tutorial hours.** ECEG-6303 has **no scheduled tutorial** (unlike
ECEG-6423). The **3-hour laboratory** each week substitutes for a tutorial:
the lab session begins with a short problem-solving segment covering the
current chapter's worked examples, followed by the scheduled hands-on lab.

---

## Course description

This course extends the classical digital signal processing sequence into
**statistical** and **adaptive** signal processing. It develops the
mathematical foundation for analysing random signals in discrete time, for
estimating their spectra from finite data, for predicting and filtering
them optimally, and for adapting filter coefficients to track changing
statistics.

**Three threads run through the course.**

1. **Statistical.** Random signals, correlations, spectra, and how to
   estimate them from data.
2. **Estimation and prediction.** ML and MMSE criteria, Wiener and Kalman
   filters, linear prediction, Levinson recursion, lattice filters.
3. **Adaptation.** LMS, RLS, and their applications in equalisation, echo
   cancellation, and noise cancellation — plus the multirate tools
   (polyphase, filter banks) that make modern DSP hardware efficient.

The course is the **first of two DSP courses** in the MSc programme. Its
companion, **ECEG-6423 Multimedia Signal Processing**, applies these tools
to images, speech, and video. ECEG-6303 is a prerequisite for the more
advanced parts of ECEG-6423.

---

## Course objectives and competences

### Objectives

On successful completion of this course, the student will be able to:

- Comprehend advanced signal-processing techniques.
- Carry out **spectrum estimation and prediction**.
- Design **adaptive filters**.
- Carry out **multirate digital signal processing**.

### Competences to be acquired

- Apply digital signal-processing algorithms and **design digital signal
  processors** such as adaptive digital filters.

---

## Learning outcomes

On successful completion, students will be able to:

| # | Learning outcome | Covered in | Assessed by |
|:---:|---|---|---|
| **LO1** | Analyse discrete random processes: expectation, variance, covariance, energy, Parseval's theorem, Wiener–Khinchine, PSD, periodogram, sample autocorrelation, Sum Decomposition, spectral factorization | Ch 1 | Test 1, Final |
| **LO2** | Compute the output statistics of an LTI system driven by a WSS process; simulate white noise; simulate coloured processes by filtering white noise | Ch 1 | Lab 1, Final |
| **LO3** | Estimate spectra using nonparametric (periodogram, Bartlett, Welch) and parametric (AR, MA, ARMA) methods; analyse bias and variance of each | Ch 2 | Test 2, Lab 2, Final |
| **LO4** | Solve the Yule–Walker equations using Levinson–Durbin; select model order via AIC/MDL/FPE | Ch 2 | Lab 2, Final |
| **LO5** | Apply the maximum likelihood and MMSE criteria; derive and use the Wiener filter; derive and use the Kalman filter; derive linear predictors | Ch 3 | Test 2, Lab 3, Final |
| **LO6** | Design a Lattice filter and a Levinson-based predictor; derive the prediction-error filter as a whitening filter | Ch 3 | Lab 3, Final |
| **LO7** | Design LMS and NLMS adaptive filters; state convergence conditions, time constants, and misadjustment; apply to equalisation, echo cancellation, noise cancellation | Ch 4 | Test 3, Lab 4, Final |
| **LO8** | Derive and implement RLS (exponentially weighted and sliding window); compare RLS and LMS | Ch 4 | Lab 4, Final |
| **LO9** | Apply decimation, interpolation, and polyphase decomposition; design a rational $L/M$ resampler; state Noble identities | Ch 5 | Test 3, Lab 5, Final |
| **LO10** | Design a two-channel QMF filter bank; state perfect-reconstruction conditions; describe the DWT as iterated filter banks | Ch 5 | Lab 5, Project, Final |

---

## Teaching and learning methods

- **Lectures (2 h/week).** Motivate and derive each topic; work through at
  least one example per method; demonstrate algorithms on the laptop.
- **Laboratory (3 h/week).** Hands-on work in **Python / MATLAB**:
  implementing each chapter's core algorithms on synthetic and real data,
  measuring empirical bias/variance, convergence rates, and steady-state
  behaviour against theory.
- **Project (weeks 10–15).** A self-contained implementation of one algorithm
  family (LMS, RLS, spectrum estimation, or filter banks) applied to a real
  dataset, with report and code.
- **Independent study (5 h/week).** Reading, tutorial problems, lab prep,
  project work.

**Primary language:** Python. MATLAB code is provided alongside for students
whose lab is MATLAB-based.

---

## Chapter-by-chapter content

### Chapter 1 — Discrete Random Signal Processing

**Topics.** Discrete random processes; expectation, variance, covariance;
scalar product and energy of discrete signals; **Parseval's theorem**;
**Wiener–Khinchine relation**; **power spectral density**; **periodogram**;
**sample autocorrelation**; **Sum Decomposition Theorem**; **Spectral
Factorization Theorem**; discrete random signal processing by linear systems;
simulation of white noise; lowpass filtering of white noise.

**Key outcomes.** LO1, LO2. **Delivered week 1–4. Assessed in Test 1 and the
final examination.**

### Chapter 2 — Spectrum Estimation

**Topics.** **Nonparametric methods**: correlation method, covariance
estimator, performance analysis of estimators, unbiased and consistent
estimators, periodogram, **Bartlett** spectrum estimation, **Welch**
estimation. **Model-based approach**: AR, MA, ARMA signal modelling,
parameter estimation using the **Yule–Walker method**, **Levinson–Durbin**
recursion, model-order selection (AIC, MDL, FPE).

**Key outcomes.** LO3, LO4. **Delivered week 5–7. Assessed in Test 2, Lab 2,
and the final.**

### Chapter 3 — Linear Estimation and Prediction

**Topics.** **Maximum likelihood** criterion; efficiency of estimators;
**least-mean-squared-error** criterion; **Wiener filter**; discrete
**Wiener–Hopf** equations; recursive estimators; **Kalman filter**; linear
prediction; prediction error; **whitening filter**; **inverse filter**;
**Levinson recursion**; **lattice realization**; Levinson recursion for
solving Toeplitz systems.

**Key outcomes.** LO5, LO6. **Delivered week 8–10. Assessed in Test 2, Lab 3,
and the final.**

### Chapter 4 — Adaptive Filters

**Topics.** FIR adaptive filters; Newton's steepest-descent method; adaptive
filters based on steepest descent; **Widrow–Hoff LMS adaptive algorithm**;
**adaptive channel equalization**; **adaptive echo canceller**; **adaptive
noise cancellation**; **RLS adaptive filters**; exponentially weighted RLS;
sliding-window RLS; simplified IIR LMS adaptive filter; NLMS.

**Key outcomes.** LO7, LO8. **Delivered week 11–13. Assessed in Test 3, Lab
4, and the final.**

### Chapter 5 — Multirate Digital Signal Processing

**Topics.** Mathematical description of **sampling-rate change**;
**interpolation** and **decimation**; continuous-time model; direct
digital-domain approach; decimation by integer factor; interpolation by
integer factor; single and multistage realization; **polyphase realization**;
**Noble identities**; applications to **subband coding**; **wavelet
transform** and filter-bank implementation of wavelet expansion of signals;
two-channel QMF filter banks.

**Key outcomes.** LO9, LO10. **Delivered week 14–16. Assessed in Test 3, Lab
5, the project, and the final.**

---

## Assessment overview

| Component | Weight | Timing | Structure |
|---|:---:|---|---|
| **Tests (3)** | **20%** | Weeks 4, 8, 14 | 45-min in-class tests, closed book |
| **Project and labs** | **30%** | Project out week 10, due week 15 | Seven lab reports (15%) + project (15%) |
| **Final examination** | **50%** | Week 18 | 3 h, closed book |

See the [Assessment page](assessment.md) for detailed rubrics, submission
rules, and policies.

### Mapping of assessments to learning outcomes

| Assessment | Covers | Weight |
|---|---|:---:|
| **Test 1** (week 4) | LO1, LO2 | 5% |
| **Test 2** (week 8) | LO3, LO4, LO5 | 7% |
| **Test 3** (week 14) | LO7, LO8, LO9 | 8% |
| **Lab 1** (week 2) | LO1, LO2 | ~2% |
| **Lab 2** (week 5) | LO3, LO4 | ~2% |
| **Lab 3** (week 9) | LO5, LO6 | ~2% |
| **Lab 4** (week 12) | LO7, LO8 | ~2% |
| **Lab 5** (week 15) | LO9, LO10 | ~2% |
| **Lab reports (all)** | All | 15% total |
| **Project** | LO5–LO10 (student choice) | 15% |
| **Final exam** | All | 50% |

Per-test weighting 5% / 7% / 8% increases with difficulty. Lab per-report
weighting is roughly equal across the five labs. The project must address at
least two of LO5–LO10.

---

## Textbooks and references

### Textbooks

1. **J. G. Proakis & D. G. Manolakis**, *Digital Signal Processing: Principles,
   Algorithms, and Applications*, Pearson Education / PHI, 2007.
2. **M. H. Hayes**, *Statistical Digital Signal Processing and Modeling*,
   John Wiley & Sons, New York, 1996. *(Primary reference for Chapters 1–3
   and 5.)*
3. **A. V. Oppenheim & R. W. Schafer**, *Discrete-Time Signal Processing*, PHI.
4. **B. Venkataramani & M. Bhaskar**, *Digital Signal Processors — Architecture,
   Programming and Applications*, Tata McGraw Hill, 2002.

### Reference books

1. **A. Antoniou**, *Digital Signal Processing*, Tata McGraw Hill, 2006.
2. **M. H. Hayes**, *Schaum's Outlines: Digital Signal Processing*, Tata McGraw
   Hill, 2007.
3. **C. B. Rorabaugh**, *DSP Primer*, Tata McGraw Hill, 2005.
4. **R. J. Schilling & S. L. Harris**, *Fundamentals of Digital Signal
   Processing using MATLAB*, Thomson, 2007.
5. **A. V. Oppenheim & R. W. Schafer**, *Digital Signal Processing*, PHI, 2006.

Full references, free online versions, and tool recommendations are on the
[Resources](resources.md) page.

---

## Course policy

- All students are expected to abide by the **Code of Conduct for Students**
  in the **Senate Legislation** of the University throughout this course.
- **80% lecture attendance is mandatory.**
- **100% laboratory attendance is mandatory.** Missing more than one lab
  without documented cause results in loss of the lab-report mark for the
  affected session.
- **Academic integrity.** All submitted work must be your own. Plagiarism
  and unauthorised AI-generated work result in a zero on the affected
  component and a report to the Dean's office. See the
  [Assessment](assessment.md) page for the full academic-integrity policy.

---

## Where to find more

| Resource | Link |
|---|---|
| Chapter notes, worked examples, code | [Chapters](chapters/01-discrete-random-signal-processing.md) |
| Solutions to tutorial problems | [Solutions](solutions/index.md) |
| Exam bank (100 marks of practice questions) | [Exam Bank](exam-bank.md) |
| Week-by-week schedule | [Schedule](schedule.md) |
| Detailed assessment rubric | [Assessment](assessment.md) |
| Textbooks, tools, datasets | [Resources](resources.md) |
| Shared DSP notation, transforms, references | [DSP Primer](../shared/dsp-primer.md) |
| Shared software setup | [Software](../shared/software.md) |
| Shared reference shelf | [References](../shared/references.md) |
