# Resources — ECEG-6303

## Core textbooks

| Book | Chapters | Use |
|---|---|---|
| Hayes, *Statistical DSP and Modeling*, Wiley, 1996 | 1–6, 9–11 | The primary reference — matches the syllabus almost exactly |
| Proakis & Manolakis, *DSP: Principles, Algorithms, Applications*, Pearson, 2007 | 2, 10, 12, 13 | Classical DSP background; multirate in Chapter 12 |
| Oppenheim & Schafer, *Discrete-Time Signal Processing*, PHI | 2, 4, 6, 10, 12, 13 | Alternative rigorous treatment; multirate and wavelets in Chapter 13 |
| Antoniou, *Digital Signal Processing*, Tata McGraw Hill, 2006 | 7–9, 11–12 | Adaptive filters |
| Haykin, *Adaptive Filter Theory*, 5th ed., Pearson, 2014 | 1–5, 8–9 | The definitive reference for LMS and RLS |

## Free online materials

- **MIT 6.341/6.342** — *Discrete-Time Signal Processing* (OCW) — for review of
  the deterministic foundations.
- **MIT 6.438** — *Algorithms for Inference* (OCW) — for the Kalman and Wiener
  background.
- **Boyd, Stanford EE263** — *Introduction to Linear Dynamical Systems* —
  excellent for Kalman filters.
- **Hayes's DSP course notes** — free online at Georgia Tech.

## Software

| Language | Libraries | Purpose |
|---|---|---|
| Python | `numpy`, `scipy.signal`, `matplotlib`, `statsmodels.tsa`, `pyroomacoustics` | All labs |
| MATLAB | `filter`, `pwelch`, `ar`, `aryule`, `dsp.LMSFilter`, `dsp.RLSFilter`, `dsp.KalmanFilter` | Reference implementations |
| Java | Apache Commons Math, JTransforms, `dsp-collection` | Optional |
| C# | Math.NET Numerics | Optional |

Install Python:

```bash
pip install numpy scipy matplotlib statsmodels pyroomacoustics
```

## Datasets for project and labs

- **NOIZEUS speech corpus** — for noise-cancellation lab
- **MIT-BIH Arrhythmia Database** — for Kalman-filter lab (ECG)
- **Million Song Dataset / FMA** — for multimedia processing (shared with ECEG-6423)
- **TIMIT** — speech corpus for LPC/Levinson labs