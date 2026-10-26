# ECEG-6303 — Advanced Digital Signal Processing

**MSc in Communication Engineering · Year 1, Semester I**

This is the course website for **Advanced Digital Signal Processing**. It
extends the classical DSP sequence into statistical and adaptive signal
processing: random-signal models, spectrum estimation, Wiener and Kalman
filtering, adaptive filters (LMS, RLS), and multirate techniques.

!!! info "📖 This course is part of the MSc DSP track"
    The companion course [ECEG-6423 — Multimedia Signal Processing](../../ecg-6423/index.md)
    shares the same [DSP Primer](../../shared/dsp-primer.md) and
    [software setup](../../shared/software.md).

---

## Course at a glance

| | |
|---|---|
| **Course code** | ECEG-6303 |
| **Credit hours (ECTS)** | 3 (5) |
| **Semester** | Year 1, Semester I |
| **Status** | Compulsory Core |
| **Pre-requisites** | None |
| **Instructor** | Dr. Tesfamichael Gebregziabher |
| **Contact** | *(office, email, office hours)* |

**Student workload (hours/week)**

| Lecture | Tutorial | Practice / Laboratory | Home study |
|:---:|:---:|:---:|:---:|
| 2 | – | 3 | 5 |

---

## Learning outcomes

By the end of this course you will be able to:

1. **Analyze discrete random processes**: expectations, variance, covariance,
   energy, Parseval's theorem, Wiener–Khinchine, power spectral density,
   periodogram, sample autocorrelation, spectral factorization.
2. **Estimate spectra** using nonparametric (correlation, periodogram, Bartlett,
   Welch) and parametric (AR, MA, ARMA) methods, and evaluate their bias and
   variance.
3. **Design linear estimators and predictors**: Wiener filter, Kalman filter,
   Levinson recursion, lattice realization, whitening filter.
4. **Design adaptive filters**: steepest-descent LMS, Widrow–Hoff algorithm,
   adaptive equalization, echo cancellation, noise cancellation, RLS.
5. **Apply multirate DSP techniques**: decimation, interpolation, polyphase
   realization, filter banks, wavelet transforms.

---

## Chapter overview

<div class="grid cards" markdown>

- :material-numeric-1-circle: **Discrete Random Signal Processing**

    ---

    Discrete random processes, moments, Parseval, Wiener–Khinchine, PSD,
    periodogram, sample autocorrelation, sum decomposition, spectral factorization,
    linear filtering, simulation of white noise, lowpass filtering of white noise.

    [:octicons-arrow-right-24: Chapter 1](chapters/01-discrete-random-signal-processing.md)

- :material-numeric-2-circle: **Spectrum Estimation**

    ---

    Nonparametric methods (correlation, covariance, periodogram, Bartlett, Welch)
    and parametric model-based methods (AR, MA, ARMA, Yule–Walker).

    [:octicons-arrow-right-24: Chapter 2](chapters/02-spectrum-estimation.md)

- :material-numeric-3-circle: **Linear Estimation and Prediction**

    ---

    ML criterion, efficiency, MSE, Wiener filter, Wiener–Hopf equations,
    Kalman filter, linear prediction, Levinson recursion, lattice realization.

    [:octicons-arrow-right-24: Chapter 3](chapters/03-linear-estimation-prediction.md)

- :material-numeric-4-circle: **Adaptive Filters**

    ---

    Steepest descent, LMS (Widrow–Hoff), adaptive equalization, echo cancellation,
    noise cancellation, RLS (exponentially weighted and sliding window).

    [:octicons-arrow-right-24: Chapter 4](chapters/04-adaptive-filters.md)

- :material-numeric-5-circle: **Multirate DSP**

    ---

    Sampling-rate change, decimation, interpolation, polyphase realization,
    filter banks, wavelet transform.

    [:octicons-arrow-right-24: Chapter 5](chapters/05-multirate-dsp.md)

</div>

---

## Assessment

| Component | Weight |
|---|:---:|
| Tests | **20%** |
| Project and labs | **30%** |
| Final examination | **50%** |
| **Total** | **100%** |

See [Assessment](assessment.md) for the detailed rubric and policies.

---

## Recommended preparation

- **Signals and systems** — convolution, z-transform, DFT, sampling theory
- **Probability** — random variables, expectation, variance, distributions
- **Linear algebra** — matrix algebra, eigenvalues, positive-definiteness
- **Programming** — Python (NumPy/SciPy) for the labs; MATLAB acceptable

If any of these feels weak, start with the [DSP Primer](../../shared/dsp-primer.md)
before Chapter 1.

---

## Course policy

- All students are expected to abide by the code of conduct of the University's
  Senate Legislation throughout this course.
- **80% lecture attendance and 100% lab attendance are mandatory.**