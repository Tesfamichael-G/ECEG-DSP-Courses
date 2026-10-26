# Chapter 3 — Linear Estimation and Prediction

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** maximum likelihood criterion · efficiency of estimator ·
    least-mean-squared-error criterion · Wiener filter · discrete Wiener–Hopf
    equations · recursive estimators · Kalman filter · linear prediction ·
    prediction error · whitening filter · inverse filter · Levinson recursion ·
    lattice realization · Levinson recursion for Toeplitz systems.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. State and apply the **maximum likelihood (ML)** criterion; derive the ML
   estimator of a scalar parameter in Gaussian noise; compute the
   **Cramér–Rao lower bound (CRLB)** and check whether an estimator is
   **efficient**.
2. State the **least-mean-squared-error (MMSE)** criterion; distinguish the
   **optimal (possibly nonlinear) MMSE** from the **linear MMSE (LMMSE)**.
3. Derive the **orthogonality principle** and use it to solve LMMSE problems
   for **random vectors** and **random processes**.
4. Formulate and solve the **Wiener–Hopf** equations for the **FIR** Wiener
   filter; derive the **causal IIR** Wiener filter via **spectral factorization**.
5. Derive and apply the **discrete Kalman filter** (scalar and vector form);
   interpret the **predict-update** recursion and the **Riccati equation**.
6. Derive the **forward and backward linear prediction** equations; interpret
   **prediction error** and the **prediction-error filter**.
7. Design the **whitening filter** and the **inverse filter** as an application
   of spectral factorization and prediction-error filtering.
8. Apply **Levinson–Durbin** (recap from Chapter 2) to solve the **Toeplitz
   normal equations** of linear prediction; state its $O(p^2)$ complexity.
9. Derive the **lattice filter** from the reflection coefficients and explain
   its numerical and structural advantages (stability, modularity, order
   adaptation).
10. Implement all of the above in **Python, Java, and C#**.

---

## 3.1 Overview — what estimation problems look like

Every estimation problem has three ingredients:

- A **parameter** or **signal** to be estimated (scalar $\theta$, vector
  $\boldsymbol\theta$, or a random process $S[n]$).
- **Observations** $Y$ (a vector of samples, or a full observation record).
- A **criterion** to optimise — what makes one estimate "better" than another.

The three criteria we will develop, in increasing generality:

| Criterion | Optimises | Requires |
|---|---|---|
| **Maximum likelihood** | Probability of observations | Full probabilistic model of $Y \mid \theta$ |
| **MMSE (optimal)** | $\mathbb E[(S - \hat S)^2]$ over *all* estimators | Joint distribution of $(S, Y)$ |
| **LMMSE** | $\mathbb E[(S - \hat S)^2]$ over *linear* estimators | Only second-order statistics (means and covariances) |

The LMMSE criterion is the workhorse of digital signal processing: it requires only
first and second moments, it is linear — hence computationally tractable — and
for **jointly Gaussian** signals it **coincides** with the true MMSE.

The Wiener filter is the LMMSE estimator for **stationary processes** (frequency-domain
solution); the Kalman filter is the LMMSE estimator for **state-space models**
(recursive time-domain solution). They are two views of the same problem.

---

## 3.2 Maximum likelihood and efficiency

### 3.2.1 The ML criterion

Given $N$ observations $Y = (y[0], \dots, y[N-1])^\top$ with joint density
$p_Y(y; \theta)$ indexed by a parameter $\theta$, the **maximum likelihood estimator (MLE)**
is

$$
\boxed{\;\hat\theta_{\text{ML}} = \arg\max_\theta p_Y(Y; \theta)\;}
$$

equivalently, maximising the **log-likelihood** $\ell(\theta) = \ln p_Y(Y; \theta)$:

$$
\hat\theta_{\text{ML}} = \arg\max_\theta\, \ell(\theta).
$$

For an i.i.d. sample, $\ell(\theta) = \sum_n \ln p_y(y[n]; \theta)$.

### 3.2.2 Properties of the MLE (under regularity conditions)

| Property | Statement |
|---|---|
| **Consistency** | $\hat\theta_{\text{ML}} \to \theta$ in probability as $N \to \infty$ |
| **Asymptotic normality** | $\sqrt N (\hat\theta_{\text{ML}} - \theta) \to \mathcal N(0, I^{-1}(\theta))$ |
| **Asymptotic efficiency** | Attains the CRLB asymptotically (smallest possible variance) |
| **Invariance** | If $\hat\theta$ is MLE of $\theta$, then $g(\hat\theta)$ is MLE of $g(\theta)$ |

### 3.2.3 The Cramér–Rao lower bound (CRLB)

For any **unbiased** estimator $\hat\theta$ of a scalar parameter $\theta$,

$$
\boxed{\;\operatorname{Var}[\hat\theta] \ge \frac{1}{I(\theta)},\qquad
I(\theta) = -\mathbb E\Bigl[\frac{\partial^2 \ell(\theta)}{\partial \theta^2}\Bigr] = \mathbb E\Bigl[\Bigl(\frac{\partial \ell}{\partial \theta}\Bigr)^2\Bigr]\;}
$$

where $I(\theta)$ is the **Fisher information**. The bound is **achieved** (equality)
iff the log-likelihood has the exponential-family form
$\ell(\theta) = a(\theta)\hat\theta + b(\hat\theta) + c(\theta) + d$.

An estimator is **efficient** if it achieves the CRLB at finite $N$; the MLE is
asymptotically efficient under regularity.

### Worked example 3.1 — estimating amplitude in Gaussian noise

Observations $y[n] = A + w[n]$, $n = 0, \dots, N-1$, with $w[n] \sim \mathcal N(0, \sigma^2)$
i.i.d. Estimate the constant amplitude $A$.

**Log-likelihood.**

$$
\ell(A) = -\frac{N}{2}\ln(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum_{n=0}^{N-1}(y[n] - A)^2.
$$

**MLE.**

$$
\frac{d\ell}{dA} = \frac{1}{\sigma^2}\sum_n (y[n] - A) = 0
\;\Longrightarrow\;
\boxed{\;\hat A_{\text{ML}} = \frac{1}{N}\sum_{n=0}^{N-1} y[n] = \bar y\;}
$$

the **sample mean** — as expected for a Gaussian.

**CRLB.**

$$
\frac{d^2\ell}{dA^2} = -\frac{N}{\sigma^2} \;\Longrightarrow\; I(A) = \frac{N}{\sigma^2}.
$$

So $\operatorname{Var}[\hat A] \ge \sigma^2/N$.

**Efficiency.** $\operatorname{Var}[\hat A_{\text{ML}}] = \operatorname{Var}[\bar y] = \sigma^2/N$.
The sample mean **attains** the CRLB for every finite $N$ — it is an **efficient**
estimator. This is the finite-$N$ exception; in general only asymptotically efficient.

### Worked example 3.2 — estimating frequency in complex Gaussian noise

Observations $y[n] = A e^{j(\omega_0 n + \theta)} + w[n]$, $w$ complex white
Gaussian with variance $\sigma^2$ per real dimension, estimate $\omega_0$ with
$A$, $\theta$ known.

**Log-likelihood.** The dominant term for large $N$ is proportional to the real
part of the **correlation** between observations and the assumed signal:

$$
\ell(\omega) \approx \frac{2 A}{\sigma^2}\,\Re\Bigl\{e^{-j\theta}\sum_n y[n]e^{-j\omega n}\Bigr\}.
$$

**MLE.** Maximising $\ell(\omega)$ is equivalent to maximising
$\bigl|\sum_n y[n]e^{-j\omega n}\bigr|^2$ — the **periodogram**! The MLE of a
sinusoidal frequency is the **peak of the periodogram**. This connects Chapter 2
(spectrum estimation) and Chapter 3 (ML estimation) into one framework.

**CRLB.** The Fisher information is

$$
I(\omega_0) = \frac{2 A^2}{\sigma^2}\sum_n n^2 \sim \frac{2 A^2 N^3}{3\sigma^2},
$$

so $\operatorname{Var}[\hat\omega_0] \ge \frac{3\sigma^2}{2 A^2 N^3} \propto 1/N^3$ —
**faster** than the usual $1/N$ rate. Frequency estimation is unusually accurate.

The CRLB is only tight in the **high-SNR** regime
($\text{SNR} = A^2 N / \sigma^2 \gg 1$); below the **threshold SNR** (typically
~10 dB for periodogram-frequency estimation), the MLE loses lock and the actual
variance is much larger than the CRLB.

---

## 3.3 Least-mean-squared-error (MMSE) criterion

### 3.3.1 Problem statement

Given observations $Y$ and a target $S$, the **MMSE estimator** is

$$
\boxed{\;\hat S_{\text{MMSE}} = \arg\min_{\hat S(Y)}\, \mathbb E\bigl[(S - \hat S)^2\bigr]\;}
$$

where the minimisation is over **all functions** $\hat S : \mathcal Y \to \mathbb R$.
The optimal solution is the **conditional mean**:

$$
\boxed{\;\hat S_{\text{MMSE}}(Y) = \mathbb E[S \mid Y]\;}
$$

with minimum MSE $\operatorname{Var}[S \mid Y]$ (the conditional variance).

**Consequences.**

- The optimal estimator is generally **nonlinear** (except when $(S, Y)$ are jointly Gaussian).
- It depends on the **full joint distribution** of $(S, Y)$ — much stronger information than the second moments alone.
- In practice, the conditional mean is often hard to compute, motivating the **linear** restriction.

### 3.3.2 Linear MMSE (LMMSE)

Restrict the estimator to **affine** functions of $Y$:

$$
\hat S_{\text{LMMSE}} = a + b^\top Y.
$$

Then the solution is given by the **orthogonality principle**:

!!! abstract "Orthogonality principle"
    The LMMSE estimator $\hat S$ satisfies
    $$
    \boxed{\;\mathbb E\bigl[(S - \hat S)\, Y^\top\bigr] = 0\;}
    $$
    and (by the affine structure) $\mathbb E[S - \hat S] = 0$. The **estimation error
    is orthogonal to the observation space**.

### 3.3.3 Closed form

Let $\mu_S = \mathbb E[S]$, $\mu_Y = \mathbb E[Y]$, $\Sigma_{SY} = \mathbb E[(S - \mu_S)(Y - \mu_Y)^\top]$, $\Sigma_{YY} = \mathbb E[(Y - \mu_Y)(Y - \mu_Y)^\top]$. Then

$$
\boxed{\;\hat S_{\text{LMMSE}} = \mu_S + \Sigma_{SY}\,\Sigma_{YY}^{-1}\,(Y - \mu_Y)\;}
$$

with minimum MSE

$$
\boxed{\;\operatorname{MSE}_{\text{min}} = \Sigma_{SS} - \Sigma_{SY}\,\Sigma_{YY}^{-1}\,\Sigma_{YS}\;}
$$

where $\Sigma_{SS} = \operatorname{Var}[S]$ (scalar) and $\Sigma_{YS} = \Sigma_{SY}^\top$.

### 3.3.4 Gaussian case

!!! abstract "Linear = optimal for Gaussian"
    If $(S, Y)$ are **jointly Gaussian**, then the conditional mean is **affine**
    in $Y$, so the LMMSE estimator **is** the MMSE estimator. This is the
    fundamental reason Kalman filtering is optimal for linear Gaussian systems,
    and why LMMSE is so useful in practice — Gaussianity makes it exact.

### Worked example 3.3 — LMMSE of a random signal in additive noise

Let $S \sim \mathcal N(0, \sigma_S^2)$ and $V \sim \mathcal N(0, \sigma_V^2)$ be
independent, with $Y = S + V$. Estimate $S$ from $Y$.

**Second moments.** $\Sigma_{SY} = \mathbb E[SY] = \mathbb E[S(S+V)] = \sigma_S^2$;
$\Sigma_{YY} = \sigma_S^2 + \sigma_V^2$.

**LMMSE.**

$$
\hat S = 0 + \frac{\sigma_S^2}{\sigma_S^2 + \sigma_V^2}\, Y = \frac{\text{SNR}}{1 + \text{SNR}}\, Y,
\quad \text{SNR} = \sigma_S^2/\sigma_V^2.
$$

**Interpretation.**

- Large SNR ($\sigma_V^2 \to 0$): $\hat S \to Y$ — trust the observation.
- Small SNR ($\sigma_V^2 \to \infty$): $\hat S \to 0$ — trust the prior mean (zero).

**Minimum MSE.**

$$
\operatorname{MSE}_{\text{min}} = \sigma_S^2 - \frac{\sigma_S^4}{\sigma_S^2 + \sigma_V^2} = \frac{\sigma_S^2 \sigma_V^2}{\sigma_S^2 + \sigma_V^2} = \frac{1}{\sigma_S^{-2} + \sigma_V^{-2}}.
$$

Notice the **harmonic sum** form: the posterior precision is the sum of the
prior precision ($1/\sigma_S^2$) and the measurement precision ($1/\sigma_V^2$).
This is exactly the **information form** of the Kalman filter update — a preview
of §3.7.

Since the model is jointly Gaussian, this LMMSE **equals** the MMSE. No
nonlinear estimator can do better.

---

## 3.4 The Wiener filter — non-causal

### 3.4.1 Problem statement

Let $X[n]$ be a WSS desired signal, $Y[n]$ an observed WSS process, both jointly
WSS. We want a **linear time-invariant** filter $h[n]$ such that

$$
\hat X[n] = \sum_{k=-\infty}^{\infty} h[k]\, Y[n-k]
$$

minimises $\mathbb E[(X[n] - \hat X[n])^2]$.

**No causality requirement** — the filter may use future observations. This
is the **non-causal** (or *smoothing*) Wiener filter.

### 3.4.2 Orthogonality in the time domain

By the orthogonality principle applied at every lag $m$:

$$
\mathbb E\bigl[(X[n] - \hat X[n])Y^*[n-m]\bigr] = 0,\quad \forall m.
$$

This gives the **normal equations** (Wiener–Hopf equations in infinite form):

$$
\boxed{\;\sum_{k=-\infty}^{\infty} h[k]\, R_Y[m-k] = R_{XY}[m],\quad \forall m \in \mathbb Z.\;}
$$

### 3.4.3 Solution in the frequency domain

Take the DTFT of both sides. The convolution becomes a product:

$$
H(e^{j\omega})\, S_Y(e^{j\omega}) = S_{XY}(e^{j\omega}),
\qquad
\boxed{\;H_{\text{nc}}(e^{j\omega}) = \frac{S_{XY}(e^{j\omega})}{S_Y(e^{j\omega})}\;}
$$

where $S_{XY}$ is the **cross-power spectral density**, $S_{XY}(e^{j\omega}) = \sum_k R_{XY}[k]e^{-j\omega k}$.

**For the signal-in-noise case** with $Y[n] = X[n] + V[n]$, $X \perp V$:

$$
S_{XY} = S_{XX},\qquad S_Y = S_{XX} + S_{VV}
\;\Longrightarrow\;
\boxed{\;H_{\text{nc}}(e^{j\omega}) = \frac{S_{XX}(e^{j\omega})}{S_{XX}(e^{j\omega}) + S_{VV}(e^{j\omega})}\;}
$$

### 3.4.4 Interpretation

The non-causal Wiener filter shapes the passband according to the **local SNR**:

- Where $S_{XX} \gg S_{VV}$ (high SNR): $H \approx 1$ — pass the signal.
- Where $S_{XX} \ll S_{VV}$ (low SNR): $H \approx 0$ — reject the noise.
- Where $S_{XX} = S_{VV}$: $H \approx 1/2$ — half pass, half reject.

The **minimum MSE** is

$$
\operatorname{MSE}_{\text{min}} = \frac{1}{2\pi}\int_{-\pi}^{\pi} S_{XX}(e^{j\omega})\bigl(1 - H_{\text{nc}}(e^{j\omega})\bigr)\, d\omega = \frac{1}{2\pi}\int_{-\pi}^{\pi} \frac{S_{XX} S_{VV}}{S_{XX} + S_{VV}}\, d\omega.
$$

### Worked example 3.4 — non-causal Wiener filter for AR(1) + white noise

$X[n]$ is AR(1): $X[n] = a X[n-1] + W[n]$, $|a| < 1$, $\sigma_W^2 = 1$.
$V[n]$ is white, $\sigma_V^2 = 0.5$, independent of $X$.

**Spectra.**

$$
S_{XX}(e^{j\omega}) = \frac{1}{|1 - a e^{-j\omega}|^2}
= \frac{1}{1 - 2a\cos\omega + a^2},
\quad S_{VV} = 0.5.
$$

**Wiener filter.**

$$
H_{\text{nc}}(e^{j\omega}) = \frac{1}{1 - 2a\cos\omega + a^2 + 0.5(1 - 2a\cos\omega + a^2)}
= \frac{1}{1.5(1 - 2a\cos\omega + a^2)}.
$$

Not a low-pass — it is a **notch-rejecting** filter that attenuates the frequencies
where $S_{XX}$ is small (high frequencies when $a > 0$).

**MSE reduction.** For $a = 0.7$, numeric integration gives
$\operatorname{MSE}_{\text{min}} / \operatorname{Var}[X] \approx 0.31$ — the
Wiener filter eliminates ~70% of the estimation error, using **future observations**.

The **causal** filter (§3.6) cannot use future observations, so it does
**not** achieve the same MSE. The gap between the non-causal and causal MSEs
quantifies the "cost of real-time".

---

## 3.5 The finite-impulse-response (FIR) Wiener–Hopf equations

### 3.5.1 Setup

In practice, $h[n]$ is finite-length: $h[0], \dots, h[M-1]$. The estimation is

$$
\hat X[n] = \sum_{k=0}^{M-1} h[k]\, Y[n-k],
$$

and the orthogonality principle gives **$M$ linear equations**:

$$
\boxed{\;\sum_{k=0}^{M-1} h[k]\, R_Y[m-k] = R_{XY}[m],\quad m = 0, \dots, M-1.\;}
$$

Matrix form (Toeplitz):

$$
\boxed{\;\mathbf R_Y\, \mathbf h = \mathbf r_{XY}\;}
$$

where $\mathbf R_Y$ is the $M \times M$ Toeplitz matrix with $[\mathbf R_Y]_{ij} = R_Y[i-j]$, and $\mathbf r_{XY} = (R_{XY}[0], \dots, R_{XY}[M-1])^\top$.

The solution is $\mathbf h = \mathbf R_Y^{-1}\mathbf r_{XY}$, computed via
**Levinson–Durbin** if $\mathbf R_Y$ is Toeplitz (which it is, for WSS data).

### Worked example 3.5 — 3-tap Wiener filter

$Y[n] = X[n] + V[n]$, $X \perp V$. Known:

$$
R_X[0] = 4,\quad R_X[1] = 2,\quad R_X[k] = 0\ (|k| \ge 2);\qquad
R_V[0] = 1,\quad R_V[k] = 0\ (k \ne 0).
$$

Design a **3-tap** Wiener filter $h[0], h[1], h[2]$ to estimate $X[n]$ from $Y$.

**Correlations.**

$$
R_Y[k] = R_X[k] + R_V[k] \quad\Rightarrow\quad
R_Y[0] = 5,\ R_Y[\pm 1] = 2,\ R_Y[k] = 0\ (|k| \ge 2).
$$

$$
R_{XY}[m] = \mathbb E[X[n]Y[n-m]] = R_X[m].
$$

**Toeplitz system.**

$$
\begin{bmatrix} R_Y[0] & R_Y[1] & R_Y[2] \\ R_Y[1] & R_Y[0] & R_Y[1] \\ R_Y[2] & R_Y[1] & R_Y[0]\end{bmatrix}\begin{bmatrix} h_0 \\ h_1 \\ h_2\end{bmatrix}
= \begin{bmatrix} R_X[0] \\ R_X[1] \\ R_X[2]\end{bmatrix}
$$

$$
\begin{bmatrix} 5 & 2 & 0 \\ 2 & 5 & 2 \\ 0 & 2 & 5\end{bmatrix}\begin{bmatrix} h_0 \\ h_1 \\ h_2\end{bmatrix} = \begin{bmatrix} 4 \\ 2 \\ 0\end{bmatrix}.
$$

**Solve.** From row 3: $2h_1 + 5h_2 = 0 \Rightarrow h_2 = -0.4h_1$.
Row 1: $5h_0 + 2h_1 = 4 \Rightarrow h_0 = 0.8 - 0.4h_1$.
Row 2: $2h_0 + 5h_1 + 2h_2 = 2 \Rightarrow 2(0.8 - 0.4h_1) + 5h_1 + 2(-0.4h_1) = 2$.

$$
1.6 - 0.8h_1 + 5h_1 - 0.8h_1 = 2 \Rightarrow 1.6 + 3.4h_1 = 2 \Rightarrow h_1 = 0.4/3.4 = 0.1176.
$$

Then $h_0 = 0.8 - 0.4\cdot 0.1176 = 0.7529$, $h_2 = -0.0471$.

**Filter.** $\mathbf h = (0.7529,\ 0.1176,\ -0.0471)$.

**Sanity check.** The centre tap dominates (as expected for a lightly-smoothed
signal), the lags are small, and $h_0 + h_1 + h_2 = 0.8234$ — the filter's DC gain
is less than 1 because it needs to reject noise at DC as well.

**Minimum MSE.** $\sigma_X^2 - \mathbf r_{XY}^\top \mathbf h = 4 - (4\cdot 0.7529 + 2\cdot 0.1176 + 0\cdot (-0.0471)) = 4 - 3.2468 = 0.7532$.

Compare with the "pass-through" MSE $\sigma_V^2 = 1$: the Wiener filter reduces
estimation error by ~25%.

### 3.5.2 Numerical conditioning

The Wiener–Hopf matrix $\mathbf R_Y$ is **Toeplitz** and **positive-definite**
(from $R_Y$ being a valid autocorrelation). Its condition number $\kappa(\mathbf R_Y)$
depends on the SNR and the autocorrelation length:

- **High SNR, highly coloured** (autocorrelation decays slowly): $\kappa$ large — the filter is effectively ill-conditioned; small errors in $R_Y$ estimates produce large filter errors.
- **Low SNR, white**: $\kappa$ near 1 — well-conditioned, robust filter.

**Practical guidance.** Estimate $R_Y$, $R_{XY}$ from data using the **biased**
sample autocorrelation (positive-definiteness guaranteed by Chapter 1 §1.4.3);
solve with **Levinson–Durbin** for numerical stability. Avoid naive Gaussian
elimination when $\kappa$ exceeds $10^4$.

---

## 3.6 The causal Wiener filter

### 3.6.1 Problem statement

Same as the non-causal case, but with the **causality** constraint
$h[n] = 0$ for $n < 0$. This is the **real-time** filter; it cannot use
future observations.

The orthogonality principle still holds, but **only for non-negative lags**:

$$
\mathbb E\bigl[(X[n] - \hat X[n])Y^*[n-m]\bigr] = 0,\quad m \ge 0.
$$

The **Wiener–Hopf equations** are now **one-sided**:

$$
\boxed{\;\sum_{k=0}^{\infty} h[k]\, R_Y[m-k] = R_{XY}[m],\quad m \ge 0.\;}
$$

This is the **discrete Wiener–Hopf equation** — harder than the two-sided (non-causal)
case because a DTFT of both sides introduces unknown terms for $m < 0$.

### 3.6.2 Solution via spectral factorization

The classical solution uses **spectral factorization** (Chapter 1 §1.8) plus a
causal-part operation. Step-by-step:

**Step 1 — Factor $S_Y$.** Find the **causal, monic, minimum-phase** $L(z)$ and
$\sigma_e^2$ such that

$$
S_Y(z) = \sigma_e^2\, L(z)\, L(z^{-1}).
$$

**Step 2 — Compute the causal part.** Form

$$
G(z) = \frac{S_{XY}(z)}{\sigma_e^2\, L(z^{-1})},
$$

and take the **causal part** $[G(z)]_+$ — the part of the Laurent expansion in
non-negative powers of $z^{-1}$ (i.e. terms $z^0, z^{-1}, z^{-2}, \dots$).

**Step 3 — Compose the filter.**

$$
\boxed{\;H_{\text{causal}}(z) = \frac{[G(z)]_+}{L(z)}\;}
$$

The numerator is realisable (causal) and the denominator is causal-minimum-phase,
so $H_{\text{causal}}$ is **stable and causal**.

### 3.6.3 The "smoothing" correction

An equivalent statement: the causal Wiener filter is

$$
H_{\text{causal}}(z) = \frac{1}{\sigma_e^2 L(z)}\Bigl[\frac{S_{XY}(z)}{L(z^{-1})}\Bigr]_+.
$$

The operation $[\cdot]_+$ is what distinguishes causal from non-causal:
without it, $H = S_{XY}/(\sigma_e^2 L L^{-1}) = S_{XY}/S_Y$ — the non-causal filter
of §3.4.

**Interpretation.** The causal filter is the non-causal filter **projected onto
the causal subspace** in a particular weighted inner product (weighted by the
whitening filter $1/L$). This projection is what makes the causal filter optimal
*among causal filters*.

### Worked example 3.6 — causal Wiener filter for AR(1) + white noise

Same setup as Example 3.4: $X$ AR(1) with $a = 0.7$, $\sigma_W = 1$, $V$ white
$\sigma_V = 0.5$.

**Step 1 — Factor $S_Y$.** From Example 3.4:

$$
S_Y(e^{j\omega}) = S_{XX} + S_{VV} = \frac{1 + 0.5(1 - 2a\cos\omega + a^2)}{1 - 2a\cos\omega + a^2}.
$$

Write as $S_Y(z) = \sigma_e^2 L(z)L(z^{-1})$ with $L(z) = (1 + bz^{-1})/(1 - az^{-1})$
by factoring the numerator and denominator polynomials in $z$. Matching
coefficients gives $b \approx 0.2648$ (a small positive root), and
$\sigma_e^2 \approx 0.5$.

**Step 2 — Causal part.** $S_{XY}(z) = S_{XX}(z) = \frac{1}{(1-az^{-1})(1-az)}$. Then

$$
G(z) = \frac{S_{XX}(z)}{\sigma_e^2 L(z^{-1})} = \frac{1}{\sigma_e^2(1-az)} \cdot \frac{1-az}{1+bz} \cdot \frac{1}{1-az^{-1}} = \frac{1}{\sigma_e^2}\cdot\frac{1}{(1+bz)(1-az^{-1})}.
$$

Partial-fraction expansion with causal-part extraction gives

$$
[G(z)]_+ = \frac{A}{1-az^{-1}},\quad A = \frac{1}{\sigma_e^2(1+ab)}.
$$

**Step 3 — Compose.**

$$
H_{\text{causal}}(z) = \frac{A/(1-az^{-1})}{(1+bz^{-1})/(1-az^{-1})} = \frac{A}{1+bz^{-1}}.
$$

**Result.** The causal Wiener filter is a **first-order IIR filter** with pole at
$-b \approx -0.2648$ and gain $A$. Compare with the non-causal filter, which is
not rational — a direct consequence of the causality constraint.

**MSE penalty.** Numerically, the causal MSE is $\approx 0.42$ versus the non-causal
$\approx 0.31$ — about **35% more error**, the price of causality.

### 3.6.4 Comparison table

| Filter | Formula | Requires | Optimal among |
|---|---|---|---|
| Non-causal Wiener | $S_{XY}/S_Y$ | Full PSDs | **All** linear filters |
| Causal Wiener | $[S_{XY}/(\sigma_e^2 L^{-1})]_+ / L$ | PSDs + spectral factorization | **Causal** linear filters |
| Kalman (§3.7) | Recursive; steady state equals causal Wiener for stationary $S_Y$ | State-space model | **Causal** linear filters (equivalently) |

Kalman and causal Wiener give the **same** filter in steady state when the process
is stationary and both are causal. The Kalman formulation has computational
advantages (recursive, no large-matrix inversion, handles non-stationarity and
multi-channel observations naturally).

---

## 3.7 The Kalman filter

### 3.7.1 State-space model

The Kalman filter solves the LMMSE problem for **state-space models**:

$$
\boxed{\;x[n] = A x[n-1] + B u[n]\;}\quad\text{(state equation)}
$$

$$
\boxed{\;y[n] = C x[n] + v[n]\;}\quad\text{(observation equation)}
$$

with:

- $x[n] \in \mathbb R^d$ — state vector (unknown, to be estimated).
- $y[n] \in \mathbb R^m$ — observation vector.
- $u[n] \sim \mathcal N(0, Q)$ — **process noise** (white, independent across $n$).
- $v[n] \sim \mathcal N(0, R)$ — **measurement noise** (white, independent of $u$ and $x[0]$).
- $x[0] \sim \mathcal N(\hat x_{0|0}, P_{0|0})$ — initial state prior.

The matrices $A, B, C, Q, R$ are known. The goal is the recursive LMMSE estimate
$\hat x_{n|k} = \mathbb E[x[n] \mid y[1], \dots, y[k]]$.

Two cases:
- **Filtering** ($k = n$): estimate the current state.
- **Prediction** ($k < n$): forecast future states (used in control and tracking).
- **Smoothing** ($k > n$): estimate past states given all future observations (used in offline estimation).

The **Kalman filter** is the **filtering** case; it is a **predict-update** recursion.

### 3.7.2 The recursion

Initialise $\hat x_{0|0}$, $P_{0|0}$. For $n = 1, 2, \dots$:

**Predict step** (extrapolate the state and its covariance):

$$
\boxed{\;\hat x_{n|n-1} = A\,\hat x_{n-1|n-1}\;}
$$

$$
\boxed{\;P_{n|n-1} = A P_{n-1|n-1} A^\top + B Q B^\top\;}
$$

**Kalman gain** (how much to trust the observation):

$$
\boxed{\;K_n = P_{n|n-1} C^\top \bigl(C P_{n|n-1} C^\top + R\bigr)^{-1}\;}
$$

**Update step** (correct the prediction using the observation):

$$
\boxed{\;\hat x_{n|n} = \hat x_{n|n-1} + K_n \bigl(y[n] - C \hat x_{n|n-1}\bigr)\;}
$$

$$
\boxed{\;P_{n|n} = (I - K_n C) P_{n|n-1}\;}
$$

The term $\nu_n = y[n] - C\hat x_{n|n-1}$ is the **innovation** (or prediction
error) — the part of the observation not predicted by the model. If the model is
correct, $\nu_n$ is white with covariance $C P_{n|n-1} C^\top + R$.

### 3.7.3 Interpretation

| Quantity | Meaning |
|---|---|
| $\hat x_{n\|n-1}$ | Prior mean — belief about $x[n]$ **before** seeing $y[n]$ |
| $P_{n\|n-1}$ | Prior covariance — uncertainty in the prior |
| $K_n$ | **Weight** on the innovation — the optimal compromise between prior and observation |
| $y[n] - C\hat x_{n\|n-1}$ | Innovation — new information in $y[n]$ |
| $\hat x_{n\|n}$ | Posterior mean — belief **after** seeing $y[n]$ |
| $P_{n\|n}$ | Posterior covariance — reduced uncertainty |

**Limit cases of the Kalman gain.**

- $R \to 0$ (noiseless observations): $K_n \to C^{-1}$, so $\hat x_{n|n} \to C^{-1}y[n]$ — trust the observation completely.
- $R \to \infty$ (useless observations): $K_n \to 0$, so $\hat x_{n|n} \to \hat x_{n|n-1}$ — ignore the observation.
- $P_{n|n-1} \to 0$ (perfect prior): $K_n \to 0$ — trust the prior.

The Kalman filter **interpolates** between these extremes according to the
**relative uncertainty**.

### Worked example 3.7 — scalar Kalman filter, two steps

Model: $x[n] = 0.9\, x[n-1] + u[n]$, $u[n] \sim \mathcal N(0, 0.1)$;
$y[n] = x[n] + v[n]$, $v[n] \sim \mathcal N(0, 1)$. Initial: $\hat x_{0|0} = 0$, $P_{0|0} = 1$.
Observations: $y[1] = 2$, $y[2] = 1.5$.

**Step 1 — predict.**

$$
\hat x_{1|0} = 0.9 \cdot 0 = 0,\qquad P_{1|0} = 0.9^2 \cdot 1 + 0.1 = 0.91.
$$

**Kalman gain.**

$$
K_1 = P_{1|0}(P_{1|0} + R)^{-1} = 0.91/(0.91 + 1) = 0.4764.
$$

**Update.**

$$
\hat x_{1|1} = 0 + 0.4764\,(2 - 0) = 0.9529,\qquad
P_{1|1} = (1 - 0.4764) \cdot 0.91 = 0.4764.
$$

**Step 2 — predict.**

$$
\hat x_{2|1} = 0.9 \cdot 0.9529 = 0.8576,\qquad P_{2|1} = 0.9^2 \cdot 0.4764 + 0.1 = 0.486.
$$

**Kalman gain.**

$$
K_2 = 0.486/(0.486 + 1) = 0.327.
$$

**Update.**

$$
\hat x_{2|2} = 0.8576 + 0.327\,(1.5 - 0.8576) = 1.067,\qquad P_{2|2} = (1 - 0.327)\cdot 0.486 = 0.327.
$$

**Interpretation.** The first observation ($y[1] = 2$) pulls the estimate up to
$0.95$, but $K_1 < 1$ means we don't fully trust it (noise variance = 1, comparable
to prior variance = 0.91). The second observation ($y[2] = 1.5$) is *lower* than
predicted ($0.858$), so it nudges the estimate to $1.067$. The posterior variance
shrinks from $1 \to 0.476 \to 0.327$ — uncertainty **decreases** with more data.

### Worked example 3.8 — steady-state behaviour

Repeat the same scalar model for many steps. The Riccati recursion
$P_{n|n} = f(P_{n-1|n-1})$ converges to a **fixed point** $P^\star$ satisfying

$$
P^\star = \frac{(A^2 P^\star + Q) R}{A^2 P^\star + Q + R}.
$$

For $A = 0.9$, $Q = 0.1$, $R = 1$:

$$
P^\star = \frac{(0.81 P^\star + 0.1) \cdot 1}{0.81 P^\star + 0.1 + 1}
\;\Rightarrow\;
P^\star \approx 0.327.
$$

The steady-state Kalman gain is $K^\star = P^\star/(P^\star + 1) \approx 0.246$.

**After convergence**, the Kalman filter is a **linear time-invariant filter**:

$$
\hat x_{n|n} = A(1 - K^\star)\hat x_{n-1|n-1} + K^\star y[n] \quad\text{(scalar case)}
= 0.9\cdot 0.754\,\hat x_{n-1|n-1} + 0.246\, y[n].
$$

This LTI recursion is exactly what the **causal Wiener filter** of §3.6 produces
for the same process. The Kalman filter is the **state-space implementation** of
the causal Wiener filter, and it converges to it exponentially fast.

### 3.7.4 The Riccati equation

The **discrete algebraic Riccati equation (DARE)** for the steady-state covariance is

$$
\boxed{\;P = A P A^\top + B Q B^\top - A P C^\top (C P C^\top + R)^{-1} C P A^\top.\;}
$$

**Existence and uniqueness.** Under **stabilisability** (unstable modes of $A$
are affected by the process noise) and **detectability** (unstable modes are
observable through $C$), the DARE has a **unique positive-definite solution**
$P^\star > 0$. Under these conditions, the Kalman filter converges to the
steady-state gain $K^\star$ from any initial $P_{0|0} \ge 0$.

**Solving the DARE in practice.** Iterating the Riccati recursion (the Kalman
recursion for $P$ alone) converges in $O(\log(1/\epsilon))$ iterations —
quadratic convergence for small $Q/R$. For very large $d$, use the **structure-preserving
methods** (Potter, Joseph form, square-root filters) to avoid numerical loss of
positive-definiteness.

### 3.7.5 Why Kalman is more general than Wiener

The Wiener filter assumes:
- **Stationarity** of all processes involved.
- **Known** second-order statistics.
- **Single-channel** observations.
- **Time-invariant** filter (in steady state).

The Kalman filter handles:
- **Non-stationary** processes (time-varying $A, C, Q, R$).
- **Multi-channel** observations with cross-correlated noise.
- **Unknown initial conditions** — the covariance recursion handles it naturally.
- **Real-time** recursive computation — $O(d^3)$ per step, not a large-matrix inversion.

When the Kalman filter is used for a stationary process, **it converges to the causal Wiener filter** in steady state.

---

## 3.8 Linear prediction

### 3.8.1 The one-step forward predictor

Predict $x[n]$ from past samples $x[n-1], \dots, x[n-p]$ using a linear predictor:

$$
\boxed{\;\hat x[n] = -\sum_{k=1}^{p} a_k\, x[n-k]\;}
$$

The **prediction error** (also called the **residual** or **innovation**) is

$$
\boxed{\;e[n] = x[n] - \hat x[n] = x[n] + \sum_{k=1}^{p} a_k x[n-k]\;}
$$

and the **predictor coefficients** $a_k$ are chosen to **minimise the mean-square
prediction error** $\mathbb E[e[n]^2]$.

### 3.8.2 Normal equations

Apply the orthogonality principle to $e[n]$ and each past sample:

$$
\mathbb E\bigl[e[n] x[n-k]\bigr] = 0,\quad k = 1, \dots, p.
$$

Expanding:

$$
R_X[k] + \sum_{m=1}^{p} a_m R_X[k-m] = 0,\quad k = 1, \dots, p.
$$

These are **identical** to the Yule–Walker equations from Chapter 2 (§2.10)! The
**optimal linear predictor** for a WSS process is the **AR coefficient vector** of
the process. The **prediction error** is exactly the **innovations**.

$$
\boxed{\;\text{MMSE linear prediction} \;\equiv\; \text{AR}(p)\ \text{model estimation}.\;}
$$

This identity is fundamental: it says that **the Wold decomposition** of a WSS
process into AR$(\infty)$ + innovations **is** the optimal linear predictor.

### 3.8.3 Minimum prediction error

The MMSE is

$$
\boxed{\;P_p = \mathbb E[e^2[n]] = R_X[0] + \sum_{k=1}^{p} a_k R_X[-k] = \sigma_{w,p}^2.\;}
$$

This is the **prediction-error variance** — the same quantity as $\sigma_p^2$
in Levinson–Durbin (§2.10).

For white noise, $R_X[k] = 0$ for $k \ne 0$, so all $a_k = 0$ and $P_p = R_X[0]$:
no predictability.

For a perfectly predictable signal (e.g. a sinusoid at a known frequency), the
prediction error can go to zero with sufficient $p$ — a **line spectral** process.

### Worked example 3.9 — one-step predictor for AR(1)

$x[n] = a x[n-1] + w[n]$, $|a| < 1$, $\sigma_w^2$ given. $R_X[k] = \frac{\sigma_w^2}{1-a^2}a^{|k|}$.

**AR(1) prediction:** $\hat x[n] = -a_1 x[n-1]$ with $a_1 = -a$ (Yule–Walker for
$p=1$: $R_X[1] + a_1 R_X[0] = 0 \Rightarrow a_1 = -a$).

**MMSE:** $P_1 = R_X[0](1 - \rho^2)$ with $\rho = R_X[1]/R_X[0] = a$. So
$P_1 = \sigma_w^2/(1-a^2) \cdot (1 - a^2) = \sigma_w^2$.

The prediction error is exactly $w[n]$ — the innovations of the AR(1). This
matches the identity above.

**Prediction gain.** $P_1/\sigma_X^2 = 1 - a^2$. For $a = 0.9$, the prediction
error is only 19% of the signal variance; 81% of the signal is predictable
from the previous sample. For $a = 0.5$, only 25% is predictable.

### 3.8.4 Order-$p$ prediction and the lattice

For order $p > 1$, the predictor is obtained via Levinson–Durbin (Chapter 2
§2.10). The intermediate coefficients $a^{(m)}$ at each order give the
**order-$m$ predictor**; the reflection coefficients $k_m$ are the **partial
autocorrelations** and satisfy $|k_m| < 1$ for a stable process.

The **minimum prediction error** at order $m$ is

$$
\boxed{\;P_m = P_{m-1}(1 - k_m^2)\;}
$$

with $P_0 = R_X[0]$. This monotone decreasing sequence converges to $\sigma_w^2$
(the innovations variance) as $m \to \infty$. The **reflection coefficient**
$k_m$ measures the *fractional* improvement in prediction error gained by adding
lag $m$ to the predictor. A magnitude near 1 means a large improvement; near 0
means lag $m$ carries little new information beyond lags $1$ to $m-1$.

### 3.8.5 Backward prediction

Similarly, predict $x[n-p]$ from the *future* samples $x[n-p+1], \dots, x[n]$:

$$
\hat x[n-p] = -\sum_{k=1}^{p} c_k\, x[n-p+k].
$$

For a WSS process, the backward predictor coefficients $c_k$ equal the **reversed**
forward predictor coefficients: $c_k = a_{p+1-k}$. This symmetry is the structural
basis of the lattice filter (§3.11).

---

## 3.9 Whitening filter and inverse filter

### 3.9.1 The whitening filter

Given a WSS process $X[n]$ with rational PSD $S_X(e^{j\omega})$, the **whitening
filter** is

$$
\boxed{\;H_w(e^{j\omega}) = \frac{1}{H(e^{j\omega})},\qquad
S_X = \sigma_e^2 |H|^2.\;}
$$

Its output $E[n] = h_w * X[n]$ is **white** of variance $\sigma_e^2$.

**Interpretation.** The whitening filter is the **prediction-error filter**: it
subtracts the best linear prediction of $X[n]$ from the past, leaving only the
unpredictable part — the innovations.

**For AR($p$):** the whitening filter is the **FIR** filter $1 + \sum_{k=1}^p a_k z^{-k}$,
so

$$
e[n] = x[n] + \sum_{k=1}^p a_k x[n-k].
$$

This is exactly the AR prediction error — the input to the AR model. The
whitening filter for an AR process is **FIR** because AR is an all-pole model
(its inverse is all-zero).

**For MA($q$):** the whitening filter is IIR (the inverse of an all-zero filter),
with coefficients determined by spectral factorization.

**For ARMA($p, q$):** the whitening filter is IIR (ratio of two FIR filters) and
requires factoring the MA polynomial.

### 3.9.2 The inverse filter

The **inverse filter** restores the original process from the white innovations:

$$
X[n] = h * E[n],\qquad h = 1/h_w.
$$

For an AR process, $H(z) = 1/(1 + \sum a_k z^{-k})$ — an IIR filter with feedback.
The pair $(H, H_w)$ is a **synthesis/analysis** pair, used throughout
communications and audio:

- **Synthesis:** generate $X$ from white noise $E$ via $H$ (e.g., speech synthesis).
- **Analysis:** recover $E$ from $X$ via $H_w$ (e.g., speech recognition, LPC).

### Worked example 3.10 — speech analysis and synthesis

The **source-filter model** of speech:

- **Source:** the glottal pulse train (voiced) or turbulent noise (unvoiced) — modeled as white noise.
- **Filter:** the vocal tract — modeled as an AR($p$) filter with $p \approx 10$ for $f_s = 8$ kHz.

**Analysis (encoder):** apply the whitening filter $1/A(z)$ to the speech
$x[n]$, obtain the prediction-error (residual) $e[n]$. The residuals contain
the pitch information and the gain; the AR coefficients $a_k$ contain the
spectral envelope (formants).

**Synthesis (decoder):** drive the AR filter $1/A(z)$ with the residual $e[n]$
to reconstruct $\tilde x[n]$.

**Compression.** For $p = 10$ at $f_s = 8$ kHz and 20 ms frames, the
parameters per frame are: $10$ AR coefficients (~30 bits) + pitch (~7 bits) + gain
(~5 bits) + residual (~0–160 bits depending on codec) — total ~50–200 bits per
20 ms ≈ 2.5–10 kbps. Compare to raw PCM at 128 kbps — a **10–50× compression**,
the basis of **LPC-10**, **CELP**, **MP3**, and modern speech codecs.

The identity **prediction error = innovations** is what makes this possible: the
residual is white and can be efficiently encoded with scalar quantisation or
vector quantisation.

---

## 3.10 Levinson recursion for Toeplitz systems

### 3.10.1 Recap

From Chapter 2 §2.10, the **Levinson–Durbin** recursion solves the Toeplitz system

$$
\mathbf R_p\, \mathbf a = -\mathbf r_p
$$

in $O(p^2)$ operations, producing the AR coefficients at each order $m = 1, \dots, p$:

$$
k_m = -\frac{r_m + \sum_{i=1}^{m-1} a_i^{(m-1)} r_{m-i}}{\sigma_{m-1}^2},\qquad
a_i^{(m)} = a_i^{(m-1)} + k_m a_{m-i}^{(m-1)},\qquad
\sigma_m^2 = (1-k_m^2)\sigma_{m-1}^2.
$$

### 3.10.2 Interpretation

Levinson–Durbin builds the **predictor of order $p$** incrementally:

1. Start with order 1: $a_1^{(1)} = -r_1/r_0$ — the single-lag predictor.
2. Add lag 2: update $a_1^{(1)} \to a_1^{(2)}$ and add $a_2^{(2)}$ — the two-lag predictor.
3. Continue to order $p$.

**Each step adds one lag and updates the previous coefficients via the reflection
coefficient $k_m$.** This is not just a numerical trick — it exposes the
**hierarchical structure** of the predictor, which the lattice filter (§3.11)
exploits directly.

### 3.10.3 Split-Levinson and other variants

For very long $p$ (e.g., $p > 1000$ in audio), the standard recursion is $O(p^2)$.
The **split-Levinson** algorithm reduces to $O(p \log^2 p)$ and is used in
real-time high-order predictors. For most signal-processing applications
($p \le 50$), standard Levinson is fast enough.

**Bareiss and Schur algorithms** are alternative $O(p^2)$ Toeplitz solvers with
better numerical stability; **Levinson is exact in exact arithmetic**, but the
recursion can lose accuracy for ill-conditioned matrices (large condition number,
close to singular). In practice, the **biased** sample autocorrelation guarantees
positive-definiteness, so the recursion is stable.

### 3.10.4 Complexity comparison

| Method | Complexity | Notes |
|---|---|---|
| Gaussian elimination | $O(p^3)$ | General matrices; does not exploit Toeplitz |
| Levinson–Durbin | $O(p^2)$ | Toeplitz, symmetric, positive-definite |
| Bareiss | $O(p^2)$ | Better numerical stability |
| Split-Levinson | $O(p \log^2 p)$ | Very long $p$ |
| FFT-based | $O(p \log p)$ | Circulant embedding (used in practice) |

For $p = 50$: Gaussian $= 125{,}000$ ops, Levinson $= 2{,}500$ ops — a 50×
speedup. For $p = 500$: $5 \times 10^7$ vs $5 \times 10^5$ — a 100× speedup.
For real-time speech coding at $p = 10$, either is fine; for audio at $p = 1000$,
the difference is the feasibility boundary.

---

## 3.11 Lattice realization

### 3.11.1 The lattice structure

The **lattice filter** implements the prediction-error filter using the
**reflection coefficients** $k_1, \dots, k_p$ as its parameters, in a cascade
of **identical 2-port stages**:

$$
\boxed{\;f_m[n] = f_{m-1}[n] + k_m\, g_{m-1}[n-1],\qquad
g_m[n] = g_{m-1}[n-1] + k_m\, f_{m-1}[n]\;}
$$

with initial conditions $f_0[n] = g_0[n] = x[n]$.

**Here:**

- $f_m[n]$ = **forward prediction error** at order $m$ — the error when predicting $x[n]$ from $m$ past samples.
- $g_m[n]$ = **backward prediction error** at order $m$ — the error when predicting $x[n-m]$ from the next $m$ samples.

The forward error at order $p$ is the output: $e[n] = f_p[n]$.

### 3.11.2 Properties

!!! tip "Why lattice filters are preferred in some applications"
    **1. Stability is guaranteed by $|k_m| < 1$.** Any set of reflection coefficients
    within $(-1, 1)$ produces a stable synthesis filter. This is not true of the
    direct-form $a_k$ representation, where stability depends on all roots of
    $1 + \sum a_k z^{-k}$ being inside the unit circle.

    **2. Modular and order-recursive.** Each stage is identical and can be added
    or removed without recomputing the others. This is useful for **order-adaptive
    filtering** (choose $p$ at runtime).

    **3. Numerically robust under quantization.** Coefficients are bounded in
    magnitude, so fixed-point implementations do not overflow or lose stability.

    **4. Provides both forward and backward errors.** Useful in two-sided
    estimation problems (smoothing, Burg spectral estimation).

    **5. Matched to the Levinson recursion.** The $k_m$ are computed by Levinson
    as a by-product; no extra work is needed to obtain the lattice form.

### 3.11.3 From lattice back to direct form

To convert lattice to direct form, apply the **step-up** recursion (the inverse of
Levinson–Durbin):

$$
a_m^{(m)} = k_m,\qquad
a_i^{(m)} = a_i^{(m-1)} + k_m a_{m-i}^{(m-1)},\quad i = 1, \dots, m-1.
$$

To convert direct form to lattice, use the **step-down** recursion (Levinson itself).

### Worked example 3.11 — lattice for AR(2)

From Example 2.8 in Chapter 2: AR(2) with $r_0 = 6$, $r_1 = 4$, $r_2 = 1$ gave
reflection coefficients $k_1 = -2/3$, $k_2 = 1/2$ and AR coefficients
$a_1 = -1$, $a_2 = 0.5$.

**Stage 1 ($m = 1$).**

$$
f_1[n] = x[n] + k_1\, g_0[n-1] = x[n] - \tfrac{2}{3}x[n-1],
$$

$$
g_1[n] = x[n-1] + k_1\, x[n] = x[n-1] - \tfrac{2}{3}x[n].
$$

**Stage 2 ($m = 2$).**

$$
f_2[n] = f_1[n] + k_2\, g_1[n-1] = \bigl(x[n] - \tfrac{2}{3}x[n-1]\bigr) + \tfrac{1}{2}\bigl(x[n-2] - \tfrac{2}{3}x[n-1]\bigr)
$$

$$
= x[n] - \tfrac{2}{3}x[n-1] - \tfrac{1}{3}x[n-1] + \tfrac{1}{2}x[n-2] = x[n] - x[n-1] + \tfrac{1}{2}x[n-2].
$$

So $f_2[n] = x[n] + a_1 x[n-1] + a_2 x[n-2] = x[n] - x[n-1] + 0.5 x[n-2]$ — matching
the direct-form predictor from Example 2.7. ✓

**Interpretation.** The lattice adds one reflection coefficient per stage. Each
$|k_m| < 1$ guarantees a stable all-pole synthesis filter of order $m$ — the
**hierarchical stability** property.

---

## 3.12 Comparison and summary

### 3.12.1 Summary of the estimators

| Estimator | Criterion | Model | Form | Complexity |
|---|---|---|---|---|
| **MLE** | Max likelihood | Full pdf | $\arg\max p(Y;\theta)$ | $O(N)$ – $O(N^3)$ |
| **MMSE** | Mean-square | Joint pdf | $\mathbb E[S \mid Y]$ | Usually intractable |
| **LMMSE (vector)** | Mean-square, linear | Covariance only | $\Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$ | $O(m^3)$ matrix inverse |
| **Wiener (non-causal)** | LMMSE, stationary | PSDs | $S_{XY}/S_Y$ | Closed form |
| **Wiener (FIR)** | LMMSE, finite-length | Autocorr | Solve Toeplitz system | $O(p^2)$ Levinson |
| **Causal Wiener** | LMMSE, causal | PSDs + spectral factorization | $[S_{XY}/L]^{+}/L$ | $O(p \log p)$ |
| **Kalman** | LMMSE, state-space | State model | Predict-update recursion | $O(d^3)$ per step |
| **Linear predictor** | LMMSE, past only | Autocorr | Solve Yule–Walker | $O(p^2)$ Levinson |

### 3.12.2 When to use what

| Situation | Use |
|---|---|
| **Full pdf known, non-Gaussian** | MLE (or exact MMSE if tractable) |
| **Jointly Gaussian, stationary, offline** | Non-causal Wiener filter |
| **Jointly Gaussian, stationary, real-time** | Causal Wiener filter (via factorization) |
| **Jointly Gaussian, real-time, state-space** | Kalman filter |
| **Jointly Gaussian, FIR wanted** | FIR Wiener (Toeplitz solve) |
| **Signal modeling, compression, speech** | Linear prediction (Levinson) |
| **Whitening, coding** | Prediction-error filter (= spectral factor inverse) |

### 3.12.3 Chapter summary

1. **ML criterion** maximises the likelihood of the data; efficient under regularity; the **CRLB** is the variance floor for unbiased estimators.
2. **MMSE** (optimal, possibly nonlinear) equals the **conditional mean** $\mathbb E[S \mid Y]$; **LMMSE** restricts to linear estimators and equals the optimal **only for jointly Gaussian** problems.
3. **Orthogonality principle** — $\mathbb E[(S - \hat S)Y^\top] = 0$ — is the operational condition for LMMSE, giving the closed form $\hat S = \mu_S + \Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$.
4. **Wiener filter** is the LMMSE for stationary processes. **Non-causal:** $H = S_{XY}/S_Y$ — simple, closed form. **Causal:** requires spectral factorization; solved via one-sided Wiener–Hopf.
5. **FIR Wiener** reduces to a **Toeplitz** linear system, solved by **Levinson–Durbin** in $O(p^2)$.
6. **Kalman filter** is the LMMSE for state-space models; **predict-update** recursion; **Riccati** recursion for covariance; steady state = causal Wiener for stationary processes.
7. **Linear prediction** is the identity: optimal predictor = AR coefficients; prediction error = innovations; MMSE = innovations variance.
8. **Whitening filter** removes correlation, producing white innovations; its **inverse** reconstructs the signal. Basis of **LPC speech coding**.
9. **Levinson–Durbin** solves Toeplitz systems in $O(p^2)$; the **reflection coefficients** $k_m$ encode stability and the hierarchical structure.
10. **Lattice filter** implements prediction-error filtering using $k_m$; provides guaranteed stability, modularity, and numerical robustness.

---

## 3.13 Quiz

### A. Multiple choice (choose the single best answer)

1. The **maximum likelihood** estimator maximises:
   (a) the prior probability (b) the data likelihood (c) the conditional mean (d) the posterior probability

2. The **CRLB** gives:
   (a) a lower bound on the bias of an estimator (b) a lower bound on the variance of any unbiased estimator (c) the optimal MSE (d) the Fisher information

3. For i.i.d. Gaussian data $y[n] = A + w[n]$, the ML estimator of $A$ is:
   (a) median (b) sample mean (c) sample variance (d) RMS

4. The optimal **MMSE** estimator $\mathbb E[(S - \hat S)^2]$ is:
   (a) always linear (b) the conditional mean $\mathbb E[S \mid Y]$ (c) always Gaussian (d) the sample mean

5. The **LMMSE** estimator is:
   (a) $\mathbb E[S \mid Y]$ always (b) $\mu_S + \Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$ (c) $\mu_S$ (d) $\Sigma_{SY}\Sigma_{YY}^{-1}$

6. The **orthogonality principle** requires:
   (a) $\mathbb E[(S - \hat S)Y] = 0$ (b) $\mathbb E[(S - \hat S)Y^\top] = 0$ (c) $\mathbb E[SY] = 0$ (d) $\mathbb E[(S - \hat S)^2] = 0$

7. The **non-causal Wiener filter** is:
   (a) $S_{XY}/S_Y$ (b) $S_Y/S_{XY}$ (c) $S_{XY}S_Y$ (d) $\sigma_e^2 |H|^2$

8. The **causal Wiener filter** requires:
   (a) solving a Toeplitz system (b) spectral factorization and causal-part extraction (c) ML (d) Kalman recursion

9. The **Kalman gain** $K_n$ becomes small when:
   (a) observation noise $R$ is small (b) observation noise $R$ is large (c) prior covariance $P_{n|n-1}$ is small (d) both (b) and (c)

10. The **prediction error** $e[n] = x[n] - \hat x[n]$ for the optimal linear predictor is:
    (a) always zero (b) the innovations process (white) (c) Gaussian always (d) the autocorrelation of $x$

11. The **whitening filter** for an AR($p$) process is:
    (a) an FIR filter with coefficients $a_k$ (b) an IIR filter (c) an all-pass filter (d) a Kalman filter

12. The **lattice filter** stability condition is:
    (a) $|a_k| < 1$ for all $k$ (b) $|k_m| < 1$ for all $m$ (c) all poles inside the unit circle (d) it is always unstable

### B. Fill in the blank

1. The MLE of a scalar parameter $\theta$ is $\arg\max$\_\_\_\_\_ where $\ell$ is the log-likelihood.
2. The CRLB says $\operatorname{Var}[\hat\theta] \ge 1/$\_\_\_\_ where $I$ is the Fisher information.
3. The **efficient** estimator achieves the \_\_\_\_ at finite $N$.
4. The MMSE estimator is the conditional mean $\mathbb E[S \mid \_\_\_\_]$.
5. For jointly Gaussian signals, the LMMSE estimator equals the \_\_\_\_ estimator.
6. The **non-causal** Wiener filter is $H_{\text{nc}} = $\_\_\_\_\_.
7. The **causal** Wiener filter requires \_\_\_\_ factorization of $S_Y$.
8. The Kalman filter's steady-state behaviour matches the causal \_\_\_\_ filter for stationary signals.
9. The **whitening filter** produces a **white** output from a coloured process; the prediction-error filter is the \_\_\_\_ filter for an AR($p$) signal.
10. The lattice filter is parameterized by the \_\_\_\_\_\_ coefficients $k_1, \dots, k_p$.

### C. Define / explain

1. State the **ML criterion** and derive the MLE of the noise variance $\sigma^2$ for i.i.d. Gaussian data $y[n] \sim \mathcal N(0, \sigma^2)$. State whether the MLE is unbiased, and if not, propose the unbiased version.
2. State the **CRLB** and derive the Fisher information for the Gaussian-mean problem in Example 3.1. Show that the sample mean is efficient.
3. Distinguish **MMSE** and **LMMSE**. Under what conditions do they coincide? Give an example where they differ.
4. State the **orthogonality principle** and derive the LMMSE closed form for random vectors. Interpret each matrix in the formula.
5. Derive the **Wiener–Hopf equations** for the FIR Wiener filter from the orthogonality principle. Show that they reduce to a Toeplitz linear system.
6. Derive the **causal Wiener filter** formula via spectral factorization. Explain the role of the causal-part operation $[\cdot]_+$.
7. State the **Kalman filter recursion** (predict and update) and interpret each step. Why does the Kalman gain converge as $n \to \infty$?
8. Derive the **Yule–Walker** equations for linear prediction from the orthogonality principle. Show that they are identical to the AR parameter equations from Chapter 2.
9. Define the **whitening filter** and the **prediction-error filter**. Explain why they are identical for AR processes but different for ARMA.
10. Describe the **lattice filter** structure and state its key advantage over the direct-form filter. Why is $|k_m| < 1$ the stability condition?

---

## 3.14 Answer key

??? success "A. Multiple choice"
    1. **(b)** — the data likelihood.
    2. **(b)** — lower bound on the variance of any unbiased estimator.
    3. **(b)** — sample mean.
    4. **(b)** — conditional mean.
    5. **(b)** — $\mu_S + \Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$.
    6. **(b)** — $\mathbb E[(S - \hat S)Y^\top] = 0$.
    7. **(a)** — $S_{XY}/S_Y$.
    8. **(b)** — spectral factorization + causal part.
    9. **(d)** — both large $R$ and small $P_{n|n-1}$ reduce $K_n$.
    10. **(b)** — innovations (white).
    11. **(a)** — FIR with AR coefficients $a_k$.
    12. **(b)** — $|k_m| < 1$ for all $m$.

??? success "B. Fill in the blank"
    1. $\ell(\theta)$ (the log-likelihood)
    2. $I(\theta)$ (Fisher information)
    3. CRLB
    4. $Y$ (observations)
    5. MMSE
    6. $S_{XY}/S_Y$
    7. spectral
    8. Wiener
    9. whitening (or prediction-error)
    10. reflection

??? success "C. Define / explain — model answers"
    1. *ML criterion:* $\hat\theta_{\text{ML}} = \arg\max_\theta p_Y(Y; \theta)$; equivalently $\arg\max \ell(\theta) = \log p_Y(Y; \theta)$. *Derivation for variance:* $\ell(\sigma^2) = -\frac{N}{2}\ln(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum_n y[n]^2$. Setting $d\ell/d\sigma^2 = 0$: $-\frac{N}{2\sigma^2} + \frac{1}{2\sigma^4}\sum y[n]^2 = 0$ ⇒ $\hat\sigma^2_{\text{ML}} = \frac{1}{N}\sum y[n]^2$. *Bias:* $\mathbb E[\hat\sigma^2_{\text{ML}}] = \sigma^2$ if the data are zero-mean, so it is unbiased for zero-mean Gaussian (but biased by a factor $N/(N-1)$ when the mean is estimated as well). *Unbiased version (mean unknown):* $s^2 = \frac{1}{N-1}\sum(y[n] - \bar y)^2$ — the sample variance with $N-1$ denominator.
    2. *CRLB:* for any unbiased estimator $\hat\theta$, $\operatorname{Var}[\hat\theta] \ge 1/I(\theta)$ with $I(\theta) = -\mathbb E[\partial^2\ell/\partial\theta^2]$. *Fisher information for the Gaussian-mean problem:* $\ell(A) = -\frac{N}{2}\ln(2\pi\sigma^2) - \frac{1}{2\sigma^2}\sum(y[n]-A)^2$; $\partial^2\ell/\partial A^2 = -N/\sigma^2$, so $I(A) = N/\sigma^2$. *Sample mean:* $\operatorname{Var}[\bar y] = \sigma^2/N = 1/I(A)$ — equality holds, so the sample mean is **efficient** (attains the CRLB) for every $N$.
    3. *MMSE:* minimises $\mathbb E[(S-\hat S)^2]$ over **all** estimators $\hat S(Y)$; solution is the conditional mean $\mathbb E[S \mid Y]$. *LMMSE:* minimises the same criterion over **linear** estimators $\hat S = a + b^\top Y$; solution is $\mu_S + \Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$. *Equivalence:* for **jointly Gaussian** $(S, Y)$, the conditional mean is affine in $Y$, so the two coincide. *Difference example:* $S$ uniform on $\{-1, +1\}$, $Y = S + V$ with $V$ uniform on $[-0.5, 0.5]$. The MMSE is $\mathrm{sign}(Y)$ — a nonlinear estimator. The LMMSE is a linear function, which has a larger MSE. The difference is small for small noise, but exact-MMSE beats LMMSE when the prior is non-Gaussian.
    4. *Orthogonality:* $\mathbb E[(S - \hat S)Y^\top] = 0$ and $\mathbb E[S - \hat S] = 0$. *Derivation:* write $\hat S = a + b^\top Y$; the orthogonality of the error to $1$ gives $a = \mu_S - b^\top\mu_Y$; the orthogonality to $Y$ gives $\Sigma_{SY} - b^\top\Sigma_{YY} = 0$, so $b = \Sigma_{YY}^{-1}\Sigma_{SY}^\top = \Sigma_{YY}^{-1}\Sigma_{YS}$. *Interpretation:* $\Sigma_{SY}$ = cross-covariance (how correlated the target is with the observations); $\Sigma_{YY}^{-1}$ = inverse observation covariance (how much each observation is worth); the product weights each observation by its information and correlation to the target.
    5. *Orthogonality:* $\mathbb E[(x[n] - \hat x[n])y[n-m]] = 0$ for $m = 0, \dots, M-1$. *Expansion:* substituting $\hat x[n] = \sum_{k=0}^{M-1}h[k]y[n-k]$ and taking expectations gives $\sum_k h[k] R_Y[m-k] = R_{XY}[m]$. *Toeplitz form:* the left side is a convolution of $h$ with $R_Y$, evaluated at $m$. In matrix form with $\mathbf h = (h[0],\dots,h[M-1])^\top$ and Toeplitz $\mathbf R_Y$ with $[\mathbf R_Y]_{ij} = R_Y[i-j]$, the system is $\mathbf R_Y\mathbf h = \mathbf r_{XY}$ where $\mathbf r_{XY} = (R_{XY}[0],\dots,R_{XY}[M-1])^\top$. Solve by Levinson–Durbin in $O(M^2)$.
    6. *Formula:* $H_{\text{cas}}(z) = [S_{XY}(z)/(\sigma_e^2 L(z^{-1}))]_+ / L(z)$ where $S_Y = \sigma_e^2 L L^{-1}$ is the spectral factorization. *Role of $[\cdot]_+$:* the one-sided Wiener–Hopf has orthogonality only for $m \ge 0$. In the $z$-domain, this means we can only use the causal part of the Laurent expansion — terms $z^0, z^{-1}, z^{-2}, \dots$. The operation $[\cdot]_+$ takes exactly those terms, discarding the anti-causal part. Without it, we would recover the non-causal filter $S_{XY}/S_Y$, which uses future observations and is not realisable.
    7. *Recursion:* predict step $\hat x_{n|n-1} = A\hat x_{n-1|n-1}$, $P_{n|n-1} = AP_{n-1|n-1}A^\top + BQB^\top$; Kalman gain $K_n = P_{n|n-1}C^\top(CP_{n|n-1}C^\top + R)^{-1}$; update step $\hat x_{n|n} = \hat x_{n|n-1} + K_n(y[n] - C\hat x_{n|n-1})$, $P_{n|n} = (I - K_nC)P_{n|n-1}$. *Convergence:* the Riccati recursion $P_{n|n} = f(P_{n-1|n-1})$ is a contraction under stabilisability and detectability, so $P_{n|n} \to P^\star$ and $K_n \to K^\star$, independent of $P_{0|0}$. The steady-state filter is LTI, matching the causal Wiener filter.
    8. *Orthogonality:* $\mathbb E[e[n]x[n-k]] = 0$ for $k = 1, \dots, p$, where $e[n] = x[n] + \sum_{m=1}^p a_m x[n-m]$. *Expansion:* $\mathbb E[x[n]x[n-k]] + \sum_m a_m \mathbb E[x[n-m]x[n-k]] = 0$ gives $R_X[k] + \sum_m a_m R_X[k-m] = 0$ for $k = 1, \dots, p$. *Yule–Walker:* exactly the AR model equations from §2.10. This shows the optimal linear predictor is the **AR coefficient vector** — the same vector that defines the AR model — and the MMSE prediction error equals the innovations variance $\sigma_w^2$.
    9. *Whitening filter:* an LTI filter $H_w = 1/H$ such that $H_w * X$ is white. *Prediction-error filter:* the filter whose output is the prediction error of the optimal linear predictor. *Identity for AR:* for AR($p$), the AR polynomial $1 + \sum a_k z^{-k}$ is precisely the prediction-error filter; its output is the input to the AR model, i.e., the innovations — so $H_w = 1 + \sum a_k z^{-k}$ is FIR and $H_w$ = prediction-error filter. *Difference for ARMA:* for an ARMA($p,q$), the optimal predictor is still IIR from the AR part, but the **whitening filter** must also invert the MA part $1 + \sum b_k z^{-k}$, which is **IIR** (a non-trivial all-pole filter). So the whitening filter is IIR, not the same as the (FIR) prediction-error filter that only handles the AR part.
    10. *Structure:* a cascade of $p$ stages, each parameterised by a single reflection coefficient $k_m$. Stage $m$ takes inputs $(f_{m-1}, g_{m-1})$ and outputs $(f_m, g_m)$ via $f_m = f_{m-1} + k_m g_{m-1}z^{-1}$, $g_m = g_{m-1}z^{-1} + k_m f_{m-1}$. *Key advantages:* (i) **guaranteed stability** iff $|k_m| < 1$ for all $m$ — a simple, local condition, unlike the global root condition for direct-form $a_k$; (ii) **modular** — stages can be added/removed without recomputing others; (iii) **robust under quantization** — coefficients bounded in magnitude; (iv) **order-adaptive** — easily grow/shrink $p$ at runtime; (v) **provides both forward and backward errors** — useful for smoothing and Burg methods. *Stability:* each stage is an all-pass in the reverse direction; if $|k_m| < 1$, the reflected energy is less than the incident energy, so the cascade cannot have poles outside the unit circle. This gives the **hierarchical stability** property: every lower-order filter in the cascade is also stable.

---

## Further reading

- Hayes, *Statistical Digital Signal Processing and Modeling*, Wiley, 1996 — **Chapters 4–7 and 9**. The closest match to this chapter.
- Kay, *Fundamentals of Statistical Signal Processing: Estimation Theory*, Prentice-Hall, 1993 — the definitive reference for ML, CRLB, and MMSE. Chapters 3–7.
- Kay, *Modern Spectral Estimation*, Prentice-Hall, 1988 — for the linear-prediction and lattice chapters.
- Kailath, Sayed & Hassibi, *Linear Estimation*, Prentice-Hall, 2000 — the unified framework for LMMSE, Wiener, and Kalman. Rigorous, comprehensive.
- Anderson & Moore, *Optimal Filtering*, Dover, 2005 — the classical Kalman filter reference.
- Simon, *Optimal State Estimation*, Wiley, 2006 — modern unified treatment, with practical implementation advice.
- Haykin, *Adaptive Filter Theory*, Pearson, 5th ed., 2014 — Chapter 1 for prediction and lattice; Chapters 2–3 for Wiener; Chapters 4–9 for adaptive (Chapter 4 of this course).
- Makhoul, "Linear Prediction: A Tutorial Review", *Proc. IEEE*, 63(4), 1975 — the classic pedagogical reference for LPC and lattice.
- Markel & Gray, *Linear Prediction of Speech*, Springer, 1976 — the definitive source for the source-filter model and LPC.
- SciPy `signal.lfilter`, `statsmodels.tsa`, `filterpy.kalman`; MATLAB `dsp.KalmanFilter`, `levinson`, `aryule`, `lpc`, `dsp.LatticeFilter`.

---

## Suggested tutorial and lab problems for this chapter

1. **CRLB exercise.** $y[n] = \theta s[n] + w[n]$ for known $s[n]$, $w$ white Gaussian. Derive the MLE of $\theta$ and its CRLB. Show the MLE is $\hat\theta = \sum y[n]s[n]/\sum s^2[n]$, and its variance equals the CRLB.
2. **LMMSE with correlated observations.** $Y_1 = S + V_1$, $Y_2 = S + V_2$, with $S \sim \mathcal N(0, 1)$, $V_i \sim \mathcal N(0, 0.5)$, $V_1 \perp V_2$. Compute the LMMSE of $S$ from $(Y_1, Y_2)$ and compare with the LMMSE from $Y_1$ alone.
3. **Non-causal Wiener filter.** For $Y[n] = X[n] + V[n]$ with $X \sim$ AR(1) ($a = 0.8$, $\sigma_W = 1$), $V$ white ($\sigma_V = 0.3$), compute $H_{\text{nc}}(e^{j\omega})$. Numerically integrate to find the MSE.
4. **Causal Wiener via factorization.** Same setup as Problem 3, compute the causal Wiener filter using spectral factorization (following Example 3.6). Compare MSE with the non-causal case.
5. **Kalman filter scalar.** Model $x[n] = 0.95 x[n-1] + u[n]$, $u \sim \mathcal N(0, 0.05)$; $y[n] = 2x[n] + v[n]$, $v \sim \mathcal N(0, 1)$. Simulate 100 steps, run the Kalman filter, and plot the true state, observation, and filtered estimate. Verify the Riccati recursion converges.
6. **Kalman vs. Wiener.** For the stationary model of Problem 5, run the Kalman filter for 1000 steps and compare the steady-state gain with the causal Wiener filter designed from the AR spectrum. Verify they match to $10^{-4}$.
7. **Linear predictor.** For $R_X[0] = 10$, $R_X[1] = 8$, $R_X[2] = 6$, $R_X[k] = 0$ for $|k| > 2$, find the optimal order-2 predictor and the minimum prediction error. Use Levinson–Durbin.
8. **Whitening.** Design the whitening filter for the AR(2) process of Example 2.8. Verify that passing the AR(2) process through the whitening filter yields a white residual (autocorrelation within $\pm 3/\sqrt N$ for $k \ne 0$).
9. **Lattice conversion.** Given the AR(3) coefficients $a = (-1.5, 1.1, -0.3)$, compute the reflection coefficients using step-down (inverse Levinson). Implement the lattice filter, verify that it produces the same prediction error as the direct-form filter.
10. **Speech LPC.** Record a segment of voiced speech (~2 s at 16 kHz). Implement a pitch-synchronous LPC analyser with $p = 12$. Plot the LPC spectrum overlaid on the Welch spectrum. Verify that the peaks of the LPC spectrum correspond to the formants (F1, F2, F3).
11. **Riccati convergence.** For the scalar model $A = 1.0$, $Q = 0.1$, $R = 1$, plot $P_{n|n}$ vs. $n$ starting from $P_{0|0} = 100$. Verify convergence to the fixed point in under 20 iterations. Repeat with $Q = 0.001$; comment on the convergence rate.
12. **Real-time consideration.** In the causal Wiener filter of Problem 4, compare its performance on a **finite** sample (N = 200) to the non-causal filter. Plot the two filters' outputs on the same axes and compute the MSE. Explain the intuition behind the ~35% error penalty.