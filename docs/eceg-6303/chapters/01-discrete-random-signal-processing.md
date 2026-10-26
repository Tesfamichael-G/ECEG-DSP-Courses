# Chapter 1 — Discrete Random Signal Processing

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** discrete random processes · expectation, variance,
    covariance · scalar product and energy of discrete signals · Parseval's
    theorem · Wiener–Khinchine relation · power spectral density · periodogram ·
    sample autocorrelation · Sum Decomposition Theorem · Spectral Factorization
    Theorem · discrete random signal processing by linear systems · simulation
    of white noise · lowpass filtering of white noise.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Define a **discrete-time random process** and compute its **mean, variance,
   autocovariance, and autocorrelation** from the ensemble definition.
2. State and apply **Parseval's theorem** for finite-length and for
   infinite-energy discrete signals.
3. Derive the **Wiener–Khinchine relation** and compute the **power spectral
   density (PSD)** of a WSS sequence from its autocorrelation.
4. Compute and interpret the **periodogram**, and state its bias and variance
   properties.
5. Compute the **sample autocorrelation** in both biased and unbiased forms.
6. Apply the **Sum Decomposition Theorem** (discrete Karhunen–Loève) to
   decompose a finite-length random vector into uncorrelated components.
7. State the **Spectral Factorization Theorem** and factor a given rational
   PSD into causal minimum-phase factors.
8. Compute the **output PSD** and output autocorrelation of an LTI system
   driven by a WSS input.
9. **Simulate white noise** (uniform, Gaussian, arbitrary distribution) and
   **colour it** by lowpass filtering.

---

## 1.1 Discrete random processes

A **discrete-time random process** (also *random sequence* or *stochastic
sequence*) is a family of random variables $\{X[n] : n \in \mathbb{Z}\}$
indexed by integer time. For each fixed outcome $\omega$ in the underlying
probability space, the sequence $n \mapsto X[n, \omega]$ is a **realization**
or **sample function** of the process; the collection of all realizations is
the **ensemble**.

!!! note "Conventions"
    We use **$n$** for the discrete-time index, **$[\cdot]$** for discrete-time
    sequences, and $\mathbb{Z}$ for the set of all integers. Lowercase $x[n]$
    denotes a *specific realization*; uppercase $X[n]$ the *random variable* at
    time $n$.

For **real-valued** processes we usually work with the first and second moments:

$$
\mu_X[n] = \mathbb{E}[X[n]],\qquad
\sigma_X^2[n] = \operatorname{Var}[X[n]] = \mathbb{E}\bigl[(X[n] - \mu_X[n])^2\bigr],
$$

$$
R_X[m, n] = \mathbb{E}\bigl[X[m] X[n]\bigr],\qquad
C_X[m, n] = R_X[m,n] - \mu_X[m]\mu_X[n].
$$

For **complex-valued** processes, $X[n]$ is complex and $R_X[m,n] = \mathbb{E}[X[m]X^*[n]]$, where $*$ denotes the complex conjugate.

### 1.1.1 Wide-sense stationarity (WSS)

A discrete random process is **wide-sense stationary** if:

1. $\mu_X[n] = \mu_X$ is constant (independent of $n$), and
2. $R_X[m, n] = R_X[m - n]$ depends only on the **lag** $k = m - n$.

For a WSS process, we write $R_X[k] = \mathbb{E}[X[n+k] X^*[n]]$ (the convention is
$X[n+k]X^*[n]$, so the autocorrelation is Hermitian-symmetric: $R_X[-k] = R_X^*[k]$;
for real processes, $R_X[-k] = R_X[k]$).

**Properties of $R_X[k]$ (real WSS):**

| Property | Interpretation |
|---|---|
| $R_X[0] = \mathbb{E}[X^2[n]] = \operatorname{Var}[X] + \mu_X^2$ | Power (or energy per sample, if zero-mean) |
| $R_X[-k] = R_X[k]$ | Symmetry |
| $\lvert R_X[k]\rvert \le R_X[0]$ | Peak at lag zero |
| $R_X[k]$ positive-definite | For any coefficients $a_i$: $\sum_{i,j} a_i R_X[i-j] a_j^* \ge 0$ |

**WSS in the frequency domain** — the PSD $S_X(e^{j\omega})$ is the
DTFT of $R_X[k]$ (§1.7), and it inherits non-negativity from the
positive-definiteness of $R_X$.

---

## 1.2 Scalar product and energy of discrete signals

For two deterministic finite-length sequences $x[n], y[n]$, $n = 0, \dots, N-1$, the **inner product** is

$$
\langle x, y\rangle = \sum_{n=0}^{N-1} x[n] y^*[n].
$$

This inner product makes $\mathbb{C}^N$ a **Hilbert space**. The **norm** (or **energy**) of $x$ is

$$
\|x\|^2 = \langle x, x\rangle = \sum_{n=0}^{N-1} |x[n]|^2.
$$

**For infinite-length sequences**, the energy is

$$
E_x = \sum_{n=-\infty}^{\infty} |x[n]|^2,
$$

finite for *energy signals*; for *power signals* the relevant quantity is the average power,

$$
P_x = \lim_{N\to\infty} \frac{1}{2N+1} \sum_{n=-N}^{N} |x[n]|^2.
$$

For a **WSS random process**, the average power is $R_X[0] = \mathbb{E}[|X[n]|^2]$ — the fundamental bridge between the deterministic and the statistical view.

---

## 1.3 Parseval's theorem

For a finite-length signal $x[n]$ with $N$-point DFT $X[k] = \sum_{n=0}^{N-1} x[n] e^{-j2\pi kn/N}$:

$$
\boxed{\;\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |X[k]|^2\;}
$$

The factor $1/N$ appears because the DFT convention absorbs $N$ into the inverse transform. For the **DTFT** of an infinite-energy sequence, the form is

$$
\boxed{\;\sum_{n=-\infty}^{\infty} |x[n]|^2 = \frac{1}{2\pi} \int_{-\pi}^{\pi} |X(e^{j\omega})|^2\, d\omega\;}
$$

and for a **WSS random process**, replacing $|X|^2$ by the PSD:

$$
\boxed{\;R_X[0] = \frac{1}{2\pi}\int_{-\pi}^{\pi} S_X(e^{j\omega})\, d\omega\;}
$$

The last form is the statistical Parseval — the total power equals the integral of the PSD.

### Worked example 1.1 — Parseval for a rectangular pulse

Let $x[n] = 1$ for $n = 0, 1, \dots, N-1$ and zero otherwise.

**Time domain energy.** $\sum_{n=0}^{N-1} |x[n]|^2 = N$.

**Frequency domain.** The $N$-point DFT of $x$ is

$$
X[k] = \sum_{n=0}^{N-1} e^{-j2\pi kn/N} = \begin{cases} N, & k = 0, \\ 0, & k \ne 0. \end{cases}
$$

**Parseval's check.** $\frac{1}{N}\sum_{k=0}^{N-1} |X[k]|^2 = \frac{1}{N}\cdot N^2 = N$. ✓

**DTFT check.** $X(e^{j\omega}) = \frac{\sin(N\omega/2)}{\sin(\omega/2)} e^{-j(N-1)\omega/2}$, so $|X|^2 = \sin^2(N\omega/2)/\sin^2(\omega/2)$. Then $\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{\sin^2(N\omega/2)}{\sin^2(\omega/2)}\,d\omega = N$ ✓.

### Worked example 1.2 — Statistical Parseval for an AR(1)

Let $X[n] = a X[n-1] + W[n]$ with $|a| < 1$ and $W$ white of variance $\sigma_W^2$. From Chapter 5 of ECEG-6201, we have

$$
R_X[k] = \frac{\sigma_W^2}{1-a^2}\, a^{|k|},\qquad
S_X(e^{j\omega}) = \frac{\sigma_W^2}{|1 - a e^{-j\omega}|^2}.
$$

**Statistical Parseval.** The right side is

$$
\frac{1}{2\pi}\int_{-\pi}^{\pi} \frac{\sigma_W^2}{|1-ae^{-j\omega}|^2}\, d\omega = \frac{\sigma_W^2}{2\pi}\int_{-\pi}^{\pi} \frac{1}{1-2a\cos\omega + a^2}\, d\omega.
$$

The integral evaluates to $\frac{2\pi}{1-a^2}$ (standard contour-integral result), so the RHS is $\frac{\sigma_W^2}{1-a^2}$. ✓

This matches $R_X[0] = \frac{\sigma_W^2}{1-a^2}$. ✓

---

## 1.4 Autocorrelation and autocovariance

### 1.4.1 Definitions

For a WSS process $X[n]$:

$$
R_X[k] = \mathbb{E}\bigl[X[n+k]X^*[n]\bigr],\qquad
C_X[k] = R_X[k] - |\mu_X|^2.
$$

The autocorrelation measures **linear dependence** between samples $k$ apart; the autocovariance subtracts the mean so that it vanishes when $k \to \infty$ (for mixing processes).

### 1.4.2 Autocorrelation of common processes

| Process | $R_X[k]$ |
|---|---|
| White noise (variance $\sigma^2$) | $\sigma^2 \delta[k]$ |
| AR(1) ($a$, $\sigma_W^2$) | $\frac{\sigma_W^2}{1-a^2} a^{|k|}$ |
| Complex exponential $e^{j\omega_0 n}$ (deterministic) | $e^{j\omega_0 k}$ |
| Random-phase sinusoid | $\frac{A^2}{2}\cos(\omega_0 k)$ |
| Moving average MA(1): $X[n] = W[n] + b W[n-1]$ | $\sigma_W^2(1 + b^2)\delta[k] + b\sigma_W^2(\delta[k-1] + \delta[k+1])$ |

### 1.4.3 Sample autocorrelation

Given a finite realization $x[0], \dots, x[N-1]$, two estimators of $R_X[k]$:

**Biased estimator:**

$$
\boxed{\;\hat{R}[k] = \frac{1}{N}\sum_{n=0}^{N-1-|k|} x[n+|k|] x^*[n]\;,\quad |k| \le N-1.\;}
$$

**Unbiased estimator:**

$$
\hat{R}_u[k] = \frac{1}{N - |k|}\sum_{n=0}^{N-1-|k|} x[n+|k|]x^*[n].
$$

**Trade-offs.**

- **Bias.** $\hat{R}_u$ is unbiased; $\hat{R}$ is biased by a factor of $(N-|k|)/N$, which is $O(1/N)$ for $|k| \ll N$.
- **Variance.** $\hat{R}_u$ has *larger* variance (its denominator $N-|k|$ becomes small as $|k| \to N$), so $\hat{R}$ has *lower* mean-square error at large lags.
- **Positive-definiteness.** $\hat{R}$ is guaranteed positive-definite (needed for spectral estimation); $\hat{R}_u$ is not.

**Why the biased version is standard.** For spectral estimation, we need the estimated autocorrelation to be a valid PSD generator — i.e., its DTFT must be non-negative. Only the biased $\hat{R}$ guarantees this. So we use $\hat{R}$, accepting its $O(1/N)$ bias in return for a well-behaved spectrum.

---

## 1.5 Wiener–Khinchine relation and power spectral density

For a WSS process $X[n]$ with $\sum_{k} |R_X[k]| < \infty$:

$$
\boxed{\;S_X(e^{j\omega}) = \sum_{k=-\infty}^{\infty} R_X[k]\, e^{-j\omega k}\;}
$$

with inverse

$$
R_X[k] = \frac{1}{2\pi}\int_{-\pi}^{\pi} S_X(e^{j\omega})\, e^{j\omega k}\, d\omega.
$$

**Properties of the PSD.**

| Property | Reason |
|---|---|
| $S_X(e^{j\omega}) \ge 0$ | $R_X$ is positive-definite |
| $S_X(e^{j\omega})$ is $2\pi$-periodic | Definition as a DTFT |
| $S_X(e^{j\omega}) = S_X(e^{-j\omega})$ (real $X$) | $R_X[k] = R_X[-k]$ |
| $\frac{1}{2\pi}\int_{-\pi}^{\pi} S_X(e^{j\omega})\, d\omega = R_X[0]$ | Statistical Parseval |

When the condition $\sum |R_X[k]| < \infty$ is violated (e.g., a pure sinusoid), $S_X$ contains impulses (Dirac deltas) at the corresponding frequencies, and the integral is interpreted distributionally.

### Worked example 1.3 — PSD of the moving average MA(1)

Let $X[n] = W[n] + b W[n-1]$ with $W$ white of variance $\sigma_W^2$.

**Autocorrelation.** $\mathbb{E}[X[n]X[n+k]] = \sigma_W^2(\delta[k] + b\delta[k-1] + b\delta[k+1] + b^2 \delta[k])$, so

$$
R_X[0] = \sigma_W^2(1 + b^2),\quad R_X[\pm 1] = b\sigma_W^2,\quad R_X[k] = 0\ (|k| \ge 2).
$$

**PSD.**

$$
S_X(e^{j\omega}) = \sigma_W^2[(1+b^2) + b e^{-j\omega} + b e^{j\omega}] = \sigma_W^2[(1+b^2) + 2b\cos\omega].
$$

**Check.** $\frac{1}{2\pi}\int_{-\pi}^\pi \sigma_W^2[(1+b^2) + 2b\cos\omega]\, d\omega = \sigma_W^2(1+b^2) = R_X[0]$ ✓.

**Non-negativity constraint.** For $S_X$ to be a valid PSD we need $1 + b^2 + 2b\cos\omega \ge 0$ for all $\omega$; the minimum is $1 + b^2 - 2|b| = (1-|b|)^2 \ge 0$, always satisfied. So MA(1) is always a valid PSD generator.

---

## 1.6 The periodogram

For a finite-length realization $x[0], \dots, x[N-1]$, the **periodogram** is

$$
\boxed{\;\hat{S}_X(e^{j\omega}) = \frac{1}{N}\Bigl|\sum_{n=0}^{N-1} x[n]\, e^{-j\omega n}\Bigr|^2\;}
$$

Evaluated at the DFT frequencies $\omega_k = 2\pi k/N$, this equals $\frac{1}{N}|X[k]|^2$, where $X[k]$ is the $N$-point DFT.

**Two properties (developed in Chapter 2):**

- **Asymptotically unbiased:** $\mathbb{E}[\hat{S}_X(e^{j\omega})] \to S_X(e^{j\omega})$ as $N \to \infty$, but the bias is nonzero for finite $N$ (the periodogram is a smoothed version of $S_X$ with the Fejér kernel).
- **Inconsistent:** $\operatorname{Var}[\hat{S}_X(e^{j\omega})] \approx S_X^2(e^{j\omega})$ regardless of $N$. The estimate fluctuates as much as its mean; averaging over independent segments is required to reduce the variance (Bartlett, Welch; §2.1).

### Worked example 1.4 — periodogram of white noise

Let $x[n]$ be white Gaussian of variance $\sigma^2$, $N$ samples.

**Bias.** $\mathbb{E}[\hat{S}_X(e^{j\omega_k})] = \sigma^2$ (unbiased at DFT frequencies when $x$ is white) — for general $S_X$ the expectation is a smoothed version.

**Variance.** $\operatorname{Var}[\hat{S}_X(e^{j\omega_k})] = \sigma^4$, so the standard deviation of $\hat{S}_X$ equals its mean. This is the fundamental inconsistency: averaging over the same $N$-sample realization does **not** reduce variance; only averaging over **independent realisations** does.

**Distribution.** $\hat{S}_X(e^{j\omega_k})/\sigma^2 \sim \chi^2_2/2 = \mathrm{Exp}(1)$, independent (asymptotically) across bins. A plot of the periodogram of white noise therefore looks like exponentially-distributed noise, not a flat line — even though its expectation is flat.

**Reducing variance.** Welch's method averages periodograms over $K$ segments; the variance drops to $\sigma^4/K$. See Chapter 2 for details.

---

## 1.7 Sum Decomposition Theorem

For a length-$N$ real-valued zero-mean random vector $\mathbf{X} = (X[0], \dots, X[N-1])^\top$ with covariance matrix $\mathbf{R} = \mathbb{E}[\mathbf{X}\mathbf{X}^\top] \in \mathbb{R}^{N\times N}$, the **Sum Decomposition Theorem** (also called the **discrete Karhunen–Loève** or **principal-component** decomposition) is the eigen-decomposition

$$
\boxed{\;\mathbf{R} = \sum_{k=1}^{N} \lambda_k\, \mathbf{v}_k\, \mathbf{v}_k^\top\;}
$$

with $\mathbf{v}_k$ orthonormal eigenvectors and $\lambda_1 \ge \lambda_2 \ge \dots \ge \lambda_N \ge 0$ eigenvalues. The corresponding expansion

$$
\boxed{\;\mathbf{X} = \sum_{k=1}^{N} Z_k\, \mathbf{v}_k,\qquad Z_k = \mathbf{v}_k^\top \mathbf{X}\;}
$$

has the property that the coefficients $Z_k$ are **uncorrelated**: $\mathbb{E}[Z_i Z_j] = \lambda_i \delta_{ij}$.

### Worked example 1.5 — Sum decomposition for an AR(1)

Let $\mathbf{X} = (X[0], \dots, X[N-1])^\top$ be a length-$N$ block of an AR(1)
process with $R_X[k] = \rho^{|k|}$ (unit variance). The covariance matrix is
**Toeplitz**:

$$
\mathbf{R} = \begin{bmatrix} 1 & \rho & \rho^2 & \cdots \\ \rho & 1 & \rho & \cdots \\ \rho^2 & \rho & 1 & \cdots \\ \vdots & & & \ddots \end{bmatrix}.
$$

For $N = 2$: $\mathbf{R} = \begin{bmatrix} 1 & \rho \\ \rho & 1 \end{bmatrix}$. Eigenvalues $\lambda_1 = 1 + \rho$, $\lambda_2 = 1 - \rho$; eigenvectors $\mathbf{v}_1 = (1, 1)/\sqrt 2$, $\mathbf{v}_2 = (1, -1)/\sqrt 2$.

For $\rho = 0.8$: $\lambda_1 = 1.8$, $\lambda_2 = 0.2$. Almost all the energy is in the "symmetric" direction; the "antisymmetric" direction is nearly noise. This is why a **first-order linear predictor** captures most of the signal in a highly-correlated process.

For general $N$ and $\rho$, the eigenvectors are **discrete cosine transform (DCT)** basis vectors only in the limit $\rho \to 1$; for finite $\rho$ they are slightly different. The **eigenvalues decay** as a geometric-like sequence, and the first few principal components capture most of the variance.

**Compression implication.** Retaining the first $M$ components of the expansion — $\hat{\mathbf{X}}_M = \sum_{k=1}^{M} Z_k \mathbf{v}_k$ — minimises the mean-square reconstruction error over all rank-$M$ linear projections. The residual energy is $\sum_{k=M+1}^{N}\lambda_k$.

---

## 1.8 Spectral Factorization Theorem

The **Spectral Factorization Theorem** (also called the **Fejér–Riesz theorem** for trigonometric polynomials, and the **Kolmogorov–Wiener** or **Wiener–Hopf** factorization in the stochastic-process setting) is:

!!! abstract "Spectral Factorization"
    Let $S_X(z)$ be a **rational** power spectral density that is **positive on
    the unit circle** ($S_X(e^{j\omega}) > 0$ for all $\omega$) and satisfies the
    **Paley–Wiener condition** $\int_{-\pi}^{\pi}\log S_X(e^{j\omega})\,d\omega > -\infty$.
    Then there exists a **unique** factorization

    $$
    \boxed{\;S_X(z) = \sigma_e^2\, H(z)\, H(z^{-1})\;}
    $$

    where $H(z)$ is **causal** ($H(z) = \sum_{k\ge 0} h_k z^{-k}$ with $h_0 = 1$),
    **minimum-phase** (all poles and zeros strictly inside $|z| < 1$), and
    $\sigma_e^2 = \exp\bigl(\frac{1}{2\pi}\int\log S_X(e^{j\omega})d\omega\bigr) > 0$
    is the **innovations variance**.

**Why rational is not a serious restriction.** Any PSD that arises from a finite-dimensional linear system driven by white noise has rational $z$-transform; this covers every AR, MA, ARMA, and state-space model. Non-rational PSDs are handled by truncation or by the general Paley–Wiener theorem.

**Why minimum-phase matters.** Among all causal filters with a given magnitude response, the minimum-phase factor $H(z)$ has the **smallest group delay** and its inverse $H^{-1}(z)$ is also causal — so the **whitening filter** $1/H(z)$ is realisable in real time. Any other factorization would give a causal $H$ with a **larger** $\sigma_e^2$ (non-minimum-phase zeros waste innovation energy).

### How to factorize a rational PSD

Given $S_X(z)$ as a ratio of polynomials in $z$:

1. **Find the roots of the numerator and denominator polynomials** in $z$.
2. **Split each root pair** $\{z_0, 1/z_0^*\}$: put $z_0$ in $H(z^{-1})$ (its zero is at $z_0$, outside the unit circle) and $1/z_0^*$ in $H(z)$ (its zero is at $1/z_0^*$, inside the unit circle) — or vice versa, depending on the convention. The **minimum-phase** choice puts all interior roots in $H(z)$.
3. **Normalize** so that $H(0) = 1$ (i.e., $h_0 = 1$).
4. **Compute** $\sigma_e^2$ from the leading coefficients.

### Worked example 1.6 — factorize an ARMA(1,1) PSD

$$
S_X(z) = \frac{(1 - 0.5 z^{-1})(1 - 0.5 z)}{(1 - 0.8 z^{-1})(1 - 0.8 z)}.
$$

**Roots.** Numerator: $z = 0.5$ (from $1 - 0.5z^{-1}$, zero at $z = 0.5$) and $z = 2$ (from $1 - 0.5z$, zero at $z = 2$). Denominator: $z = 0.8$ and $z = 1.25$.

**Minimum-phase factor.** Put the interior roots (0.5, 0.8) in $H(z)$:

$$
H(z) = \frac{1 - 0.5 z^{-1}}{1 - 0.8 z^{-1}},\qquad \sigma_e^2 = 1.
$$

**Verify.**

$$
\sigma_e^2 H(z) H(z^{-1}) = \frac{(1 - 0.5z^{-1})(1 - 0.5z)}{(1 - 0.8z^{-1})(1 - 0.8z)} = S_X(z). ✓
$$

**Whitening filter.** $H^{-1}(z) = \frac{1 - 0.8z^{-1}}{1 - 0.5z^{-1}}$ — causal, stable (both poles and zeros inside the unit circle), so realisable.

**Non-minimum-phase alternative.** Using the outer roots:

$$
H_{np}(z) = \frac{1 - 2 z^{-1}}{1 - 1.25 z^{-1}}\cdot\text{(constant)}.
$$

This has the same $|H(e^{j\omega})|^2$ but is non-minimum-phase — its inverse is not causal. Avoid.

### Worked example 1.7 — factorize an AR(2) PSD

$$
S_X(z) = \frac{\sigma_w^2}{(1 - 0.6 z^{-1})(1 - 0.7 z^{-1})(1 - 0.6 z)(1 - 0.7 z)}.
$$

Minimum-phase factor:

$$
H(z) = \frac{1}{(1 - 0.6 z^{-1})(1 - 0.7 z^{-1})},\qquad \sigma_e^2 = \sigma_w^2.
$$

The **innovations** are exactly the driving white noise $w[n]$, as expected — this is an AR(2) process.

The **predictor** is $a_1 = 0.6 + 0.7 = 1.3$, $a_2 = -0.6\cdot 0.7 = -0.42$, so $\hat x[n] = 1.3 x[n-1] - 0.42 x[n-2]$, and the innovations (prediction errors) are white.

---

## 1.9 Discrete random signal processing by linear systems

Let $X[n]$ be WSS input to an LTI system with impulse response $h[n]$ and frequency response $H(e^{j\omega})$. The output is

$$
Y[n] = \sum_{k=-\infty}^{\infty} h[k] X[n-k].
$$

**Output statistics.**

$$
\boxed{\;\mu_Y = \mu_X \cdot H(e^{j0}) = \mu_X \sum_k h[k]\;}
$$

$$
\boxed{\;S_Y(e^{j\omega}) = |H(e^{j\omega})|^2\, S_X(e^{j\omega})\;}
$$

$$
\boxed{\;R_Y[k] = R_X[k] * \bigl(h * h^\dagger\bigr)[k]\;}
$$

where $h^\dagger[k] = h^*[-k]$ is the **matched (time-reversed conjugate)** filter, and $*$ denotes convolution.

**Interpretation of the three formulas.**

- The **mean** is filtered by $H(e^{j0})$: only the DC component passes; if $H$ has a zero at DC, the output mean is zero.
- The **PSD** is *multiplied* by $|H|^2$: each frequency is attenuated or amplified independently.
- The **autocorrelation** is $R_X$ convolved with the deterministic autocorrelation $h * h^\dagger$ of the filter — a generalisation of the deterministic formula $y = h * x \Rightarrow r_{yy} = r_{xx} * (h*h^\dagger)$.

The **output variance** is

$$
\sigma_Y^2 = \frac{1}{2\pi}\int_{-\pi}^{\pi} |H(e^{j\omega})|^2 S_X(e^{j\omega})\, d\omega,
$$

the weighted average of $S_X$ over $|H|^2$.

### Worked example 1.8 — white noise through a first-order lowpass

$X[n]$ white of variance $\sigma_w^2$, filtered by $H(z) = \frac{1-a}{1 - a z^{-1}}$ with $|a| < 1$ (DC gain $H(1) = 1$).

**Output PSD.**

$$
S_Y(e^{j\omega}) = \sigma_w^2 \cdot \frac{(1-a)^2}{(1 - a e^{-j\omega})(1 - a e^{j\omega})} = \frac{\sigma_w^2(1-a)^2}{1 - 2a\cos\omega + a^2}.
$$

**Output autocorrelation** (inverse DTFT of $S_Y$):

$$
R_Y[k] = \sigma_w^2\, \frac{1-a}{1+a}\, a^{|k|}.
$$

**Check variance.** $R_Y[0] = \sigma_w^2 \frac{1-a}{1+a}$. For $a = 0$: $R_Y[0] = \sigma_w^2$ ✓ (filter is identity, output = input). For $a \to 1$: $R_Y[0] \to 0$ ✓ (filter becomes an ideal DC selector; white noise has zero DC; output is zero).

**Output mean.** $\mu_Y = \mu_X \cdot H(1) = \mu_X \cdot 1 = \mu_X$. For zero-mean white noise, output is zero-mean.

**Correlation time.** The decay rate of $R_Y[k]$ is $a^{|k|}$; the **integrated autocorrelation time** is $\tau_{int} = \sum_k R_Y[k]/R_Y[0] = (1+a)/(1-a)$. As $a \to 1$, the output is highly correlated; as $a \to 0$, it approaches white noise.

---

## 1.10 Simulation of white noise

To simulate a **white noise sequence** $X[n]$ (IID with a given distribution):

### 1.10.1 Uniform noise

Any PRNG (Mersenne Twister, xorshift, PCG) produces approximately uniform $\mathcal U[0,1)$ samples. In Python: `np.random.rand()`. In Java: `java.util.Random.nextDouble()`. In C#: `System.Random.NextDouble()`.

For $\mathcal U[a, b]$: $X = a + (b-a)U$.

### 1.10.2 Gaussian noise — Box–Muller transform

If $U_1, U_2$ are independent $\mathcal U(0,1)$ samples, then

$$
\boxed{\;Z_1 = \sqrt{-2\ln U_1}\,\cos(2\pi U_2),\qquad Z_2 = \sqrt{-2\ln U_1}\,\sin(2\pi U_2)\;}
$$

are independent $\mathcal N(0,1)$ samples. This is the **Box–Muller transform**.

**Derivation sketch.** Transform $(U_1, U_2)$ to polar $(R, \Theta)$: $R = \sqrt{-2\ln U_1}$ and $\Theta = 2\pi U_2$. The Jacobian gives the joint density

$$
f_{R,\Theta}(r,\theta) = \frac{r}{2\pi} e^{-r^2/2},\quad r \ge 0,\ 0 \le \theta < 2\pi,
$$

which factors into $R \sim \mathrm{Rayleigh}(1)$ and $\Theta \sim \mathcal U[0, 2\pi)$. Then $Z_1 = R\cos\Theta$, $Z_2 = R\sin\Theta$ are the Cartesian coordinates of a 2-D standard Gaussian. ✓

### 1.10.3 Arbitrary distribution — inverse CDF

If $F$ is the CDF of the desired distribution and $U \sim \mathcal U(0,1)$, then

$$
X = F^{-1}(U)
$$

has CDF $F$. This is the **inverse transform sampling** method. Works for any strictly increasing $F$; for distributions with atoms, handle them by intervals of $U$.

### 1.10.4 Sample autocorrelation of simulated white noise

For $N$ samples of white noise with variance $\sigma^2$:

$$
\mathbb{E}[\hat{R}[k]] = \begin{cases} \sigma^2 & k = 0 \\ 0 & k \ne 0 \end{cases},\qquad
\operatorname{Var}[\hat{R}[k]] \approx \frac{\sigma^4}{N}.
$$

So $|\hat{R}[k]| \le 3\sigma^2/\sqrt N$ for $k \ne 0$ with 99% confidence. This is the **white-noise confidence band** — any empirical autocorrelation that stays within $\pm 3\sigma/\sqrt N$ for $k \ne 0$ is consistent with white noise.

### Worked example 1.9 — Box–Muller in Python

```python
import numpy as np

def box_muller(n, seed=0):
    rng = np.random.default_rng(seed)
    U1 = rng.uniform(0, 1, size=n)
    U2 = rng.uniform(0, 1, size=n)
    R = np.sqrt(-2 * np.log(U1))
    Z1 = R * np.cos(2*np.pi*U2)
    Z2 = R * np.sin(2*np.pi*U2)
    return np.concatenate([Z1, Z2])

z = box_muller(100_000)
print(f"mean: {z.mean():.4f}, var: {z.var():.4f}")  # ≈ 0, 1
```

This produces $2n$ samples from a single pair of length-$n$ uniform arrays. Empirical mean $\approx 0$, variance $\approx 1$ — verifying the transform.

---

## 1.11 Lowpass filtering of white noise — colouring

**The point.** White noise is an abstraction — no physical process has infinite bandwidth. Every real signal is **white noise filtered** by some LTI system (the source of the "colour"). Conversely, simulating a coloured process with a given PSD $S_X$ amounts to

1. Generating white noise $W[n]$ of unit variance,
2. Filtering it with the **spectral factor** $H(z)$ from §1.8: $X[n] = h * W[n]$.

This works because $S_X = |H|^2$ when $W$ has unit PSD.

### Worked example 1.10 — AR(1) as filtered white noise

$X[n] = a X[n-1] + W[n]$, $|a| < 1$, $W$ white of variance $\sigma_W^2$. The filter is $H(z) = 1/(1-az^{-1})$; output PSD $S_X = \sigma_W^2/(1-2a\cos\omega + a^2)$.

**Simulation (Python):**

```python
import numpy as np
from scipy.signal import lfilter

def simulate_ar1(a, sigma_w, N, seed=0):
    rng = np.random.default_rng(seed)
    w = rng.normal(0, sigma_w, size=N + 200)     # burn-in
    x = lfilter([1], [1, -a], w)[200:]           # discard transient
    return x

x = simulate_ar1(a=0.9, sigma_w=1.0, N=100_000)
print(f"empirical mean: {x.mean():.4f}")                 # ≈ 0
print(f"empirical variance: {x.var():.4f}")              # ≈ 1/(1-a²) = 5.263
print(f"empirical R_x[1]: {np.corrcoef(x[:-1], x[1:])[0,1]:.4f}")  # ≈ a = 0.9
```

### Worked example 1.11 — ideal lowpass (brick-wall) filter

Ideal lowpass with cutoff $\omega_c$: $|H(e^{j\omega})|^2 = 1$ for $|\omega| \le \omega_c$, zero otherwise.

**Output PSD.** $S_Y(e^{j\omega}) = \sigma_w^2$ for $|\omega| \le \omega_c$, zero otherwise.

**Output autocorrelation.** Inverse DTFT:

$$
R_Y[k] = \frac{1}{2\pi}\int_{-\omega_c}^{\omega_c}\sigma_w^2\, e^{j\omega k}\,d\omega = \sigma_w^2\,\frac{\sin(\omega_c k)}{\pi k}.
$$

**Output variance.** $R_Y[0] = \sigma_w^2 \omega_c/\pi$ — proportional to bandwidth. The **fractional bandwidth** $\omega_c/\pi$ is the fraction of white-noise power that survives.

The ideal lowpass is **not causal** (its impulse response is a sinc), so it can only be approximated in real time. This is why **finite-length FIR approximations** are used in practice, and why the AR(1) filter (causal, recursive, infinite-order) is so popular as a "practical lowpass" — one coefficient controls the cutoff.

### Worked example 1.12 — whitening a coloured process

Given a coloured process $Y[n]$ with PSD $S_Y = \sigma_e^2 |H|^2$, the **whitening filter** is $G(z) = 1/H(z)$. Its output has PSD

$$
S_{Y_{\text{white}}}(e^{j\omega}) = |G(e^{j\omega})|^2 S_Y(e^{j\omega}) = \sigma_e^2
$$

— flat, i.e., white noise of variance $\sigma_e^2$. This is the **inverse filter** of the spectral factor.

For AR(1) with $a = 0.9$: $G(z) = 1 - 0.9 z^{-1}$, i.e., $Y_{white}[n] = Y[n] - 0.9 Y[n-1]$. The output is white with variance $\sigma_e^2 = \sigma_W^2$. **Practical implication:** the first-order difference $Y[n] - 0.9 Y[n-1]$ is a whitening filter — a fact used throughout speech processing (where speech is approximated as AR) and Kalman filtering (where the innovations are the whitened observations).

---

## 1.12 Chapter summary

1. **Discrete random processes** are families of random variables indexed by integer $n$. Wide-sense stationarity means constant mean and autocorrelation depending only on lag.
2. **Parseval's theorem** — $\sum_n |x[n]|^2 = \frac{1}{2\pi}\int_{-\pi}^{\pi}|X(e^{j\omega})|^2 d\omega$ — extends to WSS processes with $R_X[0] = \frac{1}{2\pi}\int S_X\,d\omega$.
3. **Wiener–Khinchine** makes the PSD the DTFT of the autocorrelation; it is non-negative, $2\pi$-periodic, and integrates to $R_X[0]$.
4. **The periodogram** is the natural estimator; asymptotically unbiased but inconsistent ($\mathrm{Var} \approx S_X^2$).
5. **Sample autocorrelation** — use the biased form to guarantee a valid PSD; use the unbiased form only for descriptive statistics.
6. **Sum Decomposition** = eigen-decomposition of the covariance matrix; the KL / PCA expansion in finite dimension.
7. **Spectral Factorization** — every rational PSD with $\int\log S > -\infty$ has a unique causal, minimum-phase factor $S_X = \sigma_e^2|H|^2$; the innovations variance $\sigma_e^2$ is the geometric mean of $S_X$.
8. **Linear systems** — $S_Y = |H|^2 S_X$; output autocorrelation is $R_X * (h*h^\dagger)$; the output mean is $\mu_X H(1)$.
9. **White noise simulation** — uniform from any PRNG, Gaussian from Box–Muller, arbitrary from inverse CDF.
10. **Coloured noise simulation** — filter white noise with the spectral factor; **whitening** is the inverse filter.

---

## 1.13 Quiz

### A. Multiple choice (choose the single best answer)

1. A WSS sequence has:
   (a) constant variance, arbitrary mean (b) constant mean, autocorrelation depending only on lag (c) Gaussian distribution (d) finite energy

2. The energy of a sequence $x[n]$ is $\sum |x[n]|^2$. Parseval's theorem for finite-length sequences with the $N$-point DFT expresses this as:
   (a) $\sum |X[k]|^2$ (b) $\frac{1}{N}\sum |X[k]|^2$ (c) $N\sum |X[k]|^2$ (d) $\frac{1}{N^2}\sum |X[k]|^2$

3. The Wiener–Khinchine relation states that:
   (a) the PSD equals the DTFT of the autocorrelation (b) the autocorrelation is the DFT of the signal (c) the signal is the PSD (d) the noise is white

4. For a WSS process, $R_X[0]$ equals:
   (a) the mean (b) the variance (c) the total power (variance plus mean squared) (d) zero

5. The **biased** sample autocorrelation is preferred over the **unbiased** form because:
   (a) it is unbiased (b) it has lower variance and is positive-definite (c) it is easier to compute (d) it is a better estimator of the mean

6. The **periodogram** of white noise:
   (a) is a flat line (b) has mean $S_X$ but standard deviation equal to the mean (c) is Gaussian (d) is zero at all frequencies

7. The Sum Decomposition Theorem is equivalent to:
   (a) the FFT (b) the eigen-decomposition of the covariance matrix (c) the Z-transform (d) the DFT

8. The Spectral Factorization Theorem says every rational PSD with $\int\log S > -\infty$ has a unique causal, minimum-phase factor $H(z)$ with:
   (a) $S_X = |H|^2$ (b) $S_X = \sigma_e^2 |H|^2$ (c) $S_X = H$ (d) $S_X = 1/H$

9. For white noise through an LTI filter with response $H$:
   (a) output PSD = input PSD (b) output PSD = $|H|^2 \times$ input PSD (c) output is white (d) output is non-stationary

10. The Box–Muller transform generates:
    (a) uniform samples (b) exponential samples (c) Gaussian samples (d) Poisson samples

11. The **output variance** of white noise of variance $\sigma_w^2$ through a filter with impulse response $h$ is:
    (a) $\sigma_w^2$ (b) $\sigma_w^2 \sum_k |h[k]|^2$ (c) $\sigma_w^2 \sum_k h[k]$ (d) $0$

12. The **whitening filter** for a process with spectral factor $H(z)$ is:
    (a) $H(z)$ (b) $H(z^{-1})$ (c) $1/H(z)$ (d) $H(z) H(z^{-1})$

### B. Fill in the blank

1. A WSS process has constant mean and autocorrelation depending only on \_\_\_\_.
2. $R_X[0]$ equals $\sigma_X^2 + $\_\_\_\_ for a WSS process.
3. Parseval's theorem for finite-length sequences is $\sum|x[n]|^2 = \frac{1}{N}\sum |X[k]|^2$ (with the $N$-point DFT); this holds because the DFT is \_\_\_\_-preserving.
4. The Wiener–Khinchine relation is $S_X(e^{j\omega}) = \sum_{k} $\_\_\_\_ $\,e^{-j\omega k}$.
5. The periodogram is $\hat{S}_X = \frac{1}{N}\bigl|\sum_{n}$\_\_\_\_$\bigr|^2$.
6. The Sum Decomposition Theorem expresses a random vector as $\mathbf{X} = \sum_k Z_k \mathbf{v}_k$ with $\mathbb{E}[Z_i Z_j] = $\_\_\_\_ $\delta_{ij}$.
7. The Spectral Factorization Theorem gives $S_X(z) = $\_\_\_\_ $H(z) H(z^{-1})$ with $H$ minimum-phase.
8. For white noise of variance $\sigma_w^2$ through a filter with impulse response $h$, the output variance is $\sigma_w^2 \sum_k$ \_\_\_\_.
9. The whitening filter for a process with spectral factor $H$ is \_\_\_\_.
10. Box–Muller transforms two uniform samples into two \_\_\_\_ samples.

### C. Define / explain

1. Define a **discrete-time random process** and state the conditions for **wide-sense stationarity**. Give one example that is WSS but not strictly stationary.
2. State **Parseval's theorem** for finite-length sequences and for WSS processes. What is the interpretation in each case?
3. State the **Wiener–Khinchine relation**. What is the physical interpretation of the PSD's non-negativity?
4. Define the **sample autocorrelation** in both biased and unbiased forms. Why is the biased form preferred in spectral estimation?
5. Define the **periodogram** and discuss its bias and variance as $N \to \infty$. Why is it inconsistent?
6. State the **Sum Decomposition Theorem** and interpret it as a data-compression scheme. What is the truncated approximation error?
7. State the **Spectral Factorization Theorem**. Why is the minimum-phase factor unique and important?
8. For WSS input to an LTI system with response $H$, derive the output PSD, output autocorrelation, and output mean.
9. Describe the **Box–Muller transform** and state what it is used for. Give the two formulas for $Z_1, Z_2$.
10. Explain what it means to **colour** white noise, and how the **spectral factor** is used in simulation.

---

## 1.14 Answer key

??? success "A. Multiple choice"
    1. **(b)** — constant mean + autocorrelation depends only on lag.
    2. **(b)** — $\frac{1}{N}\sum|X[k]|^2$.
    3. **(a)** — PSD = DTFT of autocorrelation.
    4. **(c)** — $R_X[0] = \sigma_X^2 + |\mu_X|^2$ (total power).
    5. **(b)** — lower variance, positive-definite.
    6. **(b)** — mean = $S_X$, std = mean.
    7. **(b)** — eigen-decomposition of covariance.
    8. **(b)** — $S_X = \sigma_e^2 |H|^2$.
    9. **(b)** — $|H|^2 \times$ input PSD.
    10. **(c)** — Gaussian samples.
    11. **(b)** — $\sigma_w^2 \sum |h[k]|^2$.
    12. **(c)** — $1/H(z)$.

??? success "B. Fill in the blank"
    1. lag
    2. $|\mu_X|^2$ (or mean squared)
    3. energy
    4. $R_X[k]$
    5. $x[n] e^{-j\omega n}$
    6. $\lambda_i$ (eigenvalues)
    7. $\sigma_e^2$
    8. $|h[k]|^2$
    9. $1/H(z)$ (or $H^{-1}(z)$)
    10. Gaussian (independent standard normal)

??? success "C. Define / explain — model answers"
    1. *Discrete-time random process:* family $\{X[n]\}_{n\in\mathbb Z}$ of random variables on a common probability space. *WSS:* (i) $\mu_X[n] = \mu_X$ constant; (ii) $R_X[m,n] = \mathbb{E}[X[m]X^*[n]]$ depends only on $m - n$. *Example WSS but not SSS:* $X[n] = A\cos(\omega_0 n + \Theta)$ with $A$ a non-constant random amplitude (e.g., $A \sim$ Rayleigh) independent of $\Theta \sim \mathcal U(0,2\pi)$. Mean is zero, $R_X[k] = \frac{\mathbb{E}[A^2]}{2}\cos(\omega_0 k)$ — WSS. But $\mathbb{E}[X[n]^3] = \frac{\mathbb{E}[A^3]}{4}\cos(3\omega_0 n + \cdots)$ depends on $n$, so the third-order moments are time-varying — not strictly stationary.
    2. *Finite-length:* $\sum_{n=0}^{N-1}|x[n]|^2 = \frac{1}{N}\sum_{k=0}^{N-1}|X[k]|^2$ where $X[k]$ is the $N$-point DFT. Interpretation: energy in time = energy in frequency (the DFT is a unitary transform up to the factor $N$). *WSS:* $R_X[0] = \frac{1}{2\pi}\int_{-\pi}^{\pi}S_X(e^{j\omega})d\omega$. Interpretation: total power = integral of PSD.
    3. $S_X(e^{j\omega}) = \sum_k R_X[k]e^{-j\omega k}$. *Non-negativity interpretation:* $S_X \ge 0$ means that the autocorrelation is positive-definite, which reflects that power at any frequency cannot be negative — a physical necessity. Equivalently, the PSD is a "weight" for the spectral content of the process, so its negativity would imply negative energy at some frequency, which is impossible.
    4. *Biased:* $\hat R[k] = \frac{1}{N}\sum_{n=0}^{N-1-|k|} x[n+|k|]x^*[n]$. *Unbiased:* $\hat R_u[k] = \frac{1}{N-|k|}\sum x[n+|k|]x^*[n]$. *Why biased:* (i) it is positive-definite — its DTFT is a valid PSD; (ii) it has lower mean-square error at large lags (the variance of $\hat R_u$ blows up as $|k| \to N$); (iii) its $O(1/N)$ bias vanishes asymptotically. For **spectral estimation** it is therefore the standard choice.
    5. $\hat S_X(e^{j\omega}) = \frac{1}{N}|\sum x[n]e^{-j\omega n}|^2$. *Bias:* $\mathbb{E}[\hat S_X] \to S_X$ as $N \to \infty$, with $O(1/N)$ bias for finite $N$. *Variance:* $\operatorname{Var}[\hat S_X] \approx S_X^2$ for all $N$ — independent of record length. *Why inconsistent:* the variance does not shrink to zero, so increasing $N$ at fixed realization does not improve accuracy. The fix is averaging over independent **segments** (Bartlett, Welch), or over independent realizations.
    6. Sum Decomposition: $\mathbf{R} = \sum_k \lambda_k \mathbf{v}_k \mathbf{v}_k^\top$; the covariance matrix eigen-decomposition. *Data-compression interpretation:* the truncated expansion $\hat{\mathbf{X}}_M = \sum_{k=1}^M Z_k \mathbf{v}_k$ retains the $M$ largest-eigenvalue directions and discards the rest. *Truncation error:* $\mathbb{E}[\|\mathbf{X} - \hat{\mathbf{X}}_M\|^2] = \sum_{k=M+1}^N \lambda_k$. This is the optimal rank-$M$ linear compression in mean-square sense — the theoretical basis of PCA / KL / DCT-based compression.
    7. *Statement:* for a rational PSD positive on $|z|=1$ with $\int\log S > -\infty$, there exists a unique factorization $S_X(z) = \sigma_e^2 H(z)H(z^{-1})$ with $H$ causal, monic, minimum-phase. *Uniqueness:* among all factors with the same magnitude response, exactly one is minimum-phase (all zeros inside the unit circle). *Importance:* (i) $H^{-1}$ is causal — the whitening filter is realizable; (ii) it gives the smallest innovations variance $\sigma_e^2$ — the most efficient representation; (iii) it is the basis of every AR/ARMA model and every innovations representation.
    8. *Output PSD:* $S_Y(e^{j\omega}) = |H(e^{j\omega})|^2 S_X(e^{j\omega})$. Derivation: $\mathbb{E}[Y[n+k]Y^*[n]] = \sum_{i,j} h[i]h^*[j] R_X[k + i - j]$; taking the DTFT and reordering gives the product form. *Output autocorrelation:* $R_Y[k] = R_X[k] * (h*h^\dagger)[k]$ — convolution of input autocorrelation with the filter's deterministic autocorrelation. *Output mean:* $\mu_Y = \mu_X H(e^{j0}) = \mu_X \sum_k h[k]$ — the filter's DC gain.
    9. *Box–Muller:* generates two independent $\mathcal N(0,1)$ samples from two independent $\mathcal U(0,1)$ samples. $Z_1 = \sqrt{-2\ln U_1}\cos(2\pi U_2)$, $Z_2 = \sqrt{-2\ln U_1}\sin(2\pi U_2)$. *Use:* simulating Gaussian noise when only a uniform PRNG is available; also gives insight into the polar decomposition of the 2-D Gaussian.
    10. *Colouring* = passing white noise through a filter so the output has the desired PSD. If the target PSD is $S_X = \sigma_w^2 |H|^2$, filtering unit-variance white noise through $H(z)$ gives output with PSD $S_X$. The spectral factor $H$ is obtained from §1.8. The **simulation recipe:** (i) generate white noise of unit variance; (ii) apply $H$ via `scipy.signal.lfilter` or its Java/C# equivalent; (iii) the output is the simulated process. Conversely, **whitening** applies $H^{-1}$ to remove the correlation and produce white noise.

---

## Further reading

- Hayes, *Statistical Digital Signal Processing and Modeling*, Wiley, 1996 — Chapters 1–3, the exact match to this chapter.
- Proakis & Manolakis, *Digital Signal Processing*, Pearson, 2007 — Chapter 2 for review; Chapter 14 for spectral estimation.
- Oppenheim & Schafer, *Discrete-Time Signal Processing*, PHI — Chapters 2, 4, 10 for the deterministic foundations.
- Kay, *Modern Spectral Estimation*, Prentice-Hall, 1988 — the reference for periodogram analysis and the Sum Decomposition.
- Papoulis & Pillai, *Probability, Random Variables, and Stochastic Processes*, 4th ed. — Chapters 9–10 for the statistical background.
- Stoica & Moses, *Spectral Analysis of Signals*, Prentice-Hall, 2005 — for the general spectral factorization theory.

---

## Suggested tutorial and lab problems for this chapter

1. **Parseval check.** For $x[n] = a^n u[n]$ with $|a| < 1$, verify Parseval's theorem for infinite-length energy signals: $\sum |x[n]|^2 = \frac{1}{2\pi}\int |X(e^{j\omega})|^2 d\omega$.
2. **Autocorrelation derivation.** Derive $R_X[k]$ for $X[n] = A\cos(\omega_0 n + \Theta)$ with $\Theta$ uniform and $A$ fixed. Compare with the continuous-time result.
3. **Wiener–Khinchine inverse.** Given $S_X(e^{j\omega}) = \frac{1 - a^2}{|1 - a e^{-j\omega}|^2}$ with $|a| < 1$, compute $R_X[k]$ by inverse DTFT. Verify the answer is $\rho^{|k|}$ with $\rho = a$.
4. **Periodogram variance study.** Simulate $N = 256$ samples of white Gaussian noise, compute the periodogram, and plot its histogram over frequencies. Verify that $\hat S_X/\sigma^2 \sim \mathrm{Exp}(1)$. Repeat for $N = 4096$ and $N = 65536$ — does the variance change?
5. **Sum decomposition.** For an AR(1) with $\rho = 0.9$ and $N = 16$, compute the eigen-decomposition of the covariance matrix numerically. Plot the eigenvalues. How many components capture 99% of the variance? Compare with $\rho = 0.5$.
6. **Spectral factorization exercise.** Factorize $S_X(z) = \frac{(1 - 0.3z^{-1})(1 - 0.3z)}{(1 - 0.9z^{-1})(1 - 0.9z)}$ into a causal, minimum-phase $H(z)$ with $\sigma_e^2$. Compute the whitening filter $1/H$.
7. **Filtering of white noise.** White noise of variance $1$ passes through $H(z) = 1 + 0.5 z^{-1}$. Compute the output PSD, autocorrelation, and variance. Verify by simulation.
8. **Box–Muller statistics.** Generate 100,000 Gaussian samples using Box–Muller. Verify the empirical mean, variance, skewness, and kurtosis against the theoretical $\mathcal N(0,1)$ values. Plot the histogram and compare to the Gaussian PDF.
9. **Colouring exercise.** Generate white noise, apply an AR(2) filter with $a_1 = 1.3$, $a_2 = -0.42$. Plot the empirical PSD (using Welch) and compare against the theoretical $S_X = \sigma_w^2/|1 - 1.3e^{-j\omega} + 0.42e^{-2j\omega}|^2$.
10. **Whitening speech.** Download a segment of voiced speech (or use `pyroomacoustics` sample). Fit an AR(10) model, apply the whitening filter, and verify the residual is close to white (its autocorrelation outside lag 0 is within $\pm 3/\sqrt N$). Interpret the result in terms of the **source-filter** model of speech.