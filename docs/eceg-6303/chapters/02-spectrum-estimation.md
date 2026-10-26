# Chapter 2 — Spectrum Estimation

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** nonparametric methods (correlation method,
    covariance estimator, unbiased/consistent estimators, periodogram, Bartlett,
    Welch) · performance analysis of estimators · model-based approach (AR, MA,
    ARMA signal modelling) · parameter estimation using the Yule–Walker method.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Explain the difference between **nonparametric** and **model-based** spectrum
   estimation, and choose between them based on data length, prior knowledge,
   and resolution requirements.
2. State the **bias, variance, MSE, and consistency** of an estimator, and prove
   that the periodogram is asymptotically unbiased but **inconsistent**.
3. Compute and interpret the **Bartlett** and **Welch** spectrum estimators, and
   quantify their variance reduction relative to the raw periodogram.
4. Apply the **correlation method** (Blackman–Tukey) using a windowed
   autocorrelation estimate, and analyse the trade-off between bias (window
   choice) and variance (window length).
5. Distinguish the **autocorrelation** and **covariance** variants of the
   correlation method, and explain when each is preferred.
6. Model a WSS process as **AR, MA, or ARMA**; estimate parameters via the
   **Yule–Walker** equations; solve them efficiently using **Levinson–Durbin**.
7. Compute the **AR / MA / ARMA PSD** from the estimated model.
8. Select **model order** using AIC, MDL, and FPE; interpret the trade-off
   between fit and complexity.
9. Implement all of the above in **Python, Java, and C#**.

---

## 2.1 The spectrum-estimation problem

Given $N$ samples $x[0], \dots, x[N-1]$ of a WSS process $X[n]$, estimate the
**power spectral density** $S_X(e^{j\omega})$.

Two broad families:

| Family | Approach | Assumptions | Good when |
|---|---|---|---|
| **Nonparametric** | Directly estimate $S_X$ from data using DFT, autocorrelation, or windows | Only WSS | Little prior knowledge; long records; spectra without sharp peaks |
| **Model-based (parametric)** | Fit a model (AR, MA, ARMA) with $p$ parameters; compute $S_X$ analytically from the model | Model class is correct | Short records; sharp peaks; prior knowledge of process type |

**Two competing goals.**

- **Resolution.** The ability to distinguish two close spectral peaks. Nonparametric resolution is $\sim 1/N$ in normalised frequency. Model-based can exceed this if the model is correct.
- **Variance.** The statistical fluctuation of the estimator around its mean. Nonparametric estimators need **averaging** (Bartlett, Welch) to reduce variance, which trades resolution.

There is no universally best estimator. The choice depends on **data length**, **prior knowledge**, and **what you want to know** (peak locations? total power? smooth envelope?).

---

## 2.2 Performance criteria for estimators

An estimator $\hat\theta$ of a quantity $\theta$ (here $\theta = S_X(e^{j\omega})$ at a
given frequency) is characterised by:

$$
\boxed{\;\operatorname{Bias}[\hat\theta] = \mathbb{E}[\hat\theta] - \theta\;}
$$

$$
\boxed{\;\operatorname{Var}[\hat\theta] = \mathbb{E}\bigl[(\hat\theta - \mathbb{E}[\hat\theta])^2\bigr]\;}
$$

$$
\boxed{\;\operatorname{MSE}[\hat\theta] = \operatorname{Bias}^2 + \operatorname{Var} = \mathbb{E}\bigl[(\hat\theta - \theta)^2\bigr]\;}
$$

!!! abstract "Definitions"
    - **Unbiased:** $\operatorname{Bias} = 0$ for every $N$.
    - **Asymptotically unbiased:** $\operatorname{Bias} \to 0$ as $N \to \infty$.
    - **Consistent:** $\hat\theta \to \theta$ in probability as $N \to \infty$ —
      equivalently, **both** bias $\to 0$ and variance $\to 0$.
    - **Efficient:** among all unbiased estimators, has the smallest variance
      (the Cramér–Rao lower bound).

The key insight for spectrum estimation: **the periodogram is asymptotically
unbiased but NOT consistent** — its variance does not decrease with $N$. This is
the central motivation for Bartlett and Welch averaging, and for model-based
approaches.

### Worked example 2.1 — checking consistency of a simple estimator

Let $x[n] = \theta + e[n]$ with $e[n]$ zero-mean white noise of variance
$\sigma_e^2$, $\theta$ constant. Consider the mean estimator
$\hat\theta_N = \frac{1}{N}\sum_{n=0}^{N-1} x[n]$.

- **Bias.** $\mathbb{E}[\hat\theta_N] = \theta$, so bias $= 0$ — unbiased.
- **Variance.** $\operatorname{Var}[\hat\theta_N] = \sigma_e^2/N \to 0$.
- **Consistency.** Both bias and variance vanish ⇒ **consistent**. ✓

Now consider $\hat\theta_N' = x[0]$ — a "one-sample" estimator.

- Bias $= 0$, but variance $= \sigma_e^2$ — **does not** vanish.

So $\hat\theta_N'$ is unbiased but **inconsistent**. Consistency is not implied by
unbiasedness — it requires the variance to also go to zero.

**Application to spectrum.** The periodogram is the analogue of $\hat\theta_N'$:
unbiased (asymptotically) but with variance that does **not** shrink with $N$.
Bartlett and Welch are the analogues of $\hat\theta_N$ — they reduce variance
at the cost of a small bias.

---

## 2.3 The periodogram — bias and variance revisited

From §1.6 of Chapter 1, the periodogram is

$$
\hat{S}_X(e^{j\omega}) = \frac{1}{N}\Bigl|\sum_{n=0}^{N-1} x[n]\, e^{-j\omega n}\Bigr|^2.
$$

### 2.3.1 Bias

Taking the expectation over realisations:

$$
\mathbb{E}\bigl[\hat{S}_X(e^{j\omega})\bigr] = \frac{1}{2\pi}\int_{-\pi}^{\pi} S_X(e^{j\theta})\, W_N(\omega - \theta)\, d\theta,
$$

where the **spectral window** is the squared Dirichlet (Fejér) kernel:

$$
\boxed{\;W_N(\omega) = \frac{1}{N}\cdot\frac{\sin^2(N\omega/2)}{\sin^2(\omega/2)}\;}
$$

**Interpretation.** The expected periodogram is the true PSD **smoothed** by
$W_N$ — a window of effective width $\sim 2\pi/N$. This gives three consequences:

1. The periodogram is **asymptotically unbiased**: $W_N \to 2\pi\delta(\omega)$ as $N \to \infty$, so the smoothing vanishes.
2. For finite $N$, sharp spectral features (narrow peaks, line spectra) are **spread** over $2\pi/N$ — the *spectral leakage* and the **resolution limit**.
3. The smoothing is a **bias** in the MSE sense — it distorts the shape even where it does not change the mean of the peak value.

### 2.3.2 Variance

For a **Gaussian** WSS process, the periodogram at any single frequency is
asymptotically exponentially distributed with mean $S_X$:

$$
\boxed{\;\frac{\hat{S}_X(e^{j\omega})}{S_X(e^{j\omega})} \xrightarrow{\text{distribution}} \frac{\chi^2_2}{2} = \mathrm{Exp}(1)\;}
$$

Hence

$$
\boxed{\;\operatorname{Var}\bigl[\hat{S}_X(e^{j\omega})\bigr] \approx S_X^2(e^{j\omega})\quad \text{for all } N.\;}
$$

The standard deviation equals the mean, regardless of record length.
This is the **inconsistency** of the periodogram — increasing $N$ reduces the bias
but never the variance at any given frequency.

### Worked example 2.2 — verifying periodogram variance

Simulate $N = 4096$ samples of white Gaussian noise with $\sigma^2 = 1$;
compute the periodogram at $N/2 + 1 = 2049$ frequency bins; report mean and
standard deviation across bins.

```
mean(periodogram) = 1.00  (expected: σ² = 1)
std(periodogram)  = 1.00  (expected: σ² = 1)
```

Double $N$ to $8192$: mean stays at $1.00$, **std stays at $1.00$**.
This is what inconsistent looks like — variance does not shrink as $N$
increases. Averaging over $K$ **independent** realisations would reduce
the variance to $\sigma^4/K$.

---

## 2.4 Bartlett's method — averaging segments

### Idea

Split the $N$-sample record into $K$ **non-overlapping** segments of length
$L = \lfloor N/K \rfloor$ each. Compute the periodogram of each segment and
**average**:

$$
\boxed{\;\hat{S}^{\text{Bartlett}}_X(e^{j\omega}) = \frac{1}{K}\sum_{k=1}^{K}\hat{S}^{(k)}_X(e^{j\omega})\;}
$$

where $\hat{S}^{(k)}_X$ is the periodogram of segment $k$ using $L$ samples.

### Bias / variance trade-off

- **Bias.** Each segment's periodogram has the same asymptotic bias as a length-$L$ periodogram — i.e., resolution $\sim 2\pi/L = K \cdot 2\pi/N$. **The bias increases** by the factor $K$ compared to the full-record periodogram.
- **Variance.** If the segments are **independent** (or approximately so), the average of $K$ independent periodograms has $\operatorname{Var} = S_X^2/K$ — a $1/K$ reduction.

**Formally:**

$$
\operatorname{Bias}[\hat S^{\text{Bartlett}}] = O(L^{-1}) = O(K/N),\qquad
\operatorname{Var}[\hat S^{\text{Bartlett}}] = O(1/K).
$$

As $N \to \infty$ with $K \to \infty$ and $L \to \infty$ (so that both $K, L \to \infty$), both bias and variance vanish ⇒ **Bartlett is consistent**.

### Worked example 2.3 — Bartlett variance reduction

$N = 1024$ samples, $K = 4$ segments of length $L = 256$ each.

| Quantity | Raw periodogram | Bartlett (K = 4) |
|---|---|---|
| Resolution | $2\pi/1024$ | $2\pi/256 = 4\times$ coarser |
| Variance | $S_X^2$ | $S_X^2/4$ (÷4) |
| Standard deviation | $S_X$ | $S_X/2$ |

Trade-off: 4× worse resolution, 2× lower standard deviation. For a smooth PSD
(no narrow peaks), Bartlett is clearly better. For a narrow line, the bias
dominates — the peak is spread over a wider band.

!!! tip "Choosing $K$"
    Rule of thumb: **$K = 4$ to $8$**, giving variance reduction factors of
    $4$ to $8$. Much larger $K$ improves variance but degrades resolution so
    badly that narrow peaks become unrecognisable. Choose $K$ based on what you
    care about — if the goal is the *envelope* of a smooth spectrum, larger $K$
    is fine; if it is *peak location*, use the largest $L$ you can afford
    (which means small $K$ or Welch's overlap trick).

---

## 2.5 Welch's method — overlap and windowing

### Idea

Welch's method improves on Bartlett by:

1. **Overlap** between segments (typically 50%): more segments → lower variance
   for the same total record.
2. **Window** each segment (Hann, Hamming, or Blackman) before computing the
   periodogram: reduces **spectral leakage** from strong out-of-band components.

### The estimator

For each of $K$ segments of length $L$, overlapping by $L/2$:

$$
\boxed{\;\hat{S}^{\text{Welch}}_X(e^{j\omega}) = \frac{1}{K U}\sum_{k=1}^{K}\Bigl|\sum_{n=0}^{L-1} w[n]\, x[n + k L/2]\, e^{-j\omega n}\Bigr|^2\;}
$$

where $w[n]$ is the analysis window and $U = \frac{1}{L}\sum_n |w[n]|^2$ is a
normalisation constant that removes the window's effect on the average level.

### Bias / variance

- **Bias.** Each segmental periodogram has resolution $\sim 2\pi/L$, as in
  Bartlett, plus a small additional bias from the window's mainlobe (which
  broadens the effective smoothing by the window's coherent gain).
- **Variance.** Reduced by the number of **effectively independent** segments.
  With 50% overlap, the segments are **correlated** (they share half their
  samples), so the effective reduction is less than $1/K$. For a Hann window
  with 50% overlap, the effective number is roughly **$0.75K$**, so
  $\operatorname{Var} \approx S_X^2/(0.75K)$.

**Net effect.** Compared to Bartlett with the same total record:

- Same resolution (segment length $L$ is the same).
- **Lower variance** (more segments from overlap).
- **Reduced leakage** (window tapers the segment edges).

Welch is the default practical choice for nonparametric spectrum estimation,
and is what `scipy.signal.welch`, MATLAB's `pwelch`, and every spectrum analyser
implements.

### Worked example 2.4 — Welch vs. Bartlett on a signal with a weak tone

Signal: $x[n] = \cos(2\pi f_0 n) + w[n]$ where $f_0 = 0.2$ (normalised frequency,
inside $[-0.5, 0.5]$), $w$ white Gaussian with variance $1$ (SNR $= 0$ dB),
$N = 4096$ samples.

| Method | Parameters | Peak SNR vs. noise floor | Resolution | Variance |
|---|---|---|---|---|
| Periodogram | — | moderate | best | highest |
| Bartlett | $K=8$, $L=512$ | moderate | 8× worse | ÷8 |
| Welch (rect) | $K=15$, $L=512$, 50% overlap | moderate | 8× worse | ÷~10 |
| Welch (Hann) | $K=15$, $L=512$, 50% overlap | **best** | 8× worse | ÷~11 |

The **Hann window** reduces leakage, so the tone is confined to fewer bins —
the peak is more prominent. The **overlap** gives lower variance than Bartlett
at the same resolution.

The **periodogram** has the best resolution and the peak location is exact,
but the noise floor fluctuates so wildly that a **weak** second tone would be
hidden. Bartlett/Welch reveal weak tones because the noise floor is smoothed.

---

## 2.6 The correlation method (Blackman–Tukey)

### Idea

Estimate the autocorrelation $\hat{R}_X[k]$ from data, apply a **lag window**
$w[k]$, and take the DTFT:

$$
\boxed{\;\hat S^{\text{BT}}_X(e^{j\omega}) = \sum_{k=-(M-1)}^{M-1} w[k]\, \hat{R}_X[k]\, e^{-j\omega k}\;}
$$

where $\hat R_X[k]$ is the biased sample autocorrelation from Chapter 1 (§1.4.3),
and $w[k]$ is a window of half-length $M \le N$. This is the **Blackman–Tukey
(BT) estimator**.

### Why window the autocorrelation?

Without the lag window, $\hat S^{BT}_X$ reduces to the periodogram (a classical
identity: the DTFT of the biased sample autocorrelation equals the periodogram).
The window $w[k]$ **smooths** the periodogram — precisely what Bartlett and
Welch do by averaging, but in the lag domain instead of the time domain.

**Effect of the lag window:**

- **Narrow $w$ (small $M$)** → heavily smoothed spectrum (low variance, high bias/resolution loss).
- **Wide $w$ (large $M$)** → less smoothing (high variance, low bias) → approaches the periodogram as $M \to N$.

### Bias / variance trade-off

The bias depends on the **lag window's** DTFT (called the **spectral window**),
and the variance on the **data window's** equivalent (the **equivalent
autocorrelation window** $w^2[k]$). For the standard windows (Bartlett, Hann,
Hamming, Blackman), the bias-variance trade-off is controlled by $M$ alone:

$$
\operatorname{Bias} \sim \frac{1}{M},\qquad
\operatorname{Var} \sim \frac{M}{N}.
$$

**MSE is minimised** at $M \sim N^{1/3}$ for standard windows — a slow
optimum that makes the choice of $M$ less critical than it looks.

### Worked example 2.5 — Blackman–Tukey with a Bartlett window

Signal $x[n] = \cos(2\pi \cdot 0.2 n) + w[n]$, $w$ white, $N = 1024$.

Choose the **triangular (Bartlett) lag window**:

$$
w[k] = \begin{cases} 1 - |k|/M & |k| \le M \\ 0 & |k| > M \end{cases}
$$

with $M = 64$.

- **Effective smoothing.** The spectral window is the Fejér (squared Dirichlet) kernel of width $\sim 1/M$ — comparable to Bartlett with $L \approx 64$.
- **Variance.** Roughly $\operatorname{Var} \approx S_X^2 \cdot (M/N) \approx S_X^2 \cdot 64/1024 = S_X^2/16$.
- **Peak.** Broadened to $\sim 1/64$ in normalised frequency — a width of ~3 frequency bins of $1/N$.

Sweeping $M \in \{16, 32, 64, 128, 256\}$ traces out the bias-variance
trade-off; $M \approx 100$ is near the MSE optimum for $N = 1024$.

---

## 2.7 The covariance method (unbiased autocorrelation)

The **covariance method** (in the spectrum-estimation sense; distinct from the
covariance method for AR parameter estimation, §2.11) uses the **unbiased**
sample autocorrelation rather than the biased one:

$$
\hat R_u[k] = \frac{1}{N - |k|}\sum_{n=0}^{N-1-|k|} x[n+|k|]\, x^*[n],
\qquad |k| \le M.
$$

The estimator $\hat S_X^{cov}$ = DTFT of a windowed $\hat R_u[k]$.

### Why / why not

| Method | Bias | Variance | Positive-definite? |
|---|---|---|---|
| **Autocorrelation** (biased $\hat R$) | $O(1/N)$ | $O(M/N)$ | **Yes** — valid PSD |
| **Covariance** (unbiased $\hat R_u$) | $O(1/M)$ | $O(M/(N-M))$ | **No** — may go negative |

The covariance method has **lower bias** at each lag (unbiased for each $k$).
But at large $M$ (close to $N$) the variance blows up because the denominator
$N - |k|$ becomes small. More importantly, $\hat S_X^{cov}$ can go **negative**
at some frequencies — a mathematically impossible PSD — because the unbiased
estimate is not positive-definite.

**Practical use:** the covariance method is often preferred when the **autocorrelation
structure is the primary interest** (and you want unbiased lag estimates), not
when a valid PSD is required. For spectral estimation proper, the **biased**
form is the safer default.

### Worked example 2.6 — biased vs. unbiased autocorrelation

White noise of variance $\sigma^2$, $N = 100$.

| Lag $k$ | Biased $\mathbb E[\hat R[k]]$ | Unbiased $\mathbb E[\hat R_u[k]]$ |
|---|---|---|
| 0 | $\sigma^2$ | $\sigma^2$ |
| 10 | $\sigma^2 \cdot 90/100$ | $\sigma^2$ |
| 50 | $\sigma^2 \cdot 50/100$ | $\sigma^2$ |
| 90 | $\sigma^2 \cdot 10/100$ | $\sigma^2$ |

The biased form *deliberately* underestimates — the taper (a triangle window of
width $N$) smoothly rolls off the far-lag autocorrelation, which is exactly what
prevents $\hat S$ from going negative. This is why the biased estimator is
positive-definite: its triangular taper acts as a Bartlett lag window.

---

## 2.8 Comparing the nonparametric methods

| Method | Resolution | Variance | Leakage | Best for |
|---|---|---|---|---|
| Periodogram | $2\pi/N$ | $S_X^2$ (inconsistent) | worst | peak location, very long records |
| Bartlett | $2\pi/L = 2\pi K/N$ | $S_X^2/K$ | reduced | smooth envelopes, medium records |
| Welch | $2\pi/L$ | $S_X^2/(K\cdot\kappa)$ | **lowest** (window) | general purpose; any smooth spectrum |
| Blackman–Tukey | $\sim 2\pi/M$ | $\sim S_X^2 \cdot M/N$ | depends on $w$ | explicit bias-variance control via $M$ |
| Covariance method | $\sim 2\pi/M$ | larger than BT | depends on $w$ | unbiased lag estimates |

$K$ = segments; $L$ = segment length; $M$ = lag-window half-length;
$\kappa \approx 0.75$ for 50% overlap with Hann window.

**Rule of thumb.**

- **Welch with Hann window, 50% overlap, $K \approx 8$** — the default "first try" for any new signal.
- **Periodogram** — only when $N$ is huge ($> 10^5$) and the variance is tolerable.
- **Blackman–Tukey** — when you explicitly want to trade resolution for smoothness and can tune $M$ by hand.
- **Covariance method** — when lag estimates are the primary output (not the spectrum).

---

## 2.9 Model-based spectrum estimation — AR, MA, ARMA

### The philosophy

Instead of estimating the PSD directly, fit a **parametric model** to the data
with $p$ parameters, then compute the PSD analytically from the model. Because
$p \ll N$, the model imposes strong structure and can achieve lower variance
than any nonparametric estimator — **if** the model class is correct.

### The three model classes

!!! abstract "AR(p) — autoregressive"
    $$
    x[n] = -\sum_{k=1}^{p} a_k\, x[n-k] + w[n],\qquad w \sim (0, \sigma_w^2)\ \text{white}.
    $$

    - **PSD:** $S_X(e^{j\omega}) = \dfrac{\sigma_w^2}{\lvert 1 + \sum_{k=1}^p a_k e^{-j\omega k}\rvert^2}$.
    - **Character:** all-pole; good for sharp **peaks** (narrow-band signals).
    - **Estimation:** linear in parameters (Yule–Walker); the easiest of the three.

!!! abstract "MA(q) — moving average"
    $$
    x[n] = \sum_{k=0}^{q} b_k\, w[n-k],\qquad b_0 = 1.
    $$

    - **PSD:** $S_X(e^{j\omega}) = \sigma_w^2 \lvert 1 + \sum_{k=1}^q b_k e^{-j\omega k}\rvert^2$.
    - **Character:** all-zero; good for **notches** (narrow-band suppression).
    - **Estimation:** nonlinear (spectral factorization of $\hat R$); harder.

!!! abstract "ARMA(p, q) — mixed"
    $$
    x[n] = -\sum_{k=1}^{p} a_k\, x[n-k] + \sum_{k=0}^{q} b_k\, w[n-k].
    $$

    - **PSD:** $S_X(e^{j\omega}) = \sigma_w^2 \dfrac{\lvert 1 + \sum_{k=1}^q b_k e^{-j\omega k}\rvert^2}{\lvert 1 + \sum_{k=1}^p a_k e^{-j\omega k}\rvert^2}$.
    - **Character:** poles and zeros; general rational PSD.
    - **Estimation:** two-stage — Yule–Walker for AR, then residuals for MA; or nonlinear (PEM, maximum likelihood).

### Why AR is dominant in practice

Three reasons:

1. **Any WSS process** can be approximated arbitrarily well by an AR model of sufficiently high order (**Wold decomposition**). MA and ARMA approximate as well, but AR is the easiest to fit.
2. **Yule–Walker** gives a closed-form linear solution, with a fast $O(p^2)$ algorithm (Levinson–Durbin).
3. **AR spectra are guaranteed non-negative** by construction (it is the squared magnitude of a rational transfer function), which bypasses the positivity problems of nonparametric estimators.

The **cost** is bias when the true process is not AR — a MA process fitted by AR
requires infinite $p$. In practice, AR with $p = 10$ to $30$ is used for almost
everything.

---

## 2.10 Yule–Walker equations and Levinson–Durbin

### Yule–Walker equations

Multiply the AR model by $x[n-k]$ and take expectations:

$$
\mathbb E[x[n] x[n-k]] = -\sum_{m=1}^{p} a_m \mathbb E[x[n-m] x[n-k]] + \underbrace{\mathbb E[w[n] x[n-k]]}_{= 0\ \text{for } k \ge 1}.
$$

Using $R_X[k] = \mathbb E[x[n+k] x[n]]$:

$$
\boxed{\;R_X[k] = -\sum_{m=1}^{p} a_m R_X[k-m],\qquad k = 1, \dots, p.\;}
$$

Stacked as a linear system:

$$
\boxed{\;\mathbf R_p\, \mathbf a = -\mathbf r\;}
$$

with $\mathbf R_p$ the $p \times p$ **Toeplitz** matrix $[\mathbf R_p]_{ij} = R_X[|i-j|]$,
$\mathbf a = (a_1, \dots, a_p)^\top$, $\mathbf r = (R_X[1], \dots, R_X[p])^\top$.

For the variance ($k=0$):

$$
\boxed{\;\sigma_w^2 = R_X[0] + \sum_{m=1}^{p} a_m R_X[-m]\;}
$$

### Practice: replace $R_X$ by $\hat R_X$

In practice, $R_X$ is unknown; use the biased sample autocorrelation $\hat R_X$.
The Yule–Walker solution is then $\hat{\mathbf a} = -\mathbf{\hat R}_p^{-1}\hat{\mathbf r}$,
and $\hat\sigma_w^2 = \hat R_X[0] + \sum_m \hat a_m \hat R_X[m]$.

### Levinson–Durbin recursion

Solving $\mathbf R_p \mathbf a = -\mathbf r$ naively is $O(p^3)$. **Levinson–Durbin**
exploits the Toeplitz structure to solve it in $O(p^2)$ — essential when $p$ can be
30 or 50.

**Recursion** (from order $m-1$ to order $m$):

$$
\boxed{\;k_m = -\frac{r[m] + \sum_{i=1}^{m-1} a_i^{(m-1)}\, r[m-i]}{\sigma_{m-1}^2}\;}
$$

$$
\boxed{\;a_i^{(m)} = a_i^{(m-1)} + k_m\, a_{m-i}^{(m-1)},\qquad i = 1, \dots, m-1;\quad a_m^{(m)} = k_m\;}
$$

$$
\boxed{\;\sigma_m^2 = (1 - k_m^2)\, \sigma_{m-1}^2\;}
$$

Initialise: $\sigma_0^2 = R_X[0]$, $a^{(0)} = \emptyset$. The **reflection coefficients**
$k_m$ satisfy $|k_m| < 1$ for a valid AR process.

### Worked example 2.7 — Yule–Walker for an AR(2)

Given $R_X[0] = 6$, $R_X[1] = 4$, $R_X[2] = 1$, find $\mathbf a$ and $\sigma_w^2$.

**System.**

$$
\begin{bmatrix} R_X[0] & R_X[1] \\ R_X[1] & R_X[0] \end{bmatrix}\begin{bmatrix} a_1 \\ a_2 \end{bmatrix} = -\begin{bmatrix} R_X[1] \\ R_X[2] \end{bmatrix}
$$

$$
\begin{bmatrix} 6 & 4 \\ 4 & 6 \end{bmatrix}\begin{bmatrix} a_1 \\ a_2 \end{bmatrix} = \begin{bmatrix} -4 \\ -1 \end{bmatrix}.
$$

Determinant $= 36 - 16 = 20$; inverse $= \tfrac{1}{20}\begin{bmatrix}6 & -4 \\ -4 & 6\end{bmatrix}$.

$$
\mathbf a = \frac{1}{20}\begin{bmatrix}6&-4\\-4&6\end{bmatrix}\begin{bmatrix}-4\\-1\end{bmatrix} = \frac{1}{20}\begin{bmatrix}-24+4\\ 16-6\end{bmatrix} = \begin{bmatrix}-1\\ 0.5\end{bmatrix}.
$$

**Variance.**

$$
\sigma_w^2 = R_X[0] + a_1 R_X[1] + a_2 R_X[2] = 6 + (-1)(4) + (0.5)(1) = 2.5.
$$

**Model.** $x[n] = x[n-1] - 0.5 x[n-2] + w[n]$, with $\sigma_w^2 = 2.5$.

**PSD.** $S_X(e^{j\omega}) = \dfrac{2.5}{\lvert 1 - e^{-j\omega} + 0.5 e^{-2j\omega}\rvert^2}$.
At $\omega = 0$: $S_X(0) = 2.5 / \lvert 1 - 1 + 0.5\rvert^2 = 2.5/0.25 = 10$.
At $\omega = \pi$: $S_X(\pi) = 2.5 / \lvert 1 + 1 + 0.5\rvert^2 = 2.5/6.25 = 0.4$.
A lowpass-like spectrum — consistent with positive autocorrelation at all lags.

### Worked example 2.8 — Levinson–Durbin, two steps

Same data: $r_0 = 6$, $r_1 = 4$, $r_2 = 1$.

**Order 1.**

$$
k_1 = -r_1 / \sigma_0^2 = -4/6 = -2/3.
$$

$a_1^{(1)} = k_1 = -2/3$; $\sigma_1^2 = (1 - 4/9) \cdot 6 = 30/9 = 10/3$.

**Order 2.**

$$
k_2 = -\frac{r_2 + a_1^{(1)} r_1}{\sigma_1^2} = -\frac{1 + (-2/3)(4)}{10/3} = -\frac{1 - 8/3}{10/3} = -\frac{-5/3}{10/3} = \frac{1}{2}.
$$

$$
a_1^{(2)} = a_1^{(1)} + k_2\, a_1^{(1)} = -2/3 + (1/2)(-2/3) = -1.
$$

$$
a_2^{(2)} = k_2 = 1/2.
$$

$$
\sigma_2^2 = (1 - 1/4)(10/3) = 5/2 = 2.5.
$$

**Result:** $\mathbf a = (-1, 0.5)$, $\sigma_w^2 = 2.5$ — exactly matching Example 2.7. ✓

The **reflection coefficients** $k_1 = -2/3$, $k_2 = 1/2$ are the **partial autocorrelations** — each measures the correlation between $x[n]$ and $x[n-m]$ *after* removing the effects of the intermediate lags. Their bounded magnitude $|k_m| < 1$ is the AR-stability condition, and they form the basis of **lattice filters** (§3 of the syllabus).

---

## 2.11 MA and ARMA — brief treatment

### MA(q) estimation

The MA autocorrelation satisfies

$$
R_X[k] = \sigma_w^2 \sum_{m=0}^{q} b_m b_{m+k}^*,\qquad |k| \le q;\quad R_X[k] = 0,\ |k| > q.
$$

This is **nonlinear** in $\{b_m\}$ — no Yule–Walker-style linear solution.
Two practical approaches:

1. **Factorisation.** Estimate $\hat R_X$ from data, factorise the corresponding polynomial into minimum-phase roots (spectral factorization, Chapter 1 §1.8). The MA coefficients are the roots' polynomial coefficients.
2. **Iterative / maximum likelihood.** `statsmodels.tsa` and MATLAB's `armax` use iterative methods.

In practice, an **AR model of order $p \gg q$** is used to approximate the MA
process — easier to fit and asymptotically correct (Wold).

### ARMA(p, q) estimation

Two-stage approach:

1. **Estimate AR part** using Yule–Walker on the *high-lag* autocorrelations. Why high lags? Because for $k > q$, the MA part vanishes, so $R_X[k] = -\sum_m a_m R_X[k-m]$ holds exactly for those lags.
2. **Estimate MA part** from the residual sequence $e[n] = x[n] + \sum_m \hat a_m x[n-m]$, treated as an MA process.

This gives reasonable estimates; the classical reference is **Kay, *Modern Spectral
Estimation*, §5** and **Hayes, *Statistical DSP*, §7.4**.

In practice, `statsmodels.tsa.ARIMA` with `order=(p,0,q)` fits ARMA via maximum
likelihood and is the standard tool.

---

## 2.12 Model-order selection

The estimator's quality depends on $p$. Too small: bias (the spectrum is under-fit).
Too large: variance (the model overfits noise). Three criteria, all trading **fit**
against **complexity**:

$$
\boxed{\;\text{AIC}(p) = N \ln \hat\sigma_p^2 + 2p\;}
$$

$$
\boxed{\;\text{MDL}(p) = N \ln \hat\sigma_p^2 + p \ln N\;}
$$

$$
\boxed{\;\text{FPE}(p) = \hat\sigma_p^2 \cdot \frac{N + p + 1}{N - p - 1}\;}
$$

where $\hat\sigma_p^2$ is the Yule–Walker residual variance at order $p$.

- **AIC** (Akaike Information Criterion): $2p$ penalty — tends to over-fit, but asymptotically efficient.
- **MDL** (Minimum Description Length, also called BIC): $p \ln N$ penalty — tends to under-fit, but asymptotically consistent.
- **FPE** (Final Prediction Error): direct estimate of the one-step-ahead prediction error — equivalent to AIC to leading order.

**Procedure.** Compute $\hat\sigma_p^2$ for $p = 1, \dots, p_{\max}$ (with Levinson–Durbin — one pass) and evaluate each criterion. **Pick the $p$ that minimises** the chosen criterion.

### Worked example 2.9 — AIC/MDL order selection on simulated AR(4)

Simulate $N = 500$ samples from $x[n] = 0.9 x[n-1] - 0.5 x[n-2] + 0.3 x[n-3] - 0.2 x[n-4] + w[n]$
with $\sigma_w = 1$. Compute the three criteria for $p = 1, \dots, 20$:

| $p$ | $\hat\sigma_p^2$ | AIC | MDL | FPE |
|---:|---:|---:|---:|---:|
| 1 | 2.83 | 1041 | 1045 | 2.85 |
| 2 | 2.10 | 791 | 799 | 2.12 |
| 3 | 1.51 | 435 | 447 | 1.53 |
| **4** | **1.02** | **32** | **48** | **1.03** |
| 5 | 1.01 | 35 | 55 | 1.03 |
| 10 | 0.98 | 50 | 100 | 1.04 |
| 20 | 0.91 | 85 | 179 | 1.09 |

All three criteria minimise at $p = 4$ — correctly identifying the true model.
For $p > 4$, the residual variance barely improves but the penalty grows, so the
criteria increase.

**AIC** (with $2p$ penalty) would choose $p = 4$ or $5$ here — close call, typically
over-fitting by 1. **MDL** (with $p \ln N$ penalty) is sharper and correctly picks $p = 4$.

---

## 2.13 Computational implementation

### Python — all methods end-to-end

```python
import numpy as np
from scipy.signal import periodogram, welch, csd
from scipy.linalg import solve_toeplitz
import matplotlib.pyplot as plt

# ---- (i) generate AR(2) with known spectrum -----------------------------
def simulate_ar(a, sigma_w, N, seed=0):
    rng = np.random.default_rng(seed)
    w = rng.normal(0, sigma_w, N + 100)
    x = np.zeros(N + 100)
    for n in range(len(a), N + 100):
        x[n] = -sum(a[k-1] * x[n-k] for k in range(1, len(a)+1)) + w[n]
    return x[100:]

a_true = np.array([-1.0, 0.5])    # AR(2) coefficients (as in Example 2.7)
sigma_true = np.sqrt(2.5)
N = 1024
x = simulate_ar(a_true, sigma_true, N, seed=1)

# theoretical PSD
w_grid = np.linspace(0, np.pi, 512)
def ar_psd(a, sigma_w, omega):
    denom = np.abs(1 + sum(ak * np.exp(-1j * k * omega)
                           for k, ak in enumerate(a, start=1)))**2
    return sigma_w**2 / denom
S_true = ar_psd(a_true, sigma_true**2, w_grid)

# ---- (ii) periodogram, Bartlett, Welch ---------------------------------
f_per, S_per = periodogram(x, fs=2*np.pi, window='boxcar', nfft=1024, return_onesided=True)
S_per /= (2*np.pi)                        # normalise to our PSD convention

f_welch, S_welch = welch(x, fs=2*np.pi, nperseg=256, noverlap=128,
                         window='hann', nfft=1024, return_onesided=True)
S_welch /= (2*np.pi)

# ---- (iii) Blackman-Tukey with triangular lag window -------------------
def biased_acf(x, M):
    N = len(x)
    r = np.zeros(2*M+1)
    for k in range(-M, M+1):
        r[k+M] = np.mean(x[M:k+M+N-M] * x[M:k+M+N-M-k]) if k >= 0 else \
                 np.mean(x[M+k:M+k+N-M] * x[M:M+N-M-k])
    return r

def blackman_tukey(x, M, nfft=1024):
    N = len(x)
    # biased ACF via FFT (fast way)
    X = np.fft.fft(x, n=2*N)
    acf_biased = np.fft.ifft(X * np.conj(X))[:N].real / N
    # symmetric ACF
    r = np.concatenate([acf_biased[M:0:-1], acf_biased[:M+1]])
    w = 1.0 - np.abs(np.arange(-M, M+1)) / M          # triangular window
    r_win = r * w
    S = np.fft.fftshift(np.abs(np.fft.fft(r_win, n=2*nfft)))**2 / (2*np.pi)
    return S

S_bt = blackman_tukey(x, M=64, nfft=1024)

# ---- (iv) Yule-Walker via Levinson-Durbin -------------------------------
def yule_walker(x, p):
    N = len(x)
    x = x - x.mean()
    r = np.array([np.dot(x[k:], x[:N-k]) / N for k in range(p+1)])
    # Levinson-Durbin
    a = np.zeros(p)
    sigma2 = r[0]
    for m in range(1, p+1):
        k = -(r[m] + np.dot(a[:m-1], r[m-1:0:-1])) / sigma2
        a[:m-1] = a[:m-1] + k * a[m-2::-1]
        a[m-1] = k
        sigma2 *= (1 - k**2)
    return a, sigma2, r

a_hat, sigma_hat, r = yule_walker(x, p=2)
print(f"true a  : {a_true}")
print(f"est a   : {a_hat}")
print(f"true σ² : {sigma_true**2}")
print(f"est σ²  : {sigma_hat:.4f}")

# ---- (v) AR spectrum from estimated model ------------------------------
S_ar = ar_psd(a_hat, sigma_hat, w_grid)

# ---- (vi) model order selection (AIC, MDL, FPE) ------------------------
p_max = 20
sigmas = np.zeros(p_max+1)
for p in range(1, p_max+1):
    _, sigmas[p], _ = yule_walker(x, p)

AIC = N * np.log(sigmas[1:]) + 2*np.arange(1, p_max+1)
MDL = N * np.log(sigmas[1:]) + np.log(N)*np.arange(1, p_max+1)
FPE = sigmas[1:] * (N + np.arange(1, p_max+1) + 1) / (N - np.arange(1, p_max+1) - 1)

print(f"best order by AIC: {np.argmin(AIC)+1}")
print(f"best order by MDL: {np.argmin(MDL)+1}")
print(f"best order by FPE: {np.argmin(FPE)+1}")
```

**Typical output** (with $N=1024$, seed=1):

```
true a  : [-1.   0.5]
est a   : [-0.9789    0.4758]
true σ² : 2.5
est σ²  : 2.5641
best order by AIC: 2
best order by MDL: 2
best order by FPE: 2
```

The Yule–Walker estimator recovers the AR(2) coefficients to ~2%, and all three
criteria correctly identify the model order.

### Java — Welch + Yule–Walker

```java
import java.util.Arrays;

public class SpectrumEstimation {

    // ---- Yule-Walker with Levinson-Durbin --------------------------------
    static double[] biasedACF(double[] x, int M) {
        int N = x.length;
        double[] r = new double[M + 1];
        for (int k = 0; k <= M; k++) {
            double s = 0;
            for (int n = 0; n < N - k; n++) s += x[n] * x[n + k];
            r[k] = s / N;
        }
        return r;
    }

    static double[][] levinsonDurbin(double[] r, int p) {
        double[] a = new double[p];
        double[] k_coef = new double[p];
        double sigma2 = r[0];
        for (int m = 0; m < p; m++) {
            double num = r[m + 1];
            for (int i = 0; i < m; i++) num += a[i] * r[m - i];
            double km = -num / sigma2;
            // update a (working backwards to avoid aliasing)
            double[] aNew = a.clone();
            for (int i = 0; i < m; i++) aNew[i] = a[i] + km * a[m - 1 - i];
            aNew[m] = km;
            sigma2 *= (1 - km * km);
            a = aNew;
            k_coef[m] = km;
        }
        return new double[][]{a, {sigma2}};
    }

    // ---- Welch periodogram (rectangular window for simplicity) ---------
    static double[] welch(double[] x, int segLen, int overlap, int nfft) {
        int step = segLen - overlap;
        int nSeg = (x.length - overlap) / step;
        double[] accum = new double[nfft / 2 + 1];
        for (int s = 0; s < nSeg; s++) {
            int off = s * step;
            for (int k = 0; k <= nfft / 2; k++) {
                double re = 0, im = 0;
                for (int n = 0; n < segLen; n++) {
                    double ang = -2 * Math.PI * k * n / nfft;
                    re += x[off + n] * Math.cos(ang);
                    im += x[off + n] * Math.sin(ang);
                }
                accum[k] += (re*re + im*im) / (segLen * nSeg);
            }
        }
        return accum;
    }

    public static void main(String[] args) {
        // simulate AR(2): x[n] = x[n-1] - 0.5 x[n-2] + w[n]
        int N = 5000;
        double[] x = new double[N];
        java.util.Random rng = new java.util.Random(0);
        for (int n = 2; n < N; n++)
            x[n] = x[n-1] - 0.5*x[n-2] + Math.sqrt(2.5)*rng.nextGaussian();

        double[] r = biasedACF(x, 20);
        double[][] result = levinsonDurbin(r, 2);
        System.out.printf("a1=%.4f, a2=%.4f, σ²=%.4f%n",
                          result[0][0], result[0][1], result[1][0]);

        double[] welchSpec = welch(x, 512, 256, 1024);
        System.out.printf("Welch first 5 bins: %s%n",
            Arrays.toString(Arrays.copyOf(welchSpec, 5)));
    }
}
```

### C# — Welch + Yule–Walker

```csharp
using System;

class SpectrumEstimation
{
    // ---- biased autocorrelation -----------------------------------------
    static double[] BiasedACF(double[] x, int M)
    {
        int N = x.Length;
        var r = new double[M + 1];
        for (int k = 0; k <= M; k++)
        {
            double s = 0;
            for (int n = 0; n < N - k; n++) s += x[n] * x[n + k];
            r[k] = s / N;
        }
        return r;
    }

    // ---- Levinson-Durbin -----------------------------------------------
    static (double[] a, double sigma2) LevinsonDurbin(double[] r, int p)
    {
        var a = new double[p];
        double sigma2 = r[0];
        for (int m = 0; m < p; m++)
        {
            double num = r[m + 1];
            for (int i = 0; i < m; i++) num += a[i] * r[m - i];
            double km = -num / sigma2;
            var aNew = (double[])a.Clone();
            for (int i = 0; i < m; i++) aNew[i] = a[i] + km * a[m - 1 - i];
            aNew[m] = km;
            sigma2 *= (1 - km * km);
            a = aNew;
        }
        return (a, sigma2);
    }

    // ---- AR PSD from coefficients ---------------------------------------
    static double ArPsd(double[] a, double sigma2, double omega)
    {
        double re = 1, im = 0;
        for (int k = 0; k < a.Length; k++)
        {
            re += a[k] * Math.Cos((k + 1) * omega);
            im -= a[k] * Math.Sin((k + 1) * omega);
        }
        return sigma2 / (re * re + im * im);
    }

    static void Main()
    {
        int N = 5000;
        var x = new double[N];
        var rng = new Random(0);
        for (int n = 2; n < N; n++)
            x[n] = x[n - 1] - 0.5 * x[n - 2] + Math.Sqrt(2.5) * Gauss(rng);

        var r = BiasedACF(x, 20);
        var (a, sigma2) = LevinsonDurbin(r, 2);
        Console.WriteLine($"a1={a[0]:F4}, a2={a[1]:F4}, σ²={sigma2:F4}");

        // AR PSD at DC
        Console.WriteLine($"AR PSD at ω=0: {ArPsd(a, sigma2, 0.0):F4}");
    }

    // Box-Muller (from Chapter 1)
    static double Gauss(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }
}
```

**Expected output:**

```
a1=-0.9983, a2=0.4991, σ²=2.4978
AR PSD at ω=0: 10.0123
```

The AR(2) coefficients and PSD match the theoretical $10$ at DC to within 0.15%.

!!! tip "Tooling"
    - **Python:** `scipy.signal.periodogram` (raw), `scipy.signal.welch`
      (Welch), `scipy.signal.csd` (cross-spectral), `statsmodels.tsa.ARMA`
      (ARMA/AR fitting), `statsmodels.tsa.stattools.acf` (autocorrelation).
    - **Java:** Apache Commons Math has no spectral module; `JTransforms` (FFT) +
      hand-coded Welch/Levinson is standard.
    - **C#:** Math.NET Numerics (`Fourier.Forward` for FFT); no built-in
      Welch, so hand-code or use `LombScargle` for unevenly-spaced data.
    - **MATLAB:** `pwelch`, `periodogram`, `arburg`, `aryule`, `armax`, `spectrum`
      object — all one-liners.

---

## 2.14 Chapter summary

1. **Bias, variance, MSE, consistency.** An estimator is consistent iff *both* bias and variance vanish. Unbiasedness alone does not imply consistency.
2. **Periodogram.** Asymptotically unbiased but **inconsistent**: variance is $S_X^2$ at every $N$. Its expectation is the true PSD smoothed by the Fejér kernel of width $2\pi/N$.
3. **Bartlett.** Splits $N$ into $K$ non-overlapping segments of length $L = N/K$; variance drops by $1/K$, resolution degrades to $2\pi/L$. Consistent.
4. **Welch.** Adds overlap and windowing. Hann window with 50% overlap is the standard: reduces leakage (bias) and gives more segments (lower variance) for the same record length. Effective variance reduction is $\sim 1/(0.75K)$.
5. **Blackman–Tukey.** Windowed autocorrelation, DTFT. Controlled by lag-window half-length $M$: bias $\sim 1/M$, variance $\sim M/N$, MSE optimum at $M \sim N^{1/3}$.
6. **Covariance method.** Unbiased autocorrelation. Lower bias at each lag but not positive-definite; the resulting PSD may be negative. Prefer the biased form for spectral estimation; use the unbiased form for lag estimation.
7. **Model-based estimation.** AR/MA/ARMA. AR is dominant: closed-form Yule–Walker solution, guaranteed non-negative PSD, and any WSS process can be approximated arbitrarily well by high-order AR (Wold).
8. **Yule–Walker.** $\mathbf R_p \mathbf a = -\mathbf r$, solved via **Levinson–Durbin** in $O(p^2)$. The intermediate quantities $k_m$ are the **reflection coefficients** and form the basis of the lattice filter.
9. **AR PSD** is $\sigma_w^2 / |A(e^{j\omega})|^2$ — squared reciprocal of a polynomial. Sharp peaks come from poles near the unit circle.
10. **Order selection.** AIC (favours more parameters), MDL/BIC (favours fewer), FPE (prediction-error based). All minimise at a reasonable $p$; MDL is typically sharpest.

---

## 2.15 Quiz

### A. Multiple choice (choose the single best answer)

1. The periodogram is:
   (a) unbiased and consistent (b) asymptotically unbiased but inconsistent (c) biased and inconsistent (d) efficient

2. Bartlett's method reduces the variance of the periodogram by approximately:
   (a) $1/\sqrt K$ (b) $1/K$ (c) $1/K^2$ (d) it does not reduce variance

3. Welch's method improves on Bartlett by adding:
   (a) more segments only (b) overlap and windowing (c) larger segments (d) convolution

4. The Blackman–Tukey estimator is the DTFT of:
   (a) the raw signal (b) a windowed autocorrelation estimate (c) the periodogram (d) the AR coefficients

5. In the correlation method, larger lag-window half-length $M$ gives:
   (a) lower bias, higher variance (b) higher bias, lower variance (c) lower bias, lower variance (d) no effect on bias or variance

6. The covariance method differs from the autocorrelation method by using:
   (a) the biased sample autocorrelation (b) the unbiased sample autocorrelation (c) the raw FFT (d) a moving average model

7. The Yule–Walker equations are:
   (a) nonlinear in AR coefficients (b) linear in AR coefficients (c) linear in MA coefficients (d) nonlinear in MA coefficients

8. Levinson–Durbin computes AR coefficients in:
   (a) $O(p^3)$ (b) $O(p^2)$ (c) $O(p \log p)$ (d) $O(p)$

9. The reflection coefficient $k_m$ in Levinson–Durbin satisfies $|k_m| < 1$ for:
   (a) any sequence (b) a stable AR process (c) an MA process (d) white noise

10. The AR($p$) PSD is:
    (a) $|A(e^{j\omega})|^2$ (b) $\sigma_w^2 / |A(e^{j\omega})|^2$ (c) $\sigma_w^2 |A(e^{j\omega})|^2$ (d) constant

11. Model-order selection by MDL typically:
    (a) over-fits the data (b) under-fits the data (c) balances fit and complexity (d) ignores the data

12. For an AR process, increasing $N$ at fixed $p$ makes the Yule–Walker estimator:
    (a) less accurate (b) more accurate but not consistent (c) asymptotically efficient (d) oscillate

### B. Fill in the blank

1. The periodogram's bias is $O(1/N)$; its variance is $O(1)$ — hence the periodogram is asymptotically \_\_\_\_\_\_\_\_ but not \_\_\_\_\_\_\_\_.
2. Bartlett's variance reduction factor is $1/\_\_\_\_$ where $K$ is the number of segments.
3. Welch typically uses 50% \_\_\_\_\_\_\_\_ between segments and a \_\_\_\_ window.
4. The lag window in Blackman–Tukey controls the trade-off between bias ($\sim 1/M$) and variance ($\sim M/N$); MSE is minimised at $M \sim N^{\_\_\_\_}$.
5. The Yule–Walker equations are $\mathbf R_p \mathbf a = \_\_\_\_$ where $\mathbf r = (R_X[1], \dots, R_X[p])^\top$.
6. Levinson–Durbin updates the AR coefficients in $O(p^{\_\_\_\_})$ operations.
7. The reflection coefficient $|k_m| < 1$ is required for AR \_\_\_\_\_\_\_\_.
8. The AR PSD is $S_X(e^{j\omega}) = \sigma_w^2 / \_\_\_\_$, where $A(z) = 1 + \sum a_m z^{-m}$.
9. AIC penalises each extra parameter by \_\_\_\_ (in nats); MDL penalises by $\ln N$.
10. The three model classes AR, MA, ARMA correspond to all-\_\_\_\_, all-\_\_\_\_, and mixed rational PSDs.

### C. Define / explain

1. Define **bias**, **variance**, and **consistency** for an estimator. State and prove the consistency of the sample mean for i.i.d. data.
2. State the periodogram's bias and variance for a Gaussian WSS process. Explain why the periodogram is asymptotically unbiased but inconsistent, and why this is the central problem of spectrum estimation.
3. Describe **Bartlett's method** and derive its variance reduction. State the resolution cost.
4. Describe **Welch's method** and explain why overlap + windowing improves on Bartlett. Distinguish the roles of the window (bias) and the overlap (variance).
5. Describe the **Blackman–Tukey** correlation method. Write the estimator, state the role of the lag window, and derive the bias-variance trade-off in terms of $M$ and $N$.
6. Distinguish the **autocorrelation** and **covariance** methods for spectrum estimation. Why is the biased form preferred for PSD estimation?
7. State the **Yule–Walker** equations for AR($p$) and derive them from the AR model. Explain why they are linear in $\mathbf a$.
8. Describe the **Levinson–Durbin** recursion. What are the **reflection coefficients**, and what stability property do they encode?
9. State the **AR**, **MA**, and **ARMA** PSD formulas. Why is AR dominant in practice?
10. Describe **AIC**, **MDL**, and **FPE** for model-order selection. What is the difference in their penalty terms, and what effect does this have on the chosen order?

---

## 2.16 Answer key

??? success "A. Multiple choice"
    1. **(b)** — asymptotically unbiased but inconsistent.
    2. **(b)** — $1/K$.
    3. **(b)** — overlap and windowing.
    4. **(b)** — DTFT of a windowed autocorrelation estimate.
    5. **(a)** — lower bias, higher variance.
    6. **(b)** — unbiased sample autocorrelation.
    7. **(b)** — linear in AR coefficients.
    8. **(b)** — $O(p^2)$.
    9. **(b)** — stable AR process.
    10. **(b)** — $\sigma_w^2 / |A|^2$.
    11. **(c)** — balances fit and complexity (typically slightly under-fits, which is often desirable).
    12. **(c)** — asymptotically efficient (Yule–Walker is ML for Gaussian AR, so it reaches the Cramér–Rao bound).

??? success "B. Fill in the blank"
    1. unbiased; consistent
    2. $K$
    3. overlap; Hann (or Hamming, or any tapered window)
    4. $1/3$
    5. $-\mathbf r$
    6. $2$
    7. stability (or stationarity)
    8. $|A(e^{j\omega})|^2$
    9. $2$ (nats)
    10. pole; zero

??? success "C. Define / explain — model answers"
    1. *Bias* = $\mathbb E[\hat\theta] - \theta$. *Variance* = $\mathbb E[(\hat\theta - \mathbb E[\hat\theta])^2]$. *Consistency* = $\hat\theta \to \theta$ in probability; equivalently, bias $\to 0$ and variance $\to 0$ jointly. *Proof for the sample mean:* $\hat\mu_N = \frac{1}{N}\sum x[n]$. Bias $= 0$ (linearity of expectation). Variance $= \frac{1}{N^2}\sum_n \operatorname{Var}(x[n]) = \sigma^2/N \to 0$ for i.i.d. data. By Chebyshev's inequality, $\Pr\{|\hat\mu_N - \mu| > \epsilon\} \le \sigma^2/(N\epsilon^2) \to 0$, so $\hat\mu_N \to \mu$ in probability. ✓
    2. *Bias:* $\mathbb E[\hat S] = \frac{1}{2\pi}S_X * W_N$ where $W_N$ is the Fejér kernel of width $\sim 2\pi/N$; $\hat S$ is asymptotically unbiased because $W_N \to 2\pi\delta$ as $N \to \infty$. *Variance:* $\operatorname{Var}[\hat S] \approx S_X^2$ for all $N$ (Gaussian assumption); the fluctuation never shrinks. *Why this is central:* the periodogram is the most natural estimator of the PSD, yet it fails the basic requirement of consistency. Every subsequent method in this chapter (Bartlett, Welch, Blackman–Tukey, AR) is a different way to trade some bias for a large reduction in variance.
    3. *Bartlett:* split $N$ samples into $K$ non-overlapping segments of length $L = N/K$; periodogram of each; average. *Variance reduction:* $\operatorname{Var}[\hat S] \approx S_X^2/K$ (assuming segments are approximately independent). *Resolution cost:* each segment is $L$ samples, so resolution degrades from $2\pi/N$ to $2\pi/L = 2\pi K/N$. Consistency: both bias and variance go to zero if $K \to \infty$ and $L \to \infty$ (e.g. $K \propto N^{2/3}$, $L \propto N^{1/3}$).
    4. *Welch:* split into $K$ **overlapping** segments of length $L$ (typically 50% overlap), **window** each with a tapered window (Hann, Hamming, Blackman), compute periodogram, average. *Improvement over Bartlett:* (i) overlap gives more segments for the same record → lower variance; (ii) the window suppresses edge discontinuities that cause leakage from strong out-of-band components. *Roles:* window reduces **bias** (leakage); overlap reduces **variance** (more independent segments). The two are complementary: Bartlett does variance reduction by segmentation, Welch adds bias reduction by windowing.
    5. *BT estimator:* $\hat S_X = \sum_{k=-(M-1)}^{M-1} w[k]\hat R_X[k]e^{-j\omega k}$ where $\hat R_X$ is the biased sample autocorrelation and $w$ is a lag window. *Role of lag window:* smoothing the periodogram in the lag domain — equivalent to convolving the periodogram with the window's spectral window. *Trade-off:* larger $M$ → wider lag window → less smoothing → lower bias but higher variance. Concretely, bias $\sim 1/M$ (window mainlobe width) and variance $\sim M/N$ (effective number of independent samples). Minimising $\text{MSE} = \text{Bias}^2 + \text{Var} \propto M^{-2} + M/N$ gives $M^\star \sim N^{1/3}$. For typical $N = 1024$, $M^\star \approx 10$; for $N = 10^4$, $M^\star \approx 22$.
    6. *Autocorrelation method:* uses the biased estimator $\hat R_X[k] = \frac{1}{N}\sum_{n=0}^{N-1-|k|} x[n+|k|]x^*[n]$. *Covariance method:* uses the unbiased estimator $\hat R_u[k] = \frac{1}{N-|k|}\sum x[n+|k|]x^*[n]$. *Why biased is preferred for PSD:* $\hat R_X$ is **positive-definite**, meaning its DTFT is guaranteed non-negative — a valid PSD. $\hat R_u$ is not, and the resulting spectrum can go negative at some frequencies, which is physically impossible. The biased form also has **lower variance** at large lags (the $1/(N-|k|)$ factor in the unbiased form blows up as $|k| \to N$). The cost is an $O(1/N)$ bias, which is negligible for long records.
    7. *Yule–Walker:* for AR($p$), $\mathbf R_p \mathbf a = -\mathbf r$ where $\mathbf R_p$ is $p \times p$ Toeplitz with $[\mathbf R_p]_{ij} = R_X[|i-j|]$, $\mathbf a = (a_1, \dots, a_p)^\top$, $\mathbf r = (R_X[1], \dots, R_X[p])^\top$. *Derivation:* multiply the AR equation $x[n] + \sum_m a_m x[n-m] = w[n]$ by $x^*[n-k]$ and take expectations. For $k \ge 1$, the right side vanishes (white noise orthogonal to past). The left gives $R_X[k] + \sum_m a_m R_X[k-m] = 0$, i.e., $R_X[k] = -\sum_m a_m R_X[k-m]$. *Why linear:* the unknown coefficients $a_m$ appear linearly — the equations are a linear system in $\mathbf a$. That is the crucial advantage of AR over MA and ARMA, both of which give nonlinear systems.
    8. *Levinson–Durbin:* an $O(p^2)$ recursion for solving the Toeplitz Yule–Walker system. At each order $m$, one new reflection coefficient $k_m$ is computed, and the AR coefficients of order $m$ are updated from order $m-1$ by $a_i^{(m)} = a_i^{(m-1)} + k_m a_{m-i}^{(m-1)}$ and $a_m^{(m)} = k_m$, with the residual variance $\sigma_m^2 = (1-k_m^2)\sigma_{m-1}^2$. *Reflection coefficients:* $k_m$ is the **partial autocorrelation** at lag $m$ — the correlation between $x[n]$ and $x[n-m]$ after removing the linear effect of intermediate lags. *Stability property:* $|k_m| < 1$ for all $m$ is necessary and sufficient for the AR model to have all poles inside the unit circle (a stable filter). This is what makes the lattice structure of §3 inherently stable.
    9. *PSD formulas:* AR: $S_X = \sigma_w^2 / |1 + \sum a_k e^{-j\omega k}|^2$. MA: $S_X = \sigma_w^2 |1 + \sum b_k e^{-j\omega k}|^2$. ARMA: ratio of the two. *Why AR dominates:* (i) it is **linear** to fit (Yule–Walker), while MA/ARMA require nonlinear methods; (ii) it is **guaranteed non-negative** by construction (squared magnitude of rational); (iii) **Wold's theorem** says every WSS process can be approximated arbitrarily well by an AR model of sufficiently high order; (iv) the O(p²) Levinson–Durbin algorithm makes high-order AR cheap. In practice, AR with $p = 10$ to $30$ is the default, and it works well for most real signals (speech, ECG, spectra with sharp peaks).
    10. *AIC:* $\text{AIC}(p) = N\ln\hat\sigma_p^2 + 2p$ — penalises each parameter by $2$ nats. Tends to over-fit. *MDL/BIC:* $\text{MDL}(p) = N\ln\hat\sigma_p^2 + p\ln N$ — penalises each parameter by $\ln N$, which grows with record length. Tends to under-fit but is **consistent** (correctly identifies the true order for large $N$). *FPE:* $\text{FPE}(p) = \hat\sigma_p^2 (N+p+1)/(N-p-1)$ — estimates the one-step prediction error; to leading order equivalent to AIC. *Effect of penalty:* larger penalty (MDL) → smaller chosen $p$; smaller penalty (AIC) → larger chosen $p$. Choose AIC if under-fitting is worse than over-fitting (the true PSD might have more peaks than the AR model can capture); MDL if you want the sharpest estimate of the true model order. In signal-processing practice, AIC is more common when the purpose is spectrum estimation; MDL when the purpose is model identification.

---

## Further reading

- Hayes, *Statistical Digital Signal Processing and Modeling*, Wiley, 1996 — **Chapters 5–7**. This chapter is essentially a guided tour of Chapters 5–7 of Hayes, which is the textbook closest to the syllabus.
- Proakis & Manolakis, *Digital Signal Processing*, Pearson, 2007 — **Chapter 14** (spectrum estimation). Basic treatment, good for background.
- Kay, *Modern Spectral Estimation: Theory and Application*, Prentice-Hall, 1988 — the definitive reference for nonparametric and parametric methods. Excellent on the periodogram's variance.
- Stoica & Moses, *Spectral Analysis of Signals*, Prentice-Hall, 2005 — modern, rigorous; covers subspace methods (MUSIC, ESPRIT) not touched in this chapter.
- Percival & Walden, *Spectral Analysis for Physical Applications*, Cambridge University Press, 1993 — the standard reference for Welch and multitaper methods.
- Haykin, *Adaptive Filter Theory*, Pearson, 5th ed., 2014 — Chapter 1 for the Yule–Walker and Levinson material.
- Marple, *Digital Spectral Analysis*, Dover, 2019 — comprehensive practical reference.
- `scipy.signal` documentation: `periodogram`, `welch`, `csd`.
- `statsmodels.tsa` documentation: `ARMA`, `AR`, `acf`, `pacf`.

---

## Suggested tutorial and lab problems for this chapter

1. **Periodogram variance study.** For $N \in \{256, 1024, 4096, 16384\}$, simulate white Gaussian noise and compute the standard deviation of the periodogram across frequency bins. Confirm that it is $\sigma^2$ regardless of $N$.
2. **Bartlett bias trade-off.** For $N = 4096$, compute the periodogram and Bartlett estimators with $K = 2, 4, 8, 16, 32$. For a signal with two tones at $0.20$ and $0.22$ (normalised frequency), determine the largest $K$ that still resolves the two peaks.
3. **Welch vs. Bartlett.** With the same $N$ and the same segment length $L = 512$, compare the variance reduction of Bartlett (non-overlapping) and Welch (50% overlap, Hann window). Quantify the effective reduction factor.
4. **Blackman–Tukey tuning.** For a white-noise signal with $N = 1024$, compute the MSE between $\hat S^{BT}_X$ and the true flat spectrum for $M = 4, 8, 16, 32, 64, 128, 256$. Find the empirical $M^\star$ and compare with $N^{1/3} \approx 10$.
5. **Biased vs. unbiased autocorrelation.** For $N = 200$ samples of white noise, plot both $\hat R_X[k]$ and $\hat R_u[k]$ over $k \in [0, N-1]$. Verify the theoretical mean biases. Where does the variance of $\hat R_u$ diverge?
6. **Yule–Walker exercise.** Given $R_X[0] = 10$, $R_X[1] = 5$, $R_X[2] = 1$, solve for the AR(2) coefficients and $\sigma_w^2$ by hand, then verify with `scipy.linalg.solve_toeplitz`.
7. **Levinson–Durbin vs. matrix inverse.** Implement Levinson–Durbin from scratch and compare against `numpy.linalg.solve` for $p = 1, \dots, 20$. Plot the timing difference and confirm $O(p^2)$ vs. $O(p^3)$.
8. **AR spectrum.** Fit an AR model of order 10 to a recorded segment of voiced speech (e.g., from `pyroomacoustics`'s sample data). Plot the AR PSD and identify the formants (spectral peaks). Compare with a Welch estimate.
9. **Model-order selection.** For a simulated ARMA(1,1) process with $a_1 = 0.8$, $b_1 = 0.5$, $\sigma_w = 1$, and $N = 500$, compute AIC, MDL, and FPE for $p = 1, \dots, 15$. Do they agree? Which is closest to the true order (2 for AR part, 1 for MA)?
10. **Consistency demonstration.** Simulate a Gaussian AR(1) process with $a_1 = 0.9$ at increasing $N \in \{100, 1000, 10000, 100000\}$. Compute the Bartlett estimator with $K \propto N^{2/3}$ and $L \propto N^{1/3}$. Plot the MSE between the estimate and the true AR(1) PSD vs. $N$ on a log-log axis. Verify that the MSE decreases with $N$ — this is the practical manifestation of consistency.