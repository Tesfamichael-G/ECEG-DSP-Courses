# Exam Bank — ECEG-6303 Advanced Digital Signal Processing

40 questions spanning all five chapters. Use for self-assessment and exam
preparation. **Answers and grading rubric at the end.**

**Instructions (for a mock exam):**

- Time: **3 hours**
- Closed book; **one A4 handwritten formula sheet** allowed
- Non-programmable calculator permitted
- Total marks: **100**
- Marks per question in brackets

---

## Section A — Discrete Random Signal Processing (Ch 1) — 20 marks

### A1. Multiple choice [2]

For a finite-length sequence, **Parseval's theorem** with the $N$-point DFT states:

(a) $\sum |x[n]|^2 = \sum |X[k]|^2$
(b) $\sum |x[n]|^2 = \frac{1}{N} \sum |X[k]|^2$
(c) $\sum |x[n]|^2 = N \sum |X[k]|^2$
(d) $\sum |x[n]|^2 = \frac{1}{N^2} \sum |X[k]|^2$

### A2. Multiple choice [2]

The **Wiener–Khinchine relation** for a WSS discrete process states that:

(a) The PSD is the DFT of the signal
(b) The PSD is the DTFT of the autocorrelation
(c) The autocorrelation is the DFT of the signal
(d) The PSD is the mean of the signal

### A3. Short answer [3]

For a finite-length sample $x[0], \ldots, x[N-1]$, state the **biased** and
**unbiased** sample autocorrelation estimators. Which is preferred for
spectral estimation and why?

### A4. Short answer [3]

State the **Sum Decomposition Theorem** (discrete Karhunen–Loève). What is the
truncation error when keeping only $M$ of $N$ components?

### A5. Derivation [5]

State the **Spectral Factorization Theorem** for a rational PSD $S_X(z)$.
Show that if $S_X = \sigma_e^2 H(z)H(z^{-1})$ with $H$ causal and minimum-phase,
then $H^{-1}(z)$ is a causal whitening filter.

### A6. Short answer [5]

For a WSS input $X$ to an LTI system with impulse response $h$:
(i) state the output PSD, output autocorrelation, and output mean;
(ii) explain why the mean is filtered by $H(e^{j0})$ but the PSD is multiplied by $|H(e^{j\omega})|^2$.

---

## Section B — Spectrum Estimation (Ch 2) — 20 marks

### B1. Multiple choice [2]

The **periodogram** of a Gaussian WSS process is:

(a) Consistent (b) Inconsistent — variance does not shrink with $N$
(c) Biased and inconsistent (d) Exact

### B2. Multiple choice [2]

**Welch's method** improves on Bartlett's by adding:

(a) More segments only (b) Overlap and windowing
(c) Larger segments (d) Higher FFT resolution

### B3. Short answer [3]

State the **bias-variance trade-off** of the periodogram. Why is it a problem
in practice, and how do Bartlett and Welch address it?

### B4. Short answer [3]

Describe the **Blackman–Tukey** estimator. State the role of the lag window
and the MSE-optimal choice of its half-length $M$.

### B5. Derivation [5]

For an AR($p$) process with autocorrelation $R[k]$, derive the **Yule–Walker
equations** by minimizing the prediction error $\mathcal E = \sum e[n]^2$.
State the resulting linear system and the prediction-error variance formula.

### B6. Short answer [5]

State the three **model-order selection criteria** AIC, MDL, and FPE. For
$N = 500$ training samples and a simulated AR(4) process, which criterion
would you trust most, and why?

---

## Section C — Linear Estimation and Prediction (Ch 3) — 20 marks

### C1. Multiple choice [2]

The **Cramér–Rao lower bound** gives:

(a) A lower bound on the bias of an estimator
(b) A lower bound on the variance of any unbiased estimator
(c) The minimum MSE for nonlinear estimators
(d) The Fisher information

### C2. Multiple choice [2]

For **jointly Gaussian** signals, the LMMSE estimator:

(a) Is worse than MMSE (b) Is identical to the MMSE
(c) Is nonlinear (d) Does not exist

### C3. Short answer [3]

State the **orthogonality principle** for LMMSE. Derive the closed-form
solution $\hat S = \mu_S + \Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$.

### C4. Derivation [5]

State the **Kalman filter** predict and update equations. Show that the
steady-state behaviour (via the discrete algebraic Riccati equation) matches
the causal Wiener filter for a stationary process.

### C5. Short answer [5]

State the **linear prediction** problem for a WSS sequence. Show that the
optimal predictor coefficients satisfy the **Yule–Walker** equations, and
that the prediction error is the **innovations** process.

### C6. Short answer [3]

Describe the **lattice filter** structure. State the stability condition
$|k_m| < 1$ and explain why it is a local (per-stage) condition.

---

## Section D — Adaptive Filters (Ch 4) — 20 marks

### D1. Multiple choice [2]

The **LMS** algorithm minimizes at each step:

(a) The global MSE (b) The instantaneous squared error $e^2(n)$
(c) The weighted least-squares cost (d) The maximum likelihood

### D2. Multiple choice [2]

**RLS** convergence rate is:

(a) Strongly dependent on eigenvalue spread
(b) Independent of eigenvalue spread
(c) Always slower than LMS
(d) Unbounded

### D3. Short answer [3]

State the **LMS convergence condition** (mean and mean-square), the **time
constant** of the slowest mode, and the **misadjustment** formula.

### D4. Derivation [5]

Derive the **LMS algorithm** from steepest descent on the MSE surface.
Replace the true gradient with the instantaneous estimate and derive the
update rule $\mathbf w(n+1) = \mathbf w(n) + \mu\, e(n)\,\mathbf x(n)$.

### D5. Short answer [4]

State the **NLMS** update. Explain why NLMS converges at a rate independent
of input power, and give one application where this matters.

### D6. Short answer [4]

Describe **adaptive noise cancellation**. Why must the reference noise be
correlated with the noise in the primary input? State the achievable SNR
improvement in terms of correlation coefficient $\rho$.

---

## Section E — Multirate DSP (Ch 5) — 20 marks

### E1. Multiple choice [2]

For **decimation by $M$**, the required anti-aliasing filter cutoff is:

(a) $\pi$ (b) $\pi/M$ (c) $\pi/(2M)$ (d) $\pi L$

### E2. Multiple choice [2]

The **Noble identity** for decimation states:

(a) $H(z) \downarrow M = \downarrow M \, H(z^M)$
(b) $H(z) \downarrow M = \downarrow M \, H(z)$
(c) $H(z) \downarrow M = M H(z/M)$
(d) $H(z) \downarrow M = \downarrow M \, H(-z)$

### E3. Short answer [3]

Describe the **polyphase decomposition** of a length-$N_h$ FIR filter for
$M$ channels. State the computational advantage of the polyphase decimator
over the direct form.

### E4. Derivation [5]

Derive the **polyphase decomposition** $H(z) = \sum_{r=0}^{M-1} z^{-r} E_r(z^M)$
and the corresponding polyphase decimator structure. Show the number of
multiplications per output sample.

### E5. Short answer [4]

For a two-channel **QMF** filter bank with $H_1(z) = H_0(-z)$, $F_0(z) = H_0(z)$,
$F_1(z) = -H_1(z)$, state the **perfect reconstruction** conditions and
verify that the aliasing-cancellation condition holds.

### E6. Short answer [4]

Describe the **discrete wavelet transform** (DWT) as an iterated filter bank.
Why does the DWT provide multi-resolution analysis, and what is the role of
**vanishing moments**?

---

## Answer key

??? success "Click to reveal answers"

    ### Section A

    **A1.** (b) — $\frac{1}{N}$ factor with the DFT convention where inverse has $1/N$.

    **A2.** (b).

    **A3.** **Biased:** $\hat R[k] = \frac{1}{N}\sum_{n=0}^{N-1-|k|} x[n+|k|]x^*[n]$. **Unbiased:** $\hat R_u[k] = \frac{1}{N-|k|}\sum x[n+|k|]x^*[n]$. **Biased preferred for spectral estimation** because (i) it is positive-definite (its DTFT is a valid non-negative PSD); (ii) it has lower variance at large lags; (iii) its $O(1/N)$ bias vanishes asymptotically.

    **A4.** For zero-mean random vector $\mathbf X$ of length $N$ with covariance $\mathbf R$, decompose $\mathbf R = \sum_{k=1}^N \lambda_k \mathbf v_k \mathbf v_k^\top$ and $\mathbf X = \sum_{k=1}^N Z_k \mathbf v_k$ with $\mathbb E[Z_i Z_j] = \lambda_i \delta_{ij}$. **Truncation error:** $E[\|\mathbf X - \hat{\mathbf X}_M\|^2] = \sum_{k=M+1}^N \lambda_k$ — the residual energy.

    **A5.** **Spectral Factorization Theorem:** any rational PSD $S_X(z)$ with no zeros on the unit circle and $\int\log S_X > -\infty$ can be factored uniquely as $S_X(z) = \sigma_e^2 L(z) L(z^{-1})$ with $L$ causal, monic, and minimum-phase. **Whitening:** $H^{-1}(z) = 1/L(z)$ is causal because $L$ has all poles and zeros inside the unit circle (minimum-phase = all zeros inside; causality = all poles inside). Its output is white of variance $\sigma_e^2$.

    **A6.** (i) $\mu_Y = \mu_X H(e^{j0})$, $S_Y = |H|^2 S_X$, $R_Y = R_X * (h*h^\dagger)$ where $h^\dagger[k] = h^*[-k]$. (ii) **Mean:** the DC component of $X$ is what contributes to $\mu_Y$; only the filter's DC gain $H(e^{j0})$ matters because a constant passes a filter with only its DC gain. **PSD:** each frequency $\omega$ of the input is scaled by the filter's magnitude response at that frequency — hence the pointwise multiplication by $|H(e^{j\omega})|^2$.

    ### Section B

    **B1.** (b) — inconsistent.

    **B2.** (b) — overlap and windowing.

    **B3.** Asymptotically unbiased but **inconsistent** — its variance is $S_X^2$ at every $N$. Practically, the periodogram is a "spiky noise" that never smooths out with more data. **Bartlett** averages over non-overlapping segments to reduce variance (÷$K$); **Welch** adds overlap (more segments) and windowing (lower bias). Both achieve consistency at the cost of resolution.

    **B4.** **BT:** $\hat S_X(\omega) = \sum_{k=-(M-1)}^{M-1} w[k] \hat R_X[k] e^{-j\omega k}$. Lag window $w[k]$ smooths the periodogram — larger $M$ = less smoothing = lower bias but higher variance. **Bias-variance trade-off:** Bias $\sim 1/M$, variance $\sim M/N$; MSE-optimal at $M \sim N^{1/3}$.

    **B5.** Model: $x[n] = -\sum_{k=1}^p a_k x[n-k] + e[n]$. Minimize $\mathcal E = \sum_n e[n]^2 = \sum_n (x[n] + \sum a_k x[n-k])^2$. $\partial \mathcal E/\partial a_k = 0$ gives $\sum_n (x[n] + \sum a_l x[n-l])x[n-k] = 0$, i.e., $R[k] + \sum_l a_l R[k-l] = 0$, so $\sum_l a_l R[k-l] = -R[k]$ for $k=1,\ldots,p$. Matrix form: $\mathbf R \mathbf a = -\mathbf r$ with $\mathbf R$ Toeplitz. **Variance:** $\sigma_e^2 = R[0] + \sum_k a_k R[k]$.

    **B6.** **AIC:** $N\log\hat\sigma_p^2 + 2p$ — tends to over-fit. **MDL:** $N\log\hat\sigma_p^2 + p\log N$ — consistent (correct for large $N$), tends to under-fit. **FPE:** $\hat\sigma_p^2 (N+p+1)/(N-p-1)$ — direct estimate of one-step prediction error. **For $N=500$, AR(4):** all three typically agree near $p=4$. MDL is the most reliable for identifying the true order because its penalty $\log N$ grows with data, giving asymptotic consistency. AIC's $2p$ penalty is fixed, so it tends to slightly over-fit — but at moderate $N$, both work.

    ### Section C

    **C1.** (b).

    **C2.** (b).

    **C3.** **Orthogonality principle:** the LMMSE estimator $\hat S$ of $S$ from $Y$ satisfies $\mathbb E[(S - \hat S)Y^\top] = 0$ and $\mathbb E[S - \hat S] = 0$. Write $\hat S = a + b^\top Y$; $\mathbb E[S - \hat S] = 0$ gives $a = \mu_S - b^\top \mu_Y$; $\mathbb E[(S - \hat S)Y^\top] = 0$ gives $\Sigma_{SY} - b^\top \Sigma_{YY} = 0$, so $b = \Sigma_{YY}^{-1}\Sigma_{SY}^\top = \Sigma_{YY}^{-1}\Sigma_{YS}$. Substituting: $\hat S = \mu_S + \Sigma_{SY}\Sigma_{YY}^{-1}(Y - \mu_Y)$.

    **C4.** **Predict:** $\hat x_{k|k-1} = A\hat x_{k-1|k-1}$, $P_{k|k-1} = AP_{k-1|k-1}A^\top + Q$. **Update:** $K_k = P_{k|k-1}H^\top(HP_{k|k-1}H^\top + R)^{-1}$; $\hat x_{k|k} = \hat x_{k|k-1} + K_k(z_k - H\hat x_{k|k-1})$; $P_{k|k} = (I - K_kH)P_{k|k-1}$. **Steady state:** the Riccati equation converges to $P^\star$, giving fixed $K^\star$. The steady-state filter is LTI with transfer function $H_{KF}(z) = K^\star/(1 - (A - K^\star H A)z^{-1})$ — matching the causal Wiener filter for the same signal and noise models.

    **C5.** Predictor: $\hat x[n] = -\sum_{k=1}^p a_k x[n-k]$. Error: $e[n] = x[n] + \sum a_k x[n-k]$. Minimize $\mathbb E[e^2]$ over $a_k$: $\partial/\partial a_k = 2 \mathbb E[e[n] x[n-k]] = 0$ → orthogonality. Expanding gives Yule–Walker: $R[k] + \sum_m a_m R[k-m] = 0$. **Prediction error = innovations** because the optimal predictor from the past leaves only the unpredictable part — exactly the innovations process of the Wold decomposition.

    **C6.** **Lattice:** cascade of $p$ stages, each with forward error $f_m[n] = f_{m-1}[n] + k_m g_{m-1}[n-1]$ and backward error $g_m[n] = g_{m-1}[n-1] + k_m f_{m-1}[n]$. **Stability:** $|k_m| < 1$ for all $m$ — each stage is an all-pass-like reflector with less than unit gain. **Why local:** each stage's stability depends only on that stage's $k_m$; adding a stage cannot destabilise the earlier ones. The direct-form $a_k$ coefficients, by contrast, require solving a root-finding problem to check stability globally.

    ### Section D

    **D1.** (b).

    **D2.** (b).

    **D3.** **Convergence:** $0 < \mu < 2/\lambda_{\max}$ (mean and mean-square). **Time constant:** $\tau_k = 1/(\mu\lambda_k)$; slowest mode $\tau_{\max} = 1/(\mu\lambda_{\min})$. **Misadjustment:** $M \approx \mu \operatorname{tr}(\mathbf R)/2$.

    **D4.** Steepest descent: $\mathbf w(n+1) = \mathbf w(n) - \mu \nabla J$ with $\nabla J = -2\mathbf p + 2\mathbf R\mathbf w$. Replace true gradient with instantaneous estimate $\widehat{\nabla J} = -2e(n)\mathbf x(n)$ (the gradient of $e^2(n)$). Then $\mathbf w(n+1) = \mathbf w(n) - \mu \cdot \frac{1}{2} \cdot (-2e(n)\mathbf x(n)) = \mathbf w(n) + \mu e(n)\mathbf x(n)$.

    **D5.** $\mathbf w(n+1) = \mathbf w(n) + \frac{\tilde\mu}{\|\mathbf x(n)\|^2 + \delta} e(n) \mathbf x(n)$. **Rate independence:** the effective step size is $\tilde\mu$, independent of input power $\|\mathbf x\|^2$. Plain LMS's effective step is $\mu\|\mathbf x\|^2$, which varies with input amplitude. **Application:** speech echo cancellation — far-end speech power varies by 40+ dB; NLMS converges stably throughout, plain LMS diverges at loud peaks or is too slow at quiet troughs.

    **D6.** Primary $d = s + v_0$ (signal + noise); reference $x = v_1$ (correlated noise). Adaptive filter $w$ predicts $v_0$ from $x$; output $e = d - w*x \approx s$. **Correlation requirement:** if $x$ and $v_0$ are uncorrelated, no linear filter of $x$ predicts $v_0$; adaptation has no useful gradient. **SNR improvement:** $1/(1-\rho^2)$ where $\rho$ is the correlation between $v_0$ and $v_1$. At $\rho = 0.9$: 7 dB; at $\rho = 0.99$: 17 dB; at $\rho = 0$: no improvement.

    ### Section E

    **E1.** (b) — $\pi/M$.

    **E2.** (a) — $H(z)\downarrow M = \downarrow M \, H(z^M)$.

    **E3.** Split $h[n]$ into $M$ polyphase subsequences $E_r$ (every $M$-th tap starting at offset $r$), so $H(z) = \sum_{r=0}^{M-1} z^{-r} E_r(z^M)$. **Advantage:** the polyphase decimator runs each $E_r$ at the **lower rate** $f_s/M$ — no wasted computation, no high-rate memory access. Multiplications per output same as direct form, but throughput 2× (or $M$×) faster because all computation is at the lower rate.

    **E4.** $H(z) = \sum_{k} h[k] z^{-k} = \sum_{r=0}^{M-1} z^{-r} \sum_m h[mM+r] z^{-mM} = \sum_r z^{-r} E_r(z^M)$ with $E_r(z) = \sum_m h[mM+r] z^{-m}$. **Polyphase decimator:** $y[n] = \sum_r (E_r * x_r)[n]$ where $x_r[n] = x[nM - r]$. Multiplications per output: $M \cdot (N_h/M) = N_h$, all at rate $f_s/M$; per input sample equivalent $N_h/M$.

    **E5.** **PR conditions:** aliasing cancellation $F_0 H_0(-z) + F_1 H_1(-z) = 0$; distortion elimination $F_0 H_0 + F_1 H_1 = 2 z^{-d}$. **Check with QMF choice:** $F_0 H_0(-z) + F_1 H_1(-z) = H_0(z)H_0(-z) + (-H_0(-z))H_0(z) = 0$ ✓. The distortion condition becomes $H_0^2(z) - H_0^2(-z) = 2 z^{-d}$, the half-band condition (power-complementary in frequency).

    **E6.** **DWT as iterated filter bank:** apply a 2-channel QMF to the signal, then apply the same QMF to the lowpass output, recursively $J$ times. Produces $J$ detail subbands (highpass at each level) plus 1 approximation subband. **Multi-resolution:** short wavelets (high frequency) resolve transient structure in time; long wavelets (low frequency) resolve fine spectral structure. **Vanishing moments:** a wavelet with $p$ vanishing moments has zero inner product with any polynomial of degree $<p$, so the detail coefficients of smooth signals are near zero — sparse representation → compression and denoising.

---

## Grading rubric (for instructors)

| Band | Marks | Criteria |
|---|---|---|
| **A** | 90–100 | All correct or with trivial slips; rigorous derivations; clear interpretation |
| **B** | 75–89 | Most correct; occasional gaps; derivations mostly sound |
| **C** | 60–74 | Basic understanding; some errors in derivations or interpretation |
| **F** | < 60 | Major conceptual gaps; significant derivations missing |

**Common error checks (and deductions):**

- **Missing $1/N$ in Parseval** (DFT convention): −1.
- **Claiming periodogram is consistent:** −2.
- **Yule–Walker sign convention** ($\mathbf R\mathbf a = +\mathbf r$ instead of $-\mathbf r$): −1 to −2.
- **Kalman gain missing matrix inverse:** −2.
- **LMS convergence bound** $0 < \mu < 2/\lambda_{\min}$ instead of $2/\lambda_{\max}$: −2.
- **Noble identity direction** reversed: −2.
- **QMF aliasing condition** missing one filter sign: −1.
- **$|k_m| < 1$ interpreted as global stability without justification:** −1.
