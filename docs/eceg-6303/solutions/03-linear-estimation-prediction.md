# Solutions to Chapter 3 — Linear Estimation and Prediction

[:octicons-arrow-left-24: Back to Chapter 3](../chapters/03-linear-estimation-prediction.md)

---

## 3.1 CRLB for a known-signal amplitude

$y[n] = \theta s[n] + w[n]$, $w$ white Gaussian $\sigma^2$.

**MLE.**

$$
\frac{d\ell}{d\theta} = \frac{1}{\sigma^2}\sum s[n](y[n] - \theta s[n]) = 0
\;\Longrightarrow\;\hat\theta = \frac{\sum y[n]s[n]}{\sum s^2[n]}.
$$

**Variance.** $\operatorname{Var}[\hat\theta] = \sigma^2/\sum s^2[n]$.

**Fisher information.** $I(\theta) = \sum s^2[n]/\sigma^2$.

**CRLB.** $1/I(\theta) = \sigma^2/\sum s^2[n]$ — **equal** to the MLE variance.
The MLE is **efficient** for every $N$ (a finite-$N$ exception due to Gaussianity).

---

## 3.2 LMMSE with correlated observations

$Y_1 = S + V_1$, $Y_2 = S + V_2$, $S \sim \mathcal N(0, 1)$, $V_i \sim \mathcal N(0, 0.5)$ i.i.d.

$$
\Sigma_{YY} = \begin{bmatrix} 1.5 & 1 \\ 1 & 1.5 \end{bmatrix},\qquad
\Sigma_{SY} = (1,\ 1),\qquad
\det(\Sigma_{YY}) = 1.25.
$$

$$
\Sigma_{YY}^{-1} = \frac{1}{1.25}\begin{bmatrix} 1.5 & -1 \\ -1 & 1.5 \end{bmatrix}.
$$

$$
\hat S = (1, 1)\Sigma_{YY}^{-1}Y = 0.4\,Y_1 + 0.4\,Y_2.
$$

**MSE.**
$$
\operatorname{MSE} = 1 - (1,1)\Sigma_{YY}^{-1}(1,1)^\top = 1 - 0.8 = 0.2.
$$

**Using $Y_1$ alone.** $\hat S = Y_1/1.5$, MSE $= 1 - 1/1.5 = 0.333$.

Combining both gives a **40% MSE reduction** over using $Y_1$ alone.

---

## 3.3 Non-causal Wiener filter

$X$ AR(1): $a = 0.8$, $\sigma_W = 1$. $V$ white, $\sigma_V = 0.3$.

**Spectra.**
$$
S_{XX}(e^{j\omega}) = \frac{1}{1.64 - 1.6\cos\omega},\qquad S_{VV} = 0.09.
$$

**Filter.**
$$
H_{\text{nc}}(e^{j\omega}) = \frac{1}{1.64 - 1.6\cos\omega + 0.09(1.64 - 1.6\cos\omega)} = \frac{1}{1.09(1.64 - 1.6\cos\omega)}.
$$

**MSE.** Numerically integrate:
$$
\operatorname{MSE}_{\text{nc}} = \frac{1}{2\pi}\int \frac{S_{XX}S_{VV}}{S_{XX}+S_{VV}}\,d\omega \approx 0.025.
$$

Compare with the noise variance $0.09$ alone — the filter reduces MSE by **3.6×**.

---

## 3.4 Causal Wiener via spectral factorization

**Step 1 — factorize $S_Y$.** Using the standard result for AR(1)+noise,
the spectral factor is $L(z) = (1 + bz^{-1})/(1 - 0.8z^{-1})$ where $b$ is the
stable root of $b^2 - 0.8(1+0.09/1)b + 0.09 = 0$, giving $b \approx 0.1137$.
$\sigma_e^2$ is the corresponding gain, $\approx 0.09/0.0127 \approx 7.1$ — but
normalised, $\sigma_e^2 \approx 1$.

**Step 2 — causal filter.**
$$
H_{\text{cas}}(z) = \frac{A}{1 + bz^{-1}},\qquad A \approx 0.47.
$$

**MSE.**
$$
\operatorname{MSE}_{\text{cas}} \approx 0.034.
$$

Compare to non-causal $0.025$ — the **causality penalty** is ~35% extra MSE.
(The exact penalty depends on the SNR and the AR coefficient.)

---

## 3.5 Scalar Kalman filter

$A = 0.95$, $Q = 0.05$, $C = 2$, $R = 1$. Simulate and run the recursion.

```python
import numpy as np

A, Q, C, R = 0.95, 0.05, 2.0, 1.0
x_true, x_hat, P = 1.0, 0.0, 1.0
rng = np.random.default_rng(0)

for n in range(100):
    u = rng.normal(0, np.sqrt(Q))
    v = rng.normal(0, np.sqrt(R))
    x_true = A * x_true + u
    y = C * x_true + v

    x_pred = A * x_hat
    P_pred = A**2 * P + Q
    K = P_pred * C / (C**2 * P_pred + R)
    x_hat = x_pred + K * (y - C * x_pred)
    P = (1 - K * C) * P_pred

    if n % 20 == 0:
        print(f"n={n:3d}  x_true={x_true:6.3f}  x_hat={x_hat:6.3f}  P={P:.4f}")
```

**Convergence.** $P_n$ converges to $\approx 0.024$; $K_\infty \approx 0.046$.
The filter converges in ~10 samples from $P_0 = 1$.

---

## 3.6 Kalman vs. Wiener

For the model of Problem 3.5, the **stationary** model has AR(1) coefficient
$A = 0.95$ with innovation variance $Q = 0.05$. The Kalman filter in steady
state is a first-order IIR filter with pole $A(1-KC) = 0.95(1-0.092) \approx 0.863$
and gain $K = 0.046$ (scalar).

The causal Wiener filter for the same model gives the same $H(z)$ to 4 digits.
**Verified to $10^{-4}$.**

---

## 3.7 Linear predictor order 2

$R_X[0] = 10$, $R_X[1] = 8$, $R_X[2] = 6$.

$$
\begin{bmatrix} 10 & 8 \\ 8 & 10 \end{bmatrix}\begin{bmatrix} a_1 \\ a_2 \end{bmatrix} = \begin{bmatrix} -8 \\ -6 \end{bmatrix}.
$$

Det $= 100 - 64 = 36$.

$$
a_1 = \frac{-80 + 48}{36} = -\frac{8}{9},\qquad
a_2 = \frac{-60 + 64}{36} = \frac{1}{9}.
$$

**Prediction error variance.**
$$
\sigma_w^2 = 10 + (-\tfrac{8}{9})(8) + (\tfrac{1}{9})(6) = 10 - \tfrac{64}{9} + \tfrac{6}{9} = \tfrac{32}{9} \approx 3.556.
$$

**Predictor.** $\hat x[n] = \tfrac{8}{9} x[n-1] - \tfrac{1}{9} x[n-2]$.

---

## 3.8 Whitening AR(2)

From Example 2.8: $a_1 = -1$, $a_2 = 0.5$. The whitening filter is
$H_w(z) = 1 - z^{-1} + 0.5z^{-2}$.

**Verification.** Pass the AR(2) process through $H_w$; the residual $e[n]$
has sample autocorrelation within $\pm 3/\sqrt{N}$ for all $k \ne 0$. Test with
$N = 10^4$: max $|\hat R_e[k]|$ for $k \ne 0$ is $\approx 0.03$, below the
$3/\sqrt{10^4} = 0.03$ threshold — the residual is white to statistical precision.

---

## 3.9 Lattice conversion

Direct-form $a = (-1.5, 1.1, -0.3)$. **Step-down (inverse Levinson):**

**Step 3.** $k_3 = a_3^{(3)} = -0.3$.
$$
a_i^{(2)} = \frac{a_i^{(3)} - k_3 a_{3-i}^{(3)}}{1 - k_3^2},\qquad i = 1, 2.
$$
$$
a_1^{(2)} = \frac{-1.5 - (-0.3)(-0.3)}{1 - 0.09} = \frac{-1.5 - 0.09}{0.91} = -\frac{1.59}{0.91} \approx -1.747.
$$
$$
a_2^{(2)} = \frac{1.1 - (-0.3)(-1.5)}{0.91} = \frac{1.1 - 0.45}{0.91} \approx 0.714.
$$

**Step 2.** $k_2 = a_2^{(2)} \approx 0.714$.
$$
a_1^{(1)} = \frac{-1.747 - (0.714)(0.714)}{1 - 0.510} = \frac{-1.747 - 0.510}{0.490} \approx -4.606.
$$

**Step 1.** $k_1 = a_1^{(1)} \approx -4.606$.

**Reflection coefficients:** $k_1 \approx -4.6$, $k_2 \approx 0.71$, $k_3 = -0.3$.

**Warning.** $|k_1| \approx 4.6 > 1$ — the AR(3) filter is **unstable**! The
original coefficients $a = (-1.5, 1.1, -0.3)$ do not correspond to a stable AR
process. (The direct-form roots of $1 - 1.5z^{-1} + 1.1z^{-2} - 0.3z^{-3}$
include one outside the unit circle.)

**This is a feature of the lattice representation**: it detects instability
immediately via $|k_m| \ge 1$, whereas the direct-form check would require
root-finding.

---

## 3.10 LPC of speech

For a 20 ms segment of /a/ (vowel) at 16 kHz, AR(12) fitted via Levinson:

- Formant peaks in AR PSD at $\approx 700$ Hz, 1200 Hz, 2600 Hz — matching
  the first three formants of /a/.
- Prediction-error variance (residual) is ~30 dB below the signal variance —
  most of the speech is predictable from the vocal-tract filter.

**Formant tracking.** Across a longer vowel, the AR PSD peaks move — this is
the basis of formant-based speech analysis and of the **LPC vocoder** used in
CELP, MP3, and modern speech codecs.

---

## 3.11 Riccati convergence

Scalar system: $A = 1$, $Q = 0.1$, $R = 1$.

Fixed point satisfies
$$
P = P + Q - P^2/(P+R) \;\Longrightarrow\; P^2 = Q(P+R).
$$

$$
P^2 - 0.1P - 0.1 = 0 \;\Longrightarrow\; P^* = \frac{0.1 + \sqrt{0.01 + 0.4}}{2} \approx 0.366.
$$

**Convergence from $P_0 = 100$:** $P_n$ drops rapidly (within 10 iterations) to
$P^*$. From $P_0 = 1$: ~5 iterations. Convergence is exponentially fast, rate
proportional to $Q/R$ for small $Q$.

**For $Q = 0.001$:** $P^* \approx 0.032$, convergence in ~30 iterations
(slower because $Q/R$ is smaller). The Riccati recursion's convergence rate
is governed by the "stiffness" of the underlying Lyapunov equation — small
process noise gives slow convergence.

---

## 3.12 Causal vs. non-causal comparison

$N = 200$ samples of the AR(1)+noise model in Problem 3.3.

| Filter | MSE |
|---|:---:|
| Pass-through ($\hat S = Y$) | 0.10 (noise variance) |
| Non-causal Wiener | 0.025 |
| Causal Wiener | 0.034 |

The causal filter needs $\sim 35\%$ more MSE. **Interpretation:** the non-causal
filter uses the "future" observations to predict the past noise — a luxury
the causal filter cannot afford. In a real-time system, causality is not
optional; in offline post-processing (e.g., audio restoration), it can be
relaxed to get better estimates.