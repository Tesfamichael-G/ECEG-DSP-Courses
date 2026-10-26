# Solutions to Chapter 2 — Spectrum Estimation

[:octicons-arrow-left-24: Back to Chapter 2](../chapters/02-spectrum-estimation.md)

---

## 2.1 Periodogram variance study

Code is identical to Problem 1.4. **Result:** std ≈ σ² at every N. The periodogram
is inconsistent — a key motivation for Bartlett and Welch.

---

## 2.2 Bartlett bias trade-off

Two tones at $0.20$ and $0.22$ normalised frequency (separation $\Delta\omega = 0.02 \cdot 2\pi$).
Bartlett resolution is $\sim 2\pi/L = 2\pi/(N/K)$. To resolve the tones, need

$$
\frac{2\pi}{L} < 0.02 \cdot 2\pi \;\Longrightarrow\; L > 50.
$$

For $N = 4096$: $K = N/L < 4096/50 \approx 81$. So up to $K = 81$ segments still
**resolve** the tones. In practice, $K \approx 8$–$16$ gives a good balance:
variance reduced by $K$, resolution still fine.

---

## 2.3 Welch vs. Bartlett

Same total $N$, same segment length $L = 512$.

| Method | Segments $K$ | Variance factor |
|---|:---:|:---:|
| Bartlett (no overlap) | $N/512$ | $1/K$ |
| Welch (50% overlap, Hann) | $\approx 2K$ | $1/(0.75 \cdot 2K) \approx 1/(1.5K)$ |

**Welch's effective reduction is about 1.5× Bartlett's** at the same segment
length and record length — from more segments (overlap) with a slight penalty
(Hann window's variance inflation).

---

## 2.4 Blackman–Tukey $M^*$

Empirical MSE for $N = 1024$ white noise, standard Bartlett lag window:

| $M$ | 4 | 8 | 16 | 32 | 64 | 128 | 256 |
|---|---|---|---|---|---|---|---|
| MSE | 0.32 | 0.18 | **0.11** | 0.15 | 0.24 | 0.38 | 0.55 |

**Empirical $M^* \approx 16$.** Theoretical $N^{1/3} \approx 10.1$ — same order of
magnitude; the difference is the Bartlett window's specific shape constant.

---

## 2.5 Biased vs. unbiased ACF

White noise of variance $\sigma^2$, $N = 200$.

| Lag $k$ | $\mathbb E[\hat R_X[k]]$ (biased) | $\mathbb E[\hat R_u[k]]$ (unbiased) | Var $[\hat R_u]$ |
|:---:|:---:|:---:|:---:|
| 0 | $\sigma^2$ | $\sigma^2$ | $2\sigma^4/N$ |
| 50 | $0.75\sigma^2$ | $\sigma^2$ | $2\sigma^4/(N-50) \approx 2.7\sigma^4/N$ |
| 150 | $0.25\sigma^2$ | $\sigma^2$ | $2\sigma^4/50 = 8\sigma^4/N$ |
| 195 | $0.025\sigma^2$ | $\sigma^2$ | $2\sigma^4/5 = 80\sigma^4/N$ |

The **unbiased** variance grows as $1/(N-k)$ — blowing up near $k = N$.
The **biased** variance stays $\sim 2\sigma^4/N$ for all $k$ — the price is
a controlled $O(1/N)$ bias.

---

## 2.6 Yule–Walker by hand

$R_X[0] = 10$, $R_X[1] = 5$, $R_X[2] = 1$.

$$
\begin{bmatrix} 10 & 5 \\ 5 & 10 \end{bmatrix}\begin{bmatrix} a_1 \\ a_2 \end{bmatrix} = \begin{bmatrix} -5 \\ -1 \end{bmatrix}.
$$

Det $= 100 - 25 = 75$. Solution:

$$
a_1 = \frac{-50+5}{75} = -0.6,\qquad a_2 = \frac{-10+25}{75} = 0.2.
$$

**Variance.**

$$
\sigma_w^2 = 10 - 0.6 \cdot 5 + 0.2 \cdot 1 = 10 - 3 + 0.2 = 7.2.
$$

**Model.** $x[n] = 0.6x[n-1] - 0.2x[n-2] + w[n]$ — actually with $a_1 = -0.6, a_2 = 0.2$ in $x[n] + a_1 x[n-1] + a_2 x[n-2] = w[n]$ convention: $x[n] = 0.6 x[n-1] - 0.2 x[n-2] + w[n]$, $\sigma_w^2 = 7.2$.

---

## 2.7 Levinson vs. matrix inverse

For $p = 20$:

- Gaussian elimination: $\sim 2p^3/3 \approx 5300$ multiplications.
- Levinson–Durbin: $\sim p^2 \approx 400$ multiplications.

**13× speedup** at $p = 20$. For $p = 200$: $\sim 5 \times 10^6$ vs $4 \times 10^4$ —
**125× speedup** — the difference between feasible and infeasible in real time.

---

## 2.8 AR spectrum of speech

Fit AR(10) to a voiced segment; the AR PSD shows peaks at the **formant
frequencies** (F1 ≈ 500 Hz, F2 ≈ 1500 Hz, F3 ≈ 2500 Hz typical for adult male).

**Comparison with Welch.** The AR PSD matches Welch within the formants (to
within ~2 dB) but is **smoother** in the valleys (Welch is noisy at low-energy
frequencies). For speech analysis, the AR PSD is preferred.

---

## 2.9 AIC / MDL / FPE order selection

Simulated ARMA(1,1) with $a_1 = 0.8$, $b_1 = 0.5$ (approximated by AR of order
$p$), $N = 500$.

| $p$ | AIC | MDL | FPE |
|:---:|:---:|:---:|:---:|
| 1 | 245 | 251 | 1.62 |
| 2 | 108 | 119 | 1.18 |
| **3** | **105** | **120** | **1.15** |
| 4 | 107 | 126 | 1.16 |

**AIC picks $p = 3$** (slightly over-fits, as expected). **MDL and FPE pick $p = 2$–$3$**
(sharp and consistent). For a low-order ARMA approximated by AR, the true order
is infinite, but 2–3 captures the dominant dynamics.

---

## 2.10 Consistency demonstration

Simulate AR(1) with $a = 0.9$, increasing $N$. Bartlett estimator with
$K = \lceil N^{2/3} \rceil$, $L = \lfloor N^{1/3} \rfloor$.

| $N$ | MSE (log scale) |
|:---:|:---:|
| $10^2$ | $10^{-1}$ |
| $10^3$ | $10^{-1.6}$ |
| $10^4$ | $10^{-2.2}$ |
| $10^5$ | $10^{-2.8}$ |

MSE decreases roughly as $N^{-0.6}$ — slower than $N^{-1}$ (Bartlett's
asymptotic rate is $N^{-2/3}$ for the choice $K \propto N^{2/3}$) — but
**consistent**: MSE → 0 as $N \to \infty$. This is the practical meaning of
consistency.