# Solutions to Chapter 1 — Discrete Random Signal Processing

[:octicons-arrow-left-24: Back to Chapter 1](../chapters/01-discrete-random-signal-processing.md)

---

## 1.1 Parseval for a geometric sequence

$x[n] = a^n u[n]$, $|a| < 1$.

**Energy.** $E = \sum |a|^{2n} = \frac{1}{1 - a^2}$.

**DTFT.** $X(e^{j\omega}) = \frac{1}{1 - ae^{-j\omega}}$, so $|X|^2 = \frac{1}{(1-ae^{-j\omega})(1-ae^{j\omega})} = \frac{1}{1-2a\cos\omega+a^2}$.

**Parseval check.**
$$
\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{1}{1-2a\cos\omega+a^2}\,d\omega = \frac{1}{1-a^2}
$$
(standard Poisson-kernel identity). ✓

---

## 1.2 Autocorrelation of random-phase sinusoid

$X[n] = A\cos(\omega_0 n + \Theta)$, $\Theta \sim U[0, 2\pi)$.

$$
R_X[k] = \mathbb E[A\cos(\omega_0 n + \Theta)A\cos(\omega_0(n+k) + \Theta)] = \tfrac{A^2}{2}\cos(\omega_0 k).
$$

**Same form as the continuous-time case** — sampling does not change the
autocorrelation structure, only the parameter (continuous $t$ → discrete $k$).

---

## 1.3 Wiener–Khinchine inverse

$S_X(e^{j\omega}) = \dfrac{1-a^2}{|1-ae^{-j\omega}|^2}$, $|a| < 1$.

Inverse DTFT:
$$
R_X[k] = \frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{1-a^2}{1-2a\cos\omega+a^2}e^{j\omega k}\,d\omega = \rho^{|k|}
$$
with $\rho = a$. Verified by direct contour integral.

---

## 1.4 Periodogram variance study

```python
import numpy as np
from scipy.signal import periodogram

rng = np.random.default_rng(0)
sigma2 = 1.0
for N in [256, 1024, 4096, 16384]:
    x = rng.normal(0, np.sqrt(sigma2), N)
    _, S = periodogram(x, fs=2*np.pi, window='boxcar', nfft=N)
    S = S / (2*np.pi)                      # normalise to rad/sample
    print(f"N={N:6d}  mean={S.mean():.3f}  std={S.std():.3f}")
```

**Expected output:** mean ≈ 1, std ≈ 1 **for every N**. The std does **not**
shrink — the periodogram is inconsistent.

---

## 1.5 Sum decomposition for AR(1)

$R_X[k] = \rho^{|k|}$, $N = 16$; compute eigenvectors of the $16\times 16$ Toeplitz
covariance matrix.

**Eigenvalue decay:**

| $\rho$ | $\lambda_1$ | $\lambda_2$ | $\lambda_5$ | Components for 99% |
|:---:|:---:|:---:|:---:|:---:|
| 0.9 | ≈ 13.0 | ≈ 2.3 | ≈ 0.12 | ~5 |
| 0.5 | ≈ 2.0 | ≈ 1.5 | ≈ 0.5 | ~14 |

For $\rho = 0.9$ (highly correlated), the first 5 components capture 99% of
the variance — a 16:5 compression ratio. For $\rho = 0.5$ (weakly correlated),
almost all 16 components are needed.

---

## 1.6 Spectral factorization

$$
S_X(z) = \frac{(1-0.3z^{-1})(1-0.3z)}{(1-0.9z^{-1})(1-0.9z)}.
$$

**Minimum-phase factor.** Interior roots at $z = 0.3$ (numerator) and $z = 0.9$
(denominator):

$$
H(z) = \frac{1-0.3z^{-1}}{1-0.9z^{-1}},\qquad \sigma_e^2 = 1.
$$

**Whitening filter.** $1/H(z) = \dfrac{1-0.9z^{-1}}{1-0.3z^{-1}}$ — causal, stable.

---

## 1.7 White noise through $H(z) = 1 + 0.5z^{-1}$

$X$ white of variance 1.

**Output PSD.**
$$
S_Y(e^{j\omega}) = |1 + 0.5e^{-j\omega}|^2 = 1.25 + \cos\omega.
$$

**Autocorrelation.**
$$
R_Y[0] = 1.25,\quad R_Y[\pm 1] = 0.5,\quad R_Y[k] = 0\ (|k| \ge 2).
$$

**Variance:** $R_Y[0] = 1.25$ — the filter amplifies noise power by $1 + 0.5^2 = 1.25$.

---

## 1.8 Box–Muller statistics

```python
import numpy as np
from scipy import stats

def box_muller(n, seed=0):
    rng = np.random.default_rng(seed)
    U1, U2 = rng.uniform(0, 1, n), rng.uniform(0, 1, n)
    R = np.sqrt(-2 * np.log(U1))
    return np.concatenate([R*np.cos(2*np.pi*U2), R*np.sin(2*np.pi*U2)])

z = box_muller(100_000)
print(f"mean={z.mean():.4f}  var={z.var():.4f}")
print(f"skew={stats.skew(z):.4f}  kurt={stats.kurtosis(z, fisher=False):.4f}")
```

**Expected:** mean ≈ 0, var ≈ 1, skew ≈ 0, kurtosis ≈ 3 — the four moments of
$\mathcal N(0,1)$.

---

## 1.9 Colouring exercise

AR(2) filter $H(z) = 1/(1 - 1.3z^{-1} + 0.42z^{-2})$.

**Theoretical PSD.**
$$
S_X(e^{j\omega}) = \frac{\sigma_w^2}{|1 - 1.3e^{-j\omega} + 0.42e^{-2j\omega}|^2}.
$$

**Empirical** — Welch on simulated output, 4 segments of 256, Hann window. The
empirical curve overlays the theoretical (to within 1 dB) when the simulation
length is $\gg 1/(1 - r^2)$ where $r = 0.7$ is the dominant pole.

---

## 1.10 Whitening speech

Fit AR(10) to a voiced-speech segment at 16 kHz. The AR coefficients
$\{a_k\}$ give the whitening filter $1 + \sum a_k z^{-k}$.

**Test of whiteness:** compute sample autocorrelation of the residual
$e[n] = x[n] + \sum a_k x[n-k]$. For $k \ne 0$, $|\hat R_e[k]|$ should be below
$3/\sqrt{N}$ for 99% of lags — consistent with white noise.

**Interpretation.** The residual is the **excitation** of the source-filter
model — the pulse train (voiced) or noise (unvoiced). The AR coefficients
capture the vocal-tract filter, whose spectral envelope is the formant
structure. This is the basis of LPC analysis.