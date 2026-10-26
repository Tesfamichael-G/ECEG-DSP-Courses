# DSP Primer

(Placeholder — this page will hold shared terminology, notation, and
transform tables for ECEG-6303 and ECEG-6423.)# DSP Primer

Shared reference for **ECEG-6303 (Advanced DSP)** and **ECEG-6423 (Multimedia
Signal Processing)**. Notation conventions, transform definitions, and key
formulas used throughout both courses.

---

## 1. Notation

| Symbol | Meaning |
|---|---|
| $x[n]$ | Discrete-time signal, $n \in \mathbb Z$ |
| $x(t)$ | Continuous-time signal, $t \in \mathbb R$ |
| $X(e^{j\omega})$ | DTFT of $x[n]$, $\omega \in [-\pi, \pi)$ |
| $X[k]$ | $N$-point DFT of $x[n]$, $k = 0, \ldots, N-1$ |
| $X(z)$ | Z-transform of $x[n]$ |
| $X(s)$ | Laplace transform of $x(t)$ |
| $f_s$ | Sampling frequency (Hz) |
| $T_s = 1/f_s$ | Sampling interval |
| $\omega = 2\pi f T_s$ | Normalised angular frequency |
| $h[n], H(z)$ | Impulse response, transfer function |
| $R_X[k]$ | Autocorrelation of a WSS sequence |
| $S_X(e^{j\omega})$ | Power spectral density |
| $*$ | Convolution |
| $\otimes$ | Kronecker product |

**Conventions.**

- **Discrete time** uses square brackets: $x[n]$.
- **Continuous time** uses parentheses: $x(t)$.
- **Frequency** in Hz is $f$; in rad/sample is $\omega$; $f = \omega/(2\pi T_s)$.
- **Fourier transform** uses $e^{-j2\pi ft}$ (negative exponent for forward).

---

## 2. Fourier transforms

### Continuous-time Fourier transform (CTFT)

$$
\boxed{\;X(f) = \int_{-\infty}^{\infty} x(t)\, e^{-j2\pi ft}\, dt,\qquad x(t) = \int_{-\infty}^{\infty} X(f)\, e^{j2\pi ft}\, df.\;}
$$

### Discrete-time Fourier transform (DTFT)

$$
\boxed{\;X(e^{j\omega}) = \sum_{n=-\infty}^{\infty} x[n]\, e^{-j\omega n},\qquad x[n] = \frac{1}{2\pi}\int_{-\pi}^{\pi} X(e^{j\omega})\, e^{j\omega n}\, d\omega.\;}
$$

$X(e^{j\omega})$ is $2\pi$-periodic in $\omega$.

### $N$-point discrete Fourier transform (DFT)

$$
\boxed{\;X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-j2\pi kn/N},\qquad x[n] = \frac{1}{N}\sum_{k=0}^{N-1} X[k]\, e^{j2\pi kn/N}.\;}
$$

The DFT is the DTFT **sampled at** $\omega_k = 2\pi k/N$.

---

## 3. Z-transform and Laplace

### Z-transform

$$
\boxed{\;X(z) = \sum_{n=-\infty}^{\infty} x[n]\, z^{-n},\qquad z \in \mathbb C.\;}
$$

**ROC (region of convergence)** is the annulus in the $z$-plane where the sum converges. For causal stable systems, all poles must lie **inside** the unit circle $|z| < 1$; the ROC must contain $|z| = 1$.

### Laplace transform

$$
\boxed{\;X(s) = \int_{-\infty}^{\infty} x(t)\, e^{-st}\, dt,\qquad s = \sigma + j\Omega.\;}
$$

ROC is a vertical strip in the $s$-plane; for causal stable systems, it must contain the imaginary axis.

**Relationship:** $z = e^{sT_s}$, so a pole at $s = p$ maps to $z = e^{pT_s}$. The left half-plane (LHP) maps to the interior of the unit disk.

---

## 4. Sampling theorem

A band-limited signal $x(t)$ with $X(f) = 0$ for $|f| > W$ is exactly recoverable from samples $x(n T_s)$ iff

$$
\boxed{\;T_s < \frac{1}{2W} \quad \Longleftrightarrow \quad f_s > 2W.\;}
$$

Recovery by sinc interpolation:

$$
\boxed{\;x(t) = \sum_{n=-\infty}^{\infty} x(nT_s)\, \operatorname{sinc}\!\left(\frac{t - nT_s}{T_s}\right),\quad \operatorname{sinc}(u) = \frac{\sin(\pi u)}{\pi u}.\;}
$$

**2-D:** apply independently to $x$ and $y$. **3-D:** apply independently to $x, y, t$.

---

## 5. Common transform pairs

| $x(t)$ or $x[n]$ | $X(f)$ or $X(e^{j\omega})$ |
|---|---|
| $\delta(t)$ | $1$ |
| $\delta[n]$ | $1$ |
| $u(t)$ | $\frac{1}{j2\pi f} + \frac{1}{2}\delta(f)$ |
| $u[n]$ | $\frac{1}{1 - e^{-j\omega}} + \pi\delta(\omega)$ (for convergence) |
| $e^{j\omega_0 n}$ | $2\pi \sum_k \delta(\omega - \omega_0 - 2\pi k)$ |
| $\cos(2\pi f_0 t)$ | $\frac{1}{2}[\delta(f - f_0) + \delta(f + f_0)]$ |
| $e^{-at}u(t)$, $a > 0$ | $\frac{1}{a + j2\pi f}$ |
| $a^n u[n]$, $|a| < 1$ | $\frac{1}{1 - a e^{-j\omega}}$ |
| Rectangular pulse | Sinc |
| Gaussian | Gaussian (self-dual) |
| $\frac{\sin(W t)}{\pi t}$ | Rectangular (band-limited) |

---

## 6. Convolution

### Continuous

$$
\boxed{\;(x * y)(t) = \int_{-\infty}^{\infty} x(\tau)\, y(t - \tau)\, d\tau,\qquad \mathcal F\{x * y\} = X(f) \cdot Y(f).\;}
$$

### Discrete

$$
\boxed{\;(x * y)[n] = \sum_{k=-\infty}^{\infty} x[k]\, y[n - k],\qquad \mathcal F\{x * y\} = X(e^{j\omega}) \cdot Y(e^{j\omega}).\;}
$$

### Circular (DFT domain)

$$
\boxed{\;(x \circledast y)[n] = \sum_{k=0}^{N-1} x[k]\, y[(n - k) \bmod N].\;}
$$

Circular convolution equals linear convolution when both sequences have length $\le N/2$ and are zero-padded to length $N$. Otherwise there is **time-domain aliasing**.

---

## 7. Key theorems

| Theorem | Statement |
|---|---|
| **Linearity** | $\mathcal T\{ax_1 + bx_2\} = a\mathcal T\{x_1\} + b\mathcal T\{x_2\}$ |
| **Time shift** | $\mathcal F\{x(t - t_0)\} = X(f)\, e^{-j2\pi f t_0}$ |
| **Frequency shift** | $\mathcal F\{x(t)\, e^{j2\pi f_0 t}\} = X(f - f_0)$ |
| **Scaling** | $\mathcal F\{x(at)\} = \frac{1}{a}X(f/a)$ |
| **Parseval (CT)** | $\int |x(t)|^2 dt = \int |X(f)|^2 df$ |
| **Parseval (DT)** | $\sum |x[n]|^2 = \frac{1}{2\pi}\int |X(e^{j\omega})|^2 d\omega$ |
| **Convolution** | $\mathcal F\{x * y\} = X \cdot Y$ |
| **Multiplication** | $\mathcal F\{x \cdot y\} = X * Y$ (convolution in frequency) |
| **Wiener–Khinchine** | $\mathcal F\{R_X\} = S_X$ (PSD is DTFT of autocorrelation) |

---

## 8. Filters

### FIR — finite impulse response

$$
\boxed{\;y[n] = \sum_{k=0}^{N-1} h[k]\, x[n-k],\qquad H(z) = \sum_{k=0}^{N-1} h[k]\, z^{-k}.\;}
$$

**Properties.** Always stable (all poles at $z=0$); linear phase possible (symmetric $h$); no feedback.

**Design.** Window method (Hamming, Kaiser); Parks–McClellan (equiripple optimal); frequency sampling.

### IIR — infinite impulse response

$$
\boxed{\;\sum_{k=0}^{M} a_k\, y[n-k] = \sum_{k=0}^{N} b_k\, x[n-k],\qquad H(z) = \frac{\sum_{k=0}^{N} b_k z^{-k}}{\sum_{k=0}^{M} a_k z^{-k}}.\;}
$$

**Properties.** More efficient (fewer coefficients for a given spec); can be unstable; nonlinear phase (unless specially designed).

**Design.** Analog prototype (Butterworth, Chebyshev, elliptic) + bilinear transform or impulse-invariant.

**Stability.** All poles inside $|z| < 1$.

---

## 9. Random signals

| Quantity | Definition |
|---|---|
| **Mean** | $\mu_X = \mathbb E[X[n]]$ |
| **Variance** | $\sigma_X^2 = \mathbb E[(X[n] - \mu_X)^2]$ |
| **Autocorrelation** | $R_X[k] = \mathbb E[X[n+k] X^*[n]]$ |
| **Autocovariance** | $C_X[k] = R_X[k] - |\mu_X|^2$ |
| **WSS** | $\mu_X$ constant, $R_X[k]$ depends only on lag |
| **PSD** | $S_X(e^{j\omega}) = \sum R_X[k] e^{-j\omega k}$ |
| **Cross-correlation** | $R_{XY}[k] = \mathbb E[X[n+k] Y^*[n]]$ |

**Key properties.**

- $R_X[0] = \mathbb E[|X|^2]$ — total power.
- $R_X[-k] = R_X^*[k]$ — Hermitian symmetry.
- $|R_X[k]| \le R_X[0]$ — max at lag 0.
- $S_X \ge 0$ — PSD is non-negative.

---

## 10. Useful constants

| Constant | Value |
|---|---|
| $\pi$ | 3.14159… |
| $e$ | 2.71828… |
| Speed of sound in air (20°C) | ≈ 343 m/s |
| Speed of sound in warm humid air (voice) | ≈ 350 m/s |
| Speed of light | $3 \times 10^8$ m/s |
| Human hearing range | 20 Hz – 20 kHz |
| Human speech bandwidth | 50 Hz – 8 kHz (most energy < 4 kHz) |
| Human visual acuity | ~1 arcmin |
| Human flicker fusion | ~50–60 Hz |
| Cinema frame rate | 24 fps |
| PAL / NTSC | 25 / 29.97 fps |

---

## 11. Common pitfalls

- **Confusing DFT frequency indices with Hz.** DFT bin $k$ corresponds to $f = k f_s/N$ Hz. Normalised frequency $\omega_k = 2\pi k/N$.
- **Forgetting the $1/N$ in the inverse DFT.** Or using the wrong convention.
- **Ignoring the ROC in Z-transform.** A transfer function without the ROC does not uniquely determine the time-domain sequence.
- **Confusing power and energy.** Energy signals have finite $\sum |x|^2$; power signals have finite $\lim \frac{1}{N}\sum |x|^2$. WSS processes are power signals.
- **Mixing continuous and discrete frequency.** $\omega$ is rad/sample; $\Omega = 2\pi f$ is rad/s. They differ by factor $T_s$.
- **Forgetting that the FFT is not the same as the DFT.** The DFT is the mathematical operation; the FFT is an algorithm for computing it.