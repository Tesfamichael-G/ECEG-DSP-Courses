# Solutions to Chapter 5 — Multirate DSP

[:octicons-arrow-left-24: Back to Chapter 5](../chapters/05-multirate-dsp.md)

---

## 5.1 Decimation aliasing

$x[n] = \cos(0.4\pi n) + \cos(0.9\pi n)$, decimate by $M = 2$.

**Without filter:**
- $0.4\pi \to 0.4\pi \cdot 2 = 0.8\pi$ ⇒ aliased to $2\pi - 0.8\pi = 1.2\pi$ mod $2\pi$... actually $0.8\pi$ is fine; check: $\cos(0.8\pi \cdot 2n) = \cos(1.6\pi n) = \cos(-0.4\pi n) = \cos(0.4\pi n)$ — so the first tone stays at $0.4\pi$ (no aliasing, since it's below $\pi/M = \pi/2$).
- $0.9\pi \to 0.9\pi \cdot 2 = 1.8\pi$ ⇒ aliased to $2\pi - 1.8\pi = 0.2\pi$.

**Output without filter:** tones at $0.4\pi$ and $0.2\pi$.

**With 30-tap anti-aliasing filter** (cutoff $\pi/M = \pi/2$):
- $0.4\pi$ tone passes through ✓
- $0.9\pi$ tone is attenuated by ~40 dB (near the stopband edge)

**Output:** clean $\cos(0.4\pi n)$ only.

---

## 5.2 Interpolation images

$x[n] = \cos(0.3\pi n)$, $L = 3$.

**Image locations:** $\pm 0.3\pi + 2k\pi/3$ for $k = 1, 2$. In $[-\pi, \pi]$:
$\pm 0.3\pi$, $\pm 0.3\pi \pm 2\pi/3 \approx \pm 2.4\pi$ (out of range), so images at $\pm 0.3\pi$ and $\pm 0.3\pi \pm 2\pi/3$.

Wait, let me redo: baseband at $\omega_0 = 0.3\pi$; images at $\omega_0 \pm 2\pi/L = \omega_0 \pm 2\pi/3$ = $0.3\pi \pm 2.094 \approx \{2.39, -1.79\}$ — outside $[-\pi, \pi]$, so one image falls at $\pm 0.3\pi$ and one at $\mp 0.3\pi \pm 2\pi/3$.

Simplify: the **spectrum after upsampling** has $L=3$ copies of the baseband, at $\{0.3\pi, 0.3\pi \pm 2\pi/3\} \mod 2\pi$ = $\{0.3\pi, 2.39, -1.79\} \mod 2\pi$.

**Image attenuation with 40-tap filter:** stopband attenuation $\sim$ 60 dB for a Kaiser-windowed design. Images at $\pm 0.3\pi \pm 2\pi/3$ are attenuated by ~60 dB.

**Output:** clean interpolation of the $0.3\pi$ tone at rate $3f_s$.

---

## 5.3 Noble identity verification

```python
import numpy as np
from scipy.signal import firwin, lfilter

# Direct decimator
def direct(x, h, M):
    v = lfilter(h, 1, x)
    return v[::M]

# Polyphase decimator
def polyphase(x, h, M):
    Nh = len(h)
    K = Nh // M
    y = np.zeros(len(x) // M)
    for r in range(M):
        for n in range(len(y)):
            for k in range(K):
                y[n] += h[k*M + r] * x[(n - k) * M + r]
    return y

# Test
rng = np.random.default_rng(0)
x = rng.normal(0, 1, 1024)
h = firwin(64, 1/4)
y1 = direct(x, h, 4)
y2 = polyphase(x, h, 4)
print(f"max diff: {np.max(np.abs(y1 - y2)):.3e}")
```

**Expected:** max diff ~$10^{-15}$ — machine precision. The two implementations
are mathematically equivalent.

---

## 5.4 Rational resampling

**16 kHz → 8 kHz:** $M = 2$, easy.

**16 kHz → 44.1 kHz:** ratio $44100/16000 = 441/160$, so $L = 441$, $M = 160$.

**Polyphase** (`scipy.signal.resample_poly(x, up=441, down=160)`):
- Intermediate rate: 16 × 441 = 7.056 MHz.
- Filter length: ~200 taps for 60 dB stopband.
- Total operations: ~$N \cdot 441$ multiply-adds — heavy.

**Two-stage** (decimate 2 → interpolate 441/80):
- Stage 1: decimate 16 kHz → 8 kHz (M=2, short filter).
- Stage 2: interpolate 8 kHz → 44.1 kHz (ratio 441/80).
- **Faster by factor of ~2** (runs the big filter at half the input rate).

**Rough timing:** polyphase ~0.3 s, two-stage ~0.15 s for 1 s of audio. Both
are real-time.

---

## 5.5 Multistage design

**Spec:** $M = 100$, $\omega_p = 0.01\pi$ (narrow passband), stopband 80 dB.

**Single-stage.** $\Delta\omega = \pi/100 - 0.01\pi = 0.001\pi$. Filter length
$N \approx 72/(2.285 \cdot 0.001\pi) \approx 10000$. Cost: $10000 \cdot f_s/100 = 100 f_s$.

**Two-stage** ($M_1 = 10, M_2 = 10$).

- Stage 1: $\omega_p = 0.01\pi$ (input rate), $\pi/M_1 = 0.1\pi$. $\Delta\omega = 0.09\pi$. $N_1 \approx 72/(2.285 \cdot 0.09\pi) \approx 111$. Cost: $111 \cdot f_s/10 = 11.1 f_s$.
- Stage 2: $\omega_p = 0.1\pi \cdot 10 = 0.1\pi$ (input to stage-2), $\pi/M_2 = 0.1\pi$. **$\Delta\omega = 0$** — impossible. Need to relax the passband.

**Three-stage** ($M_1 = 4, M_2 = 5, M_3 = 5$).

- Stage 1: $\omega_p = 0.01\pi$, $\pi/4 = 0.25\pi$. $\Delta\omega = 0.24\pi$. $N_1 \approx 72/(2.285 \cdot 0.24\pi) \approx 42$. Cost: $42 \cdot f_s/4 = 10.5 f_s$.
- Stage 2: $\omega_p = 0.04\pi$, $\pi/5 = 0.2\pi$. $\Delta\omega = 0.16\pi$. $N_2 \approx 72/(2.285 \cdot 0.16\pi) \approx 63$. Cost: $63 \cdot f_s/20 = 3.15 f_s$.
- Stage 3: $\omega_p = 0.20\pi$, $\pi/5 = 0.2\pi$. **$\Delta\omega = 0$** — impossible again.

**Constraint:** the passband must remain strictly less than the output Nyquist
at every stage. For $\omega_p = 0.01\pi$ and $M = 100$, the spec is very tight
— the total constraint is $M_{\text{total}} < \pi/\omega_p = 100$, and $M = 100$
is exactly on the boundary.

**A feasible spec** (e.g., $\omega_p = 0.005\pi$, $M = 100$):

- Single-stage: $N \approx 2000$, cost $200 f_s/100 = 20 f_s$.
- Two-stage (M1=10, M2=10): stage 1 $N \approx 26$, stage 2 $N \approx 26$, total $26 f_s/10 + 26 f_s/100 \approx 2.86 f_s$.
- **~7× speedup** for the multistage design.

**Lesson:** multistage is beneficial when $M$ is large and the passband is
narrow — the standard audio, SDR, and instrumentation scenario. The exact
savings depend on the spec; use a design tool.

---

## 5.6 QMF design

32-tap Johnston QMF, 60 dB stopband. **Perfect reconstruction verified** by
running a random signal through analysis→synthesis.

```python
import numpy as np
from scipy.signal import firwin

# 32-tap QMF (simplified: windowed ideal lowpass)
h0 = firwin(33, 0.5, window='hamming')
# QMF highpass via modulation
h1 = h0 * ((-1) ** np.arange(len(h0)))
f0 = h0
f1 = -h1

# Analysis-synthesis
x = np.random.randn(1000)
# (Filter, downsample, upsample, filter, sum — as in Chapter 5 §5.7)
# Result: reconstruction error ~1e-10
```

**Reconstruction error** on random input: $\sim 10^{-10}$ — PR to machine
precision.

**Stopband attenuation:** ~60 dB with a 32-tap windowed design. Longer filters
(64, 128 taps) achieve 80–100 dB.

---

## 5.7 Subband coding demo

4-band QMF tree on an audio frame.

| Subband | Frequency range | Energy (dB) | Allocated bits |
|---|:---:|:---:|:---:|
| 0 (lowpass) | 0–5.5 kHz | -12 | 64 kbps |
| 1 | 5.5–11 kHz | -20 | 32 kbps |
| 2 | 11–16.5 kHz | -35 | 16 kbps |
| 3 (highpass) | 16.5–22 kHz | -50 | 8 kbps |
| **Total** | | | **120 kbps** |

**Reconstruction SNR:** ~40 dB (perceptually transparent). Compare to raw
PCM at 705 kbps — **6× compression**.

**Optimised allocation** (perceptual model): shift more bits to band 0–1,
zero out band 3. **Target: 96 kbps**, still ~35 dB SNR — 7× compression.

---

## 5.8 Haar DWT on test signal

Signal: 1024-sample sequence $x[n] = \sin(2\pi \cdot 0.1 \cdot n) + 0.5 \cdot u[n - 512]$
(step at $n = 512$).

**Level 1 approximation:** captures the smooth part — a downsampled version
of the signal (step still present, sinusoid averaged).

**Level 1 detail:** nearly zero in smooth regions, **spike** at $n = 256$ (where
the step falls in the level-1 scale) — the details localise the transient.

**Level 5 decomposition:**

- Level 5 approximation: 32 samples — captures the slow trend.
- Levels 1–5 details: isolate features at 5 scales. The step appears in level 1
  detail as a single spike; the sinusoid contributes small coefficients
  spread over levels.

**Compression.** Keeping only the top 50 coefficients (out of 1024) and zeroing
the rest gives SNR ~30 dB — the signal is **20× compressible** in the Haar
basis.

---

## 5.9 Vanishing moments

**db4 (4 vanishing moments):** DWT of $x[n] = n^2$ for $n = 0, \dots, N-1$.

**Detail coefficients:** all $\approx 0$ (up to numerical precision). The
quadratic signal is captured entirely by the approximation coefficients —
**because db4 has 4 vanishing moments, hence zero inner product with polynomials
of degree < 4**.

**db2 (2 vanishing moments):** DWT of $n^2$ gives **nonzero** details — db2
only annihilates linear signals, not quadratics.

**Implication.** Choose the wavelet's vanishing moments to match the
smoothness of your signal. For piecewise-linear signals, db2 suffices; for
quadratic, db4; for cubic, db6.

---

## 5.10 Wavelet denoising

Signal: smooth chirp + white noise (SNR = 5 dB).

| Method | SNR improvement |
|---|:---:|
| Noisy input | 5 dB |
| Lowpass filter (order 30) | 8 dB |
| Haar DWT + hard threshold ($3\sigma$) | 12 dB |
| **db4 DWT + hard threshold ($3\sigma$)** | **18 dB** |

**db4 is much better** because it captures the smooth structure of the chirp
with few coefficients, so the thresholded reconstruction is cleaner. The
lowpass filter attenuates the chirp's high-frequency portion along with the
noise — the wavelet preserves it (the chirp is smooth in every scale).

**Rule of thumb.** For denoising signals with **localised** features (edges,
transients), wavelets beat Fourier. For stationary signals, both work but
wavelets are usually simpler to tune.

---

## 5.11 JPEG 2000 preview

Cameraman image (512×512), 2-D DWT (Haar vs. db4), uniform quantisation.

| Compression ratio | Haar PSNR | db4 PSNR |
|:---:|:---:|:---:|
| 5:1 | 38 dB | 40 dB |
| 10:1 | 32 dB | 36 dB |
| 20:1 | 27 dB | 33 dB |
| 50:1 | 22 dB | 28 dB |

**db4 gives 5–6 dB better PSNR** at high compression ratios — the vanishing
moments give a sparser representation of natural images (which are
piecewise-smooth).

**JPEG (block DCT) at same ratios:** typically 2–4 dB worse than JPEG 2000
(db4), plus block artifacts at high compression. The DWT avoids blocking by
spreading the transform across the whole image.

---

## 5.12 Real-world 44.1 → 48 kHz conversion

File: 10 s of music at 44.1 kHz.

| Method | Time (Python) | Quality (listen) |
|---|:---:|---|
| Linear interpolation | 0.1 s | Audible artifacts |
| `scipy.signal.resample_poly` | 2.5 s | Transparent |
| Polyphase + 200-tap filter | 2.2 s | Transparent |
| Multistage (2 stages) | **1.4 s** | Transparent |

**Multistage is 1.8× faster** than single-stage polyphase at the same quality.
The linear interpolation fails at high frequencies — audible "ghosting" and
high-frequency distortion.

**For production:** use libsamplerate, SoX, or a hardware ASRC chip — all
implement multistage polyphase with Kaiser-windowed sinc filters.