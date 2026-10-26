# Chapter 4 — Time and Frequency Domain Speech Processing Methods

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 2 tutorial + 1 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** significance of short-time analysis · time-domain
    parameters (energy, magnitude, zero-crossing rate, autocorrelation function,
    pitch estimation by parallel processing) · frequency-domain methods (short-time
    Fourier analysis, filter bank analysis, spectrographic analysis, formant
    extraction, pitch extraction, analysis–synthesis systems) · homomorphic
    analysis and LPC of speech (cepstral analysis, formant and pitch estimation,
    formulation of linear prediction problem in time domain).

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Explain **why** speech must be analyzed in short frames and derive the
   windowing / framing scheme for the standard frame rate (50–100 frames/s).
2. Compute the **short-time energy**, **magnitude**, **zero-crossing rate**,
   and **autocorrelation function**; use them for voiced/unvoiced detection,
   endpoint detection, and pitch estimation.
3. Derive and compute the **short-time Fourier transform (STFT)** and the
   **short-time autocorrelation**; interpret the STFT as a filter bank with
   narrowband and wideband trade-offs.
4. Read and interpret a **spectrogram** — formants, pitch, fricatives, plosives.
5. Design a **filter-bank** front end (uniform or mel-spaced) and compute the
   **sub-band energies**.
6. Extract **formants** by peak picking or LPC root finding; estimate **pitch**
   from the **autocorrelation**, **cepstrum**, and **harmonic product spectrum**.
7. Describe a **channel vocoder** (analysis–synthesis system) and implement one.
8. Derive the **real cepstrum** and use it to separate **pitch** (low-quefrency)
   from **formants** (high-quefrency) in speech.
9. Formulate the **linear prediction (LPC)** problem in the time domain, derive
   the **Yule–Walker equations**, and solve them by **Levinson–Durbin** in
   $O(p^2)$.
10. Extract **formants** and **pitch** from LPC coefficients.
11. Implement all of the above in **Python, Java, and C#**.

---

## 4.1 Significance of short-time analysis

### 4.1.1 Why frame-based?

Speech is **non-stationary** — its statistics change continuously as the
articulators move. But over a **short window**, it is **quasi-stationary**:

| Window length | Assumption | Trade-off |
|---|---|---|
| < 5 ms | One pitch period | Too short for reliable spectrum |
| **10–30 ms** | **Quasi-stationary** | **Standard choice** |
| > 50 ms | Two-three pitch periods | Smears transitions between phonemes |

**The rule of thumb.** Use a window of **20–30 ms** (e.g., 256 or 400 samples
at 16 kHz). This covers **2–3 pitch periods** for a typical male voice
($F_0 = 100$ Hz, $T_0 = 10$ ms) and **4–8 periods** for a female voice
($F_0 = 200$ Hz, $T_0 = 5$ ms). Long enough for a stable spectrum, short
enough to freeze the articulators.

**Overlap.** Advance the frame by **half its length** (50% overlap). This gives
a **frame rate** of $f_s / (N/2)$ — typically 50–100 frames/s. Successive frames
share half their samples, which:
- **Smooths transitions** — the feature stream is continuous.
- **Reduces quantization artifacts** — the vocoder reconstruction has fewer
  frame-boundary glitches.
- **Matches perception** — the ear integrates over ~20 ms.

### 4.1.2 Framing and windowing

**Framing** — split the signal into short segments. With frame length $N$ and
hop size $H$ (usually $N/2$), the $m$-th frame is

$$
\boxed{\;x_m[n] = x[n + mH],\quad n = 0, \dots, N-1.\;}
$$

**Windowing** — multiply each frame by a smooth window to reduce **spectral
leakage** at the frame boundaries.

$$
\boxed{\;x_w[m, n] = x_m[n]\, w[n],\quad n = 0, \dots, N-1.\;}
$$

**Common windows.**

| Window | $w[n]$ | Mainlobe width | First sidelobe |
|---|---|:---:|:---:|
| **Rectangular** | $1$ | $4\pi/N$ | $-13$ dB |
| **Hamming** | $0.54 - 0.46\cos(2\pi n/(N-1))$ | $8\pi/N$ | **$-43$ dB** |
| **Hann** | $0.5 - 0.5\cos(2\pi n/(N-1))$ | $8\pi/N$ | $-31$ dB |
| **Blackman** | $0.42 - 0.5\cos(2\pi n/(N-1)) + 0.08\cos(4\pi n/(N-1))$ | $12\pi/N$ | $-58$ dB |

!!! note "Which window?"
    **Hamming** is the standard for STFT-based analysis and vocoders. It has
    $-43$ dB first sidelobe (no visible leakage in most analyses) with mainlobe
    only 2× wider than rectangular. **Hann** is also common and slightly
    better for Welch spectrum estimation. **Rectangular** is only used when
    the frame is exactly an integer number of pitch periods (rare).

### Worked example 4.1 — frame count and coverage

Speech sampled at $f_s = 16$ kHz. Frame $N = 400$ samples (25 ms), hop $H = 200$
samples (12.5 ms, 50% overlap).

**Frame rate:** $f_s/H = 16000/200 = 80$ frames/s.

**Samples per second of speech used (after windowing):** $N \times f_s/H = 400 \times 80 = 32000$ samples — 2× the input samples, because of the 50% overlap.

**A 5-second utterance:** $5 \times 80 = 400$ frames.

**Total features per second (12 LPC coefficients per frame):** $80 \times 12 = 960$ coefficients/s. At 4 bytes/coefficient = 3.84 kB/s — compact for storage.

---

## 4.2 Time-domain parameters

Four time-domain parameters are computed **directly from the waveform** of each
frame, without any Fourier transform. They are **cheap** and useful for
**voiced/unvoiced** detection, **endpoint** detection, and **pitch**.

### 4.2.1 Short-time energy

$$
\boxed{\;E(m) = \sum_{n=0}^{N-1} x_m[n]^2.\;}
$$

**Interpretation.** The **energy** of frame $m$ — how loud the signal is in that
20–30 ms window.

**Properties.**

- **Voiced speech:** high energy (especially vowels).
- **Unvoiced fricatives:** low-to-moderate energy.
- **Silence:** near zero.
- **Plosive bursts:** brief energy spikes.

**Uses.** **Endpoint detection** (find speech vs. silence), **voiced/unvoiced
discrimination** (with a threshold), **gain control**, **speech enhancement**.

### 4.2.2 Short-time magnitude

$$
\boxed{\;M(m) = \sum_{n=0}^{N-1} |x_m[n]|.\;}
$$

Cheaper than energy (no multiplication). Similar properties. Used in hardware
where multiplication cost was significant (historical) or for simplicity.

### 4.2.3 Short-time zero-crossing rate (ZCR)

$$
\boxed{\;Z(m) = \frac{1}{2}\sum_{n=1}^{N-1} \bigl|\operatorname{sgn}(x_m[n]) - \operatorname{sgn}(x_m[n-1])\bigr|.\;}
$$

Alternatively, the **normalized** form $z(m) = Z(m)/N \in [0, 1]$.

**Interpretation.** The **number of times** the signal crosses zero per frame
(per sample, if normalized). Measures **high-frequency content**.

**Properties.**

- **Voiced speech:** low ZCR (smooth, periodic — crosses zero at the
  fundamental).
- **Unvoiced fricatives:** high ZCR (noise-like, rapid sign changes).
- **Silence:** intermediate (from ambient noise).

**Uses.** **Voiced/unvoiced discrimination**, **endpoint detection** (combined
with energy), **speech/music discrimination**, **noise detection**.

### Worked example 4.2 — energy, magnitude, and ZCR by hand

Frame $x = (1, -1, 2, -2, -1, 1, 1, -1)$.

**Energy.** $1 + 1 + 4 + 4 + 1 + 1 + 1 + 1 = 14$.

**Magnitude.** $1 + 1 + 2 + 2 + 1 + 1 + 1 + 1 = 10$.

**Sign changes.** From $+1$ to $-1$ (1), $-1$ to $+2$ (2), $+2$ to $-2$ (3),
$-2$ to $-1$ (4), $-1$ to $+1$ (5), $+1$ to $+1$ (no), $+1$ to $-1$ (6).

$Z = \tfrac12 \cdot 6 = 3$; normalized $z = 3/8 = 0.375$.

**Classification guess:** energy 14 is moderate, ZCR 0.375 is high → likely
**voiced fricative** or **mixed**. (A pure vowel would have lower ZCR; a pure
silence would have near-zero energy.)

### 4.2.4 Short-time autocorrelation function

$$
\boxed{\;R_m[k] = \sum_{n=0}^{N-1-k} x_m[n]\, x_m[n+k],\quad k = 0, 1, \ldots, K.\;}
$$

**Interpretation.** Measure of the **correlation** between the signal and a
shifted version of itself by $k$ samples. Peaks where the signal is
quasi-periodic.

**Properties.**

- $R_m[0] = E(m)$ — the zero-lag autocorrelation is the energy.
- **Voiced speech:** strong peak at $k = T_0$ (the pitch period in samples).
- **Unvoiced speech:** no strong peak; autocorrelation decays rapidly.

**Uses.** **Pitch estimation** (find the lag of the first major peak, excluding
$k=0$), **voiced/unvoiced discrimination** (peak height vs. energy), **LPC**
(the Yule–Walker equations use $R_X[k]$).

**Computational cost.** Naive: $O(NK)$ per frame. **FFT-based:** $O(N \log N)$
for $R[0]$ through $R[N-1]$, using the Wiener–Khinchin relation:
$R[k] = \text{IDFT}\{|\text{DFT}\{x\}|^2\}$. For $N = 400$, $K = 200$:
naive = $80{,}000$ ops, FFT = $\sim 4{,}000$ ops — a 20× speedup.

### Worked example 4.3 — autocorrelation for pitch detection

Frame at 16 kHz with true pitch period $T_0 = 100$ samples (so $F_0 = 160$ Hz).
Compute $R_m[k]$ for $k = 0, \ldots, 200$.

The autocorrelation will show:
- $R_m[0]$: maximum (energy).
- **Local minima** near $k = 50$ (half period).
- **Local maximum** near $k = 100$ (pitch period) — the **key peak**.
- Another local maximum near $k = 200$ (two pitch periods) — but a weaker
  one due to the finite window.

**Peak picking algorithm.**
1. Compute $R_m[k]$ for $k \in [K_{\min}, K_{\max}]$ — typical $F_0$ range
   limits: 60–400 Hz → $K \in [f_s/400, f_s/60] = [40, 267]$ at 16 kHz.
2. Find the **maximum** in that range.
3. If the peak value exceeds a threshold (e.g., 0.3·$R_m[0]$), declare **voiced** with $T_0 = $ argmax; otherwise declare **unvoiced**.

**Peak at $k = 100$** with $R_m[100]/R_m[0] \approx 0.85$ → confident voiced pitch at $F_0 = 16000/100 = 160$ Hz. ✓

**Pitfall.** **Octave errors** — the autocorrelation has peaks at $k = T_0$
**and** $k = 2T_0, 3T_0, \ldots$, plus at $T_0/2$ if the waveform is
sub-periodic. Picking the "wrong" peak gives pitch doubling or halving.
The **parallel processing** method of §4.3 addresses this.

---

## 4.3 Pitch estimation — parallel processing

### 4.3.1 The parallel processing method

The **parallel processing** pitch detector computes the autocorrelation of the
frame at multiple **downsampled** versions and combines the results:

**Step 1.** Lowpass filter the frame and downsample by a factor $D_1$; compute
the autocorrelation. The peak at $k = T_0/D_1$ gives a coarse pitch estimate.

**Step 2.** Upsample the estimate (multiply by $D_1$) to get the true lag.

**Why?** Downsampling reduces the sampling rate → reduces the number of lags
to compute → cheaper. The pitch estimate at the coarse rate is refined by
combining multiple channels.

**The classic Rabiner–Gold algorithm** uses a **parallel bank of autocorrelators**,
each operating at a different decimation factor and lag range. By choosing
decimations that cover different pitch ranges, the peak at each channel
confirms or rejects the pitch candidate from the others — this **vote-based**
approach resolves octave errors.

### 4.3.2 Other pitch detection methods (brief)

| Method | Domain | Pros | Cons |
|---|---|---|---|
| **Autocorrelation** | Time | Simple, fast | Octave errors |
| **AMDF (average magnitude difference)** | Time | Simpler than ACF; robust to noise | Slight bias |
| **Parallel processing** | Time | Resolves octave errors | More complex |
| **Cepstrum** | Frequency (via IFFT of log-spectrum) | Robust, separates pitch from formants | FFT cost |
| **YIN** | Time | Excellent accuracy; open-source standard | — |
| **HPS (harmonic product spectrum)** | Frequency | Simple, robust | FFT cost |
| **SWIPE** | Frequency | State-of-the-art accuracy | Cost |

The **YIN algorithm** (de Cheveigné & Kawahara, 2002) is the modern standard for
autocorrelation-based pitch detection — it uses the **cumulative mean normalized
difference function** to avoid octave errors and bias.

### Worked example 4.4 — parallel processing by two decimators

True $F_0 = 150$ Hz at $f_s = 16$ kHz → $T_0 = 106.7$ samples.

**Channel 1 — decimate by $D_1 = 4$.** Downsampled rate $f_{s1} = 4$ kHz. The
pitch period in the downsampled signal is $T_0/4 \approx 26.7$ samples. The
autocorrelation peak is at $k_1 \approx 27$.

**Channel 2 — decimate by $D_2 = 8$.** $f_{s2} = 2$ kHz. Period $T_0/8 \approx 13.3$
samples. Peak at $k_2 \approx 13$.

**Consistency check.** $k_1 \cdot D_1 = 27 \cdot 4 = 108$ samples; $k_2 \cdot D_2 = 13 \cdot 8 = 104$
samples. Both estimate the true period ($\sim$107) — the vote confirms the pitch
is $F_0 \approx 16000/106 \approx 151$ Hz. ✓

**Octave-error detection.** If Channel 1 reports $k_1 = 13$ (a half-period
peak), the estimated pitch would be $32000/13 \approx 2462$ Hz — way outside
the plausible range. The parallel processing detects this by disagreeing with
Channel 2; the vote corrects it.

### 4.3.3 AMDF — Average Magnitude Difference Function

A cheaper alternative to autocorrelation:

$$
\boxed{\;D_m[k] = \sum_{n=0}^{N-1-k} \bigl|x_m[n] - x_m[n+k]\bigr|.\;}
$$

**Properties.**

- $D_m[0] = 0$ (trivially — no shift).
- Minimum at $k = T_0$ (quasi-periodic).
- Minimum at $k = 2T_0, 3T_0, \ldots$.
- **Cheaper than ACF** — no multiplications, just subtractions and absolute
  values.

**Detection.** Find the **first local minimum** (after the trivial $k=0$
minimum) — the argmin over $k \in [K_{\min}, K_{\max}]$.

---

## 4.4 Short-time Fourier transform (STFT)

### 4.4.1 Definition

The **short-time Fourier transform** of a frame $x_m[n]$ is

$$
\boxed{\;X_m(e^{j\omega}) = \sum_{n=0}^{N-1} x_m[n]\, w[n]\, e^{-j\omega n}.\;}
$$

In practice, evaluated at the DFT frequencies $\omega_k = 2\pi k/N$:

$$
\boxed{\;X_m[k] = \sum_{n=0}^{N-1} x_m[n]\, w[n]\, e^{-j2\pi k n / N},\quad k = 0, \dots, N-1.\;}
$$

This is the **N-point FFT** of the windowed frame — one FFT per frame.

**Computational cost.** $O(N \log N)$ per frame, $O(N \log N \cdot f_s / H)$
per second. For $N = 400$ at $F_s = 16$ kHz with $H = 200$: $80$ frames/s ×
$400\log_2 400 \approx 3460$ ops/frame = ~$277{,}000$ ops/s. Well within
real-time on any modern CPU.

### 4.4.2 Time-frequency resolution trade-off

The **window length** $N$ determines the trade-off:

- **Wide window** (long $N$): **high frequency resolution** ($\Delta f = f_s/N$), **low time resolution** (windows are long).
- **Narrow window** (short $N$): **high time resolution**, **low frequency resolution**.

This is the **Heisenberg uncertainty principle** for signals:

$$
\boxed{\;\Delta t \cdot \Delta f \ge \frac{1}{4\pi}\;}
$$

where $\Delta t$ and $\Delta f$ are the standard deviations of the window and
its Fourier transform. You cannot have arbitrarily good time AND frequency
resolution simultaneously.

**In speech processing, the standard choices are:**

| Window length | $\Delta t$ | $\Delta f$ (at 16 kHz) | Purpose |
|---|---|---|---|
| 10 ms | 10 ms | 100 Hz | **Wideband** spectrogram — good time resolution, pitch not visible |
| 25 ms | 25 ms | 40 Hz | Standard STFT frame |
| 40 ms | 40 ms | 25 Hz | **Narrowband** spectrogram — pitch harmonics visible |

These are the **two classical spectrogram views**, and every speech tool
offers both:

- **Narrowband spectrogram:** long window (40 ms). Harmonics at $F_0$ are
  resolved; formants appear as bands. Used for **pitch analysis** and **formant
  tracking**.
- **Wideband spectrogram:** short window (10 ms). Harmonics are NOT resolved;
  each pitch pulse is a vertical stripe. Used for **voicing** and **timing**
  analysis.

### Worked example 4.5 — narrowband vs. wideband resolution

Male speech, $F_0 = 100$ Hz, $f_s = 16$ kHz.

**Narrowband** (40 ms = 640 samples): $\Delta f = 16000/640 = 25$ Hz. Harmonic
spacing is 100 Hz — resolved (25 < 100/2). Formants appear as separate bands.
Pitch visible as horizontal striations.

**Wideband** (10 ms = 160 samples): $\Delta f = 16000/160 = 100$ Hz. Harmonic
spacing = 100 Hz — **just at the resolution limit**, so individual harmonics
blur together. Formants still visible. Pitch pulses appear as vertical
striations (one per 10 ms pitch period).

**Rule of thumb.** For pitch analysis, use a window that is **at least 2–3
pitch periods long** ($\ge 2T_0$). For male voices ($T_0 \approx 10$ ms), that's
20–30 ms — the narrowband setting. For female voices ($T_0 \approx 5$ ms),
10–15 ms suffices.

### 4.4.3 The STFT as a filter bank

There is a **dual interpretation** of the STFT:

**Frame interpretation.** For each frame, compute the spectrum (standard DFT of
a windowed segment).

**Filter bank interpretation.** For each frequency $\omega_k$, the STFT
$\{X_m[k]\}_{m}$ is a **narrowband signal** at frequency $\omega_k$ — the output
of a **bandpass filter** centered at $\omega_k$ with bandwidth $\sim 1/N$.
This is exactly a **filter bank** with $N$ channels.

The two views are equivalent — a math theorem (the **filter-bank summation
theorem**) states that reconstructing $x$ from $X_m[k]$ (summing over $k$) is
identical to reconstructing from $x_m[n]$ (overlapping frames).

**Practical implication.** Speech can be modelled as a sum of **narrowband
envelopes** — each $\{X_m[k]\}_m$ is a slow-varying (envelope-like) signal at
the frame rate. This is why the **channel vocoder** works (§4.9): it transmits
only these slowly-varying envelopes plus a pitch / voicing parameter.

---

## 4.5 Filter bank analysis

### 4.5.1 Uniform filter bank

A **uniform filter bank** splits the spectrum into equal-bandwidth channels.
For $M$ channels, each covers $\pi/M$ of the frequency axis.

**Implementation.** Two equivalent options:

1. **FFT-based:** compute the N-point FFT of each frame, then **sum energy** in each channel's frequency bins. Cheapest.
2. **Filter-based:** apply $M$ bandpass FIR filters in the time domain, decimate each to rate $\pi/M$. Standard for hardware.

**Energy in channel $c$ of frame $m$:**

$$
\boxed{\;E_m(c) = \sum_{k \in \text{channel } c} |X_m[k]|^2.\;}
$$

The sequence $\{E_m(c)\}_m$ is the **channel envelope** — used directly in
vocoders, speech recognition, and speech enhancement.

### 4.5.2 Mel filter bank

Speech perception is **logarithmic** in frequency — the ear is much more
sensitive to low frequencies than high. The **mel scale** approximates this:

$$
\boxed{\;f_{\text{mel}} = 2595 \log_{10}\!\left(1 + \frac{f_{\text{Hz}}}{700}\right).\;}
$$

A **mel filter bank** places **more filters** at low frequencies and **fewer**
at high frequencies. This is the front end of every modern speech-recognition
system (MFCC features).

**Typical parameters.** 20–40 triangular filters covering 0–8000 Hz.
Lowest filter covers 0–100 Hz; highest covers 7–8 kHz.

**MFCC pipeline.**

1. Frame + window (25 ms Hamming, hop 10 ms).
2. FFT of each frame.
3. **Mel filter bank** — project onto 20–40 triangular filters.
4. **Log** — $L_m(c) = \log E_m(c)$.
5. **Discrete cosine transform** — take first 12–13 DCT coefficients. These are the **MFCCs**.
6. (Optional) **Delta / delta-delta** — first and second time derivatives, to capture dynamics.

The **DCT** decorrelates the log-mel energies and compacts them into a few
coefficients — the same energy-compaction principle as image DCT (ECEG-6423
Chapter 2), applied to the mel-spectrum.

### Worked example 4.6 — Mel filter bank spacing

At $f_s = 16$ kHz, design a 24-filter mel bank covering 0–8000 Hz.

**Mel limits.** $f_{\text{mel}}(0) = 0$; $f_{\text{mel}}(8000) = 2595 \log_{10}(1 + 8000/700) = 2595 \log_{10}(12.43) \approx 2840$ mels.

**Mel step.** $2840/25 \approx 113.6$ mels per filter boundary.

**Filter centers in Hz.**

| Filter | Mel | $f$ (Hz) |
|:---:|:---:|:---:|
| 1 | 113 | ≈ 76 |
| 2 | 227 | ≈ 168 |
| 5 | 568 | ≈ 500 |
| 10 | 1136 | ≈ 1200 |
| 15 | 1704 | ≈ 2400 |
| 20 | 2272 | ≈ 4200 |
| 24 | 2726 | ≈ 7000 |

**Observations.** The first 5 filters span 0–500 Hz (a 500 Hz range) — resolution
~100 Hz per filter. The last 5 span 4200–8000 Hz (a 3800 Hz range) — resolution
~760 Hz per filter. **The mel scale gives 7–8× more resolution to low
frequencies** — exactly matching the ear's sensitivity.

The number of **filters per octave** is roughly constant in the mel domain
(~6–10), which is what the perception literature says humans hear.

---

## 4.6 Spectrographic analysis

### 4.6.1 Definition

The **spectrogram** is a 2-D visualisation of the STFT:

- **X-axis:** time (frames).
- **Y-axis:** frequency (FFT bin or Hz).
- **Color / intensity:** magnitude (usually in dB), from light (low energy) to
  dark (high energy).

The **short-time magnitude** is typically plotted:

$$
S[m, k] = 20 \log_{10}|X_m[k]| + \text{offset (dB)}.
$$

**Two spectrogram types** (from §4.4.2):

| Type | Window | Shows | Use |
|---|---|---|---|
| **Narrowband** | 40 ms | Harmonics, formants as horizontal bands | Pitch analysis, formant tracking |
| **Wideband** | 10 ms | Pitch pulses as vertical striations, formant transitions | Timing, voicing, phonetics |

Speech processing tools (Praat, Wavesurfer, MATLAB's `spectrogram`, Python's
`matplotlib.specgram`) support both. **The choice depends on what you want to
see.**

### 4.6.2 Reading a spectrogram

| Visible feature | Interpretation |
|---|---|
| **Horizontal dark bands** | Formants — vowel identity |
| **Closely spaced parallel stripes** (narrowband) | Pitch harmonics — voiced |
| **Vertical stripes** (wideband) | Pitch pulses |
| **Diffuse high-frequency energy** | Fricative (/s/, /ʃ/, /f/) |
| **Short vertical spike** | Plosive burst |
| **Dark low-frequency band, no high harmonics** | Nasal |
| **Formant transitions** | Coarticulation between phonemes |
| **Silent regions** | Silence or stop closure |

### Worked example 4.7 — reading a spectrogram by hand

Consider the spectrogram of the word "banana" /bəˈnænə/, wideband (10 ms):

| Time | Pattern | Phone |
|---|---|---|
| 0–50 ms | Blank (silence for closure) | /b/ closure |
| 50–60 ms | Short vertical burst | /b/ release |
| 60–150 ms | Two horizontal bands near 500, 1500 Hz | /ə/ (schwa) |
| 150–200 ms | Bands near 500, 1500 Hz, muted | /n/ (nasal — damped) |
| 200–300 ms | Bands near 660, 1720 Hz (different) | /æ/ |
| 300–350 ms | Bands again near 500, 1500 Hz | /n/ |
| 350–450 ms | Bands near 500, 1500 Hz | /ə/ |

**Formant tracks:** $F_1$ moves 500→660→500; $F_2$ moves 1500→1720→1500. The
**$F_1$–$F_2$ trajectory** is the direct fingerprint of the vowel sequence.

**Voicing:** vertical striations throughout the voiced segments show pitch pulses
at $F_0 \approx 100$ Hz (visible as ~10 stripes per 100 ms).

**Nasals:** visible as a **reduction** in the higher formants and an extra
**low-frequency resonance** (the nasal murmur) around 300 Hz — the /m/, /n/,
/ŋ/ signature.

**Training note.** Phoneticians train for months to read spectrograms reliably.
For most signal processing applications, you don't need to read them — you just
need to know what features exist so you can design algorithms to extract them.

---

## 4.7 Formant extraction

Formants are the vocal-tract resonances. Two practical approaches:

### 4.7.1 Spectral peak picking

The naive method: find local maxima in the smoothed magnitude spectrum.

1. Compute the **spectral envelope** (from the LPC spectrum or by lowpass-filtering the DFT magnitude).
2. Find local maxima above a threshold.
3. Sort by frequency; keep the first $K$ (typically 3–5).

**Challenges.**

- **Pitch harmonics** create spurious peaks — smoothing is essential.
- **Closely spaced formants** (e.g., $F_2$ and $F_3$ close) merge into one peak.
- **Narrow peaks** may be missed by smoothing.

### 4.7.2 LPC root finding

The LPC polynomial

$$
A(z) = 1 + \sum_{k=1}^{p} a_k z^{-k}
$$

has $p$ roots, typically in complex-conjugate pairs. Each **pair** corresponds to
one formant. The formant frequency is the root's **angle**:

$$
\boxed{\;F_i = \frac{f_s}{2\pi}\,\angle z_i,\quad \text{BW}_i = -\frac{f_s}{\pi}\ln|z_i|.\;}
$$

**Advantages.**

- **Each pole is a formant** — no peak picking required.
- **Bandwidth** comes out naturally (from the pole radius).
- **$p = 2K$ gives $K$ formants** — pick $p = 10$ for 5 formants, $p = 12$–$16$ for 6–8.

**Disadvantages.**

- Requires estimating LPC coefficients (solve Yule–Walker).
- Root-finding is $O(p^3)$ — moderate for $p \le 20$.
- **Unstable or spurious roots** appear at low energy or during nasal sounds.

**Practical.** Use LPC roots for the main formants; supplement with spectral
peaks when roots are ambiguous.

### Worked example 4.8 — formant from LPC pole

LPC(2) with $a_1 = -1.5$, $a_2 = 0.7$.

**Roots.** $z^2 - 1.5z + 0.7 = 0 \Rightarrow z = 0.75 \pm j 0.37$ (from Chapter 3
Example 3.5).

**Formant frequency.** $|z| \approx 0.837$, $\angle z \approx 0.457$ rad.

$$
F = \frac{f_s}{2\pi} \cdot 0.457 \approx \frac{16000}{6.283} \cdot 0.457 \approx 1164\ \text{Hz}.
$$

Wait — Chapter 3's example used $f_s = 8000$: $F \approx 582$ Hz. With
$f_s = 16$ kHz: $F \approx 1164$ Hz. So the formant scales with sampling
rate — always specify $f_s$.

**Bandwidth.** $\text{BW} = -(f_s/\pi) \ln(0.837) \approx (16000/\pi) \cdot 0.178 \approx 906$ Hz
at 16 kHz (≈453 Hz at 8 kHz).

**Quality factor.** $Q = F/\text{BW} \approx 1164/906 \approx 1.28$ — a broad,
low-Q formant. Higher-$Q$ formants (narrower bandwidths) come from higher-$p$
models with poles closer to the unit circle.

**For a real vowel** ($p = 10$), you'd get 5 formants with decreasing bandwidth
and increasing frequency. Typical: $F_1 \approx 500$ Hz (BW 60), $F_2 \approx 1500$ Hz (BW 80),
$F_3 \approx 2500$ Hz (BW 100), $F_4 \approx 3500$ Hz (BW 150), $F_5 \approx 4500$ Hz (BW 200).

---

## 4.8 Pitch extraction — frequency domain

### 4.8.1 Harmonic product spectrum (HPS)

For a **voiced** frame, the spectrum has peaks at $F_0, 2F_0, 3F_0, \ldots$. The
**HPS** exploits this by multiplying downsampled versions of the magnitude
spectrum:

$$
\boxed{\;P[k] = \prod_{r=1}^{R} |X_r[k]|,\quad X_r[k] = |X[r \cdot k]|.\;}
$$

where $X_r$ is the spectrum **decimated by $r$**. Where the true pitch
$F_0$ lies, **all harmonics** align → the product has a peak. Where a
non-harmonic frequency lies, only one or two harmonics contribute → the product
is small.

**Algorithm.**

1. Compute the FFT of the frame.
2. Compute the magnitude spectrum $|X|$.
3. Form the HPS: $P = |X| \cdot |X|_{\downarrow 2} \cdot |X|_{\downarrow 3} \cdot \ldots$.
4. Find the peak of $P$ within the plausible pitch range (60–400 Hz).
5. The peak location gives $F_0$.

**Advantages.** Simple, cheap (one FFT + $R$ decimations + $R$ multiplications).

**Disadvantages.** Fails for weak fundamental (some speakers); limited by FFT
resolution ($\Delta f = f_s/N$). For $N = 400$ at 16 kHz, resolution is 40 Hz —
too coarse for 60 Hz pitch. Use **longer frames** (40+ ms) or **interpolation** of
the peak.

### 4.8.2 Cepstrum-based pitch

See §4.10 — the cepstrum is the standard frequency-domain pitch method.

### Worked example 4.9 — HPS with 3 factors

Frame has $F_0 = 200$ Hz, $f_s = 16$ kHz, $N = 512$ (FFT resolution 31.25 Hz).
Spectrum peaks at bins $k = 200/31.25 = 6.4$, $12.8$, $19.2$, $25.6$, $\ldots$
(the true pitch is at bin 6.4 — between bins, so the peak "leaks" to bins 6 and 7).

**HPS with $R = 3$.**

At bin $k = 6$: $|X[6]| \cdot |X[12]| \cdot |X[18]| \approx$ all moderate (all near
the true peak, close to bins 6.4, 12.8, 19.2). $P[6]$ is high.

At bin $k = 12$: $|X[12]| \cdot |X[24]| \cdot |X[36]|$ — first factor matches
$2F_0$, but the others don't align. Product lower.

At bin $k = 6.5$ (interpolated): all three factors peak → maximum of $P$.

**Interpolation matters.** With FFT resolution 31.25 Hz, the peak at 200 Hz is
at bin 6.4 — picking bin 6 or 7 gives 187 or 219 Hz, both 6–10% off. Parabolic
interpolation of the HPS peak recovers 200 Hz to ~1 Hz accuracy.

---

## 4.9 Analysis–synthesis systems — the channel vocoder

### 4.9.1 The channel vocoder

The **channel vocoder** (Dudley, 1939) is the first speech analysis–synthesis
system. It splits the speech into $C$ frequency channels, extracts the **energy
envelope** of each, and reconstructs by modulating a noise or pulse-train source
with the envelopes.

**Analysis.**

1. Filter into $C$ bandpass channels (typically 16–20, uniform or log-spaced).
2. Rectify and lowpass each channel → **envelope** $E_c(t)$.
3. Estimate **pitch** $F_0$ and **voicing** (voiced/unvoiced).

**Synthesis.**

1. Generate a source: **pulse train** at $F_0$ (voiced) or **white noise** (unvoiced).
2. Modulate each channel's envelope $E_c$ by the source filtered through the
   corresponding bandpass filter.
3. Sum the channels.

**Output.** A speech signal whose **envelopes match** the original but whose
**fine structure** is regenerated from the source model — hence "vocoder"
(**vo**ice **cod**er).

**Bit rate.** With 20 channels × 50 frames/s × ~5 bits/channel = 5 kbps, plus
$F_0$ (6 bits) and voicing (1 bit) — ~5–6 kbps. Older analog vocoders used
~2.4 kbps.

**Sound quality.** Recognisable but **robotic** — the phase information is lost
(envelopes only). This is the famous "vocoder" sound, used for effect in
electronic music (e.g., Daft Punk, Kraftwerk).

### 4.9.2 Modern analysis–synthesis

Modern speech codecs improve on the channel vocoder by:

- **LPC-based analysis** (Chapter 3) — the all-pole model captures formants
  more compactly than 20 filter-bank envelopes.
- **Residual coding** — transmit the prediction error (not just the
  source model), preserving waveform fine structure.
- **CELP (Code-Excited Linear Prediction)** — model the residual as a
  codebook lookup, giving near-toll quality at 4.8–16 kbps.

This is the family that includes G.729, AMR, and Opus.

### Worked example 4.10 — channel vocoder by hand

Speech frame 25 ms, filter bank of 4 channels: 0–500, 500–1500, 1500–3000,
3000–8000 Hz. Observed channel energies (in dB): $(-20, -15, -10, -30)$.
Voicing: voiced, $F_0 = 150$ Hz.

**Reconstruction.**

1. Source: pulse train at 150 Hz.
2. Bandpass into the 4 channels.
3. Multiply each by the envelope $10^{E_c/20}$ → factors $0.1, 0.178, 0.316, 0.0316$.
4. Sum.

The output sounds like the original vowel (same 4 band energies, same pitch)
but with **no phase detail** — the harmonic phases are determined by the source
pulse train, not the original. **5–10× compression** vs. raw PCM, at the cost
of "robotic" timbre.

**Enhancement.** Use 16–20 channels (instead of 4) and finer frame rate (100/s);
the robotic effect reduces but doesn't disappear.

---

## 4.10 Homomorphic analysis and the cepstrum

### 4.10.1 The homomorphic idea

In Chapter 2, homomorphic filtering was used to separate **multiplicative** image
components (illumination × reflectance) via **log**. In speech, the source-filter
model gives a similar multiplicative structure in the spectrum:

$$
\boxed{\;|X(f)| = |G(f)| \cdot |H(f)|\;}
$$

where $|G|$ is the **source** spectrum (glottal pulse train, $F_0$-periodic) and
$|H|$ is the **vocal-tract** spectrum (slowly-varying formant envelope).

**Taking the log** turns multiplication into addition:

$$
\log|X(f)| = \log|G(f)| + \log|H(f)|.
$$

Now the source and filter are **additive** — separable by linear filtering.

### 4.10.2 The real cepstrum

The **real cepstrum** is the **inverse Fourier transform of the log-magnitude
spectrum**:

$$
\boxed{\;c[n] = \mathcal F^{-1}\{\log|X(f)|\}.\;}
$$

In practice, for a framed signal:

$$
c_m[n] = \text{IFFT}\{\log|\text{FFT}\{x_m\}|\}.
$$

**Naming.** "Cepstrum" is a play on "spectrum" (invented by Bogert et al., 1963);
the independent variable is the **quefrency** (a play on "frequency"), in
**seconds or samples**, and the "low-quefrency" / "high-quefrency" regions
correspond to the **envelope** and **fine structure** respectively.

### 4.10.3 Separating source from filter

The key insight: **the source and the filter live at different quefrencies.**

- **Source** (pitch pulses) → periodic structure with period $T_0$ → in the
  cepstrum, a **peak at quefrency $n = T_0$** (or $f_s T_0$ samples). This is
  the **pitch peak**.
- **Filter** (formants) → slow variation in log-spectrum → **low quefrency**
  components ($n$ near 0).

**Liftering** (filtering in the cepstral domain):

- **Low-pass lifter** — keep $c[n]$ for $|n| < n_0$ → recovers the **spectral
  envelope** (formants). Inverse-transform → **smooth spectral envelope**.
- **High-pass lifter** — keep $c[n]$ for $|n| > n_0$ → recovers the **fine
  structure** (harmonics).

**Quefrency cutoff.** Choose $n_0$ between the formant time scale (~1–2 ms) and
the pitch period (~5–10 ms). Typical: $n_0 \approx 3$–$4$ ms, i.e., $n_0 \approx 50$
samples at 16 kHz.

### 4.10.4 Cepstral pitch estimation

**Algorithm.**

1. Frame, window, compute FFT.
2. Take log-magnitude: $L = \log|X|$.
3. IFFT: $c = \text{IFFT}\{L\}$ — the real cepstrum.
4. Find the **peak** of $c[n]$ in the **pitch quefrency range**:
   $n \in [f_s/400, f_s/60] = [40, 267]$ samples at 16 kHz (pitch 60–400 Hz).
5. Peak position $n_{\text{peak}}$ gives $T_0 = n_{\text{peak}}$ samples, $F_0 = f_s/n_{\text{peak}}$.
6. If the peak exceeds a threshold relative to the low-quefrency energy → voiced.

**Advantages over autocorrelation.**

- **Robust to formants** — formants live at low quefrency, don't interfere with
  pitch peak.
- **Deals with weak fundamentals** — even if $|X(f_0)|$ is small, the log-spectrum
  still has the pitch-period structure at $n = T_0$.
- **Cheap** — one FFT + one IFFT (or two FFTs total).

**Disadvantages.**

- **Phase ignored** — log-magnitude spectrum loses phase; slight accuracy
  penalty.
- **Requires long frames** (40+ ms) so the pitch period fits comfortably at
  moderate quefrency.
- **Sensitive to DC** — subtract the mean before FFT to avoid huge low-quefrency
  peak.

### Worked example 4.11 — cepstrum peak for $F_0 = 200$ Hz

$f_s = 16$ kHz, $T_0 = 80$ samples ($F_0 = 200$ Hz).

**Log-magnitude spectrum.** Harmonic peaks every $f_s/T_0 = 200$ Hz.

**Real cepstrum.** IFFT of this log-spectrum. Structure:
- **Low quefrency** ($n < 50$): the smooth formant envelope.
- **Peak at $n = 80$**: the pitch period — a strong peak (the IFFT of a harmonic
  comb is a comb at the pitch period).
- **Peak at $n = 160$**: two pitch periods (weaker, from second-harmonic
  interaction).

**Detection.** Peak at $n = 80$ → $F_0 = 16000/80 = 200$ Hz. ✓

**Threshold.** Peak height relative to the surrounding quefrency: if
$\max_{n \in [40, 267]} c[n] > \alpha \cdot c[0]$ for some $\alpha \approx 0.1$,
declare voiced. The exact threshold depends on the recording SNR and can be
tuned.

### 4.10.5 Mel-frequency cepstral coefficients (MFCC)

If we replace the **DFT frequencies** in the cepstrum computation with **mel
frequencies** (via the mel filter bank from §4.5.2), we get **MFCCs**:

**MFCC pipeline.**

1. **Frame** (25 ms) and **window** (Hamming).
2. **FFT** → magnitude spectrum $|X[k]|$.
3. **Mel filter bank** → energies $E_c$ in $C$ mel-spaced bands.
4. **Log** → $L_c = \log E_c$.
5. **DCT** → $c_1, c_2, \dots, c_J$ (first $J \approx 12$–13 coefficients).

The **DCT** replaces the IFFT — it decorrelates the log-mel energies and packs
them into the first few coefficients. **MFCCs** are the **standard feature
vector** for speech recognition, speaker recognition, and many speech
processing tasks.

**Why MFCC not raw cepstrum?** The **mel warping** matches human perception — it
puts the spectral resolution where it matters (low frequencies), which improves
recognition accuracy by 5–10% over uniform cepstrum.

---

## 4.11 Linear prediction of speech (LPC)

### 4.11.1 The LPC model

**LPC** models each speech sample as a **linear combination** of the previous $p$
samples plus a **prediction error**:

$$
\boxed{\;x[n] = -\sum_{k=1}^{p} a_k\, x[n-k] + e[n].\;}
$$

**Reinterpreted:** $e[n] = x[n] + \sum_{k=1}^{p} a_k x[n-k] = \sum_{k=0}^{p} a_k x[n-k]$ (with $a_0 = 1$).

**Model.** The signal $x[n]$ is the output of an **all-pole filter** $H(z) = 1/A(z)$
driven by the **prediction error** $e[n]$:

$$
\boxed{\;H(z) = \frac{1}{1 + \sum_{k=1}^{p} a_k z^{-k}}.\;}
$$

**Physical interpretation.** $H(z)$ is the **vocal-tract filter**; $e[n]$ is the
**glottal source**. This is the **source-filter model** realized as an AR(p) process.

### 4.11.2 Why LPC of speech?

**Three properties of the vocal tract make LPC a natural fit:**

1. **Vocal tract ≈ all-pole.** Each formant is a 2nd-order resonance (pole pair). $K$ formants → $2K$ poles → AR(2K). For $K = 5$ (F1–F5), $p = 10$ is enough.
2. **Glottal source ≈ white noise** (unvoiced) or a **pulse train** (voiced). Either way, it's broadband and not captured by the filter — makes $e[n]$ approximately white.
3. **Formants are acoustically important** — they determine vowel identity. LPC coefficients concentrate the formant information in a compact vector.

### 4.11.3 Derivation — Yule–Walker equations

Minimise the prediction error **energy** over a frame:

$$
\boxed{\;\mathcal E = \sum_{n} e[n]^2 = \sum_n \Bigl(x[n] + \sum_{k=1}^{p} a_k x[n-k]\Bigr)^2.\;}
$$

**Take partials and set to zero.**

$$
\frac{\partial \mathcal E}{\partial a_k} = 2\sum_n \Bigl(x[n] + \sum_{l=1}^{p} a_l x[n-l]\Bigr) x[n-k] = 0,\quad k = 1, \ldots, p.
$$

Rearranged:

$$
\boxed{\;\sum_{l=1}^{p} a_l \underbrace{\Bigl(\sum_n x[n-l] x[n-k]\Bigr)}_{= R[k-l]} = -\underbrace{\sum_n x[n] x[n-k]}_{= R[k]}.\;}
$$

**The Yule–Walker equations:**

$$
\boxed{\;\mathbf R\, \mathbf a = -\mathbf r,\quad
R_{ij} = R[|i-j|],\quad r_k = R[k].\;}
$$

with $\mathbf a = (a_1, \ldots, a_p)^\top$ and $\mathbf R$ the $p \times p$
**Toeplitz** autocorrelation matrix.

**Prediction error variance.**

$$
\boxed{\;\sigma_e^2 = R[0] + \sum_{k=1}^{p} a_k R[k].\;}
$$

### 4.11.4 Solving the Yule–Walker equations — Levinson–Durbin

**Direct inversion:** $O(p^3)$.

**Levinson–Durbin** (from ECEG-6303 Chapter 2 §2.10): exploits the Toeplitz
structure to solve in $O(p^2)$.

**Recursion** (from order $m-1$ to order $m$):

$$
\boxed{\;k_m = -\frac{R[m] + \sum_{i=1}^{m-1} a_i^{(m-1)} R[m-i]}{\sigma_{m-1}^2}\;}
$$

$$
\boxed{\;a_i^{(m)} = a_i^{(m-1)} + k_m a_{m-i}^{(m-1)},\quad i = 1, \ldots, m-1;\quad a_m^{(m)} = k_m.\;}
$$

$$
\boxed{\;\sigma_m^2 = (1 - k_m^2)\sigma_{m-1}^2.\;}
$$

Initialise $\sigma_0^2 = R[0]$, $a^{(0)} = \emptyset$.

**Reflection coefficients** $k_m$: partial autocorrelation at lag $m$; satisfy
$|k_m| < 1$ for a stable AR filter.

**Why "prediction error" and "reflection coefficients" coincide.** The physical
interpretation: $k_m$ is the fraction of energy **removed** by adding lag $m$
to the predictor. Magnitude near 1 = lag $m$ is very informative; near 0 = lag $m$
adds little. This is precisely what "reflection coefficient" means in the
acoustic-tube model (Chapter 3 §3.5.3).

### Worked example 4.12 — LPC(2) by Yule–Walker

Frame autocorrelation: $R[0] = 6$, $R[1] = 4$, $R[2] = 1$.

**Yule–Walker system** (from Chapter 2 Example 2.7):

$$
\begin{bmatrix} R[0] & R[1] \\ R[1] & R[0] \end{bmatrix}\begin{bmatrix} a_1 \\ a_2 \end{bmatrix} = -\begin{bmatrix} R[1] \\ R[2] \end{bmatrix}
$$

$$
\begin{bmatrix} 6 & 4 \\ 4 & 6 \end{bmatrix}\begin{bmatrix} a_1 \\ a_2 \end{bmatrix} = \begin{bmatrix} -4 \\ -1 \end{bmatrix}.
$$

Solution: $a_1 = -1$, $a_2 = 0.5$.

**Prediction error variance.**

$$
\sigma_e^2 = R[0] + a_1 R[1] + a_2 R[2] = 6 + (-1)(4) + (0.5)(1) = 2.5.
$$

**LPC filter.**

$$
H(z) = \frac{1}{1 - z^{-1} + 0.5 z^{-2}}.
$$

**Roots.** $1 - z^{-1} + 0.5z^{-2} = 0 \Rightarrow z^2 - z + 0.5 = 0 \Rightarrow z = 0.5 \pm j0.5$.

$|z| = 1/\sqrt 2 \approx 0.707$, $\angle z = \pi/4$. **Formant at $F = (f_s/2\pi)(\pi/4) = f_s/8$.**
At $f_s = 8000$: $F = 1000$ Hz. **Bandwidth:** $-\ln|z| \cdot f_s/\pi = -\ln(1/\sqrt 2) \cdot 8000/\pi \approx 1760$ Hz. Very wide
(a heavily damped resonance).

**Interpretation.** The two LPC coefficients capture **one** broad formant
because $p=2$ is enough for only one pole pair. Real speech requires $p = 10$
for 5 formants. This example just illustrates the mechanics.

---

## 4.12 LPC-based formant and pitch estimation

### 4.12.1 Formants from LPC roots

For an LPC(p) model with roots $z_i = r_i e^{j\theta_i}$:

$$
\boxed{\;F_i = \frac{f_s}{2\pi}\,\theta_i,\quad \text{BW}_i = -\frac{f_s}{\pi}\,\ln r_i.\;}
$$

**Extracting formants.**

1. Solve Yule–Walker for $A(z) = 1 + \sum a_k z^{-k}$.
2. Find all $p$ roots $z_i$.
3. **Keep complex-conjugate pairs** (real roots are not formants — they're DC or
   Nyquist poles).
4. Compute $F_i$ and $\text{BW}_i$ for each pair.
5. **Sort by frequency**, keep the first 3–5 (the lowest formants).
6. **Filter by bandwidth** — reject pairs with BW > ~500 Hz (typically spurious).

**Typical results** for a vowel with $p = 10$:
- 5 complex-conjugate pairs → 5 formants.
- Some may be **too wide** (spurious) or **too high** (above Nyquist/2).

**The pole-pair formants** are the LPC coefficients' physical meaning, and this
is the basis of the **formant vocoder** and **speech synthesis** systems.

### 4.12.2 Pitch from LPC residual

The prediction **error** $e[n] = x[n] + \sum a_k x[n-k]$ — the LPC **residual** —
is the **glottal source**. For voiced speech, $e[n]$ contains the pitch pulses.
Apply pitch detection (autocorrelation, AMDF, or cepstrum) to $e[n]$ to estimate
$F_0$.

**Advantages.**

- **Whiter residual** — the LPC filter has removed the formants, so pitch peaks
  in $e[n]$ are clearer than in $x[n]$.
- **Especially useful** for voiced speech with strong formants — the ACF of
  $x[n]$ has peaks at $T_0$, but also peaks from formant periodicity; the
  residual's ACF has only the pitch peak.

**Practical.** For high-quality pitch tracking, use the **SIFT algorithm**
(Markel, 1972) — LPC residual + autocorrelation peak picking. It's the
industry standard for pitch detection in LPC-based codecs.

### Worked example 4.13 — formant extraction from LPC(10)

LPC(10) of a vowel frame at 16 kHz. The 10 roots (5 conjugate pairs):

| Pair | Root magnitude | Root angle | $F$ (Hz) | BW (Hz) |
|:---:|:---:|:---:|:---:|:---:|
| 1 | 0.98 | 0.20 | 510 | 100 |
| 2 | 0.96 | 0.59 | 1500 | 200 |
| 3 | 0.94 | 0.98 | 2500 | 305 |
| 4 | 0.91 | 1.38 | 3510 | 460 |
| 5 | 0.87 | 1.77 | 4500 | 660 |

**Formants:** $F_1 \approx 510$ Hz, $F_2 \approx 1500$ Hz, $F_3 \approx 2500$ Hz.

**Interpretation.** $(F_1, F_2) = (510, 1500)$ is near the **schwa /ə/** region
of the vowel chart. The bandwidths increase with frequency (typical — high
formants are more damped).

**Filtering.** Pair 5 has BW = 660 Hz, above the typical threshold (500 Hz) —
it's likely **spurious** (a spurious root from the finite-window analysis). Keep
the first 3–4 formants; ignore the rest.

---

## 4.13 Summary of methods — choosing your tool

| Task | Time-domain | Frequency-domain | LPC |
|---|---|---|---|
| **Voiced/unvoiced** | Energy + ZCR (cheap) | Spectral flatness | LPC residual + autocorrelation |
| **Pitch** | Autocorrelation, AMDF, parallel processing | Cepstrum, HPS | Residual + ACF (SIFT) |
| **Formants** | — | Spectral peaks | LPC roots |
| **Envelope** | — | Smoothed FFT magnitude | LPC spectrum $1/|A(e^{j\omega})|^2$ |
| **Channel energies** | — | Filter bank / mel | — |
| **Feature vectors** | — | MFCC | LPC, LPCC, PLP |
| **Coding** | DPCM | Channel vocoder | CELP, G.729 |

**Rules of thumb.**

- **Cheapest:** time-domain (energy, ZCR, AMDF) — for voiced/unvoiced, endpoint.
- **Most accurate pitch:** cepstrum (frequency) or SIFT (LPC residual).
- **Best formants:** LPC roots or smoothed spectral peaks.
- **Standard recognition features:** MFCC (or PLP).
- **Standard coding:** LPC + residual (CELP family).

---

## 4.14 Computational implementation

### Python — full pipeline

```python
import numpy as np
from scipy.signal import get_window
from scipy.fft import dct

def frame_signal(x, frame_length, hop_length):
    """Split signal into overlapping frames."""
    n_frames = 1 + (len(x) - frame_length) // hop_length
    frames = np.zeros((n_frames, frame_length))
    for m in range(n_frames):
        frames[m] = x[m*hop_length : m*hop_length + frame_length]
    return frames

def short_time_energy(frames):
    return np.sum(frames**2, axis=1)

def short_time_zcr(frames):
    return 0.5 * np.mean(np.abs(np.diff(np.sign(frames), axis=1)), axis=1)

def short_time_autocorr(frame, max_lag):
    """Biased (for stability) autocorrelation up to max_lag."""
    n = len(frame)
    # FFT-based for efficiency
    X = np.fft.rfft(frame, n=2*n)
    r = np.fft.irfft(np.abs(X)**2, n=2*n)[:max_lag + 1]
    return r / n

def pitch_from_acf(frame, fs, f_min=60, f_max=400):
    """Autocorrelation-based pitch detection."""
    K_min = int(fs / f_max)
    K_max = int(fs / f_min)
    r = short_time_autocorr(frame, K_max)
    # Peak in [K_min, K_max]
    k_peak = K_min + np.argmax(r[K_min:K_max + 1])
    # Voicing decision
    if r[k_peak] / r[0] > 0.3:
        return fs / k_peak, True
    return 0.0, False

def stft(x, frame_length=400, hop_length=200, window='hamming'):
    """Short-time Fourier transform (one-sided)."""
    w = get_window(window, frame_length)
    frames = frame_signal(x, frame_length, hop_length) * w
    return np.fft.rfft(frames, axis=1)

def real_cepstrum(frame, n_fft=1024):
    """Real cepstrum of a single frame."""
    X = np.fft.rfft(frame, n=n_fft)
    log_mag = np.log(np.abs(X) + 1e-10)
    # Compute real cepstrum via IFFT of the log spectrum
    # Full spectrum: mirror
    log_mag_full = np.concatenate([log_mag, log_mag[-2:0:-1]])
    c = np.fft.irfft(log_mag_full, n=n_fft)
    return c[:n_fft//2]

def pitch_from_cepstrum(frame, fs, f_min=60, f_max=400):
    """Cepstrum-based pitch detection."""
    c = real_cepstrum(frame)
    n_min = int(fs / f_max)
    n_max = int(fs / f_min)
    n_peak = n_min + np.argmax(c[n_min:n_max + 1])
    if c[n_peak] > 0.1 * c[0]:
        return fs / n_peak, True
    return 0.0, False

def lpc(frame, p, fs):
    """LPC coefficients via Levinson-Durbin."""
    # Autocorrelation
    N = len(frame)
    r = np.correlate(frame, frame, mode='full')[N-1:N+p]
    # Levinson-Durbin
    a = np.zeros(p + 1)
    a[0] = 1.0
    E = r[0]
    for m in range(1, p + 1):
        k = -(r[m] + np.dot(a[1:m], r[m-1:0:-1])) / E
        # Update
        a_new = a.copy()
        a_new[m] = k
        for i in range(1, m):
            a_new[i] = a[i] + k * a[m - i]
        a = a_new
        E *= (1 - k**2)
    return a, E

def formants_from_lpc(a, fs):
    """Formants from LPC polynomial roots."""
    roots = np.roots(a)
    roots = roots[np.imag(roots) > 0]  # upper half plane
    formants = []
    for z in roots:
        f = fs / (2 * np.pi) * np.angle(z)
        bw = -fs / np.pi * np.log(np.abs(z))
        if 50 < f < fs / 2 and bw < 500:
            formants.append((f, bw))
    formants.sort()
    return formants

# ---- demo -------------------------------------------------------------
if __name__ == "__main__":
    fs = 16000
    # Simulate a synthetic voiced signal: pitch 150 Hz, formants 500, 1500, 2500
    t = np.arange(int(0.5 * fs)) / fs
    x = sum(np.sin(2*np.pi*f*t) for f in [150, 300, 450, 600])
    # Add some noise
    x += 0.1 * np.random.randn(len(x))

    # Frame-based features
    frames = frame_signal(x, frame_length=400, hop_length=200)
    E = short_time_energy(frames)
    Z = short_time_zcr(frames)
    print(f"Frames: {len(frames)}, mean E = {E.mean():.2f}, mean ZCR = {Z.mean():.4f}")

    # Pitch from ACF and cepstrum on one frame
    frame = frames[10]
    f_acf, v_acf = pitch_from_acf(frame, fs)
    f_cep, v_cep = pitch_from_cepstrum(frame, fs)
    print(f"Pitch: ACF = {f_acf:.1f} Hz (voiced={v_acf}), "
          f"cepstrum = {f_cep:.1f} Hz (voiced={v_cep})")

    # LPC and formants
    a, E_err = lpc(frame, p=10, fs=fs)
    formants = formants_from_lpc(a, fs)
    print(f"Prediction error variance: {E_err:.4f}")
    print(f"Formants: {[(round(f, 1), round(bw, 1)) for f, bw in formants[:5]]}")
```

### Java — ZCR + short-time energy

```java
public class SpeechParameters {

    public static double shortTimeEnergy(double[] frame) {
        double E = 0;
        for (double x : frame) E += x * x;
        return E;
    }

    public static double shortTimeZCR(double[] frame) {
        int crossings = 0;
        for (int n = 1; n < frame.length; n++) {
            if ((frame[n - 1] >= 0) != (frame[n] >= 0)) crossings++;
        }
        return (double) crossings / (frame.length - 1);
    }

    public static double[] autocorrelation(double[] frame, int maxLag) {
        int N = frame.length;
        double[] R = new double[maxLag + 1];
        for (int k = 0; k <= maxLag; k++) {
            double s = 0;
            for (int n = 0; n < N - k; n++) s += frame[n] * frame[n + k];
            R[k] = s;
        }
        return R;
    }

    public static double[] pitchFromACF(double[] frame, double fs,
                                         double fMin, double fMax) {
        int kMin = (int) (fs / fMax);
        int kMax = (int) (fs / fMin);
        double[] R = autocorrelation(frame, kMax);
        int kPeak = kMin;
        for (int k = kMin; k <= kMax; k++)
            if (R[k] > R[kPeak]) kPeak = k;
        boolean voiced = R[kPeak] / R[0] > 0.3;
        return new double[]{ voiced ? fs / kPeak : 0.0, voiced ? 1.0 : 0.0 };
    }

    public static void main(String[] args) {
        int N = 400;
        double[] frame = new double[N];
        for (int n = 0; n < N; n++)
            frame[n] = Math.sin(2 * Math.PI * 150 * n / 16000.0);

        System.out.printf("Energy: %.2f%n", shortTimeEnergy(frame));
        System.out.printf("ZCR:    %.4f%n", shortTimeZCR(frame));
        double[] p = pitchFromACF(frame, 16000, 60, 400);
        System.out.printf("Pitch:  %.1f Hz (voiced=%b)%n", p[0], p[1] == 1.0);
    }
}
```

### C# — Mel filter bank

```csharp
using System;

class MelFilterBank
{
    static double HzToMel(double f) => 2595.0 * Math.Log10(1 + f / 700.0);
    static double MelToHz(double m) => 700.0 * (Math.Pow(10, m / 2595.0) - 1);

    // Build triangular filters evenly spaced in mel space
    static double[][] BuildMelFilters(int nFilters, int nFft, double fs,
                                       double fMin = 0, double fMax = 8000)
    {
        double melMin = HzToMel(fMin), melMax = HzToMel(fMax);
        double melStep = (melMax - melMin) / (nFilters + 1);
        var filters = new double[nFilters][];
        for (int m = 0; m < nFilters; m++)
        {
            double centerMel = melMin + (m + 1) * melStep;
            double leftMel = centerMel - melStep;
            double rightMel = centerMel + melStep;
            double centerHz = MelToHz(centerMel);
            double leftHz = MelToHz(leftMel);
            double rightHz = MelToHz(rightMel);
            int binLeft = (int)(leftHz / fs * nFft);
            int binCenter = (int)(centerHz / fs * nFft);
            int binRight = (int)(rightHz / fs * nFft);
            filters[m] = new double[nFft / 2 + 1];
            for (int k = binLeft; k < binCenter; k++)
                filters[m][k] = (double)(k - binLeft) / (binCenter - binLeft);
            for (int k = binCenter; k < binRight && k < filters[m].Length; k++)
                filters[m][k] = (double)(binRight - k) / (binRight - binCenter);
        }
        return filters;
    }

    static void Main()
    {
        int nFft = 512;
        double fs = 16000;
        int nFilters = 24;
        var filters = BuildMelFilters(nFilters, nFft, fs);
        for (int m = 0; m < nFilters; m += 4)
        {
            int center = Array.IndexOf(filters[m], filters[m].Max());
            Console.WriteLine($"Filter {m,2}: peak at bin {center,3} = {center * fs / nFft:F0} Hz");
        }
    }
}
```

**Expected output for the mel bank.**

```
Filter  0: peak at bin   1 = 31 Hz
Filter  4: peak at bin   5 = 156 Hz
Filter  8: peak at bin  14 = 438 Hz
Filter 12: peak at bin  30 = 938 Hz
Filter 16: peak at bin  58 = 1813 Hz
Filter 20: peak at bin 105 = 3281 Hz
```

The filters are **denser at low frequency** — the mel spacing in action.

!!! tip "Tooling"
    - **Python:** `librosa` for STFT, mel spectrogram, MFCC, pitch (`pyin`); `scipy.signal` for filter design and windows; `python_speech_features` for classic MFCC; `praat-parselmouth` for phonetics-grade formant/pitch tracking.
    - **Java:** JTransforms (FFT); TarsosDSP (`PitchDetection`, `MFCC`); no standard LPC library — implement Levinson manually.
    - **C#:** Math.NET Numerics (FFT); NAudio (audio I/O); hand-roll MFCC/LPC.
    - **MATLAB:** `spectrogram`, `mfcc`, `lpc`, `levinson`, `pitchest`, `findpeaks` — the reference.
    - **Praat:** the standard tool for phonetic analysis — formants, pitch, spectrogram.

---

## 4.15 Chapter summary

1. **Short-time analysis** is the foundation — speech is quasi-stationary over 10–30 ms. Frame + window (Hamming, 25 ms, 50% overlap).
2. **Time-domain parameters:** energy, magnitude, ZCR, autocorrelation — cheap and useful for voiced/unvoiced, endpoint, and (with ACF) pitch.
3. **Parallel processing pitch detection** — autocorrelation of downsampled channels, combined by voting — resolves octave errors.
4. **STFT** is the frequency-domain workhorse; window length controls the **time-frequency resolution trade-off** (Heisenberg); narrowband = 40 ms (pitch), wideband = 10 ms (timing).
5. **Filter banks** — uniform, or **mel-spaced** for perceptual matching. The mel bank is the front end of MFCC.
6. **Spectrograms** — narrowband shows harmonics (pitch); wideband shows pitch pulses; formants are visible in both. Phoneticians read them for phone sequences.
7. **Formant extraction** — spectral peak picking (with smoothing) or LPC root finding. LPC gives formant **frequency and bandwidth** from the pole's angle and magnitude.
8. **Pitch extraction** — ACF (time), cepstrum (frequency), HPS (harmonic product spectrum), LPC residual + SIFT. Each has trade-offs; cepstrum and SIFT are the most robust.
9. **Analysis–synthesis** — the **channel vocoder** is the classic example; modern codecs (CELP) improve on it by transmitting the LPC residual rather than regenerating it from a source model.
10. **Cepstrum** separates source (high quefrency, pitch peak) from filter (low quefrency, envelope). It's the basis for **MFCC**, the standard recognition feature.
11. **LPC** — AR(p) all-pole model of the vocal tract; **Yule–Walker equations**; **Levinson–Durbin** solves them in $O(p^2)$; **reflection coefficients** $k_m$ characterize stability.
12. **LPC formant extraction** — pole angle = formant frequency; pole radius = bandwidth. Filter $F < f_s/2$ and BW < 500 Hz to reject spurious roots.

---

## 4.16 Quiz

### A. Multiple choice

1. The standard **frame length** for short-time speech analysis is:
   (a) 1 ms (b) 10 ms (c) 25 ms (d) 100 ms

2. **Short-time zero-crossing rate** is highest for:
   (a) vowels (b) silence (c) unvoiced fricatives (d) nasals

3. In the **autocorrelation** function, the peak at $k = T_0$ indicates:
   (a) silence (b) voiced speech with pitch period $T_0$ (c) unvoiced noise (d) DC

4. **Narrowband spectrograms** use a window of approximately:
   (a) 5 ms (b) 10 ms (c) 40 ms (d) 100 ms

5. A **mel filter bank** places more filters:
   (a) uniformly (b) at high frequencies (c) at low frequencies (d) at random

6. **Formants** can be extracted most reliably from:
   (a) the raw waveform (b) energy only (c) LPC roots or smoothed spectral peaks (d) ZCR

7. The **cepstrum** is the:
   (a) FFT of the spectrum (b) IFFT of the log-magnitude spectrum (c) squared magnitude of the spectrum (d) DFT of the autocorrelation

8. In the cepstrum, the **pitch peak** appears at quefrency equal to:
   (a) $F_0$ (b) $T_0$ (samples) (c) half the frame length (d) $N/2$

9. The **Yule–Walker** equations solve for:
   (a) pitch (b) LPC coefficients (c) cepstrum (d) spectrogram

10. **Levinson–Durbin** solves the Yule–Walker system in:
    (a) $O(p)$ (b) $O(p^2)$ (c) $O(p^3)$ (d) $O(2^p)$

11. The **channel vocoder** transmits:
    (a) the raw waveform (b) filter-bank envelopes plus pitch/voicing (c) LPC coefficients only (d) the cepstrum

12. For $p = 10$ LPC analysis of a vowel, the number of formants is typically:
    (a) 1 (b) 2 (c) 5 (d) 10

### B. Fill in the blank

1. The standard frame length for short-time analysis is \_\_\_\_ ms.
2. Window length controls the trade-off between \_\_\_\_ and \_\_\_\_ resolution.
3. The \_\_\_\_ window is standard for STFT analysis (43 dB sidelobe).
4. Autocorrelation peak at lag $k$ indicates periodicity with period $k$ — used for \_\_\_\_ estimation.
5. Mel filter banks place \_\_\_\_ filters at low frequencies.
6. The \_\_\_\_ is the IFFT of the log-magnitude spectrum.
7. In the cepstrum, the pitch peak is at quefrency \_\_\_\_ (in samples).
8. LPC models the vocal tract as an \_\_\_\_ filter.
9. Yule–Walker equations are solved by \_\_\_\_–\_\_\_\_ in $O(p^2)$.
10. LPC formant frequency is the \_\_\_\_ of the LPC pole; bandwidth is the \_\_\_\_ of its magnitude.

### C. Define / explain

1. Explain why **short-time analysis** is necessary for speech. What is the standard frame length and hop size, and why?
2. Compute the **short-time energy**, **magnitude**, and **zero-crossing rate** of $x = (2, -1, 1, -3, -2, 1, -1, 2)$.
3. Define the **short-time autocorrelation** and explain how it is used for **pitch detection**. Why is **parallel processing** needed?
4. State the **short-time Fourier transform** and derive the **time-frequency resolution trade-off**. What window length is used for narrowband vs. wideband spectrograms?
5. Describe the **mel filter bank**. Why is the mel scale used, and what is the typical number of filters?
6. Interpret a **spectrogram** of a voiced vowel, an unvoiced fricative, and a plosive. What features distinguish them?
7. Define the **cepstrum** and explain the **source-filter separation** in the cepstral domain. Why does the pitch peak appear at high quefrency and the formants at low quefrency?
8. State the **LPC model** of speech. Derive the **Yule–Walker equations** from the minimum-MSE criterion.
9. State the **Levinson–Durbin recursion**. What are the **reflection coefficients** and what do they mean physically?
10. Describe **LPC-based formant extraction**. How are formants extracted from the LPC polynomial roots?
11. Compare **autocorrelation**, **cepstrum**, **HPS**, and **LPC residual (SIFT)** for pitch detection. Which is best for which situation?

---

## 4.17 Answer key

??? success "A. Multiple choice"
    1. **(c)** — 25 ms.
    2. **(c)** — unvoiced fricatives.
    3. **(b)** — voiced, pitch period $T_0$.
    4. **(c)** — 40 ms.
    5. **(c)** — low frequencies.
    6. **(c)** — LPC roots or smoothed spectral peaks.
    7. **(b)** — IFFT of log-magnitude.
    8. **(b)** — $T_0$ samples.
    9. **(b)** — LPC coefficients.
    10. **(b)** — $O(p^2)$.
    11. **(b)** — envelopes + pitch/voicing.
    12. **(c)** — 5 (one per pole pair).

??? success "B. Fill in the blank"
    1. 20–30 (typically 25)
    2. time; frequency
    3. Hamming
    4. pitch
    5. more
    6. cepstrum
    7. $T_0$ (pitch period in samples)
    8. all-pole (AR)
    9. Levinson–Durbin
    10. angle; negative log

??? success "C. Define / explain — model answers"
    1. *Why short-time:* speech is non-stationary (vocal tract moves continuously), but quasi-stationary over 10–30 ms — long enough to contain 2–3 pitch periods and thus support a stable spectral estimate, short enough that the articulator configuration is roughly frozen. *Standard parameters:* frame length 25 ms (400 samples at 16 kHz), hop size 10 ms (160 samples) or 12.5 ms (200 samples, 50% overlap). The 50% overlap smooths the transition between frames (continuous feature stream), matches auditory integration (~20 ms), and improves vocoder quality.
    2. *Computation for $x = (2, -1, 1, -3, -2, 1, -1, 2)$:*
       Energy = $4+1+1+9+4+1+1+4 = 25$.
       Magnitude = $2+1+1+3+2+1+1+2 = 13$.
       Sign changes: $(+,-)$ (1), $(-,+)$ (2), $(+,-)$ (3), $(-,-)$ no, $(-,+)$ (4), $(+,-)$ (5), $(-,+)$ (6). ZCR = $6/2 = 3$; normalized $3/7 \approx 0.43$ (using $N-1 = 7$ as denominator).
    3. *Short-time autocorrelation:* $R_m[k] = \sum_n x_m[n] x_m[n+k]$ — measure of self-similarity at lag $k$. *Pitch detection:* for voiced speech, $R_m[k]$ peaks at $k = T_0$, the pitch period; find the first peak above a threshold in the plausible pitch range ($f_s/400$ to $f_s/60$ samples). *Parallel processing:* the ACF has peaks at $T_0, 2T_0, 3T_0, \ldots$ and (for some signals) at $T_0/2$ — leading to **octave errors** (pitch doubled or halved). Parallel processing computes ACF at multiple decimation factors and combines results by **voting**: a candidate pitch is accepted only if confirmed by several channels. This dramatically reduces octave errors at the cost of more computation.
    4. *STFT:* $X_m(e^{j\omega}) = \sum_n x_m[n] w[n] e^{-j\omega n}$ — Fourier transform of each windowed frame. *Trade-off:* window length $N$ determines both frequency resolution ($\Delta f = f_s/N$) and time resolution ($\Delta t \approx N/f_s$). Their product is bounded by the Heisenberg uncertainty principle: $\Delta t \cdot \Delta f \ge 1/(4\pi)$. *Narrowband:* 40 ms (long window, $\Delta f \approx 25$ Hz) — resolves pitch harmonics, useful for formant tracking. *Wideband:* 10 ms (short window, $\Delta f \approx 100$ Hz) — resolves pitch pulses as vertical stripes, useful for voicing and timing.
    5. *Mel filter bank:* a set of triangular bandpass filters whose **center frequencies** are equally spaced in **mel**, not Hz. The mel scale $f_{\text{mel}} = 2595\log_{10}(1 + f_{\text{Hz}}/700)$ approximates the ear's nonlinear frequency sensitivity. Because the ear has finer resolution at low frequencies (up to ~500 Hz) and coarser at high (~8 kHz), the mel bank places more filters where perceptual resolution is higher. *Typical:* 20–40 filters covering 0–8000 Hz. The output energies (often log-transformed and DCT'ed) form the basis of MFCC features.
    6. *Voiced vowel:* horizontal dark formant bands (e.g., $F_1 \approx 500$ Hz, $F_2 \approx 1500$ Hz, $F_3 \approx 2500$ Hz); vertical striations showing pitch pulses; energy concentrated below 4 kHz. *Unvoiced fricative:* diffuse broadband high-frequency energy (above 3 kHz), no formants, no pitch striations. *Plosive:* short vertical burst (30–50 ms) at a specific frequency region, preceded by a silent closure. *Distinguishing features:* formant bands (vowel), voicing striations (voiced vs. unvoiced), burst/closure (plosive).
    7. *Cepstrum:* $c[n] = \text{IFFT}\{\log|X[k]|\}$ — inverse Fourier transform of the log-magnitude spectrum. *Source-filter separation:* in the log spectrum, the source (periodic pulse train) and the filter (vocal-tract formants) are **additive**. In the cepstral domain, they occupy **different quefrency ranges**: the filter (slow variation) lives at **low quefrency** (near $n=0$); the source (rapid harmonic structure) lives at **high quefrency** (at the pitch period $n = T_0$). *Why:* the log-magnitude spectrum of the source has comb structure with period $F_0$; the IFFT of a comb with period $F_0$ is a peak at quefrency $T_0 = 1/F_0$. The formant envelope is smooth, so its IFFT is concentrated at low quefrency.
    8. *LPC model:* $x[n] = -\sum_{k=1}^p a_k x[n-k] + e[n]$ where $e$ is a white driving source. Equivalently, $x$ is the output of an all-pole AR filter $H(z) = 1/A(z)$ with $A(z) = 1 + \sum a_k z^{-k}$, driven by $e$. *Derivation:* minimise $\mathcal E = \sum e[n]^2$ with $e[n] = x[n] + \sum a_k x[n-k]$. Take partial derivatives $\partial\mathcal E/\partial a_k = 0$: $2\sum e[n] x[n-k] = 0$ for $k = 1, \ldots, p$. Substitute $e[n] = x[n] + \sum_l a_l x[n-l]$: $R[k] + \sum_l a_l R[k-l] = 0$, i.e., $\sum_l a_l R[k-l] = -R[k]$. In matrix form: $\mathbf R \mathbf a = -\mathbf r$ where $R_{ij} = R[|i-j|]$. These are the Yule–Walker equations.
    9. *Levinson–Durbin:* an $O(p^2)$ recursion for solving the Toeplitz system $\mathbf R \mathbf a = -\mathbf r$. At each order $m$ it computes a **reflection coefficient** $k_m$ and updates the AR coefficients of order $m$ from order $m-1$: $k_m = -(R[m] + \sum a_i^{(m-1)} R[m-i])/\sigma_{m-1}^2$; $a_i^{(m)} = a_i^{(m-1)} + k_m a_{m-i}^{(m-1)}$; $\sigma_m^2 = (1 - k_m^2)\sigma_{m-1}^2$. *Reflection coefficients:* $k_m$ is the **partial autocorrelation** at lag $m$ — the correlation between $x[n]$ and $x[n-m]$ after removing the linear effect of lags $1, \ldots, m-1$. Physically, it measures how much energy is removed by adding lag $m$ to the predictor. The condition $|k_m| < 1$ for all $m$ is necessary and sufficient for the AR filter to be stable (all poles inside the unit circle).
    10. *LPC formant extraction:* given LPC coefficients $A(z) = 1 + \sum_{k=1}^p a_k z^{-k}$, find its $p$ roots $z_i$ (numerically, via polynomial root-finding). Complex roots come in conjugate pairs, one pair per formant. For each pair:
        $F = (f_s/2\pi) \cdot \angle z_i$ (formant frequency);
        $\text{BW} = -(f_s/\pi) \cdot \ln|z_i|$ (formant bandwidth).
        Filter the results: keep $F$ in $[50, f_s/2]$, reject pairs with BW > 500 Hz (spurious), sort by frequency, keep the first 3–5 formants. Typical: $p = 10$ gives 5 formants; only the first 3–4 are reliably meaningful.
    11. *Pitch detection comparison:*
        - **ACF (autocorrelation):** cheap, works for most signals, but **octave errors** with some waveforms. Good baseline.
        - **Cepstrum:** robust to formants (filter lives at low quefrency), handles weak fundamentals, but requires longer frames (40+ ms) and FFT cost. Good for many speakers.
        - **HPS (harmonic product spectrum):** simple, no octave errors in principle, but limited by FFT resolution. Good for clean periodic signals.
        - **SIFT (LPC residual + ACF):** best accuracy for speech, because the LPC filter whitens the residual and leaves only the pitch peak. Standard in codecs and modern analysis.
        *Which to use:* **SIFT** for coding and highest accuracy; **cepstrum** for analysis with clear pitch; **ACF** for embedded/cheap; **HPS** rarely (mostly historical). **YIN** (2002) is the modern open-source standard, better than all of these in practice.

---

## Further reading

- Rabiner & Schafer, *Digital Processing of Speech Signals*, Prentice Hall, 1978 — **the classical reference** for Chapters 3–4 of this course. Sections 4.2–4.5 for STFT and spectrograms; 4.3 for pitch; 5.4 for cepstrum; 6.3 for LPC.
- Deller, Proakis & Hansen, *Discrete-Time Processing of Speech Signals*, IEEE Press, 2000 — the modern comprehensive reference.
- Quatieri, *Discrete-Time Speech Signal Processing*, Prentice Hall, 2001 — excellent on homomorphic methods, cepstrum, and vocoders.
- O'Shaughnessy, *Speech Communications: Human and Machine*, IEEE Press, 2nd ed., 2000 — good on production and perception.
- Makhoul, "Linear Prediction: A Tutorial Review", *Proc. IEEE*, 63(4), 1975 — the definitive pedagogical reference for LPC.
- Markel & Gray, *Linear Prediction of Speech*, Springer, 1976 — the classical LPC monograph (including SIFT).
- Rabiner, "On the Use of Autocorrelation Analysis for Pitch Detection", *IEEE Trans. ASSP*, 25(1), 1977 — the parallel-processing paper.
- Hess, *Pitch Determination of Speech Signals*, Springer, 1983 — the classic pitch reference.
- de Cheveigné & Kawahara, "YIN, a Fundamental Frequency Estimator for Speech and Music", *JASA*, 111(4), 2002 — the modern standard.
- Childers, Skinner & Kemerait, "The Cepstrum: A Guide to Processing", *Proc. IEEE*, 65(10), 1977 — the classic cepstrum reference.
- Davis & Mermelstein, "Comparison of Parametric Representations for Monosyllabic Word Recognition in Continuously Spoken Sentences", *IEEE Trans. ASSP*, 28(4), 1980 — the original MFCC paper.
- Hermansky, "Perceptual Linear Predictive (PLP) Analysis of Speech", *JASA*, 87(4), 1990 — the PLP alternative to MFCC.
- Mermelstein, "Distance Measures for Speech Recognition, Psychological and Instrumental", *Pattern Recognition and Artificial Intelligence*, 1976 — the mel scale in speech.
- Gold & Morgan, *Speech and Audio Signal Processing*, Wiley, 2000 — general textbook.
- `librosa` documentation (Python); `python_speech_features`; `praat`; `HTK`; `kaldi`; MATLAB Audio Toolbox.

---

## Suggested tutorial and lab problems for this chapter

1. **Frame-based features.** Record a sentence with mixed voiced/unvoiced segments. Compute short-time energy and ZCR for each frame. Plot both on the same time axis — can you see the voiced/unvoiced distinction clearly? Where does each feature fail?

2. **Pitch from ACF.** Implement autocorrelation-based pitch detection on a sustained vowel. Compare with a reference pitch tracker (e.g., `librosa.pyin`). Measure the **octave error rate** — on how many frames does simple ACF give a wrong result?

3. **Parallel-processing pitch.** Implement a two-channel (decimate-by-2 and decimate-by-4) pitch detector. Verify it rejects octave errors on the same recording.

4. **Narrowband vs. wideband spectrogram.** Compute both spectrograms for the same sentence. Identify features visible in one but not the other. Save both plots side by side.

5. **Mel filter bank.** Implement a 24-filter mel bank, apply to 5 frames of different vowels, and report the 24 energy values. Verify that low-frequency channels differ more between vowels than high-frequency channels.

6. **Formant extraction.** Record vowels /a/, /e/, /i/, /o/, /u/. Extract formants by (a) LPC root finding with $p = 10$, (b) smoothed-spectrum peak picking. Compare. Plot on an $(F_1, F_2)$ chart — do the five vowels cluster as expected?

7. **Cepstral pitch.** Compute the real cepstrum of a voiced frame. Locate the pitch peak. Apply a **low-pass lifter** (keep $n < 50$); inverse-transform → the smooth spectral envelope. Overlay on the original log spectrum — the envelope should be smooth and follow the formants.

8. **MFCC implementation.** Implement the full MFCC pipeline: pre-emphasis → framing → Hamming → FFT → mel bank → log → DCT. Extract 13 MFCCs per frame for a sentence. Visualize as a 13 × N feature matrix. Compare with `librosa.feature.mfcc`.

9. **LPC estimation.** Implement Levinson–Durbin from scratch. Verify against `scipy.signal.levinson` and `statsmodels.tsa` on a random Toeplitz system. Report timing.

10. **LPC spectrum.** For a vowel frame, compute LPC with $p = 10$ and plot the LPC spectral envelope $1/|A(e^{j\omega})|^2$ overlaid on the FFT magnitude spectrum. The peaks of the envelope should align with the formants.

11. **Pitch from LPC residual.** Compute the LPC residual for a voiced frame. Plot the residual and its autocorrelation. Verify the pitch peak is sharper than in the original signal. Compare pitch estimates from ACF (original) vs. ACF (residual).

12. **Channel vocoder.** Implement a 16-channel vocoder with pitch/voicing extraction and re-synthesis. Compare the output with the original. The output should sound "robotic" but intelligible — this is the classic vocoder sound.

13. **Spectral envelope from cepstrum.** Compute the cepstrum of a mixed voiced/unvoiced frame. Apply both low-pass and high-pass lifters. Reconstruct separately. Verify that low-pass gives the formant envelope and high-pass gives the pitch pulses.

14. **Real-time analysis.** Process a live microphone stream: compute the spectrogram and track pitch in real time. Report the latency (should be < 50 ms for interactive use). Use `sounddevice` or `pyaudio`.