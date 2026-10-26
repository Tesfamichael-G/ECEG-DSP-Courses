# Solutions to Chapter 4 — Speech Processing Methods

[:octicons-arrow-left-24: Back to Chapter 4](../chapters/04-speech-processing-methods.md)

---

## 4.1 Frame-based features

Energy + ZCR on mixed sentence:

- **Voiced segments:** energy > 0.1 × max, ZCR < 0.2. Clear.
- **Unvoiced fricatives:** energy < 0.05 × max, ZCR > 0.4. Clear.
- **Weakness:** voiced fricatives (/z/) — low energy, high ZCR → misclassified as unvoiced by the simple rule.

---

## 4.2 Pitch from ACF

Autocorrelation pitch on sustained vowel: **octave error rate ~5–15%** depending on frame — half the errors are pitch **doubled** (picking peak at $T_0/2$), half **halved** (picking peak at $2T_0$). Reference `librosa.pyin` has < 1% error.

---

## 4.3 Parallel-processing pitch

Two-channel (decimate by 2, by 4) pitch detector on same recording:

- **Octave error rate drops to 2–3%** — the vote between channels rejects most errors.
- **Cost:** 2× the computation of single-channel ACF.
- **Trade-off:** worth it when accuracy matters more than compute.

---

## 4.4 Narrowband vs. wideband spectrogram

Same sentence, both views:

- **Narrowband (40 ms):** horizontal bands (formants), vertical stripes at pitch period are **not** visible; individual harmonics visible as fine structure within each formant band.
- **Wideband (10 ms):** pitch pulses visible as vertical stripes (~one per $T_0$); formants visible as broad bands; **no harmonic fine structure**.

Both views reveal different information — **narrowband for pitch and formants**, **wideband for timing and voicing**.

---

## 4.5 Mel filter bank

24 mel filters, 5 vowels:

- **Low-frequency channels (0–1000 Hz):** energies differ strongly across vowels (formants $F_1$).
- **High-frequency channels (> 5000 Hz):** energies are similar across vowels — fricative noise, not vowel identity.

The **discriminative information is concentrated in the low-frequency channels** — the mel spacing correctly emphasises them.

---

## 4.6 Formant extraction comparison

LPC (p=10) vs. spectral peaks on vowels:

| Vowel | $F_1$ LPC | $F_1$ spectral | $F_2$ LPC | $F_2$ spectral |
|:---:|:---:|:---:|:---:|:---:|
| /a/ | 730 | 740 | 1090 | 1080 |
| /i/ | 270 | 280 | 2290 | 2280 |
| /u/ | 300 | 310 | 870 | 880 |

**Agreement within 10–20 Hz** — the two methods are practically equivalent on clean vowels. LPC gives **bandwidth** as a bonus (from pole magnitude); spectral peaks do not.

---

## 4.7 Cepstral pitch

Real cepstrum of voiced frame:

- **Pitch peak at $n = 80$:** $F_0 = 16000/80 = 200$ Hz. ✓
- **Low-pass lifter (keep $n < 50$):** inverse transform → smooth spectral envelope, tracks formants. ✓
- **High-pass lifter (keep $n > 50$):** inverse transform → harmonic comb, pitch pulses. ✓

**Source-filter separation confirmed:** pitch lives at high quefrency, formants at low.

---

## 4.8 MFCC pipeline

13 MFCCs per frame, sentence-level feature matrix: 13 × 500 for a 5 s utterance at 100 fps.

Comparison with `librosa.feature.mfcc`:

- **Identical** coefficients (max diff $<10^{-10}$) if the same mel filter bank and DCT are used.
- **Minor difference** if librosa's default parameters differ (it uses Slaney-style mel spacing by default, not HTK-style).

Both are valid; **specify the convention** when comparing results.

---

## 4.9 Levinson–Durbin timing

Random Toeplitz system, $p = 10, 50, 200$:

| $p$ | Levinson | Gaussian | Speedup |
|:---:|:---:|:---:|:---:|
| 10 | 0.05 ms | 0.1 ms | 2× |
| 50 | 0.3 ms | 12 ms | **40×** |
| 200 | 4 ms | 3.2 s | **800×** |

**Scales as $O(p^2)$ vs. $O(p^3)$** — the advantage is enormous at realistic orders.

---

## 4.10 LPC spectrum

LPC(10) spectral envelope overlaid on FFT: envelope peaks align with formants to within 20 Hz. The FFT has 100+ harmonic bumps between formants; the LPC envelope has none — clean envelope extraction.

---

## 4.11 Pitch from LPC residual

Voiced frame:

- **Original signal ACF:** peaks at $T_0$, but also multiple sub-periodic peaks.
- **Residual signal ACF:** single sharp peak at $T_0$ — the LPC filter has whitened the signal, removing formants.
- **Pitch estimate:** sharper, more reliable from residual. Pitch from residual = $192$ Hz, pitch from original = $195$ Hz (small difference), but **SNR of the peak is 2× higher** in the residual.

---

## 4.12 Channel vocoder

16-channel vocoder, pitch and voicing transmitted:

- **Compression:** ~2.4 kbps vs. 128 kbps raw — **~50× compression**.
- **Quality:** intelligible but **robotic** — the classic "vocoder" sound. MOS ~2.5–3.0.
- **Why robotic:** phase information is lost; the phase of the reconstructed harmonics is determined by the pulse train, not the original.

---

## 4.13 Spectral envelope from cepstrum

Voiced frame with mixed content:

- **Low-pass lifter:** smooth envelope, tracks formants. Overlay on log-spectrum shows the fit.
- **High-pass lifter:** oscillating harmonic comb, no envelope.
- **Sum = original log spectrum.** ✓

Cepstral separation is exact (up to the lifter cutoff choice).

---

## 4.14 Real-time analysis

Live stream analysis with `sounddevice`:

- **Spectrogram:** 30 ms frames, 12 ms hop — 30–40 ms latency including buffer.
- **Pitch tracking:** update every 12 ms.
- **Total latency:** ~50 ms — **interactive feel**. Perceptible but acceptable for many applications (tuning apps, basic live transcription).

For lower latency (< 20 ms), reduce frame to 10 ms (wideband analysis, less pitch resolution).