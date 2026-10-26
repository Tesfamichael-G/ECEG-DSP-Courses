# Solutions to Chapter 3 — Speech Processing Introduction

[:octicons-arrow-left-24: Back to Chapter 3](../chapters/03-speech-processing-intro.md)

---

## 3.1 Pitch period measurement

Sustained /a/ at 16 kHz. In a 100 ms window, count ~14 peaks → $F_0 = 140$ Hz, $T_0 = 7.14$ ms = 114 samples. `librosa.pyin` returns ~137–142 Hz for the same frame — agreement within 2%.

---

## 3.2 Formant estimation

For adult male vowels, spectral peaks (LPC, p=10):

| Vowel | $F_1$ (Hz) | $F_2$ (Hz) |
|:---:|:---:|:---:|
| /i/ (beat) | 270 | 2290 |
| /e/ (bet) | 530 | 1840 |
| /a/ (father) | 730 | 1090 |
| /o/ (bought) | 570 | 840 |
| /u/ (boot) | 300 | 870 |

Plotted in $(F_1, F_2)$ space, they trace the classic **vowel quadrilateral** — high $F_1$ = low vowel, high $F_2$ = front vowel.

---

## 3.3 Voiced/unvoiced classification

Energy + ZCR rule on mixed recording:

- **Voiced (vowels, nasals):** energy > 10% of max, ZCR < 0.15.
- **Unvoiced (fricatives):** energy < 10%, ZCR > 0.3.
- **Plosives:** brief energy spike + ZCR fluctuation.

Rule-based classifier achieves **~90% accuracy**. Errors are mostly on voiced fricatives (/z/, /v/) that have both high energy and moderate ZCR.

---

## 3.4 Acoustic tube model

Kelly–Lochbaum: 10 sections, area vector $A_k$. Neutral tract ($A_k$ = const): formants at 500, 1500, 2500, 3500 Hz. Moving tongue forward (narrow $A_{front}$): $F_2$ rises to ~2000 Hz, $F_1$ drops — matching /i/. Moving tongue back: $F_2$ drops to ~900 Hz, $F_1$ rises — matching /u/. The model **reproduces the vowel quadrilateral** from pure physics.

---

## 3.5 LPC analysis

10th-order LPC spectrum overlaid on FFT: the LPC envelope traces the formant peaks to within 20–50 Hz, with no harmonic ripple. **LPC = smoothed spectral envelope** in physical terms.

---

## 3.6 Pre-emphasis

Applying $1 - 0.95 z^{-1}$ flattens the −12 dB/octave glottal roll-off. Effect on spectrum:

- **Before:** energy concentrated at low frequencies, high formants hard to see.
- **After:** flatter spectrum, all 4–5 formants visible as distinct peaks.

---

## 3.7 Framing and windowing

25 ms frames, 50% overlap. With **no windowing:** the spectrogram shows horizontal bands — the frame boundaries create artificial discontinuities, which appear as extra frequency content. With **Hamming window:** clean spectrogram, no frame-boundary artifacts.

---

## 3.8 Data rate

| Source | Duration | G.711 (64 kbps) | G.729 (8 kbps) | Opus (12 kbps) |
|---|:---:|:---:|:---:|:---:|
| Phone call | 1 min | 480 KB | 60 KB | 90 KB |
| Podcast | 30 min | 14.4 MB | 1.8 MB | 2.7 MB |
| Audiobook | 90 min | 43.2 MB | 5.4 MB | 8.1 MB |

For a mobile download, G.729 (5.4 MB) is best; for a music-quality download, G.711 is unacceptable (voice only) — use Opus at higher rate.

---

## 3.9 Speech bandwidth

Downsample 16 kHz → 8 kHz → 4 kHz:

- **16 kHz:** full quality.
- **8 kHz:** "telephone-like", still fully intelligible.
- **4 kHz (no anti-alias):** severe aliasing, sibilants mangled, but vowels still intelligible.

With anti-alias filter at 4 kHz: intelligible but muffled. **The Nyquist rule** matters — without the filter, aliasing destroys the signal.

---

## 3.10 Speaker features

For two speakers saying the same sentence:

- **MFCC:** 12 coefficients/frame; speakers differ most in coefficients 2–5 (formant spacing).
- **$F_0$:** male ~120 Hz, female ~200 Hz — huge difference, but speaker-individual variation is 30–50 Hz so it only partially discriminates.
- **Speaking rate:** 3–5 syllables/s varies by speaker.
- **GMM on MFCC:** ~95% speaker classification accuracy for 2 speakers, drops to ~70% for 20.

---

## 3.11 Codec comparison

Same sentence through G.711, G.729, Opus at equivalent rates:

- **G.711:** indistinguishable from source.
- **G.729:** slight buzzy quality, especially on fricatives, but fully intelligible.
- **Opus @ 12 kbps:** comparable to G.729 but smoother — modern codec design.

MOS ≈ 4.5 / 4.0 / 4.2 respectively.

---

## 3.12 Reading machine

Typical TTS rate: 150–200 words/min. Human reading: ~150 words/min. TTS at normal speed matches human reading; **accelerated TTS** (2–4×) is intelligible with training — many blind users read at 300+ words/min, far faster than sighted readers.