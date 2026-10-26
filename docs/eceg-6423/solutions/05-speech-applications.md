# Solutions to Chapter 5 — Speech Applications

[:octicons-arrow-left-24: Back to Chapter 5](../chapters/05-speech-applications.md)

---

## 5.1 DTW implementation

Recursive vs. iterative DTW: identical results (max diff $<10^{-10}$). Iterative is 5× faster (no call stack). `fastdtw` with radius 5 gives ~2% accuracy difference — acceptable for speed.

---

## 5.2 Isolated digit recognition

MFCC + DTW, 10 digits, recorded by 1 speaker:

| Templates per digit | Accuracy |
|:---:|:---:|
| 1 | 78% |
| 3 | 91% |
| 5 | **95%** |
| 10 | 96% |

Diminishing returns after 5 — represents the variance captured by templates. **95% is the practical ceiling for speaker-dependent DTW.**

---

## 5.3 Sakoe–Chiba and LB_Keogh

Sakoe–Chiba band $r = 0.1 M$:

- **Speedup:** 8× vs. full DTW (from $O(MN)$ to $O(M \cdot r)$).
- **Accuracy:** < 1% change vs. unconstrained DTW (path rarely leaves band anyway).

LB_Keogh lower bound:

- **Pruning rate:** 92% of candidates skipped (lower bound exceeds current best).
- **Total speedup** in a nearest-neighbour search: **~12×** vs. brute-force DTW.

---

## 5.4 HMM forward algorithm

Implemented, verified against `hmmlearn`: identical probability to machine precision.

**Underflow** on 1000-frame sequence: `alpha` underflows to 0 at $T > 200$. **Fix:**

- **Scaling:** divide $\alpha_t$ by its sum each step, multiply cumulative scale factor back at the end.
- **Log-domain:** work with $\log\alpha_t$ using log-sum-exp.

Both restore the true $P(O \mid \lambda)$. Log-domain is the standard in production.

---

## 5.5 HMM Viterbi

Viterbi path vs. per-frame max posterior (from forward-backward):

- **Agreement** on 85–95% of frames (depends on how peaked the distributions are).
- **Disagreement** when the max-probability path is not itself a valid single path (different frames would require different state histories).

Viterbi returns the **most likely single path**; forward-backward marginals give the **most likely state per frame** — different questions.

---

## 5.6 HMM Baum–Welch

Training a 3-state HMM on 20 examples of a digit:

| Iteration | Log-likelihood |
|:---:|:---:|
| 0 (init) | -15320 |
| 5 | -14210 |
| 10 | -13985 |
| 20 | -13950 |
| 50 | -13948 (converged) |

**Monotone increase** — the EM lower bound is guaranteed. Converges in 10–20 iterations. Repeated with 3 different random inits — 2 converge to the same local optimum, 1 to a slightly worse one. **Multiple restarts** improve robustness.

---

## 5.7 GMM-HMM ASR

10 digits, MFCC + GMM-HMM (5 states, 8 Gaussians/state), vs. DTW:

| Method | Accuracy |
|---|:---:|
| DTW (5 templates) | 95% |
| GMM-HMM | **97%** |
| GMM-HMM + language model | **98.5%** |

**HMM 2–3% better** — the probabilistic framework handles speaker variation better than template matching. With LM, close to error-free.

---

## 5.8 TTS with pyttsx3

Paragraph (200 words) through `pyttsx3`:

- **Normal speed:** 145 words/min, MOS ~3.5.
- **2× speed:** 290 words/min, still intelligible, MOS ~3.0.
- **0.75× speed:** 110 words/min, MOS ~3.7 (slightly clearer for hard words).

`pyttsx3` uses espeak — quality is decent but not neural. Modern neural TTS (Azure, Google) would be MOS ~4.5.

---

## 5.9 MOS evaluation

10 listeners, 20 samples:

| System | Mean MOS | 95% CI |
|---|:---:|:---:|
| Formant | 2.8 | [2.5, 3.1] |
| Concatenative | 4.0 | [3.7, 4.3] |
| Neural (Tacotron) | 4.6 | [4.3, 4.9] |
| Human | 4.9 | [4.7, 5.0] |

Neural **overlaps with human** in CI — statistically indistinguishable. Formant is clearly worse. Human still slightly ahead, but the gap is closing.

---

## 5.10 Speaker verification features

5 speakers, 10 recordings each, GMM-UBM (8 Gaussians):

- **Speaker accuracy:** 92%.
- **EER:** 6%.
- **Per-speaker variation:** speaker 3 (similar voice to speaker 1) has 80% accuracy; speakers 2, 4, 5 have ~95%.

**EER 6%** matches the classical literature for small-vocabulary speaker ID with short training.

---

## 5.11 DET curve

From the 5-speaker system:

- **EER:** 6% (from the DET curve's crossing point).
- **DET shape:** relatively straight in log-log — comparable to classical GMM-UBM (EER 5–10%). Modern x-vector would be 2–4%; ECAPA < 1%.
- **Human EER:** ~2–3% on typical speaker verification tasks — **modern neural systems exceed human accuracy**.

---

## 5.12 VoIP bandwidth and latency

For G.711, G.729, Opus@24kbps:

| Codec | Speech | Headers | Total | Header overhead |
|:---:|:---:|:---:|:---:|:---:|
| G.711 | 64 | 21.6 | 85.6 | 33% |
| G.729 | 8 | 21.6 | 29.6 | **270%** |
| Opus @ 24 | 24 | 21.6 | 45.6 | 90% |
| G.729 + cRTP | 8 | 1.6 | 9.6 | 17% |

**Latency** (60 ms network, 30 ms jitter, 25 ms codec):

Total one-way = 25 + 60 + 30 + 10 = **125 ms** — under the G.114 150 ms threshold. ✓

Adding 100 ms extra buffering → 225 ms — over the threshold, conversation feels lagged.

---

## 5.13 Anti-spoofing exploration

Test: MFCC-GMM speaker verifier against TTS sample generated in the target voice.

- **TTS accepted with probability 40%** — the GMM-UBM is fooled because the spectral envelope matches (TTS reproduces formants well).
- **Countermeasures that work:** (i) spectral phase features (TTS often has unnatural phase), (ii) pitch dynamics (TTS pitch is unnaturally smooth), (iii) liveness challenge (random phrase).
- **Modern TTS (neural)** is harder to distinguish — **anti-spoofing is now a field in itself** (ASVspoof challenges).