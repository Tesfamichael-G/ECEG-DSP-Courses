# Exam Bank — ECEG-6423 Multimedia Signal Processing

40 questions spanning all six chapters. Use for self-assessment and
exam preparation. **Answers and grading rubric at the end.**

**Instructions (for a mock exam):**

- Time: **3 hours**
- Format: **40 questions**, mixed MCQ, short answer, and derivation
- Closed book; **one A4 handwritten formula sheet** allowed
- Total marks: **100** (weighted by question)
- Marks per question in brackets

---

## Section A — Image Fundamentals (Ch 1) — 25 marks

### A1. Multiple choice [2]

The **Michelson contrast** of a pattern with luminances 20 and 120 cd/m² is:

(a) 0.71 (b) 0.83 (c) 5.0 (d) 0.17

### A2. Multiple choice [2]

The **contrast sensitivity function** of the human visual system peaks at
approximately:

(a) 0 cycles/degree (b) 4–8 cycles/degree (c) 40 cycles/degree (d) 200 cycles/degree

### A3. Multiple choice [2]

For an image sampled at pixel spacing $\Delta x$, the **Nyquist frequency** is:

(a) $\Delta x$ (b) $1/\Delta x$ (c) $1/(2\Delta x)$ (d) $2\Delta x$

### A4. Short answer [4]

State the **2-D sampling theorem** and give the condition for no aliasing
in both spatial directions.

### A5. Short answer [4]

Explain why the **YCbCr** colour space is used for compression instead of
RGB. Give two specific reasons.

### A6. Derivation [5]

For the uniform quantizer with step size $\Delta$, derive the mean-square
quantization error assuming the input signal is uniformly distributed
within each quantization bin. What is the resulting SNR in dB?

### A7. Short answer [3]

Define the **Kronecker product** $A \otimes B$. State the vec identity
$\text{vec}(AXB) = ?$ and explain why it is useful in image processing.

### A8. Short answer [3]

A smooth gradient is quantized to 4 levels without dither. What artifact
appears? Why does adding dither improve perceptual quality despite raising
the numerical MSE?

---

## Section B — Image Processing Techniques (Ch 2) — 25 marks

### B1. Multiple choice [2]

The **DCT** is preferred over the **DFT** for image compression because:

(a) It is faster (b) It has real coefficients and better energy compaction
(c) It requires no quantization (d) It preserves phase perfectly

### B2. Multiple choice [2]

The **KL transform** is theoretically optimal but rarely used because:

(a) It requires too much memory (b) It is data-dependent
(c) It produces complex coefficients (d) It is slow

### B3. Short answer [4]

Compare **median filtering** and **mean filtering** for removing salt-and-pepper
noise. Explain why one works better than the other.

### B4. Short answer [4]

Derive the **histogram equalization** mapping $s = (L-1) C(r)$ where
$C(r)$ is the input CDF. Show that the output has approximately uniform
histogram.

### B5. Derivation [5]

For an image degraded by $g = h * f + \nu$, derive the **Wiener filter**
$W(u, v)$ that minimizes $\mathbb E[(f - \hat f)^2]$. State the assumptions.

### B6. Short answer [4]

Describe the **JPEG compression pipeline** (8 steps from block splitting to
entropy coding). What are the three types of redundancy JPEG exploits?

### B7. Short answer [4]

State the **Radon transform** and the **Fourier slice theorem**. Use them
to derive the **filtered back-projection** algorithm.

---

## Section C — Speech Fundamentals (Ch 3) — 10 marks

### C1. Multiple choice [2]

The **fundamental frequency** $F_0$ for adult female speech is typically:

(a) 30–60 Hz (b) 80–150 Hz (c) 150–250 Hz (d) 500+ Hz

### C2. Short answer [4]

Describe the **speech production mechanism**. Name the four subsystems and
their acoustic roles.

### C3. Short answer [4]

Explain the **source-filter model**. Why does this model explain the
intelligibility of whispered speech?

---

## Section D — Speech Methods (Ch 4) — 20 marks

### D1. Multiple choice [2]

For voiced speech, the **short-time zero-crossing rate** is typically:

(a) Higher than 0.5 (b) Lower than 0.15 (c) Exactly 1 (d) Equal to the pitch

### D2. Short answer [3]

State the **time-frequency resolution trade-off** in the STFT. What window
length is used for **narrowband** vs. **wideband** spectrograms, and what
does each reveal?

### D3. Derivation [5]

State the **LPC model** of speech. Derive the **Yule–Walker equations** from
the minimum-MSE prediction criterion.

### D4. Short answer [4]

Describe the **cepstrum** $c[n] = \text{IFFT}\{\log|X[k]|\}$. Explain how it
separates source (pitch) from filter (formants). Why does the pitch peak
appear at high quefrency?

### D5. Short answer [3]

State the **Levinson–Durbin recursion**. What are the **reflection coefficients**,
and what do they mean physically?

### D6. Short answer [3]

Compare **three pitch-detection methods**: autocorrelation, cepstrum, and
LPC residual. Which is best for which situation?

---

## Section E — Speech Applications (Ch 5) — 10 marks

### E1. Multiple choice [2]

The three fundamental **HMM problems** are:

(a) Evaluation, decoding, training (b) Input, output, state
(c) Filter, predict, smooth (d) Fit, test, retrain

### E2. Short answer [3]

State the **speech-recognition problem** in terms of the Bayes factorization
$\Pr(W \mid X)$. What are the acoustic model and the language model, and
why is the factorization useful?

### E3. Short answer [3]

Describe the **DTW algorithm**. What is the recursion? What is the role
of the **Sakoe–Chiba band**?

### E4. Short answer [2]

Compare **speaker verification** and **speaker identification**. What are
their main features?

---

## Section F — Video Processing (Ch 6) — 10 marks

### F1. Multiple choice [2]

For a video sampled at frame rate $f_s$, the **wagon-wheel effect** occurs
when:

(a) $f_s$ is too high (b) $f_s$ is below 2× the wheel's rotation frequency
(c) The wheel is black (d) The wheel is transparent

### F2. Short answer [3]

State the **optical flow constraint equation** $I_x u + I_y v + I_t = 0$.
Explain the **aperture problem** and how **Lucas–Kanade** handles it.

### F3. Short answer [3]

Describe the **hybrid video codec**. What role does motion-compensated
prediction play?

### F4. Short answer [2]

Compare **H.264**, **HEVC**, and **AV1** in terms of compression efficiency
and typical application. Which would you choose for a 4K streaming service
and why?

---

## Answer key

??? success "Click to reveal answers"

    ### Section A — Image Fundamentals

    **A1.** (b) 0.83 — $C_M = (120-20)/(120+20) = 100/140 \approx 0.714$, wait: 100/140 = 0.714. Hmm, let me recompute: (120-20)/(120+20) = 100/140 = 0.714. Option (a) is 0.71. So **A1 = (a) 0.71**.

    **A2.** (b) 4–8 cycles/degree.

    **A3.** (c) $1/(2\Delta x)$.

    **A4.** 2-D sampling theorem: a band-limited image with spectrum zero outside $[-F_x, F_x] \times [-F_y, F_y]$ can be exactly reconstructed from samples on a lattice with spacing $\Delta x, \Delta y$ iff $1/\Delta x > 2F_x$ and $1/\Delta y > 2F_y$. Equivalently, sampling frequencies exceed twice the highest spatial frequencies in both directions.

    **A5.** Two reasons: (i) YCbCr separates **luminance** (perceptually important, high spatial resolution) from **chroma** (perceptually less important, can be subsampled); (ii) human vision has much lower spatial resolution for chroma than for luminance, so 4:2:0 subsampling (halving chroma resolution in both directions) reduces data by 50% with almost no perceptual loss.

    **A6.** For uniform quantizer with step $\Delta = 2A/L$, the error $e = Q(f) - f$ is approximately uniform on $[-\Delta/2, \Delta/2]$ if the signal is dense within each bin. Variance: $\sigma_e^2 = \frac{1}{\Delta}\int_{-\Delta/2}^{\Delta/2} e^2\, de = \Delta^2/12$. If signal has variance $\sigma_f^2$ and peak-to-RMS ratio $A/\sigma_f$, then $L = A/(\sqrt{12}\sigma_e)$, so $\text{SNR} = \sigma_f^2/\sigma_e^2 = 12 \sigma_f^2/\Delta^2 = \frac{L^2\sigma_f^2}{A^2} \cdot 3 = 3\cdot L^2/(A^2/\sigma_f^2)$. In dB: $\text{SNR}_{dB} \approx 6 b + 4.77 - 20\log_{10}(A/\sigma_f)$ where $b$ = bit depth.

    **A7.** Kronecker product $A \otimes B$ is the block matrix with blocks $a_{ij} B$; result size $mp \times nq$ for $A$ $m \times n$ and $B$ $p \times q$. Vec identity: $\text{vec}(AXB) = (B^\top \otimes A)\text{vec}(X)$. Useful because it converts 2-D operations into 1-D matrix-vector products, enabling fast algorithms (separable filtering, 2-D FFT) and analysis (eigenvalues of $A \otimes B$ are products).

    **A8.** Without dither: **false contours** (visible stair-step bands) at quantization boundaries. Adding dither converts the structured quantization error into random noise that the eye averages out — perceptually smoother even though numerical MSE is higher.

    ### Section B — Image Processing Techniques

    **B1.** (b) Real coefficients and better energy compaction.

    **B2.** (b) Data-dependent (basis depends on the covariance, must be transmitted).

    **B3.** **Median filtering:** an isolated outlier (salt or pepper pixel) is not the median of its neighborhood, so it is removed. **Mean filtering:** the outlier contaminates all pixels in the neighborhood (spreads the noise). Median **preserves edges** (median of a step neighborhood lies on the step); mean blurs edges. **Median is preferred for salt-and-pepper; mean is preferred for Gaussian noise.**

    **B4.** Continuous: input $r$ has PDF $p_r(r)$ and CDF $C(r) = \int_0^r p_r(s)ds$. Output $s = (L-1)C(r)$. Then $p_s(s) = p_r(r)/|ds/dr| = p_r(r)/((L-1)p_r(r)) = 1/(L-1)$ on $[0, L-1]$ — uniform. Discrete version uses empirical CDF, giving approximate uniformity.

    **B5.** **Wiener derivation:** minimize $\mathbb E[|F - W G|^2]$ over linear filters $W$ where $G = HF + N$. Optimal $W = \mathbb E[FG^*]/\mathbb E[|G|^2] = H^* S_f/(|H|^2 S_f + S_\nu)$. Assumptions: (i) $f, \nu$ zero-mean, uncorrelated; (ii) $S_f, S_\nu$ known; (iii) LSI filter; (iv) linear minimum-MSE criterion.

    **B6.** **JPEG pipeline:** (1) split into 8×8 blocks; (2) RGB→YCbCr, subsample chroma 4:2:0; (3) 8×8 DCT per block; (4) quantize with a table scaled by quality factor; (5) zigzag scan; (6) Huffman entropy coding; (7) (decoder side) inverse steps. Three redundancies exploited: **coding** (Huffman), **interpixel** (DCT), **psychovisual** (quantization table).

    **B7.** **Radon transform:** $R(\rho, \theta) = \int\int f(x,y)\delta(x\cos\theta + y\sin\theta - \rho)dxdy$ — line integrals at angle $\theta$. **Fourier slice theorem:** $\widehat R(\omega, \theta) = F(\omega\cos\theta, \omega\sin\theta)$ — 1-D FT of projection = 2-D FT of $f$ along the slice. **FBP:** (i) FT each projection; (ii) multiply by ramp $|\omega|$; (iii) IFT; (iv) back-project over all angles. The ramp corrects the $1/r$ blur of plain back-projection.

    ### Section C — Speech Fundamentals

    **C1.** (c) 150–250 Hz.

    **C2.** Four subsystems: **lungs** (power source — airflow), **larynx / vocal folds** (source — quasi-periodic pulse train for voiced, noise for unvoiced), **vocal tract** (filter — formants from acoustic resonances), **articulators** (tongue, lips, velum — configure the tract shape). Together they produce a source-filter cascade that generates the speech waveform.

    **C3.** The speech spectrum is $|X(f)| = |G(f)| \cdot |H(f)|$ — glottal source times vocal-tract filter. **Whispered speech:** the source is non-periodic (vocal folds open), so pitch is absent. But the vocal tract is unchanged, so the **formants are identical to normal speech** — the filter alone carries the linguistic information. This is why whispered speech is fully intelligible.

    ### Section D — Speech Methods

    **D1.** (b) Lower than 0.15.

    **D2.** **Trade-off:** longer window → better frequency resolution, worse time resolution (Heisenberg: $\Delta t \cdot \Delta f \ge 1/(4\pi)$). **Narrowband** (40 ms): resolves pitch harmonics and shows them as fine horizontal bands; good for pitch/formant analysis. **Wideband** (10 ms): resolves pitch pulses as vertical stripes; good for voicing and timing analysis.

    **D3.** **LPC model:** $x[n] = -\sum_{k=1}^p a_k x[n-k] + e[n]$. Minimize $\mathbb E[e[n]^2]$: $\partial/\partial a_k = 2\sum e[n]x[n-k] = 0$. Substituting $e[n] = x[n] + \sum_l a_l x[n-l]$: $\sum_l a_l R[k-l] = -R[k]$ for $k = 1, \ldots, p$. Matrix: $\mathbf R\mathbf a = -\mathbf r$ with $\mathbf R$ Toeplitz.

    **D4.** **Cepstrum:** $c[n] = \text{IFFT}\{\log|X[k]|\}$. In log-spectrum, source (pitch) and filter (formants) are additive. In cepstral domain, they occupy different quefrency ranges — filter (smooth envelope) at **low quefrency** near 0; pitch (harmonic comb) at **high quefrency** at $n = T_0$ (pitch period in samples). **Pitch at high quefrency** because the IFFT of a comb with spacing $F_0$ is a peak at period $T_0$.

    **D5.** **Levinson–Durbin:** recursion from order $m-1$ to $m$: $k_m = -(R[m] + \sum_{i=1}^{m-1}a_i^{(m-1)}R[m-i])/\sigma_{m-1}^2$; $a_i^{(m)} = a_i^{(m-1)} + k_m a_{m-i}^{(m-1)}$; $\sigma_m^2 = (1 - k_m^2)\sigma_{m-1}^2$. **Reflection coefficients:** $k_m$ = partial autocorrelation at lag $m$ — correlation between $x[n]$ and $x[n-m]$ after removing intermediate lags. $|k_m| < 1$ is the necessary and sufficient condition for AR stability.

    **D6.** **Autocorrelation:** cheap, works for most; prone to **octave errors**. **Cepstrum:** robust to formants, handles weak fundamentals; requires long frames. **LPC residual (SIFT):** best accuracy — the LPC filter whitens residual, leaving a sharp pitch peak. **SIFT for coding; cepstrum for analysis; ACF for embedded/cheap.**

    ### Section E — Speech Applications

    **E1.** (a) Evaluation, decoding, training.

    **E2.** $\hat W = \arg\max_W \Pr(W \mid X) \propto \arg\max_W \Pr(X \mid W)\Pr(W)$. **Acoustic model:** $\Pr(X \mid W)$ — likelihood of audio given word sequence. **Language model:** $\Pr(W)$ — prior probability of word sequence. **Useful:** enables independent training on different data (acoustic on transcribed speech, language on large text corpora), and lets context resolve acoustic ambiguity.

    **E3.** **DTW:** optimal monotone path aligning two sequences. Recursion: $C(i,j) = d(x_i, y_j) + \min(C(i-1,j), C(i,j-1), C(i-1,j-1))$. **Sakoe–Chiba band** restricts path to within $r$ of the diagonal — reduces cost from $O(MN)$ to $O(Mr)$, regularizes to prevent pathological warping.

    **E4.** **Verification:** 1:1 — "Is this the claimed speaker?" (authentication). **Identification:** 1:N — "Which enrolled speaker is this?" (forensics, personalization). Features: MFCCs, x-vectors, pitch statistics, speaking rate, idiolect.

    ### Section F — Video Processing

    **F1.** (b) $f_s$ below 2× the wheel's rotation frequency.

    **F2.** $I_x u + I_y v + I_t = 0$ from brightness constancy + first-order Taylor. **Aperture:** one equation, two unknowns — motion underdetermined. **LK solution:** assume constant $u, v$ over a window, giving $N$ equations, then least-squares: $\mathbf u = (A^\top A)^{-1}A^\top\mathbf b$. Aperture resolved by having enough independent constraints in the window; fails at flat regions (both $I_x, I_y$ small).

    **F3.** Pipeline: I/P/B classification → motion estimation → motion-compensated prediction → residual → DCT → quantize → entropy-code. **Motion-compensated prediction:** removes temporal redundancy — the residual (aligned difference between current and predicted frame) is much sparser than the frame itself, and the DCT compacts it further. This is the main source of compression in inter-frame coding.

    **F4.** **Efficiency (typical 1080p at PSNR 35 dB):** H.264 = 6 Mbps, HEVC = 3.8 Mbps (−38%), AV1 = 3.4 Mbps (−43%). **For 4K streaming:** **HEVC** if HEVC hardware decode is universal on the target devices (most modern TVs, phones, browsers); **AV1** if compression is more important than decode compatibility (Netflix, YouTube use AV1; older devices need fallback). **Recommendation: AV1 with H.264 fallback** for maximum efficiency; **HEVC only** if you need universal hardware decode without fallback.

---

## Grading rubric (for instructors)

| Band | Mark | Criteria |
|---|---|---|
| **A** | 90–100 | All questions correct or with minor arithmetic slips; derivations rigorous; interpretation insightful |
| **B** | 75–89 | Most questions correct; occasional conceptual gaps; derivations mostly sound |
| **C** | 60–74 | Basic understanding demonstrated; some errors in derivations or interpretation |
| **F** | < 60 | Significant conceptual gaps; major derivations missing or wrong |

**Point allocations:**

- MCQs: all-or-nothing.
- Short answers: partial credit for correct key steps.
- Derivations: partial credit per correct line.

**Common errors to watch for (and their marks):**

- Confusing **necessary** and **sufficient** conditions (deduct 1–2).
- **Wrong sign** in a Lagrange/KKT multiplier (deduct 1).
- **Missing or incorrect assumption** in a derivation (deduct 1–2).
- **Numerical error** in a final answer without correct method (deduct 2–3).
- **Correct answer with no working** (deduct full marks).
- **Copying a formula** without stating when it applies (deduct 1–2).