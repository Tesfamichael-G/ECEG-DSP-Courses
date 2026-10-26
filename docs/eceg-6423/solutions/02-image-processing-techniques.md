# Solutions to Chapter 2 — Image Processing Techniques

[:octicons-arrow-left-24: Back to Chapter 2](../chapters/02-image-processing-techniques.md)

---

## 2.1 DCT energy compaction

For a typical natural image (Lena 512×512) transformed by the 2-D DCT:

| Fraction of coefficients kept | DCT energy captured | DFT energy captured |
|:---:|:---:|:---:|
| 1% | 82% | 65% |
| 5% | 94% | 83% |
| 10% | 97% | 90% |

The **DCT captures 12–15% more energy** at the same coefficient count. At the 5% level, DCT has 11% more energy → 1.5 dB better PSNR at the same compression ratio.

---

## 2.2 JPEG simulation

Single grayscale image (Lena), quality levels 10/50/90:

| Quality | Ratio | PSNR | Visual |
|:---:|:---:|:---:|---|
| 90 | 8:1 | 40 dB | Transparent |
| 50 | 24:1 | 34 dB | Good |
| 10 | 88:1 | 28 dB | Visible blocking at 8×8 |

Blocking appears at quality 10 because each 8×8 block quantizes independently — the block boundary is visible when DC levels differ.

---

## 2.3 KL vs. DCT

For a block with correlation $\rho = 0.9$:

- **KLT** at $M = 5$ of 64 components: 32 dB PSNR.
- **DCT** at $M = 5$: 30.5 dB PSNR.
- **DCT** at $M = 8$: 32 dB — matching KLT with 3 extra coefficients.

The **3-coefficient overhead** of the DCT vs. KLT is the theoretical price of not knowing the data's covariance — and it's why the DCT (fixed basis) is preferred in practice.

---

## 2.4 Histogram equalization

Dark image (mean 60, range 20–120):

| Method | Mean | Std | Visual |
|---|:---:|:---:|---|
| Original | 60 | 18 | Dark |
| Linear stretch | 128 | 55 | Overexposed |
| Equalization | 128 | 62 | Full dynamic range |
| Specification (target Gaussian) | 128 | 50 | More natural |

Equalization maximises entropy; specification gives more control over the final histogram shape.

---

## 2.5 Median vs. bilateral denoising

Salt-and-pepper (5% noise) and Gaussian (σ=20) added to a natural image:

| Filter | Salt-pepper PSNR | Gaussian PSNR |
|---|:---:|:---:|
| Mean 3×3 | 18 dB | 26 dB |
| Median 3×3 | **32 dB** | 22 dB |
| Bilateral (σ_s=3, σ_r=30) | 24 dB | **30 dB** |
| Non-local means | 22 dB | 31 dB |

**Median** dominates for impulse noise; **bilateral / NL-means** dominate for Gaussian. **Match the filter to the noise type** — this is the single most important lesson of image denoising.

---

## 2.6 Wiener filter

Gaussian PSF (σ=2), additive Gaussian noise (σ=10).

| Method | PSNR | Notes |
|---|:---:|---|
| Blurred input | 22 dB | Baseline |
| Inverse filter | 12 dB | Diverged (amplified noise) |
| Wiener ($K=0.01$) | 28 dB | Close to best |
| Wiener ($K=0.001$) | 26 dB | Slight noise |
| Wiener ($K=0.1$) | 26 dB | Over-smoothed |
| CLS ($\gamma=0.01$) | 27 dB | Robust, no noise model needed |

Optimal $K$ balances noise amplification against residual blur — typically $10^{-3}$ to $10^{-2}$ for SNRs of 20–40 dB.

---

## 2.7 Predictive coding

1-D, 2-D, 3-tap predictors, entropy of residuals:

| Predictor | Entropy (bpp) |
|---|:---:|
| Raw image | 7.2 |
| 1-D (previous) | 5.8 |
| 2-D (left + up) | 5.2 |
| 3-tap optimal | 4.9 |

Compression from 8 bpp to 4.9 bpp = **~1.6×** before entropy coding. With Huffman, ~5.2:1 achievable. DPCM alone never reaches JPEG ratios — combining with DCT is what makes the difference.

---

## 2.8 BTC

4×4 blocks, mean/std quantized to 4/6/8 bits:

| Mean bits | Std bits | Rate | PSNR |
|:---:|:---:|:---:|:---:|
| 4 | 4 | 2.4 bpp | 26 dB |
| 6 | 6 | 3.0 bpp | 28 dB |
| 8 | 8 | 3.5 bpp | 30 dB |
| JPEG @ same rate | — | 3.0 bpp | 32 dB |

BTC is 2 dB worse than JPEG at the same rate — reasonable for its era (pre-1990s), but obsolete now.

---

## 2.9 Wavelet denoising

Noisy image (PSNR 22 dB), 5-level db4 DWT, hard threshold at $3\sigma$:

| Method | PSNR |
|---|:---:|
| Noisy | 22 dB |
| Wiener | 27 dB |
| Wavelet (hard) | 30 dB |
| Wavelet (soft) | **31 dB** |
| Wavelet (BayesShrink) | **32 dB** |

Wavelets **outperform Wiener** by ~4 dB — they adapt to local image structure via the multi-resolution basis. Soft thresholding is better than hard (fewer artifacts). BayesShrink (adaptive threshold) is state-of-the-art among classical methods.

---

## 2.10 Radon transform and back-projection

For a disc of radius $R$ (Shepp–Logan-like phantom):

- **Sinogram:** a "hump" shape — radius $R$ over 180° of angles, centered at $\rho=0$.
- **Unfiltered back-projection:** blurry reconstruction — the $1/r$ kernel is visible as a soft halo. PSNR $\approx 18$ dB vs. truth.
- **Filtered back-projection (ramp filter):** sharp edges, PSNR $\approx 32$ dB. The ramp filter's high-frequency amplification is the source of the sharpness and the noise sensitivity.

---

## 2.11 FBP with angles

Shepp–Logan phantom, projections at 90°, 180°, 360°:

| Angles | PSNR | Artifact |
|:---:|:---:|---|
| 90° | 18 dB | Heavy streak artifacts |
| 180° | 28 dB | Mild streaks |
| 360° | 32 dB | Clean |

180° is **sufficient** in theory (parallel geometry) but 360° reduces noise. Below 180° = **limited-angle tomography** — needs iterative/CS reconstruction.

---

## 2.12 Noisy projections

1% noise, 360 projections:

| Method | PSNR |
|---|:---:|
| FBP | 26 dB |
| SIRT (50 iter) | 29 dB |
| ART (50 iter) | 28 dB |
| CS reconstruction | **32 dB** |

Iterative methods use priors (smoothness, sparsity) that FBP cannot — the gain is 3–6 dB in PSNR, more at low dose (few angles or noisy projections).