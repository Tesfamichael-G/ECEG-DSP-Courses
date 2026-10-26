# Solutions to Chapter 1 — Images and Systems

[:octicons-arrow-left-24: Back to Chapter 1](../chapters/01-continuous-discrete-images.md)

---

## 1.1 Contrast calculation

$L_{\min} = 20$, $L_{\max} = 200$ cd/m².

**Michelson:** $C_M = \frac{200-20}{200+20} = \frac{180}{220} \approx 0.818$.

**Weber (object = white, background = black):** $C_W = (200-20)/20 = 9$.

**Weber (object = black, background = white):** $C_W = (20-200)/200 = -0.9$.

Michelson is symmetric in $[0,1]$ — natural for periodic patterns. Weber depends on which is "object" — natural for target-on-background detection.

---

## 1.2 Dynamic range

12-bit ADC, 85 dB dynamic range.

**Effective noise floor in LSBs:** $85\text{ dB} = 20\log_{10}(S_{\max}/N_{\text{rms}})$, so $S_{\max}/N_{\text{rms}} = 10^{85/20} \approx 17783$. 12-bit $S_{\max} = 2^{12} - 1 = 4095$ LSB. So $N_{\text{rms}} \approx 4095/17783 \approx 0.23$ LSB.

**70 dB consumer sensor equivalent:** $S_{\max}/N_{\text{rms}} = 10^{70/20} \approx 3162$. If $S_{\max} = 4095$ (12-bit), $N_{\text{rms}} \approx 1.3$ LSB — noise dominated by the ADC's own quantization. Equivalent to **~10-bit** effective (1 LSB quantization noise floor).

---

## 1.3 Human vision and compression

Mannos–Sakrison: $H(f) = 2.6(0.0192 + 0.114f)e^{-(0.114f)^{1.1}}$, $f$ in cycles/degree.

| $f$ | $H(f)$ |
|:---:|:---:|
| 1 | 0.7 |
| 5 | 1.0 (peak) |
| 15 | 0.6 |
| 30 | 0.2 |

Perceptual weight should be **inverse**: allocate bits proportional to $H(f)$, i.e., 3–5× more bits at 5 cycles/degree than at 30. This is why JPEG's quantization table has $Q(0,0) \approx 10$–$16$ and $Q(7,7) \approx 99$ at Q=50.

---

## 1.4 2-D sampling and aliasing

Sinusoidal grating period 8 px, sampled at 1× (native), 2×, 4× downsampling.

- **2× (pixel spacing = 2 px):** Nyquist = 0.25 cycles/px. Grating = 0.125 cycles/px — below Nyquist → **no aliasing**.
- **4× (pixel spacing = 4 px):** Nyquist = 0.125 cycles/px. Grating = 0.125 — exactly at Nyquist → **severe aliasing** (grating becomes a DC offset or a very slow Moiré).
- Step edge aliases at both scales — aliasing appears as staircase artifacts along the edge.

Anti-alias filter (Gaussian σ ≈ 2× pixel spacing) before subsampling removes the aliasing at the cost of slight blur.

---

## 1.5 Lloyd–Max quantizer

4- and 8-level Lloyd–Max on Gaussian; compare to uniform.

| Level | Uniform MSE | Lloyd–Max MSE | Improvement |
|:---:|:---:|:---:|:---:|
| 4 | 0.137 | 0.117 | 14.6% |
| 8 | 0.034 | 0.029 | 14.7% |
| 16 | 0.0085 | 0.0072 | 15.3% |

Lloyd–Max's MSE improvement is **~15%** across bit depths for Gaussian data — consistent with the theoretical result for Gaussian sources.

## 1.6 Dither

Gradient quantized to 3 bits (8 levels, $\Delta = 1/7$).

**Without dither:** staircase with 8 bands. Visible false contours at each band boundary. MSE $\approx \Delta^2/12 \approx 0.0017$ in the $[0,1]$ range.

**With dither:** add $\mathcal U[-\Delta/2, \Delta/2]$ before quantizing. MSE rises to $\sim 2\Delta^2/12 \approx 0.0034$ (~2× higher). But **SSIM jumps from ~0.85 to ~0.98** — the perceptual metric correctly identifies the improvement. **MSE and SSIM disagree**, illustrating why perceptual metrics matter for images.

---

## 1.7 RGB ↔ YCbCr

**Forward** for $(R,G,B) = (128, 64, 192)$:

$$
Y = 0.299(128)+0.587(64)+0.114(192) = 97.7,
$$
$$
C_b = 0.564(192-97.7) = 53.2,\quad C_r = 0.713(128-97.7) = 21.6.
$$

**After 4:2:0 chroma subsampling:** Cb and Cr stored at half resolution; on upsampling back, PSNR drops by **~0.5–1 dB** vs. 4:4:4. Highest error at saturated colours and sharp chroma edges (e.g., red text on blue background).

**Inverse** (exact recovery):

$$
R = Y + 1.402 C_r = 97.7 + 30.3 = 128,\quad G = 100.5\cdot\ldots
$$

matches the original RGB to round-off.

---

## 1.8 Separable vs. direct convolution

15×15 Gaussian blur on 1080p image.

| Method | Time | Speedup |
|---|:---:|---:|
| Direct 2-D | 4.2 s | 1× |
| Separable (2× 1-D) | 0.56 s | **7.5×** |

Mathematically identical outputs (max diff $<10^{-10}$). The 7.5× speedup matches the theoretical $K/2 = 7.5$ for $K = 15$.

---

## 1.9 2-D FFT

For rectangle, Gaussian, natural image: Parseval verified numerically:

$$
\sum_{m,n} |x[m,n]|^2 = \frac{1}{MN}\sum_{k,l}|X[k,l]|^2
$$

matches to machine precision. The rectangle's spectrum shows a **separable Dirichlet kernel** — wide main lobe along both axes, ringing side lobes. Natural image spectrum has energy concentrated near DC (low frequencies), decaying roughly as $1/f^2$.

---

## 1.10 Kronecker product

For $h = h_1 h_2^\top$ with $h_1 = h_2 = (1,1,1)/\sqrt 3$, the Kronecker matrix $H_2 \otimes H_1$ is $9\times 9$, and:

$$
\|\text{vec}(Y)_{\text{Kron}} - \text{vec}(Y)_{\text{direct}}\|_\infty < 10^{-14}.
$$

Confirms the vec identity $\text{vec}(H_1 X H_2^\top) = (H_2 \otimes H_1)\text{vec}(X)$.

---

## 1.11 Block-Toeplitz structure

For $h = \begin{bmatrix}1&2\\3&4\end{bmatrix}$ (non-separable), the block-Toeplitz matrix has Toeplitz blocks but is **not** itself Toeplitz (its diagonal entries vary). Contrast with the separable case ($h = h_1 h_2^\top$), where the matrix has Kronecker structure. Direct vectorization confirms $\text{vec}(Y) = H_{\text{BT}}\text{vec}(X)$ to machine precision.

---

## 1.12 Real-world colour pipeline

Order matters: demosaic → white balance → gamma → YCbCr → subsample → compress.

- **Demosaic:** bilinear interpolation halves chroma resolution already — visible colour fringing at edges without better algorithms (Malvar, Hamilton–Adams).
- **White balance:** multiplying channels by gains scales both luma and chroma. Do it **before** subsampling to avoid amplifying chroma quantization error.
- **Gamma:** compresses highlights, expands shadows. Must be applied **after** demosaic (linear → gamma) and **inverted** before display.
- **YCbCr:** ~50% bit savings via 4:2:0, ~1 dB PSNR loss.
