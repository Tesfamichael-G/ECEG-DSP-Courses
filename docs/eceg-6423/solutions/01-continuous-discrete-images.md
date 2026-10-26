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