# Chapter 1 — Continuous and Discrete Images and Systems

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 2 tutorial + 1 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** light, luminance, brightness, contrast · the eye ·
    the monochrome vision model · image processing problems and applications ·
    vision camera · digital processing system · 2-D sampling theory · aliasing ·
    image quantization · Lloyd–Max quantizer · dither · color images · linear
    systems and shift invariance · Fourier transform, Z-transform · matrix
    theory results · block matrices and Kronecker products.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Distinguish **radiometric** from **photometric** quantities; compute
   **Weber** and **Michelson** contrasts for simple patterns.
2. Describe the **human visual system** and state the **contrast sensitivity
   function** — the physiological basis for perceptual compression.
3. Model a **digital camera** pipeline from lens to ADC, and recognize the
   trade-offs between sensor size, dynamic range, and noise.
4. State the **2-D sampling theorem** and derive the **anti-aliasing** condition
   for images. Explain Moiré patterns.
5. Design **uniform** and **Lloyd–Max** quantizers; compute their MSE;
   understand the role of **dither** in visual quantization.
6. Convert between **RGB**, **YCbCr**, and **HSV** colour spaces; explain why
   colour compression uses luma-chroma separation.
7. Model a 2-D **linear shift-invariant (LSI)** system as a **convolution**;
   compute 2-D Fourier and Z transforms of simple images.
8. Use **Kronecker products** and **block matrices** to represent 2-D linear
   filtering as 1-D matrix algebra — the key algebraic tool for image coding
   and filtering.
9. Implement all of the above in **Python, Java, and C#**.

---

## 1.1 Light, luminance, brightness, contrast

### 1.1.1 Radiometry vs. photometry

**Radiometry** measures light in physical units of energy. **Photometry**
weights light by the sensitivity of the human eye.

| Quantity | Unit | Weighted by | Meaning |
|---|---|---|---|
| **Radiant flux** (power) | W (watt) | — | Physical energy per second |
| **Luminous flux** | lm (lumen) | V(λ) | Perceptually weighted power |
| **Radiance** | W/(sr·m²) | — | Power per unit area per solid angle |
| **Luminance** $L$ | cd/m² | V(λ) | Perceived brightness (photometric) |
| **Illuminance** | lux | V(λ) | Light falling on a surface |
| **Reflectance** | — | — | Ratio of reflected to incident |

The **luminous efficiency function** $V(\lambda)$ peaks at $\lambda \approx 555$ nm
(green-yellow) and falls to near zero at the UV (400 nm) and IR (700 nm)
extremes — human vision occupies only a small portion of the EM spectrum.

!!! note "Brightness"
    *Luminance* is an objective, measurable quantity (cd/m²). **Brightness** is
    the subjective perception of luminance, and is **nonlinearly** related to it
    — roughly by a power law (Stevens' law: $B \propto L^{0.33}$). This is why
    images are gamma-corrected before display: the file stores an approximately
    perceptual encoding of luminance, not a linear one.

### 1.1.2 Contrast

Two common definitions of contrast for a region of interest:

**Weber contrast** — relative to the background:

$$
\boxed{\;C_W = \frac{L_{\text{obj}} - L_{\text{bg}}}{L_{\text{bg}}}\;}
$$

**Michelson contrast** — relative to the average:

$$
\boxed{\;C_M = \frac{L_{\text{max}} - L_{\text{min}}}{L_{\text{max}} + L_{\text{min}}}\;}
$$

Michelson contrast lies in $[0, 1]$ and is symmetric in the two extremes;
Weber contrast can be any value, and is meaningful only when the background is
uniform and the object is small.

!!! tip "Why two definitions?"
    Michelson is used for **periodic patterns** (gratings, sinusoids) — the
    contrast of a sinusoid $L + A\sin$ is $A/L$ (Weber) or $A/L$ (Michelson,
    for zero-mean patterns).

    Weber is used for **targets on backgrounds** — e.g., a defect on a uniform
    surface, a letter on a page, a lesion on a medical image.

### Worked example 1.1 — contrast of a checkerboard

A 2-D checkerboard alternates between black squares of luminance $L_{\min} = 5$
cd/m² and white squares of luminance $L_{\max} = 100$ cd/m².

**Michelson contrast.**

$$
C_M = \frac{100 - 5}{100 + 5} = \frac{95}{105} \approx 0.905.
$$

**Weber contrast** (relative to black background, object = white square):

$$
C_W = \frac{100 - 5}{5} = 19.
$$

**Weber contrast** (relative to white background, object = black square):

$$
C_W = \frac{5 - 100}{100} = -0.95.
$$

Michelson is the natural measure for the whole pattern; Weber depends on which
square is "the object". Both give useful information — Michelson for a whole
image, Weber for a target detection task.

### Worked example 1.2 — dynamic range and display

A display can show luminances from $L_{\min} = 0.5$ cd/m² to $L_{\max} = 500$ cd/m².

**Dynamic range in stops** (log₂ of ratio):

$$
\text{DR} = \log_2\frac{L_{\max}}{L_{\min}} = \log_2(1000) \approx 9.97\ \text{stops}.
$$

**In decibels:** $\text{DR} = 20\log_{10}(L_{\max}/L_{\min}) \approx 60$ dB.

**8-bit integer representation** gives 256 levels — a dynamic range of only
$\sim$48 dB, less than the display. This is the reason **HDR** displays
(10–12 bits) and **gamma correction** exist: to map a wider range of physical
luminance into a narrower perceptual scale without banding.

---

## 1.2 The eye and the monochrome vision model

### 1.2.1 The human eye

The retina has two types of photoreceptors:

- **Rods** (~$10^8$) — sensitive to low light, no colour discrimination, peak
  sensitivity at $\lambda \approx 500$ nm (green-blue). Dominate in **scotopic**
  (dark) vision.
- **Cones** (~$6 \times 10^6$) — need bright light, three types (S/M/L — short /
  medium / long wavelength), peak sensitivity at 420 / 534 / 564 nm. Dominate
  in **photopic** (bright) vision.

**Distribution.** Cones are concentrated in the **fovea** (centre of gaze),
which subtends $\approx 2°$ of visual angle; rods are more uniformly distributed.
This is why we see fine detail only in the centre of gaze, and why peripheral
vision is primarily for motion and low light.

**Optical resolution.** The fovea can resolve $\approx 1$ arcminute of visual
angle — about $1/60°$. At a typical reading distance (40 cm), this corresponds
to $\approx 0.1$ mm spatial resolution, or $\approx 250$ pixels per inch.

### 1.2.2 The contrast sensitivity function (CSF)

The **CSF** describes the inverse of the minimum detectable contrast as a
function of **spatial frequency** (in cycles per degree of visual angle):

$$
\boxed{\;\text{CSF}(f) = \frac{1}{C_{\min}(f)}\;}
$$

**Shape.** Band-pass — peaks at $\sim$ 4–8 cycles/degree, falls off at both
low frequencies (lateral inhibition flattens large-scale gradients) and high
frequencies (optical and neural limits).

**Weber's law** for the low-frequency limit:

$$
\boxed{\;\frac{\Delta L}{L} = \text{constant} \approx 0.01 \ (1\%)\;}
$$

A 1% luminance change is just barely detectable under normal viewing conditions.

**Implication for compression.** The CSF tells you **which frequencies you can
sacrifice**. JPEG and JPEG 2000 quantize their DCT/wavelet coefficients more
aggressively at high frequencies because the eye is less sensitive there —
this is the physiological basis of perceptual coding.

### Worked example 1.3 — perceptual quantization mask

The CSF peaks at ~8 cycles/degree. For an image displayed at 30 pixels per
degree (typical mobile viewing), this corresponds to a **spatial frequency** of

$$
f_{\text{peak}} = 8\ \text{cycles/degree} / 30\ \text{pixels/degree} = 0.267\ \text{cycles/pixel}.
$$

Above ~$0.35$ cycles/pixel (corresponding to ~10.5 cycles/degree, past the
high-frequency cutoff), the eye's sensitivity drops sharply — the *human visual
cutoff* is at $\approx 40$–$60$ cycles/degree, corresponding to roughly
$0.25$–$0.35$ cycles/pixel at this viewing distance.

**Compression insight.** An image can be **downsampled** to ~2 pixels per cycle
of the highest visible frequency — i.e., to $\sim$ 60–70 pixels/degree of
viewing — without perceptual loss. This is why JPEG's "quality" setting maps to
a quantization matrix that coarsens at high frequencies, and why JPEG 2000's
wavelet decomposition can be truncated in the finest subbands with minimal
perceptual cost.

### 1.2.3 The monochrome (achromatic) vision model

For grayscale images, the HVS is well modelled as a **linear system** with a
**contrast-sensitivity transfer function** $H_{\text{HVS}}(f_x, f_y)$ that
shapes the input image's spectrum. Its magnitude is the CSF.

**Practical approximations.**

- **Mannos–Sakrison (1974):** a simple exponential model
  $H_{\text{CSF}}(f) = 2.6(0.0192 + 0.114f)e^{-(0.114f)^{1.1}}$
  where $f$ is in cycles/degree.
- **Daly (1993):** the widely used **Daly visual model** for full-image HVS
  simulation.
- **Watson (1993):** the discrete cosine transform domain model, the basis of
  the **DCTune** perceptual JPEG quantizer.

The full HVS also includes **lateral inhibition** (which sharpens edges), **luminance
adaptation** (Weber's law), and **temporal response** (flicker sensitivity). For
image compression, the CSF alone captures most of the perceptual redundancy.

---

## 1.3 Image processing problems and applications

### 1.3.1 Taxonomy of image processing

Image processing is conventionally divided into four levels:

| Level | Input | Output | Examples |
|---|---|---|---|
| **Acquisition / pre-processing** | Sensor signal | Image | Sampling, quantization, denoising, colour correction |
| **Enhancement** | Image | Image | Contrast stretching, histogram equalization, sharpening |
| **Restoration** | Degraded image | Restored image | Deconvolution, denoising, inpainting |
| **Compression** | Image | Bits | JPEG, JPEG 2000, predictive, VQ, wavelet |
| **Analysis / computer vision** | Image | Features / decisions | Segmentation, edge detection, object detection, recognition |

This course covers the first four levels (Ch 1–2) for images, then moves to
speech (Ch 3–5) and video (Ch 6).

### 1.3.2 Applications

**Medical imaging.** X-ray, CT, MRI, ultrasound, PET — image restoration,
segmentation, computer-aided diagnosis. CT and PET reconstruction from
projections (Radon transform) is the direct application of §2.4.4.

**Remote sensing.** Satellite and aerial imagery — pansharpening, change
detection, land-cover classification, crop health monitoring. Requires
geometric and radiometric correction, plus compression for storage and downlink.

**Industrial inspection.** Machine vision for defect detection, quality control,
robot guidance. High-speed processing, high accuracy, controlled lighting.

**Document processing.** OCR, barcode reading, archival digitization — deskewing,
binarization, denoising, text enhancement.

**Forensics / security.** Fingerprint enhancement, face recognition, license
plate reading, video surveillance — image restoration and feature extraction.

**Entertainment / multimedia.** Photography, film, gaming, broadcast — colour
science, HDR, denoising, compression, style transfer.

**Scientific imaging.** Microscopy, astronomy, particle physics — deconvolution,
super-resolution, denoising, particle tracking.

!!! tip "The recurring theme"
    Most of these applications involve the **same core operations** — sampling,
    quantization, transforms, filtering, compression — applied to **different
    data**. The mathematical core (Ch 1) is shared across all of them.

---

## 1.4 Vision camera and digital processing system

### 1.4.1 The sensor pipeline

A digital camera pipeline converts **photons** to **bits** through a sequence
of stages:

1. **Optics** — lens focuses light onto the sensor. Aperture controls
   illuminance; diffraction limits resolution.
2. **Sensor** — CCD or CMOS converts photons to charge / voltage. Each pixel
   integrates over a small area (typically $2$–$7\ \mu$m). Sensor size
   (e.g., 1" ≈ 16 mm diagonal for a "one-inch" sensor) and pixel size determine
   sensitivity, dynamic range, and noise.
3. **Analog front end** — amplifier + **anti-alias filter** + **ADC** to
   digital samples. ADC resolution (8, 10, 12, 14 bits) sets dynamic range.
4. **Digital processing** — demosaicing, white balance, gamma correction,
   denoising, colour correction, compression.

### 1.4.2 Colour filter arrays — the Bayer pattern

Most colour sensors use a **colour filter array (CFA)** — a mosaic of red, green,
and blue filters over the pixel array. The **Bayer pattern** is the standard:

```
R G R G R G ...
G B G B G B ...
R G R G R G ...
```

with 2× green pixels per R or B (matching the eye's higher sensitivity to
green). **Demosaicing** reconstructs full RGB at each pixel from this mosaic.

### 1.4.3 Dynamic range and noise

The **dynamic range** of a sensor is the ratio of its maximum (saturation)
signal to its noise floor:

$$
\text{DR} = \frac{S_{\max}}{N_{\text{rms}}}.
$$

Typical values:

| Sensor class | Bit depth | Dynamic range |
|---|:---:|:---:|
| Consumer compact | 8–10 | 40–50 dB |
| Consumer DSLR / mirrorless | 12–14 | 70–85 dB |
| Scientific / high-end | 16–18 | 90–110 dB |

**Noise sources.** Photon shot noise (proportional to $\sqrt{N}$ photons, the
fundamental limit), readout noise (constant per pixel), fixed-pattern noise
(pixel-to-pixel variations, correctable), dark current (temperature-dependent,
correctable), quantization noise (from the ADC).

**Anti-alias filter.** A slight optical blur is introduced **before** the sensor
to suppress frequencies above the sensor's Nyquist. This reduces Moiré but
softens the image slightly — the constant trade-off between resolution and
aliasing.

!!! tip "The Nyquist rule for sensors"
    A camera with $W$ horizontal pixels can resolve at most $W/2$ cycles across
    the sensor — the **Nyquist frequency** of the sensor. The lens and the
    optical anti-alias filter must attenuate frequencies above this to avoid
    aliasing; otherwise checkerboard and line-pattern objects will show
    visible Moiré artifacts. Modern high-resolution sensors often omit the
    anti-alias filter entirely, relying on the lens to provide enough attenuation.

---

## 1.5 2-D sampling theory

### 1.5.1 Sampling an image

A **continuous image** $f(x, y)$ on $\mathbb R^2$ is sampled on a lattice

$$
x_m = m \Delta x,\qquad y_n = n \Delta y,\qquad m, n \in \mathbb Z.
$$

The **sampled image** is

$$
\boxed{\;f_s[m,n] = f(m \Delta x, n \Delta y).\;}
$$

The **sampling frequencies** are $f_{s,x} = 1/\Delta x$ and $f_{s,y} = 1/\Delta y$
(cycles per unit length).

### 1.5.2 The 2-D sampling theorem

!!! abstract "2-D sampling theorem"
    Let $f(x,y)$ be a band-limited image with **2-D Fourier transform**
    $F(f_x, f_y)$ vanishing outside $[-F_x, F_x] \times [-F_y, F_y]$.
    The image is **exactly recoverable** from its samples iff

    $$
    \boxed{\;\frac{1}{\Delta x} > 2F_x,\qquad \frac{1}{\Delta y} > 2F_y.\;}
    $$

    i.e., the sampling frequencies exceed twice the highest spatial frequencies
    of the image.

The recovery is by **ideal 2-D sinc interpolation**:

$$
\boxed{\;f(x,y) = \sum_{m,n} f_s[m,n]\,\operatorname{sinc}\!\left(\frac{x - m\Delta x}{\Delta x}\right)\operatorname{sinc}\!\left(\frac{y - n\Delta y}{\Delta y}\right)\;}
$$

where $\operatorname{sinc}(u) = \sin(\pi u)/(\pi u)$.

### 1.5.3 Frequency-domain sampling

The 2-D spectrum of the sampled image is the **periodic extension** of the
original spectrum:

$$
\boxed{\;F_s(f_x, f_y) = \frac{1}{\Delta x \Delta y}\sum_{k,l} F\!\left(f_x - \frac{k}{\Delta x},\ f_y - \frac{l}{\Delta y}\right)\;}
$$

The original spectrum is replicated at intervals of $(1/\Delta x, 1/\Delta y)$
in the frequency plane. **Aliasing** occurs when these replicas overlap.

The condition for no aliasing is the **Nyquist condition** of §1.5.2. If it is
violated, high-frequency content **folds** into lower frequencies — the
high-frequency spatial features of the image become **aliased** as false
low-frequency patterns.

### 1.5.4 Aliasing in images

Aliasing in 2-D images has two signature forms:

**Moiré patterns.** A fine periodic pattern (e.g., a picket fence, a striped
shirt, a screen) sampled at a frequency below its own fundamental produces a
**beat pattern** at the difference frequency. The resulting Moiré has a
period that is not present in the original — a **spurious low-frequency
structure**.

**Checkerboard aliasing.** A checkerboard at the pixel-sampling scale (one
pixel on, one off) can alias to a *different* checkerboard at a different
scale, or to a uniform gray field, depending on the exact sampling offset.

### Worked example 1.4 — 2-D sampling of a sinusoid

Consider a continuous sinusoid

$$
f(x, y) = \cos(2\pi(f_x x + f_y y)),\qquad f_x = f_y = 0.6\ \text{cycles/unit}.
$$

Sample at spacing $\Delta x = \Delta y = 1$ unit ⇒ sampling frequency
$f_s = 1$ cycles/unit. **Nyquist frequency is $0.5$ cycles/unit — the signal
exceeds Nyquist.**

**Aliased frequency.** The sampled sinusoid is
$$
f_s[m, n] = \cos(2\pi \cdot 0.6 \cdot (m + n)) = \cos(2\pi (0.6m + 0.6n)).
$$

Modulo 1 (unit sampling):
$$
0.6 \equiv -0.4 \pmod 1 \quad\text{and}\quad 0.6 \equiv -0.4 \pmod 1.
$$

So the aliased frequency is at $-0.4$ cycles/unit — **below Nyquist**, so the
sampled image shows a sinusoid at $0.4$ cycles/unit instead of the true $0.6$.

**Direction.** The original was along $(1,1)$; the alias is along $(1,1)$ as
well, but with the opposite apparent propagation direction (negative frequency
= reversed phase). The visual effect: a slow sinusoid in the same direction as
the original — a Moiré pattern.

**Anti-aliasing.** Sample at $\Delta x = 0.8$ unit instead ($f_s = 1.25$,
Nyquist = $0.625$): now $0.6 < 0.625$, so the sinusoid is below Nyquist and
is correctly captured. Or, keep $\Delta x = 1$ but first lowpass the signal
with cutoff $< 0.5$ — this attenuates the $0.6$-frequency content, but avoids
the false alias.

### Worked example 1.5 — subsampling a checkerboard

A continuous checkerboard with period $T = 4$ pixels (i.e., $f = 0.25$
cycles/pixel) is sampled at $1$ pixel spacing ($f_s = 1$, Nyquist = $0.5$).
The pattern is below Nyquist — no aliasing.

Now subsample the **sampled** image by $M = 3$: new sampling frequency
$f_s' = 1/3$, new Nyquist = $1/6 \approx 0.167$ cycles/pixel.

The pattern's apparent frequency (after the original sampling) is $0.25$
cycles/pixel — above the new Nyquist $0.167$. Aliased frequency:

$$
0.25 \pmod{1/3} = 0.25 - 0.333 \cdot 0 \ldots
$$

The alias lives at $0.25$ cycles/pixel, which is above the new Nyquist but the
sampled representation sees it at
$$
0.25 - 1/3 \approx -0.083 \Rightarrow \text{frequency } 0.083 \text{ cycles/pixel (alias)}.
$$

So the subsampled checkerboard now appears to have period $12$ pixels
($1/0.083$) — a classic **Moiré beat** between the original $4$-pixel period
and the $3$-pixel sampling period. This is what you see when a fine grid is
photographed at low resolution.

---

## 1.6 Image quantization

### 1.6.1 Uniform quantization

A continuous-amplitude sample $f \in [-A, A]$ is quantized to one of $L = 2^b$
levels (where $b$ is the bit depth) with step size

$$
\boxed{\;\Delta = \frac{2A}{L}.\;}
$$

The **uniform quantizer** is

$$
\boxed{\;Q(f) = \Delta\left\lfloor \frac{f + A}{\Delta} \right\rfloor + \tfrac{\Delta}{2}\;}
$$

The **quantization error** $e = Q(f) - f$ is approximately uniform on
$[-\Delta/2, \Delta/2]$, with variance

$$
\boxed{\;\sigma_e^2 \approx \frac{\Delta^2}{12}.\;}
$$

This gives the standard **6 dB per bit** rule:

$$
\boxed{\;\text{SNR}_{Q} = 10\log_{10}\!\frac{\sigma_f^2}{\sigma_e^2} \approx 6b + 4.77 - 20\log_{10}\!\frac{f_{\max}}{\sigma_f}\ \text{dB}\;}
$$

where $f_{\max} / \sigma_f$ is the peak-to-RMS ratio of the signal (often
$\sim$4 for images).

**For images** with $b$-bit uniform quantization and full-scale Gaussian-like
statistics, $\text{SNR}_Q \approx 6b - 1$ dB. An 8-bit image has ~47 dB of
quantization SNR — usually more than enough (well above the ~30 dB needed for
perceptual transparency).

### 1.6.2 Lloyd–Max quantizer

For a **known distribution** $p(f)$ of the signal, the **Lloyd–Max quantizer**
minimizes the MSE

$$
D = \sum_{i=1}^{L} \int_{t_{i-1}}^{t_i} (f - y_i)^2 p(f)\, df
$$

over **reproduction levels** $\{y_i\}$ and **decision thresholds** $\{t_i\}$.

!!! abstract "Lloyd–Max conditions (necessary for optimality)"
    At the optimum:

    **Reproduction levels are conditional means:**
    $$
    y_i = \frac{\int_{t_{i-1}}^{t_i} f\, p(f)\, df}{\int_{t_{i-1}}^{t_i} p(f)\, df} = \mathbb E[F \mid t_{i-1} < F < t_i].
    $$

    **Decision thresholds are midpoints:**
    $$
    t_i = \frac{y_i + y_{i+1}}{2}.
    $$

    These are solved **iteratively** (alternating the two updates) from an
    initial guess (e.g., uniform).

**Properties.**

- Lloyd–Max always improves on uniform quantization when $p(f)$ is non-uniform.
- The **optimality is local** — the iteration can converge to a local minimum; use multiple starts for confidence.
- For Gaussian data, the Lloyd–Max quantizer is nearly uniform with slight over-weighting in the tails.

### Worked example 1.6 — 4-level Lloyd–Max for Gaussian

Distribution $F \sim \mathcal N(0, 1)$, $L = 4$ levels.

**Initial (uniform) thresholds and levels:**

$$
t = (-0.75,\ 0,\ 0.75),\qquad y = (-1.5,\ -0.5,\ 0.5,\ 1.5)
$$

(thresholds and levels chosen so the four cells cover $[-\infty, \infty)$ with
equal-width bins near the mode).

**Iteration 1 — update levels** using the Gaussian conditional means:

$$
y_1 = \mathbb E[F \mid F < -0.75] \approx -1.53,\quad
y_2 = \mathbb E[F \mid -0.75 < F < 0] \approx -0.34,
$$
$$
y_3 = \mathbb E[F \mid 0 < F < 0.75] \approx 0.34,\quad
y_4 = \mathbb E[F \mid F > 0.75] \approx 1.53.
$$

**Update thresholds** to midpoints:

$$
t_1 = \tfrac{y_1 + y_2}{2} \approx -0.935,\quad
t_2 = \tfrac{y_2 + y_3}{2} = 0,\quad
t_3 = \tfrac{y_3 + y_4}{2} \approx 0.935.
$$

**Iteration 2 — update levels** with new thresholds:

$$
y_1 \approx -1.51,\quad y_2 \approx -0.45,\quad y_3 \approx 0.45,\quad y_4 \approx 1.51.
$$

**Iteration 3:** minimal change — converged.

**Final Lloyd–Max (Gaussian, 4 levels):**

| Index | Threshold | Level |
|:---:|:---:|:---:|
| 1 | $-\infty$ | $-1.51$ |
| 2 | $-0.98$ | $-0.45$ |
| 3 | $0$ | $0.45$ |
| 4 | $0.98$ | $1.51$ |

**MSE.** Numerically integrate: $D \approx 0.117$. **Uniform quantization**
with the same 4 levels (over $[-3, 3]$) gives $D \approx 0.137$ — Lloyd–Max is
**~15% better** on MSE for the same number of levels.

### 1.6.3 Dither

**The problem.** Uniform quantization of a smooth gradient produces **false
contours** — visible stair-steps at the quantization boundaries. The eye is
sensitive to these because they are coherent, low-frequency artifacts.

**The solution — dither.** Add **random noise** before quantization:

$$
\boxed{\;Q_{\text{dithered}}(f) = Q(f + n),\qquad n \sim \mathcal U\!\left[-\tfrac{\Delta}{2}, \tfrac{\Delta}{2}\right].\;}
$$

The quantization error becomes **decorrelated** from the signal — it looks like
white noise rather than structured contours. **Perceptually** much better, even
though **numerically** the MSE is slightly **higher** (because the dither itself
adds noise).

**Trade-off.**

| Quantity | With dither | Without dither |
|---|---|---|
| MSE | Higher ($\sim 1.5\Delta^2/12$) | Lower ($\sim\Delta^2/12$) |
| Perceived quality | **Better** (no contours) | Worse (visible bands) |
| Cost | None (add noise) | — |

Dither is **standard** in audio CD production (16-bit with dither sounds better
than 16-bit without) and in image / video codecs when bit depth is very low
(4-bit greyscale, medical display).

### Worked example 1.7 — dither on a gradient

A smooth gradient from $0$ to $1$ is quantized to 4 levels ($\Delta = 1/4$).

**Without dither** — the output is a **staircase**: four flat bands. At each
boundary, a sharp jump creates a **false edge** visible to the eye.

**With dither** — add uniform noise in $[-\Delta/2, \Delta/2]$ before quantizing.
The output is a **speckled** version of the gradient: locally noisy, but
globally the gradient appears smooth — the eye integrates the noise and
perceives the underlying gradient. This is **perceptually superior** even though
the numerical MSE is higher.

```python
import numpy as np
import matplotlib.pyplot as plt

L = 4
x = np.linspace(0, 1, 1000)
y_no_dither = np.round(x * (L - 1)) / (L - 1)
y_dither = np.round((x + np.random.uniform(-0.5/L, 0.5/L, len(x))) * (L - 1)) / (L - 1)

fig, axes = plt.subplots(2, 1, figsize=(10, 6))
axes[0].plot(x, y_no_dither); axes[0].set_title('Without dither (false contours)')
axes[1].plot(x, y_dither);    axes[1].set_title('With dither (perceptually smooth)')
plt.show()
```

The numerical MSE of the dithered version is $\approx 1.4\times$ higher, but the
false contours are gone — the perceptual gain dominates.

---

## 1.7 Colour images

### 1.7.1 The trichromatic model

Human vision uses **three cone types** (S, M, L — short, medium, long wavelength
peak sensitivities), so any colour percept is **three-dimensional**. This is the
foundation of the **trichromatic theory** (Young–Helmholtz, 19th century) and
of every colour imaging system.

**Any colour can be matched** by a weighted sum of three **primaries** — the
weights are the **tristimulus values**.

### 1.7.2 CIE XYZ — the standard colorimetric space

The **CIE 1931 XYZ** space defines three imaginary primaries such that:

- **Y** represents **luminance** directly.
- **X** and **Z** are constructed so every physical colour has non-negative XYZ values.
- All visible colours map into the **chromaticity diagram** in $(x, y)$ where
  $x = X/(X+Y+Z)$, $y = Y/(X+Y+Z)$.

The **RGB colour primaries** of a display are a specific triangle in the
chromaticity diagram (sRGB is the standard). Colours inside the triangle are
displayable; colours outside are not (this is why "wide-gamut" displays are
needed for saturated colours).

### 1.7.3 RGB, YCbCr, HSV — three colour spaces for three purposes

| Space | Components | Use |
|---|---|---|
| **RGB** | Red, Green, Blue | Display, sensor output |
| **YCbCr** | Luma (Y), blue-difference chroma (Cb), red-difference chroma (Cr) | Compression (JPEG, MPEG), broadcast |
| **HSV / HSL** | Hue, Saturation, Value / Lightness | Colour selection, computer graphics |

**YCbCr** is **the** colour space for compression. Its **luma** component Y is a
weighted sum of RGB that roughly matches the eye's luminance sensitivity:

$$
\boxed{\;Y = 0.299\, R + 0.587\, G + 0.114\, B.\;}
$$

The **chroma** components:

$$
\boxed{\;C_b = 0.564(B - Y),\qquad C_r = 0.713(R - Y).\;}
$$

(Normalized to lie in $[-0.5, 0.5]$ before 8-bit quantization.)

**Inverse:**

$$
R = Y + 1.402\, C_r,\qquad G = Y - 0.344\, C_b - 0.714\, C_r,\qquad B = Y + 1.772\, C_b.
$$

### 1.7.4 Why compression uses YCbCr

Two reasons:

1. **Luminance is more perceptually important** than colour detail. The eye's
   spatial resolution for colour is much lower than for luminance (this is
   related to the different distributions of cone types — colour is resolved
   primarily by the S-cones, which are sparse).
2. **Chroma can be subsampled.** 4:2:0 chroma subsampling stores the full Y at
   every pixel but the Cb and Cr only at every 2×2 block — **half the raw
   data** with almost no perceptual loss. This is the foundation of JPEG and
   all video codecs.

!!! tip "The 4:x:y notation"
    In video, **chroma subsampling** is written as **J:a:b**:
    
    - **4:4:4** — no subsampling (studio quality).
    - **4:2:2** — chroma at half horizontal resolution (broadcast).
    - **4:2:0** — chroma at half horizontal AND vertical (web video, Blu-ray,
      JPEG standard).
    - **4:1:1** — chroma at quarter horizontal resolution (older formats).

### Worked example 1.8 — RGB ↔ YCbCr

**Forward.** RGB pixel $(R, G, B) = (200, 100, 50)$.

$$
Y = 0.299(200) + 0.587(100) + 0.114(50) = 59.8 + 58.7 + 5.7 = 124.2,
$$
$$
C_b = 0.564(B - Y) = 0.564(50 - 124.2) = -41.85,
$$
$$
C_r = 0.713(R - Y) = 0.713(200 - 124.2) = 54.05.
$$

For 8-bit storage, add 128 to the chroma:
$C_b^{8} = 86$, $C_r^8 = 182$. Luma stays at 124.

**Inverse.** Starting from $(Y, C_b, C_r) = (124.2, -41.85, 54.05)$:

$$
R = 124.2 + 1.402(54.05) = 124.2 + 75.78 = 200.0,
$$
$$
G = 124.2 - 0.344(-41.85) - 0.714(54.05) = 124.2 + 14.40 - 38.59 = 100.0,
$$
$$
B = 124.2 + 1.772(-41.85) = 124.2 - 74.17 = 50.0.
$$

**Perfect reconstruction** (in the ideal case). In 8-bit pipelines, the
rounding of Y and chroma values introduces $\pm 1$ LSB error — visually
imperceptible.

---

## 1.8 Linear shift-invariant systems

### 1.8.1 2-D convolution

A 2-D **linear shift-invariant (LSI)** system with **impulse response**
$h[m, n]$ computes

$$
\boxed{\;y[m, n] = (h * x)[m, n] = \sum_{k} \sum_{l} h[k, l]\, x[m-k, n-l].\;}
$$

Equivalently, in the **frequency domain**:

$$
\boxed{\;Y(\omega_1, \omega_2) = H(\omega_1, \omega_2)\, X(\omega_1, \omega_2)\;}
$$

where $X, Y, H$ are the 2-D DTFTs.

**Properties.**

- **Linearity:** $h * (\alpha x_1 + \beta x_2) = \alpha(h * x_1) + \beta(h * x_2)$.
- **Shift invariance:** $h * (x \text{ shifted}) = (h * x) \text{ shifted}$.
- **Commutativity:** $h_1 * h_2 = h_2 * h_1$.
- **Associativity:** $(h_1 * h_2) * x = h_1 * (h_2 * x)$.
- **Separability:** if $h[m, n] = h_1[m]\, h_2[n]$, then $h * x$ is achieved by filtering **along rows** with $h_1$, then **along columns** with $h_2$ — two 1-D passes, cost $O(N(M_1 + M_2))$ instead of $O(NM_1 M_2)$.

### 1.8.2 Common image filters

| Name | Kernel | Purpose |
|---|---|---|
| **Box blur** | $3\times 3$ uniform $1/9$ | Smoothing, noise reduction |
| **Gaussian** | $3\times 3$ separable (e.g., $h_1 = h_2 = (1, 2, 1)/4$) | Isotropic smoothing |
| **Sobel X** | $\begin{bmatrix} -1 & 0 & 1 \\ -2 & 0 & 2 \\ -1 & 0 & 1\end{bmatrix}$ | Horizontal gradient |
| **Sobel Y** | $\begin{bmatrix} -1 & -2 & -1 \\ 0 & 0 & 0 \\ 1 & 2 & 1\end{bmatrix}$ | Vertical gradient |
| **Laplacian** | $\begin{bmatrix} 0 & -1 & 0 \\ -1 & 4 & -1 \\ 0 & -1 & 0\end{bmatrix}$ | Edge detection (2nd derivative) |
| **Unsharp mask** | $x + \alpha(x - \text{blur}(x))$ | Sharpening |

### Worked example 1.9 — 2-D Gaussian blur as separable filter

The **2-D Gaussian** kernel

$$
h[m, n] = \exp\!\left(-\frac{m^2 + n^2}{2\sigma^2}\right)
$$

**is separable**:

$$
h[m, n] = h_1[m]\, h_2[n],\qquad h_1[m] = e^{-m^2/(2\sigma^2)},\qquad h_2[n] = e^{-n^2/(2\sigma^2)}.
$$

**Direct convolution cost:** $K^2$ multiplications per pixel (for a $K \times K$ kernel).

**Separable cost:** $K + K = 2K$ multiplications per pixel — a factor
$K/2$ reduction. For a $15 \times 15$ Gaussian, **7.5× faster** — the difference
between real-time and not.

**Python implementation:**

```python
import numpy as np
from scipy.ndimage import convolve1d

def gaussian_kernel_1d(sigma):
    half = int(3 * sigma)
    x = np.arange(-half, half + 1)
    return np.exp(-x**2 / (2 * sigma**2))

def gaussian_blur_separable(image, sigma):
    k = gaussian_kernel_1d(sigma)
    k /= k.sum()
    # Filter rows, then columns
    tmp = convolve1d(image, k, axis=1)     # horizontal
    return convolve1d(tmp, k, axis=0)       # vertical
```

The two-pass approach gives the **same result** as full 2-D convolution but at
a fraction of the cost.

---

## 1.9 2-D Fourier and Z transforms

### 1.9.1 The 2-D discrete Fourier transform

For an $M \times N$ image $x[m, n]$:

$$
\boxed{\;X[k, l] = \sum_{m=0}^{M-1}\sum_{n=0}^{N-1} x[m,n]\, e^{-j2\pi\left(\frac{km}{M} + \frac{ln}{N}\right)},\quad k = 0, \dots, M-1,\ l = 0, \dots, N-1.\;}
$$

**Inverse:**

$$
\boxed{\;x[m,n] = \frac{1}{MN}\sum_{k=0}^{M-1}\sum_{l=0}^{N-1} X[k,l]\, e^{j2\pi\left(\frac{km}{M} + \frac{ln}{N}\right)}.\;}
$$

### 1.9.2 Key properties

| Property | Statement |
|---|---|
| **Separability** | 2-D DFT = 1-D DFT along rows, then along columns (or vice versa) |
| **Linearity** | $\text{DFT}(a x + b y) = a \text{DFT}(x) + b \text{DFT}(y)$ |
| **Shift** | $\text{DFT}(x[m-m_0, n-n_0]) = X[k,l] e^{-j2\pi(km_0/M + ln_0/N)}$ |
| **Conjugate symmetry** (real $x$) | $X[M-k, N-l] = X^*[k,l]$ |
| **Periodicity** | $X[k+M, l] = X[k, l]$; $X[k, l+N] = X[k, l]$ |
| **Parseval** | $\sum_{m,n} \lvert x \rvert^2 = \frac{1}{MN}\sum_{k,l}\lvert X \rvert^2$ |
| **Convolution** | $\text{DFT}(x * h) = X \cdot H$ (circular convolution) |

### 1.9.3 Interpretation

The **magnitude** $|X[k,l]|$ describes the **amplitude** of each 2-D frequency
$(k, l)$ in the image; the **phase** $\angle X[k,l]$ describes the **spatial
alignment** of that frequency's contribution.

- **Low frequencies** (near the origin $(0,0)$ in the shifted spectrum) correspond to **slow variations** — the overall shape, the smooth content of the image.
- **High frequencies** (near the edges of the $[0, M) \times [0, N)$ grid) correspond to **fast variations** — edges, textures, and noise.

**For display** it is conventional to shift the spectrum so that low frequencies
are at the centre:

```python
X_shifted = np.fft.fftshift(np.fft.fft2(image))
```

### Worked example 1.10 — 2-D DFT of a rectangle

Let $x[m,n] = 1$ for $0 \le m < P, 0 \le n < Q$, and $0$ elsewhere in an
$M \times N$ image.

**By separability:**

$$
X[k,l] = \underbrace{\left(\sum_{m=0}^{P-1} e^{-j2\pi km/M}\right)}_{X_1[k]} \cdot \underbrace{\left(\sum_{n=0}^{Q-1} e^{-j2\pi ln/N}\right)}_{X_2[l]}.
$$

Each factor is a **geometric series**:

$$
X_1[k] = e^{-j\pi k(P-1)/M}\cdot\frac{\sin(\pi k P/M)}{\sin(\pi k/M)}.
$$

Similarly for $X_2[l]$. So the 2-D DFT of a rectangle is a **product of two 1-D
sinc-like patterns** — a separable 2-D Dirichlet kernel.

**The spectrum** has a broad main lobe around the origin (DC) and ringing
side-lobes along both axes. The shape of the main lobe along axis $k$ is
inversely proportional to $P$ — a larger rectangle has a narrower main lobe,
the image-domain analogue of a longer signal having a narrower bandwidth.

### 1.9.4 The 2-D Z transform

For a 2-D sequence $x[m,n]$:

$$
\boxed{\;X(z_1, z_2) = \sum_{m}\sum_{n} x[m,n]\, z_1^{-m} z_2^{-n}.\;}
$$

**Properties** (analogous to 1-D):

- The **region of convergence** is a 2-D annulus / torus.
- A 2-D **linear shift-invariant system** has **transfer function** $H(z_1, z_2)$
  with **impulse response** $h[m,n]$.
- **Stability** requires $H(z_1, z_2)$ to have no poles on the unit torus
  $|z_1| = |z_2| = 1$.
- **Separability:** if $h[m,n] = h_1[m] h_2[n]$, then $H(z_1, z_2) = H_1(z_1) H_2(z_2)$ —
  a **separable** 2-D Z transform.

2-D Z transforms appear mostly in **filter design** and **state-space analysis**
of image systems; the 2-D DFT is more common in practice, since it operates
directly on the sampled image.

---

## 1.10 Matrix theory, block matrices, and Kronecker products

### 1.10.1 Kronecker product — definition

For matrices $A \in \mathbb R^{m \times n}$ and $B \in \mathbb R^{p \times q}$:

$$
\boxed{\;A \otimes B = \begin{bmatrix} a_{11} B & a_{12} B & \cdots & a_{1n} B \\ a_{21} B & a_{22} B & \cdots & a_{2n} B \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} B & a_{m2} B & \cdots & a_{mn} B\end{bmatrix} \in \mathbb R^{mp \times nq}.\;}
$$

**Key properties.**

| Property | Statement |
|---|---|
| **Mixed-product** | $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$ |
| **Transpose** | $(A \otimes B)^\top = A^\top \otimes B^\top$ |
| **Inverse** | $(A \otimes B)^{-1} = A^{-1} \otimes B^{-1}$ (if both invertible) |
| **Vec identity** | $\text{vec}(AXB) = (B^\top \otimes A)\text{vec}(X)$ |
| **Spectrum** | Eigenvalues of $A \otimes B$ are $\{\lambda_i(A)\lambda_j(B)\}$ |

The **vec identity** is the key to 2-D signal processing: it converts **2-D
matrix operations** into **1-D matrix-vector multiplications**, which is what
GPUs and DSPs are optimised for.

### 1.10.2 The 2-D convolution as a block Toeplitz matrix

Given a 2-D convolution

$$
y[m,n] = \sum_{k,l} h[k,l]\, x[m-k, n-l],
$$

we can **vectorize** both sides. If $x, y \in \mathbb R^{MN}$ are the
**lexicographic** stackings of the images (row-major or column-major), and $H$
is the **block Toeplitz** matrix built from $h$, then:

$$
\boxed{\;y = H x.\;}
$$

**Structure of $H$.** If $h$ is separable, $H = H_2 \otimes H_1$ (Kronecker
product of the two 1-D Toeplitz filters). If not, $H$ is a **block Toeplitz
matrix with Toeplitz blocks** — the structure exploited by 2-D convolution
algorithms.

### Worked example 1.11 — 2-D filtering as a Kronecker product

Consider a **separable** $2 \times 2$ averaging filter:

$$
h = \begin{bmatrix} 1/4 & 1/4 \\ 1/4 & 1/4 \end{bmatrix} = \begin{bmatrix} 1/2 \\ 1/2\end{bmatrix}\begin{bmatrix} 1/2 & 1/2\end{bmatrix} = h_1 h_2^\top.
$$

**Direct 2-D filtering** of a $3 \times 3$ image $X$:

$$
Y = h * X.
$$

**Kronecker form.** Build the 1-D Toeplitz matrices $H_1$ and $H_2$ from
$h_1 = (1/2, 1/2)$ and $h_2 = (1/2, 1/2)$:

$$
H_1 = \begin{bmatrix} 1/2 & 0 \\ 1/2 & 1/2 \\ 0 & 1/2\end{bmatrix},\qquad
H_2 = \begin{bmatrix} 1/2 & 0 \\ 1/2 & 1/2 \\ 0 & 1/2\end{bmatrix}.
$$

(These are the **Toeplitz** matrices representing 1-D convolution with the
1-D filters — input length 3, output length 3, "same" convolution.)

**The vec identity:**

$$
\text{vec}(Y) = \text{vec}(H_1 X H_2^\top) = (H_2 \otimes H_1)\,\text{vec}(X).
$$

The $9 \times 9$ matrix $H_2 \otimes H_1$ is the **block Toeplitz representation**
of the 2-D filter. It has the structure

$$
H_2 \otimes H_1 = \begin{bmatrix} (1/2)H_1 & 0 & 0 \\ (1/2)H_1 & (1/2)H_1 & 0 \\ 0 & (1/2)H_1 & (1/2)H_1\end{bmatrix}.
$$

**Computation.** Applying $H_2 \otimes H_1$ to a vectorized image costs
$O(M^2 N^2)$ flops if done naively — but the **separable structure** means it
can be done as

$$
\text{vec}(Y) = (H_2 \otimes I_M)\,(I_N \otimes H_1)\,\text{vec}(X),
$$

i.e., filter **columns with $H_1$**, then **rows with $H_2$** — costing
$O(MN(M+N))$ flops. For a $1000 \times 1000$ image, this is $2 \times 10^9$ vs
$10^{12}$ — a **500× speedup** from exploiting separability.

### 1.10.3 Why Kronecker products matter for image processing

Three concrete reasons:

1. **Vectorization.** Every 2-D operation can be represented as a single matrix-vector multiplication. This is how GPUs and tensor processing units (TPUs) accelerate image processing — they love matrix multiplies.
2. **Separable filters.** Many important filters (Gaussian, box, Sobel, DCT basis functions, wavelet filters) are separable — they can be factored as $h = h_1 \otimes h_2$ and implemented as two 1-D passes.
3. **Algebraic analysis.** Eigenvalues of $A \otimes B$ are products $\lambda_i \lambda_j$; singular values of $A \otimes B$ are products $\sigma_i \sigma_j$. This makes it possible to derive closed-form properties of 2-D filters from 1-D properties — the basis of every fast algorithm.

### Worked example 1.12 — 2-D DFT as a Kronecker product

The 2-D DFT of an $N \times N$ image can be written as

$$
\boxed{\;F_{\text{2D}} = F_N \otimes F_N,\;}
$$

where $F_N$ is the $N$-point 1-D DFT matrix.

**Proof.** From the 2-D DFT definition, using separability of the exponential
kernel:

$$
X[k,l] = \sum_m\sum_n x[m,n] e^{-j2\pi km/N} e^{-j2\pi ln/N}.
$$

This is exactly the matrix-vector identity $(F_N \otimes F_N)\,\text{vec}(x)$
with $F_N[k,m] = e^{-j2\pi km/N}$ and the vec convention appropriate for the
image stacking.

**Implication.** **The 2-D FFT is faster than the 1-D FFT** would suggest:
$O(N^2 \log N)$ via separability, not $O(N^4)$ via the naive 2-D DFT
formula. Every image-processing library uses separable FFT — it's the reason
2-D filtering in the frequency domain is practical.

**Java implementation** (for the matrix structure):

```java
public class Kronecker2DDFT {
    public static double[][] kron(double[][] A, double[][] B) {
        int m = A.length, n = A[0].length;
        int p = B.length, q = B[0].length;
        double[][] C = new double[m * p][n * q];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                for (int k = 0; k < p; k++)
                    for (int l = 0; l < q; l++)
                        C[i*p + k][j*q + l] = A[i][j] * B[k][l];
        return C;
    }

    public static double[][] dftMatrix(int N) {
        double[][] F = new double[N][N];
        for (int k = 0; k < N; k++)
            for (int n = 0; n < N; n++) {
                double ang = -2 * Math.PI * k * n / N;
                F[k][n] = Math.cos(ang);  // real part for simplicity
            }
        return F;
    }

    public static void main(String[] args) {
        int N = 4;
        double[][] F = dftMatrix(N);
        double[][] F2D = kron(F, F);
        System.out.println("2-D DFT matrix (N=4) size: "
                          + F2D.length + " x " + F2D[0].length);
        // For an 8x8 image, F2D is 64 x 64.
    }
}
```

---

## 1.11 Chapter summary

1. **Photometry vs. radiometry** — human vision weights the physical light by $V(\lambda)$, peaking at ~555 nm. **Luminance** (cd/m²) is the physical quantity; **brightness** is the subjective percept.
2. **Contrast** — Weber $C_W = (L_{\text{obj}} - L_{\text{bg}})/L_{\text{bg}}$ for targets; Michelson $C_M = (L_{\max} - L_{\min})/(L_{\max} + L_{\min})$ for patterns.
3. **The eye** — rods (low light) and cones (colour, fine detail, concentrated in the fovea). The **CSF** is band-pass, peaking at 4–8 cycles/degree, and drops off sharply at high frequencies — this is the physiological basis of perceptual compression.
4. **Camera pipeline** — optics → sensor → ADC → digital processing. The **Bayer CFA** requires demosaicing; the **anti-alias filter** trades resolution for avoiding Moiré.
5. **2-D sampling theorem** — sampling must exceed twice the highest spatial frequency in both directions. Aliasing produces Moiré patterns at the beat frequencies between the image content and the sampling lattice.
6. **Quantization** — uniform quantizer MSE $\Delta^2/12$, giving the **6 dB per bit** rule. **Lloyd–Max** improves MSE for non-uniform distributions, typically by 10–20% for Gaussian data. **Dither** trades numerical MSE for perceptual smoothness — essential for low-bit-depth displays.
7. **Colour** — **YCbCr** separates luma (perceptually important) from chroma (subsamplable); this is why compression uses it. **4:2:0** halves chroma bandwidth with almost no perceptual cost.
8. **LSI systems** — 2-D convolution; **separable filters** give massive computational savings ($O(M+N)$ per pixel instead of $O(MN)$). Gaussian, box, and Sobel are all separable.
9. **2-D Fourier** — separable DFT, obtained as row-DFT then column-DFT. **Z-transform** generalises to 2-D rational transfer functions; separability of the transfer function mirrors separability of the filter.
10. **Kronecker products** — every 2-D operation can be written as matrix-vector multiplication using $A \otimes B$. The **vec identity** $\text{vec}(AXB) = (B^\top \otimes A)\text{vec}(X)$ is the bridge between 2-D convolution and 1-D matrix algebra. Separable filters are $h = h_1 \otimes h_2$ — the computational shortcut that makes 2-D processing feasible.

---

## 1.12 Quiz

### A. Multiple choice

1. The **Michelson contrast** of a pattern with luminances 10 and 50 cd/m² is:
   (a) 0.4 (b) 0.67 (c) 4 (d) 0.33

2. Weber's law states that the just-noticeable relative luminance change is:
   (a) constant (~1%) (b) linear in $L$ (c) quadratic in $L$ (d) zero

3. The **contrast sensitivity function** peaks at approximately:
   (a) 0 cycles/degree (b) 4–8 cycles/degree (c) 50 cycles/degree (d) infinite

4. The **Nyquist frequency** for a sampled image with pixel spacing $\Delta x$ is:
   (a) $\Delta x$ (b) $1/\Delta x$ (c) $1/(2\Delta x)$ (d) $2\Delta x$

5. Aliasing in 2-D images produces:
   (a) no visible effect (b) Moiré patterns (c) increased resolution (d) lower noise

6. Uniform quantization with step $\Delta$ has MSE approximately:
   (a) $\Delta^2$ (b) $\Delta^2/2$ (c) $\Delta^2/12$ (d) $\Delta/12$

7. The Lloyd–Max quantizer:
   (a) is uniform (b) minimises MSE for a known distribution (c) has no reproduction levels (d) is unstable

8. Adding **dither** before quantization:
   (a) increases MSE slightly, improves perceptual quality (b) decreases MSE (c) has no effect (d) always worse

9. In **YCbCr**, the **Y** component represents:
   (a) blue channel (b) luminance (c) chroma (d) colour saturation

10. 4:2:0 chroma subsampling:
    (a) keeps all chroma (b) halves chroma in both directions (c) removes chroma entirely (d) reduces luma

11. A **separable** filter in 2-D:
    (a) cannot exist (b) factors as $h[m,n] = h_1[m]\,h_2[n]$ (c) is non-linear (d) is always a Gaussian

12. The **Kronecker product** $A \otimes B$ of an $m \times n$ and a $p \times q$ matrix is:
    (a) $m \times n$ (b) $p \times q$ (c) $mp \times nq$ (d) $(m+p) \times (n+q)$

### B. Fill in the blank

1. Weber contrast for an object on a background is $C_W = (L_{\text{obj}} - L_{\text{bg}})/$\_\_\_\_.
2. The luminous efficiency function $V(\lambda)$ peaks at $\lambda \approx$ \_\_\_\_ nm.
3. A camera sensor with $W$ pixels can resolve at most $W/$\_\_\_\_ cycles.
4. The 2-D sampling theorem requires $1/\Delta x > $\_\_\_\_ (twice the highest spatial frequency).
5. The variance of uniform quantization error is $\sigma_e^2 = \Delta^2/$\_\_\_\_.
6. The Lloyd–Max reproduction level is the \_\_\_\_ of the corresponding decision interval.
7. Dither is added \_\_\_\_ quantization.
8. In YCbCr, the Y component is a weighted sum of RGB: $Y = 0.299R + 0.587G + $\_\_\_\_.
9. A separable filter $h[m,n] = h_1[m]h_2[n]$ has Z-transform $H(z_1,z_2) = $\_\_\_\_.
10. The vec identity for $Y = AXB$ is $\text{vec}(Y) = (B^\top \otimes $\_\_\_\_$) \text{vec}(X)$.

### C. Define / explain

1. Distinguish **radiometry** from **photometry**. Why is the human luminous efficiency function important for display engineering?
2. Define **Weber** and **Michelson** contrast. Give one application where each is the natural measure.
3. State the **contrast sensitivity function** and explain its role in JPEG and JPEG 2000 compression.
4. Describe the **camera pipeline** from optics to digital output. What role does the anti-alias filter play?
5. State the **2-D sampling theorem**. Derive the aliasing condition and explain what Moiré patterns are.
6. Define the **uniform quantizer** and derive its MSE approximation. What is the **6 dB per bit** rule?
7. Derive the **Lloyd–Max** optimality conditions (both reproduction and threshold) and describe the iterative algorithm.
8. Explain **dither**. Why does adding noise improve perceptual quality even though it increases numerical MSE?
9. Convert an RGB triplet $(R, G, B) = (128, 64, 192)$ to YCbCr. Explain the perceptual reason YCbCr is used for compression.
10. State the **2-D convolution** and **2-D DFT**. Show that a separable filter can be implemented as two 1-D operations.
11. Define the **Kronecker product**. State and interpret the **vec identity**. Give one application where it is essential.

---

## 1.13 Answer key

??? success "A. Multiple choice"
    1. **(b)** — $(50-10)/(50+10) = 40/60 = 0.667$.
    2. **(a)** — constant (~1%).
    3. **(b)** — 4–8 cycles/degree.
    4. **(c)** — $1/(2\Delta x)$.
    5. **(b)** — Moiré patterns.
    6. **(c)** — $\Delta^2/12$.
    7. **(b)** — minimises MSE for known distribution.
    8. **(a)** — slightly higher MSE, much better perceptual quality.
    9. **(b)** — luminance.
    10. **(b)** — halves chroma in both directions.
    11. **(b)** — factors as $h_1[m]h_2[n]$.
    12. **(c)** — $mp \times nq$.

??? success "B. Fill in the blank"
    1. $L_{\text{bg}}$ (background luminance)
    2. $555$
    3. $2$
    4. $2F_x$ (twice the highest spatial frequency)
    5. $12$
    6. centroid (or conditional mean)
    7. before
    8. $0.114B$
    9. $H_1(z_1)\,H_2(z_2)$
    10. $A$

??? success "C. Define / explain — model answers"
    1. *Radiometry* measures light in physical units of energy (watt, W/m²/sr). *Photometry* weights light by the eye's sensitivity $V(\lambda)$, giving **luminous flux** (lumen), **luminance** (cd/m²), etc. The luminous efficiency function is important for display engineering because it determines what luminance range is perceivable and how display gamma / bit depth should be chosen. It is also the basis of perceptual image and video compression — you allocate bits to frequencies where the eye is most sensitive (CSF).
    2. *Weber contrast:* $C_W = (L_{\text{obj}} - L_{\text{bg}})/L_{\text{bg}}$ — relative to background. Suitable for **target detection** (a defect on a uniform background, a lesion on tissue). *Michelson contrast:* $C_M = (L_{\max} - L_{\min})/(L_{\max} + L_{\min})$ — relative to average. Suitable for **periodic patterns** (gratings, checkerboards), where both extremes are "objects". Weber can be unbounded (and even negative); Michelson is in $[0,1]$.
    3. *CSF:* inverse of minimum detectable contrast at a given spatial frequency. Band-pass, peak at 4–8 cycles/degree. In **JPEG**, the CSF-derived **quantization table** coarsens high-frequency DCT coefficients (which the eye cannot resolve anyway) while preserving low and mid frequencies. In **JPEG 2000**, the wavelet coefficients are similarly weighted — finer subbands are quantized more aggressively. The perceptual gain is 5–15 dB in PSNR at the same bit rate, compared with uniform quantization.
    4. *Pipeline:* optics (lens, aperture, diffraction limit) → sensor (photons to charge / voltage, pixel size, CFA) → analog front end (amplifier, anti-alias filter, ADC) → digital processing (demosaic, white balance, gamma, denoise, compress). *Anti-alias filter role:* a slight optical lowpass, attenuates frequencies above the sensor's Nyquist — suppresses Moiré patterns from high-frequency content. Cost: softens the image slightly. Modern high-resolution sensors often omit it (aliasing is rare at very high pixel density, and the lens provides enough rolloff).
    5. *2-D sampling theorem:* a band-limited image with 2-D spectrum zero outside $[-F_x, F_x] \times [-F_y, F_y]$ can be reconstructed exactly from samples on a lattice with spacing $\Delta x, \Delta y$ iff $1/\Delta x > 2F_x$ and $1/\Delta y > 2F_y$. *Aliasing condition:* if the sampling frequencies are below $2F_x, 2F_y$, the periodic extensions of the spectrum overlap and the high-frequency content **folds** into lower frequencies. *Moiré:* the visible manifestation — beat patterns produced by the difference between the image's own frequencies and the sampling frequencies. Common in photographs of fine grids, screens, or fabrics.
    6. *Uniform quantizer:* $Q(f) = \Delta \lfloor (f+A)/\Delta \rfloor + \Delta/2$ for input $f \in [-A, A]$; step $\Delta = 2A/L$ where $L = 2^b$. *MSE:* the quantization error is approximately uniform on $[-\Delta/2, \Delta/2]$, giving variance $\Delta^2/12$. *6 dB per bit:* for full-scale signal statistics, $\text{SNR}_Q \approx 6b + 4.77 - 20\log_{10}(f_{\max}/\sigma_f)$ dB — each additional bit improves SNR by 6 dB. For images, $\text{SNR}_Q \approx 6b - 1$ dB.
    7. *Lloyd–Max conditions:* (i) **reproduction level** $y_i$ is the **centroid** of the probability mass in the corresponding decision interval: $y_i = \mathbb E[F \mid t_{i-1} < F < t_i]$; (ii) **decision threshold** $t_i$ is the **midpoint** between adjacent reproduction levels: $t_i = (y_i + y_{i+1})/2$. *Algorithm:* start with an initial guess (e.g., uniform thresholds and levels), then alternate (i) and (ii) until convergence. Each step **strictly decreases** the MSE; convergence is guaranteed (to a local minimum) by the alternating-projection property. Multiple restarts are needed if the distribution is multimodal.
    8. *Dither:* add random noise $n \sim \mathcal U[-\Delta/2, \Delta/2]$ before quantization. The quantization error becomes decorrelated from the signal — the output is no longer structured (no false contours). *Why perceptual > numerical:* the human visual system is much more sensitive to **coherent low-frequency structure** (bands, edges) than to **white noise** (which the eye averages out). A high-MSE dithered image looks smoother than a low-MSE banded image. This is an example of **MSE not being the right fidelity metric** for images — a well-known limitation that motivated SSIM, MS-SSIM, and other perceptual metrics.
    9. *RGB → YCbCr* for $(R,G,B) = (128, 64, 192)$:
       $Y = 0.299(128) + 0.587(64) + 0.114(192) = 38.27 + 37.57 + 21.89 = 97.73$.
       $C_b = 0.564(B - Y) = 0.564(192 - 97.73) = 53.14$.
       $C_r = 0.713(R - Y) = 0.713(128 - 97.73) = 21.58$.
       With 8-bit normalization (add 128): $Y = 98$, $C_b = 181$, $C_r = 150$.
       *Why YCbCr for compression:* (i) the eye has much lower spatial resolution for colour than for luminance; (ii) chroma components can therefore be **subsampled** (4:2:0 halves chroma bandwidth) with almost no perceptual loss; (iii) Y is dominated by the green channel, matching the eye's sensitivity peak. Result: ~50% data reduction with no visible artefact — the foundation of JPEG and video compression.
    10. *2-D convolution:* $y[m,n] = \sum_{k,l} h[k,l] x[m-k,n-l]$. *2-D DFT:* $X[k,l] = \sum_{m,n} x[m,n] e^{-j2\pi(km/M + ln/N)}$. *Separability:* if $h[m,n] = h_1[m] h_2[n]$, then convolution is $\sum_k h_1[k] \left(\sum_l h_2[l] x[m-k, n-l]\right) = $ filter rows by $h_2$, then columns by $h_1$ (order commutes). Computational cost is $O(N(M_1 + M_2))$ per pixel instead of $O(N M_1 M_2)$ — a factor of $\min(M_1, M_2)$ savings. Gaussian, box, Sobel, and DCT basis functions are all separable.
    11. *Kronecker product:* $A \otimes B$ is the block matrix with blocks $a_{ij} B$. Result size: $mp \times nq$ for $A$ $m \times n$ and $B$ $p \times q$. *Vec identity:* $\text{vec}(AXB) = (B^\top \otimes A)\,\text{vec}(X)$ — **crucial** because it converts 2-D matrix operations into 1-D vector operations. *Applications:* (a) 2-D convolution as a matrix-vector product with a block-Toeplitz (Kronecker-structured) matrix; (b) 2-D DFT as $F_N \otimes F_N$ — the basis of the separable 2-D FFT; (c) **separable filtering** — a separable 2-D filter is a single Kronecker product, and the filter is applied as two 1-D passes; (d) fast algorithms for DCT, wavelet transforms, and 2-D linear systems.

---

## Further reading

- Jain, *Fundamentals of Digital Image Processing*, PHI, 1995 — Chapters 1–3 are the direct source for this chapter.
- Gonzalez & Woods, *Digital Image Processing*, 4th ed., Pearson, 2018 — Chapters 2–3 for the perceptual and sampling material.
- Pratt, *Digital Image Processing*, 4th ed., Wiley, 2007 — the encyclopaedic reference, especially for colour spaces and quantization.
- Oppenheim & Schafer, *Discrete-Time Signal Processing*, PHI — Chapter 13 for 2-D sampling and filtering.
- Van Loan, *Computational Frameworks for the Fast Fourier Transform*, SIAM, 1992 — the definitive reference for Kronecker-structured algorithms.
- Brewster, *Kronecker Products and Matrix Calculus*, Dover — a compact algebra reference.
- Sharma & Trussell, "Digital Color Imaging", *IEEE Trans. Image Process.*, 6(7), 1997 — the standard colour-space reference.
- Poynton, *Digital Video and HD*, 2nd ed., Morgan Kaufmann, 2012 — the definitive colour and gamma reference.
- Mannos & Sakrison, "The Effects of a Visual Fidelity Criterion on the Encoding of Images", *IEEE Trans. Inf. Theory*, 20(4), 1974 — the original CSF-based compression work.
- `scipy.ndimage`, `scikit-image`, `Pillow`, `OpenCV` documentation.

---

## Suggested tutorial and lab problems for this chapter

1. **Contrast calculation.** For a pattern with luminances $\{L_{\min}, L_{\max}\} = \{20, 200\}$ cd/m², compute both Weber (for each as background) and Michelson contrasts. In what context does each definition give the more useful number?

2. **Dynamic range.** A camera has 12-bit ADC and 85 dB dynamic range. What is the effective noise floor in LSBs? At what bit depth would a 70 dB consumer sensor be equivalent?

3. **Human vision and compression.** Using the Mannos–Sakrison CSF model, compute the sensitivity at 1, 5, 15, 30 cycles/degree. By what factor should you weight the DCT coefficients at these frequencies to achieve perceptual transparency?

4. **2-D sampling and aliasing.** A test image contains a sinusoidal grating with period 8 pixels and a step edge. Sample it at 1×, 2×, 4× subsampling. Identify where Moiré appears in the subsampled versions. Repeat with anti-alias filtering before subsampling.

5. **Lloyd–Max quantizer.** Implement the Lloyd–Max algorithm for 4- and 8-level quantization of a Gaussian distribution. Plot the MSE vs. bit depth and compare with uniform quantization.

6. **Dither.** Quantize a gradient to 3 bits with and without dither. Compare the visual result. Compute both MSE and SSIM for each version — comment on which metric matches perception better.

7. **RGB ↔ YCbCr.** Take a colour image, convert to YCbCr, subsample the chroma to 4:2:0, upsample back with bilinear interpolation, and convert to RGB. Compute the PSNR vs. original. Comment on which image regions show the largest errors (edges, saturated colours, etc.).

8. **Separable vs. direct convolution.** Implement a 2-D Gaussian blur with (a) direct convolution with a 15×15 kernel, (b) separable two-pass convolution. Verify identical outputs; measure and compare execution times.

9. **2-D FFT.** Compute the 2-D DFT of a synthetic image (rectangle, Gaussian, and a natural image). Plot the magnitude spectrum. Verify Parseval's theorem numerically.

10. **Kronecker product.** For a $3\times 3$ averaging filter $h = h_1 h_2^\top$ with $h_1 = h_2 = (1,1,1)/\sqrt{3}$, build the 9×9 Kronecker matrix and verify it produces the same output as direct 2-D convolution on a small test image.

11. **Block-Toeplitz structure.** For a non-separable $2\times 2$ filter $h = \begin{bmatrix} 1 & 2 \\ 3 & 4\end{bmatrix}$, build the block-Toeplitz matrix by lexicographic vectorization. Verify that applying it to a vectorized image matches the 2-D convolution. Show that the matrix has Toeplitz blocks but is not itself Toeplitz (unless the filter is separable).

12. **Real-world colour pipeline.** Take a raw camera file (or simulate one) — demosaic, white balance, gamma-correct, convert to YCbCr, and visualize each stage. Comment on the visual and numerical changes at each step.