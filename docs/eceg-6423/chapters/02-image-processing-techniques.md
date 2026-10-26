# Chapter 2 — Image Processing Techniques

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 2 tutorial + 1 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** image transforms (2-D orthogonal and unitary, DFT, DCT,
    DST, Walsh, Hadamard, Haar, Slant, KL, SVD) · image enhancement (point
    operations, histograms, spatial filtering, edge enhancement, colour
    enhancement) · image restoration (observation models, inverse and Wiener
    filtering, geometric mean, nonlinear filters, splines, constrained least
    squares) · image compression and reconstruction from projections (data
    rates, predictive, transform, vector, BTC, wavelet coding, JPEG/JPEG 2000,
    Radon transform, back-projection, fan-beam, algebraic reconstruction).

This is the longest chapter in the course — it covers the four pillars of
classical image processing: **transform**, **enhance**, **restore**, **compress**.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Apply **2-D orthogonal and unitary transforms** to images, and state the
   energy-compaction properties of the DCT, KL, and DWT.
2. Choose between **DFT, DCT, Walsh–Hadamard, Haar, Slant, KL, and SVD** based
   on their properties (energy compaction, computational cost, basis structure,
   data dependence).
3. Perform **point operations** (contrast stretching, clipping, thresholding,
   density slicing) and **histogram-based** processing (equalization,
   modification, specification).
4. Design **spatial filters** (averaging, LP, HP, BP, directional, median) and
   **edge-enhancement filters** (FIR and IIR, Laplacian, Sobel, unsharp mask).
5. Model image **degradation** and design **restoration filters**: inverse,
   Wiener (minimum MSE), geometric mean, constrained least squares.
6. Design **predictive, transform, vector, block-truncation, and wavelet** image
   coders; state the JPEG and JPEG 2000 pipelines.
7. Derive and apply the **Radon transform**, **back-projection**, **filtered
   back-projection**, **fan-beam**, and **algebraic reconstruction** algorithms.
8. Implement all of the above in **Python, Java, and C#**.

---

## 2.1 The four pillars

Image processing techniques fall into four families:

| Family | Input → Output | Question it answers | Representative methods |
|---|---|---|---|
| **Transform** | Image → Coefficients | "What is the image made of?" | DFT, DCT, DWT, KL, SVD |
| **Enhancement** | Image → Better image | "How do I make it look better?" | Histogram, filtering, sharpening |
| **Restoration** | Degraded image → Original estimate | "How do I undo the damage?" | Inverse, Wiener, CLS |
| **Compression** | Image → Bits | "How do I represent it compactly?" | DPCM, DCT, VQ, wavelet |

**The transform** is the foundation — every compression and some enhancement and
restoration methods operate **in a transform domain**. The transform is chosen so
that the image's energy concentrates in a few coefficients (energy compaction),
which is what makes compression possible.

---

## 2.2 Image transforms

### 2.2.1 Orthogonal and unitary transforms

A **linear transform** of a vector $\mathbf x \in \mathbb R^N$ is a matrix-vector
product $\mathbf y = A \mathbf x$. The transform is **unitary** if

$$
\boxed{\;A\, A^H = A^H A = I \quad \Longleftrightarrow \quad A^{-1} = A^H.\;}
$$

For real matrices ($A^H = A^\top$), it is **orthogonal**: $A^{-1} = A^\top$.

**Consequences of orthogonality/unitarity.**

- **Parseval (energy preservation):** $\|\mathbf y\|^2 = \|\mathbf x\|^2$.
- **Inverse is cheap:** $\mathbf x = A^H \mathbf y$ — no matrix inverse required.
- **Condition number 1** — numerically perfectly conditioned.

**2-D extension.** For an $N \times N$ image $X$, the 2-D transform is

$$
\boxed{\;Y = A\, X\, B^\top\;}
$$

with inverse

$$
\boxed{\;X = A^H\, Y\, B.\;}
$$

If $A = B$ (symmetric transform), this simplifies to $Y = A X A^\top$.

**Separability.** The 2-D transform is **separable**: apply the 1-D transform $B$
to the columns of $X$, then the 1-D transform $A$ to the rows. Cost: $O(N^2 (A_{\text{cost}} + B_{\text{cost}}))$,
much less than the $O(N^4)$ of direct 2-D multiplication.

### 2.2.2 2-D DFT

Already introduced in Chapter 1 §1.9.1. Recap:

$$
\boxed{\;X[k,l] = \sum_{m=0}^{M-1}\sum_{n=0}^{N-1} x[m,n]\, e^{-j2\pi(km/M + ln/N)},\quad k,l \in [0, M-1] \times [0, N-1].\;}
$$

**Properties.**

- **Basis:** complex exponentials (sines and cosines) — pure frequencies.
- **Inverse:** $x[m,n] = \frac{1}{MN}\sum_{k,l} X[k,l] e^{j2\pi(km/M + ln/N)}$.
- **Separable:** 1-D DFT on rows, then on columns.
- **Fast:** FFT gives $O(MN\log MN)$ instead of $O(M^2 N^2)$.
- **Not optimal for compression:** the DFT of a typical image has energy spread over many coefficients (beyond the DC), because most images are not periodic — boundary discontinuities create high-frequency content.
- **Two main drawbacks for compression:** (i) it produces **complex** coefficients (need to encode magnitude and phase); (ii) the **discontinuities** at block boundaries waste energy.

The **DCT** fixes both problems.

### 2.2.3 The 2-D discrete cosine transform (DCT)

$$
\boxed{\;X[k,l] = c_k c_l \sum_{m=0}^{N-1}\sum_{n=0}^{N-1} x[m,n]\, \cos\frac{(2m+1)k\pi}{2N} \cos\frac{(2n+1)l\pi}{2N}\;}
$$

with normalization

$$
c_0 = \sqrt{1/N},\qquad c_k = \sqrt{2/N}\quad (k = 1, \ldots, N-1).
$$

**Inverse.**

$$
\boxed{\;x[m,n] = \sum_{k=0}^{N-1}\sum_{l=0}^{N-1} c_k c_l\, X[k,l]\, \cos\frac{(2m+1)k\pi}{2N}\cos\frac{(2n+1)l\pi}{2N}.\;}
$$

**Key properties.**

- **Real and orthogonal** — no complex arithmetic, coefficients are real.
- **Even-symmetric extension** — the DCT of an $N$-point vector is equivalent to the DFT of its **even-symmetric** extension to $2N$ points. This avoids the boundary discontinuity that plagues the DFT.
- **Excellent energy compaction** — for typical natural images, the DCT of an $8\times 8$ block concentrates 80–95% of the energy in the first 5–10 coefficients. This is why **JPEG** uses the DCT (not the DFT).
- **Computable in $O(N\log N)$** — the **fast DCT** (FDCT) is a minor variant of the FFT.

**Why the DCT and not the DFT for compression.** Two reasons:

1. **Boundary discontinuity.** An $8\times 8$ image block is not periodic, so the DFT of the block has energy in high frequencies due to the artificial jump from the last pixel to the first. The DCT of the same block has no such jump (even symmetric extension is smooth), so its high-frequency coefficients are small. Fewer high-frequency coefficients to encode = better compression.
2. **Real coefficients.** The DCT of a real image is real; the DFT is complex. Encoding real numbers is simpler and, for the DCT, more efficient in practice.

**The 8×8 DCT basis.** Each of the 64 coefficients $X[k,l]$ corresponds to a 2-D cosine basis function. The DC coefficient $X[0,0]$ is the mean; higher $(k,l)$ correspond to higher spatial frequencies in each direction. Visualize $k,l$ as "horizontal frequency" and "vertical frequency" indices — the coefficients form a 2-D array from low-frequency (top-left) to high-frequency (bottom-right).

### 2.2.4 Other sinusoidal transforms

- **DST (discrete sine transform):** uses odd-symmetric extension. Similar energy compaction to DCT; used in some audio coding.
- **Hartley:** real-valued, both DFT and DCT are special cases. Rarely used in practice.

These are primarily of historical interest; the DCT is the practical choice.

### 2.2.5 Walsh–Hadamard transform (WHT)

The **Walsh–Hadamard transform** uses **±1 basis functions** (Walsh functions).
For $N = 2^n$:

$$
\boxed{\;H_N = \frac{1}{\sqrt N}\begin{bmatrix} H_{N/2} & H_{N/2} \\ H_{N/2} & -H_{N/2}\end{bmatrix},\qquad H_2 = \frac{1}{\sqrt 2}\begin{bmatrix} 1 & 1 \\ 1 & -1\end{bmatrix}.\;}
$$

This is a **recursive** construction — $H_N$ is built by doubling $H_{N/2}$ with
sign changes. The transform is $Y = H_N X H_N^\top$ (2-D).

**Properties.**

- **Basis:** $\pm 1$ sequences — no multiplication needed for the transform itself (just sign flips and sums).
- **Fast:** $O(N\log N)$ via a **fast Walsh–Hadamard transform** (essentially the FFT but with ±1 instead of complex exponentials).
- **Not optimal for compression:** the WHT's energy compaction is worse than the DCT and KL for natural images. Its basis functions are square-wave, which don't match the smooth structure of images.
- **Used historically** in some early image codecs and as a building block for the **Haar** wavelet.

**Practical use.** The WHT is used mainly in **digital watermarking** and **cryptography** because its ±1 basis is easy to compute and invert.

### 2.2.6 The Haar transform

The **Haar transform** is the simplest **wavelet transform**. For $N = 2^n$:

$$
\boxed{\;h_0 = \frac{1}{\sqrt 2}(1, 1),\qquad h_1 = \frac{1}{\sqrt 2}(1, -1).\;}
$$

The 1-level Haar transform of an $N$-point signal applies $h_0$ and $h_1$
successively to adjacent pairs; multi-level DWT iterates on the low-pass
output. (See ECEG-6303 Chapter 5 §5.9 for the full DWT treatment.)

**2-D Haar transform.** Apply 1-D Haar to each row, then to each column. The
result is a 4-subband decomposition: **LL** (low-low), **LH** (low-high),
**HL** (high-low), **HH** (high-high) — the 2-D analogue of subband coding.

**Properties.**

- **Basis:** piecewise constant — **compact support** (2-tap filters).
- **Fast:** $O(N)$ for a 1-level transform; $O(N\log N)$ for a full multi-level DWT.
- **Decent energy compaction for piecewise-smooth images** — better than WHT, worse than DCT for natural images with strong edges.
- **Used in JPEG 2000** for its biorthogonal **CDF 9/7** wavelet — a much smoother variant of Haar.

### 2.2.7 The Slant transform

A variant of the Walsh–Hadamard that adds a **linear ramp** as its second basis
function (after the constant). The slant basis is a good match for images
with **linear gradients**.

$$
\text{Slant}_2 = \frac{1}{\sqrt 2}\begin{bmatrix} 1 & 1 \\ 1 & -1\end{bmatrix} \quad\text{(same as Haar)}, \qquad
\text{Slant}_4 = \frac{1}{2}\begin{bmatrix} 1 & 1 & 1 & 1 \\ 1 & 1 & -1 & -1 \\ 1 & -1 & -1 & 1 \\ 1 & -1 & 1 & -1\end{bmatrix}.
$$

The Slant transform is **rarely used today** — the DCT and DWT dominate. It
was historically used in some image codecs.

### 2.2.8 The Karhunen–Loève transform (KLT)

The **KL transform** (also **Hotelling transform**, or **principal component
analysis**) is the **optimal** transform in the mean-square-error sense. It
depends on the **signal statistics** — the covariance matrix $C_x$.

!!! abstract "KL transform"
    Let $C_x = \mathbb E[(\mathbf x - \mu)(\mathbf x - \mu)^\top]$ be the
    covariance of the vectorized image. Let $C_x = U \Lambda U^\top$ be its
    eigen-decomposition. Then the **KL basis** is $U$, and the transform is

    $$
    \boxed{\;\mathbf y = U^\top(\mathbf x - \mu).\;}
    $$

    The transformed coefficients have covariance $\Lambda$ — **diagonal**. The
    coefficients are **uncorrelated**, and their variances are the eigenvalues
    $\lambda_1 \ge \lambda_2 \ge \cdots \ge \lambda_N$.

**Optimality.** Among **all** linear transforms (orthogonal or not), the KL
transform minimises the MSE when the transformed coefficients are truncated
to $M < N$ components. The retained energy is $\sum_{k=1}^M \lambda_k$; the
discarded energy is $\sum_{k=M+1}^N \lambda_k$.

**The catch.** The KL basis depends on the **covariance**, which is **different
for each image** (and must be estimated from a training set). So the transform
is **data-dependent** — you have to transmit the basis to the decoder, which
defeats the purpose of compression. For this reason the KLT is **not used in
practical codecs**; it is the theoretical benchmark against which other
transforms are measured.

**Approximation.** For a **stationary** process with separable exponential
covariance $R[m,n] = \rho^{|m-n|}$ (a **1-D AR(1)** process), the KL basis
asymptotically becomes the **DCT** as $\rho \to 1$. This is the theoretical
reason the DCT is so good for natural images: natural-image autocorrelation is
roughly exponential with $\rho \approx 0.9$–$0.95$, so the DCT is nearly optimal.

### 2.2.9 Singular value decomposition (SVD)

For an $M \times N$ image matrix $A$ (rank $r \le \min(M,N)$):

$$
\boxed{\;A = U \Sigma V^\top = \sum_{k=1}^r \sigma_k\, \mathbf u_k\, \mathbf v_k^\top\;}
$$

with $U$ ($M\times M$) and $V$ ($N\times N$) orthogonal, $\Sigma$ diagonal with
**singular values** $\sigma_1 \ge \sigma_2 \ge \cdots \ge \sigma_r > 0$.

**Interpretation.** Each term $\sigma_k \mathbf u_k \mathbf v_k^\top$ is a
**rank-1 "layer"** of the image, ordered by importance ($\sigma_1$ is the most
important). The truncated SVD

$$
A_M = \sum_{k=1}^M \sigma_k \mathbf u_k \mathbf v_k^\top,\quad M < r
$$

is the **best rank-$M$ approximation** of $A$ in both the Frobenius norm and
the 2-norm (Eckart–Young theorem).

**Compression potential.** To represent $A$ by its rank-$M$ approximation,
store:

- $U_M$: $M \times M$ matrix (or just the first $M$ columns).
- $\Sigma_M$: $M$ singular values.
- $V_M$: $M \times N$ matrix (first $M$ columns of $V$).

Total parameters: $M(M + 1 + N)$ vs original $MN$. For $M \ll r$, this is
much smaller. **But** — like the KLT, the SVD is **data-dependent** (the $U, V$
factors must be transmitted), so it is **not used in practical codecs**. It is
excellent for **image denoising** and **low-rank modelling**, where the storage
overhead is acceptable.

### 2.2.10 Comparison and practical choice

| Transform | Basis | Data-dependent? | Energy compaction | Cost | Use |
|---|---|:---:|---|---|---|
| **DFT** | Complex exponentials | No | Moderate | $O(N\log N)$ | Filtering, analysis |
| **DCT** | Real cosines | No | **Excellent** | $O(N\log N)$ | **JPEG**, video codecs |
| **DST** | Real sines | No | Good | $O(N\log N)$ | Audio (rarely) |
| **WHT** | ±1 square waves | No | Moderate | $O(N\log N)$, ±1 only | Watermarking, education |
| **Haar** | Piecewise constant | No | Good for edges | $O(N\log N)$ | Fast wavelet (JPEG 2000 uses CDF 9/7) |
| **Slant** | Linear ramps | No | Moderate | $O(N\log N)$ | Historical |
| **KLT** | Eigenvectors of covariance | **Yes** | **Optimal** | $O(N^2)$ per image | Benchmark only |
| **SVD** | Singular vectors | **Yes** | **Optimal** | $O(N^3)$ | Denoising, low-rank |

**Practical choice:**

- **JPEG, MPEG, H.264** → 8×8 DCT.
- **JPEG 2000** → CDF 9/7 DWT.
- **HEVC, AV1** → DCT (and DST in some modes).
- **Denoising / low-rank restoration** → SVD or KLT (used offline).

### Worked example 2.1 — 2-D DFT of a small image and Parseval

Image $f = \begin{bmatrix} 1 & 2 \\ 3 & 4\end{bmatrix}$.

$$
F[k,l] = \sum_{m=0}^1 \sum_{n=0}^1 f[m,n] e^{-j2\pi(km/2 + ln/2)}.
$$

- $F[0,0] = 1 + 2 + 3 + 4 = 10$
- $F[0,1] = 1 + 2(-1) + 3 + 4(-1) = 1 - 2 + 3 - 4 = -2$
- $F[1,0] = 1 + 2 + 3(-1) + 4(-1) = 1 + 2 - 3 - 4 = -4$
- $F[1,1] = 1 - 2 - 3 + 4 = 0$

**Parseval.** $\frac{1}{4}(|10|^2 + |-2|^2 + |-4|^2 + |0|^2) = \frac{1}{4}(100+4+16) = 30$.
$\sum f^2 = 1 + 4 + 9 + 16 = 30$. ✓

**Interpretation.** The DC coefficient $F[0,0] = 10$ is $100/4 = 25$ times larger
than the energy of any other coefficient — this **DC concentration** is what
makes even the DFT useful for compression of smooth images.

### Worked example 2.2 — DCT of a constant 8×8 block

Block $f[m,n] = 128$ for all $m,n$.

For $k = l = 0$: $c_0 c_0 = 1/8$.

$$
X[0,0] = \frac{1}{8}\sum_{m,n} 128 \cdot 1 \cdot 1 = \frac{1}{8}\cdot 64\cdot 128 = 1024.
$$

For any $(k,l) \ne (0,0)$: one of the sums $\sum_m \cos(\ldots)$ or
$\sum_n \cos(\ldots)$ vanishes by orthogonality, so $X[k,l] = 0$.

**Result:** the DCT of a constant block is a **single nonzero coefficient** —
all 63 AC coefficients are exactly zero. This is why uniform regions compress
to almost nothing with JPEG.

The **checkerboard** block $f[m,n] = 128\cdot(-1)^{m+n}$ has the opposite
extreme: energy concentrated in $X[7,7]$ — the highest-frequency coefficient.

### Worked example 2.3 — Walsh–Hadamard transform

$x = (1, 2, 3, 4)$, $N = 4$.

$$
H_4 = \frac{1}{2}\begin{bmatrix} 1&1&1&1\\ 1&1&-1&-1\\ 1&-1&-1&1\\ 1&-1&1&-1\end{bmatrix}.
$$

$y = H_4 x$:

$$
y_0 = \tfrac12(1+2+3+4) = 5,\quad
y_1 = \tfrac12(1+2-3-4) = -2,
$$
$$
y_2 = \tfrac12(1-2-3+4) = 0,\quad
y_3 = \tfrac12(1-2+3-4) = -1.
$$

**Energy check.** $25 + 4 + 0 + 1 = 30$. Original: $1 + 4 + 9 + 16 = 30$. ✓

The signal (increasing ramp) concentrates in $y_0$ — the WHT doesn't capture
the smooth trend as well as the DCT would, but it is **cheap** (all operations
are adds/subtracts).

### Worked example 2.4 — Haar transform (1-D, 2-level)

$x = (4, 6, 10, 12)$.

**Level 1.**

$$
a_1 = \frac{1}{\sqrt 2}(x[0]+x[1],\ x[2]+x[3]) = (10/\sqrt 2,\ 22/\sqrt 2) \approx (7.07,\ 15.56),
$$
$$
d_1 = \frac{1}{\sqrt 2}(x[0]-x[1],\ x[2]-x[3]) = (-2/\sqrt 2,\ -2/\sqrt 2) \approx (-1.41,\ -1.41).
$$

**Level 2.**

$$
a_2 = \frac{1}{\sqrt 2}(a_1[0] + a_1[1]) = \frac{10+22}{2} = 16,
$$
$$
d_2 = \frac{1}{\sqrt 2}(a_1[0] - a_1[1]) = \frac{10-22}{2} = -6.
$$

**Final DWT coefficients:** $(a_2, d_2, d_1[0], d_1[1]) = (16, -6, -1.41, -1.41)$.

**Reconstruction:**

$$
a_1 = \frac{1}{\sqrt 2}(a_2 + d_2,\ a_2 - d_2) = (10/\sqrt 2,\ 22/\sqrt 2)
\;\rightarrow\;
x[0] = \frac{a_1[0]+d_1[0]}{\sqrt 2} = 4, \ldots
$$

Exact. ✓

**Energy compaction.** The largest coefficient is $a_2 = 16$ (energy $256$);
original total energy is $4^2+6^2+10^2+12^2 = 296$. So **86% of the energy** is
in **one** coefficient.

### Worked example 2.5 — KL transform on a 2×2 correlated signal

Covariance $C_x = \begin{bmatrix} 1 & 0.9 \\ 0.9 & 1\end{bmatrix}$.

Eigen-decomposition:
$$
C_x = \underbrace{\begin{bmatrix} 1/\sqrt 2 & -1/\sqrt 2 \\ 1/\sqrt 2 & 1/\sqrt 2\end{bmatrix}}_{U}
\underbrace{\begin{bmatrix} 1.9 & 0 \\ 0 & 0.1\end{bmatrix}}_{\Lambda}
\underbrace{\begin{bmatrix} 1/\sqrt 2 & 1/\sqrt 2 \\ -1/\sqrt 2 & 1/\sqrt 2\end{bmatrix}}_{U^\top}.
$$

**KL coefficients' variances:** $\lambda_1 = 1.9$, $\lambda_2 = 0.1$.

**Energy compaction.** $\frac{1.9}{1.9+0.1} = 95\%$ of the variance is in the
first component. Keeping only $y_1$ captures 95% of the signal energy — a 2:1
compression with only 5% MSE loss.

**Comparison with the DCT.** The DCT of the same signal (orthonormal DCT-2)
also concentrates energy in the first coefficient but less efficiently — for
$\rho = 0.9$, the DCT captures about 92% (versus 95% for the optimal KLT). The
gap shrinks as $\rho \to 1$; this is the theoretical basis for DCT's dominance.

### Worked example 2.6 — SVD of a small image

$A = \begin{bmatrix} 1 & 2 & 3 \\ 4 & 5 & 6 \\ 7 & 8 & 9\end{bmatrix}$.

Compute $A^\top A = \begin{bmatrix} 66 & 78 & 90 \\ 78 & 93 & 108 \\ 90 & 108 & 126\end{bmatrix}$.

Eigenvalues: $272.6, 12.4, 0$ (rank 2). Singular values:
$\sigma_1 \approx 16.51, \sigma_2 \approx 3.52, \sigma_3 = 0$.

**Best rank-1 approximation.**

$$
A_1 = \sigma_1 \mathbf u_1 \mathbf v_1^\top
$$

Numerically, $\mathbf u_1 \approx (0.23, 0.53, 0.81)$, $\mathbf v_1 \approx (0.48, 0.55, 0.62)$.

$$
A_1 \approx 16.51\begin{bmatrix} 0.11 & 0.13 & 0.14 \\ 0.26 & 0.29 & 0.33 \\ 0.39 & 0.44 & 0.50\end{bmatrix}
\approx \begin{bmatrix} 1.82 & 2.10 & 2.36 \\ 4.28 & 4.82 & 5.42 \\ 6.72 & 7.53 & 8.49\end{bmatrix}.
$$

**Approximation error:** $\|A - A_1\|_F = \sqrt{\sigma_2^2 + \sigma_3^2} = 3.52$.
Relative error: $3.52/\sqrt{272.6+12.4} = 3.52/16.88 = 20.8\%$.

**Energy captured by 1 coefficient:** $\sigma_1^2/\sum \sigma_i^2 = 272.6/285 = 95.6\%$.
The rank-1 approximation captures 95.6% of the Frobenius energy with 1/3 of the
rank — a strong "compression" but, as noted, it requires transmitting $\mathbf u_1, \mathbf v_1$.

---

## 2.3 Image enhancement

Enhancement is **subjective** — "better" is defined by the observer or the task,
not by an exact mathematical criterion. The goal is to make features of interest
more visible: contrast, edges, structure.

### 2.3.1 Point operations

A **point operation** modifies each pixel based solely on its value:
$g[m,n] = T(f[m,n])$.

| Operation | Formula | Effect |
|---|---|---|
| **Identity** | $g = f$ | — |
| **Brightness offset** | $g = f + b$ | Shifts the intensity |
| **Contrast scaling** | $g = a f$ | Amplifies / reduces contrast |
| **Clipping** | $g = \min(255, \max(0, f))$ | Prevents overflow |
| **Threshold** | $g = 255$ if $f > T$, else $0$ | Binarization |
| **Contrast stretch** | $g = 255 \cdot \frac{f - f_{\min}}{f_{\max} - f_{\min}}$ | Full dynamic range |
| **Log** | $g = c \log(1 + f)$ | Compresses dynamic range |
| **Gamma** | $g = 255 (f/255)^\gamma$ | $\gamma < 1$ brightens, $\gamma > 1$ darkens |
| **Density slicing** | $g$ piecewise-linear in $f$ | Highlights a specific intensity range |

**Contrast stretching** is the simplest useful operation — extend the intensity
range so it occupies the full $[0, L-1]$ interval.

### 2.3.2 Histograms

The **histogram** of an image is the frequency of each intensity level:
$h(k) = $ number of pixels with value $k$. The normalized histogram
$p(k) = h(k)/N$ is the **empirical PDF**.

**Histogram equalization.** Choose a point operation that makes the histogram
**approximately uniform** — maximising visual contrast without information
loss.

!!! abstract "Histogram equalization"
    For an $L$-level image with CDF

    $$
    C(k) = \sum_{j=0}^{k} p(j)
    $$

    the equalization map is

    $$
    \boxed{\;s_k = \lfloor (L-1)\, C(k) \rfloor.\;}
    $$

    This is the **discrete** version of the continuous mapping $s = (L-1) C(r)$.
    It is a **monotone** transformation (preserves pixel ordering).

**Properties.**

- **Automatic** — no parameters to tune.
- **Global** — uses the full-image histogram. Fails when the image has locally uniform regions (part of the image may be over- or under-enhanced).
- **Produces the most uniform histogram** possible with discrete levels (not perfectly uniform because of quantization).

**Histogram modification.** Replace the target uniform distribution with an
arbitrary target PDF $p_{\text{target}}(k)$. Compute the target CDF $C_{\text{target}}$
and map:

$$
s = C_{\text{target}}^{-1}\bigl(C_{\text{in}}(r)\bigr).
$$

This is **histogram specification** (or **matching**) — the most general
histogram-based point operation.

**Applications.**

- **Medical imaging:** enhance contrast in X-ray / CT without windowing.
- **Photography:** recover detail in under-/over-exposed images (equalization is the global version of the "Auto-Levels" tool).
- **Remote sensing:** balance illumination across a satellite scene.
- **Video:** stabilise brightness across frames (per-frame equalization is often too aggressive; use adaptive / CLAHE).

### Worked example 2.7 — Histogram equalization by hand

$L = 4$ gray levels. Histogram of a 10-pixel image: $n_0 = 4, n_1 = 3, n_2 = 2, n_3 = 1$.

**CDF.**

$$
C(0) = 0.4,\ C(1) = 0.7,\ C(2) = 0.9,\ C(3) = 1.0.
$$

**Mapping** $s_k = \lfloor 3 \cdot C(k) \rfloor$:

$$
s_0 = \lfloor 1.2 \rfloor = 1,\quad
s_1 = \lfloor 2.1 \rfloor = 2,\quad
s_2 = \lfloor 2.7 \rfloor = 2,\quad
s_3 = \lfloor 3.0 \rfloor = 3.
$$

**New histogram:**

| New level | Pixels | Count |
|:---:|---|:---:|
| 0 | — | 0 |
| 1 | from level 0 | 4 |
| 2 | from levels 1 and 2 | 5 |
| 3 | from level 3 | 1 |

The new histogram is roughly uniform except for level 0 (empty). The
**dynamic range is spread** from 4 levels to 3 occupied levels — for this
small example, the entropy has increased from $\sim 1.85$ bits to $\sim 1.36$ bits
(the original was more concentrated than the new one, hmm — let me recheck:
original entropy = $-0.4\log 0.4 - 0.3\log 0.3 - 0.2\log 0.2 - 0.1\log 0.1 ≈ 1.85$;
new = $-0.4\log 0.4 - 0.5\log 0.5 - 0.1\log 0.1 ≈ 1.36$). Interesting — the new
histogram has **lower** entropy in this case because the mapping collapsed two
levels into one.

**In practice**, histograms of real images are much more skewed (heavily
peaked near dark or bright), and equalization increases entropy. The small
example above is atypical.

### Worked example 2.8 — Histogram specification

Input $L = 4$ with histogram $n = (4, 3, 2, 1)$. Target PDF: $p_{\text{target}} = (0.25, 0.25, 0.25, 0.25)$.

Target CDF: $C_{\text{target}} = (0.25, 0.50, 0.75, 1.0)$.

Inverse: $C_{\text{target}}^{-1}(0.25) = 0$, $C^{-1}(0.50) = 1$, $C^{-1}(0.75) = 2$, $C^{-1}(1.0) = 3$.

Input CDF: $C_{\text{in}} = (0.4, 0.7, 0.9, 1.0)$.

**Mapping:** $s_k = C_{\text{target}}^{-1}(C_{\text{in}}(k))$.

$$
s_0 = C_{\text{target}}^{-1}(0.4) = 1,\quad
s_1 = C^{-1}(0.7) = 2,\quad
s_2 = C^{-1}(0.9) = 3,\quad
s_3 = C^{-1}(1.0) = 3.
$$

New histogram: level 0 (0), level 1 (4), level 2 (3), level 3 (3). More
uniform than the original (4, 3, 2, 1), but not exactly the target (2.5, 2.5,
2.5, 2.5) — the discrete quantization prevents perfection. The algorithm is
still the right one.

### 2.3.3 Spatial filtering

**Spatial filtering** is convolution with a small kernel:

$$
\boxed{\;g[m,n] = \sum_{k=-K}^{K}\sum_{l=-L}^{L} w[k,l]\, f[m-k, n-l].\;}
$$

Categories:

| Type | Kernel | Effect |
|---|---|---|
| **Averaging (box)** | Uniform $1/9$ (3×3) | Blurs, reduces noise |
| **Gaussian** | $\propto e^{-(m^2+n^2)/(2\sigma^2)}$ | Blurs, smoother than box |
| **Lowpass** | Any kernel with positive sum | Blurs |
| **Highpass** | $\delta - \text{lowpass}$ | Sharpens, enhances edges |
| **Bandpass** | Difference of two lowpass | Enhances textures in a frequency range |
| **Directional** | Oriented kernel | Enhances edges in a specific direction |
| **Laplacian** | $\begin{bmatrix} 0&-1&0\\-1&4&-1\\0&-1&0\end{bmatrix}$ | Detects edges (2nd derivative) |
| **Sobel X** | $\begin{bmatrix}-1&0&1\\-2&0&2\\-1&0&1\end{bmatrix}$ | Horizontal gradient |
| **Sobel Y** | $\begin{bmatrix}-1&-2&-1\\0&0&0\\1&2&1\end{bmatrix}$ | Vertical gradient |

**Boundary handling.** Padding choices:

- **Zero padding:** assume $f = 0$ outside. Introduces dark borders.
- **Replicate:** extend with edge values. Preserves brightness.
- **Mirror:** reflect at boundaries. Smooth.
- **Wrap:** circular. Good for FFT-based filtering.

### 2.3.4 Median filtering

Replace each pixel with the **median** of its $K \times K$ neighborhood.
**Nonlinear** — does not commute with sum, cannot be written as a convolution.

**Properties.**

- **Excellent for impulse (salt-and-pepper) noise** — an isolated outlier is
  eliminated (it's not the median unless the majority of the neighborhood is
  also outliers).
- **Preserves edges** — because the median of a step neighbourhood still lies
  on the step, whereas averaging would blur it.
- **Slower than convolution** — requires sorting (or selection) per pixel.
- **Not a linear filter** — frequency-domain analysis does not apply directly.

### Worked example 2.9 — median vs mean on impulse noise

3×3 patch centered on a salt pixel ($255$):

$$
P = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 255 & 0 \\ 0 & 0 & 0\end{bmatrix}.
$$

**Mean filter (3×3 uniform):** output = $\frac{1}{9}(8 \cdot 0 + 255) = 28.3$.

**Median filter:** sort the 9 values: $\{0,0,0,0,0,0,0,0,255\}$. Median = 0.

**Conclusion.** The mean filter leaves a residual "halo" of $28.3$ at each
salt pixel — the noise is spread over the whole neighbourhood. The median
filter **exactly removes** isolated impulses — the center pixel is restored to 0.

This is why **median filtering is the default for impulse-noise removal**;
averaging is the default for Gaussian-noise removal.

**Weaknesses.** Median filtering destroys fine textures and thin lines (line
thinner than $K/2$ pixels is removed). For these, use **weighted median**,
**stack filters**, or **bilateral filtering** (edge-preserving smoothing).

### 2.3.5 Edge enhancement — 2-D FIR and IIR

**The Laplacian.** The 2-D second derivative:

$$
\boxed{\;\nabla^2 f = \frac{\partial^2 f}{\partial x^2} + \frac{\partial^2 f}{\partial y^2} \approx f[m+1,n] + f[m-1,n] + f[m,n+1] + f[m,n-1] - 4 f[m,n].\;}
$$

Kernel: $\begin{bmatrix}0&1&0\\1&-4&1\\0&1&0\end{bmatrix}$. Zero response on constant
regions; peaks at edges. **Edge enhancement** by adding the Laplacian:

$$
\boxed{\;g = f - \nabla^2 f\quad\text{(2nd-derivative sharpening)}.\;}
$$

**Unsharp masking.** A more robust edge enhancement:

$$
\boxed{\;g = f + k\,(f - \text{lowpass}(f)) = (1+k) f - k \cdot \text{lowpass}(f).\;}
$$

The **highpass** component $f - \text{lowpass}(f)$ contains the edges; adding it
amplifies them. The parameter $k$ controls the amount of sharpening.

**High-boost filtering.** A generalisation:

$$
g = A f - \text{lowpass}(f),\quad A \ge 1.
$$

$A = 1$ gives unsharp masking; $A > 1$ adds the original image back — deeper
shadows / brighter highlights.

**2-D FIR filters.** Small kernels (3×3, 5×5) implemented as direct
convolution. Common edge detectors:

- **Sobel** — good gradient approximation, less noisy than simple differences.
- **Prewitt** — same structure, different weighting.
- **Roberts** — 2×2 diagonal kernels; noisier but faster.
- **Canny** — Gaussian smoothing + gradient + non-max suppression + hysteresis; the best general-purpose edge detector (not a simple convolution).

**2-D IIR filters.** Larger kernels implemented as recursive filters —
e.g., a Gaussian approximated by a cascade of first-order IIR filters. Used
when the kernel is too large for direct FIR (e.g., $\sigma \ge 5$).

### Worked example 2.10 — Unsharp masking by hand

1-D signal $f = (10, 10, 10, 50, 50, 50)$ — step edge.

3-point box lowpass: $f_{\text{lp}} = (10, 10, 23.3, 36.7, 50, 50)$.

Highpass $f - f_{\text{lp}} = (0, 0, -13.3, 13.3, 0, 0)$.

**Unsharp masking** with $k = 1$:

$$
g = f + (f - f_{\text{lp}}) = (10, 10, -3.3, 63.3, 50, 50).
$$

Clipped to $[0, 255]$: same (no clipping needed). The **edges are now overshoot
and undershoot** — the rise from 10 to 50 is amplified to a rise from
$-3.3$ to $63.3$, appearing **sharper** to the eye because the derivative is
larger at the edge.

**Exaggerating $k = 2$:**

$$
g = (10, 10, -16.7, 76.7, 50, 50).
$$

Overshoot is now $-16.7$ and $76.7$ — more dramatic, but risks clipping and
visible "halo" artifacts.

### 2.3.6 Homomorphic filtering and cepstrum

**The problem.** Illumination is multiplicative: $f[m,n] = i[m,n] \cdot r[m,n]$
where $i$ is the illumination (slowly varying, low-frequency) and $r$ is the
reflectance (rapidly varying, high-frequency). Standard linear filters
operate on **additive** signals — they distort multiplicative ones.

**The homomorphic trick.**

1. **Log:** $\log f = \log i + \log r$ — converts multiplicative to additive.
2. **Linear filter** on $\log f$: highpass to boost reflectance and suppress
   illumination variations.
3. **Exp:** $\hat f = e^{\text{filtered}}$.

The filter is typically a **highpass** (attenuates the low-frequency
illumination, passes the high-frequency reflectance) — the opposite of
lowpass smoothing. This is **homomorphic filtering**.

**Cepstrum.** Take the Fourier transform of the log:

$$
\boxed{\;\text{cepstrum}(x) = \mathcal F^{-1}\{\log|\mathcal F\{x\}|\}.\;}
$$

The cepstrum separates the **source** (periodic) from the **filter**
(smooth envelope) — used in speech pitch estimation (ECEG-6303 Ch 4) and
image deconvolution (restoring an image blurred by an unknown PSF, when the
PSF is "narrow" in the cepstral domain).

### 2.3.7 Colour image enhancement

Colour enhancement is trickier because it must preserve **hue** while modifying
**brightness** or **saturation**. Two strategies:

**Per-channel processing.** Apply the enhancement to R, G, B separately.
Pros: simple. Cons: changes hue — a bright red pixel may become pink.

**Decoupled processing.** Convert to YCbCr (or HSV), enhance only Y (luminance),
convert back. Pros: preserves hue. Cons: does not adjust saturation.

For most applications, **Y-enhancement** (YCbCr or HSV-V) is the standard
approach.

---

## 2.4 Image restoration

Restoration is **objective** — the goal is to recover the **original** image
from a known or estimated degradation. Contrast with enhancement, which is
subjective.

### 2.4.1 Degradation model

$$
\boxed{\;g[m,n] = (h * f)[m,n] + \nu[m,n]\;}
$$

where $f$ is the original image, $h$ is the **point-spread function (PSF)**
(blur), $\nu$ is additive noise, $g$ is the observed image.

In the frequency domain:

$$
\boxed{\;G(u, v) = H(u,v)\, F(u,v) + N(u,v).\;}
$$

**Common degradation sources.**

- **Motion blur:** $h$ is a 1-D box along the motion direction (or a triangle for accelerating motion).
- **Out-of-focus blur:** $h$ is a disk of radius $R$ — a "pillbox" PSF.
- **Atmospheric turbulence:** $h$ is a Gaussian with $\sigma \propto$ distance.
- **Sensor noise:** additive Gaussian, Poisson (shot noise), or salt-and-pepper.

### 2.4.2 Inverse filter

The naive restoration: divide by the PSF:

$$
\boxed{\;\hat F(u,v) = \frac{G(u,v)}{H(u,v)} = F(u,v) + \frac{N(u,v)}{H(u,v)}.\;}
$$

**The fatal flaw.** Wherever $H(u,v)$ is small or zero, the noise term
$N/H$ **blows up**. Even a tiny amount of noise (which is always present)
amplifies enormously in those frequency bins — the restored image is dominated
by noise. **Inverse filtering fails in practice** unless (i) the noise is zero
(unrealistic) or (ii) $H$ is bounded away from zero and flat (rare).

### Worked example 2.11 — Inverse filter failure

Blur $H(f) = e^{-(f/f_0)^2}$ (Gaussian). Noise $\nu \sim \mathcal N(0, \sigma^2)$,
$\sigma = 10^{-3}$, image amplitude $\sim 1$. So **input SNR = 60 dB** — excellent.

Inverse filter at $f = 2f_0$: $H = e^{-4} \approx 0.018$. Amplified noise:
$\sigma/H = 10^{-3}/0.018 \approx 0.056$ — **noise is now 5% of the signal,
SNR drops to 25 dB**.

At $f = 3f_0$: $H = e^{-9} \approx 10^{-4}$. Noise $10^{-3}/10^{-4} = 10$ —
**noise is 10× the signal. Image unusable.**

**The Gaussian's tails** are what kill the inverse filter — they decay
faster than the noise floor. Any real blur has tails, and any real noise has a
finite level; the ratio is a disaster in the far tails.

### 2.4.3 Wiener filter (minimum MSE)

The **Wiener filter** minimises the **mean-square error** between $\hat f$ and
the true $f$, given the statistics of the signal and noise:

$$
\boxed{\;W(u,v) = \frac{H^*(u,v)\, S_f(u,v)}{|H(u,v)|^2 S_f(u,v) + S_\nu(u,v)}.\;}
$$

Equivalently:

$$
\boxed{\;W(u,v) = \frac{H^*(u,v)}{|H(u,v)|^2 + S_\nu(u,v)/S_f(u,v)}.\;}
$$

Here $S_f = |F|^2$ is the **signal power spectrum** and $S_\nu = |N|^2$ is the
**noise power spectrum**.

**Interpretation.**

- Where SNR is **high** ($S_f \gg S_\nu$): $W \approx 1/H$ — **inverse filter**, perfect restoration.
- Where SNR is **low** ($S_f \ll S_\nu$): $W \approx 0$ — **attenuate**, avoid amplifying noise.
- The filter **smoothly interpolates** between these regimes.

**Practical form.** Replace $S_f, S_\nu$ by constants (often $S_\nu/S_f \approx K$
estimated from a noise-only region) to get the **approximate Wiener filter**:

$$
W(u,v) = \frac{H^*(u,v)}{|H(u,v)|^2 + K}.
$$

The scalar $K$ is the **regularisation parameter** — it prevents division by
small $H$.

### Worked example 2.12 — Wiener filter

Blur $H(f) = e^{-(f/f_0)^2}$ (same as 2.11). Signal PSD $S_f = |F(f)|^2$ —
typical natural-image power spectrum $S_f \propto 1/f^2$. Noise PSD
$S_\nu = 10^{-6}$ (60 dB SNR).

**At low frequency** ($f = 0.1 f_0$): $H \approx 1$, SNR $\approx 10^6$,
$W \approx 1$. Noise not amplified.

**At high frequency** ($f = 3 f_0$): $H = 10^{-4}$, $H^2 = 10^{-8}$,
$S_f \propto 1/9 \approx 0.11$, $S_\nu/S_f = 10^{-6}/0.11 \approx 10^{-5}$.
Then $W = H/(H^2 + 10^{-5}) = 10^{-4}/(10^{-8} + 10^{-5}) \approx 10^{-4}/10^{-5} = 10$.

**Output at high frequency:** $W \cdot G = 10 \cdot (H F + N) = 10 \cdot (10^{-4} F + N) = 10^{-3} F + 10 N$.
The noise contribution is amplified by $10$ (vs. $10^4$ for inverse) — **1000×
better**. But the signal is still heavily attenuated; the Wiener filter
**accepts** this because the SNR there is too low to recover.

**MSE.** The Wiener filter's MSE is the minimum among all linear filters for
this problem — the theoretical optimum.

### 2.4.4 Geometric mean filter

Interpolates between the **inverse** and **Wiener** filters:

$$
\boxed{\;W_{\alpha}(u,v) = \frac{H^*(u,v)}{|H(u,v)|^2 + \alpha\, S_\nu/S_f},\quad \alpha \in [0,1].\;}
$$

- $\alpha = 1$ → Wiener filter (minimum MSE).
- $\alpha = 0$ → inverse filter (perfect restoration, worst noise).
- $\alpha \in (0,1)$ → trade-off.

The name comes from the **geometric mean** form used in some texts: $W_\alpha = \left[\frac{1}{H}\right]^\alpha \left[\frac{H^*}{|H|^2 + S_\nu/S_f}\right]^{1-\alpha}$
(though the exact form varies). The general principle is the same:
**one knob trades restoration fidelity against noise amplification**.

### 2.4.5 Constrained least squares (CLS)

When the noise PSD is unknown, use a **smoothness prior** instead of the
Wiener's explicit noise model:

$$
\boxed{\;\hat F = \arg\min_F\, \|G - H F\|^2 + \gamma\, \|P F\|^2.\;}
$$

Here $P$ is a **highpass operator** (e.g., the Laplacian), so $\|PF\|^2$
measures the **roughness** of the restored image. The solution:

$$
\boxed{\;\hat F = \frac{H^*}{|H|^2 + \gamma |P|^2}\, G.\;}
$$

**Comparison with Wiener.**

- Wiener uses **statistical** priors ($S_f, S_\nu$).
- CLS uses **deterministic** priors (smoothness).
- Both have a **regularisation parameter** ($K$ or $\gamma$) that controls the
  fidelity/smoothness trade-off.
- **CLS is often preferred** when the noise level is unknown but smoothness is
  a reasonable prior — as in astronomical imaging.

### Worked example 2.13 — CLS restoration

Same blur $H(f) = e^{-(f/f_0)^2}$. Use Laplacian $P(f) = -4\pi^2 f^2$ (in the
Fourier domain), $\gamma = 0.01$.

**At low frequency** ($f = 0.1 f_0$): $P \approx 0.4 f_0^2$, $\gamma P^2 \approx 0.0016 f_0^4$.
This is small compared to $H^2 \approx 1$, so $W \approx H^* \approx 1$ — full restoration.

**At high frequency** ($f = 3 f_0$): $H^2 \approx 10^{-8}$, $P^2 \approx 36 f_0^4$,
$\gamma P^2 = 0.36 f_0^4$. This **dominates** $H^2$, so $W \approx H^*/(0.36 f_0^4) \approx 3 \cdot 10^{-4}$
— the filter **cuts off** the high frequencies before they can amplify noise.

**MSE vs. Wiener.** CLS is slightly worse than Wiener when the true
$S_f, S_\nu$ are known, but **more robust** when they are not. The regularization
parameter can be tuned to match the visual quality.

### 2.4.6 Nonlinear filters

For non-Gaussian noise or edge-preserving smoothing:

- **Median filter** — best for salt-and-pepper.
- **Bilateral filter** — edge-preserving smoothing: averages nearby pixels weighted by **both spatial distance and intensity similarity**. Preserves edges while removing noise.
- **Non-local means** — averages similar patches anywhere in the image. Excellent for Gaussian noise; slow.
- **BM3D** — block-matching + 3D filtering; the state-of-the-art denoiser for many years.
- **Total variation** — minimises the total variation $\int |\nabla f|$ subject to data fidelity. Excellent for piecewise-smooth images (medical, astronomical).

### 2.4.7 Smoothing splines and interpolation

For **inpainting** (filling missing regions) or **super-resolution**:

- **Bicubic interpolation** — smooth, 3rd-order polynomial, 16-tap. The default for image resizing.
- **Spline interpolation** — piecewise polynomial with continuous derivatives. MATLAB `spline`, `interp2`.
- **Thin-plate spline** — 2-D analogue of the cubic spline; minimises bending energy. Used in registration.
- **Smoothing spline** — fits a spline that trades off data fidelity vs. curvature penalty. A 2-D version of CLS.

---

## 2.5 Image compression

### 2.5.1 Data rates and redundancy

**Raw image data rate** for an $N \times M$ image with $b$ bits per channel,
$c$ channels, $f$ frames/second:

$$
\text{bitrate} = N \cdot M \cdot b \cdot c \cdot f\quad\text{bits/s.}
$$

**Example.** 1080p video (1920×1080), 24 fps, 8-bit RGB = $1920 \times 1080 \times 24 \times 24 \approx 1.2$ Gbps raw. Broadcast needs $\approx 20$ Mbps — a **60:1 compression**.

**Three types of redundancy.**

| Redundancy | Explanation | Reduction method |
|---|---|---|
| **Coding** | Fixed-length codes inefficient for non-uniform symbols | Huffman, arithmetic coding |
| **Spatial** | Neighbouring pixels correlated | Prediction, transform coding |
| **Psychovisual** | Eye is insensitive to certain detail | Perceptual quantization (quantization matrix) |

### 2.5.2 Predictive coding (DPCM)

**Differential pulse-code modulation (DPCM).** Predict each pixel from its
neighbours, code only the prediction error:

$$
\boxed{\;d[m,n] = f[m,n] - \hat f[m,n],\quad \hat f[m,n] = \sum_{(k,l) \in \mathcal N} a_{k,l} f[m-k, n-l].\;}
$$

The **prediction error** $d$ has a much smaller range and variance than $f$, so
it can be encoded with fewer bits. The decoder reconstructs:

$$
f[m,n] = d[m,n] + \hat f[m,n].
$$

**Predictors.**

- **1-D (previous pixel):** $\hat f = f[m-1, n]$.
- **2-D (causal neighbours):** $\hat f = a_1 f[m-1,n] + a_2 f[m,n-1] + a_3 f[m-1,n-1]$.
- **Optimal linear predictor:** solve the Yule–Walker equations for the
  autocovariance (ECEG-6303 Ch 3 §3.8) — but for 2-D, the predictor becomes a
  small FIR filter.

**Limitations.**

- **Error propagation:** the decoder reconstructs $f[m,n]$ from past **reconstructed** pixels — a single error can propagate.
- **Fixed code length** for each residual — no entropy coding (a simple DPCM system).
- **Poor for textures:** textures have little spatial correlation → prediction gain small.

**Practical use.** DPCM is used in **lossless JPEG** and as the **DC coefficient coder** in JPEG (adjacent blocks have correlated DCs).

### Worked example 2.14 — 1-D DPCM

$x = (10, 12, 14, 15, 15, 14, 12, 10)$ — a smooth ramp up and down.

**Predictor:** $\hat x[n] = x[n-1]$.

**Residuals:** $d[0] = 10$, $d[1] = 2$, $d[2] = 2$, $d[3] = 1$, $d[4] = 0$, $d[5] = -1$, $d[6] = -2$, $d[7] = -2$.

**Variance reduction:**

- Signal variance: $\approx 4.2$.
- Residual variance (excluding $d[0]$): $\approx 2.3$.

**Entropy.** Signal range 10–15 (6 levels): $\log_2 6 \approx 2.6$ bits. Residual range $-2$ to $+2$ (5 levels): $\log_2 5 \approx 2.3$ bits. Modest reduction for this example; in real images with larger dynamic range, the reduction is much larger (8 bits → 3–5 bits typical).

**Better predictor** (2-tap): $\hat x[n] = 2x[n-1] - x[n-2]$ (linear extrapolation).

Residuals: $d[0]=10$, $d[1]=2$, $d[2]=0$, $d[3]=-1$, $d[4]=-1$, $d[5]=-1$, $d[6]=-1$, $d[7]=0$. Even smaller — but slightly larger than one might hope because the corners of the ramp break the linear model.

### 2.5.3 Transform coding

**The idea.** Apply an invertible transform $T$ to a block, quantize the coefficients
(non-uniformly, more heavily for high frequency where the eye is less sensitive),
and entropy-code. The block-based **JPEG** uses an 8×8 DCT.

**JPEG pipeline** (for reference; JPEG is covered fully in most textbooks):

1. Split into 8×8 blocks.
2. Convert RGB → YCbCr, subsample chroma 4:2:0.
3. Compute 8×8 DCT of each block.
4. **Quantize** using the standard 8×8 quantization table, scaled by a quality factor $Q \in [1, 100]$:
   $$
   \hat X[k,l] = \text{round}\!\left(\frac{X[k,l]}{Q(k,l)}\right).
   $$
   The matrix $Q(k,l)$ is **larger at high frequencies** (more aggressive quantization) — perceptual.
5. **Zigzag scan** the coefficients (low to high frequency) → long runs of zeros at the end.
6. **Entropy code** (Huffman) the run-length encoded zeros + the nonzero coefficients.

**Decoding:** reverse each step. Because quantization is lossy, the reconstruction
is **not identical** to the original — but at $Q = 75$, the result is
perceptually transparent for most images.

**JPEG quality levels.**

| Quality | Compression | Visual quality | Artifacts |
|:---:|:---:|---|---|
| 95 | 10:1 | Transparent | None |
| 75 | 20:1 | Good | None visible |
| 50 | 40:1 | Fair | Slight blocking |
| 25 | 80:1 | Poor | Visible blocking, ringing |
| 10 | 200:1 | Very poor | Severe blocking, blur |

### Worked example 2.15 — DCT block coding

An 8×8 block has DCT coefficients concentrated at (0,0) = 1000, (0,1) = 500,
(1,0) = 300, (1,1) = 200, all others near zero.

**Quantization** with a Q matrix that has $Q(0,0) = 16$ (low frequency, gentle),
$Q(0,1) = Q(1,0) = 32$, $Q(1,1) = 64$, and $Q(k,l) = 128$ for the rest.

Quantized values: $1000/16 = 63$, $500/32 = 16$, $300/32 = 9$, $200/64 = 3$,
all others 0 (since the original was near zero and the quantizer rounds to 0).

**Stored: 4 nonzero coefficients** + run-length encoding that covers the
remaining 60 zeros in ~3–4 bytes. Total: ~12 bytes per block — a **42:1
compression** from 8×8 8-bit = 64 bytes.

**Reconstruction:** dequantize (multiply by Q) and inverse DCT. The block is
recovered to within the quantization error of the 4 largest coefficients —
perceptually identical to the original for a smooth block.

### 2.5.4 Vector quantization (VQ)

**Idea.** Group pixels into **vectors** (e.g., 4×4 blocks) and encode each vector
by the **index** of its nearest codeword in a **codebook**.

$$
\boxed{\;\text{codebook} = \{\mathbf c_1, \dots, \mathbf c_K\},\quad \text{encode}(\mathbf x) = \arg\min_k \|\mathbf x - \mathbf c_k\|.\;}
$$

**Compression.** For 4×4 = 16-pixel vectors and a codebook of $K = 256$ codewords:
$16 \times 8 = 128$ bits per block → 8 bits per block = **16:1** compression.

**Codebook design.** **LBG algorithm** (Linde–Buzo–Gray, 1980) — iteratively
refine centroids, like k-means. Optimality = **Lloyd–Max** condition in vector
form: each codeword is the centroid of its Voronoi region.

**Properties.**

- **Optimal** among all fixed-rate vector coders (for the training distribution).
- **Asymptotically** approaches the rate-distortion bound.
- **The catch:** the codebook must be transmitted. For $K = 256$ codewords and 16 pixels: $256 \times 16$ bytes = 4 KB overhead — negligible for large images, significant for small ones.

**Not used in modern codecs** for images (transform coding dominates), but used in:
- **Speech coding** (as the "excitation" quantizer in CELP).
- **Deep learning** (as a product quantizer in face recognition).
- **Image retrieval** (as a compressed feature descriptor).

### 2.5.5 Block truncation coding (BTC)

A fast, low-complexity block-based coder. For each 4×4 or 8×8 block:

1. Compute the **mean** $\bar x$ and **standard deviation** $\sigma_x$ of the block.
2. **Binarize** at $\bar x$: pixel is 1 if above mean, 0 if below.
3. **Two reconstruction levels** $a, b$ chosen so the block mean and variance are preserved:
   $$
   a = \bar x - \sigma_x \sqrt{q/m},\quad b = \bar x + \sigma_x \sqrt{m/q},
   $$
   where $m$ = number of 0s, $q$ = number of 1s, $m + q = N$.
4. Store $\bar x$, $\sigma_x$, and the **mask** (1 bit per pixel).

**Rate.** For a 4×4 block: 4×4 = 16 mask bits + 8 mean + 8 std = 32 bits. Original
= $16 \times 8 = 128$ bits. Compression = **4:1** (with quantization of mean
and std, better).

**Properties.**

- **Very fast** — no transform, no floating point (with quantized mean/std).
- **Decent quality** on natural images — comparable to low-quality JPEG.
- **Historical:** pre-JPEG era; still used in some embedded hardware.

### Worked example 2.16 — BTC on a 4×4 block

Block values (8-bit): $\{100, 120, 140, 160\}$ in a checkerboard pattern.

Actually, do a simpler one: block with 8 pixels = 100, 8 pixels = 200.

**Mean:** $\bar x = 150$. **Variance:** $\sigma_x = 50$. **Mask:** 1 for the 200s, 0 for the 100s.

**Reconstruction levels.** $m = q = 8$, $N = 16$. So $\sqrt{q/m} = 1$, $\sqrt{m/q} = 1$.
$a = 150 - 50 = 100$, $b = 150 + 50 = 200$. Exact! (Because the block has only two distinct values, BTC reconstructs perfectly.)

**Bits stored:** 16 mask + 8 mean + 8 std = 32 bits (vs. 16 pixels × 8 bits = 128
bits). **4:1 compression, lossless for this block.**

**Real images:** blocks with smooth structure have mean/std that vary slowly; BTC
compresses well but is no match for transform coding.

### 2.5.6 Wavelet coding

The **wavelet transform** (ECEG-6303 Ch 5 §5.9) decomposes an image into
multi-resolution coefficients. Because natural images are **piecewise smooth**,
most wavelet coefficients are **near zero** — high sparsity. Wavelet coding
exploits this:

1. **2-D DWT** (typically 5 levels with the CDF 9/7 biorthogonal wavelet).
2. **Quantize** coefficients — high-frequency subbands coarser (they capture less energy in natural images).
3. **Entropy-code** with context modelling (neighbouring coefficients are correlated — code them together).

**JPEG 2000** is the standard wavelet codec. Advantages over JPEG:

- **No blocking artifacts** — the DWT is global, not block-based.
- **Progressive transmission** — better rate-distortion trade-off (truncate the bitstream at any point).
- **Lossless and lossy** in one framework.
- **Better compression** at high PSNR (2–4 dB better than JPEG).
- **Region of interest** coding.
- **Higher computational cost** — this is why JPEG dominates the consumer market despite JPEG 2000's technical superiority.

**HEVC (H.265) and AV1** use **DCT-based** transform coding, not wavelets — the
transform tools have evolved but wavelets were not adopted due to hardware cost.

### 2.5.7 Colour image coding

Colour adds three times the raw data, but **the eye is less sensitive to chroma**.
Compression exploits this in two ways:

1. **Colour-space transform:** RGB → YCbCr (Chapter 1 §1.7).
2. **Chroma subsampling:** 4:4:4 → 4:2:0 → half the chroma bytes for ~0.1 dB PSNR loss.

Combined, YCbCr + 4:2:0 halves the raw image size (before any DCT, prediction,
or entropy coding). This is why colour compression achieves nearly 3× the ratio
of grayscale at the same perceptual quality.

### 2.5.8 JPEG vs. JPEG 2000 vs. HEIC vs. AVIF

| Codec | Transform | Block | Lossless | Progressiveness | Compression at Q=30 |
|---|---|---|:---:|:---:|---|
| JPEG | 8×8 DCT | Yes | No (except JPEG-LS) | Yes | 25:1 |
| JPEG 2000 | DWT (CDF 9/7) | No | Yes | Yes | 30:1 |
| WebP | VP8 (DCT+intra) | Yes | Yes (lossless mode) | Yes | 30:1 |
| HEIC (HEVC) | DCT | Yes | Yes | Yes | 35:1 |
| AVIF (AV1) | DCT + intra | Yes | Yes | Yes | 35:1 |

JPEG 2000 is technically superior in most respects, but JPEG's ubiquity and
its adequacy for most web use kept it dominant until newer codecs (WebP,
HEIC, AVIF) recently displaced it in browsers and mobile.

---

## 2.6 Reconstruction from projections

### 2.6.1 The Radon transform

For a 2-D function $f(x,y)$, its **Radon transform** is the set of **line integrals**
at each projection angle $\theta$ and offset $\rho$:

$$
\boxed{\;R(\rho, \theta) = \int_{-\infty}^{\infty}\int_{-\infty}^{\infty} f(x,y)\, \delta(x\cos\theta + y\sin\theta - \rho)\, dx\, dy.\;}
$$

**Geometric interpretation.** The line $x\cos\theta + y\sin\theta = \rho$ is a
straight line at distance $\rho$ from the origin, oriented at angle $\theta + 90°$.
The Radon transform integrates $f$ along that line.

**Physical meaning.** If $f$ is the attenuation coefficient of a body (as in CT),
$R(\rho, \theta)$ is the X-ray intensity loss along a ray at distance $\rho$
from the origin, oriented at angle $\theta$ — the **projection** of the body onto
a line perpendicular to the ray.

**Sinogram.** The 2-D plot of $R(\rho, \theta)$ vs. $(\rho, \theta)$ is called a
**sinogram**. A point source produces a sinusoid in the sinogram; a line
produces a horizontal line; an extended object produces an elliptical smear.

**Applications.**

- **CT (computed tomography):** reconstruct 3-D density from X-ray projections.
- **PET, SPECT:** reconstruct tracer distribution from gamma emissions.
- **Radar, sonar:** reconstruct object from scattered returns.
- **Electron microscopy:** 3-D structure of macromolecules.

### 2.6.2 Back-projection (unfiltered)

The naive inverse: smearing each projection back along its ray.

$$
\boxed{\;f_{\text{bp}}(x,y) = \int_0^\pi R(x\cos\theta + y\sin\theta,\ \theta)\, d\theta.\;}
$$

**The problem.** Back-projection is **not** the inverse Radon transform — it is
the original $f$ **convolved with $1/r$**:

$$
f_{\text{bp}}(x,y) \propto f(x,y) * \frac{1}{\sqrt{x^2 + y^2}}.
$$

So the reconstruction is **blurred** by a $1/r$ PSF — low frequencies correct,
high frequencies attenuated. Back-projection alone is never used in practice.

### 2.6.3 Filtered back-projection (FBP)

The **filtered back-projection** filter corrects this. It is based on the
**Fourier slice theorem**: the Fourier transform of the projection $R(\rho, \theta)$
at frequency $\omega$ equals the 2-D Fourier transform of $f$ along the line
through the origin at angle $\theta$:

$$
\boxed{\;\widehat{R}(\omega, \theta) = F(\omega\cos\theta, \omega\sin\theta).\;}
$$

This lets us reconstruct $f$ from all angles by inverse 2-D FT:

1. **Fourier transform** each projection: $\widehat{R}(\omega, \theta)$.
2. **Filter** in the Fourier domain with the **ramp filter** $|\omega|$ (or a smoothed version, the **Shepp–Logan** filter $|\omega| \cdot |\text{window}|$).
3. **Inverse FT** to get the filtered projection $\tilde R(\rho, \theta)$.
4. **Back-project** $\tilde R$ over all $\theta$:
   $$
   \hat f(x,y) = \int_0^\pi \tilde R(x\cos\theta + y\sin\theta,\ \theta)\, d\theta.
   $$

**The ramp filter** $|\omega|$ compensates for the $1/r$ blur of unfiltered
back-projection. In practice, a **Hamming** or **Shepp–Logan** window is applied
to suppress the noise amplification at high frequencies (the ramp is
unbounded).

**Reconstruction quality.**

- **Exact** (in the continuous, noiseless case) for full $180°$ coverage.
- **Excellent** in practice with 500–1000 projections and a good window.
- **Standard in commercial CT** — every clinical CT scanner uses FBP or a close variant (iterative methods are increasingly replacing it, but FBP is still the baseline).

### Worked example 2.17 — Radon transform of a point object

Point object at $(x_0, y_0)$: $f(x,y) = \delta(x - x_0, y - y_0)$.

$$
R(\rho, \theta) = \int\int \delta(x-x_0, y-y_0)\, \delta(x\cos\theta + y\sin\theta - \rho)\, dx\, dy
= \delta(x_0\cos\theta + y_0\sin\theta - \rho).
$$

For $\rho = x_0\cos\theta + y_0\sin\theta$, the projection is nonzero. The locus
of $\rho$ over $\theta$ is a **sinusoid** of amplitude $\sqrt{x_0^2 + y_0^2}$
and phase $\arctan(y_0/x_0)$ — the "sinogram of a point" is literally a sine
wave, which is why it is called a sinogram.

**Filtered back-projection** of this sinogram gives a sharp point (up to
finite-angle blur). **Unfiltered back-projection** gives the $1/r$ blur.

### 2.6.4 Fan-beam geometry

**Parallel-beam** is the idealised geometry: all rays at angle $\theta$ are
parallel. **Fan-beam** (more common in CT scanners) uses a **point source** and
a detector arc. The rays diverge from the source — each angle $\beta$ of the
source and each detector element $\gamma$ define a ray.

**Fan-beam → parallel-beam.** A fan-beam projection can be re-sorted to
parallel-beam by appropriate coordinate transformation (fan-to-parallel rebinning),
then the parallel FBP applies. Alternatively, **fan-beam FBP** has its own
reconstruction formula with a modified filter.

**Cone-beam.** The 3-D extension (each source position produces a cone of rays)
is the geometry of **CBCT** (cone-beam CT) and modern helical scanners.
Reconstruction uses **FDK** (Feldkamp–Davis–Kress) — the 3-D FBP approximation.

### 2.6.5 Algebraic reconstruction techniques (ART)

An alternative to FBP when the data is sparse (few angles) or noisy. **ART
formulates the reconstruction as a linear system** and solves it iteratively.

**Discretize.** The image is $N \times N$ pixels, unknown. Each projection ray
gives one equation:

$$
\boxed{\;\sum_{i,j} a_{ki,j}\, f_{i,j} = p_k,\quad a_{ki,j} = \text{length of ray } k \text{ through pixel } (i,j).\;}
$$

If there are $M$ rays, this is a linear system $A \mathbf f = \mathbf p$
with $M$ equations and $N^2$ unknowns.

**Iterative solution.**

- **Kaczmarz (ART):** cycle through the equations, projecting the current estimate onto each hyperplane:
  $$
  \mathbf f^{(m+1)} = \mathbf f^{(m)} + \lambda \frac{p_k - \mathbf a_k^\top \mathbf f^{(m)}}{\|\mathbf a_k\|^2} \mathbf a_k.
  $$
- **SART (simultaneous ART):** update using a whole projection at once — faster convergence, smoother.
- **SIRT (simultaneous iterative reconstruction):** weighted average of all projections.

**Pros.** Works with few angles; handles noise via regularization.

**Cons.** Slow (many iterations); needs careful preconditioning.

**Modern use.** Iterative reconstruction (IR) is standard on premium CT
scanners — it reduces X-ray dose by 30–50% compared to FBP while maintaining
image quality. **Compressed sensing** reconstruction (which exploits image
sparsity in a transform domain) goes further, using as few as 20–30
projections — used in some research and emerging clinical systems.

---

## 2.7 Chapter summary

1. **Transforms** are the foundation — orthogonal/unitary transforms preserve energy and invert cheaply.
2. **DCT** is the practical transform of choice for compression: real, fast, excellent energy compaction, near-KLT-optimal for natural images.
3. **KLT** is theoretically optimal but data-dependent (not usable in practice). **SVD** is similarly optimal but data-dependent (used for denoising and low-rank problems).
4. **Walsh–Hadamard, Haar, Slant** are historical — mainly of interest for their simple structure (Haar feeds into wavelets).
5. **Enhancement** is subjective — point operations (contrast stretch, gamma), histograms (equalization, specification), spatial filtering (blur, sharpen), median filtering (impulse noise), unsharp masking (edges).
6. **Restoration** is objective — inverse (fails with noise), **Wiener** (optimal MSE with known PSDs), **CLS** (robust with smoothness prior).
7. **DPCM** is simple predictive coding; **DCT coding** (JPEG) is the practical standard; **wavelet coding** (JPEG 2000) is superior technically but less widespread; **BTC** is fast but low quality.
8. **Radon transform** gives projections; **back-projection** gives blurred reconstruction; **filtered back-projection** (with the ramp filter) gives sharp reconstruction; **ART** handles few-angle / noisy cases.

---

## 2.8 Quiz

### A. Multiple choice

1. The **DCT** is preferred over the **DFT** for image compression because:
   (a) it is faster (b) its coefficients are real and it has better energy compaction (c) it is more accurate (d) it requires no quantization

2. The **KL transform** is theoretically optimal but not used in practice because:
   (a) it is too slow (b) it is data-dependent (c) it produces complex coefficients (d) it is not orthogonal

3. The **SVD** of an image is used primarily for:
   (a) compression standards (b) denoising and low-rank approximation (c) JPEG encoding (d) histogram equalization

4. **Histogram equalization** transforms the image so its histogram is:
   (a) Gaussian (b) approximately uniform (c) exponential (d) identical to the original

5. **Median filtering** is preferred over mean filtering for:
   (a) Gaussian noise (b) impulse (salt-and-pepper) noise (c) blur (d) sharpening

6. **Unsharp masking** enhances edges by:
   (a) averaging (b) subtracting a lowpass version from the image and adding it back (c) histogram equalization (d) thresholding

7. **Inverse filtering** fails in practice because:
   (a) it is too slow (b) it amplifies noise at high frequencies (c) it needs a prior (d) it is not linear

8. The **Wiener filter** minimizes:
   (a) peak error (b) mean-square error given signal and noise PSDs (c) mutual information (d) total variation

9. **JPEG** compression uses which transform?
   (a) DWT (b) 8×8 DCT (c) KLT (d) SVD

10. **Block truncation coding** represents each block by:
    (a) DCT coefficients (b) mean, standard deviation, and 1-bit mask (c) wavelet coefficients (d) KLT basis

11. The **Radon transform** of a point object produces:
    (a) a point (b) a sinusoid in the sinogram (c) a disk (d) zero

12. **Filtered back-projection** uses which filter?
    (a) Gaussian (b) ramp $|\omega|$ (c) box (d) none

### B. Fill in the blank

1. The DCT of a constant $8\times 8$ block has only one nonzero coefficient at (0,0) — the \_\_\_\_ coefficient.
2. The KLT basis is given by the \_\_\_\_ of the covariance matrix.
3. SVD of an image is $A = U\Sigma V^\top$ with \_\_\_\_ on the diagonal.
4. Histogram equalization maps intensities by their \_\_\_\_ function.
5. The Laplacian kernel is $\begin{bmatrix}0&1&0\\1&$\_\_\_\_$&1\\0&1&0\end{bmatrix}$.
6. **Wiener filter:** $W = H^*/(|H|^2 + S_\nu/$\_\_\_\_).
7. **Inverse filter** amplifies noise as $N/$\_\_\_\_ at each frequency.
8. **DPCM** codes the prediction \_\_\_\_.
9. **BTC** stores \_\_\_\_, \_\_\_\_, and a mask per block.
10. **Filtered back-projection** filter is $|\omega|$, called the \_\_\_\_ filter.

### C. Define / explain

1. State the **DFT**, **DCT**, and **KLT** in one sentence each. Compare their energy compaction for a natural image (which is best, which is worst).
2. Derive the **2-D DCT** from the 1-D DCT. Why is the DCT asymptotically optimal for AR(1) processes with $\rho \to 1$?
3. Describe the **8×8 DCT basis functions**. What does the low-frequency corner represent? The high-frequency corner?
4. Write the **JPEG compression pipeline**. What are the three types of redundancy JPEG exploits?
5. Derive the **histogram equalization mapping** $s = (L-1)C(r)$. Show that it is monotone and produces a uniform histogram in the continuous limit.
6. Explain the difference between **contrast stretching**, **histogram equalization**, and **histogram specification**.
7. Compare **median**, **mean**, **bilateral**, and **non-local means** filters. What kind of noise / structure is each suited for?
8. Derive the **inverse filter** formula. Explain with a numerical example why it fails on real images with noise.
9. Derive the **Wiener filter** formula from the minimum-MSE criterion. How does it interpolate between inverse filtering and no filtering?
10. Derive **BTC**. What is the compression ratio for a 4×4 block? What is the reconstruction quality trade-off?
11. State the **Radon transform** and the **Fourier slice theorem**. Use the slice theorem to derive the **filtered back-projection** algorithm.

---

## 2.9 Answer key

??? success "A. Multiple choice"
    1. **(b)** — real coefficients + better energy compaction.
    2. **(b)** — data-dependent.
    3. **(b)** — denoising / low-rank.
    4. **(b)** — approximately uniform.
    5. **(b)** — impulse (salt-and-pepper).
    6. **(b)** — subtract lowpass, add back.
    7. **(b)** — amplifies noise at high frequencies.
    8. **(b)** — MMSE given PSDs.
    9. **(b)** — 8×8 DCT.
    10. **(b)** — mean, std, mask.
    11. **(b)** — sinusoid in the sinogram.
    12. **(b)** — ramp $|\omega|$.

??? success "B. Fill in the blank"
    1. DC (zero-frequency)
    2. eigenvectors
    3. singular values
    4. cumulative distribution (CDF)
    5. $-4$
    6. $S_f$ (signal PSD)
    7. $H$ (blur transfer function)
    8. error (residual)
    9. mean, standard deviation
    10. ramp

??? success "C. Define / explain — model answers"
    1. **DFT:** $X[k] = \sum x[n] e^{-j2\pi kn/N}$ — Fourier basis (complex exponentials). **DCT:** $X[k] = c_k \sum x[n] \cos((2n+1)k\pi/(2N))$ — real cosines with even symmetry. **KLT:** $y = U^\top(x - \mu)$ where $U$ = eigenvectors of $x$'s covariance. **Energy compaction:** KLT optimal > DCT (very close, 92–95% for natural-image correlation) >> DFT (80–85%, boundary leakage). For practical codecs, DCT wins because KLT requires transmitting the basis.
    2. **1-D to 2-D DCT:** apply 1-D DCT along rows, then along columns. Since the DCT basis is separable, the result is a tensor product of 1-D cosine bases. **Asymptotic optimality:** for an AR(1) process with correlation $r[k] = \rho^{|k|}$, as $\rho \to 1$ the KL eigenvectors converge to the DCT basis. Practical natural images have $\rho \approx 0.9$–0.95, so DCT is nearly optimal, explaining its dominance in JPEG, MPEG, H.264, and HEVC.
    3. **8×8 DCT basis.** Each of the 64 basis functions is a 2-D separable cosine with frequency $(k, l)$ in the two directions. The top-left (0,0) is the DC term — a uniform basis (constant over the block). Top-right (0,7) is high-frequency vertically, low-frequency horizontally. Bottom-right (7,7) is the highest-frequency 2-D checkerboard — flips sign on every pixel. Low frequencies dominate typical natural images; high frequencies capture edges and textures. **Compression intuition:** natural images have most energy in low frequencies (smooth), so most DCT energy sits in the top-left corner — the rest is small and can be coarsely quantized.
    4. **JPEG pipeline:** (i) split into 8×8 blocks; (ii) RGB → YCbCr, chroma 4:2:0; (iii) 8×8 DCT per block; (iv) quantize using a table scaled by a quality factor; (v) zigzag scan; (vi) Huffman (or arithmetic) entropy coding of the runs. **Three redundancies:** (i) **coding** — Huffman exploits uneven symbol probabilities; (ii) **interpixel (spatial)** — DCT decorrelates neighbouring pixels; (iii) **psychovisual** — the quantization table leverages the eye's lower sensitivity at high frequencies, allocating bits where they are perceptually most needed.
    5. **Histogram equalization derivation.** In the continuous limit, intensity $r$ has PDF $p_r(r)$ and CDF $C(r) = \int_0^r p_r(s) ds$. Define the output $s = (L-1) C(r)$. Then $s$ has PDF $p_s(s) = p_r(r)/((L-1) p_r(r)) = 1/(L-1)$ for $s \in [0, L-1]$ — **uniform**. Monotone because $C(r)$ is monotone (integral of a non-negative density). Discrete case uses the empirical CDF.
    6. **Contrast stretching:** linear scaling that maps the observed minimum and maximum to 0 and $L-1$. Simple but ignores the histogram shape. **Histogram equalization:** nonlinear mapping based on the cumulative distribution — maximises entropy, spreads the histogram as uniformly as possible. **Histogram specification:** nonlinear mapping to a **specific target PDF** — more flexible than equalization (target need not be uniform).
    7. **Mean filter:** linear, blurs uniformly. Best for Gaussian noise, destroys edges. **Median:** nonlinear, excellent for salt-and-pepper (impulse) noise, preserves edges but destroys fine textures and thin lines. **Bilateral:** weighted average where the weight depends on both spatial distance and intensity similarity — smooths noise while preserving edges. Best for Gaussian noise in textured images. **Non-local means:** averages patches that are similar anywhere in the image. Excellent for Gaussian noise in natural images. Slower.
    8. **Inverse filter:** in Fourier domain, $\hat F = G/H = F + N/H$. **Failure:** wherever $H(f) \to 0$ (high frequency for most blur models), $1/H \to \infty$, so the noise term $N/H$ amplifies without bound. Even a tiny noise level (e.g., 0.1% of the image RMS) becomes amplified to 10–100% of the image RMS at those frequencies. The output is dominated by noise. **Example (Worked Example 2.11):** Gaussian blur $H(f) = e^{-(f/f_0)^2}$ with noise std $10^{-3}$: at $f = 3f_0$, $H = 10^{-4}$, so noise is amplified by $10^4$ — from $10^{-3}$ to $10$, i.e., 10× the signal. Unusable.
    9. **Wiener derivation.** Minimize $\mathbb E[\|F - \hat F\|^2]$ over linear filters $\hat F = W \cdot G$. Substitute $G = HF + N$ and the signal/noise independence. The optimal $W$ is $\mathbb E[F G^*]/\mathbb E[|G|^2] = H^* S_f/(|H|^2 S_f + S_\nu)$. **Interpolation:** where $S_f \gg S_\nu$ (high SNR), $W \approx 1/H$ — inverse filtering. Where $S_f \ll S_\nu$ (low SNR), $W \approx 0$ — no restoration, no noise amplification. In between, the filter smoothly transitions based on the local SNR.
    10. **BTC.** For each 4×4 block: compute mean $\bar x$ and std $\sigma_x$; binarize at $\bar x$ (1 if above, 0 if below); set reconstruction levels $a = \bar x - \sigma_x\sqrt{q/m}$, $b = \bar x + \sigma_x\sqrt{m/q}$; store mask (16 bits) + quantized $\bar x$ (8 bits) + quantized $\sigma_x$ (8 bits) = **32 bits per block**. **Compression:** original 4×4 8-bit = 128 bits; BTC = 32 bits; **4:1**. **Quality trade-off:** the two-level approximation loses mid-range detail — textures and smooth gradients are degraded, but edges are preserved. Reconstruction is blurry in smooth areas, crisp on edges. Acceptable for low-rate applications, poor at high rate.
    11. **Radon transform:** $R(\rho, \theta) = \int\int f(x,y) \delta(x\cos\theta + y\sin\theta - \rho) dx dy$ — line integrals at angle $\theta$. **Fourier slice theorem:** $\widehat R(\omega, \theta) = F(\omega\cos\theta, \omega\sin\theta)$ — the 1-D FT of the projection equals the 2-D FT of $f$ along the slice through the origin at angle $\theta$. **FBP derivation:** the inverse 2-D FT can be written as a polar integral: $f(x,y) = \int_0^\pi \left[\int |\omega| \widehat R(\omega, \theta) e^{j2\pi\omega(x\cos\theta + y\sin\theta)} d\omega\right] d\theta$. The inner bracket is the **filtered projection** $\tilde R(\rho, \theta)$ — the projection filtered by the ramp $|\omega|$. The outer integral is the **back-projection** — smearing the filtered projection back along its ray. So FBP = filter with ramp + back-project. The ramp compensates for the $1/r$ blur of plain back-projection.

---

## Further reading

- Jain, *Fundamentals of Digital Image Processing*, PHI, 1995 — Chapters 5–7 cover enhancement, restoration, and compression.
- Gonzalez & Woods, *Digital Image Processing*, 4th ed., Pearson, 2018 — Chapters 3–6 and 10 for the full classical treatment.
- Pratt, *Digital Image Processing*, 4th ed., Wiley, 2007 — the encyclopedic reference.
- Rao & Yip, *The Transform and Data Compression Handbook*, CRC Press, 2001 — DCT, KLT, and other transforms in depth.
- Ahmed, Natarajan & Rao, "Discrete Cosine Transform", *IEEE Trans. Computers*, 1974 — the original DCT paper.
- Oppenheim & Schafer, *Discrete-Time Signal Processing*, PHI — for the transform theory.
- Lim, *Two-Dimensional Signal and Image Processing*, Prentice-Hall, 1990 — comprehensive 2-D transform and filter theory.
- Jain, "Image data compression: A review", *Proc. IEEE*, 1981 — the classical survey.
- Clarke, *Transform Coding of Images*, Academic Press, 1985 — detailed DCT analysis.
- Wallace, "The JPEG Still Picture Compression Standard", *Commun. ACM*, 1991 — the definitive JPEG reference.
- Taubman & Marcellin, *JPEG 2000: Image Compression Fundamentals, Standards and Practice*, Springer, 2002 — the definitive JPEG 2000 book.
- **Tour of Filtered Back-Projection:** Kak & Slaney, *Principles of Computerized Tomographic Imaging*, IEEE Press, 1988 (available free online).
- Natterer, *The Mathematics of Computerized Tomography*, Wiley, 2001 — rigorous FBP theory.
- Buades, Coll & Morel, "A Non-Local Algorithm for Image Denoising", *CVPR* 2005 — the non-local means paper.
- Dabov et al., "Image Denoising by Sparse 3D Transform-Domain Collaborative Filtering", *IEEE TIP*, 2007 — the BM3D paper.
- SciPy, scikit-image, OpenCV, PyWavelets documentation.

---

## Suggested tutorial and lab problems for this chapter

1. **DCT energy compaction.** Apply the 2-D DCT to a natural image (or Lena). Sort the coefficients by magnitude. What fraction of the total energy is in the top 1%, 5%, 10% of coefficients? Repeat for the DFT — quantify the difference.

2. **JPEG simulation.** Implement the JPEG pipeline (DCT, quantization, entropy coding) for a single grayscale image. Compare compression ratios and PSNR at quality levels 10, 50, 90. Visualise the blocking artifacts at low quality.

3. **KL vs. DCT.** Compute the KLT of a small correlated image block, truncate to $M$ components, and reconstruct. Do the same with the DCT. At what $M$ does each reach 30 dB PSNR?

4. **Histogram equalization.** Implement histogram equalization and histogram specification from scratch. Apply to a dark image and to an overexposed image. Compare with `cv2.equalizeHist`.

5. **Median vs. bilateral denoising.** Add salt-and-pepper and Gaussian noise to a natural image. Compare median, bilateral, and non-local means for both noise types. Which is best for which?

6. **Wiener filter with unknown noise.** Blur a clean image with a Gaussian PSF (σ = 2). Add white Gaussian noise. Apply inverse filtering, Wiener filtering, and CLS. Compare PSNR. Tune the regularisation parameter.

7. **Predictive coding.** Implement DPCM with 1-D, 2-D, and 3-tap predictors on a grayscale image. Compute the entropy of the residuals for each. Plot the rate-distortion curve.

8. **BTC implementation.** Implement BTC with block sizes 4×4 and 8×8, with quantized mean/std to 4/6/8 bits. Compare compression ratio and PSNR to JPEG at the same rate.

9. **Wavelet denoising.** Take a noisy image, apply a 5-level DWT (db4), hard- or soft-threshold the detail coefficients at $3\sigma$, and reconstruct. Compare PSNR to the Wiener filter.

10. **Radon transform.** Implement the Radon transform of a small binary shape (e.g., a disc, a rectangle). Visualise the sinogram. Reconstruct via unfiltered back-projection and filtered back-projection. Compare the two reconstructions — the blur of unfiltered BP should be visible.

11. **FBP on a Shepp–Logan phantom.** Reconstruct the Shepp–Logan phantom (standard test image for CT) from projections at 90°, 180°, and 360°. Compare quality at each. Demonstrate the failure of FBP with too few angles.

12. **Art with noisy projections.** Add 1% Gaussian noise to the projections of a Shepp–Logan phantom. Compare FBP reconstruction vs. iterative reconstruction (ART or SIRT, 50 iterations). Report PSNR and visual quality.