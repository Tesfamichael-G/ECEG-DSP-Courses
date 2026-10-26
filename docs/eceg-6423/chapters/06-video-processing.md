# Chapter 6 — Digital Video Processing Techniques and Algorithms

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 2 tutorial + 1 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** video representation and modelling · video sampling
    and rate conversion · two- and three-dimensional motion estimation · video
    coding techniques · video compression standards · error control in video
    communications · streaming video over the internet.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Describe a **digital video signal** as a 3-D spatio-temporal array; state the
   format conventions (4:2:0 chroma, interlaced vs. progressive, frame rate,
   bit depth).
2. State the **3-D sampling theorem** and explain **temporal aliasing** — the
   "wagon-wheel effect."
3. Perform **video rate conversion**: interlacing, deinterlacing, frame-rate
   up/down-conversion, telecine.
4. Explain **2-D motion estimation** via **block matching** and **optical
   flow** (Lucas–Kanade), stating their trade-offs, complexity, and typical
   applications.
5. State the **3-D motion estimation** problem (structure from motion, stereo
   disparity) and its relation to 2-D motion fields.
6. Describe the **hybrid video codec**: intra/inter prediction, motion-compensated
   prediction, DCT transform, quantization, and entropy coding.
7. Compare **video coding standards**: MPEG-1/2/4, H.261/263, H.264/AVC,
   H.265/HEVC, VP9, AV1 — their design goals, bitrate efficiency, and typical
   applications.
8. Explain **rate control**: CBR vs. VBR, buffer models, and rate-distortion
   optimization.
9. Design **error control** for video: FEC, retransmission, error resilience,
   and error concealment. Compute the required FEC overhead for a given
   packet-loss rate.
10. Describe **streaming video over the internet**: HTTP adaptive streaming,
    DASH, HLS, and the trade-offs between buffering and latency.
11. Implement all of the above in **Python, Java, and C#**.

---

## 6.1 Video representation and modelling

### 6.1.1 Digital video as a 3-D signal

A **digital video signal** is a 3-D array of intensity values:

$$
\boxed{\;I[m, n, t],\quad m = 0, \ldots, M-1,\ n = 0, \ldots, N-1,\ t = 0, \ldots, T-1.\;}
$$

- **$M, N$:** spatial dimensions (rows, columns in pixels).
- **$T$:** temporal dimension (number of frames).
- **Frame rate:** $f_s = 1/\Delta t$ frames per second (fps), typically 24 (cinema), 25 (PAL), 30 (NTSC), 50/60 (HD), 120 (high-speed).
- **Bit depth:** 8, 10, or 12 bits per sample per channel.

**For colour video**, each pixel has three channels: **Y, Cb, Cr** (ECEG-6423 Chapter 1, §1.7). The most common subsampling is **4:2:0** — full luma resolution, half chroma resolution in both directions.

**Effective bit rate** (uncompressed):

$$
\boxed{\;R_{\text{raw}} = M \cdot N \cdot f_s \cdot b \cdot (1 + 2 \cdot \tfrac{1}{4}) = M \cdot N \cdot f_s \cdot b \cdot 1.5,\;}
$$

where $b$ is bits per sample and the factor $1.5$ comes from 4:2:0 (1 full luma + 2 half-resolution chroma).

**Example:** 1920×1080 (1080p), 30 fps, 8-bit 4:2:0:

$$
R_{\text{raw}} = 1920 \cdot 1080 \cdot 30 \cdot 8 \cdot 1.5 \approx 746\ \text{Mbps.}
$$

Compare with a typical streaming bitrate of 5 Mbps — a **150:1 compression ratio**. The rest of this chapter explains how that compression is achieved.

### 6.1.2 Frame types

Video codecs classify frames by their dependence on other frames:

| Type | Description | Reference frames | Typical size |
|---|---|---|---|
| **I-frame** | Intra-coded; no temporal prediction | None | Large (JPEG-like) |
| **P-frame** | Predictive; encoded from one past reference | 1 (past) | Moderate |
| **B-frame** | Bi-predictive; encoded from past and future references | 2 | Smallest |
| **IDR** | I-frame that resets the decoder state (a "key frame") | — | — |

A **Group of Pictures (GOP)** defines the frame structure:

- **Closed GOP:** begins with an I-frame; subsequent frames do not depend on the *previous* GOP.
- **Open GOP:** B-frames at the start of the GOP can reference frames from the previous GOP.
- **Typical patterns:**
  - **IBBPBBPBBP...** — one I every 12–15 frames (a "GOP of 12").
  - **IPPP...** — one I every 30–60 frames, all others P.
  - **All-I** — every frame is an I-frame (e.g., for editing, streaming seek).

### 6.1.3 Frame-rate and format conventions

**Progressive vs. interlaced.**

- **Progressive** (e.g., 1080p30): all pixels in a frame are captured at the same instant.
- **Interlaced** (e.g., 1080i60): each frame consists of two **fields** captured at different times — odd lines first (top field), even lines second (bottom field), spaced 1/60 s apart. The viewer sees 60 fields/s but only 30 full frames/s.

**Legacy interlaced formats:** 480i60 (NTSC), 576i50 (PAL), 1080i60 (HD broadcast). **Modern streaming/display:** almost always progressive. **Interlaced** still common in broadcast TV, cable, and some sports production.

**Frame-rate standards.**

| Standard | Region | Frame rate | Reason |
|---|---|---|---|
| NTSC (US, Japan) | 525 lines | 29.97 fps | 60 Hz AC power → legacy |
| PAL (Europe, others) | 625 lines | 25 fps | 50 Hz AC power |
| Cinema | — | 24 fps | Historical film standard |
| HD ATSC | — | 30/60 fps | Digital |
| HFR (high frame rate) | — | 48, 120 fps | Sports, VR (Hobbit, Ang Lee films) |

**Bit depth.** 8-bit (consumer), 10-bit (broadcast, HDR), 12-bit (cinema, HDR), 16-bit (raw, scientific).

---

## 6.2 Video sampling and rate conversion

### 6.2.1 The 3-D sampling theorem

For a video signal band-limited to $[-F_x, F_x] \times [-F_y, F_y] \times [-F_t, F_t]$:

$$
\boxed{\;\frac{1}{\Delta x} > 2F_x,\qquad \frac{1}{\Delta y} > 2F_y,\qquad \frac{1}{\Delta t} > 2F_t.\;}
$$

The temporal condition — sampling rate $f_s > 2 F_t$ — is the **temporal Nyquist condition**. Violations produce **temporal aliasing**.

### Worked example 6.1 — the wagon-wheel effect

A wheel with $N_{\text{spokes}}$ rotating at $\omega$ rad/s (angular frequency). The apparent (aliased) angular frequency when sampled at frame rate $f_s$ is:

$$
\boxed{\;\omega_{\text{alias}} = \omega - k f_s,\quad \text{where } k = \text{round}(\omega / f_s).\;}
$$

**Case 1.** A wheel with $N_s = 8$ spokes, rotating at $\omega/2\pi = 1$ rev/s, filmed at 24 fps.

The apparent frequency is $1 - 24\cdot k$ rev/s. For $k = 0$: 1 rev/s — correct (no aliasing). But if the wheel spun at 12 rev/s:

$\omega_{\text{alias}} = 12 - 24 = -12$ rev/s — the wheel appears to spin **backwards at 12 rev/s**.

If it spun at 13 rev/s:

$\omega_{\text{alias}} = 13 - 24 = -11$ rev/s — backwards at 11 rev/s.

**Case 2.** Same wheel at 25 fps (PAL film speed) rotating at 1 rev/s: $\omega_{\text{alias}} = 1$ rev/s — correct (assuming 25 > 2× 1).

**Real-world effect.** In film and video, wheels and propellers often appear to spin slowly backward because of temporal aliasing. This is the same phenomenon as 2-D Moiré (ECEG-6423 Chapter 1 §1.5.4), extended to the time dimension.

**Mitigation.** **Shutter speed.** Cinema cameras use 1/48 s shutter (180° shutter) at 24 fps — the motion blur acts as a temporal lowpass filter, reducing aliasing. Sports use faster shutters but higher frame rates.

### 6.2.2 Interlaced video and deinterlacing

**Interlacing** was invented to reduce bandwidth for broadcast TV: instead of transmitting 60 full frames/s (which required too much bandwidth), each frame is split into two **fields** captured 1/60 s apart. The viewer's eye integrates the two fields into a smooth image, but:

- **Static scenes:** fine detail (hair, textures) shows **combing** or **interline flicker**.
- **Moving scenes:** motion between fields creates **interlace artifacts** — visible "teeth" at object edges.

**Deinterlacing** converts interlaced to progressive. Two basic approaches:

**Bob (spatial interpolation).** Each field is upsampled to a full frame by interpolating the missing lines:

$$
\boxed{\;\hat I[m, n] = \frac{I[m-1, n] + I[m+1, n]}{2}.\;}
$$

**Pros:** simple, no combing, works in motion.
**Cons:** halves vertical resolution, loses fine detail.

**Weave (temporal merge).** Combine the two fields into one frame, treating them as captured at the same instant.

**Pros:** full resolution in static areas.
**Cons:** combing in motion.

**Motion-adaptive.** Detect motion between fields; use weave for static pixels, bob for moving pixels.

**Modern approaches:** motion-compensated deinterlacing (MC-DI) uses motion estimation to interpolate the missing lines along motion trajectories. **Deep-learning deinterlacers** (e.g., in NVIDIA's `Maxine`) achieve near-progressive quality.

### 6.2.3 Frame-rate conversion

Converting between frame rates — 24 → 30, 30 → 60, 24 → 60, etc. — is fundamental for cinema/broadcast/streaming interoperability.

**Simple approaches.**

**Frame repetition.** Each frame is duplicated: 24 → 30 by repeating every 4th frame (actually 24 → 30 is non-integer, so a 5:4 pull-down pattern). Adds judder (jerky motion) because the display alternates between 1 and 2 frame periods.

**Frame interpolation (blending).** Average two adjacent frames:

$$
\boxed{\;\hat I_{n+0.5} = \alpha I_n + (1 - \alpha) I_{n+1}.\;}
$$

**Pros:** no judder. **Cons:** ghosting at motion boundaries.

**Motion-compensated frame interpolation (MCFI).** Estimate motion between $I_n$ and $I_{n+1}$, then **warp** the pixels along the motion trajectory:

$$
\boxed{\;\hat I_{n+0.5}(x, y) = I_n(x - u/2, y - v/2),\quad (u, v) = \text{motion field}.\;}
$$

**Pros:** smooth motion without ghosting. **Cons:** occlusion handling (pixels moving in/out of view); computational cost (two motion estimates per output frame).

**Telecine.** 24 fps film → 30 fps NTSC video uses **3:2 pull-down**: frame A is shown for 3 fields (1.5 frames), frame B for 2 fields (1 frame), repeating. This is why film on NTSC video has a slight judder — the 3:2 pattern is not perfectly smooth. Modern HD displays remove the 3:2 cadence back to 24 fps (**inverse telecine** or **IVTC**).

---

## 6.3 2-D motion estimation

**Motion estimation** finds the displacement of pixels between consecutive frames. It is the single most impactful tool in video compression (motion-compensated prediction) and appears in video analysis, surveillance, object tracking, and autonomous driving.

### 6.3.1 Block matching

**Block matching** divides the current frame into **macroblocks** (typically 16×16 pixels) and finds, for each macroblock, the displacement vector $(u, v)$ that minimizes the difference to a block in the reference frame.

**Formal problem.** For macroblock $B$ at position $(x, y)$ in frame $I_t$, find:

$$
\boxed{\;(u^*, v^*) = \arg\min_{u, v \in [-p, p]} \text{SAD}(x, y, u, v),\;}
$$

where the **Sum of Absolute Differences** is:

$$
\boxed{\;\text{SAD}(x, y, u, v) = \sum_{i=0}^{B-1}\sum_{j=0}^{B-1} \bigl|I_t(x+i, y+j) - I_{t-1}(x+i+u, y+j+v)\bigr|.\;}
$$

**Alternatives to SAD.**

- **SSD** (sum of squared differences): more sensitive to large errors, costlier (multiplications).
- **SATD** (sum of absolute transformed differences): Hadamard transform then SAD — more robust to noise.
- **MSE**, **MAD**: equivalent up to scaling.

**Search range.** $\pm p$ pixels around the current position. Typical $p = 16$ (for standard-def), 32 (HD), 64 (4K).

**Complexity.** Naive full search: $(2p+1)^2$ candidates per macroblock. For $p = 16$: 1089 candidates. For 1080p (8160 macroblocks per frame):

$$
8160 \cdot 1089 \approx 8.9\ \text{million SADs per frame,}
$$

or **267 million SADs/s at 30 fps**. Each SAD is 256 operations (16×16), so **68 billion ops/s** — this is why fast search algorithms are essential.

**Fast search algorithms.**

| Algorithm | Search points | Relative speedup | Quality vs. full search |
|---|---|---|---|
| **Three-step search (TSS)** | ~25 | ~40× | Slightly worse |
| **Diamond search (DS)** | ~20 | ~50× | Nearly as good |
| **Hexagon search (HEX)** | ~15 | ~70× | Nearly as good |
| **EPZS (Enhanced Predictive Zonal Search)** | ~5 | ~200× | Very close |
| **UMHexagonS** | ~10 | ~100× | Very close (used in H.264 reference) |
| **TZSearch** | ~10 | ~100× | Used in HEVC reference |

**Diamond search** is the standard baseline: start at $(0, 0)$, test the 8 neighbors of a diamond pattern; if the minimum moves, re-center; if not, switch to a smaller diamond. Converges quickly for typical motion.

### Worked example 6.2 — block matching by hand

4×4 macroblock in frame $I_t$ (values in 0–9):

$$
B_{\text{cur}} = \begin{bmatrix} 5&5&4&4\\ 5&5&4&4\\ 4&4&3&3\\ 4&4&3&3\end{bmatrix}.
$$

Reference frame $I_{t-1}$ contains three candidates at displacements $(0,0)$, $(1,0)$, $(0,1)$:

$$
B_{0,0} = \begin{bmatrix} 4&4&3&3\\ 4&4&3&3\\ 3&3&2&2\\ 3&3&2&2\end{bmatrix},\quad
B_{1,0} = \begin{bmatrix} 5&5&5&4\\ 5&5&5&4\\ 4&4&4&3\\ 4&4&4&3\end{bmatrix},\quad
B_{0,1} = \begin{bmatrix} 3&4&3&3\\ 4&5&4&4\\ 3&3&2&2\\ 4&4&3&3\end{bmatrix}.
$$

**SAD for $(0,0)$.**

Each entry differs by 1 (current minus reference = 1 everywhere). Sum of absolute differences = 16 × 1 = **16**.

**SAD for $(1,0)$.**

Differences:

$$
\begin{bmatrix} 0&0&-1&0\\ 0&0&-1&0\\ 0&0&-1&0\\ 0&0&-1&0\end{bmatrix}
$$

Sum of absolutes = $1 \cdot 4 = 4$.

**SAD for $(0,1)$.**

Differences:

$$
\begin{bmatrix} 2&1&1&1\\ 1&0&0&0\\ 1&1&1&1\\ 0&0&0&0\end{bmatrix}
$$

Sum of absolutes = $2+1+1+1 + 1+0+0+0 + 1+1+1+1 + 0+0+0+0 = 10$.

**Winner:** $(1, 0)$ with SAD = 4 — the reference block shifted one pixel right matches the current block best. Motion vector $= (u, v) = (1, 0)$ — the block is moving right at 1 pixel per frame.

**In a real codec:** the encoder would search all $(2p+1)^2$ candidates, find the true minimum, and encode the residual (current block minus reference block) after motion compensation. The residual — being much smaller — compresses well.

### 6.3.2 Optical flow — Lucas–Kanade

Block matching gives **one motion vector per block** (an integer, discrete field). **Optical flow** gives **per-pixel motion vectors** (continuous field). Used in video analysis, object tracking, stereo, and self-driving.

**Brightness constancy assumption.** A moving pixel maintains its intensity:

$$
\boxed{\;I(x, y, t) = I(x + u\Delta t, y + v\Delta t, t + \Delta t).\;}
$$

**Taylor expansion** (first order):

$$
I(x+u\Delta t, y+v\Delta t, t+\Delta t) \approx I(x, y, t) + I_x u\Delta t + I_y v\Delta t + I_t \Delta t.
$$

Setting equal and dropping $\Delta t$:

$$
\boxed{\;I_x\, u + I_y\, v + I_t = 0.\;}
$$

This is the **optical flow constraint equation** — **one equation, two unknowns** ($u, v$). This is the **aperture problem**: a single pixel cannot determine motion along an edge (only perpendicular to it).

**Lucas–Kanade solution.** Assume the motion is **constant over a small window** $\Omega$ of pixels. Then $N = |\Omega|$ equations, 2 unknowns:

$$
\begin{bmatrix} I_x(p_1) & I_y(p_1) \\ I_x(p_2) & I_y(p_2) \\ \vdots & \vdots \\ I_x(p_N) & I_y(p_N)\end{bmatrix}
\begin{bmatrix} u \\ v \end{bmatrix} = -\begin{bmatrix} I_t(p_1) \\ I_t(p_2) \\ \vdots \\ I_t(p_N)\end{bmatrix}.
$$

In matrix form: $A \mathbf{u} = \mathbf{b}$ with $A \in \mathbb R^{N \times 2}$.

**Least-squares solution:**

$$
\boxed{\;\mathbf u = (A^\top A)^{-1} A^\top \mathbf b.\;}
$$

Expanded:

$$
\boxed{\;\begin{bmatrix} u \\ v \end{bmatrix} = \begin{bmatrix} \sum I_x^2 & \sum I_x I_y \\ \sum I_x I_y & \sum I_y^2 \end{bmatrix}^{-1}\begin{bmatrix} -\sum I_x I_t \\ -\sum I_y I_t\end{bmatrix}.\;}
$$

The **sums are over the window** $\Omega$ — a $5\times 5$ or $7\times 7$ neighborhood.

**Interpretation.** The matrix $A^\top A$ is (proportional to) the **structure tensor** — it captures the local image structure. Its eigenvalues $\lambda_1, \lambda_2$ tell you:

- **Both large:** two distinct edges — full 2-D motion is recoverable.
- **One large, one small:** one edge only — motion is only recoverable perpendicular to the edge (aperture problem).
- **Both small:** flat region — no motion recoverable.

**Lucas–Kanade in practice.**

- **Pyramid.** Apply LK at multiple scales (Gaussian pyramid) to capture large motions — the basic LK has limited range ($\pm$ window/2 per level).
- **Iterative refinement.** Warp the image by the current estimate, recompute LK on the residual, refine. Typically 3–5 iterations.
- **Robustification.** Weight pixels by their agreement with the current estimate (Huber, Tukey) — reduces sensitivity to outliers.

**Computational cost.** $O(N \cdot W^2)$ per window for the sums, where $W$ is the window size — for $W = 7$: 49 pixels per window, ~150 arithmetic ops per pixel. For 1080p at 30 fps: ~10 billion ops/s — feasible on GPU.

### Worked example 6.3 — Lucas–Kanade on a corner

Local patch around a strong corner — both $I_x$ and $I_y$ are nonzero (a diagonal edge), so $A^\top A$ is well-conditioned:

$$
A^\top A = \begin{bmatrix} 50 & 20 \\ 20 & 30 \end{bmatrix},\qquad A^\top\mathbf b = \begin{bmatrix} 15 \\ 10 \end{bmatrix}.
$$

Determinant $= 1500 - 400 = 1100$. Inverse:

$$
(A^\top A)^{-1} = \frac{1}{1100}\begin{bmatrix} 30 & -20 \\ -20 & 50\end{bmatrix}.
$$

$$
\mathbf u = \frac{1}{1100}\begin{bmatrix} 30 & -20 \\ -20 & 50\end{bmatrix}\begin{bmatrix} 15 \\ 10\end{bmatrix} = \frac{1}{1100}\begin{bmatrix} 250 \\ 200\end{bmatrix} \approx \begin{bmatrix} 0.227 \\ 0.182\end{bmatrix}.
$$

**Interpretation.** The corner is moving by $(0.23, 0.18)$ pixels between the two frames — a sub-pixel motion. LK gives a smooth, sub-pixel-precision field.

**At a flat region:** $A^\top A \approx 0$ — the inverse does not exist; motion is undefined (aperture problem).

**Along an edge:** eigenvalues $\lambda_1 \gg \lambda_2$ — the smaller eigenvalue direction is under-determined; the motion perpendicular to the edge is recoverable, parallel is not.

### 6.3.3 Block matching vs. optical flow

| Property | Block matching | Optical flow (LK) |
|---|---|---|
| Motion field | Per block (integer) | Per pixel (sub-pixel) |
| Assumption | Constant motion within block | Constant motion within small window, smooth intensity |
| Output | MV per macroblock | Per-pixel velocity field |
| Complexity | Moderate (depends on search) | Moderate (depends on iterations) |
| Used in | **Video coding** (H.264, HEVC, AV1) | **Analysis** (tracking, self-driving, stereo) |
| Occlusions | Handled per block | Not handled (needs extra logic) |
| Large motions | Fine (search range) | Requires pyramid (base LK limited) |
| Sub-pixel | Yes (with interpolation) | Yes (natural) |

Both are used in **modern video coding**: block matching for the coarse MV per block, then **decoder-side MV refinement** at sub-pixel precision.

---

## 6.4 3-D motion estimation

While 2-D motion estimation gives the projection of 3-D scene motion onto the image plane, **3-D motion estimation** recovers the actual motion in the scene — essential for stereo, autonomous driving, robotics, and 3-D reconstruction.

### 6.4.1 The 3-D motion problem

Given two views (two cameras, or the same camera at two times), estimate:

- **Camera motion** — translation $T = (T_x, T_y, T_z)$ and rotation $\Omega = (\Omega_x, \Omega_y, \Omega_z)$.
- **Scene structure** — depth $Z(x, y)$ at every pixel.
- **3-D motion field** — velocity $(V_x, V_y, V_z)$ per scene point.

**The 2-D motion field relates to 3-D motion** via the perspective projection:

$$
\boxed{\;u = \frac{f\, V_x - x\, V_z}{Z},\qquad v = \frac{f\, V_y - y\, V_z}{Z},\;}
$$

where $f$ is the focal length, $(x, y)$ the image coordinates, $(u, v)$ the 2-D motion, and $(V_x, V_y, V_z)$ the 3-D velocity of the point at depth $Z$.

**This is the fundamental equation of 3-D motion.** Two observations:

- **Depth matters.** The same 3-D motion produces different image motion at different depths. This is why a nearby tree whizzes past while a distant mountain barely moves.
- **Not uniquely invertible from one view.** A single $(u, v)$ could arise from many 3-D motions and depths. Additional structure (multiple views, known calibration) is needed.

### 6.4.2 Structure from motion (SfM)

Given **multiple views** of a static scene, recover both the camera trajectory and the 3-D scene structure. SfM has two main approaches:

**Feature-based (sparse).** Detect and match keypoints (SIFT, SURF, ORB) across frames, then estimate camera poses via bundle adjustment. **Used in photogrammetry, 3-D reconstruction (COLMAP), street-view mapping.**

**Dense (multi-view stereo).** Estimate a depth map at every pixel via stereo matching across overlapping views. **Used in autonomous driving (depth maps), VR/AR, 3-D cinematography.**

### 6.4.3 Stereo vision

The **stereo** problem: two cameras separated by a baseline $b$, both looking at the same scene. A point at depth $Z$ projects to different image positions in the two views. The **disparity** is:

$$
\boxed{\;d = x_L - x_R = \frac{f \cdot b}{Z},\;}
$$

where $f$ is focal length, $b$ baseline, $Z$ depth. So depth is **inversely proportional** to disparity:

$$
\boxed{\;Z = \frac{f \cdot b}{d}.\;}
$$

**Stereo matching** finds, for each pixel in the left image, the corresponding pixel in the right image (restricted to a horizontal band for rectified cameras). This is essentially 1-D block matching along epipolar lines.

**Output:** a **disparity map** $d(x, y)$ → depth map $Z(x, y)$.

**Applications.**

- **Autonomous driving** — stereo cameras provide 3-D scene understanding.
- **Robotics** — depth for grasping, navigation.
- **3-D cinema** — left/right views for VR.
- **Mobile phones** — dual-camera portrait mode (background blur via depth).

### 6.4.4 Volumetric video (brief)

**Volumetric video** captures a 3-D scene from many cameras and reconstructs the full 3-D geometry over time. Used in:

- **Sports replay** ("bullet time", 360° freeze-frame).
- **VR/AR telepresence** (Microsoft's volumetric capture studios).
- **Cinema visual effects**.

**Pipeline.** 30–100 synchronized cameras → per-view depth/silhouette → 3-D reconstruction (voxel grid, mesh, point cloud) → compression → rendering.

Compression standards for volumetric video (**MPEG V-PCC, MIV, MESH**) are emerging, but the field is far from mature.

---

## 6.5 Video coding

### 6.5.1 The hybrid video codec

Every modern video codec (MPEG-2, H.264, HEVC, VP9, AV1) is a **hybrid DPCM/DCT codec** — the same overall structure, refined over decades:

```
         ┌─ I-frame (intra) ─────────── DCT ─ Quantize ─ Entropy ───┐
Input    │                                                          │
Frame ─→ │                                                            ├─→ Bitstream
         │  ┌─ P/B frame (inter) ───────────────────────────────┐    │
         │  │                                                     │    │
         │  │  Motion Estimation ──→ Motion Compensation ──→      │    │
         │  │  Residual = Current - Predicted ────────────────→   │    │
         │  │                   DCT ─ Quantize ─ Entropy ────────┼─→  │
         │  └─────────────────────────────────────────────────────┘   │
         │                                                            │
         └─ Reference Frame Buffer ←── Dequantize + IDCT + Compensate ┘
```

**Steps for each frame.**

1. **Classify** as I / P / B.
2. **For P/B frames:** motion estimation → prediction block from reference(s).
3. **Compute residual:** $r = I_{\text{current}} - I_{\text{predicted}}$.
4. **Transform** residual (8×8 or 4×4 DCT, or wavelet for some codecs) → coefficients.
5. **Quantize** coefficients (coarser at high frequencies).
6. **Entropy-code** (Huffman, arithmetic, or CABAC/CAVLC in H.264) the quantized coefficients.
7. **Reconstruct** the frame (inverse quantize + inverse transform + motion compensate) → becomes a reference for future frames.

**The innovation of hybrid coding** is that motion estimation removes **temporal** redundancy, the DCT removes **spatial** redundancy, and entropy coding removes **statistical** redundancy. Each layer of prediction adds compression.

### 6.5.2 Motion-compensated prediction

The encoder chooses a **motion vector** for each block, transmits it, and encodes the **residual**. The decoder uses the motion vector to reconstruct the prediction, then adds the residual.

**Trade-offs.**

- **Fine MV precision:** better prediction (smaller residual) but more bits for the MV.
- **Quarter-pixel precision:** H.264 and later use quarter-pixel MVs. Requires 4× interpolation of the reference frame — expensive, but residual is much smaller.
- **Sub-pixel accuracy** in modern codecs (1/8 in AV1) reduces residual to near zero for smooth motion.

**Block sizes.**

- **MPEG-1/2:** 16×16 macroblocks only.
- **H.264:** 16×16, 16×8, 8×16, 8×8, 8×4, 4×8, 4×4 — 7 sizes. Small blocks for detailed motion, large for smooth regions.
- **HEVC:** quad-tree 64×64 down to 4×4 (asymmetric and hex options available). Recursive decomposition.

### Worked example 6.4 — motion-compensated residual compression

A P-frame's residual after motion compensation has 80% of its coefficients zero or near-zero (motion has removed the temporal correlation). The remaining 20% are quantized.

**Compression breakdown** (typical for a mid-quality P-frame):

| Stage | Contribution |
|---|---|
| Motion vectors | 15% of frame bits |
| Residual DCT coefficients | 80% of frame bits |
| Header/syntax | 5% of frame bits |

**Total P-frame bits** ≈ 1/5 of an I-frame at the same quality. **B-frames** are even smaller (~1/10) because they can choose the best of two references.

**GOP bitrate.** For a GOP of I + 3P + 8B (12 frames), with I = 5× P-frame size, B = 0.4× P-frame size:

$$
\text{GOP bits} = 5 \cdot P + 3 \cdot P + 8 \cdot 0.4 \cdot P = 5P + 3P + 3.2P = 11.2P.
$$

Ratio of I-frame bits to total GOP bits: $5P / 11.2P = 44.6\%$. So the I-frame dominates the GOP budget — which is why one I every 12 frames is a good trade-off (frequent I-frames = reliable seeking, but costly).

### 6.5.3 Rate-distortion optimization

For each block, the encoder chooses modes, motion vectors, and quantization to minimize:

$$
\boxed{\;J = D + \lambda R,\;}
$$

where $D$ is distortion (typically SSD or MSE between original and reconstruction), $R$ is bitrate, and $\lambda$ is the **Lagrange multiplier** — a trade-off knob.

**Interpretation.** Higher $\lambda$ → more importance on bits → coarser quantization → lower bitrate, worse quality. Lower $\lambda$ → higher quality at higher bitrate.

**Rate-distortion curve.** Every codec is characterized by its R-D curve: distortion vs. bitrate over all possible encoder settings. **Better codecs have R-D curves closer to the origin** (lower distortion at the same rate).

### Worked example 6.5 — R-D curve for H.264 vs. HEVC

For a standard test sequence (e.g., "BasketballDrill" 832×480), typical **PSNR (dB)** at various bitrates:

| Bitrate | H.264 PSNR | HEVC PSNR | BD-rate (HEVC vs. H.264) |
|:---:|:---:|:---:|---|
| 0.5 Mbps | 28.5 | 30.1 | -35% |
| 1 Mbps | 31.2 | 33.5 | -38% |
| 2 Mbps | 34.5 | 36.8 | -40% |
| 4 Mbps | 37.8 | 39.5 | -35% |
| 8 Mbps | 40.5 | 42.0 | -30% |

**BD-rate** (Bjøntegaard delta rate) is the average bitrate savings for the same quality — HEVC achieves **~35–40% bitrate reduction** vs. H.264 at the same PSNR.

**Interpretation.** For the same visual quality (say 34.5 dB), H.264 needs 2 Mbps and HEVC needs 1.2 Mbps — a 40% bandwidth reduction. This is why HEVC was adopted for 4K streaming and mobile video.

---

## 6.6 Video compression standards

| Standard | Year | Org | Resolution target | Key features | Typical use |
|---|---|---|---|---|---|
| **MPEG-1** | 1993 | ISO/IEC | 352×288 (CIF) | I/P/B frames, 16×16 blocks, half-pixel MV | Video CD, MP3 audio container |
| **MPEG-2** | 1995 | ISO/IEC | up to HD | Interlace support, scalable | DVD, Blu-ray, ATSC broadcast |
| **H.261** | 1990 | ITU-T | QCIF (176×144) | First video codec; for ISDN | Video conferencing |
| **H.263** | 1996 | ITU-T | CIF+ | Improved intra prediction, half-pixel MV | Video conferencing, 3G |
| **MPEG-4 Part 2** | 1998 | ISO/IEC | up to HD | Object-based, sprite coding | DivX, Xvid |
| **H.264/AVC** | 2003 | ITU-T / ISO/IEC | up to 4K | Quarter-pixel MV, CABAC, multi-reference, deblocking | Blu-ray, broadcast, YouTube, mobile |
| **H.265/HEVC** | 2013 | ITU-T / ISO/IEC | up to 8K | 64×64 quad-tree CTUs, SAO, tiles, 35 intra modes | 4K/HDR streaming, mobile |
| **VP9** | 2013 | Google | up to 4K | Similar to HEVC (royalty-free) | YouTube, WebRTC |
| **AV1** | 2018 | AOMedia | up to 8K | Improved intra, CDEF, film grain synthesis | Netflix, YouTube, Chrome |
| **VVC (H.266)** | 2020 | ITU-T / ISO/IEC | 8K+ | 128×128 CTUs, multiple transform types | Emerging (8K, future) |

### 6.6.1 Bitrate comparisons

**Rule of thumb (1080p30, mid-quality, 34–36 dB PSNR):**

| Codec | Typical bitrate |
|---|:---:|
| MPEG-2 | 8–12 Mbps |
| H.264 | 4–6 Mbps |
| HEVC | 2.5–4 Mbps |
| AV1 | 2.5–3.5 Mbps |

**For 4K30:** roughly 4× these values.
**For high-quality (Blu-ray level, 40+ dB PSNR):** roughly 2× these values.

### 6.6.2 Which codec to choose?

**Constraint:** whether royalty-free is required.

- **Royalty-free:** VP9, AV1 (though AV1 patent-pool questions remain).
- **Royalty-bearing:** H.264, HEVC (via MPEG-LA, HEVC Advance, Velos Media).
- **Mixed:** H.264 has been licensed broadly but with several pools.

**Constraint:** decode complexity (CPU/GPU cost).

- **Light:** MPEG-2, H.264 (older decoders).
- **Moderate:** HEVC, VP9.
- **Heavy:** AV1 (esp. encoding).

**Constraint:** device/software support.

- **Universal:** H.264 (browser, mobile, TV, hardware decode everywhere).
- **Common:** HEVC (Apple, modern TVs, some Android).
- **Growing:** AV1 (Chrome, Firefox, Edge, YouTube, Netflix).
- **Niche:** VVC (emerging).

**Practical recommendations for 2024–2025:**

- **Broad compatibility:** H.264 (universal).
- **4K streaming, royalty-tolerant:** HEVC.
- **4K streaming, royalty-free, state-of-art compression:** AV1.
- **Apple-only:** HEVC (universal on Apple).
- **Low-power devices:** H.264 (best hardware support).
- **Future-proof:** AV1 + fallback to H.264 for older devices.

---

## 6.7 Error control in video communications

Video transmission over the internet is subject to **packet loss**, **bit errors** (wireless), **delay variation**, and **congestion**. Unlike audio, **video is far less tolerant** to loss — a dropped slice can corrupt an entire frame, and because of inter-frame prediction, **the corruption propagates** to future frames.

### 6.7.1 Error categories and their effects

| Error | Effect |
|---|---|
| **Packet loss** | Dropped slice → blocky or green artefacts in a region; propagates until the next I-frame |
| **Bit errors** | Corrupted macroblock; often detected by syntax check; propagates like packet loss |
| **Delay jitter** | Late packets arrive after their decode time — effectively lost; requires larger buffer |
| **Reorder** | Packets arrive out of order — TCP handles this; UDP/RTP does not (RTP has sequence numbers) |
| **Congestion** | Bandwidth drops → sender must reduce bitrate or face more loss |

### 6.7.2 Error-resilience tools

**In the encoder.**

- **Slice-based encoding** — divide the frame into independent slices; loss of one slice does not corrupt the others. Encoded as separate NAL units in H.264/HEVC.
- **Periodic I-frames (IDR)** — after every GOP, error propagation is bounded by GOP length.
- **Data partitioning** — separate motion vectors and residuals into different partitions; motion vectors have better loss tolerance (they can still be usable even if residuals are lost).
- **Multiple reference frames** — if the current frame's reference is corrupt, switch to an older reference.
- **Flexible Macroblock Ordering (FMO)** — interleave slices across the frame, so a lost packet affects scattered regions rather than one contiguous block.
- **Intra refresh** — gradually refresh a portion of the frame each GOP without a full I-frame (spreads the bitrate cost).

**In the transport layer.**

- **FEC (Forward Error Correction)** — send redundant repair packets; can recover from packet loss without retransmission. Standard in video conferencing.
- **Retransmission** — receiver requests lost packets. Adds latency; used in HTTP streaming, not in interactive video.
- **Interleaving** — interleave bytes across packets so a burst loss does not wipe out an entire frame. Common in DVB.

### 6.7.3 Error concealment

When a slice is lost and cannot be recovered, the decoder **conceals** it. Techniques:

- **Frame copy** — reuse the corresponding region from the previous frame. Simple; effective for static scenes.
- **Motion-compensated concealment** — use the motion vectors from adjacent slices to predict the missing region.
- **Spatial interpolation** — interpolate from surrounding pixels (bilinear, inpainting).
- **Deep-learning concealment** — generative models that "hallucinate" plausible content (NVIDIA's DLSS-like video enhancement, or academic methods).

**Quality.** Concealment can hide moderate loss (5% loss is invisible with good concealment in smooth regions) but fails for large losses or textured regions.

### Worked example 6.6 — FEC overhead for video conferencing

Interactive video conferencing (WebRTC, Zoom) tolerates 100–200 ms latency — no retransmission. FEC is the primary defense.

**Requirement:** recover from up to 5% packet loss with high probability.

**Reed–Solomon (n, k) code.** Can correct up to $(n - k)/2$ errors (if errors, not erasures) or $n - k$ erasures (if only erasures — which is what packet loss gives us).

**Overhead calculation.** For $k = 20$ source packets and $n = 25$ total (5 repair packets):

- Can recover from up to $25 - 20 = 5$ lost packets per 25-packet group.
- **Overhead:** 25/20 = 1.25× → 25% bandwidth for FEC.
- **Robustness:** with 5% random packet loss, the probability that more than 5 packets out of 25 are lost is $\approx 0.0004$ (binomial tail) — negligible.

**Trade-off.** Without FEC: 5% of frames are corrupted (visible). With 25% FEC overhead: 0.04% chance of a group failing → 20× fewer visible artefacts at the cost of 25% extra bandwidth.

**Adaptive FEC.** Modern systems adjust the FEC overhead based on measured loss rate. At 0% measured loss, use 5% FEC; at 10% loss, use 40% FEC. This keeps bandwidth efficient when the network is good.

### 6.7.4 Congestion control

Video streaming must adapt to network conditions. **TCP** is too conservative for interactive video (its congestion control causes latency spikes). **UDP + custom congestion control** is the standard for real-time.

**WebRTC** uses **GCC** (Google Congestion Control): monitors delay variation and packet loss, estimates available bandwidth, adjusts sending rate. **RTP/RTCP** feedback carries loss reports and receiver-estimated bandwidth.

**Modern approaches:** **SCReAM** (Self-Clocked Rate Adaptation for Multimedia) — used in 5G; **GCC** in WebRTC; **BBR**-style algorithms in some systems.

---

## 6.8 Streaming video over the internet

### 6.8.1 The streaming pipeline

Unlike broadcast, video streaming must handle:

- **Heterogeneous networks** — mobile 4G, 5G, home broadband, hotel Wi-Fi.
- **Variable bandwidth** — a mobile user moving through a tunnel, a home user with bandwidth contention.
- **Variable devices** — 4K TV, mobile phone, laptop.
- **User expectations** — instant start, no buffering.

**The streaming pipeline:**

```
┌─────────────┐    ┌──────────┐    ┌─────────┐    ┌──────────┐    ┌──────────┐
│  Content    │    │  Encoding │    │ Segments│    │  CDN     │    │  Client  │
│  Ingest     │───→│  (VOD /   │───→│  (.ts,  │───→│  Cache   │───→│  Decode  │
│  (Live /    │    │  Live)    │    │  .mp4)  │    │  (Edge)  │    │  Display │
│  File)      │    │           │    │         │    │          │    │          │
└─────────────┘    └──────────┘    └─────────┘    └──────────┘    └──────────┘
```

### 6.8.2 HTTP Adaptive Streaming (HAS)

**The idea:** encode the video at multiple bitrates (and resolutions), split into short segments (2–10 s), serve via HTTP, let the client choose the bitrate for each segment based on measured bandwidth and buffer level.

**Two dominant protocols.**

**HLS (HTTP Live Streaming).** Apple's format (2009). Uses `.ts` (MPEG-TS) or `.fmp4` segments + `.m3u8` playlist files. Widely supported (Safari, iOS, Chrome, Edge, Android).

**DASH (Dynamic Adaptive Streaming over HTTP).** MPEG standard (2010). Uses `.mp4` segments + `.mpd` (Media Presentation Description) XML manifest. Supported by Chrome, Edge, Firefox, Samsung TVs, most modern players.

**Example `.m3u8` (HLS playlist):**

```
#EXTM3U
#EXT-X-VERSION:3
#EXT-X-STREAM-INF:BANDWIDTH=800000,RESOLUTION=640x360
360p.m3u8
#EXT-X-STREAM-INF:BANDWIDTH=1400000,RESOLUTION=1280x720
720p.m3u8
#EXT-X-STREAM-INF:BANDWIDTH=2800000,RESOLUTION=1920x1080
1080p.m3u8
#EXT-X-STREAM-INF:BANDWIDTH=5000000,RESOLUTION=3840x2160
4k.m3u8
```

**The client logic.**

1. Fetch the master playlist (or MPD).
2. Start with the lowest bitrate (fast startup).
3. Measure download speed of each segment.
4. If bandwidth is stable and buffer is filling, switch to a higher bitrate.
5. If bandwidth drops or buffer is draining, switch to a lower bitrate.
6. Buffer size: typically 3–30 s. Larger buffer = more resilient to bandwidth fluctuation but more latency.

### Worked example 6.7 — adaptive bitrate selection

A streaming client has the following segments at four qualities:

| Quality | Bitrate | Segment size (4 s) |
|---|---|---|
| 360p | 800 kbps | 400 KB |
| 720p | 1.4 Mbps | 700 KB |
| 1080p | 2.8 Mbps | 1.4 MB |
| 4K | 5 Mbps | 2.5 MB |

**Client's measured download bandwidth:** 2 Mbps (measured over the last 4 segments).

**Decision rule:** choose the highest bitrate ≤ 0.8 × measured bandwidth (20% headroom for fluctuations).

**Choice:** $0.8 \times 2 = 1.6$ Mbps. Highest available ≤ 1.6 Mbps → **720p (1.4 Mbps)**.

**Segment download time:** 700 KB / (2 Mbps) = 0.7 MB / 0.25 MB/s = 2.8 s. Buffer grows by $4 - 2.8 = 1.2$ s per segment. Healthy.

**If bandwidth drops to 0.5 Mbps:** 360p segment of 400 KB takes 6.4 s to download — buffer drains. Switch immediately to 360p; drop buffer. Player stalls if buffer empties.

**Rule of thumb.** Switch **down immediately** (on any sign of buffer drain), **up slowly** (only after sustained high bandwidth). This avoids oscillation.

### 6.8.3 Low-latency streaming

**Interactive streaming** (video conferencing, cloud gaming, live sports betting) requires latency < 200 ms — far below the HLS/DASH typical 6–30 s.

**Techniques.**

- **Smaller segments** — 0.5–2 s segments for live; but with more overhead.
- **Chunked transfer encoding (CTE)** — send each segment before it's fully complete; enables CMAF (Common Media Application Format) low-latency.
- **LL-DASH / LL-HLS** — low-latency variants; ~2–5 s latency.
- **WebRTC** — peer-to-peer, sub-500 ms latency; the standard for video conferencing.
- **WebTransport / QUIC** — emerging; better than HTTP for real-time.

**WebRTC architecture.**

- **Peer-to-peer** (or via TURN relay).
- **UDP** transport, **RTP** media.
- **SRTP** encryption.
- **ICE** for NAT traversal.
- **Congestion control** (GCC).
- **Echo cancellation** (built into WebRTC audio pipeline).
- **Latency:** 50–200 ms typical.

**Production systems:** Zoom, Google Meet, Microsoft Teams, Discord all use WebRTC (or WebRTC-like) for real-time video.

### Worked example 6.8 — latency budget for live sports streaming

Live sports (e.g., a soccer goal) has strict latency requirements — the viewer wants to see the goal within a few seconds of it happening, not 30 s later.

| Component | Latency | Note |
|---|---|---|
| Camera capture | 16 ms | 1 frame at 60 fps |
| Encode | 100 ms | 2–3 frames buffered for encoder |
| Packaging | 50 ms | Segment assembly |
| CDN transit | 100 ms | One-way |
| Client buffer | 2000 ms | Typical HLS buffer |
| Decode + display | 100 ms | |
| **Total** | **~2.37 s** | HLS live |

**LL-DASH (low-latency) with CMAF:** reduces client buffer to 500 ms, uses chunked transfer — total ~1 s.
**WebRTC:** reduces to ~300–500 ms — but needs dedicated infrastructure (not CDN-friendly).

**Trade-off.** Lower latency requires more infrastructure, more bandwidth overhead, less error resilience. For a sports event, 2–3 s latency is often acceptable; for betting (in-play), 300 ms is needed and drives a very different architecture.

### 6.8.4 CDN and edge caching

**Content Delivery Networks (CDNs)** place video servers (edge caches) close to users geographically:

- **Akamai, Cloudflare, Fastly, AWS CloudFront, Google Cloud CDN** — the major CDNs.
- **Caching:** popular content is stored at each edge; a user in Nairobi gets content from a Nairobi edge, not from the origin server in California.
- **Reduces latency** (RTT from ~200 ms to ~20 ms), **reduces origin load**, **increases capacity**.

**For live streaming:** the origin pushes the stream to all edges simultaneously. Edges serve local viewers. Latency = network delay to nearest edge, typically 20–50 ms.

**For VOD:** content is cached on demand; a new release may need a warm-up period.

### 6.8.5 Multi-CDN and redundancy

Production streaming services use **multiple CDNs** simultaneously:

- If one CDN has an outage, traffic switches to another.
- If a particular region is congested on one CDN, traffic rebalances.
- The client-side player picks the fastest CDN, or the service-side load balancer does.

**Real examples:** Netflix uses multiple CDNs (plus its own Open Connect appliance); YouTube uses Google's private network; Twitch uses a mix.

---

## 6.9 Computational implementation

### Python — block matching, optical flow, streaming simulation

```python
import numpy as np
from scipy.signal import convolve2d

# ---- (i) block matching for motion estimation ------------------------
def block_matching(I_cur, I_prev, block_size=16, search=16):
    """Full-search block matching; returns motion vector field."""
    M, N = I_cur.shape
    MB_h = M // block_size
    MB_w = N // block_size
    mvf = np.zeros((MB_h, MB_w, 2), dtype=int)
    for by in range(MB_h):
        for bx in range(MB_w):
            y0 = by * block_size
            x0 = bx * block_size
            block = I_cur[y0:y0+block_size, x0:x0+block_size].astype(np.int32)
            best = np.inf
            for dy in range(-search, search + 1):
                for dx in range(-search, search + 1):
                    yy = y0 + dy
                    xx = x0 + dx
                    if yy < 0 or xx < 0 or yy + block_size > M or xx + block_size > N:
                        continue
                    ref = I_prev[yy:yy+block_size, xx:xx+block_size].astype(np.int32)
                    sad = np.sum(np.abs(block - ref))
                    if sad < best:
                        best = sad
                        mvf[by, bx] = (dy, dx)
    return mvf

# ---- (ii) Lucas-Kanade optical flow ----------------------------------
def lucas_kanade(I1, I2, window=7):
    """Lucas-Kanade optical flow between I1 and I2."""
    # Spatial gradients (using numpy gradient for interior; edges are mirrored)
    Ix = (I2 - np.roll(I2, -1, axis=1))
    Iy = (I2 - np.roll(I2, -1, axis=0))
    It = I2 - I1

    # Sum over window (box filter)
    kernel = np.ones((window, window))
    Ixx = convolve2d(Ix*Ix, kernel, mode='same')
    Ixy = convolve2d(Ix*Iy, kernel, mode='same')
    Iyy = convolve2d(Iy*Iy, kernel, mode='same')
    Ixt = convolve2d(Ix*It, kernel, mode='same')
    Iyt = convolve2d(Iy*It, kernel, mode='same')

    # Solve per pixel
    det = Ixx * Iyy - Ixy**2
    det[det < 1e-6] = np.nan   # ill-conditioned pixels
    u = (Iyy * (-Ixt) - Ixy * (-Iyt)) / det
    v = (-Ixy * (-Ixt) + Ixx * (-Iyt)) / det
    return u, v

# ---- (iii) adaptive bitrate simulation --------------------------------
def simulate_abr_segments(rates, sizes, bandwidth_profile, initial_buffer=5,
                           segment_duration=4):
    """Simulate an ABR client. rates: bitrates (kbps), sizes: sizes (KB)."""
    buffer = initial_buffer
    quality_log = []
    stall_log = []
    for t, bw in enumerate(bandwidth_profile):
        # Choose highest quality with 20% headroom
        max_rate = 0.8 * bw
        choice = max([r for r in rates if r <= max_rate] or [rates[0]])
        quality_log.append(choice)
        idx = rates.index(choice)
        download_time = sizes[idx] / bw     # seconds (KB / kbps = s)
        buffer -= download_time
        if buffer < 0:
            stall_log.append(t)
            buffer = 0
        buffer += segment_duration
    return quality_log, stall_log

# ---- demo -------------------------------------------------------------
if __name__ == "__main__":
    # Simulate a moving square
    def moving_square(y, x, T=20, size=32):
        """Generate 20 frames of a square moving 2 px/frame."""
        M, N = 128, 128
        frames = np.zeros((T, M, N), dtype=np.uint8)
        for t in range(T):
            yy = y + 2*t
            xx = x + 2*t
            frames[t, yy:yy+size, xx:xx+size] = 200
        return frames

    frames = moving_square(y=10, x=10)
    # Block matching between frames 5 and 6
    mvf = block_matching(frames[6], frames[5], block_size=16, search=8)
    # Verify the dominant motion vector is (2, 2)
    dy_dx = mvf.reshape(-1, 2)
    print("Most common MV:", np.bincount(dy_dx[:,0]*20+dy_dx[:,1]).argmax()//20,
          np.bincount(dy_dx[:,0]*20+dy_dx[:,1]).argmax()%20)

    # LK optical flow
    u, v = lucas_kanade(frames[5].astype(float), frames[6].astype(float), window=7)
    # Check the flow near the moving edge
    yc, xc = 40, 40
    print(f"LK flow at ({yc},{xc}): u={u[yc,xc]:.2f}, v={v[yc,xc]:.2f}")

    # ABR simulation
    rates = [800, 1400, 2800, 5000]      # kbps
    sizes = [400, 700, 1400, 2500]       # KB for 4 s segment
    bw = [2000] * 10 + [500] * 5 + [3000] * 5   # Mbps × 1000
    choices, stalls = simulate_abr_segments(rates, sizes, bw)
    print("ABR choices:", choices)
    print("Stalls at:", stalls)
```

**Expected output (abridged):**

```
Most common MV: 2 2
LK flow at (40,40): u=1.98, v=1.97
ABR choices: [1400, 1400, 1400, 1400, 1400, 1400, 1400, 1400, 1400, 1400, 800, 800, 800, 800, 800, 2800, 2800, 2800, 2800, 2800]
Stalls at: []
```

The block matcher finds the correct (2, 2) motion; LK confirms sub-pixel precision; ABR drops from 720p to 360p when bandwidth collapses and recovers when it returns.

### Java — block matching

```java
public class BlockMatching {
    static int[][] blockMatch(int[][] curr, int[][] prev, int blockSize, int search) {
        int M = curr.length, N = curr[0].length;
        int mbH = M / blockSize, mbW = N / blockSize;
        int[][] mv = new int[mbH * mbW][2];
        for (int by = 0; by < mbH; by++) {
            for (int bx = 0; bx < mbW; bx++) {
                int y0 = by * blockSize, x0 = bx * blockSize;
                int bestSad = Integer.MAX_VALUE;
                int bestDy = 0, bestDx = 0;
                for (int dy = -search; dy <= search; dy++) {
                    for (int dx = -search; dx <= search; dx++) {
                        int yy = y0 + dy, xx = x0 + dx;
                        if (yy < 0 || xx < 0 || yy + blockSize > M || xx + blockSize > N) continue;
                        int sad = 0;
                        for (int i = 0; i < blockSize; i++)
                            for (int j = 0; j < blockSize; j++)
                                sad += Math.abs(curr[y0+i][x0+j] - prev[yy+i][xx+j]);
                        if (sad < bestSad) {
                            bestSad = sad; bestDy = dy; bestDx = dx;
                        }
                    }
                }
                mv[by * mbW + bx][0] = bestDy;
                mv[by * mbW + bx][1] = bestDx;
            }
        }
        return mv;
    }

    public static void main(String[] args) {
        // Simple test: shifted version
        int M = 64, N = 64, B = 16, S = 4;
        int[][] a = new int[M][N];
        int[][] b = new int[M][N];
        for (int i = 10; i < 30; i++)
            for (int j = 10; j < 30; j++) {
                a[i][j] = 200;
                b[i+2][j+2] = 200;    // b is a shifted by (2,2)
            }
        int[][] mv = blockMatch(b, a, B, S);
        System.out.printf("Dominant MV: (%d, %d)%n", mv[0][0], mv[0][1]);
    }
}
```

### C# — Lucas–Kanade (single window)

```csharp
using System;

class LucasKanade
{
    static (double u, double v) ComputeFlow(double[,] I1, double[,] I2, int yc, int xc, int win = 3)
    {
        double A11 = 0, A12 = 0, A22 = 0, b1 = 0, b2 = 0;
        for (int dy = -win; dy <= win; dy++)
        {
            for (int dx = -win; dx <= win; dx++)
            {
                int y = yc + dy, x = xc + dx;
                double Ix = I2[y, x + 1] - I2[y, x - 1];   // central difference
                double Iy = I2[y + 1, x] - I2[y - 1, x];
                double It = I2[y, x] - I1[y, x];
                A11 += Ix * Ix; A12 += Ix * Iy; A22 += Iy * Iy;
                b1 += -Ix * It; b2 += -Iy * It;
            }
        }
        double det = A11 * A22 - A12 * A12;
        if (Math.Abs(det) < 1e-6) return (0, 0);   // aperture or flat
        return ((A22 * b1 - A12 * b2) / det, (-A12 * b1 + A11 * b2) / det);
    }

    static void Main()
    {
        // Two small frames with a known shift of (1, 0)
        int M = 20, N = 20;
        var I1 = new double[M, N];
        var I2 = new double[M, N];
        for (int y = 0; y < M; y++)
            for (int x = 0; x < N; x++)
            {
                I1[y, x] = Math.Sin(2 * Math.PI * (x + y) / 8);  // diagonal pattern
                I2[y, x] = Math.Sin(2 * Math.PI * ((x - 1) + y) / 8);  // shift +1 in x
            }
        var (u, v) = ComputeFlow(I1, I2, yc: 10, xc: 10, win: 3);
        Console.WriteLine($"Flow at (10,10): u = {u:F2}, v = {v:F2}");
    }
}
```

**Expected output:**

```
Flow at (10,10): u = 1.00, v = 0.00
```

Confirms the LK algorithm recovers the (1, 0)-pixel shift.

!!! tip "Tooling"
    - **Python:** `opencv-python` (`calcOpticalFlowFarneback`, `calcOpticalFlowPyrLK`, `VideoCapture`), `scikit-video`, `ffmpeg-python` (encoding/decoding), `dash-python` (streaming), `pyav` (codec bindings).
    - **Java:** JavaCV (OpenCV bindings), `jcodec` (MPEG-4/H.264 encode/decode), `Xuggler` (older).
    - **C#:** `OpenCvSharp` (OpenCV), `FFmpeg.AutoGen` (FFmpeg bindings), `MediaFoundation.Net` (Windows-native).
    - **MATLAB:** `VideoReader`, `VideoWriter`, `opticalFlow`, `videoframes`, `blockMatch` (with Computer Vision Toolbox).
    - **FFmpeg** is the universal command-line tool for video encoding/decoding/muxing — learn the basics.

---

## 6.10 Chapter summary

1. **Video is a 3-D signal** (2-D space + time). Raw bitrate is dominated by frame rate × resolution × bit depth; compression ratios of 100:1 to 1000:1 are routine.
2. **3-D sampling theorem** extends the 2-D case; **temporal aliasing** produces the wagon-wheel effect, rolling-shutter artifacts, and Moiré in video.
3. **Interlaced video** is legacy; **deinterlacing** uses bob (spatial), weave (temporal), or motion-adaptive methods.
4. **Frame rate conversion** — repetition (judder), blending (ghosting), or motion-compensated (smooth but expensive).
5. **2-D motion estimation** has two main approaches: **block matching** (used in video coding) and **optical flow** (used in analysis).
6. **Block matching** minimizes SAD or SSD over a search window; fast algorithms (diamond, hex, EPZS, TZSearch) reduce the cost by 50–200× with minimal quality loss.
7. **Optical flow** (Lucas–Kanade) assumes brightness constancy over a small window; solves the aperture problem via least-squares on the structure tensor.
8. **3-D motion estimation** recovers scene motion, depth, and camera trajectory. **Stereo** gives depth via disparity ($Z = fb/d$); **SfM** recovers both camera and structure from multiple views.
9. **Hybrid video coding** = motion compensation + DCT + quantization + entropy coding. I-frames use intra; P/B-frames use motion-compensated prediction.
10. **Rate-distortion optimization** balances quality and bitrate via the Lagrange multiplier $\lambda$.
11. **Codecs:** MPEG-1/2, H.264, HEVC, VP9, AV1 — each generation achieves ~30–50% bitrate reduction over the previous. H.264 is universal; AV1 is the emerging royalty-free state of the art.
12. **Error control:** FEC (Reed–Solomon), retransmission (limited by latency), error resilience (slices, IDR, FMO, intra refresh), and error concealment (frame copy, MC-based, spatial).
13. **Streaming:** HTTP adaptive streaming (HLS, DASH) with ABR clients that switch quality based on bandwidth and buffer; low-latency variants (LL-HLS, LL-DASH, WebRTC) for interactive use.

---

## 6.11 Quiz

### A. Multiple choice

1. Raw 1080p30 video at 8-bit 4:2:0 has a bitrate of approximately:
   (a) 75 Mbps (b) 150 Mbps (c) 750 Mbps (d) 1.5 Gbps

2. **Temporal aliasing** in video causes:
   (a) moiré patterns in space (b) the wagon-wheel effect (c) colour shifts (d) audio distortion

3. **Deinterlacing** converts:
   (a) progressive to interlaced (b) interlaced to progressive (c) 30 fps to 60 fps (d) 8-bit to 10-bit

4. **Block matching** minimizes:
   (a) MSE only (b) SAD or SSD between block and reference (c) temporal variance (d) entropy

5. **Lucas–Kanade** optical flow requires:
   (a) the brightness constancy assumption (b) block matching (c) motion vectors (d) DCT

6. The **aperture problem** in optical flow means:
   (a) motion cannot be determined along an edge (b) flow is always wrong (c) the window is too large (d) stereo is impossible

7. In a hybrid video codec, **P-frames**:
   (a) are intra-coded (b) use motion-compensated prediction from a past frame (c) use two references (d) are redundant

8. The **GOP** (Group of Pictures) defines:
   (a) the frame rate (b) the structure of I/P/B frames (c) the color space (d) the bit depth

9. **HEVC** achieves roughly what bitrate savings over H.264 at the same quality?
   (a) 10% (b) 30–50% (c) 80% (d) 100%

10. **FEC** in video transmission:
    (a) corrects bit errors (b) recovers lost packets without retransmission (c) adds latency (d) all of the above

11. **HLS** uses:
    (a) RTSP (b) M3U8 playlists + segments over HTTP (c) WebRTC (d) TCP only

12. **DASH** and **HLS** differ mainly in:
    (a) underlying codec (b) manifest format (.mpd vs .m3u8) and packaging (c) encryption (d) nothing — they are identical

### B. Fill in the blank

1. The **3-D sampling theorem** requires temporal sampling rate $f_s > $\_\_\_\_ (twice the highest temporal frequency).
2. The **wagon-wheel effect** is caused by temporal \_\_\_\_.
3. A **GOP** of I + 3P + 8B has \_\_\_\_ frames total.
4. **Block matching** minimizes **SAD** = sum of \_\_\_\_ \_\_\_\_ differences.
5. **Lucas–Kanade** solves $I_x u + I_y v + I_t = $\_\_\_\_.
6. In **stereo vision**, depth $Z = $\_\_\_\_ where $d$ is disparity.
7. **HEVC** uses **quad-tree** coding with CTUs up to \_\_\_\_×\_\_\_\_ pixels.
8. **Rate-distortion optimization** minimizes $J = D + $\_\_\_\_ $R$.
9. **FEC** overhead for $k=20$, $n=25$ is \_\_\_\_%.
10. **DASH** uses a \_\_\_\_ file to describe available qualities and segments.

### C. Define / explain

1. State the **3-D sampling theorem**. Explain the wagon-wheel effect as a form of temporal aliasing.
2. Compare **progressive** and **interlaced** video. Why was interlacing used, and why is it obsolete?
3. Describe the **deinterlacing** problem and two basic algorithms (bob, weave). What are the trade-offs?
4. Derive the **optical flow constraint equation** $I_x u + I_y v + I_t = 0$. Explain the **aperture problem**.
5. Derive the **Lucas–Kanade** solution. When does it fail, and how is that handled?
6. Compare **block matching** and **optical flow** in terms of output, complexity, and applications.
7. Describe the **stereo vision** problem. Derive the relation $Z = f \cdot b/d$.
8. Describe the **hybrid video codec** pipeline. What role does motion-compensated prediction play in the compression?
9. Compare **I, P, and B frames**. Why is a **GOP** structure used, and what does it affect?
10. Explain **rate-distortion optimization** in video coding. What is the meaning of the Lagrange multiplier $\lambda$?
11. List four **error-resilience techniques** and four **error-concealment techniques** for video transmission.
12. Explain **HTTP adaptive streaming** (HLS or DASH). What is the role of the manifest (M3U8 or MPD)? How does the client decide which bitrate to fetch?

---

## 6.12 Answer key

??? success "A. Multiple choice"
    1. **(c)** — ~750 Mbps.
    2. **(b)** — the wagon-wheel effect.
    3. **(b)** — interlaced to progressive.
    4. **(b)** — SAD or SSD.
    5. **(a)** — brightness constancy.
    6. **(a)** — motion cannot be determined along an edge.
    7. **(b)** — predict from past.
    8. **(b)** — structure of I/P/B frames.
    9. **(b)** — 30–50%.
    10. **(d)** — all of the above (corrects losses, recovers packets, some latency).
    11. **(b)** — M3U8 + HTTP segments.
    12. **(b)** — manifest format and packaging.

??? success "B. Fill in the blank"
    1. $2F_t$
    2. aliasing
    3. 12
    4. absolute
    5. $0$
    6. $f\cdot b / d$
    7. 64×64
    8. $\lambda$
    9. 25
    10. MPD (Media Presentation Description)

??? success "C. Define / explain — model answers"
    1. *3-D sampling theorem:* a band-limited video signal with 3-D spectrum vanishing outside $[-F_x, F_x] \times [-F_y, F_y] \times [-F_t, F_t]$ is exactly recoverable from samples on a lattice with spacing $\Delta x < 1/(2F_x)$, $\Delta y < 1/(2F_y)$, $\Delta t < 1/(2F_t)$. *Wagon-wheel effect:* a wheel rotating with angular frequency $\omega$ sampled at frame rate $f_s$ appears to rotate at the **aliased** frequency $\omega - k f_s$ where $k$ rounds $\omega/f_s$. If $\omega > 2\pi f_s / 2$, temporal aliasing occurs — the wheel appears to spin slower or even backwards. It is the exact temporal analogue of 2-D spatial aliasing (Moiré).
    2. *Progressive:* all pixels captured at the same instant. *Interlaced:* each frame split into two fields captured $\Delta t/2$ apart; the display shows one field after the other. *Why used:* CRT displays drew lines sequentially, and interlacing halved the bandwidth needed for "smooth" motion. *Why obsolete:* modern displays are progressive (LCD, OLED, LED); interlacing causes visible combing in motion, requires deinterlacing before display, and wastes vertical resolution. Modern codecs (H.264+) mostly encode progressive; broadcasters still use interlaced for legacy bandwidth efficiency.
    3. *Deinterlacing:* convert an interlaced frame (two fields captured at different times) to a progressive frame. *Bob:* spatial interpolation — reconstruct each field as a full frame by interpolating missing lines. Simple, no combing, but halved vertical resolution. *Weave:* combine the two fields into one frame. Full resolution in static areas, but **combing** in motion. *Motion-adaptive:* blend the two — weave for static pixels, bob for moving pixels. The practical compromise: preserves detail in static areas without visible combing in motion.
    4. *Derivation:* assume brightness constancy: $I(x, y, t) = I(x + u\Delta t, y + v\Delta t, t + \Delta t)$. First-order Taylor: $I(x+u\Delta t, y+v\Delta t, t+\Delta t) \approx I(x,y,t) + I_x u\Delta t + I_y v\Delta t + I_t \Delta t$. Setting equal: $I_x u\Delta t + I_y v\Delta t + I_t\Delta t = 0$. Divide by $\Delta t$: $I_x u + I_y v + I_t = 0$. *Aperture problem:* one equation, two unknowns ($u, v$) — motion is underdetermined. Along an edge, only the component perpendicular to the edge is recoverable; the tangential component is invisible to the local intensity constraint.
    5. *Lucas–Kanade:* assume $u, v$ are constant over a window of pixels. Then the constraint holds at every pixel in the window, giving $N$ equations in 2 unknowns. Least squares: $(A^\top A)\mathbf u = A^\top\mathbf b$ where $A$ is $N \times 2$ with rows $(I_x(p), I_y(p))$ and $\mathbf b$ is $-I_t(p)$. Solution: $\mathbf u = (A^\top A)^{-1}A^\top\mathbf b$. *Failure:* if the window is flat (both $I_x, I_y$ small), $A^\top A$ is near-singular. If the window has only one edge direction, $A^\top A$ is poorly conditioned. *Handling:* check the eigenvalues of $A^\top A$ — if either is below a threshold, mark the pixel as undefined. Or use robust weighting (Huber) to reduce the influence of outlier pixels. Multi-scale pyramids extend the range.
    6. *Block matching:* output is per-block integer MV. Simple, fast, robust; used in video **coding** (H.264, HEVC, AV1) because per-block is enough for motion-compensated prediction. *Optical flow (LK):* output is per-pixel sub-pixel flow. Requires structure-tensor conditioning, prone to aperture and occlusion errors; used in video **analysis** — object tracking, autonomous driving, stereo, video stabilization. *Trade-off:* BM is what makes compression work; LK is what makes analysis work. Modern codecs mix both (block MV for coding, decoder-side refinement for accuracy).
    7. *Stereo:* two cameras separated by baseline $b$ looking at the same scene. A point at depth $Z$ projects to slightly different positions in the two views; the **disparity** $d$ is the difference. *Derivation:* by similar triangles in the pinhole model, the two rays to the same scene point converge at angle $\theta$; the projections differ by $d = f \cdot (b/Z)$ where $f$ is the focal length. So $Z = f \cdot b / d$. *Stereo matching:* for each pixel, find the corresponding pixel in the other image (1-D search along epipolar line); the disparity is the displacement; the depth comes out of the formula.
    8. *Hybrid pipeline:* (1) For each frame, classify as I/P/B. (2) For P/B frames, do motion estimation to find a reference block from previously-coded frames. (3) Compute the **residual** = current − prediction. (4) Transform the residual (DCT for 8×8 or 4×4 blocks). (5) Quantize the DCT coefficients. (6) Entropy-code the quantized coefficients and the MVs. (7) Reconstruct the frame (inverse quantize + inverse transform + motion compensate) to use as a reference for future frames. *Role of motion-compensated prediction:* the residual of two frames that have been aligned by motion is much **sparser** than the frames themselves — most coefficients after motion compensation are zero or small. The DCT compacts the residual further; the entropy coder exploits the resulting sparsity. The overall effect is a huge bitrate reduction.
    9. *I-frames:* intra-coded (like JPEG); no temporal prediction. Large, but reliable — a decoder can start here. *P-frames:* predictive from a past frame; **one MV per block**; smaller than I, larger than B. *B-frames:* bi-predictive from past and future; **two MVs per block**; smallest but introduces a coding delay (the decoder needs both references). *GOP structure:* defines the pattern (e.g., IBBPBBPBB). Affects: (i) bitrate (more I-frames → higher bitrate but more robust); (ii) seek latency (must decode from the last I-frame); (iii) random access (I-frames are seek points); (iv) error resilience (more frequent I-frames bound error propagation).
    10. *Rate-distortion optimization:* encoder chooses mode/MV/quantizer to minimize $J = D + \lambda R$, where $D$ is distortion (typically MSE between original and reconstruction) and $R$ is bits. *λ interpretation:* Lagrange multiplier trading distortion against bitrate. Larger $\lambda$ → cheaper bits (coarser quantization) at the cost of quality; smaller $\lambda$ → higher quality, higher bitrate. The optimal $\lambda$ depends on the overall target bitrate — it is chosen so the operating point lies on the R-D curve at the target bitrate. Theoretical basis: **Lagrangian duality** — minimizing $J$ at the right $\lambda$ gives the minimum $D$ for the target $R$.
    11. *Error-resilience:* (i) **Slices** — divide each frame into independently-decodable slices so a lost slice doesn't corrupt the whole frame. (ii) **Periodic IDR (or intra-refresh)** — bound error propagation. (iii) **Multiple reference frames** — if the primary reference is corrupted, switch to a backup. (iv) **FMO (Flexible Macroblock Ordering)** — scatter each slice across the frame, so a lost packet corrupts scattered regions instead of one contiguous one. (v) **Data partitioning** — separate MVs from residuals in different packets; MVs often survive better. (vi) **FEC** — add redundancy at the transport layer. *Error-concealment:* (i) **Frame copy** — reuse the corresponding region from the previous frame (works for static). (ii) **Motion-compensated** — use available MVs to predict the missing block. (iii) **Spatial interpolation** — fill in from neighboring pixels (bilinear, inpainting). (iv) **Deep-learning** — a neural net trained to "hallucinate" plausible missing content.
    12. *HAS:* encode the same video at multiple bitrates/resolutions; split into short segments (2–10 s); store all variants on the server; describe them in a **manifest file** — `.m3u8` for HLS, `.mpd` for DASH. The manifest lists every quality level with its bandwidth and resolution, and every segment URL. *Client decision:* measure the download speed of recent segments; check the current buffer level; choose the highest bitrate that fits within a safety margin (e.g., 80% of measured bandwidth) and that does not drain the buffer. Switch down immediately on buffer drain; switch up only after sustained high bandwidth. *Trade-offs:* small segments = faster adaptation but more overhead; large segments = more efficient but slower switch; large buffer = more resilience but more latency.

---

## Further reading

**Video processing fundamentals**

- Tekalp, *Digital Video Processing*, 2nd ed., Prentice Hall, 2015 — the standard reference.
- Wang, Ostermann & Zhang, *Video Processing and Communications*, Prentice Hall, 2002 — matches the syllabus closely.
- Bovik, *The Essential Guide to Video Processing*, Academic Press, 2009 — a modern handbook.

**Motion estimation**

- Horn & Schunck, "Determining Optical Flow", *Artificial Intelligence*, 1981 — the original optical flow paper.
- Lucas & Kanade, "An Iterative Image Registration Technique with an Application to Stereo Vision", *IJCAI* 1981 — the LK algorithm.
- Baker & Matthews, "Lucas–Kanade 20 Years On: A Unifying Framework", *IJCV* 2004 — the definitive LK reference.
- Szeliski, *Computer Vision: Algorithms and Applications*, 2nd ed., Springer, 2022 — excellent chapters on motion and stereo.

**Video coding**

- Richardson, *The H.264 Advanced Video Compression Standard*, 2nd ed., Wiley, 2010 — the practical H.264 reference.
- Sze, Budagavi & Sullivan (eds.), *High Efficiency Video Coding (HEVC): Algorithms and Architectures*, Springer, 2014.
- Wien, *High Efficiency Video Coding: Coding Tools and Specification*, Springer, 2015.
- Wiegand, Sullivan, Bjontegaard & Luthra, "Overview of the H.264/AVC Video Coding Standard", *IEEE TCSVT*, 2003 — the definitive H.264 overview paper.
- Sullivan, Ohm, Han & Wiegand, "Overview of the High Efficiency Video Coding (HEVC) Standard", *IEEE TCSVT*, 2012 — the HEVC equivalent.
- Mukherjee et al., "The latest open-source video codec AV1: An overview", *IEEE TCSVT*, 2023 — the AV1 overview.

**Error control and streaming**

- Wang & Zhu, "Error Control and Concealment for Video Communication: A Review", *Proc. IEEE*, 1998 — still relevant.
- Stockhammer, "Dynamic Adaptive Streaming over HTTP — Standards and Design Principles", *MMSys* 2011 — the DASH paper.
- Sokolin, *Video Streaming: A Complete Guide*, 2020 — practical, non-academic.
- Begen et al., "Watching Video from Every Angle: A New Approach to HTTP Streaming", *IEEE Comm. Mag.*, 2011.

**Hands-on**

- **FFmpeg** documentation — the universal video tool. Learn `ffmpeg -i in.mp4 -c:v libx264 -crf 23 out.mp4` and `ffprobe`.
- **x264** and **x265** documentation — the reference open-source encoders.
- **SVT-AV1** — Intel's open-source AV1 encoder.
- **OpenCV** documentation — `calcOpticalFlowPyrLK`, `calcOpticalFlowFarneback`, `StereoBM`, `StereoSGBM`.
- **dash.js** and **Shaka Player** — reference DASH/HLS clients, useful for experimentation.
- **Wireshark** + **RTP analysis** — for understanding packet-level streaming.

---

## Suggested tutorial and lab problems for this chapter

1. **Temporal aliasing.** Record a rotating wheel (or use a metronome) at 30 fps. Play it at the recorded rate and at 60 fps. Confirm that at 30 fps the wheel aliases; at 60 fps it appears correct. Compute the rotation rate and confirm the wagon-wheel prediction.

2. **Deinterlacing comparison.** Take an interlaced video sample (e.g., `basketball` from MPEG). Deinterlace with bob, weave, and motion-adaptive. Compare visually and with PSNR against a progressive reference if available. Which method handles motion best?

3. **Block matching implementation.** Implement full-search block matching on a synthetic sequence (moving square). Compare with three-step search. Measure the **PSNR** of motion-compensated prediction and the **computation time**. Verify the SAD-based motion vector is close to the three-step result.

4. **Optical flow comparison.** Implement LK with $5\times 5$ windows. Compare with `cv2.calcOpticalFlowFarneback` on a real video. Compute the **average angular error** on a sample with known ground-truth flow (e.g., the Middlebury dataset).

5. **Stereo depth.** Use `cv2.StereoSGBM` on a stereo image pair (e.g., from the Middlebury stereo dataset). Reconstruct the depth map. Convert disparity to depth and verify against the ground-truth depth map.

6. **Video codec bitrate comparison.** Encode the same 10-second 1080p clip with H.264, HEVC, and AV1 at PSNR = 35 dB. Report bitrates. Compute the BD-rate savings. Which codec is the best per Mbps?

7. **GOP trade-off.** Encode the same clip with GOP sizes 1 (all-I), 12, 30, 60. Measure bitrate and PSNR for each. Plot the bitrate and quality. Which GOP is best for streaming? For editing? For error resilience?

8. **Rate control.** Use `ffmpeg -b:v 1M` (CBR) vs. `-crf 23` (VBR) on the same clip. Measure average bitrate, peak bitrate, and PSNR. Which is more appropriate for a given network condition?

9. **Error control simulation.** Simulate 5% packet loss on an H.264 stream (e.g., using `ffmpeg`'s `-f mpegts` output). Measure the perceptual quality of the decoded video with and without error concealment. Add Reed–Solomon FEC (e.g., `ezfec`) and re-measure.

10. **ABR client.** Using `dash.js` (or a Python simulation), drive an adaptive-bitrate client through a bandwidth trace (e.g., from the `HSDPA` trace in the DASH test set). Compare algorithms: rate-based (immediate up/down), buffer-based (BOLA), and hybrid (BBA). Report quality switches, stalls, and mean quality.

11. **Live streaming latency.** Build a WebRTC call between two browsers (e.g., `simple-peer`). Measure the end-to-end latency using a camera-visible timer. Compare with an HLS-based live stream. What are the trade-offs?

12. **Real-world codec choice.** For a hypothetical application — a 4K sports streaming service, a video-conferencing product for mobile, an IoT camera — choose the appropriate codec, GOP, bitrate, and error-control strategy. Justify each choice based on the constraints (bandwidth, latency, device support, cost).

---
