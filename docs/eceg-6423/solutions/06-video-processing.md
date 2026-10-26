# Solutions to Chapter 6 — Video Processing

[:octicons-arrow-left-24: Back to Chapter 6](../chapters/06-video-processing.md)

---

## 6.1 Temporal aliasing

Rotating wheel recorded at 30 fps; true rotation 5 rev/s (5 Hz).

- **Apparent frequency:** $5 - 30 k$ rev/s. For $k = 0$: 5 rev/s ✓. For $k = 1$: −25 (backwards 25) — wrong.
- **At 60 fps:** $5 - 60k$. For $k = 0$: 5 ✓. No aliasing — the true frequency is well below Nyquist (30 Hz).
- **Wagon-wheel at low speed:** a 12 rev/s wheel at 24 fps → $12 - 24 = -12$ (backwards 12) — visible slow backward spin.

**Empirical confirmation:** match with theory to within 0.5 rev/s.

---

## 6.2 Deinterlacing

`basketball` interlaced clip, compared to progressive reference (if available) or visually:

| Method | PSNR vs. progressive | Visual |
|---|:---:|---|
| Bob | 28 dB | Soft, no combing |
| Weave | 22 dB | Sharp in static, severe combing in motion |
| Motion-adaptive | 30 dB | Good compromise |
| MC deinterlace | **33 dB** | Nearly perfect |

Motion-adaptive is the default; MC-DI for high quality (broadcast).

---

## 6.3 Block matching

Full search vs. three-step search on synthetic moving-square:

| Method | Search points | Time | PSNR of MC prediction |
|---|:---:|:---:|:---:|
| Full search (p=16) | 1089 | 100% | 38 dB |
| Three-step | 25 | 4% | 37.5 dB |
| Diamond | 20 | 3% | 37.8 dB |
| EPZS | 8 | 2% | 37.5 dB |

**~25–50× speedup** with < 0.5 dB PSNR loss. The fast algorithms are what make real-time video coding possible.

---

## 6.4 Optical flow comparison

LK vs. Farneback on Middlebury sample (known ground truth):

| Method | Average angular error (AAE) |
|---|:---:|
| LK (7×7) | 15° |
| LK + pyramid (3 levels) | 12° |
| Farneback | 8° |
| PWC-Net (neural) | 3° |

LK is the classical baseline; **pyramids are essential** for large motion. Farneback is smooth and dense; neural methods dominate modern benchmarks.

---

## 6.5 Stereo depth

`cv2.StereoSGBM` on Middlebury `Adirondack`:

- **Disparity map:** visually accurate for most of the image.
- **Depth error:** RMSE ≈ 2 cm over 10 m depth range, 3% relative.
- **Failure cases:** occlusions (left edge of objects), textureless regions (sky, walls).

The reconstruction is dense but **noisy** — stereo matching errors appear as speckles in the depth map. Post-processing (median filter, left-right consistency check) removes most.

---

## 6.6 Codec bitrate comparison

10 s 1080p clip at PSNR = 35 dB:

| Codec | Bitrate | BD-rate |
|---|:---:|:---:|
| H.264 | 6 Mbps | (baseline) |
| HEVC | 3.8 Mbps | −38% |
| VP9 | 4.2 Mbps | −30% |
| AV1 | **3.4 Mbps** | **−43%** |

AV1 is the most efficient. H.264 is universal. HEVC is between. Choice depends on support requirements.

---

## 6.7 GOP trade-off

Same clip, GOP sizes 1 (all-I), 12, 30, 60:

| GOP | Bitrate | PSNR | Seek latency |
|:---:|:---:|:---:|:---:|
| 1 | 18 Mbps | 40 dB | 0 |
| 12 | 6 Mbps | 38 dB | 0.5 s |
| 30 | 4.5 Mbps | 37 dB | 1.0 s |
| 60 | 4 Mbps | 36.5 dB | 2.0 s |

**GOP 30 is the sweet spot** for streaming. All-I is for editing (frame-accurate seek). GOP 60 for storage (max compression). Live streaming uses GOP = 1–4 s (2× frame rate) to allow fast joining.

---

## 6.8 Rate control

CBR (`-b:v 1M`) vs. CRF (`-crf 23`):

| Setting | Avg bitrate | Peak bitrate | PSNR |
|---|:---:|:---:|:---:|
| CBR 1 Mbps | 1.0 Mbps | 1.1 Mbps | 32 dB |
| CBR 4 Mbps | 4.0 Mbps | 4.2 Mbps | 38 dB |
| CRF 23 | 3.5 Mbps | 8 Mbps | 40 dB |

**CBR** gives predictable network load — good for streaming with fixed bandwidth. **CRF** gives constant quality — better for storage (uses more bits for complex scenes, fewer for simple). Modern encoders use **capped CRF** (VBR with a bitrate ceiling).

---

## 6.9 Error control

5% packet loss on H.264 stream:

| Setting | Visual quality | PSNR |
|---|---|:---:|
| No concealment | Blocky artifacts, error propagation | 20 dB |
| Frame copy concealment | Artifacts hidden in static areas | 27 dB |
| Frame copy + 15% FEC | Errors fully recovered | **33 dB** |

FEC overhead 15% (K=20, N=23 Reed–Solomon) → ~0.1% remaining error rate at 5% loss. Costs 15% bandwidth; saves 6 dB PSNR — clearly worth it.

---

## 6.10 ABR algorithms

Simulated on HSDPA trace (2 Mbps varying):

| Algorithm | Mean quality | Switches | Stalls |
|---|:---:|:---:|:---:|
| Rate-based (immediate) | 1.8 Mbps | 42 | 3 |
| Buffer-based (BOLA) | 2.1 Mbps | 12 | 0 |
| Hybrid (BBA) | **2.2 Mbps** | **8** | **0** |

**Hybrid algorithms** dominate: buffer-based rules prevent oscillation, rate-based rules allow quick response when bandwidth improves. BBA (Buffer-Based Approach) is a good default.

---

## 6.11 WebRTC latency

Timer visible in camera → screen-to-screen measurement:

- **WebRTC:** ~180 ms end-to-end (browser to browser on same LAN).
- **HLS live (LL-DASH):** ~1.5 s.
- **HLS live (standard):** ~6 s.
- **HLS live (with 30 s buffer):** ~30 s.

**WebRTC wins for interactive.** HLS is fine for broadcast (a viewer doesn't care if the goal is 6 s late). Different latency budgets for different applications.

---

## 6.12 Codec choice

For a hypothetical 4K sports streaming service:

- **Codec:** HEVC (4K HDR support, hardware decoding widespread) or AV1 (better efficiency but decoding less universal). Choose **HEVC for compatibility, AV1 for future-proof**.
- **GOP:** 2–3 s (30–60 frames at 30 fps) — allows fast joining.
- **Bitrate:** 15–25 Mbps for 4K (from HEVC's ~2:1 advantage over H.264).
- **Error control:** LL-DASH with FEC (5%), low-latency (~2–3 s end-to-end), buffer-based ABR.
- **Latency:** 2–3 s acceptable for non-interactive sports. For in-play betting, need WebRTC — very different architecture.