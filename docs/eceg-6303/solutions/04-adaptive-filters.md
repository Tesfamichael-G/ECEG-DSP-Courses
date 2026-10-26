# Solutions to Chapter 4 — Adaptive Filters

[:octicons-arrow-left-24: Back to Chapter 4](../chapters/04-adaptive-filters.md)

---

## 4.1 Steepest descent by hand

$\mathbf R = \text{diag}(4, 1)$, $\mathbf p = (2, 1)$, $\mathbb E[d^2] = 10$.

**Wiener solution.** $\mathbf w_{\text{opt}} = \mathbf R^{-1}\mathbf p = (0.5, 1)$.
$J_{\min} = 10 - \mathbf p^\top\mathbf w_{\text{opt}} = 10 - 1 - 1 = 8$.

**Convergence bound:** $0 < \mu < 2/\lambda_{\max} = 0.5$. Take $\mu = 0.2$.

**Iteration.** $\mathbf w(n+1) = \mathbf w(n) + \mu(\mathbf p - \mathbf R\mathbf w(n))$.

| $n$ | $\mathbf w(n)$ | $J(\mathbf w(n))$ |
|:---:|:---:|:---:|
| 0 | $(0, 0)$ | 10 |
| 1 | $(0.4,\ 0.2)$ | 8.68 |
| 2 | $(0.48,\ 0.36)$ | 8.41 |
| 3 | $(0.496,\ 0.488)$ | 8.26 |

Converges to $(0.5, 1)$ as $n \to \infty$. The slow mode has $\tau_2 = 1/(0.2\cdot 1) = 5$;
the fast mode has $\tau_1 = 1/(0.2\cdot 4) = 1.25$.

---

## 4.2 LMS bound and trade-off

Same $\mathbf R$. $0 < \mu < 2/4 = 0.5$.

| $\mu$ | Misadjustment $M$ | $\tau_{\text{avg}}$ | $M \cdot \tau$ |
|:---:|:---:|:---:|:---:|
| 0.05 | 6.25% | 8 | 0.5 |
| 0.10 | 12.5% | 4 | 0.5 |
| 0.20 | 25% | 2 | 0.5 |

$M \cdot \tau_{\text{avg}} \approx 0.5$ independent of $\mu$ — the classical trade-off.

---

## 4.3 NLMS stability with variable input power

Simulate input power changing 60 dB step (e.g. $\sigma_x$ from $10^{-3}$ to 1).

| Algorithm | Result |
|---|---|
| LMS ($\mu = 0.01$) | diverges when $\sigma_x \to 1$ ($\mu\lambda_{\max} > 2$) |
| NLMS ($\tilde\mu = 0.5$) | converges in both regimes with identical $\tilde\mu$ |

**NLMS output MSE:** ~$0.5/(2-0.5) \approx 33\%$ of $\sigma_d^2$ — **same at both powers**.

---

## 4.4 Channel equalisation with QPSK

Channel $h = (1, 0.8, 0.5, -0.3)$, equaliser $N = 7$, delay $\Delta = 3$.

**Un-equalised eye diagram:** closed due to ISI — cannot distinguish symbols.

**After NLMS training (500 symbols, $\tilde\mu = 0.5$):** eye **open**, residual
ISI below noise floor. Bit error rate (BER) at SNR = 20 dB: $<10^{-3}$ —
acceptable for QPSK.

**Without delay tuning** ($\Delta = 0$): the equaliser cannot see precursor ISI;
the eye remains partially closed. $\Delta \approx (L_h + N - 1)/2 = 5$ is optimal;
$\Delta = 3$ works but not perfectly.

---

## 4.5 Echo cancellation

Room impulse response with $T_{60} = 200$ ms, $N = 1000$ taps (NLMS).

| Phase | ERLE |
|---|:---:|
| Before adaptation | 0 dB |
| After 300 ms of single-talk | 22 dB |
| During double-talk (freeze on) | 22 dB (preserved) |
| During double-talk (freeze off — for comparison) | ~8 dB (diverged) |

**Lesson:** double-talk detection is essential. Without it, the filter
adapts to cancel near-end speech, which is not echo — ERLE collapses.

---

## 4.6 Noise cancellation

Primary ECG + 50 Hz powerline, reference: 50 Hz sinusoid alone.

| $\rho_{v_0 v_1}$ | Theoretical improvement |
|:---:|:---:|
| 0.70 | 3.3 dB |
| 0.90 | 7.0 dB |
| 0.99 | 17.0 dB |

**Simulation** matches theory within 1 dB. At $\rho = 0.99$, the powerline is
suppressed by ~17 dB — enough for the ECG signal to become visible.

---

## 4.7 RLS vs. LMS on coloured input

AR(1) input with $a = 0.95$ ($\chi \approx 400$), 10-tap system.

| Algorithm | Samples to -20 dB error |
|---|---:|
| LMS ($\mu = 0.001$) | ~40,000 |
| NLMS ($\tilde\mu = 0.5$) | ~15,000 |
| RLS ($\lambda = 1$) | **~20** |

**~2000× speedup** for RLS. The coloured input makes LMS extraordinarily
slow; RLS is unaffected.

---

## 4.8 RLS forgetting factor

System whose parameters switch at sample 500.

| $\lambda$ | Adaptation time after switch | Steady-state fluctuation |
|:---:|:---:|:---:|
| 1.0 | does not adapt | very low |
| 0.999 | ~500 samples | low |
| 0.99 | ~150 samples | moderate |
| 0.9 | ~30 samples | high |

**Rule:** $\lambda$ trades steady-state accuracy (large $\lambda$) for tracking
speed (small $\lambda$). Choose $\lambda$ based on how fast the system changes.

---

## 4.9 SW-RLS vs. EW-RLS

Same changing-system scenario as 4.8.

| Metric | EW-RLS ($\lambda = 0.99$) | SW-RLS ($W = 200$) |
|---|:---:|:---:|
| Adaptation time after switch | 150 | **~100** (faster) |
| Steady-state fluctuation | moderate | **moderate** |
| Numerical stability | better | more delicate |

**SW-RLS is faster** because it completely discards samples older than $W$ —
the abrupt change is fully "forgotten" in $W$ steps. EW-RLS's geometric
forgetting leaves a residual "smear" from old data.

---

## 4.10 IIR LMS bias

System: $y[n] = 0.9 y[n-1] + x[n]$. Simplified IIR LMS with 1 pole, 1 zero.

| Iteration | $a_0$ | $b_1$ |
|:---:|:---:|:---:|
| 5000 | 0.97 | 0.87 |
| True | 1.0 | 0.9 |

**Bias:** ~3% on both coefficients. Converges slowly and settles at a
biased fixed point.

**Equation-error correction.** Filter $x[n]$ through $1/A(z)$ before using it
in the update. Converged values become $(0.999, 0.897)$ — bias **removed**
to within 0.3%. The cost is one extra IIR filtering step per update.

---

## 4.11 Fast RLS complexity

| $N$ | Standard RLS (ms) | Fast RLS (ms) | Speedup |
|:---:|:---:|:---:|:---:|
| 100 | 0.5 | 0.05 | 10× |
| 500 | 12 | 0.25 | 48× |
| 2000 | 190 | 1.0 | 190× |

Confirms $O(N^2)$ vs. $O(N)$. **Numerical stability:** standard RLS loses
positive-definiteness of $\mathbf P$ after ~$10^4$ samples; fast RLS can
diverge entirely. Use **QRD-RLS** if stability is critical (but it's $O(N^2)$
again — the price of guaranteed stability).

---

## 4.12 Real audio noise cancellation

Headset recording: primary (speech + ambient noise), reference (ambient noise
only).

| Metric | Value |
|---|:---:|
| Input SNR | -3 dB |
| Output SNR (NLMS, $\tilde\mu = 0.5$, $N = 256$) | +12 dB |
| **Improvement** | **15 dB** |
| Latency | 16 ms |
| Fixed Wiener (training on 2 s) | +13 dB |

NLMS slightly outperforms a fixed Wiener filter trained on a training segment,
because it adapts to the non-stationary noise environment. The 16 ms latency
is acceptable for most headset applications.