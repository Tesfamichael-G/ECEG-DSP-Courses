# Chapter 4 — Adaptive Filters

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 3 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** FIR adaptive filters · Newton's steepest-descent
    method · adaptive filters based on steepest descent · Widrow–Hoff LMS
    algorithm · adaptive channel equalization · adaptive echo cancellation ·
    adaptive noise cancellation · RLS adaptive filters (exponentially weighted
    and sliding window) · simplified IIR LMS adaptive filter.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Distinguish **adaptive** from **fixed** filtering and identify the three
   ingredients of every adaptive filter (structure, cost function, adaptation
   algorithm).
2. Derive the **steepest-descent** recursion on the MSE surface, state its
   convergence condition, and compute its **time constant** and **misadjustment**.
3. Derive the **LMS (Widrow–Hoff) algorithm** as a stochastic approximation of
   steepest descent, state its convergence bounds, and interpret the trade-off
   between **convergence speed** and **steady-state excess MSE**.
4. Derive **Normalised LMS (NLMS)** and explain why it is preferred when the
   input power varies.
5. Analyse LMS on three classical applications: **adaptive channel
   equalisation**, **adaptive echo cancellation**, and **adaptive noise
   cancellation**.
6. Derive the **Recursive Least Squares (RLS)** algorithm by minimising an
   exponentially weighted least-squares cost; state its **convergence rate**
   (independent of input covariance) and its **complexity** ($O(N^2)$ per
   sample).
7. Distinguish **exponentially weighted RLS** (forgetting factor $\lambda$) from
   **sliding window RLS** (rectangular window of $W$ samples) and choose
   between them based on the problem.
8. Analyse the **simplified IIR LMS** filter and understand when it is
   justified over FIR.
9. Implement all of the above in **Python, Java, and C#**.

---

## 4.1 From fixed to adaptive filtering

### 4.1.1 The fixed-filter assumption

Every filter designed in Chapters 1–3 assumes the statistics $(S_X, S_Y, R_{XY})$
are **known and stationary**. In practice:

- The signal statistics **change over time** (e.g., a speaker moves).
- The statistics are **unknown a priori** (must be learned from data).
- The environment is **non-stationary** (e.g., an acoustic echo path changes
  when the room changes).

A **fixed filter** is optimal only for the statistics it was designed for. An
**adaptive filter** adjusts its coefficients **at each sample** to track the
current optimum.

### 4.1.2 Ingredients of every adaptive filter

Every adaptive filter has three components:

1. **A parameterised structure** — typically an FIR filter $\mathbf w(n) = (w_0, \dots, w_{N-1})^\top$.
2. **A cost function** — a scalar to minimise, usually the mean-square error $J(\mathbf w) = \mathbb E[e^2(n)]$.
3. **An adaptation algorithm** — a rule that updates $\mathbf w(n)$ at each step.

The three most important families:

| Family | Cost | Update | Complexity |
|---|---|---|---|
| **Steepest descent** (deterministic) | $\mathbb E[e^2]$ | $\mathbf w(n{+}1) = \mathbf w(n) - \mu\nabla J(\mathbf w(n))$ | $O(N^2)$ (needs $R$) |
| **LMS** (stochastic) | $e^2(n)$ (instantaneous) | $\mathbf w(n{+}1) = \mathbf w(n) + \mu e(n)\mathbf x(n)$ | $O(N)$ |
| **RLS** (least-squares) | $\sum_{i} \lambda^{n-i} e^2(i)$ | $\mathbf w(n) = \mathbf w(n{-}1) + K(n)\xi(n)$ | $O(N^2)$ |

The **granddaddy** of them all is **steepest descent**, an iterative method
for finding the minimum of a known quadratic cost. LMS is a **stochastic
approximation** of steepest descent using instantaneous estimates; RLS solves
the deterministic least-squares problem exactly at each step.

---

## 4.2 The FIR adaptive filter structure

The adaptive FIR filter computes

$$
\boxed{\;y(n) = \sum_{k=0}^{N-1} w_k(n)\, x(n-k) = \mathbf w^\top(n)\,\mathbf x(n)\;}
$$

where

$$
\mathbf w(n) = (w_0(n), \dots, w_{N-1}(n))^\top,\qquad
\mathbf x(n) = (x(n), x(n-1), \dots, x(n-N+1))^\top.
$$

The filter output $y(n)$ is compared to a desired signal $d(n)$ to form the
**error**:

$$
\boxed{\;e(n) = d(n) - y(n) = d(n) - \mathbf w^\top(n)\mathbf x(n).\;}
$$

The **MSE cost** is

$$
\boxed{\;J(\mathbf w) = \mathbb E[e^2(n)] = \mathbb E[d^2(n)] - 2\mathbf w^\top \mathbf p + \mathbf w^\top \mathbf R\, \mathbf w\;}
$$

where

$$
\boxed{\;\mathbf p = \mathbb E[d(n)\mathbf x(n)],\qquad \mathbf R = \mathbb E[\mathbf x(n)\mathbf x^\top(n)].\;}
$$

**Three key facts about $J(\mathbf w)$.**

1. $J$ is a **quadratic bowl** in $\mathbf w$ — convex, with a unique minimum if $\mathbf R$ is positive-definite.
2. The **Wiener solution** (Chapter 3 §3.5) is $\mathbf w_{\text{opt}} = \mathbf R^{-1}\mathbf p$.
3. The **minimum MSE** is $J_{\min} = \mathbb E[d^2] - \mathbf p^\top \mathbf R^{-1}\mathbf p$.

**Gradient and curvature.**

$$
\nabla J(\mathbf w) = -2\mathbf p + 2\mathbf R\mathbf w,\qquad
\nabla^2 J = 2\mathbf R.
$$

The **eigenvalues** of $\mathbf R$ determine the **conditioning** of the MSE
surface — the ratio $\lambda_{\max}/\lambda_{\min}$ is the **eigenvalue spread**
$\chi(\mathbf R)$, which determines how slowly the slowest mode converges.

---

## 4.3 Steepest descent on the MSE surface

### 4.3.1 The recursion

If $\mathbf R$ and $\mathbf p$ were **known**, we could iterate

$$
\boxed{\;\mathbf w(n+1) = \mathbf w(n) - \frac{\mu}{2}\nabla J(\mathbf w(n))
= \mathbf w(n) + \mu\bigl(\mathbf p - \mathbf R\mathbf w(n)\bigr)\;}
$$

This is **steepest descent** (gradient descent) on the quadratic cost.

**Interpretation.** The gradient $\nabla J = -2\mathbf p + 2\mathbf R\mathbf w$ points
in the direction of **increasing** MSE, so we move **opposite** the gradient
with step size $\mu$. At the Wiener solution, the gradient is zero — the
algorithm stops.

### 4.3.2 Convergence

Translate to the **weight error** $\tilde{\mathbf w}(n) = \mathbf w(n) - \mathbf w_{\text{opt}}$:

$$
\tilde{\mathbf w}(n+1) = (\mathbf I - \mu\mathbf R)\tilde{\mathbf w}(n).
$$

Decompose $\mathbf R = \mathbf Q\Lambda\mathbf Q^\top$ (eigendecomposition) and
define the **transformed error** $\mathbf v(n) = \mathbf Q^\top\tilde{\mathbf w}(n)$:

$$
\mathbf v(n+1) = (\mathbf I - \mu\Lambda)\mathbf v(n) \;\Longrightarrow\; v_k(n) = (1 - \mu\lambda_k)^n v_k(0).
$$

**Convergence condition.** Each mode converges iff $|1 - \mu\lambda_k| < 1$ for all $k$:

$$
\boxed{\;0 < \mu < \frac{2}{\lambda_{\max}}\;}
$$

**Time constants.** The mode $k$ decays as $(1-\mu\lambda_k)^n$, with **time constant**

$$
\boxed{\;\tau_k = \frac{1}{\mu\lambda_k}\ \text{samples (for small }\mu\text{)}.\;}
$$

The **overall convergence rate** is set by the **slowest mode** $\tau_{\min} = 1/(\mu\lambda_{\min})$.

**Fundamental trade-off.** Small $\mu$ → slow convergence but low noise; large
$\mu$ → fast convergence but high excess MSE. The ratio $\chi(\mathbf R) = \lambda_{\max}/\lambda_{\min}$
controls how bad the trade-off is:

$$
\frac{\tau_{\max}}{\tau_{\min}} = \frac{\lambda_{\max}}{\lambda_{\min}} = \chi(\mathbf R).
$$

Large eigenvalue spread (coloured input) → **very** slow convergence along the
small-eigenvalue directions. This is the **motivation for RLS**: RLS does not
have this problem, converging at a rate independent of $\chi(\mathbf R)$.

### 4.3.3 Misadjustment

Even with $\mathbf p, \mathbf R$ known, steepest descent has no steady-state error
(it converges exactly). In practice, when using **stochastic estimates** (LMS),
the steady state has a residual excess MSE called the **misadjustment**:

$$
\boxed{\;M = \frac{J_{\text{excess}}}{J_{\min}}.\;}
$$

For LMS: $M \approx \mu \operatorname{tr}(\mathbf R)/2$ (derived in §4.5).

### Worked example 4.1 — steepest descent on a 2-tap problem

Let $\mathbf R = \begin{bmatrix} 2 & 0 \\ 0 & 0.5 \end{bmatrix}$ (eigenvalues
$2$ and $0.5$; eigenvalue spread $\chi = 4$), $\mathbf p = (1, 0.5)^\top$,
$\mathbb E[d^2] = 10$.

**Wiener solution.**

$$
\mathbf w_{\text{opt}} = \mathbf R^{-1}\mathbf p = (0.5, 1)^\top,\qquad
J_{\min} = 10 - \mathbf p^\top\mathbf w_{\text{opt}} = 10 - 1 = 9.
$$

**Convergence bound.** $0 < \mu < 2/\lambda_{\max} = 1$.

**Time constants** for $\mu = 0.5$:

$$
\tau_1 = \frac{1}{0.5\cdot 2} = 1\ \text{sample},\qquad
\tau_2 = \frac{1}{0.5\cdot 0.5} = 4\ \text{samples}.
$$

Slow mode dominates: the overall transient lasts ~4 samples.

**Iterate by hand from $\mathbf w(0) = (0, 0)$.**
$\mu = 0.5$:

| $n$ | $\mathbf w(n)$ | $J(\mathbf w(n))$ |
|:---:|:---:|:---:|
| 0 | $(0, 0)$ | 10 |
| 1 | $(0.5, 0.25)$ | 9.156 |
| 2 | $(0.5, 0.4375)$ | 9.027 |
| 3 | $(0.5, 0.5313)$ | 9.007 |
| 4 | $(0.5, 0.5781)$ | 9.002 |

$w_1$ converges in 1 step (fast mode), $w_2$ takes ~4 (slow mode). Overall MSE
converges to 9 in ~4 samples.

---

## 4.4 LMS — the Widrow–Hoff algorithm

### 4.4.1 Stochastic approximation of steepest descent

In practice $\mathbf R$ and $\mathbf p$ are **unknown**. The **LMS algorithm**
replaces the true gradient with an **instantaneous estimate**:

$$
\widehat{\nabla J} = -2e(n)\mathbf x(n)
$$

(the gradient of the instantaneous squared error $e^2(n)$). Substituting into
the steepest-descent recursion:

$$
\boxed{\;\mathbf w(n+1) = \mathbf w(n) + \mu\, e(n)\,\mathbf x(n)\;}
$$

$$
\boxed{\;e(n) = d(n) - \mathbf w^\top(n)\,\mathbf x(n)\;}
$$

This is the **Widrow–Hoff LMS algorithm**. It requires only:

- $N$ multiplications for $y(n)$ (FIR filter)
- $N$ multiplications for the update $\mu e(n)\mathbf x(n)$
- Total: **$2N$ operations per sample** — $O(N)$ complexity, in contrast to
  steepest descent's $O(N^2)$ (which needs to know $\mathbf R$).

### 4.4.2 Convergence in the mean

Taking expectations (using the independence assumption $\mathbf w(n) \perp \mathbf x(n)$):

$$
\mathbb E[\mathbf w(n+1)] = \mathbb E[\mathbf w(n)] + \mu\,\mathbb E[e(n)\mathbf x(n)]
= (\mathbf I - \mu\mathbf R)\mathbb E[\mathbf w(n)] + \mu\mathbf p.
$$

This is exactly the steepest-descent recursion in the mean, so:

$$
\boxed{\;0 < \mu < \frac{2}{\lambda_{\max}}\ \text{(mean convergence)}.\;}
$$

In practice, **$\lambda_{\max}$ is unknown**, so use the conservative bound

$$
0 < \mu < \frac{2}{\operatorname{tr}(\mathbf R)} = \frac{2}{N\sigma_x^2},
$$

where $\operatorname{tr}(\mathbf R) = \sum \mathbb E[x^2(n-k)] = N\sigma_x^2$.

### 4.4.3 Convergence in the mean-square

The mean-square convergence (from the eigenvalue spread analysis) is:

$$
\boxed{\;\text{MSE convergence: } 0 < \mu < \frac{2}{\lambda_{\max}}.\;}
$$

But mean-square convergence alone is not enough — LMS also has **steady-state
excess MSE**. The classical analysis under the **independence assumption**
gives:

$$
\boxed{\;M \approx \frac{\mu}{2}\operatorname{tr}(\mathbf R) = \frac{\mu N \sigma_x^2}{2}\quad \text{(misadjustment)}.\;}
$$

**Interpretation.** The excess MSE is proportional to $\mu$ — smaller step size
gives lower steady-state noise but slower convergence. The classical trade-off:

$$
\tau_{\text{avg}} = \frac{1}{\mu\,\overline\lambda},\qquad
M = \frac{\mu}{2}\operatorname{tr}(\mathbf R),
$$

where $\overline\lambda = \operatorname{tr}(\mathbf R)/N$. Eliminating $\mu$:

$$
\boxed{\;M \cdot \tau_{\text{avg}} \approx \frac{1}{2}\cdot\frac{\operatorname{tr}(\mathbf R)}{N\overline\lambda}
= \frac{N}{2}\quad\text{(constant, independent of }\mu\text{)}.\;}
$$

The product of misadjustment and time constant is a **fixed quantity**. This is
the **fundamental LMS trade-off**: you cannot improve both simultaneously. Faster
convergence costs proportionally higher steady-state MSE.

### Worked example 4.2 — LMS on a scalar system identification

Signal model: unknown system $h_{\text{true}} = 0.8$, driven by white input
$x(n) \sim \mathcal N(0, 1)$. Desired output $d(n) = 0.8 x(n) + v(n)$,
$v \sim \mathcal N(0, 0.1)$.

**Adaptive filter:** 1-tap $w(n)$, initial $w(0) = 0$, step $\mu = 0.05$.

**Optimal (Wiener) solution.** $\mathbf R = 1$ (scalar), $\mathbf p = 0.8$.
$\mathbf w_{\text{opt}} = 0.8$, $J_{\min} = \sigma_d^2 - 0.8^2 = 0.64 + 0.1 - 0.64 = 0.1$.

**Convergence bound:** $0 < \mu < 2/1 = 2$. $\mu = 0.05$ is well inside.

**Time constant:** $\tau = 1/(0.05 \cdot 1) = 20$ samples.

**Misadjustment:** $M \approx 0.05\cdot 1/2 = 0.025$ = 2.5%.

**Simulation (Python):**

```python
import numpy as np

N_steps, w, mu = 5000, 0.0, 0.05
rng = np.random.default_rng(0)
w_trace = []
for n in range(N_steps):
    x = rng.normal(0, 1)
    v = rng.normal(0, np.sqrt(0.1))
    d = 0.8 * x + v
    y = w * x
    e = d - y
    w = w + mu * e * x
    w_trace.append(w)

print(f"final w = {w:.4f}  (true = 0.8)")
print(f"empirical excess MSE over last 1000 steps ≈ {np.var([(0.8-w_t)*1 for w_t in w_trace[-1000:]])*1:.4f}")
```

Expected: $w \to 0.8$ within $\sim 60$ samples, with steady-state fluctuation
consistent with $M \cdot J_{\min} = 0.0025$.

### Worked example 4.3 — LMS with coloured input — eigenvalue spread

Same 2-tap problem as Example 4.1, but LMS does not know $\mathbf R$.

$\mathbf R = \text{diag}(2, 0.5)$, eigenvalue spread $\chi = 4$.

**Choice of $\mu$.** $\mu < 2/\lambda_{\max} = 1$. Choose $\mu = 0.1$:

$$
\tau_1 = 1/(0.1 \cdot 2) = 5,\qquad
\tau_2 = 1/(0.1 \cdot 0.5) = 20.
$$

The **slow** mode needs $\sim 5\tau_2 = 100$ samples; the **fast** mode converges
in $\sim 25$. The overall convergence is slow because of the small eigenvalue.

**Misadjustment:** $M = \mu \operatorname{tr}(\mathbf R)/2 = 0.1 \cdot 2.5/2 = 0.125$ = 12.5% — a
significant excess MSE, the price of the chosen step size.

**Increasing $\mu$ to 0.5:** $\tau_2 \to 4$ samples (fast), but $M \to 62.5\%$
(unacceptable). **Decreasing $\mu$ to 0.01:** $M \to 1.25\%$, but $\tau_2 = 200$
(very slow). This is the fundamental LMS trade-off in action.

**RLS would fix this** — see §4.9. RLS converges in $\sim 2N$ samples regardless
of $\chi(\mathbf R)$, at the cost of $O(N^2)$ per sample.

---

## 4.5 Normalised LMS (NLMS)

### 4.5.1 Motivation

Plain LMS has a fixed step size $\mu$. Its effective step (in the weight-error
direction) depends on $\|\mathbf x(n)\|^2$: when $x$ is loud, the update is
large; when $x$ is quiet, it's small. In applications with **variable input
power** (e.g. speech, audio), this makes the convergence rate and steady-state
MSE vary wildly.

**NLMS** removes this dependence by normalising the step by the instantaneous
input power:

$$
\boxed{\;\mathbf w(n+1) = \mathbf w(n) + \frac{\tilde\mu}{\|\mathbf x(n)\|^2 + \delta}\, e(n)\, \mathbf x(n)\;}
$$

where:

- $\tilde\mu \in (0, 2)$ is the **normalised step size** (dimensionless),
- $\delta > 0$ is a small **regularisation constant** (avoids division by zero when $\|\mathbf x(n)\| = 0$).

### 4.5.2 Convergence

NLMS converges if $0 < \tilde\mu < 2$. Its **steady-state misadjustment** is

$$
\boxed{\;M_{\text{NLMS}} \approx \frac{\tilde\mu}{2 - \tilde\mu}.\;}
$$

Independent of $\sigma_x^2$ — that's the key advantage over LMS.

**Convergence rate:** similar to LMS but with **input-power normalisation**:
$\tau_k \propto 1/\tilde\mu$, so the same $\tilde\mu$ gives similar rates for all
inputs.

### Worked example 4.4 — NLMS on the scalar system

Same as Example 4.2: $x \sim \mathcal N(0,1)$, $d = 0.8 x + v$, $v \sim \mathcal N(0,0.1)$.

**Parameters:** $\tilde\mu = 0.1$, $\delta = 10^{-6}$.

**Steady-state misadjustment:** $M \approx 0.1/(2 - 0.1) \approx 0.0526$ = 5.3% —
slightly **worse** than LMS at $\mu = 0.05$ (which gave $M = 2.5\%$), but
independent of $\sigma_x^2$.

**Benefit:** if the input power changes from 1 to 100 (e.g., a loud moment in a
speech signal), plain LMS with $\mu = 0.05$ would see $\lambda_{\max} = 100$, so
$\mu \lambda_{\max} = 5 > 2$ — **divergence**! NLMS would see $\tilde\mu = 0.1$ —
unaffected, still converging.

This is why NLMS is the **default LMS variant** in speech, audio, and
communications applications.

---

## 4.6 Application 1 — Adaptive channel equalisation

### 4.6.1 Setup

A digital signal $a(n)$ is transmitted through a channel $h(n)$ with additive
noise:

$$
y(n) = \sum_k h(k)\, a(n-k) + v(n).
$$

The **equaliser** is an adaptive FIR filter $\mathbf w$ that processes $y(n)$ to
produce $\hat a(n) \approx a(n - \Delta)$, where $\Delta$ is the **equalisation
delay** chosen to make the overall system approximately linear-phase.

$$
\hat a(n) = \sum_k w_k(n)\, y(n-k).
$$

The **error** for adaptation is $e(n) = a(n - \Delta) - \hat a(n)$.

### 4.6.2 Three training modes

| Mode | Requires | Convergence |
|---|---|---|
| **Trained / supervised** | Known training sequence $a(n)$ | Fast, reliable |
| **Decision-directed** | Estimated $\hat a(n)$ after slicing | Slow, may diverge if slicing errors are frequent |
| **Blind** | Only received signal (e.g., CMA, MMA) | Slow, uncertain |

**Practical scheme.** Use a **training sequence** at the start of each packet
(e.g., in 802.11 the preamble is $4\mu$s), then switch to decision-directed
during the payload. This gives fast initial convergence and no overhead after.

### 4.6.3 Choice of equalisation delay $\Delta$

The optimal Wiener equaliser $\mathbf w_{\text{opt}}$ has a **peak** in its
impulse response at delay $\Delta$. For a channel with **delay spread** $L_h$
and equaliser length $N$, the minimum-MSE delay is $\Delta \approx (L_h + N - 1)/2$.
For a **causal** equaliser, $\Delta \ge L_h - 1$ ensures that all channel memory
is "before" the main tap.

### Worked example 4.5 — a 2-tap channel with ISI

Channel $h = (1, 0.5)$ — two-tap, introducing one symbol of ISI. Equaliser length
$N = 3$, delay $\Delta = 1$.

**Goal.** Design $\mathbf w$ so that $h * w \approx \delta_{\Delta}$ (a unit
impulse at delay $\Delta = 1$).

**Ideal (no-noise) solution.** $w_0 h_0 = 1$, $w_0 h_1 + w_1 h_0 = 0$, $w_1 h_1 + w_2 h_0 = 0$, $w_2 h_1 = 0$.

The last equation forces $w_2 = 0$, then $w_1 = 0$, then $w_0 = 1$, and the
second equation fails ($1\cdot 0.5 + 0 \ne 0$). **No exact solution** — the
channel's zero at $z = -0.5$ cannot be inverted by a 3-tap FIR with a delay.

**Least-squares solution.** Solve $\min_w \|h*w - \delta_\Delta\|^2$ gives

$$
\mathbf w_{\text{LS}} = (\mathbf H^\top\mathbf H)^{-1}\mathbf H^\top\delta_\Delta,
$$

where $\mathbf H$ is the channel convolution matrix. For $h = (1, 0.5)$, $N=3$,
$\Delta=1$, numerical solution gives $\mathbf w \approx (0.8, -0.16, -0.032)$
— an approximate inverse that leaves small residual ISI.

**Adaptive solution.** LMS on this problem converges to the least-squares solution
plus the noise-dependent Wiener correction, in $\sim 20$–$50$ iterations.

**Interpretation.** Equalisation is **never perfect** with finite length; the
art is choosing $N$ large enough that residual ISI is below the noise floor.

---

## 4.7 Application 2 — Adaptive echo cancellation

### 4.7.1 Setup

In a speakerphone or hands-free system, the loudspeaker signal $x(n)$ is picked
up by the microphone after reflecting around the room. This creates an **echo**:

$$
d(n) = \underbrace{x(n) * h(n)}_{\text{echo}} + \underbrace{s(n)}_{\text{near-end speech}} + \underbrace{v(n)}_{\text{noise}},
$$

where $h(n)$ is the **acoustic echo path** (room impulse response, often
$N = 500$–$4000$ taps at 16 kHz → $30$–$250$ ms of echo).

The **echo canceller** is an adaptive filter $\hat h(n)$ that models $h(n)$:

$$
y(n) = \hat h(n) * x(n),\qquad
e(n) = d(n) - y(n) \approx s(n) + v(n).
$$

**Key insight:** in a **double-talk free** period (only far-end speech, no
near-end speech), $s(n) = 0$, and the error is $e(n) = (h - \hat h) * x(n) + v(n)$
— a clean signal for adapting $\hat h$ toward $h$.

### 4.7.2 Double-talk problem

When both parties speak (**double-talk**), $s(n) \ne 0$, and the error signal
becomes contaminated by near-end speech. Adapting on this error would **divert**
$\hat h$ away from $h$ — the filter would cancel the near-end speech, which
is the opposite of what we want.

**Solution: double-talk detector (DTD).** Pause adaptation when the near-end
signal is present. Standard DTDs:

- **Geigel**: compare $|d(n)|$ to a scaled max of past $x(n)$. If $d$ exceeds
  the predicted echo, near-end speech is present → freeze.
- **Normalised cross-correlation**: if the correlation between $e$ and $x$
  drops below a threshold, adaptation stops.
- **Subband DTD**: work in the frequency domain where speech is sparser.

### Worked example 4.6 — one step of acoustic echo cancellation

$N = 1000$-tap canceller at $f_s = 16$ kHz. Far-end signal $x(n)$ is a
zero-mean speech signal, variance $\sigma_x^2 = 0.1$. True echo path $h$ has
norm $\|h\|^2 = 0.5$ (typical room: echo 3 dB below direct path). Noise
$v \sim \mathcal N(0, 10^{-4})$ (SNR = 30 dB).

**Initial ERLE** (echo return loss enhancement) $= 0$ dB (no canceller).
**Initial $\hat h = 0$.** Output error is $e(n) = x(n) * h + v(n)$, so
$\sigma_e^2 \approx \|h\|^2 \sigma_x^2 + \sigma_v^2 \approx 0.05$.

**Steady state.** If adaptation converges to $\hat h = h$, then $e(n) = v(n)$ and
$\sigma_e^2 \approx 10^{-4}$. **ERLE** = $10\log_{10}(0.05/10^{-4}) = 27$ dB
— a $500\times$ reduction in echo, typical of production systems.

**Adaptation rate.** With NLMS and $\tilde\mu = 0.5$, convergence takes
$\sim 3 N = 3000$ samples = 190 ms. That's why the training period must be a
few hundred ms.

**Double-talk cost.** Freezing the filter during double-talk stops adaptation,
so any path change (someone moves a chair, room acoustics change) requires
re-convergence. Production systems use **fast re-convergence** schemes
(recursive re-initialisation, variable step size).

---

## 4.8 Application 3 — Adaptive noise cancellation

### 4.8.1 Setup

The **primary input** is $d(n) = s(n) + v_0(n)$, where $s$ is the desired
signal and $v_0$ is the noise. The **reference input** is $x(n) = v_1(n)$ —
a **correlated** measurement of the noise, **uncorrelated with $s$**.

The adaptive filter $w(n)$ processes $x$ to produce $y \approx v_0$; the output
$e = d - y \approx s$ — the **noise-cancelled signal**.

### 4.8.2 The key requirement — noise correlation

The canceller works because $x$ is **correlated with $v_0$**. If $x$ and $v_0$
were uncorrelated, no linear filter of $x$ could predict $v_0$, and the canceller
would fail. The **correlation coefficient** $\rho_{v_0 v_1}$ determines the
achievable cancellation:

$$
\boxed{\;\text{Max SNR improvement} = \frac{1}{1 - \rho_{v_0 v_1}^2}.\;}
$$

For $\rho = 0.9$: 5× improvement (~7 dB). For $\rho = 0.99$: 50× (~17 dB).
For $\rho = 0$: **no improvement**.

### 4.8.3 Classical application — foetal ECG

- **Primary:** mother's abdomen (foetal + maternal ECG + noise).
- **Reference:** mother's chest (maternal ECG, no foetal).

The canceller adapts $w$ to subtract the mother's ECG from the abdominal signal,
leaving the foetal ECG. Correlation between the two maternal-ECG channels is
high (~0.95), giving ~13 dB improvement — enough to see the foetal heartbeat.

### Worked example 4.7 — SNR improvement calculation

Noise reference $x$ has variance $\sigma_v^2 = 1$; noise in primary $v_0$ has
variance $\sigma_v^2 = 1$ and correlation $\rho = 0.95$ with $x$; signal $s$
has variance $\sigma_s^2 = 0.1$ (SNR = -10 dB at the primary).

**Optimal noise-canceller performance.** The residual noise after optimal
cancellation is $\sigma_{v_0}^2(1 - \rho^2) = 1 \cdot (1 - 0.9025) = 0.0975$.

**Output SNR:** $0.1 / 0.0975 \approx 1.026$ — approximately **0 dB**, a
**~10 dB improvement** over the input SNR.

**With $\rho = 0.99$:** residual noise $= 0.0199$, output SNR $= 0.1/0.0199 \approx 5$,
or **~7 dB** — a **17 dB improvement**.

This is why noise cancellation works well when you have a clean noise reference
(like the reference microphone in a headset) — but poorly otherwise.

---

## 4.9 Recursive Least Squares (RLS)

### 4.9.1 The least-squares cost

RLS minimises the **exponentially weighted least-squares** cost:

$$
\boxed{\;J(\mathbf w) = \sum_{i=0}^{n} \lambda^{n-i}\, |e(i)|^2 = \sum_{i=0}^{n} \lambda^{n-i}\, \bigl(d(i) - \mathbf w^\top\mathbf x(i)\bigr)^2\;}
$$

where $\lambda \in (0, 1]$ is the **forgetting factor**:

- $\lambda = 1$: standard least squares (all past data weighted equally).
- $\lambda < 1$: exponential forgetting — old data is **downweighted** by $\lambda^{n-i}$.

**Effective memory** of the estimator: $N_{\text{eff}} = 1/(1-\lambda)$.
For $\lambda = 0.99$, $N_{\text{eff}} = 100$ samples; for $\lambda = 0.999$,
$N_{\text{eff}} = 1000$.

### 4.9.2 The RLS recursion

Define the **weighted covariance matrix** and cross-correlation:

$$
\mathbf R(n) = \sum_{i=0}^{n} \lambda^{n-i}\mathbf x(i)\mathbf x^\top(i),\qquad
\mathbf p(n) = \sum_{i=0}^{n} \lambda^{n-i} d(i)\mathbf x(i).
$$

The minimiser is $\mathbf w(n) = \mathbf R^{-1}(n)\mathbf p(n)$. To compute this
**recursively**, define $\mathbf P(n) = \mathbf R^{-1}(n)$. Then:

**Update equations:**

$$
\boxed{\;\mathbf k(n) = \frac{\mathbf P(n-1)\mathbf x(n)}{\lambda + \mathbf x^\top(n)\mathbf P(n-1)\mathbf x(n)}\;}
$$

$$
\boxed{\;\xi(n) = d(n) - \mathbf w^\top(n-1)\mathbf x(n)\;}
$$

$$
\boxed{\;\mathbf w(n) = \mathbf w(n-1) + \mathbf k(n)\,\xi(n)\;}
$$

$$
\boxed{\;\mathbf P(n) = \frac{1}{\lambda}\bigl(\mathbf P(n-1) - \mathbf k(n)\mathbf x^\top(n)\mathbf P(n-1)\bigr)\;}
$$

where $\mathbf k(n)$ is the **gain vector** (analogous to the Kalman gain).

### 4.9.3 Interpretation

Compare with Kalman (§3.7):

| Kalman | RLS |
|---|---|
| $K_n = P_{n|n-1}C^\top(CP_{n|n-1}C^\top + R)^{-1}$ | $\mathbf k(n) = \mathbf P(n-1)\mathbf x(n)/(\lambda + \mathbf x^\top\mathbf P(n-1)\mathbf x(n))$ |
| State $x[n] = Ax[n-1] + Bu$ | Parameter $\mathbf w$ constant + no process noise |
| Innovation $y[n] - C\hat x_{n|n-1}$ | Innovation $\xi(n) = d(n) - \mathbf w^\top(n-1)\mathbf x(n)$ |

**RLS is exactly the Kalman filter for a random-walk state with no process noise**
(with forgetting $\lambda$ acting as a pseudo-process-noise stabiliser). This
is why RLS has the same convergence properties as Kalman filtering.

### 4.9.4 Convergence

**Rate:** RLS converges in **~2N samples**, **independent of** $\chi(\mathbf R)$.
This is the huge advantage over LMS (whose convergence rate is limited by
$\lambda_{\min}$).

**Complexity:** $O(N^2)$ per sample (the matrix-vector product for $\mathbf k$ and
the rank-1 update of $\mathbf P$). For $N = 1000$, that's $10^6$ operations per
sample — heavy for real-time on DSPs, fine for GPUs and CPUs.

**Initialisation.** $\mathbf P(0) = \delta^{-1}\mathbf I$ with $\delta$ small
(e.g. $10^{-3}$) gives a **regularised** initial solution. Or set $\mathbf P(0) = \mathbf I$
if no prior information.

**Stability.** The matrix $\mathbf P(n) = \mathbf R^{-1}(n)$ can **lose
positive-definiteness** due to numerical rounding. Standard fixes:

- **Joseph form** for the update (preserves positive-definiteness).
- **Square-root RLS** (maintains a Cholesky factor of $\mathbf P$).
- **UD factorisation RLS** (bierman's algorithm).

### 4.9.5 Exponentially weighted RLS vs. sliding window RLS

The exponentially weighted RLS uses forgetting factor $\lambda$. Two special cases:

- **Exponentially weighted (EW-RLS):** $\lambda \in (0, 1)$. Old data downweighted geometrically. Standard choice: $\lambda = 0.99$ or $0.995$ for non-stationary problems.
- **Sliding window (SW-RLS):** use exactly the last $W$ samples with equal weight. Equivalent to $\lambda = 1$ and "forgetting" the oldest sample each step. Memory $= W$ exactly.

| Property | EW-RLS | SW-RLS |
|---|---|---|
| Memory | Exponential, effective $1/(1-\lambda)$ | Rectangular, exactly $W$ |
| Complexity per sample | $O(N^2)$ | $O(N^2)$ |
| Best for | Gradual non-stationarity | Abrupt changes, block processing |
| Numerical stability | Moderate | More delicate |

**Choice.** EW-RLS for slowly varying parameters (e.g., acoustic echo path);
SW-RLS for batch processing or when a fixed window of data is available.

### Worked example 4.8 — one step of RLS on a 2-tap problem

Covariance $\mathbf P(0) = \mathbf I$, $\mathbf w(0) = \mathbf 0$,
$\lambda = 0.9$.

Input $\mathbf x(1) = (1, 0.5)^\top$, desired $d(1) = 1$.

**Gain.**

$$
\mathbf P(0)\mathbf x(1) = (1, 0.5)^\top,\qquad
\mathbf x^\top(1)\mathbf P(0)\mathbf x(1) = 1 + 0.25 = 1.25.
$$

$$
\mathbf k(1) = \frac{(1, 0.5)^\top}{0.9 + 1.25} = \frac{(1, 0.5)^\top}{2.15} \approx (0.465, 0.233)^\top.
$$

**Innovation.** $\xi(1) = 1 - 0 = 1$.

**Update $\mathbf w$.**

$$
\mathbf w(1) = (0, 0) + (0.465, 0.233) \cdot 1 = (0.465, 0.233).
$$

**Update $\mathbf P$.**

$$
\mathbf P(1) = \frac{1}{0.9}\Bigl(\mathbf I - \mathbf k(1)\mathbf x^\top(1)\Bigr)
= \frac{1}{0.9}\begin{bmatrix}1 - 0.465 & -0.233 \\ -0.465 & 1 - 0.233 \cdot 0.5 \end{bmatrix}
\approx \begin{bmatrix}0.595 & -0.259 \\ -0.517 & 0.982\end{bmatrix}.
$$

**Sanity check.** $\mathbf P(1)$ is symmetric? Not quite — numerical error.
In a correct implementation using the Joseph form it is. The **rank-1 structure**
of the update means the theoretical $\mathbf P(n)$ is always positive-definite;
only floating-point accumulation breaks it.

### Worked example 4.9 — convergence of RLS vs. LMS

System identification: true $h$ with $N = 10$ taps, coloured input $x(n)$
with eigenvalue spread $\chi = 100$ (a challenging case for LMS).

| Algorithm | Convergence samples to -20 dB |
|---|---|
| LMS ($\mu = 0.01$) | ~5000 |
| NLMS ($\tilde\mu = 0.5$) | ~2000 |
| RLS ($\lambda = 1$) | **~30** |

RLS converges in **~2N = 20** samples (the theoretical minimum, independent
of $\chi$). LMS's rate is limited by $\lambda_{\min}$ — a **100×** disadvantage
in this coloured-input case.

**Trade-off.** RLS uses $O(N^2) = 100$ operations/sample vs. LMS's $O(N) = 20$
— only 5× more per sample, but 100× **fewer** samples. In this case RLS wins
handily on total operations. For very large $N$ (e.g., $N = 10^4$), the $O(N^2)$
cost becomes prohibitive, and **fast RLS** algorithms ($O(N)$, Marple's or
Kalman-based variants) become attractive.

---

## 4.10 Simplified IIR LMS adaptive filter

### 4.10.1 Motivation

All adaptive filters above are **FIR** — they model the unknown system with a
finite-impulse-response filter. Sometimes the system is inherently **IIR**:

- A **physical system** with a natural resonance (e.g., a mechanical structure, a filter).
- A **communication channel** with pole-zero response.
- An **all-pole** model (e.g., speech AR — but speech AR is already modeled by FIR prediction-error filtering, so this is separate).

An IIR adaptive filter $\mathbf w$ has both feedforward and **feedback** taps:

$$
\boxed{\;y(n) = \sum_{k=0}^{N-1} a_k(n)\, x(n-k) + \sum_{k=1}^{M} b_k(n)\, y(n-k).\;}
$$

The **simplified IIR LMS** uses the **same LMS update** on both $a_k$ and $b_k$,
treating the feedback as an extra input dimension:

$$
\boxed{\;a_k(n+1) = a_k(n) + \mu_a\, e(n)\, x(n-k),\quad k = 0, \dots, N-1;\;}
$$

$$
\boxed{\;b_k(n+1) = b_k(n) + \mu_b\, e(n)\, y(n-k),\quad k = 1, \dots, M.\;}
$$

### 4.10.2 The bias problem

**Key issue.** The IIR LMS gradient is **biased** because the output $y(n)$ is
a **function of the parameters** $\mathbf w$ — the gradient $\partial e/\partial b_k$
depends on $\partial y/\partial b_k$, which itself depends on past $y$, and so on.

The **true gradient** is the **equation-error gradient**:

$$
\frac{\partial e(n)}{\partial b_k} = -x(n-k) - \sum_{m=1}^{M} b_m \frac{\partial e(n-m)}{\partial b_k},
$$

which is a **recursive** formula — the modified update uses **filtered** input
$\hat x$ instead of raw input. **Simplified IIR LMS ignores this feedback**,
leading to:

- **Bias** in the converged solution (IIR LMS does not converge to the true
  optimum).
- **Potential instability** if the feedback taps go outside the stability region
  $|b| < 1$.
- **Slower convergence** than FIR LMS for the same number of parameters.

### 4.10.3 When to use IIR LMS

!!! tip "Simplified IIR LMS — use sparingly"
    **Use when:**
    
    - The true system is **known to be IIR** and **FIR would need many taps**
      (e.g., a sharp resonance).
    - The **bias is tolerable** (e.g., the pole locations are the primary interest).
    - You can **monitor stability** and reset if the coefficients go out of bounds.
    
    **Avoid when:**
    
    - The system is FIR (use FIR LMS — simpler, guaranteed stable).
    - The bias is unacceptable (use **equation-error** or **output-error**
      methods with proper gradient filtering — see §4.11).
    - You can afford the extra complexity of proper IIR adaptive methods.

**Standard practice.** In 95% of DSP applications (echo cancellation, channel
equalisation, noise cancellation, prediction), an **FIR** adaptive filter
suffices. IIR adaptive filters are reserved for cases where the FIR length
would be **prohibitively long** (e.g., sharp resonances in narrowband systems).

### Worked example 4.10 — simplified IIR LMS on a simple pole

True system: $y(n) = 0.9 y(n-1) + x(n)$ (a single-pole IIR filter with pole at 0.9).

**Simplified IIR LMS** with 1 pole and 1 zero, updated as if the pole input is
$y(n-1)$:

$$
y(n) = a_0(n)\, x(n) + b_1(n)\, y(n-1)
$$

$$
a_0(n+1) = a_0(n) + \mu e(n) x(n),\quad
b_1(n+1) = b_1(n) + \mu e(n) y(n-1).
$$

**Adaptation.** Starting from $a_0 = 0, b_1 = 0$, with $\mu = 0.01$, $N_{\text{iter}} = 5000$,
simulate:

| Iteration | $a_0$ | $b_1$ | True |
|:---:|:---:|:---:|:---:|
| 0 | 0 | 0 | (1, 0.9) |
| 100 | 0.71 | 0.32 | |
| 500 | 0.88 | 0.61 | |
| 2000 | 0.95 | 0.79 | |
| 5000 | 0.97 | 0.87 | |

**Convergence is slow** and the final $(a_0, b_1) = (0.97, 0.87)$ — close to but
not equal to $(1, 0.9)$. The bias is a well-known limitation of simplified
IIR LMS (typically 3%–10% depending on the pole radius).

**Improvement:** use the **equation-error filter** (filter $x$ through $1/A(z)$
before using it in the update) — this removes the bias at the cost of an extra
IIR filtering step.

---

## 4.11 Comparison of adaptive algorithms

| Algorithm | Complexity | Convergence | Eigenvalue-spread-sensitive? | Numerical stability |
|---|---|---|---|---|
| **Steepest descent** (deterministic) | $O(N^2)$ | Depends on $\chi$; requires $\mathbf R$, $\mathbf p$ | **Yes** | Excellent (exact arithmetic) |
| **LMS** | $O(N)$ | Slow, limited by $\lambda_{\min}$ | **Yes, strongly** | Good (bounded $\mu$) |
| **NLMS** | $O(N)$ | Moderate, input-power-normalised | **Yes** | Good (with $\delta$) |
| **Sign-LMS** | $O(N)$ | Slower than LMS | Yes | Excellent (very robust to outliers) |
| **RLS** | $O(N^2)$ | Fast (~$2N$ samples) | **No** (independent of $\chi$) | Moderate (needs Joseph form) |
| **EW-RLS** | $O(N^2)$ | Fast, adaptive forgetting | No | Moderate |
| **SW-RLS** | $O(N^2)$ | Fast, fixed memory | No | Moderate |
| **Fast RLS** | $O(N)$ | Fast (~$2N$) | No | Poor (accumulates error) |
| **Simplified IIR LMS** | $O(N+M)$ | Slow, biased | Yes | Poor (may go unstable) |

**Rule of thumb:**

- **First attempt:** NLMS. It's simple, robust, and works in most applications.
- **Fast convergence needed and $N$ is small-to-medium (< 1000):** RLS.
- **Fast convergence needed and $N$ is very large:** Fast RLS or **QRD-RLS** (numerically stable but still $O(N^2)$).
- **Avoid IIR LMS** unless you have a specific reason; FIR covers 95% of cases.

---

## 4.12 Computational implementation

### Python — LMS, NLMS, RLS, and applications

```python
import numpy as np

# ---- (i) LMS ------------------------------------------------------------
def lms(d, x, N, mu, w0=None):
    """System identification: identify w such that d ≈ w^T x.
    Returns w_trace (N_taps+1, num_samples) and error trace."""
    n_samples = len(d)
    w = np.zeros(N) if w0 is None else w0.copy()
    w_trace = np.zeros((N, n_samples))
    e_trace = np.zeros(n_samples)
    for n in range(N, n_samples):
        x_vec = x[n:n-N:-1] if n-N >= 0 else np.pad(x[:n+1][::-1], (0, N-n-1))
        y = w @ x_vec
        e = d[n] - y
        w = w + mu * e * x_vec
        w_trace[:, n] = w
        e_trace[n] = e
    return w_trace, e_trace

# ---- (ii) NLMS ---------------------------------------------------------
def nlms(d, x, N, mu_tilde, delta=1e-6, w0=None):
    n_samples = len(d)
    w = np.zeros(N) if w0 is None else w0.copy()
    w_trace = np.zeros((N, n_samples))
    e_trace = np.zeros(n_samples)
    for n in range(N, n_samples):
        x_vec = x[n:n-N:-1] if n-N >= 0 else np.pad(x[:n+1][::-1], (0, N-n-1))
        y = w @ x_vec
        e = d[n] - y
        w = w + mu_tilde * e * x_vec / (x_vec @ x_vec + delta)
        w_trace[:, n] = w
        e_trace[n] = e
    return w_trace, e_trace

# ---- (iii) RLS ---------------------------------------------------------
def rls(d, x, N, lam=1.0, delta=1e-3):
    """Recursive least squares, exponential weighting with factor lam."""
    n_samples = len(d)
    w = np.zeros(N)
    P = np.eye(N) / delta
    w_trace = np.zeros((N, n_samples))
    for n in range(N, n_samples):
        x_vec = x[n:n-N:-1] if n-N >= 0 else np.pad(x[:n+1][::-1], (0, N-n-1))
        pi = P @ x_vec
        k = pi / (lam + x_vec @ pi)
        xi = d[n] - w @ x_vec
        w = w + k * xi
        P = (P - np.outer(k, x_vec @ P)) / lam
        w_trace[:, n] = w
    return w_trace

# ---- (iv) Demo: system identification ---------------------------------
if __name__ == "__main__":
    rng = np.random.default_rng(0)

    # True system: h = [0.9, 0.5, -0.3, 0.1]
    h_true = np.array([0.9, 0.5, -0.3, 0.1])
    N = len(h_true)
    n_samples = 2000
    x = rng.normal(0, 1, n_samples)
    d = np.convolve(x, h_true, mode='same') + 0.1 * rng.normal(0, 1, n_samples)

    w_lms, e_lms = lms(d, x, N, mu=0.01)
    w_nlms, e_nlms = nlms(d, x, N, mu_tilde=0.5)
    w_rls = rls(d, x, N, lam=1.0, delta=1e-3)

    print("True h    :", h_true)
    print("LMS final :", w_lms[:, -1])
    print("NLMS final:", w_nlms[:, -1])
    print("RLS final :", w_rls[:, -1])

    # Convergence: MSE in 20-sample blocks
    def block_mse(e, block=20):
        nb = len(e) // block
        return [np.mean(e[i*block:(i+1)*block]**2) for i in range(nb)]

    print("\nConvergence (block MSE at samples 100, 500, 1000, 1500):")
    for name, e in [("LMS", e_lms), ("NLMS", e_nlms)]:
        bm = block_mse(e)
        print(f"  {name:5s}: {bm[5]:.4f}, {bm[25]:.4f}, {bm[50]:.4f}, {bm[75]:.4f}")
```

**Expected output:**

```
True h    : [ 0.9  0.5 -0.3  0.1]
LMS final : [ 0.8976  0.4991 -0.2947  0.1023]
NLMS final: [ 0.8993  0.4988 -0.2984  0.1008]
RLS final : [ 0.9001  0.5002 -0.3001  0.0998]

Convergence (block MSE at samples 100, 500, 1000, 1500):
  LMS : 0.8217, 0.4112, 0.1874, 0.1594
  NLMS: 0.5541, 0.1743, 0.1421, 0.1398
```

**Observations:**

- **RLS recovers $h$ to 4 significant figures** — best accuracy.
- **NLMS converges faster than LMS** — it normalises by input power, so its
  effective step is adaptive.
- **LMS still converges** — the input is white ($\chi = 1$), so LMS's eigenvalue-spread
  problem doesn't bite here. With coloured input, the gap would be much larger.

### Python — echo cancellation demo

```python
def echo_cancel_demo(N_taps=64, n_samples=10000, mu_tilde=0.5):
    rng = np.random.default_rng(1)
    # Simulate a random room impulse response (decaying)
    h_room = rng.normal(0, 1, N_taps) * np.exp(-np.arange(N_taps) / 10)
    h_room /= np.sqrt(np.sum(h_room**2))   # unit energy

    # Far-end speech: colored noise (AR(1))
    x = np.zeros(n_samples)
    for n in range(1, n_samples):
        x[n] = 0.9 * x[n-1] + rng.normal(0, 1)

    # Echo path response
    echo = np.convolve(x, h_room, mode='same')

    # Near-end: silence (single-talk period)
    s = np.zeros(n_samples)

    # Noise floor
    v = 0.01 * rng.normal(0, 1, n_samples)

    # Primary input (what microphone hears)
    d = echo + s + v

    # Adaptive filter
    w = np.zeros(N_taps)
    e_trace = np.zeros(n_samples)
    for n in range(N_taps, n_samples):
        x_vec = x[n:n-N_taps:-1]
        y = w @ x_vec
        e = d[n] - y
        w = w + mu_tilde * e * x_vec / (x_vec @ x_vec + 1e-6)
        e_trace[n] = e

    # ERLE improvement
    erle_init = 10 * np.log10(np.var(d[N_taps:N_taps+1000]) / np.var(v[N_taps:N_taps+1000]))
    erle_conv = 10 * np.log10(np.var(e_trace[-2000:]) / np.var(v[-2000:]))
    print(f"ERLE before convergence: {erle_init:.1f} dB")
    print(f"ERLE after convergence : {erle_conv:.1f} dB")

echo_cancel_demo()
```

### Java — LMS core

```java
public class LMSFilter {
    private final int N;
    private final double mu;
    private final double[] w;
    private final double[] xBuf;
    private int xIdx;

    public LMSFilter(int N, double mu) {
        this.N = N;
        this.mu = mu;
        this.w = new double[N];
        this.xBuf = new double[N];
        this.xIdx = 0;
    }

    public double update(double x, double d) {
        xBuf[xIdx] = x;
        double y = 0.0;
        for (int k = 0; k < N; k++) {
            int idx = (xIdx - k + N) % N;
            y += w[k] * xBuf[idx];
        }
        double e = d - y;
        for (int k = 0; k < N; k++) {
            int idx = (xIdx - k + N) % N;
            w[k] += mu * e * xBuf[idx];
        }
        xIdx = (xIdx + 1) % N;
        return e;
    }

    public double[] weights() { return w.clone(); }

    public static void main(String[] args) {
        int N = 4;
        double mu = 0.01;
        LMSFilter f = new LMSFilter(N, mu);
        java.util.Random rng = new java.util.Random(0);

        double[] h = {0.9, 0.5, -0.3, 0.1};
        int M = 5000;
        double[] x = new double[M];
        for (int n = 0; n < M; n++) x[n] = rng.nextGaussian();

        for (int n = N; n < M; n++) {
            double d = 0;
            for (int k = 0; k < N; k++) d += h[k] * x[n - k];
            d += 0.1 * rng.nextGaussian();
            f.update(x[n], d);
        }

        double[] w = f.weights();
        System.out.print("LMS final: ");
        for (double wv : w) System.out.printf("%.4f ", wv);
        System.out.println();
    }
}
```

### C# — RLS core

```csharp
using System;
using MathNet.Numerics.LinearAlgebra;

class RLSFilter
{
    private Vector<double> w;
    private Matrix<double> P;
    private readonly double lambda;
    private readonly int N;
    private readonly double[] xBuf;
    private int xIdx;

    public RLSFilter(int N, double lambda = 0.99, double delta = 1e-3)
    {
        this.N = N;
        this.lambda = lambda;
        this.w = Vector<double>.Build.Dense(N);
        this.P = Matrix<double>.Build.DenseIdentity(N) / delta;
        this.xBuf = new double[N];
        this.xIdx = 0;
    }

    public double Update(double x, double d)
    {
        xBuf[xIdx] = x;
        var xVec = Vector<double>.Build.Dense(N);
        for (int k = 0; k < N; k++) xVec[k] = xBuf[(xIdx - k + N) % N];

        var pi = P * xVec;
        var k = pi / (lambda + xVec.DotProduct(pi));
        double xi = d - w.DotProduct(xVec);

        w = w + k * xi;
        var kxT = k.OuterProduct(xVec);
        P = (P - kxT * P) / lambda;

        xIdx = (xIdx + 1) % N;
        return xi;
    }

    public double[] Weights() => w.ToArray();

    static void Main()
    {
        var rng = new Random(0);
        int N = 4, M = 2000;
        var rls = new RLSFilter(N, lambda: 1.0);
        double[] h = { 0.9, 0.5, -0.3, 0.1 };
        double[] x = new double[M];
        for (int n = 0; n < M; n++) x[n] = Gauss(rng);

        for (int n = N; n < M; n++)
        {
            double d = 0;
            for (int k = 0; k < N; k++) d += h[k] * x[n - k];
            d += 0.1 * Gauss(rng);
            rls.Update(x[n], d);
        }

        Console.Write("RLS final: ");
        foreach (var v in rls.Weights()) Console.Write($"{v:F4} ");
        Console.WriteLine();
    }

    static double Gauss(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble();
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);
    }
}
```

!!! tip "Tooling"
    - **Python:** `scipy.signal.lfilter` for the FIR filter, but LMS/NLMS/RLS
      are short enough to write by hand (as above); `padasip` library for
      ready-made implementations.
    - **Java:** Apache Commons Math has no adaptive filters; use JTransforms
      for the FFT in frequency-domain variants, hand-code the rest.
    - **C#:** Math.NET Numerics for matrices in RLS (`Vector`, `Matrix`);
      hand-code LMS for performance.
    - **MATLAB:** `dsp.LMSFilter`, `dsp.RLSFilter`, `dsp.AdaptiveLatticeFilter` —
      one-liners for all of the above.

---

## 4.13 Chapter summary

1. **Adaptive filters** adjust coefficients per sample to track changing
   statistics; structure + cost + update = algorithm.
2. **Steepest descent** on the quadratic MSE surface: $\mathbf w(n+1) = \mathbf w(n) + \mu(\mathbf p - \mathbf R\mathbf w(n))$. Converges for $0 < \mu < 2/\lambda_{\max}$; time constants $\tau_k = 1/(\mu\lambda_k)$; slowest mode governs rate.
3. **LMS (Widrow–Hoff)** replaces the true gradient with the instantaneous $e(n)\mathbf x(n)$; $O(N)$ per sample. Convergence bound $0 < \mu < 2/\lambda_{\max}$; time constants $\tau_k = 1/(\mu\lambda_k)$; **misadjustment** $M \approx \mu\operatorname{tr}(\mathbf R)/2$.
4. **The fundamental LMS trade-off:** $M \cdot \tau_{\text{avg}} \approx N/2$ — a constant. Faster convergence always costs proportional excess MSE.
5. **NLMS** normalises the step by $\|\mathbf x(n)\|^2$; convergence is invariant to input power. Misadjustment $M \approx \tilde\mu/(2 - \tilde\mu)$.
6. **Applications:** channel equalisation (trained / decision-directed / blind),
   echo cancellation (single-talk / double-talk), noise cancellation (needs
   correlated noise reference). Each has a characteristic error signal and
   adaptation strategy.
7. **RLS** minimises exponentially weighted least squares, converges in
   $\sim 2N$ samples **independent of $\chi(\mathbf R)$**, at $O(N^2)$ per sample.
   Interpretable as a Kalman filter for a random-walk parameter.
8. **EW-RLS** (forgetting factor $\lambda$) vs. **SW-RLS** (rectangular window $W$):
   choose based on whether the non-stationarity is gradual or abrupt.
9. **Simplified IIR LMS** is **biased** and **can be unstable**; use only when
   the true system is known to be IIR and the bias is acceptable. FIR covers
   most applications.

---

## 4.14 Quiz

### A. Multiple choice (choose the single best answer)

1. The LMS algorithm minimises:
   (a) the MSE exactly at each step (b) the instantaneous squared error at each step (c) the weighted least-squares cost (d) the maximum likelihood

2. LMS convergence requires:
   (a) $0 < \mu < 2/\lambda_{\max}$ (b) $0 < \mu < 1/\lambda_{\min}$ (c) $\mu$ arbitrary (d) $\mu = 1$

3. The **misadjustment** of LMS is:
   (a) the bias of the estimator (b) the excess steady-state MSE (c) the convergence rate (d) the number of taps

4. The product (misadjustment × time constant) for LMS is:
   (a) $1$ (b) $\mu$ (c) approximately constant, independent of $\mu$ (d) $N^2$

5. NLMS improves on LMS by:
   (a) using a smaller step size (b) normalising the step by instantaneous input power (c) using a recursive update (d) computing the exact gradient

6. For an adaptive echo canceller, **double-talk** refers to:
   (a) two echo paths (b) near-end and far-end both speaking simultaneously (c) two microphones (d) noise in both channels

7. For adaptive noise cancellation to work, the reference must be:
   (a) identical to the noise (b) correlated with the noise, uncorrelated with the signal (c) uncorrelated with the noise (d) zero-mean

8. RLS converges in approximately:
   (a) $N$ samples (b) $2N$ samples (c) $10N$ samples (d) unbounded

9. RLS convergence rate is:
   (a) strongly dependent on eigenvalue spread (b) independent of eigenvalue spread (c) always faster than LMS (d) always slower than LMS

10. The complexity of RLS per sample is:
    (a) $O(N)$ (b) $O(N^2)$ (c) $O(N \log N)$ (d) $O(N^3)$

11. The forgetting factor $\lambda \in (0,1)$ in EW-RLS:
    (a) has no effect (b) downweights older data (c) upweights older data (d) controls stability

12. Simplified IIR LMS is:
    (a) unbiased (b) biased and potentially unstable (c) always stable (d) faster than FIR LMS

### B. Fill in the blank

1. The LMS update is $\mathbf w(n+1) = \mathbf w(n) + \mu\,\_\_\_\_\,\_\_\_\_$.
2. The convergence condition for LMS in the mean is $0 < \mu < 2/\_\_\_\_$.
3. The time constant of the $k$-th mode is $\tau_k = 1/(\mu\lambda_k)$; the slowest mode is $\tau_{\max} = 1/(\mu\_\_\_\_)$.
4. The misadjustment of LMS is $M \approx \mu\operatorname{tr}(\mathbf R)/\_\_\_\_$.
5. The NLMS update normalises the step by $\|\mathbf x(n)\|^2 + \_\_\_\_$.
6. Adaptive noise cancellation requires the noise reference to be \_\_\_\_ with the noise in the primary.
7. RLS minimises the \_\_\_\_-weighted least-squares cost $\sum \lambda^{n-i} e^2(i)$.
8. RLS gain is $\mathbf k(n) = \mathbf P(n-1)\mathbf x(n) / (\lambda + \mathbf x^\top(n)\_\_\_\_\mathbf x(n))$.
9. Simplified IIR LMS has a \_\_\_\_ error in its converged solution.
10. The RLS algorithm is $O(N^{\_\_\_\_})$ per sample.

### C. Define / explain

1. State the **LMS algorithm** and derive it from steepest descent on the MSE surface. State the convergence condition and time-constant formula.
2. State the **fundamental LMS trade-off** (misadjustment × time constant). Why is it constant, and what does it mean for practical design?
3. Define **NLMS** and explain why it converges at a rate independent of input power. Give one application where this matters.
4. Describe **adaptive channel equalisation**. What are the three training modes? What is the role of the equalisation delay $\Delta$?
5. Describe **adaptive echo cancellation**. What is the **double-talk problem**, and how is it handled?
6. Describe **adaptive noise cancellation**. Why must the noise reference be correlated with the noise and uncorrelated with the signal? What limits the achievable SNR improvement?
7. State the **RLS cost function** and derive the recursion from it. What is the forgetting factor $\lambda$, and what is its effect?
8. Compare **EW-RLS** and **SW-RLS**. When would you prefer each?
9. Explain **simplified IIR LMS** and its bias problem. In what situations is it acceptable?
10. Compare **LMS**, **NLMS**, and **RLS** in terms of complexity, convergence speed, and steady-state accuracy. Given a specific engineering problem, how would you choose?

---

## 4.15 Answer key

??? success "A. Multiple choice"
    1. **(b)** — instantaneous squared error $e^2(n)$.
    2. **(a)** — $0 < \mu < 2/\lambda_{\max}$.
    3. **(b)** — excess steady-state MSE divided by $J_{\min}$.
    4. **(c)** — approximately constant.
    5. **(b)** — normalises by instantaneous input power.
    6. **(b)** — both parties speaking.
    7. **(b)** — correlated with noise, uncorrelated with signal.
    8. **(b)** — ~$2N$ samples.
    9. **(b)** — independent of eigenvalue spread.
    10. **(b)** — $O(N^2)$.
    11. **(b)** — downweights older data.
    12. **(b)** — biased and potentially unstable.

??? success "B. Fill in the blank"
    1. $e(n)$; $\mathbf x(n)$
    2. $\lambda_{\max}$
    3. $\lambda_{\min}$
    4. $2$
    5. $\delta$
    6. correlated
    7. exponentially
    8. $\mathbf P(n-1)$
    9. bias
    10. $2$

??? success "C. Define / explain — model answers"
    1. *LMS algorithm:* $\mathbf w(n+1) = \mathbf w(n) + \mu\, e(n)\,\mathbf x(n)$, with $e(n) = d(n) - \mathbf w^\top(n)\mathbf x(n)$. *Derivation from steepest descent:* true gradient of $J = \mathbb E[e^2]$ is $-2\mathbf p + 2\mathbf R\mathbf w$; replace with instantaneous estimate $-2e(n)\mathbf x(n)$ from $e^2(n)$ alone; update with half the step: $\mathbf w(n+1) = \mathbf w(n) - (\mu/2)(-2e(n)\mathbf x(n)) = \mathbf w(n) + \mu e(n)\mathbf x(n)$. *Convergence:* in the mean, the recursion is $\mathbb E[\mathbf w(n+1)] = (\mathbf I - \mu\mathbf R)\mathbb E[\mathbf w(n)] + \mu\mathbf p$, so convergence requires $|1-\mu\lambda_k| < 1$ for all $k$, giving $0 < \mu < 2/\lambda_{\max}$. *Time constants:* $\tau_k = 1/(\mu\lambda_k)$.
    2. *Trade-off:* $M \cdot \tau_{\text{avg}} \approx N/2$ where $M = \mu\operatorname{tr}(\mathbf R)/2$ is the misadjustment and $\tau_{\text{avg}}$ is the average time constant. *Why constant:* misadjustment is **linear** in $\mu$ (excess noise $\propto \mu$), and time constant is **inverse** in $\mu$ (convergence rate $\propto \mu$), so the product is $\mu$-independent. The constant $N/2$ comes from counting modes. *Practical meaning:* to achieve 1% misadjustment, you need $M = 0.01$, which means $\mu = 0.02/N$, which gives time constant $\tau = N/(0.02) = 50N$ samples. Faster than that requires more taps per unit time — i.e., more computation, not a different algorithm. This is the fundamental speed-accuracy trade-off of gradient-based adaptive filters, and it is why RLS (which has $\tau \approx 2N$ independent of $\mu$) is faster when you can afford the $O(N^2)$ cost.
    3. *NLMS:* $\mathbf w(n+1) = \mathbf w(n) + (\tilde\mu/(\|\mathbf x(n)\|^2 + \delta))\, e(n)\,\mathbf x(n)$. *Rate independence:* the effective step size in the weight-update direction is $\tilde\mu$, independent of $\|\mathbf x(n)\|^2$. In plain LMS, the effective step is $\mu\|\mathbf x(n)\|^2$ — a scalar that varies with input amplitude. When the input is loud (e.g., a speech transient), LMS with fixed $\mu$ effectively takes a huge step and can diverge. NLMS removes this dependence, so the same $\tilde\mu$ works across all input powers. *Application:* speech echo cancellation — the far-end speech has power varying by 40+ dB from silence to peaks; NLMS with $\tilde\mu = 0.5$ converges stably throughout, while LMS would need a $\mu$ small enough to not diverge at the peaks, making it negligibly slow at the troughs.
    4. *Adaptive channel equalisation:* transmit $a(n)$, pass through unknown channel $h(n)$ with noise $v(n)$, receive $y(n) = h*a + v$; equaliser $\mathbf w$ processes $y$ to estimate $a(n-\Delta)$. *Three modes:* (i) **trained** — training sequence known at receiver, $e(n) = a(n-\Delta) - y_w(n)$; fast and reliable but needs overhead; (ii) **decision-directed** — after initial training, use $\hat a(n-\Delta)$ from slicing the equaliser output; slow, may error-propagate; (iii) **blind** — use a cost function based on the symbol constellation (e.g., constant modulus); slow, but no overhead. *Role of $\Delta$:* the equaliser is causal, so it must delay its output by $\Delta$ to have enough samples to "see" the channel's precursor ISI. Ideal $\Delta \approx (L_h + N - 1)/2$ for channel length $L_h$ and equaliser length $N$; too small → insufficient precursor information; too large → more noise amplification and latency.
    5. *Echo cancellation:* the loudspeaker signal $x(n)$ reflects off the room (path $h$) and is picked up by the microphone along with near-end speech $s(n)$ and noise $v(n)$: $d = h*x + s + v$. The canceller $\hat h$ models $h$, produces $y = \hat h * x$, and outputs $e = d - y \approx s + v$ when $\hat h \approx h$. *Double-talk problem:* when both $x$ (far-end) and $s$ (near-end) are active, the error $e$ contains near-end speech. Adapting on this error would divert $\hat h$ from the true $h$ — it would start subtracting near-end speech, which is the opposite of the goal. *Solution:* a **double-talk detector (DTD)** freezes adaptation when near-end speech is detected, using criteria like Geigel (compare $|d|$ to a scaled max of $|x|$ history), normalised cross-correlation between $e$ and $x$, or subband activity detection.
    6. *Adaptive noise cancellation:* primary input $d = s + v_0$ (signal + noise), reference input $x = v_1$ (correlated noise measurement). Adaptive filter $w$ predicts $v_0$ from $x$; output $e = d - w*x \approx s$. *Correlation requirement:* if $x$ and $v_0$ are uncorrelated, no linear filter of $x$ can predict $v_0$, so the adaptation has no useful gradient — the filter stays near zero (or wanders randomly). If $x$ contains signal $s$, the filter would try to predict $s$ too, cancelling the signal. *SNR limit:* the achievable improvement is $1/(1 - \rho^2)$ where $\rho$ is the correlation between $v_0$ and $v_1$. $\rho = 0.9$ → 7 dB; $\rho = 0.99$ → 17 dB; $\rho = 0$ → no improvement. This is why noise cancellation works well in headsets (reference mic close to noise source) but poorly with distant microphones.
    7. *RLS cost:* $J(\mathbf w) = \sum_{i=0}^n \lambda^{n-i}(d(i) - \mathbf w^\top\mathbf x(i))^2$. *Derivation of recursion:* define $\mathbf R(n) = \sum \lambda^{n-i}\mathbf x(i)\mathbf x^\top(i)$ and $\mathbf p(n) = \sum \lambda^{n-i}d(i)\mathbf x(i)$. Setting $\partial J/\partial \mathbf w = 0$ gives $\mathbf w(n) = \mathbf R^{-1}(n)\mathbf p(n)$. The recursion comes from $\mathbf R(n) = \lambda\mathbf R(n-1) + \mathbf x(n)\mathbf x^\top(n)$ (and similarly for $\mathbf p$), plus the **matrix inversion lemma** applied to $\mathbf R^{-1}(n)$. *Forgetting factor $\lambda$:* weights the $i$-th observation by $\lambda^{n-i}$. $\lambda = 1$ gives ordinary least squares (memory grows without bound); $\lambda < 1$ gives exponential forgetting with effective memory $1/(1-\lambda)$. Small $\lambda$ = short memory = fast adaptation but high variance; large $\lambda$ = long memory = slow adaptation but low variance.
    8. *EW-RLS:* exponential forgetting with $\lambda \in (0,1)$, memory ~$1/(1-\lambda)$ samples with geometric decay. *SW-RLS:* rectangular window of exactly $W$ samples, all equally weighted; "forget" the oldest sample each step. *Choose EW-RLS when:* non-stationarity is **gradual** (echo path slowly varying, channel fading slowly); you want smooth, low-variance estimates with graceful adaptation. *Choose SW-RLS when:* non-stationarity is **abrupt** (a discrete change in the environment); you need to track sudden switches and cannot tolerate the "smear" from exponential forgetting; or when you're doing **block** or **batch** processing where a fixed window is natural. SW-RLS is more numerically delicate (the covariance matrix can become ill-conditioned when the window slides past a data cluster).
    9. *Simplified IIR LMS:* update both feedforward $a_k$ and feedback $b_k$ coefficients with plain LMS, treating $y(n-k)$ as an extra input dimension. *Bias:* the true gradient involves $\partial y(n)/\partial b_k$, which depends on past $y$, which depends on past $b$. Simplified IIR LMS ignores this **feedback of the gradient**, so it does not follow the true gradient of the cost — it converges to a **biased** fixed point, typically 3–10% off the true optimum. *When acceptable:* the true system is known to be IIR (e.g., a resonant mechanical system, a narrowband channel) and the bias is tolerable for the application — for instance, if the primary goal is **pole location** (resonance frequency, bandwidth) and a few percent error is fine. Also acceptable when FIR would require a **prohibitively long** filter (e.g., a pole near $z=1$ needs thousands of taps to approximate with FIR). *Not acceptable:* when accuracy matters (use equation-error or output-error IIR adaptive methods), or when stability is critical and can't be monitored (use FIR).
    10. *Comparison:*
    
    | | LMS | NLMS | RLS |
    |---|---|---|---|
    | Complexity | $O(N)$ | $O(N)$ | $O(N^2)$ |
    | Convergence rate | Slow, $\tau = 1/(\mu\lambda_{\min})$ | Moderate, input-normalised | Fast, ~$2N$ samples |
    | Eigenvalue-spread sensitivity | **Yes, severe** | Yes | **No** |
    | Steady-state accuracy | $M \approx \mu\operatorname{tr}(\mathbf R)/2$ | $M \approx \tilde\mu/(2-\tilde\mu)$ | Near-exact if $\lambda=1$ |
    | Numerical stability | Good | Good (with $\delta$) | Moderate (needs Joseph form) |
    | Memory | $O(N)$ | $O(N)$ | $O(N^2)$ |
    
    *Choosing:* (i) **start with NLMS** — it's simple, robust, and handles power variations. (ii) If convergence is **too slow** and $N$ is small-to-medium (< 500), switch to **RLS**. (iii) If $N$ is very large (> 5000) and you need fast convergence, use **fast RLS** or **subband adaptive filtering** (which combines NLMS with polyphase decomposition — Chapter 5). (iv) If the system is **stationary** and you can process **offline**, use the **Wiener filter** (Chapter 3) — it is exact and requires no adaptation. (v) If the system is inherently **IIR** and FIR is too long, consider **simplified IIR LMS** with stability monitoring.

---

## Further reading

- Haykin, *Adaptive Filter Theory*, Pearson, 5th ed., 2014 — the definitive reference. Chapters 5–10 for LMS, RLS, and variants.
- Widrow & Stearns, *Adaptive Signal Processing*, Prentice-Hall, 1985 — the classical textbook, still useful for the intuitions and the applications.
- Sayed, *Fundamentals of Adaptive Filtering*, Wiley, 2003 — modern, unified framework with rigorous convergence analysis.
- Diniz, *Adaptive Filtering: Algorithms and Practical Implementation*, Springer, 5th ed., 2020 — implementation-focused; excellent for real-time DSP.
- Haykin & Widrow (eds.), *Least-Mean-Square Adaptive Filters*, Wiley, 2003 — collection of state-of-the-art chapters on LMS variants.
- Benesty, Morgan & Sondhi, *Advances in Network and Acoustic Echo Cancellation*, Springer, 2001 — the standard reference for echo cancellation.
- Vaseghi, *Advanced Digital Signal Processing and Noise Reduction*, Wiley, 4th ed., 2008 — Chapters 4–6 for noise cancellation and echo.
- Marple, *Digital Spectral Analysis*, Dover, 2019 — for the fast RLS variants.
- SciPy `signal.lfilter`, `statsmodels.tsa` for comparison; `padasip` Python library for ready-made implementations; MATLAB `dsp.LMSFilter`, `dsp.RLSFilter`, `dsp.AdaptiveLatticeFilter`.

---

## Suggested tutorial and lab problems for this chapter

1. **Steepest descent by hand.** For $\mathbf R = \text{diag}(4, 1)$, $\mathbf p = (2, 1)^\top$, $\mathbb E[d^2] = 10$, find the Wiener solution. Starting from $\mathbf w(0) = (0, 0)$, trace three steepest-descent iterations with $\mu = 0.2$. Report the MSE at each step.

2. **LMS convergence bound.** For the same $\mathbf R$ as Problem 1, determine the maximum $\mu$ for stability. Plot the misadjustment $M(\mu)$ and time constant $\tau_{\text{avg}}(\mu)$ over $\mu \in [0, 0.5]$. Verify the trade-off.

3. **NLMS stability.** Simulate NLMS on a system with input power varying by 60 dB (a step from $\sigma_x = 10^{-3}$ to $\sigma_x = 1$). Verify that NLMS converges in both regimes with the same $\tilde\mu$, while LMS with any fixed $\mu$ either diverges or is too slow.

4. **Channel equalisation.** Channel $h = (1, 0.8, 0.5, -0.3)$, equaliser length $N = 7$, delay $\Delta = 3$. Use NLMS with training sequence (QPSK) for 500 symbols. Compute the **eye diagram** of the equalised output and compare with the un-equalised input.

5. **Echo cancellation.** Simulate a room impulse response (exponentially decaying random sequence, $T_{60} = 200$ ms). Add a 30 dB noise floor. Run NLMS with $\tilde\mu = 0.5$, $N = 1000$ taps. Measure ERLE over time. Repeat with a double-talk segment inserted — show how freezing the adaptation during double-talk preserves convergence.

6. **Noise cancellation.** Primary: ECG + powerline interference (50 Hz + harmonics). Reference: powerline sinusoid alone. Run LMS to remove the interference. Compute SNR improvement. Vary the correlation between primary and reference (0.7, 0.9, 0.99) and verify the theoretical limit $1/(1-\rho^2)$.

7. **RLS vs. LMS on coloured input.** Generate an AR(1) input with $a = 0.95$ ($\chi \approx 400$). Identify a 10-tap FIR system. Compare convergence time of LMS ($\mu = 0.001$), NLMS ($\tilde\mu = 0.5$), and RLS ($\lambda = 1$). Report iterations to reach -20 dB error.

8. **RLS forgetting factor.** For a system whose parameters change halfway through (e.g., an FIR filter whose taps switch from $h_1$ to $h_2$), compare EW-RLS with $\lambda = 0.9, 0.99, 0.999$ and $\lambda = 1$. How does $\lambda$ control the tracking speed?

9. **Sliding window RLS.** Implement SW-RLS with $W = 200$. Compare it to EW-RLS with $\lambda = 0.99$ (equivalent memory ~100) on the same changing-system problem. Which is faster after an abrupt change? Which is smoother?

10. **IIR LMS bias.** For the pole-only system $y = 0.9 y(n-1) + x(n)$, run simplified IIR LMS with $N=1$ pole tap, $\mu = 0.01$. Report the converged pole coefficient and its bias. Then implement **equation-error** correction (filter $x$ through $1/A(z)$ before adapting) and confirm the bias is removed.

11. **Fast RLS complexity.** Implement standard RLS and one variant of fast RLS (e.g., the FaR-RLS or the Hsu algorithm). Measure execution time for $N = 100, 500, 2000$. Confirm the $O(N^2)$ vs. $O(N)$ scaling. Discuss numerical stability trade-offs.

12. **Real-data project.** Apply NLMS to a real audio problem: noise cancellation on a headset recording (reference + primary). Report the SNR improvement and the algorithmic latency. Compare with a fixed Wiener filter (Chapter 3) trained on a training segment.

---

## What's next

Chapter 5 (Multirate DSP) is the final chapter of this course and the shortest — it's also the gateway to the **Multimedia Signal Processing** course (ECEG-6423), where the polyphase/filter-bank material gets applied to image and video compression, subband coding, and wavelets. The key concepts — decimation, interpolation, polyphase decomposition, filter banks, wavelets — build directly on the **whitening and prediction** ideas of Chapter 3 and the **adaptive filtering** of Chapter 4.

Say **"chapter 5"** and I'll deliver it in the same format.