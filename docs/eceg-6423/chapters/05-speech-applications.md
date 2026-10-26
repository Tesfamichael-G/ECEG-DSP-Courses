# Chapter 5 — Applications of Speech Processing

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 2 tutorial + 1 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** speech recognition — specification of a typical
    application, voice response system · algorithms: spectral estimation,
    **dynamic time warping (DTW)**, **hidden Markov model (HMM)** · speech
    synthesis: text-to-speech, voice over IP · speaker verification and
    identification.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Describe the four major **application families** of speech processing, and
   identify the technology stack each uses.
2. State the **speech-recognition problem** in its three levels (isolated word,
   connected word, continuous speech) and describe the standard pipeline:
   feature extraction → acoustic model → language model → decoder.
3. Derive and apply **dynamic time warping (DTW)** for template matching; state
   its complexity, its constraint variants (Sakoe–Chiba band, Itakura
   parallelogram), and the modern lower-bound acceleration.
4. Derive and apply **hidden Markov models (HMM)**: forward algorithm,
   Viterbi decoding, and Baum–Welch training. State the three fundamental HMM
   problems and their solutions.
5. Compare **HMM-based** and **DTW-based** recognizers. Explain why HMMs
   replaced DTW and why neural networks have now displaced classical HMMs in
   production.
6. Describe the four **speech-synthesis** generations (concatenative, formant,
   HMM-based, neural) and their trade-offs.
7. State the modern **speaker verification** stack: feature extraction, model
   (GMM-UBM / i-vector / x-vector / ECAPA-TDNN), and score normalization.
   Compute EER, FAR, and FRR.
8. Describe the **Voice-over-IP** stack: codec, RTP, jitter buffer, packet-loss
   concealment, and end-to-end latency budget.
9. Implement **DTW**, the **forward algorithm**, and **Viterbi** in Python, Java,
   and C#.

---

## 5.1 Overview of the four application families

Speech processing supports four main applications, each built on the analysis
tools of Chapter 4:

| Application | Input | Output | Core technology |
|---|---|---|---|
| **Recognition** | Audio | Text / intent | HMM, DTW, deep nets |
| **Synthesis** | Text | Audio | Concatenative, neural TTS |
| **Speaker recognition** | Audio | Identity | GMM-UBM, i-vector, x-vector |
| **Coding / transmission** | Audio | Bits | LPC, CELP, Opus |

This chapter covers all four, with **recognition** receiving the most
attention — it is the hardest problem and introduces the two fundamental
algorithms (DTW and HMM) that appear throughout the field.

---

## 5.2 The speech-recognition problem

### 5.2.1 Problem statement

Given an acoustic signal $x(t)$ (typically 5–10 s of speech), find the most
likely word sequence $\hat W = (w_1, \ldots, w_M)$:

$$
\boxed{\;\hat W = \arg\max_W \Pr(W \mid X)\;}
$$

By Bayes' rule:

$$
\Pr(W \mid X) = \frac{\Pr(X \mid W)\Pr(W)}{\Pr(X)}.
$$

Since $\Pr(X)$ is independent of $W$, the argmax reduces to:

$$
\boxed{\;\hat W = \arg\max_W \underbrace{\Pr(X \mid W)}_{\text{acoustic model}} \cdot \underbrace{\Pr(W)}_{\text{language model}}\;}
$$

This **factorization** is the foundation of every ASR system. The acoustic model
scores how well the audio matches each word; the language model scores how
likely each word sequence is; the **decoder** searches for the best combined
score.

### 5.2.2 Three levels of difficulty

| Level | Task | Vocabulary | Constraint |
|---|---|---|---|
| **Isolated word** | One word at a time, with pauses | Small–medium | No context |
| **Connected word** | Fluent speech, small set of words | Small | No language model |
| **Continuous speech** | Natural speech, large vocabulary | 50k–1M+ words | Full language model |

**Production systems** today are **continuous, large-vocabulary** and usually
**speaker-independent**, with vocabulary sizes ranging from 20k (dictation) to
over a million (voice search), and language models trained on billions of
sentences.

### 5.2.3 The standard pipeline

```
Audio  →  Pre-emphasis  →  Framing  →  Windowing  →  Feature
                                                     Extraction
                                                        ↓
                                             (MFCC / PLP / filterbank)
                                                        ↓
                        +-------------------------------------------------+
                        |                                                 |
                        |  Acoustic Model (HMM/DNN): p(X | W)             |
                        |  Language Model (n-gram / RNN): p(W)            |
                        |  Lexicon (dictionary): pronunciation rules      |
                        |  Decoder (Viterbi / WFST): search best W        |
                        |                                                 |
                        +-------------------------------------------------+
                                                        ↓
                                                     Text
```

**Acoustic model.** Scores the likelihood of the acoustic observations given a
phoneme or word sequence. Historically HMM-GMM; now DNN/RNN/Transformer.

**Language model.** Scores the a priori likelihood of word sequences.
Historically n-grams; now recurrent / Transformer LMs.

**Lexicon.** Maps words to phoneme sequences (pronunciation dictionary).

**Decoder.** Searches the combined model for the most likely word sequence.
Dynamic-programming algorithm (Viterbi over the state lattice).

### 5.2.4 Specification of a typical application — voice response

A **voice response system** (e.g., an IVR menu) has much stricter constraints
than dictation:

| Requirement | Value |
|---|---|
| Vocabulary | 10–100 commands (isolated words) |
| Speaker dependency | Speaker-independent preferred |
| Noise tolerance | Must work in a noisy environment |
| Latency | < 300 ms (user expects immediate response) |
| Word accuracy | > 98% (errors are more costly in an IVR than in dictation) |
| Language | One language, trained on domain-specific corpus |

**Design choices for an IVR.**

- **Isolated word / small vocabulary** — simpler acoustic model, no language model needed (or a small grammar).
- **MFCC features + HMM** (or deep net) — the standard architecture.
- **Keyword spotting** — listen for a small set of keywords in continuous speech, no full transcription needed. Simpler than full ASR.
- **Confirm-by-repeat** — ask the user to repeat the choice for ambiguous cases (≥ 85% confidence).

The IVR is a **classic engineering application** — a small, well-defined problem
where accuracy can be pushed to 99%+ with domain-specific training.

---

## 5.3 Dynamic Time Warping (DTW)

### 5.3.1 The need for time warping

The same word, spoken twice, has different **durations** — one utterance may
take 400 ms, another 600 ms. A naive **Euclidean distance** between the
feature sequences is meaningless — the sequences have different lengths and
their features are misaligned in time.

**DTW** solves this by finding an **optimal alignment** — a monotone warping
path that matches the two sequences frame-by-frame, allowing different
speeds.

### 5.3.2 The DTW algorithm

**Problem statement.** Given two feature sequences

$$
X = (x_1, \ldots, x_M),\qquad Y = (y_1, \ldots, y_N),
$$

with $x_i, y_j \in \mathbb{R}^d$ (e.g., MFCC vectors, $d = 13$), find the
**alignment path** $\pi = ((\pi_1^X, \pi_1^Y), \ldots)$ that minimizes the
**accumulated cost**:

$$
\boxed{\;\text{DTW}(X, Y) = \min_{\pi} \sum_{(i,j) \in \pi} d(x_i, y_j)\;}
$$

subject to constraints:

- **Boundary:** path starts at $(1,1)$ and ends at $(M,N)$.
- **Monotonicity:** if $(i,j)$ follows $(i',j')$, then $i \ge i'$ and $j \ge j'$.
- **Step size:** $(i,j)$ follows $(i',j')$ with $(i-i', j-j') \in \{(0,1), (1,0), (1,1)\}$ — path can only move right, up, or diagonally.
- **Locality:** (optional) path stays within a band around the diagonal — the **Sakoe–Chiba band** $|i/M - j/N| \le r$ or the **Itakura parallelogram**.

**The DP recursion.** Compute the accumulated cost matrix:

$$
\boxed{\;C(i,j) = d(x_i, y_j) + \min\begin{cases}C(i-1,j) & \text{(insertion)} \\ C(i,j-1) & \text{(deletion)} \\ C(i-1,j-1) & \text{(match)}\end{cases}\;}
$$

Initial conditions: $C(1,1) = d(x_1,y_1)$; $C(i,1) = C(i-1,1) + d(x_i, y_1)$; $C(1,j) = C(1,j-1) + d(x_1, y_j)$.

**DTW distance:** $C(M,N)$.

**Optimal path:** trace back from $(M,N)$ to $(1,1)$ by following the argmin at each step.

**Complexity.** $O(MN)$ time and $O(MN)$ memory. With the Sakoe–Chiba band of width $r$, it reduces to $O(Mr)$ — a major saving for long sequences.

### Worked example 5.1 — DTW on two simple sequences

$X = (1, 2, 3)$, $Y = (1, 3)$, distance $d(x_i,y_j) = |x_i - y_j|$.

**Distance matrix $D$:**

| | $y=1$ | $y=3$ |
|---|:---:|:---:|
| $x=1$ | 0 | 2 |
| $x=2$ | 1 | 1 |
| $x=3$ | 2 | 0 |

**Accumulated cost matrix $C$:**

$C(1,1) = 0$, $C(1,2) = 0+2=2$, $C(2,1) = 0+1=1$, $C(3,1) = 1+2=3$.

$C(2,2) = 1 + \min(C(1,1), C(1,2), C(2,1)) = 1 + \min(0,2,1) = 1$.

$C(3,2) = 0 + \min(C(2,1), C(2,2), C(3,1)) = 0 + \min(1,1,3) = 1$.

**Final matrix:**

| | $y=1$ | $y=3$ |
|---|:---:|:---:|
| $x=1$ | 0 | 2 |
| $x=2$ | 1 | 1 |
| $x=3$ | 3 | **1** |

**DTW distance:** $C(3,2) = 1$.

**Backtracked path.** Starting at $(3,2)$: predecessor is $C(2,2)$ (min of 1,1,3). From $(2,2)$: predecessor is $C(1,1)$ (min of 0,2,1). From $(1,1)$: start.

**Path:** $(1,1) \to (2,2) \to (3,2)$. So $x_1 \leftrightarrow y_1$, $x_2 \leftrightarrow y_3$, $x_3 \leftrightarrow y_3$.

**Interpretation.** $X$ "stretches" the second sample to match both $y_2$ and $y_3$ — DTW warped the time axis optimally. Total cost 1 = $|1-1| + |2-3| + |3-3|$. ✓

### 5.3.3 Sakoe–Chiba band

The **Sakoe–Chiba band** restricts the path to within $r$ samples of the
diagonal:

$$
\boxed{\;|i \cdot \tfrac{N}{M} - j| \le r,\quad \text{or equivalently } \left|\frac{i}{M} - \frac{j}{N}\right| \le \frac{r}{\min(M,N)}.\;}
$$

**Effects.**

- **Speedup:** $O(MN) \to O(r \min(M,N))$ — for long sequences, a huge speedup.
- **Regularization:** prevents pathological warping (e.g., aligning one frame of $X$ with all of $Y$). For similar-length sequences, this is desirable.
- **Small accuracy loss** in most speech tasks — the path rarely leaves the band anyway.

Typical choices: $r = 0.05 M$–$0.15 M$ (5–15% of the sequence length).

### 5.3.4 LB_Keogh and lower bounds

For template matching over a large database, computing DTW for every candidate
is expensive. **Lower-bound distances** prune candidates without computing the
full DTW:

**LB_Keogh.** Compute the upper and lower envelopes of $Y$ within the
Sakoe–Chiba band:

$$
U_i = \max_{j \in [i-r, i+r]} y_j,\qquad L_i = \min_{j \in [i-r, i+r]} y_j.
$$

Then:

$$
\boxed{\;\text{LB\_Keogh}(X, Y) = \sum_{i=1}^M \begin{cases}(x_i - U_i)^2, & x_i > U_i \\ (x_i - L_i)^2, & x_i < L_i \\ 0, & \text{otherwise}\end{cases}\;}
$$

is a **lower bound** on the DTW distance: $\text{LB\_Keogh}(X,Y) \le \text{DTW}(X,Y)$.

**Use.** For each candidate, compute the cheap lower bound first. If it
exceeds the current best DTW distance, skip the full DTW. Typically prunes
90–99% of candidates — the difference between O(N·M²) and O(N·M) in practice.

The modern standard in time-series matching is the **UCR Suite** (Rakthanmanon
et al., 2012), which combines LB_Keogh with early abandoning, Z-normalization,
and reordering — enables DTW search over billions of subsequences.

### 5.3.5 DTW in isolated-word recognition

**Training.** Store one or more templates per word (recorded examples). Average
across speakers and conditions if available.

**Testing.** For each candidate template $W$, compute $\text{DTW}(X, W)$ on the
MFCC sequence of the test utterance. The word with minimum DTW distance is the
recognized word:

$$
\boxed{\;\hat w = \arg\min_{W \in \text{vocab}} \text{DTW}(X, W).\;}
$$

**Word accuracy.** For a 10-word isolated-digit task with 3 templates per word,
speaker-dependent DTW gives ~99% accuracy; speaker-independent ~95%. This was
the state of the art in the 1970s–80s and still works well for small tasks.

**Limitations.**

- **Template storage** — memory grows with vocabulary and number of templates.
- **No statistical modelling** — no way to represent within-word variability beyond the templates.
- **No language model** — cannot exploit word-sequence statistics.
- **Speaker-dependence** — accuracy drops sharply across speakers.

These limitations motivated the transition to HMMs.

---

## 5.4 Hidden Markov Models (HMM)

### 5.4.1 The HMM framework

An **HMM** models a sequence of observations as generated by an underlying
**Markov chain** whose states are hidden. Formally, an HMM is characterized by:

- **$N$ states** $S = \{1, \ldots, N\}$ — for speech, each phoneme gets 3–5 states (beginning / middle / end).
- **Initial distribution** $\boldsymbol\pi = (\pi_1, \ldots, \pi_N)$ — $P(q_1 = i)$.
- **Transition matrix** $A = [a_{ij}]$ where $a_{ij} = P(q_t = j \mid q_{t-1} = i)$.
- **Emission distribution** $B = \{b_j(\cdot)\}$ — probability of observing a feature vector given the state. For classical HMM-GMM, $b_j(x) = \sum_k c_{jk}\, \mathcal N(x; \mu_{jk}, \Sigma_{jk})$ (Gaussian mixture).

!!! abstract "Three fundamental HMM problems"
    1. **Evaluation.** Given observation $O = (o_1, \ldots, o_T)$ and model $\lambda = (A, B, \boldsymbol\pi)$, compute $P(O \mid \lambda)$. — **Forward algorithm.**
    2. **Decoding.** Given $O$ and $\lambda$, find the most likely hidden state sequence $\hat Q = (\hat q_1, \ldots, \hat q_T)$. — **Viterbi algorithm.**
    3. **Training.** Given one or more observation sequences and an initial model, find the parameters $\lambda = (A, B, \boldsymbol\pi)$ that maximize $P(O \mid \lambda)$. — **Baum–Welch (EM) algorithm.**

### 5.4.2 Evaluation — the forward algorithm

Define the **forward variable**

$$
\boxed{\;\alpha_t(i) = P(o_1, o_2, \ldots, o_t, q_t = i \mid \lambda)\;}
$$

— the probability of the observation prefix and being in state $i$ at time $t$.

**Recursion.**

$$
\boxed{\;\alpha_1(i) = \pi_i\, b_i(o_1)\;}
$$

$$
\boxed{\;\alpha_t(j) = \Bigl[\sum_{i=1}^{N} \alpha_{t-1}(i)\, a_{ij}\Bigr]\, b_j(o_t),\quad t = 2, \ldots, T;\;}
$$

$$
\boxed{\;P(O \mid \lambda) = \sum_{i=1}^{N} \alpha_T(i).\;}
$$

**Cost.** $O(N^2 T)$ — for $N = 3$ states, $T = 400$ frames, $\sim 3600$ multiplications. **Trivial** compared to the $O(N^T)$ naive enumeration.

### Worked example 5.2 — forward algorithm on a 2-state HMM

**Model.**

- States $S = \{1, 2\}$
- $\boldsymbol\pi = (0.6, 0.4)$
- $A = \begin{bmatrix} 0.7 & 0.3 \\ 0.4 & 0.6 \end{bmatrix}$
- Emission:
  - State 1: $b_1(A) = 0.5$, $b_1(B) = 0.5$
  - State 2: $b_2(A) = 0.1$, $b_2(B) = 0.9$
- Observation: $O = A, B, A$

**$\alpha_1$:**

$$
\alpha_1(1) = 0.6 \times 0.5 = 0.30,\quad \alpha_1(2) = 0.4 \times 0.1 = 0.04.
$$

**$\alpha_2$:**

$$
\alpha_2(1) = [0.30 \times 0.7 + 0.04 \times 0.4] \times 0.5 = [0.21 + 0.016] \times 0.5 = 0.113,
$$

$$
\alpha_2(2) = [0.30 \times 0.3 + 0.04 \times 0.6] \times 0.9 = [0.09 + 0.024] \times 0.9 = 0.1026.
$$

**$\alpha_3$:**

$$
\alpha_3(1) = [0.113 \times 0.7 + 0.1026 \times 0.4] \times 0.5 = 0.06007,
$$

$$
\alpha_3(2) = [0.113 \times 0.3 + 0.1026 \times 0.6] \times 0.1 = 0.009546.
$$

**Final.** $P(O \mid \lambda) = 0.06007 + 0.00955 = 0.06962$.

The observation sequence $A, B, A$ has ~7% probability under the model. **The
forward algorithm avoids the $N^T = 8$ full path enumerations** — with $N=3$,
$T=400$, it's the difference between $3^{400}$ (impossible) and 360,000
(instant).

### 5.4.3 Decoding — the Viterbi algorithm

The **Viterbi algorithm** finds the most likely state sequence:

$$
\boxed{\;\hat Q = \arg\max_Q P(Q, O \mid \lambda).\;}
$$

**Recursion.** Define $\delta_t(i) = \max_{q_1,\ldots,q_{t-1}} P(q_1,\ldots,q_{t-1}, q_t = i, o_1, \ldots, o_t \mid \lambda)$.

$$
\boxed{\;\delta_1(i) = \pi_i\, b_i(o_1),\qquad \psi_1(i) = 0;\;}
$$

$$
\boxed{\;\delta_t(j) = \max_{i}\, [\delta_{t-1}(i)\, a_{ij}]\, b_j(o_t),\qquad \psi_t(j) = \arg\max_i\, \delta_{t-1}(i)\, a_{ij}.\;}
$$

**Backtrack.**

$$
q_T^* = \arg\max_i \delta_T(i),\qquad q_t^* = \psi_{t+1}(q_{t+1}^*),\quad t = T-1, \ldots, 1.
$$

**Differences from forward.**

- **Forward** uses **sum** over the previous states → gives **total probability** $P(O \mid \lambda)$.
- **Viterbi** uses **max** over the previous states → gives **best single path** $\max_Q P(Q, O \mid \lambda)$.

In the **log domain**, "sum" becomes **log-sum-exp** and "max" becomes plain max
— the operations stay in the log domain for numerical stability.

### Worked example 5.3 — Viterbi on the same 2-state HMM

Same model and observation $O = A, B, A$.

**$\delta_1$:** $\delta_1(1) = 0.30$, $\delta_1(2) = 0.04$; both $\psi_1 = 0$.

**$\delta_2$:**

$$
\delta_2(1) = \max(0.30 \times 0.7, 0.04 \times 0.4) \times 0.5 = 0.21 \times 0.5 = 0.105,\quad \psi_2(1) = 1.
$$

$$
\delta_2(2) = \max(0.30 \times 0.3, 0.04 \times 0.6) \times 0.9 = 0.09 \times 0.9 = 0.081,\quad \psi_2(2) = 1.
$$

**$\delta_3$:**

$$
\delta_3(1) = \max(0.105 \times 0.7, 0.081 \times 0.4) \times 0.5 = 0.0735 \times 0.5 = 0.03675,\quad \psi_3(1) = 1.
$$

$$
\delta_3(2) = \max(0.105 \times 0.3, 0.081 \times 0.6) \times 0.1 = 0.0486 \times 0.1 = 0.00486,\quad \psi_3(2) = 2.
$$

**Termination.** $q_3^* = \arg\max(0.03675, 0.00486) = 1$.

**Backtrack.** $q_3^* = 1 \Rightarrow q_2^* = \psi_3(1) = 1 \Rightarrow q_1^* = \psi_2(1) = 1$.

**Viterbi path:** $(1, 1, 1)$ — all three states are state 1. The observation
$A, B, A$ is best explained by staying in state 1 throughout, even though
state 2 emits $B$ with high probability — because state 1 started with much
higher probability.

**Comparison with forward.** Forward sums over **all paths** (including
$(1,1,1)$, $(1,1,2)$, $(1,2,1)$, ..., $(2,2,2)$), giving the total likelihood.
Viterbi returns the **single best path** $(1,1,1)$ — the most likely
explanation, not the summed probability.

**Use in ASR.** Viterbi decodes each word's HMM, then the language model
re-ranks the candidate sequences. The combined Viterbi + language model search
is a **single DP** over the lattice of all possible word sequences.

### 5.4.4 Training — Baum–Welch (forward-backward)

The **Baum–Welch algorithm** is EM for HMMs: iteratively re-estimate the
parameters to increase $P(O \mid \lambda)$ until convergence.

**Backward variable.**

$$
\boxed{\;\beta_t(i) = P(o_{t+1}, \ldots, o_T \mid q_t = i, \lambda)\;}
$$

$$
\beta_T(i) = 1,\qquad \beta_t(i) = \sum_{j=1}^{N} a_{ij}\, b_j(o_{t+1})\, \beta_{t+1}(j).
$$

**State posterior.**

$$
\boxed{\;\gamma_t(i) = \frac{\alpha_t(i)\, \beta_t(i)}{P(O \mid \lambda)} = P(q_t = i \mid O, \lambda).\;}
$$

**Transition posterior.**

$$
\boxed{\;\xi_t(i, j) = \frac{\alpha_t(i)\, a_{ij}\, b_j(o_{t+1})\, \beta_{t+1}(j)}{P(O \mid \lambda)} = P(q_t = i, q_{t+1} = j \mid O, \lambda).\;}
$$

**M-step updates.**

$$
\boxed{\;\hat a_{ij} = \frac{\sum_{t=1}^{T-1} \xi_t(i, j)}{\sum_{t=1}^{T-1} \gamma_t(i)},\qquad
\hat b_j(k) = \frac{\sum_{t : o_t = k} \gamma_t(j)}{\sum_{t=1}^{T} \gamma_t(j)},\qquad
\hat\pi_i = \gamma_1(i).\;}
$$

For continuous observations with GMM emissions, $\hat\mu_{jk}$ and
$\hat\Sigma_{jk}$ are updated using the **state- and component-level
posteriors** — an additional responsibility variable per mixture component.

### Worked example 5.4 — one Baum–Welch update

Same 2-state model, observation $O = A, B, A$. We computed $\alpha$ in 5.2.
Now compute $\beta$ (backward) with $\beta_3(1) = \beta_3(2) = 1$.

$$
\beta_2(1) = 0.7 \times 0.5 \times \beta_3(1) + 0.3 \times 0.1 \times \beta_3(2) = 0.35 + 0.03 = 0.38,
$$

$$
\beta_2(2) = 0.4 \times 0.5 \times 1 + 0.6 \times 0.1 \times 1 = 0.20 + 0.06 = 0.26.
$$

$$
\beta_1(1) = 0.7 \times 0.5 \times 0.38 + 0.3 \times 0.9 \times 0.26 = 0.133 + 0.0702 = 0.2032,
$$

$$
\beta_1(2) = 0.4 \times 0.5 \times 0.38 + 0.6 \times 0.9 \times 0.26 = 0.076 + 0.1404 = 0.2164.
$$

**Check** (should equal $P(O)$): $\pi_1 b_1(A)\beta_1(1) + \pi_2 b_2(A)\beta_1(2) = 0.6 \times 0.5 \times 0.2032 + 0.4 \times 0.1 \times 0.2164 = 0.06096 + 0.008656 = 0.0696$. ✓

**State posteriors $\gamma$.**

$$
\gamma_1(1) = \frac{0.30 \times 0.2032}{0.0696} = 0.876,\quad
\gamma_1(2) = \frac{0.04 \times 0.2164}{0.0696} = 0.124.
$$

$$
\gamma_2(1) = \frac{0.113 \times 0.38}{0.0696} = 0.617,\quad
\gamma_2(2) = \frac{0.1026 \times 0.26}{0.0696} = 0.383.
$$

$$
\gamma_3(1) = \frac{0.06007 \times 1}{0.0696} = 0.863,\quad
\gamma_3(2) = \frac{0.009546 \times 1}{0.0696} = 0.137.
$$

**Transition posteriors $\xi$** (only $t=1,2$).

$$
\xi_1(1,1) = \frac{0.30 \times 0.7 \times 0.5 \times 0.38}{0.0696} = \frac{0.0399}{0.0696} = 0.5733,
$$

$$
\xi_2(1,1) = \frac{0.113 \times 0.7 \times 0.5 \times 1}{0.0696} = \frac{0.03955}{0.0696} = 0.568.
$$

**M-step: update $a_{11}$.**

$$
\hat a_{11} = \frac{\xi_1(1,1) + \xi_2(1,1)}{\gamma_1(1) + \gamma_2(1)} = \frac{0.5733 + 0.568}{0.876 + 0.617} = \frac{1.1413}{1.493} = 0.7644.
$$

**Comparison.** Old $a_{11} = 0.70$, new $a_{11} = 0.7644$ — the transition
probability into state 1 has **increased** after observing $A,B,A$, which
favoured staying in state 1. ✓ This is exactly the algorithm's intent: adjust
the transition matrix to better explain the observations.

**Numerical note.** The $\alpha$ and $\beta$ values **shrink exponentially**
with $T$, underflowing in floating point for $T > 100$. **Scaling** is
required (or work in the **log domain** with **log-sum-exp**). Standard
implementations use log-domain arithmetic throughout.

### 5.4.5 HMM for isolated word recognition

**Architecture.** One HMM per word in the vocabulary, each with $N = 3$–$5$
states (roughly one state per phoneme or per subphoneme segment).

**Training.** For each word, train its HMM on multiple examples (Baum–Welch).

**Testing.** For each candidate word's HMM $\lambda_w$, compute
$P(O \mid \lambda_w)$ (forward). The word with highest likelihood:

$$
\boxed{\;\hat w = \arg\max_w P(O \mid \lambda_w).\;}
$$

**Extension to continuous speech.** Connect the word HMMs into a **sentence
HMM** using the language model and a **lexicon** (pronunciation dictionary).
The decoder runs Viterbi over the combined state lattice — the **word
lattice** — to find the best word sequence.

### 5.4.6 Why HMM replaced DTW

| Aspect | DTW | HMM |
|---|---|---|
| Modelling | Template matching | Probabilistic |
| Within-class variability | Only via multiple templates | Via transition and emission distributions |
| Segment duration | Warping path | Implicit in self-transitions |
| Training | Store templates | Baum–Welch uses all data |
| Language model | Not straightforward | Natural via state topology |
| Continuous speech | Hard | Natural via word lattice |

HMM's **probabilistic framework** makes it naturally extendable to larger
vocabularies, continuous speech, and probabilistic language models — that is
why it dominated ASR for 30 years (1980s–2010s).

### 5.4.7 The modern landscape — neural ASR

Since ~2012, deep neural networks have progressively replaced the components
of the HMM pipeline:

| Component | Classical | Neural |
|---|---|---|
| Acoustic model | HMM-GMM | DNN-HMM (2012), then end-to-end (CTC, RNN-T, Transformer) |
| Language model | n-gram | RNN / Transformer LM |
| Decoder | Viterbi over WFST | Beam search over neural decoder |
| Feature extraction | MFCC | Learned (raw waveform or log-mel) |

**Modern systems** (Whisper, Conformer, wav2vec 2.0) are **end-to-end** — a
single neural network maps audio → text, with no explicit HMM or lexicon. Word
error rates on standard benchmarks (LibriSpeech) have dropped from ~10%
(classical HMM-GMM) to ~2% (best neural systems).

**HMM still matters because:**

- Its **mathematics** (forward-backward, Viterbi) underlies the neural tools too — CTC and RNN-T use forward-backward and Viterbi variants.
- Its **concepts** — states, transitions, decoding over a lattice — are reused in neural decoders.
- **Small-vocabulary, low-resource** applications still use HMM-based systems effectively.
- **Understanding HMMs** gives insight into every sequence-labelling problem.

---

## 5.5 Speech synthesis (text-to-speech)

### 5.5.1 The TTS pipeline

**Text-to-speech** converts written text to spoken audio through five stages:

1. **Text normalization** — expand abbreviations, numbers, dates, currency: "Dr." → "doctor", "2024" → "twenty twenty-four" (or "two thousand and twenty-four" depending on context).
2. **Text analysis** — sentence segmentation, part-of-speech tagging, syntactic parsing.
3. **Grapheme-to-phoneme (G2P)** — convert words to phoneme sequences using rules, a dictionary, or a neural network.
4. **Prosody prediction** — pitch contour, duration, stress, pauses, and speaking rate from syntax and semantics.
5. **Waveform generation** — synthesize the acoustic signal from the phoneme + prosody sequence.

Stages 1–4 are **linguistic**; stage 5 is **acoustic**.

### 5.5.2 Four generations of synthesis

**1. Concatenative synthesis (1990s–2000s).**

- **Unit selection** — record a large corpus (10–100 hours) covering many phonemes and contexts.
- **Synthesis** — for each target phoneme, select the best-matching recorded unit from the corpus, minimize join cost + target cost.
- **Quality:** near-human for matched conditions; degrades outside the recorded domain.
- **Storage:** large corpus (hundreds of MB).
- **Example:** AT&T Natural Voices, IBM ViaVoice, modern "unit-selection" TTS.

**2. Formant synthesis (1960s–1990s).**

- **Rule-based** — model the vocal tract as a formant filter, generate glottal pulses, modulate with formant trajectories.
- **Quality:** robotic, intelligible.
- **Storage:** kilobytes.
- **Example:** MITalk, DECtalk, Klatt synthesizer (Stephen Hawking's voice).

**3. HMM-based parametric synthesis (2000s–2010s).**

- **Statistical** — HMM models of the speech frames, trained on a corpus.
- **Synthesis** — generate a speech parameter trajectory from the HMM, then vocode.
- **Quality:** smooth but buzzy (no phase detail).
- **Storage:** small model (a few MB).
- **Example:** HTS (HMM-based TTS), MaryTTS.

**4. Neural synthesis (2016–present).**

- **End-to-end** — a neural network takes text (or phonemes) and outputs audio directly.
- **Architectures:**
  - **WaveNet** (2016) — autoregressive sample-by-sample synthesis, very high quality but slow.
  - **Tacotron 1/2** (2017–2018) — sequence-to-sequence, generates a mel spectrogram, then a vocoder.
  - **FastSpeech 1/2** (2019–2020) — non-autoregressive, much faster.
  - **VITS, NaturalSpeech, Bark, VoiceBox** — current state of the art.
- **Quality:** nearly indistinguishable from human speech in ideal conditions.
- **Storage:** 20–100 MB per voice.
- **Latency:** ~50 ms with optimized inference.
- **Examples:** Google WaveNet, Amazon Polly, Azure Neural TTS, OpenAI Voice Engine, ElevenLabs.

### 5.5.3 Modern neural TTS in detail

The **Tacotron 2 + WaveNet** pipeline is representative:

1. **Text encoder** — character or phoneme sequence → hidden representation (encoder RNN/Transformer).
2. **Attention decoder** — generates a **mel spectrogram** frame by frame (80 mel bins per 12.5 ms).
3. **Vocoder** — mel spectrogram → waveform. **WaveNet** (autoregressive) or **HiFi-GAN** (adversarial) or **WaveGlow** (flow-based).

**End-to-end training.** The whole pipeline is trained on (text, audio) pairs
without separate acoustic and prosodic labels. **Voice cloning** requires only
a few seconds of target voice audio (few-shot learning with speaker embeddings).

### Worked example 5.5 — TTS pipeline walk-through

**Input text.** "The 2024 budget is $1.2 billion."

**Stage 1 — Text normalization.**

$$
\text{"The two-thousand-and-twenty-four budget is one-point-two billion dollars."}
$$

**Stage 2 — Text analysis.** Sentence: declarative. Words: "The" (det), "two-thousand-and-twenty-four" (num), "budget" (noun), "is" (verb), "one-point-two" (num), "billion" (num), "dollars" (noun).

**Stage 3 — G2P.**

- "The" → **DH AH** (or **DH IY** before vowel)
- "two-thousand" → **T UW TH AW Z AH N D**
- "budget" → **B AH JH IH T**
- "is" → **IH Z**
- "one" → **W AH N**
- "point" → **P OY N T**
- "two" → **T UW**
- "billion" → **B IH L Y AH N**
- "dollars" → **D AA L ER Z**

**Stage 4 — Prosody.**

- **Falling intonation** at the end of a declarative sentence.
- **Stress:** "BUDget" (noun stress on first syllable), "bilLION" (numeral stress on second syllable).
- **Duration:** numerals speak slower (more processing time for listeners).
- **Pauses:** none within the sentence; potential short pause after "is".

**Stage 5 — Waveform generation.**

Neural TTS generates a mel spectrogram from the phoneme + prosody sequence,
then converts to a waveform with the vocoder. Total latency for an 8-word
sentence: ~500 ms on GPU, ~2 s on CPU (with WaveNet; faster with HiFi-GAN).

### 5.5.4 Quality evaluation — MOS

**Mean opinion score (MOS)** is the standard subjective metric for TTS quality:

- **5 = Excellent** (indistinguishable from human)
- **4 = Good** (natural, minor artifacts)
- **3 = Fair** (intelligible, noticeable artifacts)
- **2 = Poor** (intelligible with effort)
- **1 = Bad** (unintelligible)

Ratings collected from 20+ listeners per sample, over 50+ sentences, from 3+
different systems. Averaged.

**Typical MOS by generation:**

| Generation | MOS |
|---|---|
| Formant (DECtalk) | 2.5–3.0 |
| Concatenative | 4.0–4.3 |
| HTS (HMM) | 3.5–3.8 |
| Tacotron 2 + WaveNet | **4.5–4.8** |
| FastSpeech 2 + HiFi-GAN | 4.3–4.6 |
| Best neural (2023+) | 4.6–4.9 |

The jump from ~4.0 (concatenative) to ~4.7 (neural) is what made neural TTS
displace every other approach in production. It is now the default for
virtual assistants, audiobook generation, and voice interfaces.

---

## 5.6 Speaker verification and identification

### 5.6.1 Problem definitions

**Speaker identification.** Given a test utterance, determine which of $N$
**enrolled speakers** produced it:

$$
\boxed{\;\hat s = \arg\max_{s \in \{1, \ldots, N\}} P(s \mid X).\;}
$$

- **Closed-set** — the speaker is known to be enrolled; must pick the closest.
- **Open-set** — the speaker may be unknown; must first detect, then classify.

**Speaker verification.** Given a test utterance $X$ and a claimed identity $s$,
decide whether $X$ was spoken by $s$:

$$
\boxed{\;\text{score}(X, s) \gtrless \theta \;\Rightarrow\; \text{accept or reject}.\;}
$$

**Speaker diarization.** "Who spoke when?" — segment a multi-speaker recording
by speaker identity.

### 5.6.2 Features and modeling — the evolution

**Features.** MFCCs (13-dimensional, plus deltas) are the standard. Recent
systems use **larger context** (stacked frames, longer analysis windows) and
**learned features** from deep networks.

**Modelling evolution.**

**1. GMM-UBM (2000s).**

- Train a **Universal Background Model (UBM)** — a large GMM trained on many speakers, capturing speaker-independent distribution of speech.
- For each enrolled speaker, adapt the UBM using **MAP adaptation** on their audio.
- Score: **log-likelihood ratio** between speaker-adapted model and UBM.

**2. i-vector (2011).**

- Represent each utterance as a fixed-length **i-vector** (typically 400–600 dimensions):

$$
\boxed{\;M = m + T\,w,\;}
$$

where $M$ is the supervector (concatenation of GMM means), $m$ is the UBM mean supervector, $T$ is a low-rank matrix, and $w$ is the i-vector.

- **Scoring:** cosine similarity between i-vectors, with **PLDA** (probabilistic linear discriminant analysis) for probabilistic comparison.

**3. x-vector (2017).**

- **DNN-based** — a small network classifies speakers on a training set, and an internal layer (before the classifier) becomes the x-vector.
- **Extraction:** run the test utterance through the network, average the activations of the target layer over frames → x-vector.
- **Scoring:** cosine or PLDA on x-vectors.

**4. ECAPA-TDNN (2020).**

- **State of the art** — 1D-convolutional, **self-attention** pooling, channel-dependent squeeze-excitation.
- EER as low as **0.5%** on VoxCeleb (test set).

### 5.6.3 Evaluation metrics

**False Acceptance Rate (FAR).** Fraction of impostor attempts incorrectly accepted:

$$
\boxed{\;\text{FAR} = \frac{\text{false accepts}}{\text{total impostor trials}}.\;}
$$

**False Rejection Rate (FRR).** Fraction of genuine attempts incorrectly rejected:

$$
\boxed{\;\text{FRR} = \frac{\text{false rejects}}{\text{total genuine trials}}.\;}
$$

**Equal Error Rate (EER).** The threshold at which $\text{FAR} = \text{FRR}$ — a single-number summary of the system. Lower is better.

**Detection Error Tradeoff (DET) curve.** FAR vs. FRR plotted over all thresholds. The **DET curve** shape distinguishes good from poor systems (close to the origin is good).

**Typical EER values:**

| System | EER (on VoxCeleb1-O) |
|---|---|
| GMM-UBM | ~5–10% |
| i-vector + PLDA | ~2–4% |
| x-vector | ~1–2% |
| ECAPA-TDNN | ~0.5–1% |

Modern systems are **more accurate than humans** at speaker verification.

### Worked example 5.6 — EER and FAR/FRR

A speaker-verification system evaluated on 1000 trials: 500 **genuine** (same
speaker) and 500 **impostor** (different speakers). At a candidate threshold
$\theta$:

- Genuine trials: 450 accepted, 50 rejected → FRR = 50/500 = 10%.
- Impostor trials: 20 accepted, 480 rejected → FAR = 20/500 = 4%.

The system is currently at $(\text{FAR}, \text{FRR}) = (4\%, 10\%)$. By raising
the threshold, both FAR and FRR change: fewer false accepts but more false
rejects. The **EER** is the point where they cross.

**Iterating thresholds.** At another threshold $\theta'$:

- Genuine: 470 accepted, 30 rejected → FRR = 6%.
- Impostor: 40 accepted, 460 rejected → FAR = 8%.

At some intermediate threshold $\theta''$: FAR = FRR = 6% → **EER ≈ 6%**.

**Application.** For **access control** (banking), you might choose a threshold
with FAR = 0.1%, FRR = 10% — security-critical. For **forensic analysis** (where
the analyst verifies), FAR = 5%, FRR = 1% — high sensitivity, human confirmation.

### 5.6.4 Text-dependent vs. text-independent

**Text-dependent.** The user speaks a **fixed passphrase** ("my voice is my
password"). Easier because the system knows the target content. Security
benefit: an impostor who doesn't know the phrase is rejected even if their
voice matches.

**Text-independent.** The user speaks **any** content. Harder because content
varies. Better for forensic work (no constraint on suspect's speech).

**Deep-speaker verification** (e.g., for banking) often uses **text-dependent**
with **anti-spoofing** (detecting recorded replay, synthetic speech, and
voice conversion).

### 5.6.5 Anti-spoofing

Modern speaker-verification systems must defend against:

- **Replay attacks** — replay a recording of the enrolled speaker.
- **Synthetic speech** — TTS-generated audio of the target.
- **Voice conversion** — transcode a different speaker's voice to sound like the target.

**Countermeasures:**

- **Liveness detection** — challenge-response (random passphrase), lip-sync checks, pop noise / breath sounds.
- **Deepfake detectors** — trained to distinguish real from synthesized speech (spectral artifacts, phase discontinuities).
- **Random prompts** — the user must read a random sequence (not pre-recordable).

This is an **active research area** — every advance in TTS/voice conversion
raises a new class of attacks.

---

## 5.7 Voice over IP (VoIP)

### 5.7.1 The VoIP stack

**VoIP** transmits speech over IP networks (Internet, enterprise LANs, mobile
data). The full stack has five layers:

| Layer | Protocol | Role |
|---|---|---|
| **Signalling** | SIP, H.323, WebRTC | Call setup, teardown |
| **Media transport** | RTP / UDP | Audio packet delivery |
| **Codec** | G.711, G.729, Opus, AMR | Speech compression |
| **Jitter buffer** | — | Compensate network delay |
| **Quality (QoS)** | DSCP, RSVP | Priority for real-time traffic |

The **codec** is what this course has been building toward — every algorithm
from speech analysis (Chapter 4) is used in the codec.

### 5.7.2 Bandwidth budget

**Narrowband VoIP (G.711):**

- Speech: 64 kbps.
- RTP/UDP/IP headers: 54 bytes/20 ms packet = 21.6 kbps.
- **Total ~85 kbps per call** — 33% overhead just for headers.

**Compressed VoIP (G.729):**

- Speech: 8 kbps.
- RTP/UDP/IP headers: 21.6 kbps.
- **Total ~30 kbps** — headers are 2.7× the payload.

**Header compression (cRTP).** Compresses 54-byte headers to ~4 bytes:

- G.729 + cRTP: 8 + 2 = ~10 kbps per call.

Modern mobile (VoLTE) uses AMR-WB or EVS with header compression — ~15 kbps
per call is typical.

### 5.7.3 Latency budget

**End-to-end latency** for interactive voice:

| Source | Typical delay |
|---|---|
| Encoder | 10–30 ms |
| Packetization | 20 ms |
| Network transit | 50–100 ms |
| Jitter buffer | 20–60 ms |
| Decoder | 10–30 ms |
| **Total (one-way)** | **~130–250 ms** |

**ITU-T G.114** recommends **< 150 ms one-way** for toll-quality; **150–400 ms**
is acceptable; **> 400 ms** makes conversation difficult.

### 5.7.4 Jitter buffer and packet loss

**Jitter.** Packets arrive with **variable delay** — some fast, some slow.
The **jitter buffer** delays each packet by a small amount, smoothing out the
variation.

**Trade-off.** Larger buffer absorbs more jitter but adds more latency. Typical
choice: buffer size ≈ 2–4× the observed jitter.

**Packet loss.** Network packets can be dropped. Three defenses:

1. **FEC (Forward Error Correction).** Send redundant packets — each frame carries bits of the previous frames. Recovers isolated losses at the cost of extra bandwidth.
2. **PLC (Packet Loss Concealment).** Synthesize missing frames from past frames (repetition, LPC extrapolation, pitch-based). No extra bandwidth; quality drop.
3. **Retransmission.** Request retransmission of missing packets. Adds latency — not viable for interactive voice unless the buffer is large.

Modern systems (Opus, EVS) include **built-in PLC and FEC** — they can recover
from 5–10% packet loss with minimal perceptual impact.

### Worked example 5.7 — bandwidth and latency for a mobile call

**4G VoLTE call.**

- Codec: AMR-WB at 12.65 kbps.
- Packetization: 20 ms frames.
- Header compression: RoHC → 4 bytes (32 bits).
- Total bitrate: $12.65 + 32/0.020/1000 = 12.65 + 1.6 = 14.25$ kbps.

**Call with 5% packet loss:**

- **Without FEC/PLC:** 5% of frames are silent or artifact — very noticeable.
- **With Opus's built-in PLC:** most losses are concealed by extrapolation — the listener barely notices.

**Latency budget** (typical VoLTE):

- Codec: 25 ms (AMR-WB frame + algorithmic delay).
- Network: 40 ms (4G RTT one-way).
- Jitter buffer: 40 ms.
- Decoder: 5 ms.
- **Total one-way: ~110 ms** — well below 150 ms.

This is why VoLTE calls have noticeably lower latency than 3G / GSM calls
(~200 ms one-way).

---

## 5.8 Aids to the handicapped — applications in speech

(This is a brief reprise of Chapter 3 §3.7.5, with concrete TTS/ASR technology.)

| Technology | Users | Underlying tech |
|---|---|---|
| **Screen readers** (NVDA, JAWS, VoiceOver) | Blind | TTS on screen content |
| **Speech-to-text** (YouTube, Teams, Google Meet captions) | Deaf / hard-of-hearing | Neural ASR |
| **Hearing aids** | Hard-of-hearing | Speech enhancement, noise reduction, compression |
| **Cochlear implants** | Profoundly deaf | Audio filter bank + neural stimulation |
| **Voice banking** (for ALS, Parkinson's) | Progressive speech loss | Neural TTS on recorded phrases |
| **AAC (Augmentative & Alternative Communication)** | Non-verbal | Symbol-to-TTS |
| **Speech therapy apps** | Children with speech disorders | Formant / pitch visualization, pronunciation scoring |

The **social impact** of these applications is enormous — a blind user's reading
speed drops from ~20 pages/day (Braille) to ~200 pages/day (TTS), and a
hearing-aid user's quality of life improves dramatically with modern
enhancement. Speech processing is one of the highest-leverage areas of DSP
for real-world benefit.

---

## 5.9 Chapter summary

1. **Four applications:** recognition, synthesis, speaker ID, coding. Each reuses the analysis tools from Chapter 4.
2. **Speech recognition** factorizes into acoustic model $\Pr(X \mid W)$ and language model $\Pr(W)$; the decoder searches their product.
3. **DTW** optimally warps the time axis to align two sequences; $O(MN)$ with Sakoe–Chiba band reducing it to $O(Mr)$. LB_Keogh enables fast retrieval over databases.
4. **HMM** is a probabilistic sequence model. Three problems: **evaluation** (forward, $O(N^2 T)$), **decoding** (Viterbi, $O(N^2 T)$), **training** (Baum–Welch EM).
5. **HMM replaced DTW** because it handles within-class variability, extends naturally to large vocabulary and continuous speech, and integrates cleanly with language models.
6. **Neural ASR** (2012–present) uses end-to-end networks that learn features and decode jointly; word error rates have dropped to ~2% on standard benchmarks. HMM concepts still underpin the sequence models (CTC, RNN-T).
7. **TTS** has four generations: formant (robotic), concatenative (natural but rigid), HMM-based (smooth but buzzy), **neural** (nearly human). Neural TTS dominates production.
8. **Speaker verification** uses MFCC → i-vector / x-vector / ECAPA-TDNN → PLDA scoring. EERs have dropped from ~7% (GMM-UBM) to < 1% (ECAPA-TDNN) — more accurate than humans.
9. **VoIP** is the engineering application of speech processing: codec + RTP + jitter buffer + PLC + QoS. Latency must stay < 150 ms; packet loss must be handled gracefully.
10. **Aids to the handicapped** — screen readers, captioning, hearing aids, voice banking — are where speech processing has the largest real-world impact.

---

## 5.10 Quiz

### A. Multiple choice

1. Speech recognition factorizes $\Pr(W \mid X)$ into:
   (a) acoustic model only (b) acoustic × language model (c) forward × backward (d) feature × lexicon

2. **Dynamic time warping** finds the optimal:
   (a) time alignment between two sequences (b) frequency alignment (c) phoneme sequence (d) feature vector

3. DTW's time complexity for sequences of length $M$ and $N$ is:
   (a) $O(M + N)$ (b) $O(M \log N)$ (c) $O(MN)$ (d) $O(M^N)$

4. The **Sakoe–Chiba band** restricts the DTW path to:
   (a) a circle (b) a diagonal band (c) the entire matrix (d) the anti-diagonal

5. The three fundamental HMM problems are:
   (a) evaluation, decoding, training (b) input, output, state (c) filter, predict, smooth (d) fit, predict, retrain

6. The **forward algorithm** computes:
   (a) the best state sequence (b) the probability of the observation sequence (c) the transition matrix (d) the emission distribution

7. The **Viterbi algorithm** finds:
   (a) the total likelihood (b) the best state sequence (c) the parameters (d) the initial distribution

8. The **Baum–Welch** algorithm is an instance of:
   (a) gradient descent (b) expectation-maximization (c) Viterbi (d) Newton's method

9. **Neural TTS** systems like Tacotron 2 typically output:
   (a) waveform directly (b) mel spectrogram then waveform (c) phoneme sequence (d) formant tracks

10. The standard metric for TTS naturalness is:
    (a) PSNR (b) WER (c) MOS (d) EER

11. The **EER** of a speaker-verification system is the threshold where:
    (a) FAR = 1 (b) FRR = 1 (c) FAR = FRR (d) both are zero

12. Modern state-of-the-art speaker verification uses:
    (a) GMM-UBM (b) i-vector (c) x-vector (d) ECAPA-TDNN

### B. Fill in the blank

1. Speech recognition uses the **Bayes factorization** $\hat W = \arg\max \Pr(X \mid W) \cdot \Pr(\_\_\_\_)$.
2. **DTW** finds the optimal **warping path** between two sequences, with complexity $O(MN)$, reduced to $O(Mr)$ by the \_\_\_\_–\_\_\_\_ band.
3. The **LB\_Keogh** lower bound \_\_\_\_ the DTW distance.
4. In an HMM, the **forward algorithm** computes $\Pr(\_\_\_\_ \mid \lambda)$.
5. The **Viterbi** algorithm uses \_\_\_\_ (max or sum?) over previous states.
6. The **Baum–Welch** algorithm uses both **forward** $\alpha$ and **backward** $\beta$ variables to compute \_\_\_\_ posteriors.
7. **Concatenative synthesis** selects recorded \_\_\_\_ from a large corpus.
8. **Neural TTS** typically outputs a **mel spectrogram**, then a \_\_\_\_ converts it to a waveform.
9. **Speaker verification** evaluates via **FAR**, **FRR**, and \_\_\_\_.
10. **VoIP** uses **RTP** over \_\_\_\_ for real-time audio transport.

### C. Define / explain

1. State the **speech-recognition problem** in its Bayes factorization. What are the acoustic model and the language model, and why is the factorization useful?
2. Describe the **DTW algorithm**. Give its recursion, complexity, and one constraint variant (Sakoe–Chiba or Itakura).
3. Define the three **fundamental HMM problems**. For each, state the algorithm that solves it and its complexity.
4. Derive the **forward algorithm** recursion. Show that computing $P(O \mid \lambda)$ naively has cost $O(N^T)$ but the forward algorithm reduces it to $O(N^2 T)$.
5. Compare **Viterbi** and **forward** — how are they similar, how do they differ, and what does each return?
6. State the **Baum–Welch** update equations. Explain why it is an EM algorithm and why it converges to a **local** maximum.
7. Explain why **HMM replaced DTW** as the dominant ASR framework. Give three specific reasons.
8. Describe the **neural TTS pipeline** (Tacotron + WaveNet). How does it differ from concatenative synthesis?
9. Compare **GMM-UBM**, **i-vector**, **x-vector**, and **ECAPA-TDNN** for speaker verification. State the approximate EER of each.
10. State the **VoIP stack**: codec, transport, jitter buffer, QoS. What are the latency budget components for a typical call?

---

## 5.11 Answer key

??? success "A. Multiple choice"
    1. **(b)** — acoustic × language model.
    2. **(a)** — time alignment.
    3. **(c)** — $O(MN)$.
    4. **(b)** — diagonal band.
    5. **(a)** — evaluation, decoding, training.
    6. **(b)** — probability of observation sequence.
    7. **(b)** — best state sequence.
    8. **(b)** — expectation-maximization.
    9. **(b)** — mel spectrogram then waveform.
    10. **(c)** — MOS.
    11. **(c)** — FAR = FRR.
    12. **(d)** — ECAPA-TDNN.

??? success "B. Fill in the blank"
    1. $W$
    2. Sakoe–Chiba
    3. is a lower bound on (≤)
    4. $O \mid \lambda$ (observation sequence given the model)
    5. max
    6. state and transition
    7. units (of speech)
    8. vocoder
    9. EER
    10. UDP

??? success "C. Define / explain — model answers"
    1. *Problem:* $\hat W = \arg\max_W \Pr(W \mid X)$, where $X$ is the audio and $W$ is the word sequence. Bayes factorizes: $\Pr(W \mid X) \propto \Pr(X \mid W)\Pr(W)$. *Acoustic model:* $\Pr(X \mid W)$ — the probability that the audio $X$ was produced from word sequence $W$. Models the mapping from text to acoustics; learned from (audio, transcription) pairs. *Language model:* $\Pr(W)$ — the a priori probability of word sequence $W$. Captures grammar and world knowledge; learned from text corpora. *Why useful:* separating the acoustic and language parts allows independent training on different data — acoustic model on transcribed speech, language model on much larger text-only corpora. It also lets the decoder score candidates according to both factors, so ambiguous acoustics can be resolved by language context ("I scream" vs. "ice cream").
    2. *DTW:* given sequences $X = (x_1, \ldots, x_M)$ and $Y = (y_1, \ldots, y_N)$ with local distance $d$, find the optimal monotonic path $\pi$ aligning them. *Recursion:* $C(i,j) = d(x_i, y_j) + \min(C(i-1,j), C(i,j-1), C(i-1,j-1))$, with DTW distance $C(M,N)$. *Complexity:* $O(MN)$ time and space. *Sakoe–Chiba band:* restrict $|i/M - j/N| \le r$ (or a sample-count version) — reduces to $O(Mr)$ time. Prevents pathological warping (aligning one frame to many), and regularizes the solution to stay near the diagonal. Typical $r = 0.05$–$0.15$ of the sequence length.
    3. **Evaluation:** given $\lambda$ and $O$, compute $P(O \mid \lambda)$. Solved by the **forward algorithm** in $O(N^2 T)$. **Decoding:** given $\lambda$ and $O$, find $\hat Q = \arg\max_Q P(Q \mid O, \lambda)$. Solved by the **Viterbi algorithm** in $O(N^2 T)$. **Training:** given $O$, find $\lambda$ maximizing $P(O \mid \lambda)$. Solved by **Baum–Welch (EM)** in $O(N^2 T)$ per iteration (with several iterations to convergence).
    4. *Recursion:* $\alpha_1(i) = \pi_i b_i(o_1)$; $\alpha_t(j) = [\sum_{i} \alpha_{t-1}(i) a_{ij}] b_j(o_t)$; $P(O \mid \lambda) = \sum_i \alpha_T(i)$. *Complexity argument:* the naive computation of $P(O \mid \lambda)$ enumerates all $N^T$ state sequences, each requiring $T$ multiplications — total $O(N^T)$. The forward algorithm recognizes that $\alpha_t(j)$ is the **same value** regardless of which specific history led to state $j$ — so each $\alpha_t(j)$ is computed once. There are $N$ states × $T$ time steps = $NT$ alphas, each summing over $N$ predecessors → $O(N^2 T)$ total. The key insight: **dynamic programming** reuses subproblems to collapse the exponential tree into a polynomial lattice.
    5. *Similarity:* both are dynamic-programming recursions over the same $(state, time)$ lattice, with identical complexity $O(N^2 T)$. *Difference:* forward uses **sum** over predecessors (computes total probability mass reaching state $j$); Viterbi uses **max** over predecessors (finds the single best path ending at $j$). *Output:* forward returns $P(O \mid \lambda)$, the total probability of the observation under any state sequence. Viterbi returns $\hat Q = \arg\max_Q P(Q, O \mid \lambda)$ (and the corresponding max probability) — the single most likely state sequence. *Use:* forward for evaluation and training; Viterbi for decoding (in ASR, deciding which phoneme was spoken).
    6. *Baum–Welch M-step:* $\hat a_{ij} = \sum_{t} \xi_t(i,j) / \sum_t \gamma_t(i)$; $\hat b_j(k) = \sum_{t:o_t=k} \gamma_t(j) / \sum_t \gamma_t(j)$; $\hat\pi_i = \gamma_1(i)$; and for GMM emissions, $\hat\mu_{jk}$, $\hat\Sigma_{jk}$ are updated with component-level posteriors. *Why EM:* the E-step computes the **posterior state occupancies** given the current parameters (the expected counts of transitions and emissions). The M-step sets the parameters to the **normalized expectations** — i.e., the maximum-likelihood estimates given the posterior. This is exactly the EM template: iteratively refine the posterior and the parameters. *Why local max:* the E-step increases a lower bound on the log-likelihood (the ELBO); the M-step maximizes that lower bound. Both steps increase $P(O \mid \lambda)$, and the algorithm converges to a **stationary point** — typically a local maximum but not necessarily the global one. Initialization matters; **multiple restarts** and **good priors** help.
    7. *HMM replaced DTW for three reasons:* (i) **Probabilistic modelling** — HMM represents within-class variability as distributions (Gaussian emissions, transition probabilities) rather than as a fixed set of templates; requires no template selection. (ii) **Statistical training** — Baum–Welch uses **all training data** to jointly optimize all parameters, whereas DTW only stores templates (and each template captures one recording's characteristics). (iii) **Composability** — HMMs for individual words/phonemes connect via the lexicon and language model into a **sentence HMM**, enabling continuous large-vocabulary speech recognition. DTW-based systems don't compose cleanly into a language model framework. The historical transition (1980s) was driven by these three advantages, at a cost of greater complexity and training time.
    8. *Neural TTS pipeline:* (i) Text/phoneme encoder — embed the input into a hidden representation. (ii) Attention decoder — autoregressively (Tacotron 2) or in parallel (FastSpeech 2) generate a **mel spectrogram**. (iii) Neural vocoder — convert mel spectrogram to waveform (WaveNet, HiFi-GAN, WaveGlow). *Difference from concatenative:* concatenative selects and splices recorded units from a large corpus — requires storing hundreds of MB, only speaks what was recorded. Neural TTS is a **parametric model**: it learns a continuous mapping from text to waveform, is **storage-efficient** (tens of MB for the model, no corpus), can speak **arbitrary text**, and is **smoother** because it interpolates in the model's latent space rather than splicing discrete recordings.
    9. *Speaker-verification systems:* (i) **GMM-UBM** (2000s) — large GMM trained on many speakers (UBM), MAP-adapted per enrolled speaker. EER ~5–10% on VoxCeleb. (ii) **i-vector** (2011) — fixed-length low-rank representation of the GMM supervector. EER ~2–4%. (iii) **x-vector** (2017) — DNN trained for speaker classification, internal layer used as the embedding. EER ~1–2%. (iv) **ECAPA-TDNN** (2020) — 1D-conv + self-attention pooling. EER ~0.5–1%. *Trend:* increasingly deep, end-to-end, and data-hungry. Each generation roughly **halves the EER** — from ~7% to < 1% in one decade.
    10. *VoIP stack:* **codec** (G.711 at 64 kbps, G.729 at 8 kbps, AMR-WB at 12.65 kbps, Opus adjustable) → **RTP/UDP/IP transport** → **jitter buffer** (absorbs network delay variation, adds 20–60 ms) → **QoS** (DSCP priority, RSVP reservations). *Latency budget:* codec 10–30 ms, packetization 20 ms, network 50–100 ms, jitter buffer 20–60 ms, decoder 10–30 ms — total **~110–250 ms one-way**. ITU-T G.114 recommends < 150 ms for toll quality; > 400 ms makes conversation difficult. *Packet-loss handling:* FEC (redundant packets) + PLC (synthesize missing frames from past).

---

## Further reading

**Recognition**

- Rabiner & Juang, *Fundamentals of Speech Recognition*, Prentice Hall, 1993 — **the HMM reference**; Chapters 5–7 cover the material in §5.4 thoroughly.
- Huang, Acero & Hon, *Spoken Language Processing*, Prentice Hall, 2001 — comprehensive modern reference (pre-neural).
- Jurafsky & Martin, *Speech and Language Processing*, 3rd ed., 2024 — free online; the ASR chapters are excellent.
- Jelinek, *Statistical Methods for Speech Recognition*, MIT Press, 1997 — the classical statistical approach.
- Yu & Deng, *Automatic Speech Recognition: A Deep Learning Approach*, Springer, 2015 — the DNN transition.
- Graves, *Supervised Sequence Labelling with Recurrent Neural Networks*, Springer, 2012 — CTC and RNN for ASR.
- Radford et al., "Robust Speech Recognition via Large-Scale Weak Supervision", OpenAI (2022) — the Whisper paper.

**DTW and time-series matching**

- Sakoe & Chiba, "Dynamic Programming Algorithm Optimization for Spoken Word Recognition", *IEEE Trans. ASSP*, 1978 — the classic reference.
- Ratanamahatana & Keogh, "Making Time-Series Classification More Accurate Using Learned Constraints", *SDM* 2004 — DTW in data mining.
- Rakthanmanon et al., "Searching and Mining Trillions of Time Series Subsequences under Dynamic Time Warping", *KDD* 2012 — the UCR Suite.
- Müller, *Information Retrieval for Music and Motion*, Springer, 2007 — DTW with extensions (subsequence DTW, etc.).

**Synthesis**

- Taylor, *Text-to-Speech Synthesis*, Cambridge University Press, 2009 — the classical reference.
- van den Oord et al., "WaveNet: A Generative Model for Raw Audio", *arXiv* 2016.
- Wang et al., "Tacotron: Towards End-to-End Speech Synthesis", *Interspeech* 2017.
- Shen et al., "Natural TTS Synthesis by Conditioning WaveNet on Mel Spectrogram Predictions", *ICASSP* 2018 (Tacotron 2).
- Ren et al., "FastSpeech 2: Fast and High-Quality End-to-End Text to Speech", *ICLR* 2021.
- Kong, Kim & Bae, "HiFi-GAN: Generative Adversarial Networks for Efficient and High Fidelity Speech Synthesis", *NeurIPS* 2020.

**Speaker recognition**

- Kinnunen & Li, "An Overview of Text-Independent Speaker Recognition", *Speech Communication*, 52(1), 2010 — the classical survey.
- Dehak et al., "Front-End Factor Analysis for Speaker Verification", *IEEE TASLP*, 2011 — i-vector.
- Snyder et al., "X-Vectors: Robust DNN Embeddings for Speaker Recognition", *ICASSP* 2018.
- Desplanques, Thienpondt & Demuynck, "ECAPA-TDNN: Emphasized Channel Attention, Propagation and Aggregation in TDNN Based Speaker Verification", *Interspeech* 2020.
- Wu et al., "Anti-Spoofing for Text-Independent Speaker Verification: An Overview", *arXiv* 2020.

**VoIP**

- Hersent, Petit & Gurle, *Beyond VoIP Protocols*, Wiley, 2005.
- Perkins, *RTP: Audio and Video for the Internet*, Addison-Wesley, 2003.
- ITU-T Recommendations G.114 (latency), G.711, G.729, G.722.2 (AMR-WB), G.718, G.722.1.
- Valin, Vos & Terriberry, "Definition of the Opus Audio Codec", RFC 6716, IETF, 2012.

**Aids to the handicapped**

- Loizou, *Speech Enhancement: Theory and Practice*, 2nd ed., CRC Press, 2013 — cochlear implant and hearing-aid algorithms.
- Levitt, *Aural Rehabilitation*, Singular Publishing, 2007.

---

## Suggested tutorial and lab problems for this chapter

1. **DTW implementation.** Implement DTW from scratch (recursive and iterative). Test on two sequences with known ground-truth alignment. Verify against `fastdtw` library.

2. **DTW for isolated digits.** Record 10 examples of each of 10 digits (0–9). Use 5 as templates, 5 as tests. Compute MFCCs, apply DTW, and report digit recognition accuracy. How does accuracy change with 1, 3, 5 templates?

3. **Sakoe–Chiba band and LB_Keogh.** Add the Sakoe–Chiba band to your DTW. Measure the speedup. Then implement LB_Keogh and verify it prunes candidates correctly (lower bound never exceeds true DTW).

4. **HMM forward algorithm.** Implement the forward algorithm from scratch. Verify against `hmmlearn` on a small test HMM. Test the **underflow** issue on a long sequence — then fix it with **scaling** or **log-domain**.

5. **HMM Viterbi.** Implement Viterbi. Compare the path with forward's most probable state (via max posterior). Verify that for many sequences, they coincide.

6. **HMM Baum–Welch.** Train a 3-state HMM on repeated examples of a single word. Then use it for isolated word recognition (one HMM per word). Report per-word accuracy.

7. **GMM-HMM ASR.** Use `hmmlearn` or `python_speech_features` to build a small ASR for digits. Train 10 GMM-HMMs (one per digit), test on a separate set. Compare with DTW accuracy.

8. **TTS with `pyttsx3`.** Generate TTS for a paragraph. Vary speaking rate and pitch. Measure MOS informally by having 5 friends rate 10 samples each. Compare with the "expected" MOS of the engine.

9. **MOS evaluation.** Collect MOS ratings from 10 listeners on 20 TTS samples (5 from formant, 5 from concatenative, 5 from neural, 5 from human). Compute means and 95% CIs. Which pairs are statistically distinguishable?

10. **Speaker verification features.** Extract MFCCs (with deltas) for 5 speakers saying "hello" 10 times each. Train a GMM per speaker (5 Gaussians) on 8 recordings, test on 2. Compute per-speaker accuracy and EER across the 5 speakers.

11. **EER calculation.** For the speaker-verification system of Problem 10, plot the DET curve (FAR vs. FRR over all thresholds). Mark the EER point. Compare with the human EER baseline (~2–3%).

12. **VoIP bandwidth and latency.** Compute the bandwidth overhead of RTP/UDP/IP headers for a G.711 call, a G.729 call, and an Opus call at 24 kbps. Then compute the latency budget for a call with 60 ms network delay, 30 ms jitter buffer, and 25 ms codec delay. Is it under the 150 ms G.114 threshold?

13. **Anti-spoofing exploration.** Use a TTS engine to generate speech in a target voice (from a sample recording — use `pyttsx3` for a basic attempt, or a cloud API for higher quality). Test whether a simple MFCC-GMM speaker verifier accepts the synthetic sample. Discuss countermeasures.