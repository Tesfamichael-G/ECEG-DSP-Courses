# Chapter 3 — Introduction to Digital Speech Processing and Recognition

!!! abstract "Chapter at a glance"
    **Contact hours:** 2 lecture + 2 tutorial + 1 laboratory · **Home study:** 5 h
    
    **Coverage (syllabus):** digital speech processing · digital transmission and
    storage of speech · speech synthesis · speaker verification and identification ·
    speech recognition and aids to the handicapped · speech production mechanism ·
    classification of speech sounds · nature of the speech signal · models of
    speech production.

---

## Learning outcomes

By the end of this chapter you should be able to:

1. Describe the **speech production mechanism**: lungs, larynx, vocal tract, nasal
   cavity, articulators, and the acoustic consequences of each.
2. Classify speech sounds as **vowels, diphthongs, and consonants**; distinguish
   **voiced** from **unvoiced** sounds; place consonants by **place and manner**
   of articulation.
3. Analyse the **speech waveform**, spectrum, and spectrogram. Identify
   **formants**, **pitch**, and the **voiced/unvoiced distinction**.
4. State the **source-filter model** of speech production and explain how it
   underlies LPC (linear predictive coding) and formant synthesis.
5. Compute key **speech parameters**: fundamental frequency, pitch period,
   formant frequencies, **short-time energy**, **zero-crossing rate**, and speech
   bandwidth.
6. Describe the **data rates** of raw and compressed speech, and compute storage
   and transmission requirements.
7. State the four main **applications**: transmission/storage, synthesis,
   recognition, speaker verification/ID — and the technology behind each.
8. Describe **aids to the handicapped** based on speech processing (reading
   machines, prostheses, hearing aids).
9. Implement basic **speech analysis** in Python, Java, and C#.

---

## 3.1 What is speech?

**Speech** is the acoustic signal produced by the human vocal apparatus to
communicate linguistic information. It is:

- **Acoustic** — a pressure wave in air, measurable with a microphone.
- **Time-varying** — the signal changes continuously as the speaker articulates.
- **Non-stationary** — its statistical properties change over tens of milliseconds.
- **Quasi-stationary** — over short intervals (~10–30 ms), it can be
  approximated as stationary, which is the foundation of **short-time analysis**.

**Speech vs. audio.** *Speech* is human-produced, has linguistic structure
(phonemes, words, sentences), and occupies a specific bandwidth (50 Hz – 8 kHz
for most of its energy). *Audio* is a broader category (music, environmental
sounds) that shares signal-processing tools but differs in statistics and
perception.

**Language vs. speaker.** A speech signal carries two layers of information:

- **Linguistic** — what is being said (words, phonemes). Speaker-independent.
- **Paralinguistic / speaker-specific** — who is saying it (identity, emotion,
  age, gender, health). Speaker-specific.

This two-layer structure is why the same signal supports both **speech
recognition** (extract the words) and **speaker recognition** (extract the
identity) — the two tasks are theoretically separate even though they operate
on the same waveform.

---

## 3.2 The speech production mechanism

Speech is produced by the coordinated action of four subsystems, from the
bottom up:

### 3.2.1 The lungs (power source)

The lungs provide the **airflow** that drives speech. During speech, we exhale
steadily for several seconds; the airflow is regulated by the diaphragm and
intercostal muscles. The result is a quasi-steady **volume velocity** of air
through the vocal tract — typically 100–300 cm³/s during normal phonation.

Without airflow, no speech. Whispering produces speech with airflow but without
vocal-fold vibration; normal speech uses both.

### 3.2.2 The larynx and vocal folds (sound source)

The **larynx** houses the **vocal folds** (formerly "vocal cords"), two
elastic membranes that can be:

- **Open** (abducted) — air flows freely → **unvoiced** sounds (/s/, /f/, /p/, /t/, /k/).
- **Vibrating** (adducted, with tension) — air flows through in **puffs** → **voiced** sounds (/a/, /z/, /b/, /d/, /g/, and all vowels).

**The phonation mechanism.** When the vocal folds are closed and pressure builds
up beneath them, they are pushed open; the airflow then **reduces** the pressure
(Bernoulli effect), pulling them back together; the cycle repeats. The result is
a **quasi-periodic pulse train** of air puffs at the **fundamental frequency** $F_0$.

**Fundamental frequency.** The rate of vocal-fold vibration:

| Speaker | Typical $F_0$ | Pitch period $T_0 = 1/F_0$ |
|---|---|---|
| Adult male | 80–150 Hz | 6.7–12.5 ms |
| Adult female | 150–250 Hz | 4–6.7 ms |
| Child | 200–400 Hz | 2.5–5 ms |

$F_0$ is what we perceive as **pitch**. The harmonics at $k F_0$ ($k = 1, 2, 3, \ldots$)
are what give voiced speech its rich spectral structure.

**The glottal pulse.** Each cycle of vocal-fold vibration produces a **glottal
pulse** — a brief, asymmetric wave with a **fundamental period** $T_0$ and a
broadband harmonic spectrum. The glottal pulse is the **excitation** of the
speech signal.

### 3.2.3 The vocal tract (filter)

Above the larynx, the **vocal tract** — the pharynx (throat), oral cavity (mouth),
and nasal cavity (nose) — acts as an **acoustic filter** that shapes the
glottal pulse into recognizable speech sounds.

**The vocal tract as an acoustic tube.** Approximate dimensions of an adult male:

- Length: ~17 cm (a quarter-wavelength resonator → first resonance at $c/(4L) \approx 340/(4 \cdot 0.17) \approx 500$ Hz).
- Cross-section: varies from 0 to ~5 cm² along its length, controlled by the tongue, jaw, lips, and velum.
- Speed of sound in warm, humid air: $c \approx 350$ m/s.

**Formants.** The **resonances** of the vocal-tract tube are called **formants**
$F_1, F_2, F_3, \ldots$. They appear as peaks in the spectrum of voiced speech.

| Vowel | $F_1$ (Hz) | $F_2$ (Hz) | $F_3$ (Hz) |
|---|:---:|:---:|:---:|
| /i/ (beat) | 270 | 2290 | 3010 |
| /ɪ/ (bit) | 390 | 1990 | 2550 |
| /ɛ/ (bet) | 530 | 1840 | 2480 |
| /æ/ (bat) | 660 | 1720 | 2410 |
| /a/ (father) | 730 | 1090 | 2440 |
| /ɔ/ (bought) | 570 | 840 | 2410 |
| /ʊ/ (book) | 440 | 1020 | 2240 |
| /u/ (boot) | 300 | 870 | 2240 |

(Typical adult male values, varies by speaker and dialect.)

**The first two formants carry most vowel identity.** In the $(F_1, F_2)$ plane,
vowels form a "vowel quadrilateral" — the classic **vowel chart** (IPA) is
directly the image of the vowel space in the $(F_1, F_2)$ coordinates.

### 3.2.4 The articulators

The **articulators** shape the vocal tract:

- **Tongue** — the most important articulator. Its **body** (front/back) controls
  the low formants; its **tip/blade** controls high formants and constrictions
  for consonants.
- **Lips** — opening (rounded) or closing.
- **Velum (soft palate)** — closes the nasal cavity (oral sounds) or opens it
  (nasal sounds like /m/, /n/).
- **Jaw** — gross up/down movement affecting the vocal-tract length.
- **Pharynx** — constriction affects the higher formants.
- **Epiglottis** — protects the airway during swallowing.

The **configuration of the articulators at any instant** determines the **shape
of the vocal tract**, which determines the **formants**, which determines the
sound. Speech is the continuous trajectory of that configuration over time.

### Worked example 3.1 — acoustic tube resonances

Model the vocal tract as a **uniform tube** of length $L = 17$ cm, closed at the
glottis and open at the lips.

**Boundary conditions.** Closed at one end, open at the other — the **quarter-wave
resonator**.

**Resonant frequencies.**

$$
f_n = \frac{(2n - 1)\, c}{4 L},\qquad n = 1, 2, 3, \ldots
$$

With $c = 350$ m/s, $L = 0.17$ m:

$$
f_1 = \frac{350}{4 \cdot 0.17} \approx 515\ \text{Hz},
$$
$$
f_2 = \frac{3 \cdot 350}{4 \cdot 0.17} \approx 1544\ \text{Hz},
$$
$$
f_3 = \frac{5 \cdot 350}{4 \cdot 0.17} \approx 2574\ \text{Hz}.
$$

These are in the right range for a **neutral vowel** (schwa, /ə/), which is
produced with a uniform tract — the vowels /i/, /u/, /a/ all deviate from this
by articulator movement.

**Speaker variation.** A female with $L = 14$ cm has formants ~20% higher; a
child with $L = 10$ cm has formants ~70% higher. This is why men's and women's
vowel spectra look "shifted" — the same articulatory configuration gives
different formant frequencies. Recognizing vowels **speaker-independently** is
one of the hard problems of speech recognition (solved by **vocal-tract-length
normalization**).

### Worked example 3.2 — pitch period from a waveform

A voiced speech signal has an obvious periodic structure. Count the cycles
in a 100 ms window: if there are 12 cycles, then

$$
F_0 = \frac{12}{0.1\ \text{s}} = 120\ \text{Hz},\qquad
T_0 = \frac{1}{F_0} = 8.33\ \text{ms}.
$$

The pitch period $T_0$ is visible as the spacing between successive large peaks
in the waveform; $F_0$ is its reciprocal.

At a sampling rate $f_s = 16$ kHz, the pitch period is

$$
T_0 \cdot f_s = \frac{8.33\ \text{ms}}{62.5\ \mu\text{s}} \approx 133\ \text{samples}.
$$

So pitch detection in the digital domain looks for peaks separated by ~130
samples in this case. This is the **time-domain pitch estimation** covered in
Chapter 4.

---

## 3.3 Classification of speech sounds

Speech is a sequence of **phones** — minimal sound units. Each phone belongs to
one of three broad classes:

### 3.3.1 Vowels

**Vowels** are the "open" sounds: the vocal tract is **unobstructed**, airflow
is laminar, and the sound is **always voiced**. The vowel's identity is
determined by the **articulatory configuration**, which determines the
**formants** $F_1, F_2, F_3$.

**Classical articulatory description.**

- **Tongue height** (high / mid / low) ↔ $F_1$ (low $F_1$ = high vowel).
- **Tongue backness** (front / central / back) ↔ $F_2$ (high $F_2$ = front vowel).
- **Rounding** (rounded / unrounded) ↔ $F_3$ (lowers $F_3$ for rounded).

**In English**, the ~12 monophthongs occupy the vowel quadrilateral: /i, ɪ, ɛ,
æ, a, ɔ, ʊ, u, ʌ, ə, ɜ, ɑ/.

**Diphthongs** are vowels whose formants **glide** from one vowel target to
another within a syllable: /aɪ/, /aʊ/, /ɔɪ/, /eɪ/, /oʊ/. Their spectral
trajectory is what distinguishes them from the monophthongs.

### 3.3.2 Consonants

**Consonants** involve a **constriction** or **closure** of the vocal tract.
They are classified by three attributes:

**1. Voicing.** Voiced (vocal folds vibrate) vs. unvoiced (no vibration).

**2. Place of articulation.** Where the constriction occurs:

| Place | Articulators | Examples |
|---|---|---|
| Bilabial | Upper and lower lips | /p/, /b/, /m/ |
| Labiodental | Lower lip + upper teeth | /f/, /v/ |
| Dental | Tongue tip + teeth | /θ/, /ð/ |
| Alveolar | Tongue tip + alveolar ridge | /t/, /d/, /n/, /s/, /z/, /l/, /r/ |
| Postalveolar | Tongue blade + hard palate | /ʃ/, /ʒ/, /tʃ/, /dʒ/ |
| Palatal | Tongue body + hard palate | /j/ |
| Velar | Tongue back + soft palate | /k/, /g/, /ŋ/ |
| Glottal | Vocal folds | /h/, /ʔ/ |

**3. Manner of articulation.** How the airflow is modified:

| Manner | Description | Examples |
|---|---|---|
| **Plosives (stops)** | Complete closure → sudden release (burst) | /p/, /b/, /t/, /d/, /k/, /g/ |
| **Fricatives** | Narrow constriction → turbulent noise | /s/, /z/, /f/, /v/, /ʃ/, /ʒ/ |
| **Affricates** | Closure + slow release | /tʃ/, /dʒ/ |
| **Nasals** | Full oral closure + open nasal cavity | /m/, /n/, /ŋ/ |
| **Liquids (approximants)** | Partial closure, resonant | /l/, /r/ |
| **Glides (semivowels)** | Transitional, like vowels | /w/, /j/ |

### 3.3.3 Summary table

| Sound class | Voiced? | Constriction? | Pitch present? | Waveform look |
|---|:---:|:---:|:---:|---|
| **Vowels** | Yes | No (open) | Yes | Periodic, high energy |
| **Diphthongs** | Yes | Transitional | Yes | Periodic, gliding |
| **Liquids, glides** | Yes | Partial | Yes | Periodic, formant-like |
| **Nasals** | Yes | Full oral + open nasal | Yes | Periodic, damped |
| **Voiced fricatives** | Yes | Narrow | Yes + turbulent noise | Periodic + noise |
| **Unvoiced fricatives** | No | Narrow | No | Noise |
| **Voiced plosives** | Yes (closure + burst) | Full → burst | Brief | Silent closure, then burst, then periodic |
| **Unvoiced plosives** | No | Full → burst | No | Silent closure, then burst |

### Worked example 3.3 — voiced/unvoiced decision from the signal

We observe two 20 ms speech frames:

**Frame A:** periodic peaks every ~8 ms, high amplitude, smooth spectral
envelope with peaks around 700 Hz and 1200 Hz.

**Frame B:** noisy, no periodic structure, low amplitude, broadband spectrum
with energy concentrated at 4–8 kHz.

**Classification:**

- Frame A: **voiced** — periodic (vocal-fold vibration), low formants (vowel or liquid), high energy. Likely a vowel.
- Frame B: **unvoiced** — aperiodic (turbulent noise), high frequencies (fricative). Likely a fricative like /s/ or /ʃ/.

**Quantitative features** (Chapter 4):

- **Short-time energy** — high for voiced, low for unvoiced fricatives.
- **Zero-crossing rate (ZCR)** — low for voiced (smooth), high for unvoiced (rapid sign changes).
- **Autocorrelation peak** — high at $k = T_0$ for voiced, near zero for unvoiced.

These three features together give reliable voiced/unvoiced decisions — the
front end of every speech recognition system.

---

## 3.4 The nature of the speech signal

### 3.4.1 Time-domain view

The speech **waveform** $x(t)$ (or $x[n]$) has very different characteristics
depending on the phone:

- **Voiced speech:** quasi-periodic, with period $T_0 = 1/F_0$; slow envelope modulation (vocal-tract filter shaping the harmonics).
- **Unvoiced speech:** noise-like, aperiodic, low energy (except for fricatives, which are moderate).
- **Silence:** near zero (with background noise).
- **Plosives:** brief silence followed by a burst (impulse-like).

The waveform is **not directly interpretable** for linguists — this is why
the **spectrum** and **spectrogram** are more useful.

### 3.4.2 Frequency domain

The **spectrum** of a short-time speech frame (typically 10–30 ms) shows:

- **Harmonics** at $k F_0$ (from the glottal pulse) — visible as a comb of peaks.
- **Formants** — the envelope of the harmonic peaks, with resonances at $F_1, F_2, F_3$, etc.
- **High-frequency roll-off** — the glottal spectrum falls at ~-12 dB/octave.

**A typical vowel spectrum** shows the harmonics (spaced $F_0$ apart) modulated
by the **formant envelope** (broad peaks at $F_1, F_2, F_3$). The formants
determine the vowel identity; the harmonic spacing determines the pitch.

### 3.4.3 The spectrogram

The **spectrogram** is the time-frequency representation of speech — the most
important visualization in speech processing. It is computed by **short-time
Fourier analysis** (Chapter 4 §4.4) — window the signal into short frames,
compute the spectrum of each, display as an image (time on x-axis, frequency
on y-axis, magnitude as colour).

**Reading a spectrogram:**

| Visible feature | Interpretation |
|---|---|
| Horizontal bands (dark) | **Formants** — vowel identity |
| Vertical striations | **Pitch periods** — voiced speech |
| Broad diffuse energy at high frequency | **Fricatives** (/s/, /ʃ/) |
| Short vertical spike | **Plosive burst** |
| Dark region near 500 Hz + harmonics | **Nasal** (extra formant at low frequency) |
| Blank regions | Silence or very low energy |

The spectrogram is the **phonetician's primary tool** — trained phoneticians
can read the phoneme sequence from a spectrogram, though not perfectly.

### 3.4.4 Speech bandwidth and sampling

**Where is speech's energy?**

| Band | Content | Contribution |
|---|---|---|
| 50–300 Hz | Fundamental frequency, low harmonics | Pitch, voicing |
| 300–1000 Hz | First formant $F_1$ | Vowel height |
| 1000–3000 Hz | Second formant $F_2$, most consonant info | Vowel backness, consonant place |
| 3000–8000 Hz | Third formant $F_3$, fricatives, plosive bursts | Fine consonant detail |
| 8000–20000 Hz | Airborne noise, negligible speech info | Sparsely used |

**Speech bandwidths:**

- **Narrowband telephony:** 300–3400 Hz (8 kHz sampling). Intelligible but "telephone-like" — sounds like a phone call. What the PSTN uses.
- **Wideband / HD voice:** 50–7000 Hz (16 kHz sampling). Noticeably clearer; standard in VoLTE, VoIP, and modern codecs.
- **Super-wideband:** 50–14000 Hz (32 kHz sampling). Used in some conferencing systems.
- **Full-band:** 20–20000 Hz (44.1 or 48 kHz sampling). Studio quality.

**Nyquist sampling.** To capture up to $W$ Hz of speech, we need $f_s > 2W$:

| Speech bandwidth | Required sampling rate |
|---|---|
| 3.4 kHz (narrowband) | 8 kHz (standard telephony) |
| 7 kHz (wideband) | 16 kHz |
| 14 kHz (super-wideband) | 32 kHz |
| 20 kHz (full-band) | 44.1 or 48 kHz |

The 8 kHz rate for telephony is the historical legacy: 4 kHz is about the
bandwidth needed for **intelligible** speech; 8 kHz satisfies Nyquist; it's
been the standard since the 1960s and is still used in every mobile phone call.

### Worked example 3.4 — speech data rates

**Raw PCM, narrowband telephony.**

- Sampling rate: 8 kHz
- Bit depth: 8 bits (μ-law compressed) or 16 bits (linear)
- Bitrate: $8000 \times 8 = 64$ kbps (μ-law), $8000 \times 16 = 128$ kbps (linear PCM).

**Raw PCM, wideband.**

- Sampling: 16 kHz
- Bit depth: 16 bits
- Bitrate: $16000 \times 16 = 256$ kbps.

**CD-quality audio.**

- Sampling: 44.1 kHz, 16 bits, stereo
- Bitrate: $44100 \times 16 \times 2 = 1411$ kbps.

**Compressed speech codecs:**

| Codec | Rate | Speech quality |
|---|---|---|
| G.711 (μ-law PCM) | 64 kbps | Toll quality |
| G.726 (ADPCM) | 32 kbps | Toll quality |
| G.729 (CS-ACELP) | 8 kbps | Toll quality |
| AMR-NB | 4.75–12.2 kbps | Mobile quality |
| AMR-WB | 6.6–23.85 kbps | HD voice |
| Opus | 6–510 kbps | Adjustable |
| LPC-10 | 2.4 kbps | Military / very low rate |
| MELPe | 600 bps | Military / very low rate |

**Storage example.** One hour of speech:

| Codec | Bitrate | Storage for 1 hour |
|---|---|---|
| G.711 | 64 kbps | 28.8 MB |
| G.729 | 8 kbps | 3.6 MB |
| AMR-WB | 12.2 kbps | 5.5 MB |
| LPC-10 | 2.4 kbps | 1.08 MB |

The compression from G.711 (28.8 MB) to LPC-10 (1.08 MB) is ~27:1 — and the
speech is still **intelligible** (though somewhat robotic). This is the
headline of speech coding.

---

## 3.5 Models of speech production

### 3.5.1 The source-filter model

The most influential model of speech production, introduced by **Fant (1960)**:

!!! abstract "Source-filter model"
    Speech is modelled as a **glottal source** filtered by the **vocal tract**:

    $$
    X(f) = G(f)\, H(f),
    $$

    where:

    - $G(f)$ = spectrum of the **glottal pulse** (source) — broadband, with
      harmonic structure at $F_0$.
    - $H(f)$ = transfer function of the **vocal tract** (filter) — all-pole,
      with resonances at the formants.

**Consequences.**

- The **source** determines **pitch** (harmonics at $F_0$) and **voicing** (periodic
  vs. noise).
- The **filter** determines the **formants** — the vowel identity, consonant
  place of articulation.
- Source and filter are **independent** — you can whisper at the same pitch
  (source off, filter on) or speak the same vowel at different pitches (source
  changes, filter stable). This is why whispering is intelligible.

### 3.5.2 The LPC model

The **linear predictive coding (LPC)** model is a **parametric realization**
of the source-filter model. It assumes:

- **Filter:** an all-pole (AR) filter of order $p$ (typically $p = 10$ for 8 kHz
  speech, $p = 16$–$20$ for 16 kHz).
  $$
  H(z) = \frac{1}{1 + \sum_{k=1}^{p} a_k z^{-k}}
  $$
- **Source:** either a **pulse train** at pitch $F_0$ (voiced) or a **white noise**
  generator (unvoiced).

The LPC coefficients $a_k$ are estimated by **Levinson–Durbin** (Chapter 2 of
ECEG-6303) from the autocorrelation of short-time frames. The **prediction
error** $e[n]$ is the **residual** — the driving source, which is close to
the glottal pulse train (voiced) or white noise (unvoiced).

**Why LPC matters.**

- Compact: $p = 10$ coefficients + $F_0$ + gain + voicing ≈ 40–50 bits per 20 ms frame, giving 2.4 kbps.
- Fast: Levinson–Durbin is $O(p^2)$.
- Interpretable: the LPC coefficients are the AR parameters of the vocal tract; the LPC spectrum approximates the spectral envelope (formants).
- The basis of every modern speech codec — CELP, AMR, Opus, G.729.

### Worked example 3.5 — LPC spectrum

An LPC model of order $p = 2$ for a vowel has coefficients $a_1 = -1.5$,
$a_2 = 0.7$, gain $G = 1$.

**Filter.** $H(z) = 1/(1 - 1.5z^{-1} + 0.7z^{-2})$.

**Frequency response.** The denominator roots:

$$
1 - 1.5 z^{-1} + 0.7 z^{-2} = 0 \;\Longrightarrow\; z^2 - 1.5 z + 0.7 = 0 \;\Longrightarrow\; z = 0.75 \pm j 0.37.
$$

$|z| = \sqrt{0.75^2 + 0.37^2} \approx 0.837$, angle $\theta = \arctan(0.37/0.75) \approx 0.457$ rad.

**Resonant frequency.** $f = \theta/(2\pi) \cdot f_s$. At $f_s = 8$ kHz:

$$
f \approx \frac{0.457}{2\pi} \cdot 8000 \approx 582\ \text{Hz}.
$$

**Bandwidth.** Related to the pole radius: $\text{BW} \approx (-\ln|z|) \cdot f_s/\pi \approx 0.178 \cdot 8000/\pi \approx 453$ Hz.

So this 2-pole model has a resonance at **~580 Hz** with **~450 Hz bandwidth** —
a single formant. A real vowel needs $p = 10$ (5 formants × 2 poles each), giving
peaks at $F_1, F_2, F_3, F_4, F_5$.

**LPC spectrum.** $|H(e^{j\omega})|^2$ peaks at the resonance frequencies — so
plotting $1/|A(e^{j\omega})|^2$ gives the spectral envelope, from which the
formants can be read directly.

### 3.5.3 The articulatory model

A more detailed model that simulates the **physical** vocal tract: a cascade
of acoustic tubes whose cross-sectional areas vary along the tract. The
**Kelly–Lochbaum** model (1962) is the classic digital-tube model.

**Structure.** The tract is divided into $N$ short sections, each a two-way
delay line (forward and backward acoustic waves). The vocal-tract shape is
parameterized by the areas $A_1, \ldots, A_N$; changing these parameters
changes the formants — exactly like moving the articulators.

**Advantages.** Physically meaningful — the parameters are **articulatory**
(tongue position, lip rounding) rather than **acoustic** (LPC coefficients).

**Disadvantages.** More parameters (one per section), harder to estimate from
speech, computationally more expensive.

**Modern use.** Articulatory synthesis is used in research on speech production,
speech therapy (visual feedback of tongue position), and animation (talking
heads). Production TTS uses LPC or unit selection, not articulatory synthesis.

### 3.5.4 Comparison of models

| Model | Parameters | Complexity | Estimation | Use |
|---|---|:---:|---|---|
| Source-filter | $F_0$ + formants | Low | Manual or from spectrum | Speech analysis, synthesis |
| **LPC** | $p$ AR coefficients + $F_0$ + gain + voicing | Medium | Levinson–Durbin | **Standard in codecs and recognition** |
| Articulatory | Section areas (5–20) + source | High | Not typically estimated | Research, animation |

The LPC model is the **sweet spot**: compact enough for low-bit-rate coding,
physical enough to interpret, and mathematically tractable (Yule–Walker / Levinson).

---

## 3.6 Digital speech processing

### 3.6.1 From analogue to digital

Speech is a continuous-time analogue signal; digital processing requires:

1. **Anti-alias filter** — analog lowpass at $f_s/2$ (e.g., 4 kHz for 8 kHz sampling) to prevent aliasing.
2. **Sampling** at $f_s$ (e.g., 8 kHz narrowband, 16 kHz wideband).
3. **Quantization** to $b$ bits per sample (16-bit linear is standard for processing; 8-bit μ-law is standard for telephony).
4. **Encoding** (μ-law, A-law, or linear PCM) for storage / transmission.

**μ-law and A-law compression.** Human hearing is logarithmic — the perceptual
difference between quiet sounds is much smaller than the perceptual difference
between loud sounds. μ-law / A-law compress the dynamic range by 8-bit
logarithmic quantization, giving quality comparable to 12-bit linear — at
**half the bit rate** of 16-bit linear.

The standard telephone bitstream is **8 kHz × 8 bits = 64 kbps** μ-law PCM,
known as **G.711**. It is the foundation of the entire PSTN.

### 3.6.2 The processing pipeline

Typical digital speech processing pipeline:

1. **Pre-emphasis** — a first-order highpass $H(z) = 1 - \alpha z^{-1}$ with $\alpha \approx 0.95$, compensating the glottal roll-off (-6 dB/octave) and flattening the spectrum for analysis.
2. **Framing** — divide the signal into short overlapping frames (20–30 ms, 50% overlap).
3. **Windowing** — apply Hamming or Hann window to each frame.
4. **Feature extraction** — compute the features for the task (LPC, MFCC, pitch, energy, ZCR).
5. **Decision / modelling** — classification, coding, synthesis.

### 3.6.3 Pre-emphasis and framing

**Why pre-emphasis?** The glottal spectrum falls at ~-12 dB/octave. Without
compensation, the analysis emphasises low frequencies where the energy is
concentrated, and high-frequency formants (like $F_3$ and $F_4$) become hard
to resolve. The pre-emphasis filter $1 - 0.95z^{-1}$ inverts the slope,
flattening the spectrum.

**Why framing?** Speech is **quasi-stationary**: over ~20 ms, the vocal-tract
configuration is approximately constant, so the signal's statistical
properties are stationary. Framing lets us apply stationary tools (Fourier,
autocorrelation, LPC) to each frame independently.

**Why windowing?** A rectangular window has a poor sidelobe response (~-13 dB
first sidelobe), causing **spectral leakage** — energy from one frequency
bin "leaks" into neighbours. The Hamming window has -43 dB sidelobes, the Hann
window -31 dB — better for analysis.

### Worked example 3.6 — framing and overlap

Speech sampled at $f_s = 16$ kHz. Frame size $N = 400$ samples (25 ms); overlap
$= 50\%$ (200 samples).

- **Number of frames per second:** $16000 / (400 - 200) = 80$ frames/s.
- **Frame rate:** 80 Hz.
- For a 3-second utterance: 240 frames.
- **Number of features per second** (if each frame has 12 LPC coefficients): $80 \times 12 = 960$ LPC coefficients/s.

This is the standard feature rate for speech analysis — the "resolution" of any
time-varying speech model.

### 3.6.4 Speech coding (compression)

**The problem.** Raw speech at 64 kbps (G.711) is too much for mobile networks,
satellite links, and archival storage. Compression trades fidelity for bits.

**Classification of codecs.**

| Type | Approach | Bit rate | Quality |
|---|---|---|---|
| **Waveform** | Direct time-domain quantization (PCM, DPCM, ADPCM) | 64–32 kbps | Excellent |
| **Parametric** | Model-based (LPC, vocoder) | 2.4–4.8 kbps | Fair-robotic |
| **Hybrid** | Waveform + model (CELP family) | 4.8–24 kbps | Toll quality |

**The waveform / parametric trade-off.** Waveform codecs preserve the exact
waveform → high fidelity but high bit rate. Parametric codecs preserve only the
**perceptual** properties (formants, pitch) → low bit rate but lose individual
voice details (the famous "vocoder" sound). Hybrid codecs combine both —
a linear prediction model (LPC) for the formants, plus a codebook-excited
residual for detail (CELP).

**Standard codecs:**

| Standard | Technique | Rate |
|---|---|---|
| G.711 (μ-law PCM) | Waveform | 64 kbps |
| G.726 (ADPCM) | Waveform | 32 kbps |
| G.728 (LD-CELP) | Hybrid | 16 kbps |
| G.729 (CS-ACELP) | Hybrid | 8 kbps |
| AMR-NB | ACELP | 4.75–12.2 kbps |
| AMR-WB | ACELP | 6.6–23.85 kbps |
| Opus | CELP + LPC (hybrid) | 6–510 kbps |
| LPC-10 | Parametric (US DoD) | 2.4 kbps |

Modern mobile phone calls use **AMR-NB or AMR-WB** (or **EVS** for the newest
VoLTE systems), typically at 12.2–24 kbps. The PSTN at the other end uses
G.711 at 64 kbps — the two are transcoded at the network boundary.

---

## 3.7 Applications of speech processing

Speech processing has four major application families:

### 3.7.1 Digital transmission and storage

**Transmission.** Every phone call, VoIP session, and video conference carries
compressed speech. The main challenges:

- **Bandwidth efficiency** — compress to fit the channel (typically 4.8–24 kbps for mobile).
- **Robustness to errors** — channel errors must not destroy the speech (error concealment, FEC, retransmission).
- **Latency** — real-time constraints (typically < 200 ms end-to-end).
- **Jitter** — variable network delay must be buffered (de-jitter buffer).

**Storage.** Voicemail, podcasts, call archives, audiobooks, digital assistants
(Siri, Alexa) — all store compressed speech.

### 3.7.2 Speech synthesis (text-to-speech)

**Speech synthesis** converts text to spoken audio. Three generations:

**Concatenative synthesis** — splice pre-recorded speech units (phones, diphones,
syllables) from a recorded corpus. High naturalness but large storage.

**Parametric synthesis** — generate speech from a model (LPC or vocoder) with
text-derived parameters. Low storage but "robotic". The famous DECtalk used
this approach.

**Neural synthesis (modern)** — TTS systems based on deep neural networks
(WaveNet, Tacotron, Transformer-TTS). Naturalness nearly indistinguishable from
human speech. Current systems include **Google's WaveNet / Tacotron**,
**Amazon Polly**, **Microsoft Azure TTS**, **OpenAI's Voice Engine**.

**TTS pipeline.**

1. **Text normalization** — expand abbreviations, numbers, symbols.
2. **Text-to-phoneme conversion** — grapheme-to-phoneme (G2P) rules or dictionaries.
3. **Prosody prediction** — pitch, duration, stress from syntactic and semantic analysis.
4. **Waveform generation** — concatenative or neural.
5. **Post-processing** — smoothing, normalisation.

### 3.7.3 Speech recognition

**Automatic speech recognition (ASR)** converts audio to text. Three generations:

**Template matching (1950s–1980s)** — store spectral templates of words, match
the input to the closest template via **dynamic time warping (DTW)**.

**Statistical (1980s–2010s)** — model speech units with **hidden Markov models
(HMM)** and Gaussian mixtures. Train on large corpora; use Viterbi decoding.
The dominant paradigm for 30 years (Chapter 5 covers HMM in detail).

**Neural (2010s–present)** — end-to-end deep learning. **DeepSpeech**, **wav2vec**,
**Whisper**, and modern commercial systems achieve ~5% word error rate on
clean read speech. Modern ASR needs far less handcrafted feature engineering
— the neural network learns the features and the acoustic model jointly.

**ASR challenges.**

- **Speaker variability** — pitch, accent, vocal-tract length.
- **Environmental noise** — background noise, reverberation.
- **Speaking style** — casual, fast, disfluent speech.
- **Vocabulary size** — large-vocabulary continuous speech recognition (LVCSR) with 100k+ words.
- **Language model** — grammar, semantics, context.

### 3.7.4 Speaker verification and identification

**Speaker recognition** identifies the speaker from their voice.

- **Speaker identification** — "Which of the enrolled speakers is this?" (1:N classification).
- **Speaker verification** — "Is this speaker who they claim to be?" (1:1 verification).
- **Speaker diarization** — "Who spoke when?" in a multi-speaker recording.

**Features.** Speaker identity is encoded in:

- **Spectral envelope** (vocal-tract length → formant spacing).
- **Pitch statistics** (mean $F_0$, pitch dynamics).
- **Speaking rate, rhythm, prosody.**
- **Pronunciation and idiolect** (words chosen, accents).

**Classical approach.** Gaussian mixture models (GMM) on **MFCC features**, with
**universal background models (UBM)** for normalisation.

**Modern approach.** Deep neural networks (**d-vectors**, **x-vectors**,
**ECAPA-TDNN**). Very low error rates on standard benchmarks (NIST SRE).

**Application.** Authentication (banking by phone), forensic identification,
access control, personalised assistants (Alexa voice profiles).

### 3.7.5 Aids to the handicapped

Speech processing supports the disabled through:

**Reading machines for the blind** — OCR + TTS: scan text, synthesize speech.
Classic systems: Kurzweil Reading Machine (1976), modern apps: Voice Dream
Reader, Microsoft's Seeing AI.

**Hearing aids and cochlear implants** — real-time speech enhancement, noise
reduction, frequency shaping, binaural processing. Modern hearing aids use
machine-learning-based speech enhancement.

**Speech-to-text for the deaf** — ASR + display; automated captioning
(YouTube, Microsoft Teams, Google Meet).

**Text-to-speech for the mute** — "voice banking" for ALS patients (record
phrases while still able, then synthesize); Stephen Hawking's voice (DECtalk).

**Augmentative and alternative communication (AAC)** — apps that let users
select words/pictures and synthesize speech, used by non-verbal individuals.

**Speech therapy tools** — visualization of formants, pitch, loudness for
feedback; automated scoring of pronunciation.

**Parkinson's / dysarthria aids** — speech enhancement, slowed playback,
intelligibility-boosting processing.

Speech processing is one of the most impactful areas of DSP for
**human-centred technology** — the technology's value to a blind reader or a
hearing-aid user is immediate and dramatic.

### Worked example 3.7 — reading machine data rate

A blind user wants to listen to a book. The book has 400 pages, ~500 words per
page, ~5 characters per word + spaces ≈ 6 chars/word → 2400 chars/page → 960 000
characters total.

**TTS output rate:** typical English speech is ~150 words/min = 900 chars/min.

**Total playback time:** $960000 / 900 \approx 1067$ min ≈ **18 hours**.

**Storage at different rates:**

| Codec | Bitrate | Storage |
|---|---|---|
| G.711 | 64 kbps | 8.5 MB |
| G.729 | 8 kbps | 1.1 MB |
| LPC-10 | 2.4 kbps | 320 KB |
| Modern neural TTS (on-device) | typically 24–48 kbps | 3.2–6.4 MB |

For **real-time playback** the storage requirement is modest; the bigger
challenge is the TTS model itself (neural TTS models are 20–100 MB,
downloaded once).

---

## 3.8 Chapter summary

1. **Speech** is a non-stationary but **quasi-stationary** acoustic signal, carrying linguistic and speaker-specific information.
2. **Production mechanism:** lungs (power) → larynx/vocal folds (source) → vocal tract (filter) → articulators (shape).
3. **Vocal-fold vibration** at $F_0$ (80–400 Hz) produces voiced speech; vocal-fold opening produces unvoiced speech; the two are the fundamental axis of classification.
4. **Vocal-tract resonances** are the **formants** $F_1, F_2, F_3, \ldots$ — the primary determinant of vowel identity.
5. **Speech sounds** are **vowels** (open, voiced), **diphthongs** (gliding vowels), and **consonants** (constricted; classified by voicing, place, manner).
6. **Spectrogram** is the time-frequency representation, the most important visualization of speech; formants, pitch, fricatives, and plosives are all visible.
7. **Speech bandwidth:** 300–3400 Hz for telephony (8 kHz sampling); 50–7000 Hz for wideband (16 kHz); 20–20 kHz for full-band.
8. **Source-filter model:** speech = glottal source × vocal-tract filter. Explains the independence of pitch (source) and vowel identity (filter).
9. **LPC model:** all-pole vocal tract, driven by a pulse train (voiced) or noise (unvoiced); $p \approx 10$ coefficients capture the spectral envelope; estimated via Levinson–Durbin.
10. **Data rates:** raw PCM 64–256 kbps; compressed 2.4–24 kbps. Speech coding is what makes mobile telephony and VoIP possible.
11. **Four applications:** transmission/storage, synthesis (TTS), recognition (ASR), speaker ID. Modern systems use neural networks for all four.
12. **Aids to the handicapped** — reading machines, hearing aids, cochlear implants, captioning, augmentative communication, voice banking.

---

## 3.9 Quiz

### A. Multiple choice

1. **Voiced speech** is produced by:
   (a) turbulence at a constriction (b) vocal-fold vibration (c) plosive release (d) nasal resonance

2. The **fundamental frequency** $F_0$ for an adult male is typically:
   (a) 20–50 Hz (b) 80–150 Hz (c) 300–800 Hz (d) 2–5 kHz

3. **Formants** are:
   (a) harmonics of $F_0$ (b) vocal-tract resonances (c) glottal pulses (d) consonant bursts

4. The first formant $F_1$ of a **high vowel** like /i/ is:
   (a) high (b) low (c) zero (d) equal to $F_2$

5. Which consonant type is **always unvoiced**?
   (a) vowel (b) nasal (c) voiceless fricative (d) liquid

6. The **source-filter model** states that:
   (a) speech = source + filter (b) speech = source × filter (c) speech = source - filter (d) speech = filter / source

7. In the LPC model, the vocal tract is modelled as:
   (a) an FIR filter (b) an all-pole filter (c) a comb filter (d) a linear-phase filter

8. **Narrowband telephony** uses a sampling rate of:
   (a) 4 kHz (b) 8 kHz (c) 16 kHz (d) 44.1 kHz

9. **G.711** compresses speech at:
   (a) 2.4 kbps (b) 8 kbps (c) 64 kbps (d) 128 kbps

10. The **spectrogram** shows:
    (a) time only (b) frequency only (c) time-frequency magnitude (d) phase only

11. **Speaker verification** is a:
    (a) 1:N classification (b) 1:1 verification (c) N:N comparison (d) N:1 synthesis

12. The **source-filter model** implies that pitch and vowel identity are:
    (a) the same thing (b) independent (c) inversely proportional (d) uncorrelated with the signal

### B. Fill in the blank

1. The four subsystems of speech production are lungs, larynx, vocal tract, and \_\_\_\_.
2. The vocal-fold vibration rate is the \_\_\_\_ frequency $F_0$.
3. The vocal-tract resonances are called \_\_\_\_.
4. The first two formants $F_1, F_2$ determine vowel \_\_\_\_.
5. **Consonants** are classified by voicing, place, and \_\_\_\_ of articulation.
6. The **spectrogram** shows \_\_\_\_, \_\_\_\_, and \_\_\_\_ of the speech.
7. Speech telephony uses sampling at \_\_\_\_ kHz to cover 300–3400 Hz.
8. The **source-filter model** expresses the speech spectrum as $G(f) \cdot $\_\_\_\_.
9. **LPC** stands for \_\_\_\_\_\_\_\_ coding.
10. **G.729** compresses speech to \_\_\_\_ kbps.

### C. Define / explain

1. Describe the **speech production mechanism** in one paragraph. Name the four subsystems and their acoustic roles.
2. Classify vowels and consonants. For consonants, name the three articulatory dimensions and give two examples of each.
3. Draw (or describe) the **spectrogram** of a spoken word. What features correspond to formants, pitch, fricatives, plosives?
4. Explain the **source-filter model**. Why does it explain the intelligibility of whispered speech?
5. Derive the **acoustic tube resonances** for a uniform vocal tract of length $L = 17$ cm closed at the glottis and open at the lips. What are the first three resonances?
6. Define **LPC** and explain its relation to the source-filter model. What does "prediction error" correspond to physically?
7. Compute the **data rate** of raw 16-bit PCM at 16 kHz sampling. Compare with G.711, G.729, and Opus at 12 kbps. Discuss the trade-off.
8. List the four main **applications of speech processing** and give one concrete example of each.
9. Describe **aids to the handicapped** based on speech processing. Name three specific technologies and their users.
10. Compare **speaker verification** and **speaker identification**. What are the main features used for each?

---

## 3.10 Answer key

??? success "A. Multiple choice"
    1. **(b)** — vocal-fold vibration.
    2. **(b)** — 80–150 Hz.
    3. **(b)** — vocal-tract resonances.
    4. **(b)** — low (high vowel = low $F_1$).
    5. **(c)** — voiceless fricative (and voiceless plosives).
    6. **(b)** — speech = source × filter (multiplicative).
    7. **(b)** — all-pole (AR) filter.
    8. **(b)** — 8 kHz.
    9. **(c)** — 64 kbps.
    10. **(c)** — time-frequency magnitude.
    11. **(b)** — 1:1 verification.
    12. **(b)** — independent.

??? success "B. Fill in the blank"
    1. articulators (or nasal cavity)
    2. fundamental
    3. formants
    4. identity
    5. manner
    6. time; frequency; magnitude
    7. 8
    8. $H(f)$
    9. linear predictive
    10. 8

??? success "C. Define / explain — model answers"
    1. *Production mechanism:* air from the **lungs** drives the **larynx**, where the **vocal folds** either vibrate (voiced) or remain open (unvoiced). The airflow is filtered by the **vocal tract** — the pharynx, oral cavity, and nasal cavity — which resonantly shapes the spectrum according to the configuration of the **articulators** (tongue, lips, velum, jaw). The result is a speech waveform whose harmonic structure (source) is modulated by formant resonances (filter). The four subsystems: lungs (power), larynx (source), vocal tract (filter), articulators (shape control).
    2. *Vowels:* open, voiced, no constriction; identity determined by the formants $F_1, F_2$ (related to tongue height and backness). *Consonants:* constricted or closed vocal tract. Three dimensions: **(i) voicing** — voiced vs. unvoiced (e.g., /b/ vs. /p/); **(ii) place** — bilabial, alveolar, velar, etc. (e.g., /p/ bilabial, /t/ alveolar, /k/ velar); **(iii) manner** — plosive, fricative, nasal, liquid, glide (e.g., /p/ plosive, /s/ fricative, /m/ nasal).
    3. *Spectrogram of "banana" /bəˈnænə/:* (i) silent closure for /b/ — blank; (ii) burst and short silence; (iii) formant bands for /ə/ — horizontal dark bands around 500 and 1500 Hz; (iv) /n/ — formant bands shifted slightly, with a nasal murmur in the low frequency; (v) /æ/ — formants at ~660 and ~1720 Hz; (vi) /n/ again; (vii) /ə/ again. *Formants:* horizontal dark bands. *Pitch:* vertical striations within voiced segments. *Fricatives:* diffuse broadband energy at high frequency. *Plosives:* brief vertical spike at release.
    4. *Source-filter model:* the speech spectrum is the product of the glottal source spectrum (harmonic at $F_0$, with -12 dB/octave roll-off) and the vocal-tract filter spectrum (all-pole, with formant resonances). *Whispered speech:* the source is not periodic (vocal folds open) so no pitch, but the filter (vocal tract shape) is unchanged — the formants are the same, so vowel identity is preserved. The speech is intelligible because the **filter** carries the linguistic information; the **source** only carries pitch and voicing, which are prosodic not lexical.
    5. *Tube resonances:* for a tube of length $L$, closed at one end and open at the other, the resonant wavelengths are $\lambda_n = 4L/(2n-1)$ → frequencies $f_n = (2n-1) c/(4L)$. With $c = 350$ m/s, $L = 0.17$ m: $f_1 \approx 515$ Hz, $f_2 \approx 1544$ Hz, $f_3 \approx 2574$ Hz. *Interpretation:* these are the formants of a **neutral vowel** (schwa, /ə/) — the uniform tract produces equally spaced formants. Vowels like /i/ and /u/ arise from deviating the tract from this uniform shape.
    6. *LPC:* **linear predictive coding** models each speech sample as a linear combination of $p$ previous samples plus a prediction error: $x[n] = -\sum_{k=1}^p a_k x[n-k] + e[n]$. The coefficients $a_k$ are estimated to minimise the prediction-error energy (Yule–Walker equations, Levinson–Durbin). *Source-filter relation:* the AR filter $1/A(z)$ models the vocal tract; the prediction error $e[n]$ models the **glottal source** — it is close to a pulse train for voiced speech, or white noise for unvoiced. So LPC **decomposes** the speech into filter (vocal tract) and source (glottal pulse train) — a parametric realization of the source-filter model.
    7. *Raw PCM:* 16 kHz × 16 bit = 256 kbps. *G.711:* 8 kHz × 8 bit (μ-law) = 64 kbps — a 4:1 reduction with minor quality loss. *G.729:* 8 kbps — a 32:1 reduction from 256 kbps; toll quality (fine for voice but loses some nuance). *Opus at 12 kbps:* similar quality to G.729 with more modern features (variable rate, FEC, stereo). *Trade-off:* lower bit rate → more compression → more modelling assumptions → more fragile to non-speech sounds, accents, noise. Below ~2.4 kbps, only parametric vocoders work, and they sound robotic. Above ~16 kbps, waveform codecs give near-transparent quality. The "knee" of the trade-off is ~8–12 kbps.
    8. *Applications:* (i) **Transmission / storage** — compressed speech in every phone call, VoIP, audiobook, and podcast. (ii) **Speech synthesis (TTS)** — text-to-speech for virtual assistants (Siri, Alexa), screen readers, audiobook generation. (iii) **Speech recognition (ASR)** — dictation, automated captioning, voice search, virtual assistants. (iv) **Speaker recognition** — voice authentication in banking, forensic identification, personalised experiences. (Also: aids to the handicapped, which are applications of (ii) and (iii) with a social focus.)
    9. *Aids:* (i) **Reading machines** (OCR + TTS) — blind users; Kurzweil Reading Machine and modern apps (Seeing AI, Voice Dream Reader). (ii) **Hearing aids and cochlear implants** — hard-of-hearing users; real-time speech enhancement, noise reduction, frequency shaping. (iii) **Captioning / speech-to-text** — deaf users; automated captioning on YouTube, Teams, Zoom. (Additional: augmentative and alternative communication for non-verbal users; speech therapy visualisation tools; voice banking for ALS patients — Stephen Hawking's DECtalk voice.)
    10. *Verification vs. identification:* (i) **Verification** — "Is this speaker who they claim to be?" 1:1 comparison with a claimed identity. Used for **authentication** (banking by phone, access control). (ii) **Identification** — "Which enrolled speaker is this?" 1:N classification. Used for **forensics**, personalised assistants, speaker diarization. *Features:* for both, the key features are **MFCC** (mel-frequency cepstral coefficients) or **x-vectors** — short-time spectral envelopes that capture speaker-specific vocal-tract characteristics. Additional: **$F_0$ statistics** (mean and range of pitch), **speaking rate**, **prosody**, **pronunciation idiosyncrasies**. Modern systems use deep networks trained on speaker classification (x-vectors, ECAPA-TDNN) for state-of-the-art accuracy.

---

## Further reading

- Rabiner & Schafer, *Digital Processing of Speech Signals*, Prentice Hall, 1978 — the classical textbook; Chapters 1–3 cover everything in this chapter in more depth.
- O'Shaughnessy, *Speech Communications: Human and Machine*, 2nd ed., IEEE Press, 2000 — excellent on speech production and perception.
- Deller, Proakis & Hansen, *Discrete-Time Processing of Speech Signals*, IEEE Press, 2000 — the standard graduate reference.
- Quatieri, *Discrete-Time Speech Signal Processing*, Prentice Hall, 2001 — exceptional on source-filter and homomorphic methods.
- Flanagan, *Speech Analysis Synthesis and Perception*, 2nd ed., Springer, 1972 — the classic physical-acoustics reference.
- Fant, *Acoustic Theory of Speech Production*, Mouton, 1960 — the original source-filter model paper.
- Stevens, *Acoustic Phonetics*, MIT Press, 1998 — the definitive reference on the acoustic phonetics of English.
- Ladefoged, *A Course in Phonetics*, 7th ed., Cengage, 2014 — the standard introduction to phonetics.
- Huang, Acero & Hon, *Spoken Language Processing*, Prentice Hall, 2001 — modern comprehensive reference, covers all four applications.
- Rabiner & Juang, *Fundamentals of Speech Recognition*, Prentice Hall, 1993 — the HMM reference (Chapter 5).
- Taylor, *Text-to-Speech Synthesis*, Cambridge University Press, 2009 — the standard TTS reference.
- Kinnunen & Li, "An Overview of Text-Independent Speaker Recognition", *Speech Communication*, 52(1), 2010 — the speaker-recognition survey.
- Gold & Morgan, *Speech and Audio Signal Processing*, Wiley, 2000 — the general signal-processing reference.
- `librosa` (Python audio), `python_speech_features` (MFCC), `kaldi` (ASR), `HTK` (ASR), `Praat` (phonetic analysis), MATLAB Audio Toolbox.

---

## Suggested tutorial and lab problems for this chapter

1. **Pitch period measurement.** Record yourself saying a sustained /a/ for 3 seconds. Count the peaks in a 100 ms window. Compute $F_0$. Compare with a **pitch tracker** (e.g., `librosa.pyin`).

2. **Formant estimation.** Record a vowel (/i/, /u/, /a/, /e/) and plot its spectrum. Identify $F_1$ and $F_2$ from the peaks in the spectral envelope. Plot your vowels on an $(F_1, F_2)$ chart — do they cluster the way the standard vowel quadrilateral predicts?

3. **Voiced/unvoiced classification.** Implement a simple voiced/unvoiced decision based on **short-time energy** and **zero-crossing rate**. Test on a recording containing vowels and fricatives (/s/, /ʃ/). Report classification accuracy.

4. **Acoustic tube model.** Write a Python script that simulates the vocal tract as a lossless acoustic tube (Kelly–Lochbaum) with variable cross-sectional areas. Generate vowels by moving the "tongue" position (vary the area function). Confirm the formants shift the way the vowel chart predicts.

5. **LPC analysis.** Fit an LPC model of order $p = 10$ to a vowel frame. Plot the LPC spectrum overlaid on the DFT spectrum. Verify that the LPC peaks align with the formants.

6. **Pre-emphasis.** Apply the pre-emphasis filter $1 - 0.95 z^{-1}$ to a speech signal. Plot the spectra before and after. What is the effect on the harmonic envelope and formant visibility?

7. **Framing and windowing.** Implement framing (25 ms, 50% overlap) and Hamming windowing. Plot the sequence of frames as a spectrogram. Verify that with no windowing, the spectrogram shows visible discontinuities.

8. **Data rate calculation.** A 1-minute telephone call, a 30-minute podcast, and a 90-minute audiobook — compute the storage required at G.711, G.729, and Opus (12 kbps). Which would be acceptable for a mobile download? For a music album?

9. **Speech bandwidth.** Downsample a 16 kHz speech recording to 8 kHz (with anti-alias filtering) and then to 4 kHz (the latter without anti-aliasing, to demonstrate aliasing). Compare the three versions: which is intelligible? Which has audible aliasing artifacts?

10. **Speaker features.** Compute $F_0$ statistics, MFCC, and speaking rate for two speakers (e.g., yourself and a friend) saying the same sentence. Which features are most discriminative? Would a GMM on MFCC features separate them?

11. **Codec comparison.** Record the same sentence, encode it with G.711, G.729 (or use `ffmpeg` with `libopencore_amrnb` for AMR-NB), and Opus at equivalent rates. Compare the audio quality — which sounds most natural? Which has the most artifacts?

12. **Reading machine simulation.** Take a text passage, run it through a TTS engine (e.g., `pyttsx3` or Google's TTS API), and measure the output rate (words per second). Compare with human reading speed (~150 words/min) — how much faster is TTS? Does intelligibility suffer at higher speeds?