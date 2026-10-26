# Resources — ECEG-6423

## Core textbooks

**Image processing**

- Jain, *Fundamentals of Digital Image Processing*, PHI, 1995 — matches the syllabus almost exactly.
- Gonzalez & Woods, *Digital Image Processing*, 4th ed., Pearson, 2018 — modern, comprehensive, with code.
- Pratt, *Digital Image Processing*, 4th ed., Wiley, 2007 — the encyclopaedic reference.
- Sid Ahmed, *Image Processing*, McGraw-Hill, 1995.

**Speech processing**

- Rabiner & Schafer, *Digital Processing of Speech Signals*, Prentice Hall, 1978 — the classical reference; covers Ch 3–4 almost line-by-line.
- Rabiner & Juang, *Fundamentals of Speech Recognition*, Prentice Hall, 1993 — the definitive HMM reference (Ch 5).
- Huang, Acero & Hon, *Spoken Language Processing*, Prentice Hall, 2001 — modern, comprehensive.
- O'Shaughnessy, *Speech Communications: Human and Machine*, 2nd ed., IEEE Press, 2000.
- Quatieri, *Discrete-Time Speech Signal Processing*, Prentice Hall, 2001 — excellent on homomorphic and cepstral methods.

**Video processing**

- Wang, Ostermann & Zhang, *Video Processing and Communications*, Prentice Hall, 2002.
- Tekalp, *Digital Video Processing*, 2nd ed., Prentice Hall, 2015.
- Richardson, *The H.264 Advanced Video Compression Standard*, 2nd ed., Wiley, 2010.
- Sze, Budagavi & Sullivan (eds.), *High Efficiency Video Coding (HEVC)*, Springer, 2014.

## Free online materials

- **Gonzalez & Woods companion site** — [imageprocessingplace.com](https://www.imageprocessingplace.com/) — slides, MATLAB code, and problem sets for the textbook.
- **MIT 6.344 / 6.345** — *Automatic Speech Recognition* (OCW) — lectures and problem sets for Ch 5.
- **Stanford EE368** — *Digital Image Processing* (course materials online).
- **Bernd Girod's Stanford EE398B** — *Image and Video Compression* — excellent free lecture slides.
- **Hynek Hermansky's speech course notes** (Johns Hopkins) — for the front-end analysis material.
- **scikit-image tutorials** — [scikit-image.org/docs/stable/user_guide](https://scikit-image.org/docs/stable/user_guide.html).

## Software

| Language | Libraries | Purpose |
|---|---|---|
| Python | `numpy`, `scipy`, `scikit-image`, `Pillow`, `OpenCV`, `librosa`, `soundfile`, `PyWavelets` | All labs |
| MATLAB | Image Processing Toolbox, Signal Processing Toolbox, Deep Learning Toolbox | Reference implementations |
| Java | `ImageJ` (image analysis), `TarsosDSP` (audio), `JCodec` (video) | Optional |
| C# | `Math.NET Numerics`, `OpenCvSharp`, `NAudio` | Optional |

Install Python:

```bash
pip install numpy scipy scikit-image Pillow opencv-python librosa soundfile PyWavelets matplotlib jupyter
```

## Datasets

| Dataset | Use | Link |
|---|---|---|
| **USC-SIPI Image Database** | Standard image-processing benchmarks | [sipi.usc.edu/database](https://sipi.usc.edu/database/) |
| **BSDS500** | Berkeley segmentation dataset | [berkeley.edu](https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html) |
| **Kodak lossless** | 24 colour images for compression testing | [r0k.us/graphics/kodak](http://r0k.us/graphics/kodak/) |
| **TIMIT** | Speech corpus (American English) | LDC (license required) |
| **LibriSpeech** | Large open speech corpus | [openslr.org/12](https://www.openslr.org/12/) |
| **NOIZEUS** | Speech-in-noise corpus | [ecs.utdallas.edu/loizou/speech/noizeus](https://ecs.utdallas.edu/loizou/speech/noizeus/) |
| **Xiph.org video test media** | Video compression testing | [media.xiph.org](https://media.xiph.org/) |
| **UVG video sequences** | High-resolution video benchmarks | [ultravideo.fi](https://ultravideo.fi/) |