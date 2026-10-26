# Software Setup

Shared setup instructions for **ECEG-6303** and **ECEG-6423**.

---

## Python (primary language for both courses)

### Installation

1. Download **Python 3.11 or later** from [python.org](https://www.python.org/downloads/).
2. During installation, **check "Add Python to PATH"** (Windows).
3. Verify: `python --version` (or `python3 --version` on macOS/Linux).

### Virtual environment (recommended)

```bash
# Create environment
python -m venv eceg-env

# Activate
source eceg-env/bin/activate         # macOS / Linux
ecg-env\Scripts\activate             # Windows

# Later, to exit:
deactivate
```

### Required packages

```bash
pip install numpy scipy matplotlib jupyter
pip install scikit-image Pillow opencv-python
pip install librosa soundfile
pip install PyWavelets
pip install statsmodels
pip install hmmlearn
pip install simpy
pip install padasip
```

### Optional (for advanced topics)

```bash
pip install cvxpy
pip install filterpy pykalman
pip install markovchain PyDTMC
pip install python_speech_features
pip install praat-parselmouth
pip install pydub
pip install moviepy
```

### Verify installation

```python
import numpy, scipy, matplotlib, sklearn, librosa, cv2
print(f"numpy: {numpy.__version__}")
print(f"scipy: {scipy.__version__}")
print(f"matplotlib: {matplotlib.__version__}")
print(f"sklearn: {sklearn.__version__}")
print(f"librosa: {librosa.__version__}")
print(f"opencv: {cv2.__version__}")
```

---

## MATLAB (required for ECEG-6423 labs)

Available through your university's licence — contact IT services.

**Required toolboxes:**

- **Image Processing Toolbox** (Chapter 2)
- **Signal Processing Toolbox** (Chapter 4)
- **Audio Toolbox** (Chapter 4–5)
- **Computer Vision Toolbox** (Chapter 6, optional)
- **Deep Learning Toolbox** (Chapters 5–6, optional)

**Verify:**

```matlab
ver
```

lists all installed toolboxes.

---

## Java (optional)

If you prefer Java for the labs:

1. Install **JDK 17** from [adoptium.net](https://adoptium.net/).
2. Maven or Gradle for dependency management.

**Required dependencies** (`pom.xml`):

```xml
<dependency>
  <groupId>org.apache.commons</groupId>
  <artifactId>commons-math3</artifactId>
  <version>3.6.1</version>
</dependency>
<dependency>
  <groupId>com.github.wendykierp</groupId>
  <artifactId>JTransforms</artifactId>
  <version>3.1</version>
</dependency>
<dependency>
  <groupId>org.openpnp</groupId>
  <artifactId>opencv</artifactId>
  <version>4.7.0-0</version>
</dependency>
```

---

## C# (optional)

1. Install **.NET 8 SDK** from [dotnet.microsoft.com](https://dotnet.microsoft.com/).
2. Create a console project: `dotnet new console -o MyLab`.

**Required packages:**

```bash
dotnet add package MathNet.Numerics
dotnet add package OpenCvSharp4
dotnet add package OpenCvSharp4.runtime.win
dotnet add package NAudio
```

---

## Getting the course materials

Clone the course repo:

```bash
git clone https://github.com/tesfamichael-g/ECEG-DSP-Courses.git
cd ECEG-DSP-Courses
```

The `docs/ecg-6303/` and `docs/ecg-6423/` folders contain the lecture notes
and code examples. The `docs/ecg-6303/chapters/` and `docs/ecg-6423/chapters/`
folders contain the worked Python/Java/C# code snippets from each chapter.

---

## Jupyter notebooks

Both courses use Jupyter for interactive labs.

```bash
jupyter notebook
```

Then navigate to `notebooks/` in the course folder. Each chapter's lab
exercises are provided as `.ipynb` files.

**Tip:** use **JupyterLab** instead of classic Notebook for a better UI:

```bash
pip install jupyterlab
jupyter lab
```

---

## Data and audio files

- **Speech samples:** `assets/audio/` in the course folder
- **Image samples:** `assets/images/` — Kodak, USC-SIPI, BSDS500 subsets
- **Video samples:** `assets/video/` — Xiph test clips, UVG sequences

Download the full datasets from the links in the [Resources](resources.md)
page if you want more extensive testing.

---

## Recommended IDEs

| Language | IDE | Why |
|---|---|---|
| **Python** | VS Code, PyCharm, JupyterLab | Best for scientific computing |
| **MATLAB** | MATLAB IDE | Built-in |
| **Java** | IntelliJ IDEA, Eclipse | Standard |
| **C#** | Visual Studio, VS Code | Standard |

**VS Code** with the Python + Jupyter extensions is the recommended setup
for most students.

---

## Common installation issues

### `pip install` fails with permission error

Use a virtual environment (as above) or `pip install --user`.

### `numpy` import fails on macOS M1/M2

Use the arm64 Python from python.org, or `pip install --upgrade numpy`.

### `opencv-python` fails on Windows

Try `pip install opencv-python-headless` (no GUI dependencies).

### `librosa` import fails

`librosa` has many dependencies (numba, llvmlite). Sometimes `pip install librosa` fails; try installing dependencies manually:

```bash
pip install numba
pip install soundfile
pip install librosa
```

### `import` errors after installing packages

Restart the Python kernel (in Jupyter) or the terminal (in scripts). Package
imports are cached at Python startup.

---

## Testing the environment

Run this quick test to verify everything is working:

```python
import numpy as np
import scipy
import matplotlib.pyplot as plt
import librosa
import cv2
from scipy.signal import welch, firwin
from scipy.fft import dct

# Quick test: generate noise, filter, plot
x = np.random.randn(1000)
b = firwin(51, 0.1)
y = scipy.signal.lfilter(b, 1, x)

f, pxx = welch(y, fs=1000, nperseg=256)

fig, ax = plt.subplots()
ax.semilogy(f, pxx)
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("PSD")
ax.set_title("Welch PSD — test")
plt.savefig("test_plot.png")
print("✓ All packages working")
```

If this runs without errors and produces `test_plot.png`, the environment is
ready.