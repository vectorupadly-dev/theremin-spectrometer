# Theremino Spectrometer V5.0 – Cross-Platform Python Edition (Linux / macOS / Windows)

A native, cross-platform port of **Theremino Spectrometer V5.0** written in Python (**PyQt6 + OpenCV + NumPy**) for USB optical spectrometers based on a CMOS sensor and a diffraction grating (e.g. the popular AliExpress spectrometers with the **Sonix USB 2.0 Camera `0c68:6464` / `USB-ZH`** controller).

---

## ✨ Key features

1. **Cross-platform UVC camera support (Plug & Play):**
   - **Linux:** native kernel `V4L2` backend (`cv2.CAP_V4L2`, `/dev/video*` devices).
   - **macOS:** native `AVFoundation` framework (`cv2.CAP_AVFOUNDATION`).
   - **Windows:** native `DirectShow` interface (`cv2.CAP_DSHOW` – the same as in the original `Theremino_Spectrometer.exe`).
   - **Demo mode:** automatic simulation of a compact fluorescent lamp (CFL Hg/Eu) spectrum with sensor noise when no camera is connected.
2. **Faithful port of the Theremino V5.0 signal processing (DSP) pipeline:**
   - **ROI** selection in per-mille units (`Start X`, `End X`, `Start Y`, `Size Y`) and horizontal mirroring (`FlipH`).
   - Per-column luminance extraction ($0.299\,R + 0.587\,G + 0.114\,B$) with a live preview of the camera strip and automatic contrast gain (AGC).
   - **Multi-point wavelength calibration** ($\text{nm}$) – fully compatible with Theremino `Calib*.txt` files, plus a built-in calibration point editor.
   - **Temporal and spatial filters:** *Rising Speed*, *Falling Speed*, a two-pass *Spatial Averaging* smoothing filter and a non-linear scale (*Log Scale*).
   - **Relative measurements:** background subtraction (*Background / Dark Frame*) and division by a reference spectrum (*Reference*) for absorption and transmission measurements.
3. **CIE 1931 colorimetry and colour temperature (CCT in Kelvin):**
   - Calculation of **CIE 1931 $(x, y)$** chromaticity coordinates by integrating the spectrum with the $2^\circ$ standard observer functions ($\bar{x}, \bar{y}, \bar{z}$).
   - Estimation of the **correlated colour temperature CCT ($\text{K}$)** (McCamy and Hernández-Andrés formulas) and of the distance from the blackbody locus **$D_{uv}$** (Ohno method).
   - Live plot of the theoretical **Planck radiation curve** for the computed CCT, plus a preview of the source colour (sRGB swatch).
   - Support for CMOS sensor spectral sensitivity correction files (`IrradianceCoeffs`, e.g. `LaoCorrections.txt`).
4. **Photometry – illuminance measurement in lux ($\text{lx}$):**
   - Integration of the spectrum with the human photopic sensitivity curve **$V(\lambda)$** (peak at $555\text{ nm}$), with the $V(\lambda)$ curve shown on the plot.
   - Compensation for exposure time ($t_{\text{exp}}$) and gain ($\text{Gain}$), with a warning when the sensor is saturated.
   - Single-point calibration against a reference (e.g. an external lux meter), saved to `calib/lux_calib.json`.
5. **Data export:**
   - Saving the spectrum and measurement parameters (including `CCT_Kelvin`, `CIE1931_xy`, `Duv`, `Illuminance_Lux`) to a Theremino-compatible `.csv` file.

---

## 📁 Project structure

```text
theremino-spectrometer/
├── spectrometer.py                  # Main GUI application (PyQt6 + OpenCV + NumPy)
├── README.md                        # Project documentation and run instructions
├── Theremino_Spectrometer_INI.txt   # Original Theremino V5.0 configuration file
├── calib/                           # Wavelength (nm) and lux meter calibration files
│   ├── Calib1.txt                   # Default 4-point calibration (326 nm – 722 nm)
│   ├── Calib2.txt
│   ├── Calib3.txt
│   ├── Calib_WEBCAM.txt
│   └── lux_calib.json               # Lux meter calibration factor (created after calibration)
└── irradiance/                      # CMOS sensor spectral sensitivity correction files
    ├── LaoCorrections.txt           # Sensor sensitivity correction from the Theremino project
    └── Coeffs_FLAT.txt              # Flat response (no correction)
```

---

## 🚀 Installation and running on each OS

Requires **Python 3.9+** and three libraries: `opencv-python`, `numpy` and `PyQt6`.

---

### 🐧 1. Linux (Ubuntu / Debian / Mint / Fedora / Arch)

#### Step 1: Install dependencies
```bash
sudo apt update
sudo apt install python3 python3-pip python3-opencv python3-numpy python3-pyqt6 acl v4l-utils
```
*(Alternatively, in a `venv` virtual environment: `pip install opencv-python numpy PyQt6`).*

#### Step 2: USB camera permissions (`/dev/video*`)
Plug the spectrometer into a USB port. You can check that the system sees it with `lsusb` (you should see e.g. `ID 0c68:6464 Sonix Technology Co., Ltd. USB 2.0 Camera`).

To give a regular user access to `/dev/video0`:
```bash
# Grants permission immediately in the current session (no restart needed):
sudo setfacl -m u:$USER:rw /dev/video*

# Adds the user to the video group permanently (takes effect after logging in again):
sudo usermod -aG video $USER
```

#### Step 3: Run the application
```bash
python3 spectrometer.py
```

---

### 🍎 2. macOS (Intel and Apple Silicon M1 / M2 / M3 / M4)

macOS has a built-in UVC driver and detects the spectrometer as soon as it is plugged in (Plug & Play).

#### Step 1: Install Python and libraries
If you use the [Homebrew](https://brew.sh/) package manager:
```bash
brew install python
pip3 install opencv-python numpy PyQt6
```

#### Step 2: Run the application
```bash
python3 spectrometer.py
```

#### Notes for macOS:
1. **Camera access permission:** On first run, macOS will show a system dialog asking whether to allow *Terminal* (or *VS Code* / *iTerm*) to access the camera. Click **Allow**.
2. **Choosing the camera on MacBooks:** On laptops with a built-in FaceTime camera, `Camera #0` is usually the camera above the screen. After plugging in the USB spectrometer, select **`Camera #1`** (or the next one) from the drop-down list in the top-left corner of the application.

---

### 🪟 3. Windows (Windows 10 / Windows 11)

Windows automatically installs the standard UVC driver when the spectrometer is plugged in and exposes it via DirectShow.

#### Step 1: Install Python and libraries
1. Download and install Python 3 from [python.org](https://www.python.org/downloads/) (tick **"Add python.exe to PATH"** during installation) or from the Microsoft Store.
2. Open a command prompt (`cmd` or `PowerShell`) and install the required packages:
   ```powershell
   pip install opencv-python numpy PyQt6
   ```

#### Step 2: Run the application
Go to the project directory in the terminal and type:
```powershell
python spectrometer.py
```

#### Notes for Windows:
* Make sure the original `Theremino_Spectrometer.exe` (or any other application using the camera) is not running in the background – on Windows only one process at a time can use a USB camera stream.
* If your laptop has a built-in webcam, switch the source in the top bar of the application from `Camera #0` to **`Camera #1`**.

---

## 🎯 How to calibrate the spectrometer correctly

### 1. Setting the image strip (ROI)
The camera preview at the top shows a horizontal rainbow of dispersed light and an orange frame:
* Adjust **`Start Y`** and **`Size Y`** so that the selection covers the brightest central band of the spectrum (without the black background above and below it).
* If the colours run the wrong way (red on the left, violet on the right), tick or untick **`Flip X (FlipH)`** – blue/violet ($400\text{ nm}$) should be on the left and red ($700\text{ nm}$) on the right.

### 2. Wavelength calibration ($\text{nm}$) with a fluorescent lamp (CFL) or a laser
1. Point the spectrometer at a **compact fluorescent lamp (CFL)** or a surface lit by one – it has very sharp, physically fixed mercury and europium peaks:
   - **Blue mercury peak (Hg):** $435.8\text{ nm}$
   - **Green mercury peak (Hg):** $546.1\text{ nm}$
   - **Orange-red europium peak (Eu):** $611.6\text{ nm}$
2. You can also use a laser pointer (e.g. violet $405\text{ nm}$, green $532\text{ nm}$, red $650\text{ nm}$).
3. Click **`Edit calibration...`** and enter the pixel positions read for the known peaks together with their wavelengths in $\text{nm}$.

### 3. Measuring colour temperature CCT ($\text{K}$) and illuminance ($\text{lx}$)
1. **Untick `Auto Exposure`** in the *V4L2 sensor controls* panel.
2. Set the **`Exp.`** slider so that the plot of the measured light source is high (about $60\text{–}85\%$ of the scale) but **not flattened at the top** (no red sensor saturation warning).
3. Fully cover the spectrometer input and click **`Background`** to subtract dark current and sensor noise.
4. Tick **`Correction (LaoCorrections)`** to compensate for the CMOS sensor's drop in sensitivity at the edges of the visible band.
5. To scale the **lux meter ($\text{lx}$)**: place a reference lux meter (or a phone with a lux meter app) next to the spectrometer input, click **`Calibrate with reference (lx)...`** and enter the reading in lux. The factor is saved automatically to `calib/lux_calib.json`.
