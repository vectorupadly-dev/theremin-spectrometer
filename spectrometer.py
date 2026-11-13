#!/usr/bin/env python3
"""
Theremino Spectrometer V5.0 - Native Linux Port (Python + PyQt6 + OpenCV)
Supports UVC cameras (e.g. USB-ZH Sonix 0c68:6464 from AliExpress), Theremino calibration files,
spectral sensitivity correction (Irradiance Coeffs), correlated colour temperature CCT (K) estimation,
CIE 1931 (x, y) chromaticity coordinates, Duv deviation and illuminance measurement in lux (lx)
with exposure time and gain compensation plus single-point calibration.
"""

import json
import math
import os
import platform
import sys
import time
from pathlib import Path
from typing import List, Tuple, Optional

import cv2
import numpy as np
from PyQt6.QtCore import Qt, QTimer, QPointF
from PyQt6.QtGui import (
    QImage, QPixmap, QPainter, QPen, QColor, QFont, QPolygonF
)
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QSpinBox, QDoubleSpinBox, QCheckBox, QPushButton, QComboBox,
    QGroupBox, QFileDialog, QMessageBox, QDialog, QTableWidget,
    QTableWidgetItem, QHeaderView, QSlider, QFrame, QInputDialog
)


def wavelength_to_rgb(wavelength: float) -> Tuple[int, int, int]:
    """Exact port of the WavelengthToColor function from Theremino Module_Utils.vb."""
    if 380 <= wavelength < 440:
        r = -(wavelength - 440.0) / (440.0 - 380.0)
        g = 0.0
        b = 1.0
    elif 440 <= wavelength < 490:
        r = 0.0
        g = (wavelength - 440.0) / (490.0 - 440.0)
        b = 1.0
    elif 490 <= wavelength < 510:
        r = 0.0
        g = 1.0
        b = -(wavelength - 510.0) / (510.0 - 490.0)
    elif 510 <= wavelength < 580:
        r = (wavelength - 510.0) / (580.0 - 510.0)
        g = 1.0
        b = 0.0
    elif 580 <= wavelength < 645:
        r = 1.0
        g = -(wavelength - 645.0) / (645.0 - 580.0)
        b = 0.0
    elif 645 <= wavelength <= 780:
        r = 1.0
        g = 0.0
        b = 0.0
    else:
        r = g = b = 0.0

    ir2, ir1, uv2, uv1 = 780.0, 650.0, 420.0, 380.0
    if wavelength > ir2 or wavelength < uv1:
        factor = 0.0
    elif wavelength > ir1:
        factor = (ir2 - wavelength) / (ir2 - ir1)
    elif wavelength < uv2:
        factor = (wavelength - uv1) / (uv2 - uv1)
    else:
        factor = 1.0

    return (
        int(max(0, min(255, 255 * r * factor))),
        int(max(0, min(255, 255 * g * factor))),
        int(max(0, min(255, 255 * b * factor))),
    )


def piecewise_interpolate(x_arr: np.ndarray, y_arr: np.ndarray, x_vals: np.ndarray) -> np.ndarray:
    """
    Exact port of the Interpolate function from Theremino Module_Calibrations.vb
    (linear interpolation with extrapolation beyond the end points).
    """
    x_arr = np.asarray(x_arr, dtype=np.float64)
    y_arr = np.asarray(y_arr, dtype=np.float64)
    x_vals = np.asarray(x_vals, dtype=np.float64)

    if len(x_arr) < 2:
        return np.zeros_like(x_vals)

    idx = np.searchsorted(x_arr, x_vals, side="left")
    idx = np.clip(idx, 1, len(x_arr) - 1)

    x0 = x_arr[idx - 1]
    x1 = x_arr[idx]
    y0 = y_arr[idx - 1]
    y1 = y_arr[idx]

    dx = np.where(np.abs(x1 - x0) < 1e-12, 1e-12, x1 - x0)
    return y0 + ((x_vals - x0) / dx) * (y1 - y0)


def cie1931_cmf(wavelengths_nm: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Analytic approximation of the CIE 1931 2-degree standard observer colour matching
    functions (Wyman, Sloan, Shirley 2013).
    The y_bar(lambda) component is also the photopic luminosity function V(lambda), peaking at 555 nm.
    """
    w = np.asarray(wavelengths_nm, dtype=np.float64)

    def pw_gauss(x, amp, mu, s1, s2):
        s = np.where(x < mu, s1, s2)
        return amp * np.exp(-0.5 * ((x - mu) * s) ** 2)

    x_bar = (
        pw_gauss(w, 1.056, 599.8, 0.0264, 0.0323)
        + pw_gauss(w, 0.362, 442.0, 0.0624, 0.0374)
        - pw_gauss(w, 0.065, 501.1, 0.0490, 0.0382)
    )
    y_bar = (
        pw_gauss(w, 0.821, 568.8, 0.0213, 0.0247)
        + pw_gauss(w, 0.286, 530.9, 0.0613, 0.0322)
    )
    z_bar = (
        pw_gauss(w, 1.217, 437.0, 0.0845, 0.0278)
        + pw_gauss(w, 0.681, 459.0, 0.0385, 0.0725)
    )

    mask = (w >= 380.0) & (w <= 780.0)
    return (
        np.where(mask, np.maximum(0.0, x_bar), 0.0),
        np.where(mask, np.maximum(0.0, y_bar), 0.0),
        np.where(mask, np.maximum(0.0, z_bar), 0.0),
    )


def analyze_colorimetry_and_photometry(
    wavelengths_nm: np.ndarray, intensities: np.ndarray
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Tuple[int, int, int], float]:
    """
    Computes from the spectrum:
    - CCT (K): correlated colour temperature
    - CIE 1931 (x, y): chromaticity coordinates
    - Duv: distance from the Planckian locus
    - sRGB: approximate colour of the light source
    - photopic_integral: photopic integral int S(lambda) * V(lambda) d_lambda for lux (lx) measurement
    """
    if len(wavelengths_nm) < 10 or len(intensities) < 10:
        return None, None, None, None, (0, 0, 0), 0.0

    # Sort in ascending wavelength order
    if wavelengths_nm[0] > wavelengths_nm[-1]:
        w = wavelengths_nm[::-1]
        spd = intensities[::-1]
    else:
        w = wavelengths_nm
        spd = intensities

    vis_mask = (w >= 380.0) & (w <= 780.0)
    if int(np.sum(vis_mask)) < 10:
        return None, None, None, None, (0, 0, 0), 0.0

    spd_clean = np.maximum(0.0, spd)
    x_bar, y_bar, z_bar = cie1931_cmf(w)

    trapz_fn = getattr(np, "trapezoid", getattr(np, "trapz", None))
    X = float(trapz_fn(spd_clean * x_bar, w))
    Y = float(trapz_fn(spd_clean * y_bar, w))  # Photopic integral int S(lambda)*V(lambda) d_lambda
    Z = float(trapz_fn(spd_clean * z_bar, w))

    photopic_integral = max(0.0, Y)

    spd_vis = spd_clean[vis_mask]
    if float(np.max(spd_vis)) < 150.0 or (X + Y + Z) <= 1e-6 or Y <= 1e-6:
        return None, None, None, None, (0, 0, 0), photopic_integral

    total = X + Y + Z
    x = X / total
    y = Y / total

    denom = 0.1858 - y
    if abs(denom) < 1e-6:
        return None, x, y, None, (0, 0, 0), photopic_integral

    n = (x - 0.3320) / denom
    cct = 449.0 * (n ** 3) + 3525.0 * (n ** 2) + 6823.3 * n + 5520.33

    if cct > 10000.0 or cct < 1500.0:
        denom_ha = y - 0.1870
        if abs(denom_ha) > 1e-6:
            n_ha = (x - 0.3366) / denom_ha
            cct_ha = (
                -949.86315
                + 6253.80338 * math.exp(-n_ha / 0.92159)
                + 28.70599 * math.exp(-n_ha / 0.20039)
                + 0.00004 * math.exp(-n_ha / 0.07125)
            )
            if 1000.0 <= cct_ha <= 50000.0:
                cct = cct_ha

    cct = float(np.clip(cct, 1000.0, 50000.0))

    # Duv calculation (Ohno 2014)
    uv_denom = -2.0 * x + 12.0 * y + 3.0
    if abs(uv_denom) > 1e-6:
        u = 4.0 * x / uv_denom
        v = 6.0 * y / uv_denom
        lfp = math.hypot(u - 0.292, v - 0.24)
        a = math.acos(np.clip((u - 0.292) / max(1e-9, lfp), -1.0, 1.0))
        lbb = (
            -0.00616793 * (a ** 6)
            + 0.0893944 * (a ** 5)
            - 0.5179722 * (a ** 4)
            + 1.5317403 * (a ** 3)
            - 2.4243787 * (a ** 2)
            + 1.925865 * a
            - 0.471106
        )
        duv = float(lfp - lbb)
    else:
        duv = 0.0

    # Approximate sRGB colour
    Xn, Yn, Zn = X / Y, 1.0, Z / Y
    r_lin = 3.2406 * Xn - 1.5372 * Yn - 0.4986 * Zn
    g_lin = -0.9689 * Xn + 1.8758 * Yn + 0.0415 * Zn
    b_lin = 0.0557 * Xn - 0.2040 * Yn + 1.0570 * Zn

    max_rgb = max(r_lin, g_lin, b_lin, 1e-6)
    r_lin = max(0.0, r_lin / max_rgb)
    g_lin = max(0.0, g_lin / max_rgb)
    b_lin = max(0.0, b_lin / max_rgb)

    srgb = (
        int(np.clip((r_lin ** (1.0 / 2.2)) * 255.0, 0, 255)),
        int(np.clip((g_lin ** (1.0 / 2.2)) * 255.0, 0, 255)),
        int(np.clip((b_lin ** (1.0 / 2.2)) * 255.0, 0, 255)),
    )

    return cct, x, y, duv, srgb, photopic_integral


def planck_relative_spectrum(wavelengths_nm: np.ndarray, temp_k: float, target_max: float) -> np.ndarray:
    """
    Computes the theoretical blackbody radiation spectrum (Planck's law)
    for a given temperature T (in Kelvin), normalised over the visible band for plotting.
    """
    w_m = np.maximum(1e-9, wavelengths_nm * 1e-9)
    c2 = 1.438776877e-2  # hc / k_B [m * K]
    exponent = np.clip(c2 / (w_m * max(500.0, temp_k)), 1e-3, 100.0)
    b_lambda = (1.0 / (w_m ** 5)) / (np.expm1(exponent))

    vis_mask = (wavelengths_nm >= 400.0) & (wavelengths_nm <= 700.0)
    max_b = float(np.max(b_lambda[vis_mask])) if np.any(vis_mask) else float(np.max(b_lambda))
    if max_b <= 0:
        return np.zeros_like(wavelengths_nm)
    return (b_lambda / max_b) * (target_max * 0.8)


def decode_exposure_seconds(exp_val: int, is_v4l2_abs: bool = True) -> float:
    """
    Converts a UVC camera exposure value into exposure time in seconds:
    - On Linux (V4L2 Exposure Time, Absolute: 1..5000): 1 unit = 100 µs (0.0001 s),
      e.g. 157 = 15.7 ms, 5000 = 500 ms.
    - On Windows (DirectShow: -13..0): t = 2^(exp_val) s.
    """
    if is_v4l2_abs and exp_val > 0:
        return float(exp_val) * 0.0001
    return float(2.0 ** float(exp_val))


class SpectrumEngine:
    """Spectrometer signal processing engine (port of Module_Spectrometer.vb + CCT/CIE + Lux)."""

    def __init__(self):
        # ROI parameters (per mille 0..1000, defaults from Theremino_Spectrometer_INI.txt)
        self.flip_h: bool = True
        self.start_x: int = 0
        self.end_x: int = 1000
        self.start_y: int = 503
        self.size_y: int = 98

        # Filter parameters
        self.spatial_averaging: int = 0      # 0..10
        self.rising_speed: int = 100         # 1..100
        self.falling_speed: int = 100        # 1..100
        self.log_scale_param: float = 0.0    # -10..10

        # Normalised calibration (0..1 -> nm), defaults to Theremino Calib1.txt
        self.calib_bins_norm: List[float] = [0.0, 0.350217335390947, 0.517448481289962, 1.0]
        self.calib_nm: List[float] = [326.0, 456.130509259259, 524.555941358025, 722.0]

        # Sensor sensitivity correction (Irradiance Coeffs from Module_IrradianceCoeffs.vb)
        self.irradiance_enabled: bool = False
        self.irrad_nm: List[float] = [400.0, 800.0]
        self.irrad_coeffs: List[float] = [1.0, 1.0]

        # Exposure and gain parameters for lux meter compensation
        self.auto_exposure_active: bool = True
        self.exposure_step: int = 157 if platform.system() == "Linux" else -5
        self.camera_gain: int = 0            # 0..64 (Brightness / Gain)

        # Lux meter calibration factor (default estimate for slit + diffraction grating)
        self.lux_calib_factor: float = 0.00018
        self.lux_calibrated_by_user: bool = False

        # Working arrays
        self.num_samples: int = 1920
        self.up_down_filtered: np.ndarray = np.zeros(self.num_samples, dtype=np.float32)
        self.spatial_filtered: np.ndarray = np.zeros(self.num_samples, dtype=np.float32)
        self.background: Optional[np.ndarray] = None
        self.reference: Optional[np.ndarray] = None

        # Results for the visible slice
        self.visible_nm: np.ndarray = np.linspace(326.0, 722.0, self.num_samples)
        self.visible_samples: np.ndarray = np.zeros(self.num_samples, dtype=np.float32)
        self.max_peak_val: float = 0.1
        self.max_peak_nm: float = 0.0
        self.max_visible_val: float = 1000.0
        self.peak_area: float = 0.0

        # Colorimetric (CCT / CIE 1931) and photometric (Lux) results
        self.cct_kelvin: Optional[float] = None
        self.cie_x: Optional[float] = None
        self.cie_y: Optional[float] = None
        self.duv: Optional[float] = None
        self.swatch_rgb: Tuple[int, int, int] = (0, 0, 0)

        self.raw_photopic_rate: float = 0.0  # Photopic integral compensated for exposure time and gain
        self.lux_value: float = 0.0          # Illuminance [lx]
        self.sensor_saturated: bool = False  # Whether the sensor is saturated (pixels = 255)

        # Camera slice for preview (RGB)
        self.roi_preview_rgb: Optional[np.ndarray] = None
        self.src_x0: int = 0
        self.src_dx: int = 1920
        self.src_y0: int = 0
        self.src_dy: int = 100

    def load_calibration_file(self, filepath: str) -> bool:
        """Loads a Calib*.txt file in Theremino format."""
        try:
            bins_norm = None
            nms = None
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("CalibrationBins="):
                        val_str = line.split("=", 1)[1].strip()
                        bins_norm = [float(x) for x in val_str.split("|") if x.strip()]
                    elif line.startswith("CalibrationNanometers="):
                        val_str = line.split("=", 1)[1].strip()
                        nms = [float(x) for x in val_str.split("|") if x.strip()]
            if bins_norm and nms and len(bins_norm) == len(nms) and len(bins_norm) >= 2:
                self.calib_bins_norm = bins_norm
                self.calib_nm = nms
                return True
        except Exception as e:
            print(f"Error loading calibration: {e}")
        return False

    def load_irradiance_file(self, filepath: str) -> bool:
        """Loads a spectral sensitivity correction file (e.g. Theremino LaoCorrections.txt)."""
        try:
            nms = []
            coeffs = []
            with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("-") or "nm" in line.lower():
                        continue
                    parts = line.replace(",", ".").split()
                    if len(parts) >= 2:
                        try:
                            nms.append(float(parts[0]))
                            coeffs.append(float(parts[1]))
                        except ValueError:
                            pass
            if len(nms) >= 2:
                self.irrad_nm = nms
                self.irrad_coeffs = coeffs
                return True
        except Exception as e:
            print(f"Error loading irradiance correction: {e}")
        return False

    def load_lux_calibration(self, filepath: str):
        if not os.path.exists(filepath):
            return
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                data = json.load(f)
                self.lux_calib_factor = float(data.get("lux_calib_factor", self.lux_calib_factor))
                self.lux_calibrated_by_user = bool(data.get("calibrated_by_user", True))
        except Exception:
            pass

    def save_lux_calibration(self, filepath: str):
        try:
            with open(filepath, "w", encoding="utf-8") as f:
                json.dump(
                    {
                        "lux_calib_factor": self.lux_calib_factor,
                        "calibrated_by_user": self.lux_calibrated_by_user,
                    },
                    f,
                    indent=2,
                )
        except Exception as e:
            print(f"Error saving lux meter calibration: {e}")

    def bin_to_nm(self, bins: np.ndarray) -> np.ndarray:
        calib_bins = np.asarray(self.calib_bins_norm, dtype=np.float64) * self.num_samples
        return piecewise_interpolate(calib_bins, np.asarray(self.calib_nm), bins)

    def nm_to_bin(self, nms: np.ndarray) -> np.ndarray:
        calib_bins = np.asarray(self.calib_bins_norm, dtype=np.float64) * self.num_samples
        return piecewise_interpolate(np.asarray(self.calib_nm), calib_bins, nms)

    def process_frame(self, frame_bgr: np.ndarray):
        """Processes a single camera frame following the Theremino pipeline + computes CCT and Lux."""
        if self.flip_h:
            frame_bgr = cv2.flip(frame_bgr, 1)

        src_h, src_w = frame_bgr.shape[:2]
        if src_w != self.num_samples or len(self.up_down_filtered) != src_w:
            self.num_samples = src_w
            self.up_down_filtered = np.zeros(src_w, dtype=np.float32)
            self.spatial_filtered = np.zeros(src_w, dtype=np.float32)
            self.background = None
            self.reference = None

        # Compute X and Y slice (exactly as in Spectrometer_SetSourceParams)
        src_x0 = (src_w * self.start_x) // 1000
        src_dx = src_w - src_x0 + (src_w * (self.end_x - 1000)) // 1000
        if src_x0 + src_dx > src_w:
            src_dx = src_w - src_x0
        if src_dx <= 0:
            src_dx = 1
        src_x0 = max(0, min(src_x0, src_w - 1))

        src_dy = (src_h * self.size_y) // 1000
        src_y0 = src_h - (src_h * self.start_y // 1000) - src_dy
        if src_y0 + src_dy > src_h:
            src_dy = src_h - src_y0
        if src_dy <= 0:
            src_dy = 1
        src_y0 = max(0, min(src_y0, src_h - 1))

        self.src_x0, self.src_dx = src_x0, src_dx
        self.src_y0, self.src_dy = src_y0, src_dy

        # Vertical ROI slice
        roi_bgr = frame_bgr[src_y0:src_y0 + src_dy, :]
        roi_rgb = cv2.cvtColor(roi_bgr, cv2.COLOR_BGR2RGB)

        # Check whether the sensor is saturated within the ROI (> 252)
        self.sensor_saturated = bool(np.percentile(roi_rgb, 99.5) >= 252.0)

        # Compute luminance for each column (0.299 R + 0.587 G + 0.114 B)
        r_mean = roi_rgb[:, :, 0].mean(axis=0, dtype=np.float32)
        g_mean = roi_rgb[:, :, 1].mean(axis=0, dtype=np.float32)
        b_mean = roi_rgb[:, :, 2].mean(axis=0, dtype=np.float32)
        received_samples = 0.299 * r_mean + 0.587 * g_mean + 0.114 * b_mean

        # Automatic brightness gain for the preview strip (Theremino IMAGE AGC)
        max_img_val = float(np.max(received_samples)) if len(received_samples) > 0 else 1.0
        gain = min(255.0 / max(max_img_val, 0.01), 100.0)
        preview = np.clip(roi_rgb.astype(np.float32) * max(1.0, gain * 0.6), 0, 255).astype(np.uint8)
        self.roi_preview_rgb = preview

        # Scale to 0..65536 (ScaleCoeff_to65536 = 256)
        scaled = received_samples * 256.0

        # RisingSpeed / FallingSpeed temporal filter
        k_up = (self.rising_speed / 100.0) ** 1.5
        k_dw = (self.falling_speed / 100.0) ** 1.5
        diff = scaled - self.up_down_filtered
        rates = np.where(diff > 0, k_up, k_dw).astype(np.float32)
        self.up_down_filtered += diff * rates

        # Wavelength calibration (Array_UpDownFiltered_To_Array_Calibrated)
        nm_min = float(np.clip(self.bin_to_nm(np.array([0.0]))[0], 0.0, 4000.0))
        nm_max = float(np.clip(self.bin_to_nm(np.array([float(self.num_samples)]))[0], nm_min + 10.0, 4000.0))

        indices = np.arange(self.num_samples, dtype=np.float64)
        nm_grid = nm_min + (nm_max - nm_min) * indices / self.num_samples
        bins_mapped = np.clip(np.rint(self.nm_to_bin(nm_grid)).astype(np.int32), 0, self.num_samples - 1)
        calibrated = self.up_down_filtered[bins_mapped].astype(np.float32)

        # Sensor spectral sensitivity correction (CorrectForIrradiance from Module_IrradianceCoeffs.vb)
        if self.irradiance_enabled and len(self.irrad_nm) >= 2:
            coeffs = piecewise_interpolate(
                np.asarray(self.irrad_nm), np.asarray(self.irrad_coeffs), nm_grid
            ).astype(np.float32)
            calibrated = calibrated * np.maximum(0.0, coeffs)

        # Spatial filter (Add_SpatialFilter - two-pass IIR left->right and right->left)
        k_filter = 0.1 + 0.9 * ((10 - self.spatial_averaging) / 10.0)
        if k_filter >= 0.999:
            spatial = calibrated.copy()
        else:
            spatial = np.zeros_like(calibrated)
            v = 0.0
            for i in range(self.num_samples):
                v += (float(calibrated[i]) - v) * k_filter
                spatial[i] = v * 0.5
            v = 0.0
            for i in range(self.num_samples - 1, -1, -1):
                v += (float(calibrated[i]) - v) * k_filter
                spatial[i] += v * 0.5

        self.spatial_filtered = spatial.copy()

        # Background subtraction (Dark Frame / Background)
        if self.background is not None and len(self.background) == len(spatial):
            spatial = np.maximum(0.0, spatial - self.background)

        # Copy of the spectrum before dividing by the reference (for absolute lux measurement)
        spatial_abs = spatial.copy()

        # Division by the reference spectrum (Reference)
        if self.reference is not None and len(self.reference) == len(spatial):
            vis_slice = spatial[src_x0:src_x0 + src_dx]
            max_ref = max(0.1, float(np.max(vis_slice)) if len(vis_slice) > 0 else 0.1)
            ref_safe = np.where(self.reference > 0, self.reference, 1.0)
            spatial = (spatial * max_ref) / ref_safe

        # Slice shown on the plot (StartX .. EndX)
        self.visible_samples = spatial[src_x0:src_x0 + src_dx].copy()
        self.visible_nm = nm_grid[src_x0:src_x0 + src_dx].copy()

        if len(self.visible_samples) > 0:
            max_idx = int(np.argmax(self.visible_samples))
            self.max_peak_val = max(0.1, float(self.visible_samples[max_idx]))
            self.max_peak_nm = float(self.visible_nm[max_idx])
        else:
            self.max_peak_val = 0.1
            self.max_peak_nm = 0.0

        self.max_visible_val = max(1000.0, self.max_peak_val)

        # Area integral (Peak Area)
        if len(self.visible_samples) > 1:
            trapz_fn = getattr(np, "trapezoid", getattr(np, "trapz", None))
            self.peak_area = float(trapz_fn(self.visible_samples, self.visible_nm))
        else:
            self.peak_area = 0.0

        # Colorimetric analysis (CCT, CIE x,y, Duv) and photopic V(lambda) integral
        abs_visible = spatial_abs[src_x0:src_x0 + src_dx]
        cct, cx, cy, duv, srgb, photopic_int = analyze_colorimetry_and_photometry(
            self.visible_nm, abs_visible
        )

        if cct is not None and self.cct_kelvin is not None:
            self.cct_kelvin = 0.8 * self.cct_kelvin + 0.2 * cct
            self.cie_x = 0.8 * (self.cie_x or cx) + 0.2 * cx
            self.cie_y = 0.8 * (self.cie_y or cy) + 0.2 * cy
            self.duv = 0.8 * (self.duv or duv) + 0.2 * duv
        else:
            self.cct_kelvin = cct
            self.cie_x = cx
            self.cie_y = cy
            self.duv = duv
        self.swatch_rgb = srgb

        # Illuminance calculation in lux (lx)
        # Compensation for exposure time (in manual mode) and Gain
        is_linux = platform.system() == "Linux"
        if not self.auto_exposure_active:
            exp_sec = decode_exposure_seconds(self.exposure_step, is_linux)
        else:
            exp_sec = 0.0157 if is_linux else decode_exposure_seconds(-5, False)

        gain_factor = 1.0 + (max(0, self.camera_gain) / 64.0)
        current_rate = photopic_int / max(1e-6, exp_sec * gain_factor)

        self.raw_photopic_rate = 0.8 * self.raw_photopic_rate + 0.2 * current_rate
        self.lux_value = self.raw_photopic_rate * self.lux_calib_factor

    def set_background(self):
        self.background = self.spatial_filtered.copy()

    def reset_background(self):
        self.background = None

    def set_reference(self):
        ref = self.spatial_filtered.copy()
        ref[ref < 200.0] = 99999.0
        self.reference = ref

    def reset_reference(self):
        self.reference = None


class SpectrumGraphWidget(QWidget):
    """Widget drawing the spectrum plot in Theremino style + optional Planck and V(lambda) curves."""

    def __init__(self, engine: SpectrumEngine, parent=None):
        super().__init__(parent)
        self.engine = engine
        self.use_colors: bool = True
        self.show_peaks: bool = True
        self.show_dips: bool = False
        self.show_trim_points: bool = False
        self.show_planck_curve: bool = True
        self.show_vlambda_curve: bool = False
        self.cursor_pos: Optional[Tuple[int, int]] = None
        self.setMouseTracking(True)
        self.setMinimumHeight(320)

    def mouseMoveEvent(self, event):
        self.cursor_pos = (int(event.position().x()), int(event.position().y()))
        self.update()

    def leaveEvent(self, event):
        self.cursor_pos = None
        self.update()

    def _get_log_scale(self) -> float:
        log_n = self.engine.log_scale_param
        if log_n >= 0:
            return 1.0 / (1.0 + log_n / 5.0)
        return 1.6 ** (-log_n)

    def _y_from_val(self, val: float, top: float, bottom: float, max_val: float, log_scale: float) -> float:
        h = max(1.0, bottom - top)
        norm = np.clip(val / max(1.0, max_val), 0.0, 1.0)
        if abs(log_scale - 1.0) > 1e-4:
            norm = norm ** log_scale
        return bottom - norm * h

    def _val_from_y(self, y: float, top: float, bottom: float, max_val: float, log_scale: float) -> float:
        h = max(1.0, bottom - top)
        norm = np.clip((bottom - y) / h, 0.0, 1.0)
        if abs(log_scale - 1.0) > 1e-4:
            norm = norm ** (1.0 / log_scale)
        return norm * max_val

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)

        w = self.width()
        h = self.height()
        left = 48
        top = 20
        right = w - 10
        bottom = h - 25
        plot_w = max(10, right - left)
        plot_h = max(10, bottom - top)

        # AliceBlue background as in Theremino
        painter.fillRect(0, 0, w, h, QColor(240, 248, 255))

        nms = self.engine.visible_nm
        samples = self.engine.visible_samples
        if len(nms) < 2 or len(samples) < 2:
            return

        nm_start = float(nms[0])
        nm_end = float(nms[-1])
        nm_span = max(1.0, nm_end - nm_start)
        nm_coeff = plot_w / nm_span
        max_val = self.engine.max_visible_val
        log_scale = self._get_log_scale()

        # Vertical grid (wavelength nm)
        font_small = QFont("Sans", 8)
        painter.setFont(font_small)

        pen_100 = QPen(QColor(40, 40, 40), 1)
        pen_50 = QPen(QColor(110, 110, 110), 1)
        pen_10 = QPen(QColor(185, 185, 185), 1)

        start_nm_int = int(nm_start // 10) * 10
        end_nm_int = int(nm_end) + 1
        for nm_tick in range(start_nm_int, end_nm_int, 10):
            if nm_tick < nm_start:
                continue
            x = int(left + (nm_tick - nm_start) * nm_coeff)
            if x > right:
                break
            draw_txt = False
            if nm_tick % 100 == 0:
                painter.setPen(pen_100)
                painter.drawLine(x, top, x, bottom)
                draw_txt = True
            elif nm_tick % 50 == 0:
                painter.setPen(pen_50)
                painter.drawLine(x, top, x, bottom)
                if nm_coeff > 1.2:
                    draw_txt = True
            else:
                if nm_coeff > 1.5:
                    painter.setPen(pen_10)
                    painter.drawLine(x, top, x, bottom)
                if nm_coeff > 4.0:
                    draw_txt = True

            if draw_txt:
                painter.setPen(QColor(0, 0, 0))
                painter.drawText(x - 14, top - 5, f"{nm_tick}")

        # Horizontal grid (intensity)
        for pct in range(0, 105, 10):
            y = int(bottom - (pct / 100.0) * plot_h)
            painter.setPen(pen_100 if pct == 100 or pct == 0 else pen_10)
            painter.drawLine(left, y, right, y)
            if pct > 0:
                val = self._val_from_y(y, top, bottom, max_val, log_scale)
                txt = f"{val/1000:.0f}K" if val >= 10000 else f"{val:.0f}"
                painter.setPen(QColor(0, 0, 0))
                painter.drawText(4, y + 4, txt)

        # Fill with visible colours (rainbow under the curve)
        if self.use_colors:
            for px in range(left, right + 1):
                nm_here = nm_start + (px - left) / nm_coeff
                idx_float = (px - left) / plot_w * (len(samples) - 1)
                idx = int(np.clip(idx_float, 0, len(samples) - 1))
                val = float(samples[idx])
                y = int(self._y_from_val(val, top, bottom, max_val, log_scale))
                if y < bottom:
                    r, g, b = wavelength_to_rgb(nm_here)
                    if r > 0 or g > 0 or b > 0:
                        painter.setPen(QColor(r, g, b))
                        painter.drawLine(px, bottom - 1, px, y)

        # Spectrum line outline
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        pen_curve = QPen(QColor(0, 70, 0), 1.5)
        painter.setPen(pen_curve)
        poly = QPolygonF()
        step = max(1, len(samples) // plot_w)
        for i in range(0, len(samples), step):
            x = left + (nms[i] - nm_start) * nm_coeff
            y = self._y_from_val(float(samples[i]), top, bottom, max_val, log_scale)
            poly.append(QPointF(x, y))
        painter.drawPolyline(poly)

        # Optional theoretical blackbody (Planck) curve for the computed CCT
        if self.show_planck_curve and self.engine.cct_kelvin is not None:
            planck_vals = planck_relative_spectrum(nms, self.engine.cct_kelvin, max_val)
            pen_planck = QPen(QColor(180, 40, 140, 200), 1.8, Qt.PenStyle.DashLine)
            painter.setPen(pen_planck)
            poly_p = QPolygonF()
            for i in range(0, len(nms), max(1, step * 2)):
                x = left + (nms[i] - nm_start) * nm_coeff
                y = self._y_from_val(float(planck_vals[i]), top, bottom, max_val, log_scale)
                poly_p.append(QPointF(x, np.clip(y, top, bottom)))
            painter.drawPolyline(poly_p)

        # Optional human photopic sensitivity curve V(lambda) (used for lux calculation)
        if self.show_vlambda_curve:
            _, v_lambda, _ = cie1931_cmf(nms)
            v_scaled = v_lambda * (max_val * 0.85)
            pen_v = QPen(QColor(0, 130, 210, 210), 1.8, Qt.PenStyle.DotLine)
            painter.setPen(pen_v)
            poly_v = QPolygonF()
            for i in range(0, len(nms), max(1, step * 2)):
                x = left + (nms[i] - nm_start) * nm_coeff
                y = self._y_from_val(float(v_scaled[i]), top, bottom, max_val, log_scale)
                poly_v.append(QPointF(x, np.clip(y, top, bottom)))
            painter.drawPolyline(poly_v)

        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)

        # Draw calibration points (TrimScale)
        if self.show_trim_points:
            for b_norm, nm_cal in zip(self.engine.calib_bins_norm, self.engine.calib_nm):
                x = int(left + (nm_cal - nm_start) * nm_coeff)
                if left <= x <= right:
                    painter.setPen(QPen(QColor(220, 0, 0), 1, Qt.PenStyle.DashLine))
                    painter.drawLine(x, top, x, bottom)
                    painter.fillRect(x - 18, top + 2, 36, 16, QColor(255, 255, 120))
                    painter.setPen(QColor(180, 0, 0))
                    painter.drawRect(x - 18, top + 2, 36, 16)
                    painter.setPen(QColor(0, 0, 0))
                    painter.drawText(x - 15, top + 14, f"{nm_cal:.0f}")

        # Peak detection and marking (Theremino MarkAllPeaks)
        if self.show_peaks or self.show_dips:
            delta = max(4, (20 * len(samples)) // max(1, plot_w))
            for i in range(delta, len(samples) - delta):
                v = float(samples[i])
                if self.show_peaks and v * 25.0 > max_val:
                    window = samples[i - delta:i + delta + 1]
                    if v >= float(np.max(window)) and v > float(samples[i - 1]):
                        x = int(left + (nms[i] - nm_start) * nm_coeff)
                        y = int(self._y_from_val(v, top, bottom, max_val, log_scale))
                        lbl_y = max(top + 5, y - 22)
                        painter.setPen(QColor(220, 0, 0))
                        painter.drawLine(x, y, x, bottom)
                        painter.fillRect(x - 18, lbl_y, 38, 15, QColor(255, 255, 140))
                        painter.setPen(QColor(0, 140, 0))
                        painter.drawRect(x - 18, lbl_y, 38, 15)
                        painter.setPen(QColor(0, 0, 0))
                        painter.drawText(x - 15, lbl_y + 12, f"{nms[i]:.1f}")

        # Axis frame
        painter.setPen(pen_100)
        painter.drawRect(left, top, plot_w, plot_h)

        # Info line under the plot (cursor + CCT + Lux + Max)
        font_bold = QFont("Sans", 9, QFont.Weight.Bold)
        painter.setFont(font_bold)
        painter.setPen(QColor(10, 30, 90))
        cct_txt = f"CCT: {self.engine.cct_kelvin:.0f} K" if self.engine.cct_kelvin else "CCT: --- K"
        lux_txt = f"E: {self.engine.lux_value:.1f} lx"
        if self.cursor_pos and left <= self.cursor_pos[0] <= right and top <= self.cursor_pos[1] <= bottom:
            cx, cy = self.cursor_pos
            cur_nm = nm_start + (cx - left) / nm_coeff
            cur_val = self._val_from_y(cy, top, bottom, max_val, log_scale)
            painter.setPen(QPen(QColor(100, 100, 100), 1, Qt.PenStyle.DotLine))
            painter.drawLine(cx, top, cx, bottom)
            painter.drawLine(left, cy, right, cy)
            painter.setPen(QColor(10, 30, 90))
            info_str = f"Cursor: {cur_val:.0f} @ {cur_nm:.1f} nm   |   Max: {self.engine.max_peak_val:.0f} @ {self.engine.max_peak_nm:.1f} nm   |   {cct_txt}   |   {lux_txt}"
        else:
            info_str = f"Max: {self.engine.max_peak_val:.0f} @ {self.engine.max_peak_nm:.1f} nm   |   {cct_txt}   |   {lux_txt}   |   Area: {self.engine.peak_area:.0f}   |   Range: {nm_start:.0f} – {nm_end:.0f} nm"
        painter.drawText(left, h - 6, info_str)


class CameraPreviewWidget(QWidget):
    """Widget showing the camera image strip with the ROI frame (like PBox_Samples in Theremino)."""

    def __init__(self, engine: SpectrumEngine, parent=None):
        super().__init__(parent)
        self.engine = engine
        self.setMinimumHeight(80)
        self.setMaximumHeight(130)

    def paintEvent(self, event):
        painter = QPainter(self)
        w, h = self.width(), self.height()
        painter.fillRect(0, 0, w, h, QColor(20, 20, 20))

        rgb = self.engine.roi_preview_rgb
        if rgb is None or rgb.size == 0:
            painter.setPen(QColor(200, 200, 200))
            painter.drawText(15, h // 2, "No camera image (check the connection / /dev/video* permissions)")
            return

        ih, iw, _ = rgb.shape
        qimg = QImage(rgb.data, iw, ih, 3 * iw, QImage.Format.Format_RGB888)
        scaled_pix = QPixmap.fromImage(qimg).scaled(w, h, Qt.AspectRatioMode.IgnoreAspectRatio, Qt.TransformationMode.FastTransformation)
        painter.drawPixmap(0, 0, scaled_pix)

        # Orange frame of the X slice (StartX .. EndX) as in Theremino
        x0 = int(self.engine.src_x0 / max(1, iw) * w)
        dx = int(self.engine.src_dx / max(1, iw) * w)
        pen = QPen(QColor(230, 140, 0), 3)
        painter.setPen(pen)
        painter.drawRect(x0 + 1, 2, max(2, dx - 3), h - 4)


class CalibrationDialog(QDialog):
    """Calibration point editor (equivalent of Trim Scale in Theremino)."""

    def __init__(self, engine: SpectrumEngine, parent=None):
        super().__init__(parent)
        self.engine = engine
        self.setWindowTitle("Wavelength calibration (nm)")
        self.resize(420, 300)

        layout = QVBoxLayout(self)
        self.table = QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["Position (pixel 0..W)", "Wavelength (nm)"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        self._load_points()

        btn_layout = QHBoxLayout()
        btn_add = QPushButton("Add point")
        btn_del = QPushButton("Remove selected")
        btn_apply = QPushButton("Apply")
        btn_close = QPushButton("Close")

        btn_add.clicked.connect(self._add_row)
        btn_del.clicked.connect(self._del_row)
        btn_apply.clicked.connect(self._apply)
        btn_close.clicked.connect(self.accept)

        btn_layout.addWidget(btn_add)
        btn_layout.addWidget(btn_del)
        btn_layout.addWidget(btn_apply)
        btn_layout.addWidget(btn_close)
        layout.addLayout(btn_layout)

    def _load_points(self):
        self.table.setRowCount(0)
        for b_norm, nm in zip(self.engine.calib_bins_norm, self.engine.calib_nm):
            row = self.table.rowCount()
            self.table.insertRow(row)
            px = b_norm * self.engine.num_samples
            self.table.setItem(row, 0, QTableWidgetItem(f"{px:.1f}"))
            self.table.setItem(row, 1, QTableWidgetItem(f"{nm:.2f}"))

    def _add_row(self):
        row = self.table.rowCount()
        self.table.insertRow(row)
        self.table.setItem(row, 0, QTableWidgetItem("960.0"))
        self.table.setItem(row, 1, QTableWidgetItem("546.0"))

    def _del_row(self):
        if self.table.rowCount() <= 2:
            QMessageBox.warning(self, "Calibration", "At least 2 calibration points are required.")
            return
        row = self.table.currentRow()
        if row >= 0:
            self.table.removeRow(row)

    def _apply(self):
        points = []
        for r in range(self.table.rowCount()):
            try:
                px = float(self.table.item(r, 0).text().replace(",", "."))
                nm = float(self.table.item(r, 1).text().replace(",", "."))
                b_norm = np.clip(px / max(1, self.engine.num_samples), 0.0, 1.0)
                points.append((b_norm, nm))
            except Exception:
                pass
        points.sort(key=lambda x: x[0])
        if len(points) >= 2:
            self.engine.calib_bins_norm = [p[0] for p in points]
            self.engine.calib_nm = [p[1] for p in points]
            self._load_points()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Theremino Spectrometer V5.0 - Native Linux Edition (Python + CCT + Lux)")
        self.resize(1240, 760)

        self.base_dir = Path(__file__).resolve().parent
        self.lux_config_path = str(self.base_dir / "calib" / "lux_calib.json")
        self.engine = SpectrumEngine()

        # Load the default Theremino calibration file if it exists
        default_calib = self.base_dir / "calib" / "Calib1.txt"
        if default_calib.exists():
            self.engine.load_calibration_file(str(default_calib))

        # Load the LaoCorrections.txt sensor sensitivity correction file if it exists
        default_irrad = self.base_dir / "irradiance" / "LaoCorrections.txt"
        if default_irrad.exists():
            self.engine.load_irradiance_file(str(default_irrad))

        # Load the saved lux meter calibration if it exists
        self.engine.load_lux_calibration(self.lux_config_path)

        self.cap: Optional[cv2.VideoCapture] = None
        self.fps_counter: float = 0.0
        self.last_frame_time: float = time.time()

        self._build_ui()

        self.timer = QTimer(self)
        self.timer.timeout.connect(self._on_timer)
        self.timer.start(30)

        # Try to open the camera automatically
        self._open_camera(0)

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        # Top toolbar and camera selection
        top_bar = QHBoxLayout()

        top_bar.addWidget(QLabel("Camera:"))
        self.cmb_cam = QComboBox()
        sys_name = platform.system()
        for i in range(4):
            label = f"/dev/video{i}" if sys_name == "Linux" else f"Camera #{i}"
            self.cmb_cam.addItem(label, i)
        self.cmb_cam.addItem("Demo mode (Hg/CFL simulation)", -1)
        self.cmb_cam.currentIndexChanged.connect(self._on_cam_changed)
        top_bar.addWidget(self.cmb_cam)

        top_bar.addWidget(QLabel("Resolution:"))
        self.cmb_res = QComboBox()
        self.cmb_res.addItems(["1920x1080", "1280x720", "640x480"])
        self.cmb_res.currentIndexChanged.connect(self._on_res_changed)
        top_bar.addWidget(self.cmb_res)

        self.lbl_fps = QLabel("0.0 FPS")
        self.lbl_fps.setStyleSheet("font-weight: bold; color: #005500; padding: 0 8px;")
        top_bar.addWidget(self.lbl_fps)

        top_bar.addStretch()

        self.btn_bg = QPushButton("Background")
        self.btn_bg.setCheckable(True)
        self.btn_bg.clicked.connect(self._toggle_background)
        top_bar.addWidget(self.btn_bg)

        self.btn_ref = QPushButton("Reference")
        self.btn_ref.setCheckable(True)
        self.btn_ref.clicked.connect(self._toggle_reference)
        top_bar.addWidget(self.btn_ref)

        btn_load_cal = QPushButton("Load Calib...")
        btn_load_cal.clicked.connect(self._load_calib_dialog)
        top_bar.addWidget(btn_load_cal)

        btn_edit_cal = QPushButton("Edit calibration...")
        btn_edit_cal.clicked.connect(self._edit_calib_dialog)
        top_bar.addWidget(btn_edit_cal)

        btn_save_csv = QPushButton("Save spectrum CSV")
        btn_save_csv.clicked.connect(self._save_spectrum_csv)
        top_bar.addWidget(btn_save_csv)

        main_layout.addLayout(top_bar)

        # Camera strip preview (ROI)
        self.preview_widget = CameraPreviewWidget(self.engine)
        main_layout.addWidget(self.preview_widget)

        # Middle parameter panel
        params_layout = QHBoxLayout()

        # Group 1: Sensor area (ROI)
        grp_roi = QGroupBox("Image slice (ROI 0..1000)")
        grid_roi = QGridLayout(grp_roi)

        self.chk_flip_h = QCheckBox("Flip X (FlipH)")
        self.chk_flip_h.setChecked(self.engine.flip_h)
        self.chk_flip_h.toggled.connect(lambda v: setattr(self.engine, "flip_h", v))
        grid_roi.addWidget(self.chk_flip_h, 0, 0, 1, 2)

        self.spn_start_x = self._make_spin(0, 999, self.engine.start_x, lambda v: setattr(self.engine, "start_x", v))
        self.spn_end_x = self._make_spin(1, 1000, self.engine.end_x, lambda v: setattr(self.engine, "end_x", v))
        self.spn_start_y = self._make_spin(0, 999, self.engine.start_y, lambda v: setattr(self.engine, "start_y", v))
        self.spn_size_y = self._make_spin(1, 1000, self.engine.size_y, lambda v: setattr(self.engine, "size_y", v))

        grid_roi.addWidget(QLabel("Start X:"), 1, 0)
        grid_roi.addWidget(self.spn_start_x, 1, 1)
        grid_roi.addWidget(QLabel("End X:"), 1, 2)
        grid_roi.addWidget(self.spn_end_x, 1, 3)
        grid_roi.addWidget(QLabel("Start Y:"), 2, 0)
        grid_roi.addWidget(self.spn_start_y, 2, 1)
        grid_roi.addWidget(QLabel("Size Y:"), 2, 2)
        grid_roi.addWidget(self.spn_size_y, 2, 3)
        params_layout.addWidget(grp_roi)

        # Group 2: Filters (Spatial, Rising/Falling Speed, LogScale)
        grp_flt = QGroupBox("Filters and scale (Theremino DSP)")
        grid_flt = QGridLayout(grp_flt)

        self.spn_spatial = self._make_spin(0, 10, self.engine.spatial_averaging, lambda v: setattr(self.engine, "spatial_averaging", v))
        self.spn_rise = self._make_spin(1, 100, self.engine.rising_speed, lambda v: setattr(self.engine, "rising_speed", v))
        self.spn_fall = self._make_spin(1, 100, self.engine.falling_speed, lambda v: setattr(self.engine, "falling_speed", v))
        self.spn_log = QDoubleSpinBox()
        self.spn_log.setRange(-10.0, 10.0)
        self.spn_log.setSingleStep(0.5)
        self.spn_log.setValue(self.engine.log_scale_param)
        self.spn_log.valueChanged.connect(lambda v: setattr(self.engine, "log_scale_param", v))

        grid_flt.addWidget(QLabel("Smoothing (0-10):"), 0, 0)
        grid_flt.addWidget(self.spn_spatial, 0, 1)
        grid_flt.addWidget(QLabel("Log scale:"), 0, 2)
        grid_flt.addWidget(self.spn_log, 0, 3)
        grid_flt.addWidget(QLabel("Rising Speed %:"), 1, 0)
        grid_flt.addWidget(self.spn_rise, 1, 1)
        grid_flt.addWidget(QLabel("Falling Speed %:"), 1, 2)
        grid_flt.addWidget(self.spn_fall, 1, 3)
        params_layout.addWidget(grp_flt)

        # Group 3: V4L2 camera controls (Exposure / Brightness / White balance)
        grp_cam = QGroupBox("V4L2 sensor controls")
        grid_cam = QGridLayout(grp_cam)

        self.chk_auto_exp = QCheckBox("Auto Exposure")
        self.chk_auto_exp.setChecked(True)
        self.chk_auto_exp.toggled.connect(self._on_auto_exp_toggled)
        grid_cam.addWidget(self.chk_auto_exp, 0, 0, 1, 3)

        is_linux = platform.system() == "Linux"
        self.lbl_exp_val = QLabel("15.7 ms (157)" if is_linux else "31.2 ms (-5)")
        self.sld_exp = QSlider(Qt.Orientation.Horizontal)
        if is_linux:
            self.sld_exp.setRange(1, 5000)
            self.sld_exp.setValue(157)
        else:
            self.sld_exp.setRange(-13, -1)
            self.sld_exp.setValue(-5)
        self.sld_exp.setEnabled(False)
        self.sld_exp.valueChanged.connect(self._on_exp_changed)
        grid_cam.addWidget(QLabel("Exp.:"), 1, 0)
        grid_cam.addWidget(self.sld_exp, 1, 1)
        grid_cam.addWidget(self.lbl_exp_val, 1, 2)

        self.sld_gain = QSlider(Qt.Orientation.Horizontal)
        self.sld_gain.setRange(-64, 64)
        self.sld_gain.setValue(0)
        self.sld_gain.valueChanged.connect(self._on_gain_changed)
        grid_cam.addWidget(QLabel("Brightness:"), 2, 0)
        grid_cam.addWidget(self.sld_gain, 2, 1, 1, 2)

        params_layout.addWidget(grp_cam)

        # Group 4: Colour temperature CCT (Kelvin) & CIE 1931
        grp_cct = QGroupBox("Colour temperature (CIE 1931)")
        vbox_cct = QVBoxLayout(grp_cct)

        top_cct_row = QHBoxLayout()
        self.lbl_cct_big = QLabel("--- K")
        self.lbl_cct_big.setStyleSheet("font-size: 17px; font-weight: bold; color: #103070;")
        top_cct_row.addWidget(self.lbl_cct_big)

        self.frm_swatch = QFrame()
        self.frm_swatch.setFixedSize(34, 22)
        self.frm_swatch.setStyleSheet("background-color: #000000; border: 1px solid #444; border-radius: 4px;")
        top_cct_row.addWidget(self.frm_swatch)
        vbox_cct.addLayout(top_cct_row)

        self.lbl_cie_xy = QLabel("CIE x: ---   y: ---   Duv: ---")
        self.lbl_cie_xy.setStyleSheet("font-size: 11px; color: #333;")
        vbox_cct.addWidget(self.lbl_cie_xy)

        irrad_row = QHBoxLayout()
        self.chk_irrad = QCheckBox("Correction (LaoCorrections)")
        self.chk_irrad.setChecked(False)
        self.chk_irrad.toggled.connect(lambda v: setattr(self.engine, "irradiance_enabled", v))
        irrad_row.addWidget(self.chk_irrad)

        btn_irrad_file = QPushButton("...")
        btn_irrad_file.setFixedWidth(26)
        btn_irrad_file.clicked.connect(self._load_irrad_dialog)
        irrad_row.addWidget(btn_irrad_file)
        vbox_cct.addLayout(irrad_row)

        params_layout.addWidget(grp_cct)

        # Group 5: Lux meter (illuminance in lx)
        grp_lux = QGroupBox("Illuminance (Lux meter)")
        vbox_lux = QVBoxLayout(grp_lux)

        self.lbl_lux_big = QLabel("0.0 lx")
        self.lbl_lux_big.setStyleSheet("font-size: 18px; font-weight: bold; color: #006622;")
        vbox_lux.addWidget(self.lbl_lux_big)

        self.lbl_lux_status = QLabel("Estimated (Auto Exp. mode)")
        self.lbl_lux_status.setStyleSheet("font-size: 11px; color: #555;")
        vbox_lux.addWidget(self.lbl_lux_status)

        btn_cal_lux = QPushButton("Calibrate with reference (lx)...")
        btn_cal_lux.setToolTip(
            "Lets you enter a known value in lux (e.g. from a reference lux meter)\n"
            "to accurately scale the attenuation of the slit and diffraction grating."
        )
        btn_cal_lux.clicked.connect(self._calibrate_lux_dialog)
        vbox_lux.addWidget(btn_cal_lux)

        params_layout.addWidget(grp_lux)

        # Group 6: Plot options
        grp_view = QGroupBox("Plot view")
        vbox_view = QVBoxLayout(grp_view)

        self.chk_colors = QCheckBox("Spectrum colours (Colors)")
        self.chk_colors.setChecked(True)
        self.chk_colors.toggled.connect(lambda v: setattr(self.graph_widget, "use_colors", v))

        self.chk_peaks = QCheckBox("Mark peaks (Peaks)")
        self.chk_peaks.setChecked(True)
        self.chk_peaks.toggled.connect(lambda v: setattr(self.graph_widget, "show_peaks", v))

        self.chk_planck = QCheckBox("Planck curve (CCT)")
        self.chk_planck.setChecked(True)
        self.chk_planck.toggled.connect(lambda v: setattr(self.graph_widget, "show_planck_curve", v))

        self.chk_vlambda = QCheckBox("Eye sensitivity V(λ) (Lux)")
        self.chk_vlambda.setChecked(False)
        self.chk_vlambda.toggled.connect(lambda v: setattr(self.graph_widget, "show_vlambda_curve", v))

        self.chk_trim = QCheckBox("Calibration points (Trim)")
        self.chk_trim.setChecked(False)
        self.chk_trim.toggled.connect(lambda v: setattr(self.graph_widget, "show_trim_points", v))

        vbox_view.addWidget(self.chk_colors)
        vbox_view.addWidget(self.chk_peaks)
        vbox_view.addWidget(self.chk_planck)
        vbox_view.addWidget(self.chk_vlambda)
        vbox_view.addWidget(self.chk_trim)
        params_layout.addWidget(grp_view)

        main_layout.addLayout(params_layout)

        # Main spectrum plot
        self.graph_widget = SpectrumGraphWidget(self.engine)
        main_layout.addWidget(self.graph_widget, stretch=1)

    def _make_spin(self, min_v: int, max_v: int, val: int, callback) -> QSpinBox:
        spn = QSpinBox()
        spn.setRange(min_v, max_v)
        spn.setValue(val)
        spn.valueChanged.connect(callback)
        return spn

    def _open_camera(self, dev_index: int):
        if self.cap is not None:
            self.cap.release()
            self.cap = None

        if dev_index < 0:
            return

        sys_name = platform.system()
        if sys_name == "Linux":
            cap = cv2.VideoCapture(dev_index, cv2.CAP_V4L2)
        elif sys_name == "Darwin":
            cap = cv2.VideoCapture(dev_index, cv2.CAP_AVFOUNDATION)
        elif sys_name == "Windows":
            cap = cv2.VideoCapture(dev_index, cv2.CAP_DSHOW)
        else:
            cap = cv2.VideoCapture(dev_index)

        if not cap.isOpened():
            cap = cv2.VideoCapture(dev_index)

        if cap.isOpened():
            w, h = map(int, self.cmb_res.currentText().split("x"))
            cap.set(cv2.CAP_PROP_FOURCC, cv2.VideoWriter_fourcc(*"MJPG"))
            cap.set(cv2.CAP_PROP_FRAME_WIDTH, w)
            cap.set(cv2.CAP_PROP_FRAME_HEIGHT, h)
            cap.set(cv2.CAP_PROP_FPS, 15)
            self.cap = cap
        else:
            print(f"Could not open camera with index {dev_index}. Switching to Demo mode.")
            self.cmb_cam.blockSignals(True)
            self.cmb_cam.setCurrentIndex(self.cmb_cam.count() - 1)
            self.cmb_cam.blockSignals(False)

    def _on_cam_changed(self, idx: int):
        dev_id = self.cmb_cam.currentData()
        self._open_camera(dev_id)

    def _on_res_changed(self, idx: int):
        dev_id = self.cmb_cam.currentData()
        if dev_id is not None and dev_id >= 0:
            self._open_camera(dev_id)

    def _on_auto_exp_toggled(self, checked: bool):
        self.sld_exp.setEnabled(not checked)
        self.engine.auto_exposure_active = checked
        if self.cap is not None and self.cap.isOpened():
            self.cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, 3.0 if checked else 1.0)
            if not checked:
                self.cap.set(cv2.CAP_PROP_EXPOSURE, float(self.sld_exp.value()))

    def _on_exp_changed(self, val: int):
        self.engine.exposure_step = val
        is_linux = platform.system() == "Linux"
        sec = decode_exposure_seconds(val, is_linux)
        if sec >= 0.1:
            self.lbl_exp_val.setText(f"{sec*1000:.0f} ms ({val})")
        elif sec >= 0.001:
            self.lbl_exp_val.setText(f"{sec*1000:.1f} ms ({val})")
        else:
            self.lbl_exp_val.setText(f"{sec*1e6:.0f} µs ({val})")

        if self.cap is not None and self.cap.isOpened():
            self.cap.set(cv2.CAP_PROP_EXPOSURE, float(val))

    def _on_gain_changed(self, val: int):
        self.engine.camera_gain = val
        if self.cap is not None and self.cap.isOpened():
            self.cap.set(cv2.CAP_PROP_BRIGHTNESS, float(val))

    def _calibrate_lux_dialog(self):
        if self.engine.raw_photopic_rate <= 1e-3:
            QMessageBox.warning(self, "Lux meter calibration", "The light signal is too weak for calibration.")
            return

        val, ok = QInputDialog.getDouble(
            self,
            "Lux meter calibration (lx)",
            "Preferably disable 'Auto Exposure' first to lock the shutter time.\n\n"
            "Enter the current illuminance falling on the sensor in lux (lx)\n"
            "(e.g. measured with a reference lux meter placed next to the spectrometer):",
            value=max(1.0, round(self.engine.lux_value, 1)),
            min=0.1,
            max=200000.0,
            decimals=1,
        )
        if ok and val > 0:
            self.engine.lux_calib_factor = val / self.engine.raw_photopic_rate
            self.engine.lux_calibrated_by_user = True
            self.engine.lux_value = val
            self.engine.save_lux_calibration(self.lux_config_path)
            QMessageBox.information(
                self,
                "Calibration saved",
                f"Lux meter calibration factor saved to:\n{self.lux_config_path}",
            )

    def _generate_demo_frame(self) -> np.ndarray:
        """Generates a synthetic CFL spectrum image (Hg/Eu peaks: 405, 436, 546, 611 nm)."""
        w, h = 1920, 400
        frame = np.zeros((h, w, 3), dtype=np.uint8)
        nms = self.engine.bin_to_nm(np.arange(w, dtype=np.float64))

        peaks = [
            (404.7, 90.0, 2.5),
            (435.8, 190.0, 2.8),
            (487.7, 55.0, 6.0),
            (546.1, 240.0, 2.5),
            (580.0, 45.0, 8.0),
            (611.6, 210.0, 3.2),
        ]
        spec = np.full(w, 4.0, dtype=np.float32)
        for center_nm, amp, sigma in peaks:
            spec += amp * np.exp(-0.5 * ((nms - center_nm) / sigma) ** 2)

        noise = np.random.normal(0, 1.8, size=w).astype(np.float32)
        spec = np.clip(spec + noise, 0, 255)

        row_bgr = np.zeros((1, w, 3), dtype=np.uint8)
        for i in range(w):
            r, g, b = wavelength_to_rgb(float(nms[i]))
            scale = spec[i] / 255.0
            row_bgr[0, i] = [int(b * scale), int(g * scale), int(r * scale)]

        y_center = h // 2
        band_h = 80
        frame[y_center - band_h // 2:y_center + band_h // 2, :] = row_bgr
        if self.engine.flip_h:
            frame = cv2.flip(frame, 1)
        return frame

    def _on_timer(self):
        frame = None
        if self.cap is not None and self.cap.isOpened():
            ret, f = self.cap.read()
            if ret and f is not None:
                frame = f

        if frame is None:
            frame = self._generate_demo_frame()

        now = time.time()
        dt = max(1e-4, now - self.last_frame_time)
        self.last_frame_time = now
        self.fps_counter = 0.85 * self.fps_counter + 0.15 * (1.0 / dt)
        self.lbl_fps.setText(f"{self.fps_counter:.1f} FPS")

        self.engine.process_frame(frame)

        # Update the CCT / CIE 1931 panel
        if self.engine.cct_kelvin is not None:
            cct = self.engine.cct_kelvin
            if cct < 3300:
                desc = "warm"
            elif cct < 5300:
                desc = "neutral"
            else:
                desc = "cool / daylight"
            self.lbl_cct_big.setText(f"{cct:.0f} K ({desc})")
            self.lbl_cie_xy.setText(
                f"CIE x: {self.engine.cie_x:.4f}   y: {self.engine.cie_y:.4f}   Duv: {self.engine.duv:+.4f}"
            )
            r, g, b = self.engine.swatch_rgb
            self.frm_swatch.setStyleSheet(
                f"background-color: rgb({r}, {g}, {b}); border: 1px solid #444; border-radius: 4px;"
            )
        else:
            self.lbl_cct_big.setText("--- K (weak signal)")
            self.lbl_cie_xy.setText("CIE x: ---   y: ---   Duv: ---")
            self.frm_swatch.setStyleSheet("background-color: #222; border: 1px solid #444; border-radius: 4px;")

        # Update the lux meter panel (lx)
        lux = self.engine.lux_value
        if lux >= 1000.0:
            self.lbl_lux_big.setText(f"{lux:.0f} lx")
        elif lux >= 10.0:
            self.lbl_lux_big.setText(f"{lux:.1f} lx")
        else:
            self.lbl_lux_big.setText(f"{lux:.2f} lx")

        if self.engine.sensor_saturated:
            self.lbl_lux_big.setStyleSheet("font-size: 18px; font-weight: bold; color: #cc0000;")
            self.lbl_lux_status.setText("WARNING: Sensor saturated! Reduce exposure.")
            self.lbl_lux_status.setStyleSheet("font-size: 11px; font-weight: bold; color: #cc0000;")
        else:
            self.lbl_lux_big.setStyleSheet("font-size: 18px; font-weight: bold; color: #006622;")
            cal_tag = "Calibrated" if self.engine.lux_calibrated_by_user else "Estimated"
            mode_tag = "Auto Exp." if self.engine.auto_exposure_active else "Fixed Exp."
            self.lbl_lux_status.setText(f"{cal_tag} ({mode_tag})")
            self.lbl_lux_status.setStyleSheet("font-size: 11px; color: #444;")

        self.preview_widget.update()
        self.graph_widget.update()

    def _toggle_background(self, checked: bool):
        if checked:
            self.engine.set_background()
        else:
            self.engine.reset_background()

    def _toggle_reference(self, checked: bool):
        if checked:
            self.engine.set_reference()
        else:
            self.engine.reset_reference()

    def _load_calib_dialog(self):
        calib_dir = str(self.base_dir / "calib")
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load Theremino calibration file", calib_dir, "Calibration Files (*.txt);;All Files (*)"
        )
        if filepath:
            if self.engine.load_calibration_file(filepath):
                QMessageBox.information(self, "Calibration", f"Calibration loaded from:\n{Path(filepath).name}")
            else:
                QMessageBox.warning(self, "Error", "Invalid calibration file format.")

    def _load_irrad_dialog(self):
        irrad_dir = str(self.base_dir / "irradiance")
        filepath, _ = QFileDialog.getOpenFileName(
            self, "Load IrradianceCoeffs correction file", irrad_dir, "Irradiance Files (*.txt);;All Files (*)"
        )
        if filepath:
            if self.engine.load_irradiance_file(filepath):
                self.chk_irrad.setChecked(True)
                self.chk_irrad.setText(f"Correction ({Path(filepath).stem[:14]})")
                QMessageBox.information(self, "Irradiance", f"Sensitivity correction loaded:\n{Path(filepath).name}")
            else:
                QMessageBox.warning(self, "Error", "Invalid IrradianceCoeffs file format.")

    def _edit_calib_dialog(self):
        dlg = CalibrationDialog(self.engine, self)
        dlg.exec()

    def _save_spectrum_csv(self):
        filepath, _ = QFileDialog.getSaveFileName(
            self, "Save spectrum CSV (Theremino format)", str(self.base_dir / "spectrum.csv"), "CSV Files (*.csv *.txt)"
        )
        if not filepath:
            return
        try:
            with open(filepath, "w", encoding="utf-8") as f:
                f.write("Sensor;WebCam\n")
                f.write("-----------------------\n")
                f.write(f"Rec.Samples;{self.engine.num_samples}\n")
                f.write(f"FramesPerSec;{self.fps_counter:.1f} fps\n")
                f.write(f"Spatial avg.;{self.engine.spatial_averaging}\n")
                f.write(f"RisingSpeed;{self.engine.rising_speed}\n")
                f.write(f"FallingSpeed;{self.engine.falling_speed}\n")
                f.write(f"NanometersMin;{self.engine.visible_nm[0]:.2f}\n")
                f.write(f"NanometersMax;{self.engine.visible_nm[-1]:.2f}\n")
                f.write(f"Peak Area;{self.engine.peak_area:.1f}\n")
                f.write(f"Illuminance_Lux;{self.engine.lux_value:.2f} lx\n")
                if self.engine.cct_kelvin is not None:
                    f.write(f"CCT_Kelvin;{self.engine.cct_kelvin:.0f} K\n")
                    f.write(f"CIE1931_xy;{self.engine.cie_x:.4f}, {self.engine.cie_y:.4f}\n")
                    f.write(f"Duv;{self.engine.duv:+.4f}\n")
                f.write("-----------------------\n")
                f.write("Nanometers;Intensity\n")
                f.write("-----------------------\n")
                for nm, val in zip(self.engine.visible_nm, self.engine.visible_samples):
                    f.write(f"{nm:.2f};{val:12.1f}\n")
            QMessageBox.information(self, "Saved", f"Spectrum saved to file:\n{filepath}")
        except Exception as e:
            QMessageBox.critical(self, "Save error", str(e))


def main():
    app = QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
