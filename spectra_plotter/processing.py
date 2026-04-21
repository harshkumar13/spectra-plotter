"""Spectral data processing: cleaning, smoothing, binning, cosmic ray removal."""

from typing import Sequence, Tuple

import numpy as np
from scipy.signal import savgol_filter, medfilt

from .line_catalog import TELLURIC_BANDS


def clean_spectrum(wave: np.ndarray, flux: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    valid = np.isfinite(wave) & np.isfinite(flux) & (wave > 0)
    wave, flux = wave[valid], flux[valid]
    if wave.size == 0:
        return wave, flux
    order = np.argsort(wave)
    wave, flux = wave[order], flux[order]
    _, uniq = np.unique(wave, return_index=True)
    return wave[np.sort(uniq)], flux[np.sort(uniq)]


def apply_telluric_mask(
    wave_obs: np.ndarray, flux: np.ndarray, band_keys: Sequence[str],
) -> np.ndarray:
    if not band_keys:
        return flux
    masked = flux.astype(float).copy()
    for key in band_keys:
        if key not in TELLURIC_BANDS:
            continue
        _, (lo, hi) = TELLURIC_BANDS[key]
        masked[(wave_obs >= lo) & (wave_obs <= hi)] = np.nan
    return masked


def _fill_nans(y: np.ndarray) -> np.ndarray:
    y = y.astype(float).copy()
    if y.size == 0:
        return y
    x = np.arange(y.size)
    good = np.isfinite(y)
    if good.all() or not good.any():
        return y if good.any() else np.zeros_like(y)
    y[~good] = np.interp(x[~good], x[good], y[good])
    return y


def cosmic_ray_clip(
    flux: np.ndarray, kernel: int = 9, sigma: float = 8.0,
) -> np.ndarray:
    y = flux.astype(float).copy()
    nan_mask = ~np.isfinite(y)
    y_filled = _fill_nans(y)
    kernel = max(3, kernel | 1)
    baseline = medfilt(y_filled, kernel_size=kernel)
    resid = y_filled - baseline
    mad = np.median(np.abs(resid - np.median(resid)))
    if mad == 0 or not np.isfinite(mad):
        return flux
    spike = np.abs(resid) > sigma * 1.4826 * mad
    y_out = y_filled.copy()
    y_out[spike] = baseline[spike]
    y_out[nan_mask] = np.nan
    return y_out


def apply_savgol(flux: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    y = flux.astype(float).copy()
    nan_mask = ~np.isfinite(y)
    y_filled = _fill_nans(y)
    n = y_filled.size
    if n < 5 or window < 3:
        return y
    window = max(3, int(window) | 1)
    window = min(window, n if n % 2 == 1 else n - 1)
    polyorder = min(max(1, int(polyorder)), window - 1)
    y_smooth = savgol_filter(y_filled, window_length=window, polyorder=polyorder, mode="interp")
    y_smooth[nan_mask] = np.nan
    return y_smooth


def apply_binning(
    wave: np.ndarray, flux: np.ndarray, bin_size: int,
) -> Tuple[np.ndarray, np.ndarray]:
    bin_size = int(bin_size)
    if bin_size <= 1:
        return wave, flux
    m = (wave.size // bin_size) * bin_size
    if m < bin_size:
        return wave, flux
    return (
        np.nanmean(wave[:m].reshape(-1, bin_size), axis=1),
        np.nanmean(flux[:m].reshape(-1, bin_size), axis=1),
    )


def apply_gaussian_smooth(flux: np.ndarray, sigma: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter1d
    y = flux.astype(float).copy()
    nan_mask = ~np.isfinite(y)
    y_filled = _fill_nans(y)
    if y_filled.size < 5 or sigma <= 0:
        return y
    y_smooth = gaussian_filter1d(y_filled, sigma=sigma)
    y_smooth[nan_mask] = np.nan
    return y_smooth


def compute_snr(flux: np.ndarray, window: int = 50) -> float:
    y = flux[np.isfinite(flux)]
    if y.size < window * 2:
        return 0.0
    mid = y.size // 2
    segment = y[mid - window:mid + window]
    if segment.std() == 0:
        return 0.0
    return float(segment.mean() / segment.std())


def measure_ew(
    wave: np.ndarray, flux: np.ndarray, line_center: float, width: float = 20.0,
) -> float:
    mask = (wave >= line_center - width) & (wave <= line_center + width)
    if mask.sum() < 3:
        return np.nan
    w_seg, f_seg = wave[mask], flux[mask]
    continuum = np.nanmean([f_seg[:3].mean(), f_seg[-3:].mean()])
    if continuum == 0:
        return np.nan
    ew = np.trapz(1.0 - f_seg / continuum, w_seg)
    return float(ew)
