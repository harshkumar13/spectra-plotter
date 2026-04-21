"""Spectrum file I/O: reading FITS, ASCII, and CSV spectral data.

Handles a wide variety of formats encountered in transient astronomy:
- FITS binary tables with various column naming conventions
- FITS 1D spectra with WCS wavelength solution
- FITS 2D arrays [wavelength, flux]
- Multi-extension FITS files
- Whitespace-delimited ASCII (.txt, .dat, .ascii)
- Comma-separated CSV (.csv)
- Files with comment headers (# or %)
- Files with or without column headers
"""

import os
import re
from typing import List, Tuple

import numpy as np
import pandas as pd
from astropy.io import fits

SPECTRA_EXTENSIONS = (
    ".fits", ".fit", ".fits.gz", ".fit.gz", ".fz",
    ".txt", ".dat", ".ascii", ".csv",
)


def is_spectrum_file(path: str) -> bool:
    base = os.path.basename(path)
    if base.startswith("."):
        return False
    p = path.lower()
    return any(p.endswith(ext) for ext in SPECTRA_EXTENSIONS)


def list_spectra_files(folder: str) -> List[str]:
    if not folder or not os.path.isdir(folder):
        return []
    paths = []
    for name in sorted(os.listdir(folder)):
        if name.startswith("."):
            continue
        full = os.path.join(folder, name)
        if os.path.isfile(full) and is_spectrum_file(full):
            paths.append(full)
    return paths


def _find_column(names, candidates):
    names_lower = [n.lower().strip() for n in names]
    for cand in candidates:
        for i, n in enumerate(names_lower):
            if n == cand.lower():
                return names[i]
        for i, n in enumerate(names_lower):
            if cand.lower() in n:
                return names[i]
    return None


_WAVE_CANDIDATES = [
    "OPT_WAVE", "wave", "WAVE", "wavelength", "WAVELENGTH",
    "lam", "lambda", "LAMBDA", "wl", "WL", "wavel",
    "obs_wave", "rest_wave", "angstrom", "ANGSTROM",
]

_FLUX_CANDIDATES = [
    "OPT_FLUX", "OPT_FLAM", "flux", "FLUX", "flam", "FLAM",
    "counts", "COUNTS", "intensity", "INTENSITY",
    "obs_flux", "f_lam", "f_lambda",
]


def _read_fits(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    with fits.open(file_path) as hdul:
        for ext_idx in range(len(hdul) - 1, -1, -1):
            ext = hdul[ext_idx]
            if ext.data is None:
                continue

            if hasattr(ext.data, "names") and ext.data.names:
                wave_col = _find_column(list(ext.data.names), _WAVE_CANDIDATES)
                flux_col = _find_column(list(ext.data.names), _FLUX_CANDIDATES)
                if wave_col and flux_col:
                    return (
                        np.array(ext.data[wave_col], dtype=float).flatten(),
                        np.array(ext.data[flux_col], dtype=float).flatten(),
                    )

        for ext_idx in range(len(hdul)):
            ext = hdul[ext_idx]
            if ext.data is None:
                continue

            data = np.array(ext.data, dtype=float)
            header = ext.header

            if data.ndim == 1:
                n = data.size
                crval1 = header.get("CRVAL1")
                cdelt1 = header.get("CDELT1", header.get("CD1_1"))
                crpix1 = header.get("CRPIX1", 1.0)

                if crval1 is not None and cdelt1 is not None:
                    dc_flag = header.get("DC-FLAG", 0)
                    wavelength = crval1 + (np.arange(n) - (crpix1 - 1)) * cdelt1
                    if dc_flag == 1:
                        wavelength = 10.0 ** wavelength
                    return wavelength, data.flatten()

            if data.ndim == 2 and data.shape[0] >= 2:
                return data[0].flatten(), data[1].flatten()

            if data.ndim >= 3:
                squeezed = data.squeeze()
                if squeezed.ndim == 1:
                    n = squeezed.size
                    crval1 = header.get("CRVAL1")
                    cdelt1 = header.get("CDELT1", header.get("CD1_1"))
                    crpix1 = header.get("CRPIX1", 1.0)
                    if crval1 is not None and cdelt1 is not None:
                        wavelength = crval1 + (np.arange(n) - (crpix1 - 1)) * cdelt1
                        dc_flag = header.get("DC-FLAG", 0)
                        if dc_flag == 1:
                            wavelength = 10.0 ** wavelength
                        return wavelength, squeezed
                elif squeezed.ndim == 2 and squeezed.shape[0] >= 2:
                    return squeezed[0].flatten(), squeezed[1].flatten()

    raise ValueError(f"Could not extract wavelength/flux from FITS file: {file_path}")


def _read_ascii(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    with open(file_path, "r", errors="replace") as fh:
        raw_lines = fh.readlines()

    data_lines = []
    for line in raw_lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("%"):
            continue
        if re.match(r"^[a-zA-Z]", stripped) and not any(c.isdigit() for c in stripped[:20]):
            continue
        data_lines.append(stripped)

    if not data_lines:
        raise ValueError("No data lines found in ASCII file.")

    first = data_lines[0]
    if re.match(r"^[a-zA-Z]", first):
        data_lines = data_lines[1:]

    if not data_lines:
        raise ValueError("No numeric data lines found after skipping headers.")

    is_csv = "," in data_lines[0]
    is_tab = "\t" in data_lines[0] and "," not in data_lines[0]

    rows = []
    for line in data_lines:
        if is_csv:
            parts = line.split(",")
        elif is_tab:
            parts = line.split("\t")
        else:
            parts = line.split()

        nums = []
        for p in parts:
            p = p.strip()
            if not p:
                continue
            try:
                nums.append(float(p))
            except ValueError:
                continue
        if len(nums) >= 2:
            rows.append(nums[:2])

    if not rows:
        try:
            df = pd.read_csv(
                file_path, comment="#", header=None,
                sep=r"\s+|,|\t", engine="python",
            )
            if df.shape[1] >= 2:
                wave = pd.to_numeric(df.iloc[:, 0], errors="coerce").to_numpy(dtype=float)
                flux = pd.to_numeric(df.iloc[:, 1], errors="coerce").to_numpy(dtype=float)
                valid = np.isfinite(wave) & np.isfinite(flux)
                if valid.sum() >= 2:
                    return wave[valid], flux[valid]
        except Exception:
            pass
        raise ValueError("Could not parse any wavelength/flux pairs from ASCII file.")

    arr = np.array(rows, dtype=float)
    return arr[:, 0], arr[:, 1]


def read_spectrum(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    p = file_path.lower()
    if any(p.endswith(ext) for ext in (".fits", ".fit", ".fits.gz", ".fit.gz", ".fz")):
        return _read_fits(file_path)
    return _read_ascii(file_path)
