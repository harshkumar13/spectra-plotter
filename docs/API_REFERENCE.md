# SpectraPlotter -- API Reference

This document describes the public Python API for SpectraPlotter modules. All functions can be imported and used independently in scripts and notebooks.

---

## `spectra_plotter.io` -- File I/O

### `read_spectrum(file_path: str) -> Tuple[np.ndarray, np.ndarray]`

Read a spectrum file and return (wavelength, flux) arrays.

**Parameters:**
- `file_path` -- Path to a FITS, ASCII, or CSV spectrum file.

**Returns:**
- `wave` -- 1D numpy array of wavelengths (Angstroms).
- `flux` -- 1D numpy array of flux values.

**Raises:**
- `ValueError` -- If the file cannot be parsed.

**Supported formats:**
- FITS binary tables (auto-detects column names)
- FITS 1D with WCS headers (CRVAL1/CDELT1/CRPIX1, including log-lambda)
- FITS 2D arrays (row 0 = wavelength, row 1 = flux)
- ASCII/text with whitespace, tab, or comma delimiters
- CSV with or without headers

**Example:**
```python
from spectra_plotter.io import read_spectrum

wave, flux = read_spectrum("SN2024abc_20240115.fits")
wave, flux = read_spectrum("spectrum.dat")
wave, flux = read_spectrum("data.csv")
```

---

### `list_spectra_files(folder: str) -> List[str]`

List all recognized spectrum files in a directory.

**Parameters:**
- `folder` -- Directory path to scan.

**Returns:**
- Sorted list of absolute file paths with recognized spectrum extensions.

**Recognized extensions:** `.fits`, `.fit`, `.fits.gz`, `.fit.gz`, `.fz`, `.txt`, `.dat`, `.ascii`, `.csv`

**Example:**
```python
from spectra_plotter.io import list_spectra_files

files = list_spectra_files("/Users/you/data/spectra/")
# ['/Users/you/data/spectra/SN2024abc.fits', '/Users/you/data/spectra/SN2024def.dat']
```

---

### `is_spectrum_file(path: str) -> bool`

Check if a file path has a recognized spectrum extension.

---

## `spectra_plotter.processing` -- Data Processing

### `clean_spectrum(wave, flux) -> Tuple[np.ndarray, np.ndarray]`

Clean a raw spectrum by removing bad data and sorting.

**Operations performed (in order):**
1. Remove entries where wavelength or flux is NaN/Inf
2. Remove entries with wavelength <= 0
3. Sort by wavelength (ascending)
4. Remove duplicate wavelength entries (keeps first)

**Parameters:**
- `wave` -- Wavelength array.
- `flux` -- Flux array (same length as wave).

**Returns:**
- Cleaned (wave, flux) tuple.

**Example:**
```python
from spectra_plotter.processing import clean_spectrum

wave_clean, flux_clean = clean_spectrum(wave, flux)
```

---

### `apply_telluric_mask(wave_obs, flux, band_keys) -> np.ndarray`

Mask flux in telluric absorption bands by setting values to NaN.

**Parameters:**
- `wave_obs` -- Observed-frame wavelength array (before redshift correction).
- `flux` -- Flux array.
- `band_keys` -- List of telluric band identifiers. Available keys:
  - `"O2_B"` -- O2 B-band (6865-6969 A)
  - `"O2_A"` -- O2 A-band (7600-7700 A)
  - `"H2O_9300"` -- H2O (9300-9600 A)
  - `"H2O_11300"` -- H2O (11300-11600 A)
  - `"H2O_13500"` -- H2O (13500-15000 A)
  - `"H2O_18140"` -- H2O (18140-19500 A)
  - `"H2O_24500"` -- H2O (24500-27000 A)

**Returns:**
- Flux array with masked regions set to NaN (copy, does not modify input).

**Example:**
```python
from spectra_plotter.processing import apply_telluric_mask

flux_masked = apply_telluric_mask(wave_obs, flux, ["O2_A", "H2O_13500"])
```

---

### `cosmic_ray_clip(flux, kernel=9, sigma=8.0) -> np.ndarray`

Conservative cosmic ray removal using median filtering and MAD clipping.

**Algorithm:**
1. Compute a median-filtered baseline (kernel size `kernel`).
2. Compute residuals (data - baseline).
3. Compute MAD (median absolute deviation) of residuals.
4. Replace points where |residual| > sigma * 1.4826 * MAD with the baseline value.

**Parameters:**
- `flux` -- Flux array.
- `kernel` -- Median filter kernel size (default: 9, must be odd).
- `sigma` -- Clipping threshold in units of scaled MAD (default: 8.0).

**Returns:**
- Clipped flux array. NaN positions are preserved.

**Example:**
```python
from spectra_plotter.processing import cosmic_ray_clip

flux_clean = cosmic_ray_clip(flux, kernel=9, sigma=8.0)
```

---

### `apply_savgol(flux, window, polyorder) -> np.ndarray`

Apply Savitzky-Golay smoothing filter.

**Parameters:**
- `flux` -- Flux array.
- `window` -- Filter window length in pixels (must be odd; will be forced odd if even).
- `polyorder` -- Polynomial order for the filter (typically 2-4).

**Returns:**
- Smoothed flux array. NaN positions are preserved.

**Notes:**
- Window is clamped to array size if too large.
- Polyorder is clamped to window-1 if too large.
- Returns original array if fewer than 5 data points or window < 3.

**Example:**
```python
from spectra_plotter.processing import apply_savgol

flux_smooth = apply_savgol(flux, window=21, polyorder=3)
```

---

### `apply_gaussian_smooth(flux, sigma) -> np.ndarray`

Apply Gaussian smoothing (convolution with a Gaussian kernel).

**Parameters:**
- `flux` -- Flux array.
- `sigma` -- Gaussian kernel standard deviation in pixels.

**Returns:**
- Smoothed flux array. NaN positions are preserved.

**Example:**
```python
from spectra_plotter.processing import apply_gaussian_smooth

flux_smooth = apply_gaussian_smooth(flux, sigma=5.0)
```

---

### `apply_binning(wave, flux, bin_size) -> Tuple[np.ndarray, np.ndarray]`

Bin a spectrum by averaging groups of adjacent pixels.

**Parameters:**
- `wave` -- Wavelength array.
- `flux` -- Flux array.
- `bin_size` -- Number of pixels per bin. Must be >= 2; if 1, returns original arrays.

**Returns:**
- (wave_binned, flux_binned) tuple with reduced length (original_length // bin_size).

**Notes:**
- Trailing pixels that don't fill a complete bin are discarded.
- Uses `np.nanmean`, so NaN values are ignored within each bin.

**Example:**
```python
from spectra_plotter.processing import apply_binning

wave_bin, flux_bin = apply_binning(wave, flux, bin_size=5)
# If wave had 1000 points, wave_bin has 200 points
```

---

### `compute_snr(flux, window=50) -> float`

Estimate signal-to-noise ratio from a central segment of the spectrum.

**Parameters:**
- `flux` -- Flux array.
- `window` -- Half-width of the central segment (default: 50, so 100 pixels total).

**Returns:**
- SNR estimate (mean / std of the central segment). Returns 0.0 if insufficient data.

**Example:**
```python
from spectra_plotter.processing import compute_snr

snr = compute_snr(flux, window=50)
print(f"SNR ~ {snr:.1f}")
```

---

### `measure_ew(wave, flux, line_center, width=20.0) -> float`

Measure the equivalent width of a spectral line.

**Parameters:**
- `wave` -- Rest-frame wavelength array.
- `flux` -- Flux array (should be in physical or normalized units).
- `line_center` -- Rest-frame wavelength of the line center (Angstroms).
- `width` -- Half-width of the integration window (Angstroms, default: 20).

**Returns:**
- Equivalent width in Angstroms. **Positive = absorption, Negative = emission.**
- Returns NaN if fewer than 3 data points in the window or if continuum is zero.

**Algorithm:**
1. Select data within [line_center - width, line_center + width].
2. Estimate continuum as the mean of the first and last 3 points.
3. Integrate (1 - flux/continuum) using the trapezoidal rule.

**Example:**
```python
from spectra_plotter.processing import measure_ew

# Measure H-alpha equivalent width
ew = measure_ew(wave, flux, line_center=6563.0, width=20.0)
print(f"EW(H-alpha) = {ew:.2f} A")
# Positive value = absorption line
# Negative value = emission line
```

---

## `spectra_plotter.line_catalog` -- Line Database

### `SPECTRAL_LINES: dict`

Dictionary of spectral line groups. Structure:
```python
{
    "H (Balmer)": {
        "H-alpha (6563 A)": 6563.0,
        "H-beta (4861 A)": 4861.0,
        ...
    },
    "He I": { ... },
    ...
}
```

### `LINE_COLORS: dict`

Mapping from element group name to hex color string for plot display.

### `TELLURIC_BANDS: dict`

Dictionary of telluric absorption bands:
```python
{
    "O2_A": ("O2 A-band (7600-7700 A)", (7600.0, 7700.0)),
    ...
}
```

### `DEFAULT_TELLURIC_MASK: list`

Default telluric band keys to mask: `["H2O_13500", "H2O_18140", "H2O_24500"]`.

---

## `spectra_plotter.app` -- Dash Application

### `create_app(upload_folder=None) -> Dash`

Create and configure the Dash application.

**Parameters:**
- `upload_folder` -- Optional path to a folder for uploaded files. If None, uses the default internal `uploads/` directory.

**Returns:**
- Configured `Dash` application instance.

**Example:**
```python
from spectra_plotter.app import create_app

app = create_app(upload_folder="/path/to/spectra")
app.run(host="127.0.0.1", port=8050, debug=True)
```

### `main()`

CLI entry point. Parses command-line arguments (`--port`, `--host`, `--folder`, `--debug`) and launches the app.
