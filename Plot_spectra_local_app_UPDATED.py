"""Plot_spectra_local_app.py

A local Dash app to load (ASCII/FITS) spectra, apply optional processing, and
visualize multiple spectra with interactive line-identification overlays.

Key improvements:
- User-controlled Savitzky–Golay smoothing or binning.
- User-selectable telluric masking (applied in observed wavelength space).
- Redshifts persist per-spectrum and never reset unless the user edits them.
- Option to scan any local folder (server-side) in addition to uploading files.
- Improved plotting UX (labels, grid, unified hover, zoom persistence).
"""

import base64
import json
import os
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from astropy.io import fits
from scipy.signal import savgol_filter, medfilt

# Dash imports with compatibility across Dash 1.x/2.x
try:
    from dash import Dash, dcc, html, callback_context  # Dash>=2
    from dash.dependencies import Input, Output, State, ALL
except Exception:  # pragma: no cover
    import dash  # type: ignore
    import dash_core_components as dcc  # type: ignore
    import dash_html_components as html  # type: ignore
    from dash.dependencies import Input, Output, State, ALL  # type: ignore
    callback_context = dash.callback_context  # type: ignore
    Dash = dash.Dash  # type: ignore

# -----------------------------
# Local storage folder (uploads)
# -----------------------------
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

SPECTRA_EXTENSIONS = (
    ".fits",
    ".fit",
    ".fits.gz",
    ".fit.gz",
    ".fz",
    ".txt",
    ".dat",
    ".ascii",
    ".csv",
)

# -----------------------------
# Spectral lines & colors
# -----------------------------
SPECTRAL_LINES = {
    "H": {
        "Hα (6563 Å)": 6563,
        "Hβ (4861 Å)": 4861,
        "Hγ (4340 Å)": 4340,
        "P9 (9229 Å)": 9229,
        "P8 (9546 Å)": 9546,
        "P7 (10052 Å)": 10052,
        "P6 (10938 Å)": 10938,
        "P5 (12818 Å)": 12818,
        "P4 (18751 Å)": 18751,
        "Brγ (21655 Å)": 21655,
    },
    "He I": {
        "He I (3888 Å)": 3888,
        "He I (5876 Å)": 5876,
        "He I (6678 Å)": 6678,
        "He I (7045 Å)": 7045,
        "He I (10830 Å)": 10830,
        "He I (20581 Å)": 20581,
    },
    "He II": {"He II (4686 Å)": 4686},
    "C I": {
        "C I (8335 Å)": 8335,
        "C I (9111 Å)": 9111.80,
        "C I (9094 Å)": 9094.80,
        "C I (9405 Å)": 9405.73,
        "C I (9658 Å)": 9658,
        "C I (10691.25 Å)": 10691,
        "C I (9405 Å)": 9405.73,
        "C I (11754 D Å)": 11754.0,
    },
    "C II": {
        "C II (5145 Å)": 5145,
        "C II (5890 Å)": 5890,
        "C II (6580 Å)": 6580,
        "C II (7234 Å)": 7234,
        "C II (17124 Å)": 17124,
    },
    "O I": {
        "O I (7773 Å)": 7773,
        "O I (8446 Å)": 8446,
        "O I (9266 Å)": 9266,
        "O I (11287 Å)": 11287,
        "O I (13164 Å)": 13164,
        "O I (13664 Å)": 13664,
    },
    "O II": {
        "O II (3738 Å)": 3738,
        "O II (3960 Å)": 3960,
        "O II (4115 Å)": 4115,
        "O II (4358 Å)": 4358,
        "O II (4651 Å)": 4651,
        "O II (7320 Å)": 7320,
    },
    "[O II]": {"[O II] (3726 Å)": 3726, "[O II] (3729 Å)": 3729},
    "[O III]": {
        "[O III] (4363 Å)": 4363,
        "[O III] (4959 Å)": 4959,
        "[O III] (5007 Å)": 5007,
        "[O III] (17525 Å)": 17525,
    },
    "N II": {"N II (6583 Å)": 6583},
    "Si I": {"Si I (10585 Å)": 10585, "Si I (12032 Å)": 12032},
    "Si II": {
        "Si II (4128 Å)": 4128,
        "Si II (4131 Å)": 4131,
        "Si II (5958 Å)": 5958,
        "Si II (5979 Å)": 5979,
        "Si II (6347 Å)": 6347,
        "Si II (6371 Å)": 6371,
        "Si II (12600 Å)": 12600,
        "Si II (15334 Å)": 15334,
        "Si II (18165 Å)": 18165,
    },
    "S I": {"S I (10457 Å)": 10457, "S I (13809 Å)": 13809, "Si I (16040 Å)": 16040},
    "S II": {"S II (6731 Å)": 6731, "S II (9069 Å)": 9069, "S II (10287 Å)": 10287},
    "Ca II": {
        "Ca II H(3934 Å)": 3934,
        "Ca II K(3964 Å)": 3964,
        "Ca II (8498 Å)": 8498,
        "Ca II (8542 Å)": 8542,
        "Ca II (8662 Å)": 8662,
    },
    "[Ca II]": {"[Ca II] (7292 Å)": 7292, "[Ca II] (7324 Å)": 7324},
    "Mg I": {"Mg I (8807 Å)": 8807, "Mg I (11828 Å)": 11828, "Mg I (15025 Å)": 15025},
    "Mg II": {
        "Mg II (2796 Å)": 2796,
        "Mg II (2798 Å)": 2798,
        "Mg II (2803 Å)": 2803,
        "Mg II (4481 Å)": 4481,
        "Mg II (7877 Å)": 7877,
        "Mg II (7896 Å)": 7896,
        "Mg II (8214 Å)": 8214,
        "Mg II (8235 Å)": 8235,
        "Mg II (9218 Å)": 9218,
        "Mg II (9244 Å)": 9244,
        "Mg II (9244 Å)": 9632,
        "Mg II (10914 Å)": 10914,
        "Mg II (12821 Å)": 12821,
        "Mg II (14431 Å)": 14431,
        "Mg II (17110 Å)": 17110,
        "Mg II (21062 Å)": 21062,
    },
    "Na I": {
        "Na I (5890 Å)": 5890,
        "Na I (5890 Å)": 5896,
        "Na I (8183 Å)": 8183,
        "Na I (8195 Å)": 8195,
        "Na I (11382 Å)": 11382,
    },
    "Fe II": {
        "Fe II (3440 Å)": 3440.6060,
        "Fe II (3581 Å)": 3581.1931,
        "Fe II (3750 Å)": 3750,
        "Fe II (4924 Å)": 4924,
        "Fe II (5018 Å)": 5018,
        "Fe II (5169 Å)": 5169,
        "Fe II (5363 Å)": 5363,
        "Fe II (9202 Å)": 9202,
        "Fe II (9998 Å)": 9998,
        "Fe II (10501 Å)": 10501,
        "Fe II (12570 Å)": 12570,
        "Fe II (16435 Å)": 16435,
    },
    "[Ar] III": {
        "[Ar III] (7751 Å)": 7751,
        "[Ar III]  (7136 Å)": 7136,
        "[Ar III]  (5192 Å)": 5192,
    },
    "Sc II": {"Sc II (6246 Å)": 6246, "Sc II (5669 Å)": 5669},
    "Ba II": {"Ba II (6142 Å)": 6142, "Ba II (6497 Å)": 6497},
    "Fe III": {"Fe III (11250 Å)": 11250, "Fe III (12297 Å)": 12297},
    "Ni II": {"Ni II (11252 Å)": 11252, "Ni II (13558 Å)": 13558},
    "Ti II": {"Ti II (12070 Å)": 12070},
    "Cr II": {"Cr II (14500 Å)": 14500},
    "Ni II": {"Ni II (11460 Å)": 11460},
    "Telluric": {
        "O₂ A-band (6865-6969 Å)": 6900,
        "O₂ A-band (7600-7700 Å)": 7650,
        "Water band (9300-9600 Å)": 9450,
        "Water band (11300-11600 Å)": 11450,
        "Water band (13500-15000 Å)": 14250,
    },
    "[Fe II]": {
        "[Fe II] (8617 Å)": 8617,
        "[Fe II] (8892 Å)": 8892,
        "[Fe II] (8926 Å)": 8926,
    },
    "[Ni II]": {
        "[Ni II] (8795 Å)": 8795,
    },
    "Host Lines": {
        "Hγ (4341 Å)": 4341,
        "Hα (6563 Å)": 6563,
        "Hβ (4861 Å)": 4861,
        "N II (6548 Å)": 6548,
        "N II (6583 Å)": 6583,
        "[O II] (3727 Å)": 3727,
        "[O III] (4959 Å)": 4959,
        "[O III] (5007 Å)": 5007,
        "Na I D (5890 Å)": 5890,
        "Na I D (5896 Å)": 5896,
        "Mg II (2798 Å)": 2798,
        "S II (6717 Å)": 6717,
        "S II (6731 Å)": 6731,
        "Ca II H (3969 Å)": 3969,
        "Ca II K (3934 Å)": 3934,
        "[S III] (9069)Å)": 9069,
    },
    "SLSN I prepeak": {
        "O II (3134.720 Å)": 3134.7205,
        "O II (3911.95 Å)": 3911.95,
        "O II (3973.2562 Å)": 3973.2562,
        "O II (4075.862 Å)": 4075.862,
        "O II (4414.905 Å)": 4414.905,
        "O II (4641.810 Å)": 4641.810,
        "O II (4649.135 Å)": 4649.135,
    },
    "SLSN I peak": {
        "O II (3135 Å)": 3134.7205,
        "O II (3912 Å)": 3911.95,
        "O II (3973 Å)": 3973.2562,
        "O II (4076 Å)": 4075.862,
        "O II (4415 Å)": 4414.905,
        "O II (4642 Å)": 4641.810,
        "O II (4649 Å)": 4649.135,
        "O II (4115 Å)": 4115,
        "O II (4358 Å)": 4358,
        "O II (3960 Å)": 3960,
        "Fe II (3440 Å)": 3440.6060,
        "Fe II (3581 Å)": 3581.1931,
        "Fe II (3750 Å)": 3750,
        "Fe II (4924 Å)": 4924,
        "Fe II (5018 Å)": 5018,
        "Fe II (5169 Å)": 5169,
        "Fe II (5363 Å)": 5363,
        "Fe III (4397 Å)": 4397,
        "Fe III (4421 Å)": 4421,
        "Fe III (4430 Å)": 4430,
        "Fe III (5129 Å)": 5129,
        "Fe III (5158 Å)": 5158,
        "Si II (6355 Å)": 6355,
        "C II (6580 Å)": 6580,
        "C II (7234 Å)": 7234,
        "O I (7234 Å)": 7774,
    },
    "SLSN I late": {
        "Mg I] (4571 Å)": 4571,
        "[Ca II] (7291 Å)": 7291,
        "[Ca II] (7324 Å)": 7324,
        "[O I] (6300 Å)": 6300,
        "[O I] (6364 Å)": 6364,
        "O I (7774 Å)": 7774,
    },
    "SLSN I nebular": {"[O II] (7319 Å)": 7319, "[O II] (7330 Å)": 7330},
}
ELEMENT_COLORS = {
    "H": "red",
    "He I": "orange",
    "He II": "yellow",
    "C I": "lightblue",  # Added color for C I
    "C II": "lightgreen",  # Added color for C II
    "[O III]": "purple",
    "O I": "green",
    "O II": "blue",
    "[O II]": "blue",
    "O III": "purple",
    "N II": "brown",
    "Si I": "teal",
    "Si II": "cyan",
    "S I": "pink",
    "S II": "magenta",
    "Ca II": "gold",
    "[Ca II]": "gold",
    "Mg I": "lime",
    "Mg II": "olive",
    "Na I": "navy",
    "Fe II": "gray",
    "Sc II": "teal",
    "Ba II": "peru",
    "[Ar] III": "lawngreen",
    "Fe III": "silver",
    "[Fe II]": "maroon",
    "Fe II": "maroon",
    "Ni II": "violet",
    "[Ni II]": "violet",
    "Ti II": "indigo",  # Added color for Ti II
    "Cr II": "darkred",  # Added color for Cr II
    "Telluric": "black",
    "Host Lines": "grey",
    "SLSN I prepeak": "violet",
    "SLSN I peak": "green",
    "SLSN I late": "red",
    "SLSN I nebular": "blue",
}

# -----------------------------
# Telluric bands (observed frame)
# -----------------------------
# These are *observed-wavelength* ranges (Å).
TELLURIC_BANDS: Dict[str, Tuple[str, Tuple[float, float]]] = {
    "O2_B": ("O₂ B-band (6865–6969 Å)", (6865.0, 6969.0)),
    "O2_A": ("O₂ A-band (7600–7700 Å)", (7600.0, 7700.0)),
    "H2O_9300": ("H₂O (9300–9600 Å)", (9300.0, 9600.0)),
    "H2O_11300": ("H₂O (11300–11600 Å)", (11300.0, 11600.0)),
    "H2O_13500": ("H₂O (13500–15000 Å)", (13500.0, 15000.0)),
    "H2O_18140": ("H₂O (18140–19500 Å)", (18140.0, 19500.0)),
    "H2O_24500": ("H₂O (24500–27000 Å)", (24500.0, 27000.0)),
}
DEFAULT_TELLURIC_MASK = ["H2O_13500", "H2O_18140", "H2O_24500"]  # preserve old behavior


# -----------------------------
# File helpers
# -----------------------------
def _is_spectrum_file(path: str) -> bool:
    base = os.path.basename(path)
    if base.startswith("."):
        return False
    p = path.lower()
    return any(p.endswith(ext) for ext in SPECTRA_EXTENSIONS)


def list_spectra_files(folder: str) -> List[str]:
    if not folder or not os.path.isdir(folder):
        return []
    paths: List[str] = []
    for name in os.listdir(folder):
        if name in {".DS_Store"}:
            continue
        full = os.path.join(folder, name)
        if os.path.isfile(full) and _is_spectrum_file(full):
            paths.append(full)
    return sorted(paths)


def read_fits_file(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    hdul = fits.open(file_path)
    if len(hdul) > 1:
        data = hdul[1].data
        if data is None or not hasattr(data, "names"):
            raise ValueError("FITS extension 1 has no table data.")

        if "OPT_WAVE" in data.names:
            wavelength = data["OPT_WAVE"]
        elif "wave" in data.names:
            wavelength = data["wave"]
        else:
            raise ValueError("Wavelength column not found in FITS file (OPT_WAVE/wave).")

        if "OPT_FLUX" in data.names:
            flux = data["OPT_FLUX"]
        elif "OPT_FLAM" in data.names:
            flux = data["OPT_FLAM"]
        elif "flux" in data.names:
            flux = data["flux"]
        else:
            raise ValueError("Flux column not found in FITS file (OPT_FLUX/OPT_FLAM/flux).")
    else:
        data = hdul[0].data
        if data is None or len(data) < 2:
            raise ValueError("FITS primary HDU has no usable [wave, flux] arrays.")
        wavelength = data[0]
        flux = data[1]

    return np.array(wavelength, dtype=float), np.array(flux, dtype=float)


def read_ascii_spectrum(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    arr = None
    try:
        arr = np.genfromtxt(file_path, comments="#", invalid_raise=False, delimiter=None)
    except Exception:
        arr = None

    if arr is None or (isinstance(arr, np.ndarray) and (arr.ndim == 1 or arr.shape[-1] < 2)):
        arr = np.genfromtxt(file_path, comments="#", invalid_raise=False, delimiter=",")

    if not isinstance(arr, np.ndarray) or arr.ndim != 2 or arr.shape[1] < 2:
        df = pd.read_csv(
            file_path,
            comment="#",
            header=None,
            sep=r"\s+|,",
            engine="python",
        )
        if df.shape[1] < 2:
            raise ValueError("ASCII spectrum must have at least two columns: wavelength flux")
        wave = df.iloc[:, 0].to_numpy(dtype=float)
        flux = df.iloc[:, 1].to_numpy(dtype=float)
        return wave, flux

    wave = arr[:, 0].astype(float)
    flux = arr[:, 1].astype(float)
    return wave, flux


def read_spectrum(file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    p = file_path.lower()
    if any(p.endswith(ext) for ext in (".fits", ".fit", ".fits.gz", ".fit.gz", ".fz")):
        return read_fits_file(file_path)
    return read_ascii_spectrum(file_path)


def clean_spectrum(wave: np.ndarray, flux: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    m = np.isfinite(wave) & np.isfinite(flux)
    wave, flux = wave[m], flux[m]
    m2 = wave > 0
    wave, flux = wave[m2], flux[m2]
    if wave.size == 0:
        return wave, flux

    order = np.argsort(wave)
    wave, flux = wave[order], flux[order]
    _, uniq_idx = np.unique(wave, return_index=True)
    uniq_idx = np.sort(uniq_idx)
    return wave[uniq_idx], flux[uniq_idx]


# -----------------------------
# Processing helpers
# -----------------------------
def apply_telluric_mask_observed(
    wave_obs: np.ndarray,
    flux: np.ndarray,
    selected_band_keys: Sequence[str],
) -> np.ndarray:
    if not selected_band_keys:
        return flux
    masked_flux = flux.astype(float).copy()
    for key in selected_band_keys:
        if key not in TELLURIC_BANDS:
            continue
        _, (start, end) = TELLURIC_BANDS[key]
        band_mask = (wave_obs >= start) & (wave_obs <= end)
        masked_flux[band_mask] = np.nan
    return masked_flux


def fill_nans_1d(y: np.ndarray) -> np.ndarray:
    y = y.astype(float).copy()
    n = y.size
    if n == 0:
        return y
    x = np.arange(n)
    good = np.isfinite(y)
    if good.all():
        return y
    if not good.any():
        return np.zeros_like(y)
    y[~good] = np.interp(x[~good], x[good], y[good])
    return y


def robust_cosmic_ray_clip(flux: np.ndarray, kernel: int = 9, sigma: float = 8.0) -> np.ndarray:
    y = flux.astype(float).copy()
    nan_mask = ~np.isfinite(y)
    y_filled = fill_nans_1d(y)

    if kernel < 3:
        kernel = 3
    if kernel % 2 == 0:
        kernel += 1

    baseline = medfilt(y_filled, kernel_size=kernel)
    resid = y_filled - baseline
    mad = np.median(np.abs(resid - np.median(resid)))
    if mad == 0 or not np.isfinite(mad):
        return flux
    threshold = sigma * 1.4826 * mad
    spike = np.abs(resid) > threshold

    y_out = y_filled.copy()
    y_out[spike] = baseline[spike]
    y_out[nan_mask] = np.nan
    return y_out


def apply_savgol(flux: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    y = flux.astype(float).copy()
    nan_mask = ~np.isfinite(y)
    y_filled = fill_nans_1d(y)

    n = y_filled.size
    if n < 5:
        return y

    window = int(window) if window is not None else 0
    polyorder = int(polyorder) if polyorder is not None else 2

    if window < 3:
        return y
    if window % 2 == 0:
        window += 1
    window = min(window, n if n % 2 == 1 else n - 1)
    window = max(window, 3)

    polyorder = max(1, polyorder)
    polyorder = min(polyorder, window - 1)

    y_smooth = savgol_filter(y_filled, window_length=window, polyorder=polyorder, mode="interp")
    y_smooth[nan_mask] = np.nan
    return y_smooth


def apply_binning(wave: np.ndarray, flux: np.ndarray, bin_size: int) -> Tuple[np.ndarray, np.ndarray]:
    if bin_size is None or bin_size <= 1:
        return wave, flux
    bin_size = int(bin_size)
    n = wave.size
    m = (n // bin_size) * bin_size
    if m < bin_size:
        return wave, flux

    w = wave[:m].reshape(-1, bin_size)
    f = flux[:m].reshape(-1, bin_size)
    w_b = np.nanmean(w, axis=1)
    f_b = np.nanmean(f, axis=1)
    return w_b, f_b


def interp_ignore_nan(x: np.ndarray, y: np.ndarray, x_new: np.ndarray) -> np.ndarray:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    return np.interp(x_new, x[m], y[m])


# -----------------------------
# Dash app
# -----------------------------
app = Dash(__name__)
server = app.server

initial_options = [{"label": os.path.basename(p), "value": p} for p in list_spectra_files(UPLOAD_FOLDER)]

app.layout = html.Div(
    [
        dcc.Store(id="redshift-store", data={}),
        html.H1("Spectra Plotting Tool", style={"textAlign": "center", "marginBottom": "0.2rem"}),
        html.Div(
            [
                html.Div(
                    [
                        html.H3("Data source", style={"marginBottom": "0.25rem"}),
                        dcc.RadioItems(
                            id="file-source",
                            options=[
                                {"label": "Uploads folder", "value": "uploads"},
                                {"label": "Scan a local folder", "value": "folder"},
                            ],
                            value="uploads",
                            inline=True,
                        ),
                        html.Div(
                            [
                                html.Label("Folder path (on this machine):"),
                                dcc.Input(
                                    id="folder-path",
                                    type="text",
                                    value=os.path.abspath(UPLOAD_FOLDER),
                                    style={"width": "70%"},
                                ),
                                html.Button("Scan", id="scan-folder", n_clicks=0, style={"marginLeft": "0.5rem"}),
                                html.Div(
                                    "Tip: Since this app runs locally, scanning lets you plot files without uploading them.",
                                    style={"fontSize": "0.85rem", "opacity": 0.8, "marginTop": "0.25rem"},
                                ),
                            ],
                            id="folder-controls",
                            style={"display": "none", "marginTop": "0.5rem"},
                        ),
                        html.Div(
                            [
                                html.Label("Upload spectra files (ASCII/FITS):"),
                                dcc.Upload(
                                    id="upload-data",
                                    children=html.Div(["Drag & drop or ", html.Button("Browse…")]),
                                    multiple=True,
                                    style={
                                        "width": "100%",
                                        "height": "60px",
                                        "lineHeight": "60px",
                                        "borderWidth": "1px",
                                        "borderStyle": "dashed",
                                        "borderRadius": "8px",
                                        "textAlign": "center",
                                        "marginTop": "0.5rem",
                                    },
                                ),
                                html.Div(
                                    f"Uploaded files are saved to: {os.path.abspath(UPLOAD_FOLDER)}",
                                    style={"fontSize": "0.85rem", "opacity": 0.8, "marginTop": "0.25rem"},
                                ),
                            ],
                            id="upload-controls",
                            style={"marginTop": "0.5rem"},
                        ),
                    ],
                    style={"flex": "1 1 420px", "padding": "0.75rem", "border": "1px solid #ddd", "borderRadius": "12px"},
                ),
                html.Div(
                    [
                        html.H3("Display & processing", style={"marginBottom": "0.25rem"}),
                        html.Label("Select spectra to display:"),
                        dcc.Dropdown(
                            id="file-selector",
                            options=initial_options,
                            value=[],
                            multi=True,
                            placeholder="Select one or more spectra…",
                        ),
                        html.Div(id="redshift-inputs", style={"marginTop": "0.75rem"}),

                        html.Div(
                            [
                                html.Label("Velocity (km/s):"),
                                dcc.Input(id="velocity-input", type="number", value=0, step=1),
                            ],
                            style={"marginTop": "0.75rem"},
                        ),

                        html.Hr(style={"margin": "0.75rem 0"}),

                        html.Label("Telluric masking (observed wavelengths):"),
                        dcc.Checklist(
                            id="telluric-mask",
                            options=[{"label": v[0], "value": k} for k, v in TELLURIC_BANDS.items()],
                            value=DEFAULT_TELLURIC_MASK,
                        ),

                        html.Hr(style={"margin": "0.75rem 0"}),

                        html.Label("Smoothing / binning:"),
                        dcc.RadioItems(
                            id="proc-mode",
                            options=[
                                {"label": "None", "value": "none"},
                                {"label": "Savitzky–Golay smoothing", "value": "savgol"},
                                {"label": "Binning (mean)", "value": "bin"},
                            ],
                            value="none",
                        ),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.Label("Savitzky–Golay window (points, odd):"),
                                        dcc.Slider(
                                            id="savgol-window",
                                            min=3,
                                            max=301,
                                            step=2,
                                            value=11,
                                            tooltip={"placement": "bottom", "always_visible": False},
                                        ),
                                    ],
                                    style={"marginTop": "0.5rem"},
                                ),
                                html.Div(
                                    [
                                        html.Label("Savitzky–Golay polyorder:"),
                                        dcc.Slider(
                                            id="savgol-polyorder",
                                            min=1,
                                            max=6,
                                            step=1,
                                            value=3,
                                            tooltip={"placement": "bottom", "always_visible": False},
                                        ),
                                    ],
                                    style={"marginTop": "0.75rem"},
                                ),
                            ],
                            id="savgol-controls",
                            style={"display": "none"},
                        ),
                        html.Div(
                            [
                                html.Label("Bin size (points):"),
                                dcc.Input(id="bin-size", type="number", value=5, min=2, step=1),
                            ],
                            id="bin-controls",
                            style={"display": "none", "marginTop": "0.5rem"},
                        ),

                        html.Hr(style={"margin": "0.75rem 0"}),

                        html.Label("Normalization:"),
                        dcc.RadioItems(
                            id="norm-mode",
                            options=[
                                {"label": "Off", "value": "off"},
                                {"label": "Each spectrum (max)", "value": "each"},
                                {"label": "Common overlap (max)", "value": "common"},
                            ],
                            value="common",
                        ),
                        html.Div(
                            [
                                html.Label("Y-axis:"),
                                dcc.RadioItems(
                                    id="y-scale",
                                    options=[{"label": "Linear", "value": "linear"}, {"label": "Log", "value": "log"}],
                                    value="linear",
                                    inline=True,
                                ),
                            ],
                            style={"marginTop": "0.5rem"},
                        ),
                        dcc.Checklist(
                            id="cr-clip",
                            options=[{"label": "Conservative spike removal (cosmic rays)", "value": "on"}],
                            value=[],
                            style={"marginTop": "0.5rem"},
                        ),
                    ],
                    style={"flex": "1 1 520px", "padding": "0.75rem", "border": "1px solid #ddd", "borderRadius": "12px"},
                ),
            ],
            style={"display": "flex", "flexWrap": "wrap", "gap": "0.75rem", "margin": "0.75rem 0"},
        ),

        html.Div(
            [
                html.H3("Line IDs"),
                html.Div(
                    [
                        html.Div(
                            [
                                html.Label(element, style={"color": ELEMENT_COLORS.get(element, "black")}),
                                dcc.Checklist(
                                    id={"type": "element-line-selector", "index": element},
                                    options=[{"label": line, "value": line} for line in lines],
                                    value=[],
                                    inline=True,
                                ),
                            ],
                            style={"marginBottom": "0.25rem"},
                        )
                        for element, lines in SPECTRAL_LINES.items()
                    ],
                    style={"columnCount": 2, "columnGap": "1.5rem"},
                ),
            ],
            style={"padding": "0.75rem", "border": "1px solid #ddd", "borderRadius": "12px"},
        ),

        html.Div([dcc.Graph(id="spectra-plot", style={"height": "800px"})], style={"marginTop": "0.75rem"}),
    ],
    style={"maxWidth": "1600px", "margin": "0 auto", "padding": "0.75rem"},
)

@app.callback(
    Output("folder-controls", "style"),
    Output("upload-controls", "style"),
    Input("file-source", "value"),
)
def toggle_source_controls(source):
    if source == "folder":
        return {"display": "block", "marginTop": "0.5rem"}, {"display": "none"}
    return {"display": "none"}, {"display": "block", "marginTop": "0.5rem"}


@app.callback(
    Output("savgol-controls", "style"),
    Output("bin-controls", "style"),
    Input("proc-mode", "value"),
)
def toggle_processing_controls(mode):
    if mode == "savgol":
        return {"display": "block"}, {"display": "none"}
    if mode == "bin":
        return {"display": "none"}, {"display": "block", "marginTop": "0.5rem"}
    return {"display": "none"}, {"display": "none"}


@app.callback(
    Output("file-selector", "options"),
    Output("file-selector", "value"),
    Input("upload-data", "contents"),
    State("upload-data", "filename"),
    Input("scan-folder", "n_clicks"),
    State("folder-path", "value"),
    Input("file-source", "value"),
    State("file-selector", "value"),
)
def update_file_options(upload_contents, upload_filenames, scan_clicks, folder_path, source, current_selection):
    if upload_contents and upload_filenames:
        for content, filename in zip(upload_contents, upload_filenames):
            try:
                data = content.split(",", 1)[1]
                with open(os.path.join(UPLOAD_FOLDER, os.path.basename(filename)), "wb") as f:
                    f.write(base64.b64decode(data))
            except Exception:
                pass

    scan_root = UPLOAD_FOLDER if source == "uploads" else (folder_path or UPLOAD_FOLDER)
    files = list_spectra_files(scan_root)

    options = [{"label": os.path.basename(p), "value": p} for p in files]

    current_selection = current_selection or []
    selection = [p for p in current_selection if p in files]

    trig = callback_context.triggered[0]["prop_id"].split(".", 1)[0] if callback_context.triggered else ""
    if trig == "upload-data" and source == "uploads" and upload_filenames:
        new_paths = [os.path.join(UPLOAD_FOLDER, os.path.basename(fn)) for fn in upload_filenames]
        for p in new_paths:
            if p in files and p not in selection:
                selection.append(p)

    return options, selection


@app.callback(
    Output("redshift-inputs", "children"),
    Input("file-selector", "value"),
    State("redshift-store", "data"),
)
def generate_redshift_inputs(selected_files, redshift_store):
    if not selected_files:
        return []
    redshift_store = redshift_store or {}

    children = []
    for filepath in selected_files:
        z = redshift_store.get(filepath, 0.0)
        children.append(
            html.Div(
                [
                    html.Label(f"Redshift for {os.path.basename(filepath)}:"),
                    dcc.Input(
                        id={"type": "redshift-input", "index": filepath},
                        type="number",
                        value=float(z),
                        step=0.00001,
                        debounce=False,
                        persistence=True,
                        persistence_type="session",
                    ),
                ],
                style={"marginTop": "0.5rem"},
            )
        )
    return children


@app.callback(
    Output("redshift-store", "data"),
    Input({"type": "redshift-input", "index": ALL}, "value"),
    State({"type": "redshift-input", "index": ALL}, "id"),
    State("redshift-store", "data"),
)
def update_redshift_store(values, ids, store):
    store = store or {}
    if not ids:
        return store
    if not callback_context.triggered:
        return store

    trig = callback_context.triggered[0]["prop_id"].split(".", 1)[0]
    if not trig or not trig.startswith("{"):
        return store
    try:
        trig_id = json.loads(trig)
    except Exception:
        return store

    for v, _id in zip(values, ids):
        if _id == trig_id:
            if v is None:
                return store
            store[_id["index"]] = float(v)
            return store
    return store


@app.callback(
    Output("spectra-plot", "figure"),
    Input("file-selector", "value"),
    Input("redshift-store", "data"),
    Input({"type": "element-line-selector", "index": ALL}, "value"),
    Input("velocity-input", "value"),
    Input("telluric-mask", "value"),
    Input("proc-mode", "value"),
    Input("savgol-window", "value"),
    Input("savgol-polyorder", "value"),
    Input("bin-size", "value"),
    Input("norm-mode", "value"),
    Input("y-scale", "value"),
    Input("cr-clip", "value"),
    State("spectra-plot", "relayoutData"),
)
def update_plot(
    selected_files,
    redshift_store,
    selected_lines_nested,
    velocity,
    telluric_mask_keys,
    proc_mode,
    sg_window,
    sg_poly,
    bin_size,
    norm_mode,
    y_scale,
    cr_clip_value,
    relayoutData,
):
    fig = go.Figure()
    fig.update_layout(
        width=1600,
        height=800,
        hovermode="x unified",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        margin=dict(l=60, r=20, t=30, b=55),
        uirevision="spectra-uirevision",
    )
    fig.update_xaxes(title_text="Wavelength (Å)", showgrid=True, zeroline=False)
    fig.update_yaxes(title_text="Flux (normalized)" if norm_mode != "off" else "Flux", showgrid=True, zeroline=False)

    if not selected_files:
        fig.add_annotation(
            text="Select or upload spectra to begin.",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=16),
        )
        return fig

    redshift_store = redshift_store or {}
    telluric_mask_keys = telluric_mask_keys or []
    selected_lines_nested = selected_lines_nested or []
    selected_lines = [line for sublist in selected_lines_nested for line in (sublist or [])]

    velocity = float(velocity) if velocity is not None else 0.0
    do_cr_clip = (cr_clip_value or []) == ["on"] or ("on" in (cr_clip_value or []))

    spectra_data: List[Tuple[str, np.ndarray, np.ndarray]] = []

    for filepath in selected_files:
        try:
            wave_obs, flux = read_spectrum(filepath)
            wave_obs, flux = clean_spectrum(wave_obs, flux)
            if wave_obs.size < 2:
                continue

            flux = apply_telluric_mask_observed(wave_obs, flux, telluric_mask_keys)

            if do_cr_clip:
                flux = robust_cosmic_ray_clip(flux, kernel=9, sigma=8.0)

            z = float(redshift_store.get(filepath, 0.0) or 0.0)
            wave = wave_obs / (1.0 + z)

            if proc_mode == "savgol":
                flux = apply_savgol(flux, window=int(sg_window or 0), polyorder=int(sg_poly or 2))
            elif proc_mode == "bin":
                wave, flux = apply_binning(wave, flux, int(bin_size or 1))

            spectra_data.append((os.path.basename(filepath), wave, flux))
        except Exception as e:
            print(f"Error reading {filepath}: {e}")

    if not spectra_data:
        fig.add_annotation(
            text="No readable spectra found in selection.",
            x=0.5,
            y=0.5,
            xref="paper",
            yref="paper",
            showarrow=False,
            font=dict(size=16),
        )
        return fig

    # Normalization
    if norm_mode == "common" and len(spectra_data) > 1:
        mins = [np.nanmin(w) for _, w, _ in spectra_data if np.isfinite(w).any()]
        maxs = [np.nanmax(w) for _, w, _ in spectra_data if np.isfinite(w).any()]
        if mins and maxs:
            lo = float(max(mins))
            hi = float(min(maxs))
        else:
            lo, hi = np.nan, np.nan

        if np.isfinite(lo) and np.isfinite(hi) and hi > lo:
            common_wavelength = np.linspace(lo, hi, 1500)
            norms: List[float] = []
            for _, w, f in spectra_data:
                interp_f = interp_ignore_nan(w, f, common_wavelength)
                nrm = np.nanmax(interp_f)
                norms.append(float(nrm) if np.isfinite(nrm) and nrm != 0 else 1.0)

            spectra_data = [(name, w, f / nrm) for (name, w, f), nrm in zip(spectra_data, norms)]
        else:
            norm_mode = "each"

    if norm_mode == "each":
        spectra_data = [
            (name, w, f / (float(np.nanmax(f)) if np.isfinite(np.nanmax(f)) and np.nanmax(f) != 0 else 1.0))
            for name, w, f in spectra_data
        ]

    # Traces
    # If using log scale, mask non-positive values (Plotly log axes require >0)
    if y_scale == "log":
        spectra_data = [(name, w, np.where(f > 0, f, np.nan)) for name, w, f in spectra_data]

    for name, w, f in spectra_data:
        fig.add_trace(
            go.Scatter(
                x=w,
                y=f,
                mode="lines",
                name=name,
                line=dict(width=1.6),
                hovertemplate="λ=%{x:.1f} Å<br>flux=%{y:.4g}<extra></extra>",
                connectgaps=False,
            )
        )

    if y_scale == "log":
        fig.update_yaxes(type="log")

    if relayoutData:
        x0 = relayoutData.get("xaxis.range[0]")
        x1 = relayoutData.get("xaxis.range[1]")
        y0 = relayoutData.get("yaxis.range[0]")
        y1 = relayoutData.get("yaxis.range[1]")
        if x0 is not None and x1 is not None:
            fig.update_xaxes(range=[x0, x1])
        if y0 is not None and y1 is not None:
            fig.update_yaxes(range=[y0, y1])

    # Vertical lines (rest frame)
    y_max = 1.0
    try:
        y_max = max(np.nanmax(trace.y) for trace in fig.data if getattr(trace, "mode", "") == "lines")
        if not np.isfinite(y_max) or y_max <= 0:
            y_max = 1.0
    except Exception:
        y_max = 1.0

    for element, lines in SPECTRAL_LINES.items():
        for label, w0 in lines.items():
            if label not in selected_lines:
                continue
            shifted = float(w0) * (1.0 + velocity / 3e5)
            fig.add_trace(
                go.Scatter(
                    x=[shifted, shifted],
                    y=[0, y_max],
                    mode="lines",
                    line=dict(color=ELEMENT_COLORS.get(element, "black"), dash="dash"),
                    name=label,
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    return fig


if __name__ == "__main__":
    app.run_server(debug=True)
