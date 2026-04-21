# SpectraPlotter

**Interactive astronomical spectra visualization and line identification tool.**

A browser-based application for loading, processing, and visualizing astronomical spectra with real-time line identification overlays. Designed for transient astronomers working with optical and NIR spectroscopy, inspired by the plotting tools on the [Transient Name Server (TNS)](https://www.wis-tns.org/) and [WISeREP](https://www.wiserep.org/).

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Step-by-Step User Guide](#step-by-step-user-guide)
- [Supported File Formats](#supported-file-formats)
- [Spectral Line Catalog](#spectral-line-catalog)
- [Command-Line Options](#command-line-options)
- [API / Programmatic Usage](#api--programmatic-usage)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [Contributing](#contributing)
- [Contributors](#contributors)
- [License](#license)
- [Acknowledgments](#acknowledgments)

---

## Features

- **Multi-format input** -- FITS (binary tables, 1D WCS, 2D arrays, multi-extension, log-lambda), ASCII (`.txt`, `.dat`, `.ascii`), and CSV files with automatic delimiter and header detection
- **Interactive plots** -- Plotly-powered with zoom, pan, scroll-zoom, hover readouts, measurement tools, and high-res PNG export
- **200+ spectral lines** -- Organized across 35+ element/ion groups (H, He, C, N, O, Fe, Si, Ca, Mg, Na, etc.)
- **Supernova line sets** -- Pre-configured groups for SLSN-I at different phases (pre-peak, peak, late, nebular) and host galaxy lines
- **Per-spectrum redshift** -- Individual z value per file, session-persistent
- **Velocity shift** -- Apply radial velocity offsets (km/s) to all line markers
- **Telluric masking** -- Selectable atmospheric absorption bands (O2, H2O) masked in the observed frame
- **Signal processing** -- Savitzky-Golay smoothing, Gaussian smoothing, wavelength binning, cosmic ray removal
- **Normalization** -- Individual-max, common-overlap, or raw flux
- **Dark theme** -- Eye-friendly interface for observatory / nighttime use
- **CSV export** -- Download processed rest-frame spectra
- **Local-first** -- Runs on your machine; scan any folder directly

---

## Installation

### Prerequisites

- Python 3.8 or newer
- pip (Python package manager)

### Option A: Install from source (recommended for development)

```bash
git clone https://github.com/harsh-kumar-agrawal/spectra-plotter.git
cd spectra-plotter
pip install -e .
```

### Option B: Install dependencies only

```bash
git clone https://github.com/harsh-kumar-agrawal/spectra-plotter.git
cd spectra-plotter
pip install -r requirements.txt
```

### Verify installation

```bash
python -c "from spectra_plotter.app import create_app; print('OK')"
```

---

## Quick Start

```bash
# Launch the app
python -m spectra_plotter

# Or, if installed as a package:
spectra-plotter

# Point to a folder of spectra directly:
spectra-plotter --folder /path/to/your/spectra
```

Open **http://127.0.0.1:8050** in your browser. That's it.

---

## Step-by-Step User Guide

The interface has three main sections arranged from top to bottom:

```
+---------------------------------------------------------------+
|  Title bar                                                    |
+---------------------------------------------------------------+
|  [Upload/Folder] | [File selector dropdown] | [Buttons]      |  <-- Top bar
+-------------------------------+-------------------------------+
|                               |                               |
|   Interactive Plot            |   Line Identification         |
|                               |   (4-column checkbox grid)    |
|                               |                               |
+-------------------------------+-------------------------------+
|  Toolbar: Zoom | Binning | Smooth | Norm | CR clip | Log Y   |
+---------------------------------------------------------------+
|  z= inputs for each selected spectrum                         |
+---------------------------------------------------------------+
|  Status bar (monospace)                                       |
+---------------------------------------------------------------+
```

### Step 1: Load your spectra

You have two ways to get spectra into the tool:

**Option A -- Upload files**
1. The top-left shows a radio toggle: **Upload** / **Local folder**. Keep it on **Upload**.
2. Drag-and-drop your spectrum files onto the upload area (or click **browse**).
3. Files are saved to the internal `uploads/` folder.

**Option B -- Scan a local folder**
1. Switch the radio to **Local folder**.
2. A text input + **Scan** button appears. Enter the full path to a folder containing your spectra, e.g., `/Users/you/data/SN2024abc/spectra/`.
3. Click **Scan**. All recognized spectrum files in that folder appear in the dropdown.
4. This is the faster option since files are read directly without copying.

### Step 2: Select spectra to plot

1. Click the **file selector dropdown** (center of the top bar). It lists all available files.
2. Click individual files to select them, or use the **Select all** button.
3. To deselect everything, click **Clear**.
4. To re-scan the folder after adding new files, click **Reload**.

Each selected spectrum appears as a colored trace on the plot.

### Step 3: Set redshifts

Once you select spectra, a row of **z= inputs** appears below the toolbar:

```
z= [0.072] SN2024abc_20240115.fits   z= [0.072] SN2024abc_20240124.fits
```

1. Type the redshift value for each spectrum.
2. The plot immediately updates: wavelengths are corrected to the rest frame via `lambda_rest = lambda_obs / (1 + z)`.
3. Redshifts are remembered for the duration of your browser session.

### Step 4: Identify spectral lines

The **right panel** shows a 4-column grid of element groups. Each group has:
- A **checkbox next to the element name** -- toggle this to select/deselect ALL lines in that group at once
- Individual line checkboxes below it

**Typical workflow for SN classification:**
1. Check the **H (Balmer)** group checkbox to show H-alpha, H-beta, H-gamma, H-delta
2. Check **He I** to test for helium
3. For SLSN-I, check **SLSN-I Peak** to overlay the characteristic O II / Fe II / Fe III lines at once
4. Click **Clear all** (top of the panel) to remove all line markers

**Velocity shift:**
- At the top of the line ID panel, enter a **v-shift** in km/s.
- This shifts ALL line markers by `lambda_shifted = lambda_rest * (1 + v/c)`.
- Useful for measuring photospheric velocities.

**Telluric lines:**
- Below the velocity input, checkboxes for telluric bands (O2, H2O) let you show/hide these markers.
- By default, the strongest NIR water bands are masked (data set to NaN in those regions).

### Step 5: Process the data

The **toolbar below the plot** gives you quick access to processing:

| Control | What it does |
|---------|-------------|
| **Zoom Full** | Reset plot to show the full wavelength range |
| **Auto Zoom** | Auto-fit to the data extent |
| **Binning** | Enter a number > 1 to bin the spectrum (mean of N adjacent pixels). Set to 1 for no binning. |
| **Smooth** | Dropdown: `None`, `SavGol` (Savitzky-Golay), or `Gauss` (Gaussian). The number field sets the window size (SavGol) or is used to derive sigma (Gaussian). |
| **Norm** | `Off` = raw flux, `Each` = normalize each spectrum to its own max, `Common` = normalize all to their common overlap region |
| **CR clip** | Check to enable conservative cosmic ray removal (median-filter + MAD clipping) |
| **Log Y** | Check to switch y-axis to logarithmic scale |

**Example: smooth a noisy spectrum**
1. Set **Smooth** dropdown to `SavGol`
2. Set the window to `21` (must be odd; larger = smoother)
3. The plot updates in real time

**Example: bin a low-S/N NIR spectrum**
1. Set **Binning** to `5`
2. Each 5 adjacent pixels are averaged, reducing noise by ~sqrt(5)

### Step 6: Interact with the plot

The plot supports several mouse interactions:

| Action | How |
|--------|-----|
| **Zoom in** | Click-drag a rectangle on the plot |
| **Scroll zoom** | Mouse scroll wheel zooms in/out |
| **Pan** | Hold shift + drag, or use the pan tool in the toolbar |
| **Hover readout** | Hover over any point to see `wl=XXXX.X A  flux=Y.YYY` |
| **Reset zoom** | Click **Zoom Full** button, or double-click the plot |
| **Save image** | Click the camera icon in the Plotly toolbar (top-right of plot) to download a high-res PNG |
| **Draw lines** | Use the line-draw tool (in Plotly toolbar) to annotate features |

### Step 7: Export data

Click **Export CSV** in the top bar. A CSV file downloads containing:

```csv
filename,wavelength_rest_A,flux
SN2024abc_20240115.fits,3500.1234,1.23e-16
SN2024abc_20240115.fits,3501.5678,1.25e-16
...
```

The exported wavelengths are redshift-corrected (rest frame). One row per data point, multiple spectra concatenated.

### Step 8: Status bar

The bottom bar shows a monospace summary:

```
2 spectra | SN2024abc_20240115.fits z=0.0720 SNR~45.3 (3400-9200A, 4521pts) | ...
```

This gives you a quick check: number of loaded spectra, redshift, estimated signal-to-noise ratio, wavelength coverage, and number of data points.

---

## Supported File Formats

SpectraPlotter automatically detects the format based on file extension and content:

| Format | Extensions | How it's read |
|--------|-----------|---------------|
| **FITS binary table** | `.fits`, `.fit`, `.fits.gz`, `.fit.gz`, `.fz` | Searches all extensions for columns matching wavelength (`OPT_WAVE`, `wave`, `WAVE`, `wavelength`, `lam`, `lambda`, etc.) and flux (`OPT_FLUX`, `OPT_FLAM`, `flux`, `FLUX`, `flam`, `counts`, etc.). Partial name matching is supported. |
| **FITS 1D + WCS** | `.fits` | Single-array HDU with `CRVAL1`, `CDELT1`/`CD1_1`, `CRPIX1` header keywords. Supports `DC-FLAG=1` (log-lambda) spectra. |
| **FITS 2D array** | `.fits` | Primary HDU with shape `[2, N]` or `[>=2, N]` -- row 0 = wavelength, row 1 = flux. |
| **Whitespace-delimited** | `.txt`, `.dat`, `.ascii` | Two or more columns separated by spaces/tabs. Lines starting with `#` or `%` are skipped. Header rows are auto-detected and skipped. |
| **Tab-delimited** | `.txt`, `.dat` | Auto-detected when tabs are present. |
| **CSV** | `.csv` | Comma-separated, with or without a header row. |

**Tips for tricky files:**
- The parser tries multiple strategies in order: numpy `genfromtxt`, manual line-by-line parsing, then pandas as a fallback.
- Column headers (if present) are skipped automatically.
- Files with mixed comment styles (`#`, `%`) are handled.
- Wavelength must be in the first column, flux in the second.

---

## Spectral Line Catalog

The built-in catalog includes rest-frame wavelengths organized by element/ion:

| Group | Lines | Key features |
|-------|-------|-------------|
| H (Balmer) | 4 | H-alpha (6563), H-beta (4861), H-gamma (4340), H-delta (4102) |
| H (Paschen) | 6 | Pa4 through Pa9 (9229--18751 A) |
| H (Brackett) | 1 | Br-gamma (21655 A) |
| He I | 6 | 3888, 5876, 6678, 7045, 10830, 20581 A |
| He II | 1 | 4686 A |
| C I | 7 | 8335--11754 A |
| C II | 5 | 5145--17124 A |
| O I | 6 | 7773--13664 A |
| O II | 6 | 3738--7320 A |
| [O I] | 2 | 6300, 6364 A |
| [O II] | 2 | 3726, 3729 A |
| [O III] | 4 | 4363, 4959, 5007, 17525 A |
| Si II | 9 | 4128--18165 A |
| Ca II | 5 | H&K (3934, 3969) + NIR triplet (8498, 8542, 8662) |
| Fe II | 12 | 3441--16435 A |
| Fe III | 7 | 4397--12297 A |
| Na I | 5 | D doublet (5890, 5896) + NIR |
| Mg II | 15 | 2796--21062 A |
| Host Galaxy | 16 | Common nebular emission lines for redshift measurement |
| SLSN-I Pre-peak | 7 | O II lines characteristic of SLSN-I before peak |
| SLSN-I Peak | 26 | O II + Fe II + Fe III + Si II + C II lines at peak |
| SLSN-I Late | 6 | Mg I], [Ca II], [O I], O I lines in late phase |
| SLSN-I Nebular | 2 | [O II] 7319, 7330 A in nebular phase |
| Telluric | 5 | Atmospheric O2 and H2O absorption markers |

---

## Command-Line Options

```
usage: spectra-plotter [-h] [--port PORT] [--host HOST] [--folder FOLDER] [--debug]

SpectraPlotter - Astronomical spectra viewer

options:
  --port PORT      Port number (default: 8050)
  --host HOST      Host address (default: 127.0.0.1; use 0.0.0.0 for network access)
  --folder FOLDER  Default spectra folder to scan on startup
  --debug          Enable Dash debug mode (auto-reload on code changes)
```

**Examples:**

```bash
# Default
spectra-plotter

# Point to a specific data folder
spectra-plotter --folder ~/data/SN2024abc/spectra/

# Make accessible on local network
spectra-plotter --host 0.0.0.0 --port 8080

# Development mode
spectra-plotter --debug
```

---

## API / Programmatic Usage

You can also embed SpectraPlotter in your own scripts or Jupyter workflows:

```python
from spectra_plotter.app import create_app

# Create app pointing to a specific folder
app = create_app(upload_folder="/path/to/spectra")

# Run it
app.run(host="127.0.0.1", port=8050, debug=True)
```

The processing functions are available independently:

```python
from spectra_plotter.io import read_spectrum
from spectra_plotter.processing import clean_spectrum, apply_savgol, compute_snr

# Read any supported format
wave, flux = read_spectrum("spectrum.fits")
wave, flux = clean_spectrum(wave, flux)

# Smooth
flux_smooth = apply_savgol(flux, window=21, polyorder=3)

# Measure SNR
snr = compute_snr(flux)
print(f"SNR ~ {snr:.1f}")
```

---

## Project Structure

```
spectra-plotter/
├── spectra_plotter/
│   ├── __init__.py          # Package version and metadata
│   ├── __main__.py          # `python -m spectra_plotter` entry point
│   ├── app.py               # Dash application, layout, and callbacks
│   ├── io.py                # File I/O (FITS, ASCII, CSV readers)
│   ├── processing.py        # Spectral processing (clean, smooth, bin, CR clip)
│   ├── line_catalog.py      # Spectral line database, colors, telluric bands
│   └── assets/
│       └── style.css        # Dark theme CSS
├── tests/
│   ├── __init__.py
│   └── test_processing.py   # Unit tests for processing functions
├── pyproject.toml            # Package configuration (pip install -e .)
├── requirements.txt          # Dependency list
├── LICENSE                   # MIT License
├── .gitignore
└── README.md                 # This file
```

### Module overview

| Module | Purpose |
|--------|---------|
| `app.py` | Dash web application: layout (WISeREP-style), all callbacks for file selection, redshift, plotting, line overlays, export |
| `io.py` | Reads spectra from FITS (binary table, WCS, 2D), ASCII, CSV. Handles comment lines, headers, delimiters automatically. |
| `processing.py` | `clean_spectrum` (NaN/sort/dedup), `apply_telluric_mask`, `cosmic_ray_clip`, `apply_savgol`, `apply_gaussian_smooth`, `apply_binning`, `compute_snr`, `measure_ew` |
| `line_catalog.py` | `SPECTRAL_LINES` dict (200+ lines in 35+ groups), `LINE_COLORS`, `TELLURIC_BANDS` |

---

## Dependencies

| Package | Minimum version | Purpose |
|---------|----------------|---------|
| [Dash](https://dash.plotly.com/) | 2.0 | Web application framework |
| [Plotly](https://plotly.com/python/) | 5.0 | Interactive plotting |
| [NumPy](https://numpy.org/) | 1.20 | Array operations |
| [Pandas](https://pandas.pydata.org/) | 1.3 | CSV/table parsing fallback |
| [SciPy](https://scipy.org/) | 1.7 | Savitzky-Golay, Gaussian, median filters |
| [Astropy](https://www.astropy.org/) | 5.0 | FITS file I/O |

---

## Running Tests

```bash
# From the project root
python -m pytest tests/ -v

# With coverage
python -m pytest tests/ --cov=spectra_plotter --cov-report=term-missing
```

---

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Make your changes
4. Run the test suite (`python -m pytest tests/ -v`)
5. Commit (`git commit -am 'Add new feature'`)
6. Push (`git push origin feature/new-feature`)
7. Open a Pull Request

### Ideas for contributions

- Spectral template matching / cross-correlation
- Continuum estimation and subtraction
- Gaussian / Voigt line fitting
- Integration with TNS / WISeREP APIs for direct spectrum download
- Jupyter notebook widget mode

---

## Contributors

- **Harsh Kumar** (primary contributor) -- [Harvard Center for Astrophysics](https://www.cfa.harvard.edu/) / [IAIFI](https://iaifi.org/)

---

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

---

## Acknowledgments

- Inspired by the spectral plotting tools on [TNS](https://www.wis-tns.org/) and [WISeREP](https://www.wiserep.org/)
- Built with [Dash](https://dash.plotly.com/) by Plotly
- Spectral line wavelengths from [NIST Atomic Spectra Database](https://www.nist.gov/pml/atomic-spectra-database) and published supernova spectroscopy literature
