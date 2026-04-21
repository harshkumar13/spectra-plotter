# SpectraPlotter -- Detailed User Guide

This guide walks you through every feature of SpectraPlotter, with practical examples for transient astronomy workflows.

---

## 1. Installation

### 1.1 System requirements

- **Python 3.8+** (tested on 3.8, 3.9, 3.10, 3.11, 3.12)
- **OS**: macOS, Linux, Windows
- **Browser**: Chrome, Firefox, Safari, or Edge (any modern browser)
- **Disk**: ~50 MB for dependencies

### 1.2 Install steps

```bash
# Clone
git clone https://github.com/harshkumar13/spectra-plotter.git
cd spectra-plotter

# Install (editable mode -- changes to code take effect immediately)
pip install -e .

# Verify
python -c "from spectra_plotter import __version__; print(__version__)"
# Should print: 1.0.0
```

If you prefer not to install as a package:
```bash
pip install -r requirements.txt
python -m spectra_plotter
```

---

## 2. Launching the App

### 2.1 Basic launch

```bash
spectra-plotter
```

Opens at **http://127.0.0.1:8050**. Your default browser should open automatically; if not, copy-paste the URL.

### 2.2 Point to a data folder

```bash
spectra-plotter --folder /Users/you/data/SN2024abc/spectra/
```

When you launch with `--folder`, all recognized spectrum files in that directory are immediately available in the dropdown -- no need to upload or scan.

### 2.3 Network access

To let collaborators on the same network access your instance:
```bash
spectra-plotter --host 0.0.0.0 --port 8080
```

They can then open `http://your-ip:8080` in their browser.

### 2.4 Debug mode

```bash
spectra-plotter --debug
```

Enables auto-reload when you edit the source code. Useful for development.

---

## 3. Interface Layout

The interface is organized as a compact, WISeREP-inspired layout:

### 3.1 Top bar (data loading)

```
[ Upload / Local folder ]  [ File selector dropdown ........... ]  [Select all] [Clear] [Export CSV] [Reload]
```

- **Upload / Local folder** radio: switches between file upload and folder scanning
- **File selector**: multi-select dropdown listing available spectra
- **Buttons**: Select all, Clear selection, Export CSV, Reload file list

### 3.2 Main area (side-by-side)

```
+---------------------------+---------------------------+
|                           |  Line Identification      |
|   Interactive Plot        |  v-shift: [0] km/s        |
|   (Plotly)                |  Telluric: [x] O2 [x] H2O|
|                           |  ----------------------   |
|                           |  [x] H (Balmer)           |
|                           |     [ ] H-alpha            |
|                           |     [ ] H-beta             |
|                           |  [x] He I                  |
|                           |     [x] He I (5876)        |
|                           |  ...                       |
+---------------------------+---------------------------+
```

**Left**: the spectrum plot.
**Right**: the line identification panel with checkboxes.

### 3.3 Toolbar (below plot)

```
[Zoom Full] [Auto Zoom]  Binning: [1]  Smooth: [None v] [11]  |  Norm: [Common v]  [ ] CR clip  [ ] Log Y
```

All processing controls in one compact row.

### 3.4 Redshift inputs (below toolbar)

```
z= [0.072] SN2024abc_20240115.fits    z= [0.072] SN2024abc_20240124.fits
```

One `z=` input per selected spectrum, shown inline.

### 3.5 Status bar (bottom)

```
2 spectra | SN2024abc_20240115.fits z=0.0720 SNR~45.3 (3400-9200A, 4521pts) | ...
```

Monospace summary of loaded data.

---

## 4. Loading Spectra

### 4.1 Upload mode (default)

1. Make sure "Upload" is selected in the top-left radio.
2. Drag-and-drop one or more files onto the upload area.
3. Files are copied to the internal `uploads/` folder and appear in the dropdown.
4. Select them from the dropdown.

### 4.2 Local folder mode

1. Switch the radio to "Local folder".
2. Enter the full path, e.g.: `/Users/you/Desktop/spectra/`
3. Click **Scan**.
4. All `.fits`, `.dat`, `.txt`, `.csv`, `.ascii` files appear in the dropdown.
5. Files are read directly from disk -- nothing is copied.

**Tip**: If you add new files to the folder while the app is running, click **Reload** to refresh the file list.

### 4.3 Select all / clear

- **Select all**: selects every file in the current dropdown.
- **Clear**: deselects everything (the plot goes blank).

---

## 5. Supported File Formats (detailed)

### 5.1 FITS files

SpectraPlotter tries multiple strategies to extract wavelength and flux:

**Strategy 1: Binary table extensions**
- Scans all HDU extensions (starting from the last) for table data.
- Looks for wavelength columns: `OPT_WAVE`, `wave`, `WAVE`, `wavelength`, `lam`, `lambda`, `wl`, `obs_wave`, `rest_wave`, `angstrom` (case-insensitive, partial match supported).
- Looks for flux columns: `OPT_FLUX`, `OPT_FLAM`, `flux`, `FLUX`, `flam`, `counts`, `intensity`, `obs_flux`, `f_lam`, `f_lambda`.

**Strategy 2: 1D array + WCS**
- A single 1D array in any HDU, with WCS header keywords:
  - `CRVAL1` -- reference wavelength
  - `CDELT1` (or `CD1_1`) -- wavelength step
  - `CRPIX1` -- reference pixel (default 1)
  - `DC-FLAG=1` -- log-lambda mode (wavelength = 10^(CRVAL1 + ...))

**Strategy 3: 2D array**
- Shape `[N, M]` with N >= 2: row 0 = wavelength, row 1 = flux.

**Strategy 4: Squeezed 3D+**
- Higher-dimensional arrays are squeezed, then treated as 1D or 2D.

### 5.2 ASCII / text files

The parser:
1. Reads all lines from the file.
2. Strips lines starting with `#` or `%` (comments).
3. Strips lines that are purely alphabetic (column headers).
4. Auto-detects delimiter: comma, tab, or whitespace.
5. Parses each remaining line into numbers. Lines with fewer than 2 numbers are skipped.
6. Column 1 = wavelength, Column 2 = flux.

**If all of the above fail**, it falls back to `pandas.read_csv` with flexible delimiter detection.

### 5.3 Example file formats that work

**Whitespace-delimited (.dat)**
```
# Object: SN2024abc
# Date: 2024-01-15
# Instrument: EFOSC2
3500.0  1.23e-16
3501.5  1.25e-16
3503.0  1.22e-16
```

**CSV with header (.csv)**
```
wavelength,flux
3500.0,1.23e-16
3501.5,1.25e-16
```

**Tab-delimited (.txt)**
```
3500.0	1.23e-16
3501.5	1.25e-16
```

**Comment-heavy (.ascii)**
```
% Spectrum of SN 2024abc
% Observed 2024-01-15 with VLT/X-shooter
% wavelength(A)  flux(erg/s/cm2/A)
3500.0  1.23e-16
3501.5  1.25e-16
```

---

## 6. Redshift Correction

### 6.1 Setting redshifts

When spectra are selected, a `z=` input appears for each one below the toolbar. Type the redshift value. The plot updates immediately.

### 6.2 How it works

Wavelengths are corrected to the rest frame:
```
lambda_rest = lambda_observed / (1 + z)
```

Telluric masking is always applied in the **observed** frame (before redshift correction), since atmospheric absorption depends on the observed wavelength.

### 6.3 Session persistence

Redshifts are stored in your browser session. If you close the tab and reopen it, they may persist (depending on your browser settings). They are always preserved while the app is running.

---

## 7. Spectral Line Identification

### 7.1 The line ID panel

The right side of the interface shows element groups in a 4-column grid:

- **Column 1**: H (Balmer/Paschen/Brackett), He I, He II, C I, C II, N II
- **Column 2**: O I, O II, [O I], [O II], [O III], Si I, Si II, S I, S II, [S III]
- **Column 3**: Na I, Mg I, Mg I], Mg II, Ca II, [Ca II], Fe II, [Fe II], Fe III, Ni II, [Ni II]
- **Column 4**: [Ar III], Sc II, Ba II, Ti II, Cr II, Telluric, Host Galaxy, SLSN-I phases

### 7.2 Selecting lines

**Toggle an entire group**: Click the checkbox next to the element name (e.g., the checkbox next to "H (Balmer)"). This selects or deselects all lines in that group.

**Toggle individual lines**: Click individual line checkboxes (e.g., just H-alpha).

**Clear everything**: Click **Clear all** at the top of the panel.

### 7.3 How lines appear on the plot

Selected lines appear as **dotted vertical lines** at their rest-frame wavelengths. Each element group has a distinct color. The line label (e.g., "H-alpha") is shown as a small rotated annotation above the line.

### 7.4 Velocity shift

Enter a value in the **v-shift** field (km/s). All line markers shift by:
```
lambda_shifted = lambda_rest * (1 + v / c)
```
where c = 300,000 km/s.

This is useful for:
- Measuring photospheric expansion velocities in SN spectra
- Aligning absorption minima with rest-frame lines
- Checking if a feature matches a line at a particular velocity

### 7.5 Common workflows

**SN classification (Type Ia):**
1. Set redshift
2. Check **Si II** -- look for the 6355 A absorption
3. Check **Ca II** -- look for the H&K (3934/3969) and NIR triplet
4. Check **Fe II** -- multiple features in the blue

**SN classification (Type II):**
1. Check **H (Balmer)** -- strong H-alpha P-Cygni profile
2. Check **He I** for Type IIb
3. Check **Fe II** for photospheric features

**SLSN-I analysis:**
1. At early times: check **SLSN-I Pre-peak** for O II lines
2. Near peak: check **SLSN-I Peak** for O II + Fe III + Fe II
3. Late phase: check **SLSN-I Late** for [O I], [Ca II], Mg I]
4. Nebular: check **SLSN-I Nebular** for [O II]

**Host galaxy redshift:**
1. Check **Host Galaxy** group
2. Look for narrow emission lines (H-alpha, [O III] 5007, [O II] 3727)
3. Adjust redshift until the lines align

---

## 8. Processing Controls

### 8.1 Binning

Set the **Binning** field to N (integer). Each N adjacent pixels are averaged:
- wavelength_binned = mean of N wavelength values
- flux_binned = mean of N flux values

Set to **1** for no binning.

**When to use**: Low signal-to-noise spectra, especially NIR data.

### 8.2 Smoothing

Choose from the **Smooth** dropdown:

| Mode | Parameter meaning | Effect |
|------|------------------|--------|
| **None** | -- | Raw data |
| **SavGol** | Window size (odd integer, e.g., 11, 21, 51) | Savitzky-Golay filter (polynomial order 3). Preserves line shapes better than Gaussian. |
| **Gauss** | Window size (used to derive sigma = window/2) | Gaussian convolution. More aggressive smoothing. |

**When to use**:
- SavGol: general-purpose, good for preserving sharp features
- Gauss: when you want broad smoothing for continuum estimation

### 8.3 Normalization

| Mode | What it does |
|------|-------------|
| **Off** | Plot raw flux values |
| **Each** | Each spectrum is divided by its own maximum flux |
| **Common** | All spectra are normalized to the maximum flux in their common (overlapping) wavelength range |

**Common** is the default and most useful for comparing multiple spectra on the same scale.

### 8.4 Cosmic ray removal

Check **CR clip** to enable conservative cosmic ray removal:
1. A median-filtered baseline is computed (kernel size 9)
2. Residuals are calculated
3. Points exceeding 8 * MAD (median absolute deviation) above baseline are replaced with the median value

This is conservative by design -- it only clips extreme outliers.

### 8.5 Log Y

Check **Log Y** to switch the y-axis to logarithmic scale. Non-positive flux values are masked (set to NaN).

---

## 9. Telluric Masking

The toolbar includes checkboxes for atmospheric absorption bands. When checked, data in those wavelength ranges is set to NaN (appears as a gap in the plot).

| Band | Wavelength range |
|------|-----------------|
| O2 B-band | 6865--6969 A |
| O2 A-band | 7600--7700 A |
| H2O | 9300--9600 A |
| H2O | 11300--11600 A |
| H2O | 13500--15000 A |
| H2O | 18140--19500 A |
| H2O | 24500--27000 A |

**Important**: Masking is applied in the **observed** frame, before redshift correction. This is correct because telluric absorption is an artifact of Earth's atmosphere, not of the source.

---

## 10. Exporting Data

Click **Export CSV** in the top bar. The output CSV has three columns:

| Column | Description |
|--------|-------------|
| `filename` | Original filename |
| `wavelength_rest_A` | Rest-frame wavelength in Angstroms (redshift-corrected) |
| `flux` | Flux value |

The export includes all currently selected spectra. Data is cleaned (NaNs removed) but no smoothing/binning is applied to the export -- you get the redshift-corrected raw data.

---

## 11. Keyboard and Mouse Controls

### Plot interactions

| Action | How |
|--------|-----|
| Zoom in | Click-drag a rectangle |
| Scroll zoom | Mouse wheel |
| Pan | Shift + drag |
| Reset zoom | Double-click on plot, or click **Zoom Full** |
| Hover info | Move mouse over any point |
| Save PNG | Camera icon in Plotly modebar (top-right) |
| Draw line | Line-draw tool in modebar |
| Erase drawings | Eraser tool in modebar |

---

## 12. Troubleshooting

### "No readable spectra" error

- Check that your file has at least two columns (wavelength, flux)
- Check that wavelength values are positive (in Angstroms)
- For FITS files, check that the header has correct WCS keywords or that the binary table has recognizable column names
- Check the status bar for specific error messages

### Plot shows no data / flat line

- Check that the redshift is correct (a wrong z can shift data out of the visible range)
- Click **Zoom Full** to reset the view
- Check if telluric masking is removing too much data

### Lines don't align with features

- Verify the redshift value
- Try adjusting the velocity shift
- Remember that absorption minima are blueshifted from the rest wavelength

### App won't start

- Make sure all dependencies are installed: `pip install -r requirements.txt`
- Check for port conflicts: `spectra-plotter --port 8051`
- Try: `python -c "from spectra_plotter.app import create_app; print('OK')"`

---

## 13. Processing Functions API

For users who want to use the processing functions in scripts or notebooks:

### Reading spectra

```python
from spectra_plotter.io import read_spectrum, list_spectra_files

# List all spectra in a folder
files = list_spectra_files("/path/to/folder")

# Read a single file (auto-detects format)
wave, flux = read_spectrum("spectrum.fits")
wave, flux = read_spectrum("spectrum.dat")
wave, flux = read_spectrum("spectrum.csv")
```

### Cleaning

```python
from spectra_plotter.processing import clean_spectrum

# Remove NaNs, sort by wavelength, remove duplicates, remove negative wavelengths
wave, flux = clean_spectrum(wave, flux)
```

### Smoothing

```python
from spectra_plotter.processing import apply_savgol, apply_gaussian_smooth

# Savitzky-Golay (window must be odd)
flux_smooth = apply_savgol(flux, window=21, polyorder=3)

# Gaussian
flux_smooth = apply_gaussian_smooth(flux, sigma=5.0)
```

### Binning

```python
from spectra_plotter.processing import apply_binning

wave_bin, flux_bin = apply_binning(wave, flux, bin_size=5)
```

### Cosmic ray removal

```python
from spectra_plotter.processing import cosmic_ray_clip

flux_clean = cosmic_ray_clip(flux, kernel=9, sigma=8.0)
```

### Telluric masking

```python
from spectra_plotter.processing import apply_telluric_mask

flux_masked = apply_telluric_mask(wave_observed, flux, ["O2_A", "H2O_13500"])
```

### Signal-to-noise

```python
from spectra_plotter.processing import compute_snr

snr = compute_snr(flux, window=50)
```

### Equivalent width

```python
from spectra_plotter.processing import measure_ew

# Measure EW of H-alpha (positive = absorption, negative = emission)
ew = measure_ew(wave, flux, line_center=6563.0, width=20.0)
```
