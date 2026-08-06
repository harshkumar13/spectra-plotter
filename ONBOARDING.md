# Onboarding

## What this is
SpectraPlotter: a local-first Dash web app for plotting astronomical (transient)
spectra with interactive line identification, WISeREP/TNS-inspired. Author:
Harsh Kumar (CfA/IAIFI).

## Canonical codebase
The `spectra_plotter/` package is canonical.
`Plot_spectra_local_app_UPDATED.py` at the repo root is now just a shim that
launches the package app (the old standalone code, which had the redshift bug,
is in git history). Do not develop there.

## Directory map
- `spectra_plotter/app.py` — Dash app: layout + all callbacks (the frontend).
- `spectra_plotter/io.py` — spectrum readers (FITS binary table/WCS/2D, ASCII, CSV).
- `spectra_plotter/processing.py` — clean/smooth/bin/CR-clip/SNR/EW functions.
- `spectra_plotter/line_catalog.py` — `SPECTRAL_LINES` (35+ groups), `LINE_COLORS`,
  `TELLURIC_BANDS`, `DEFAULT_TELLURIC_MASK`.
- `spectra_plotter/assets/style.css` — the entire theme (dark "observatory console",
  amber accent, spectrum-gradient header strip). Prefer CSS classes here over
  inline styles in app.py.
- `spectra_plotter/uploads/` — default upload/scan folder (also `uploads/` at repo
  root holds sample spectra for manual testing).
- `tests/test_processing.py` — unit tests (processing only; no callback tests).

## Run / test
```bash
python -m spectra_plotter            # http://127.0.0.1:8050
python -m spectra_plotter --window   # local app-style window (Chrome --app mode)
python -m spectra_plotter --port 8067 --folder /path/to/spectra
python -m pytest tests/ -q
```

## Conventions and gotchas
- **Pattern-matching callback IDs embed file paths, which contain dots.** Never
  parse `callback_context.triggered[0]["prop_id"]` with `split(".")` — the JSON id
  gets truncated. Use `rsplit(".", 1)`, `ctx.triggered_id`, or (as
  `update_redshift_store` does) read the values/ids lists directly. This exact bug
  once silently killed the redshift feature.
- Rest-framing is `λ_rest = λ_obs / (1+z)`; telluric *masking* is applied in the
  observed frame before rest-framing; telluric *markers* are drawn at
  `λ/(1+z_first_spectrum)`; v-shift applies only to non-Telluric, non-Host groups.
- Plot view state: `uirevision` is keyed on (files, redshifts) — changing either
  refits the view, everything else keeps the user's zoom. The relayout re-apply
  block in `update_plot` mirrors this; keep both in sync.
- Line markers are layout *shapes* (`yref="y domain"`) + annotations, not traces.
- Element group toggling is a stateless button (`element-all-btn`): all-on if the
  group isn't fully selected, else all-off. No stored toggle state to sync.
- Traces use `Scattergl` (WebGL): headless Chrome without GPU shows "WebGL is not
  supported" — that's a test-environment limitation, not a bug.
- Callbacks can be exercised without a browser by POSTing to
  `/_dash-update-component` (see PROGRESS.md 2026-08-06 for verified examples).
- Dash ≥2.18 / Plotly ≥6 assumed (Plotly 6 serializes arrays as base64 `bdata`).
