# Progress log

## 2026-08-06 (later) — "Redshift still broken" report resolved; legacy trap removed

- Harsh reported redshift still didn't work after the fix. Verified the packaged
  app end-to-end in a **real Chrome browser** (headless + DevTools protocol,
  `scratchpad cdp_test`): typing z=0.1 rest-frames the traces (3477.56 → 3161.42 Å
  = /1.1), the axis label flips Observed→Rest, the He I name-button draws all 6
  markers, and z persists across reloads (session persistence). The app is correct.
- Root cause of the report: the **legacy standalone script at the repo root**
  (`Plot_spectra_local_app_UPDATED.py`) still had the original prop_id-split bug —
  launching it instead of the package reproduces "redshift does nothing".
  That file is now a thin shim that launches the canonical package app; the old
  code lives in git history (commit 984acdb).
- Added `--window` flag (`spectra-plotter --window`): opens the app in a local
  app-style window (Chrome `--app` mode, no tabs/URL bar; falls back to the
  default browser). Requested as a "local window app" alternative to a webpage.

## 2026-08-06 — Frontend redesign + redshift/toggle fixes

**Done:**
- **Fixed the redshift bug**: entering z never updated the plot. Root cause in
  `update_redshift_store` (app.py): the callback `prop_id` was split on the *first*
  `.`, but pattern-matching IDs embed the file path (which contains dots, e.g.
  `spec.fits`), so the JSON id was truncated and parsing failed silently. The store
  is now written directly from the displayed input values (no prop_id parsing).
- **Dynamic x-axis label**: "Observed Wavelength (Å)" when all z=0, "Rest
  Wavelength (Å)" once any nonzero z is applied. Changing z (or the file set) now
  refits the view; other tweaks (smoothing, lines, …) preserve the user's zoom
  (uirevision keyed on files+redshifts).
- **One-click element groups**: the element name itself is now a button — click
  toggles the whole group on/off (all-on if not everything selected, else all-off).
  Individual per-line checkboxes remain. "Clear all" empties every group.
- **Full dark redesign** ("night-vision observatory console"): space-black base,
  warm amber accent, monospace instrument eyebrow labels, optical-spectrum gradient
  strip under the header as the signature. All styling in
  `spectra_plotter/assets/style.css`; layout uses CSS classes instead of inline styles.
- **Line markers rebuilt**: full-height vline shapes (`yref="y domain"`) instead of
  0→ymax traces — log-scale safe, don't pollute autoscale; rotated colored labels
  with hoverable exact wavelength.
- **Physics fixes on markers**: v-shift no longer applied to Telluric or Host
  Galaxy groups; Telluric markers (observed-frame) are mapped onto the rest-frame
  axis via the first spectrum's z (noted in the marker hover).
- **New controls**: "Apply to all" redshift quick-set; per-trace vertical offset
  (Display → offset) for stacking epochs; robust-percentile Auto zoom (0.5–99.7%);
  status bar is now color-coded chips matched to trace colors; z inputs debounced
  (update on Enter/blur, not per keystroke).
- Cosmetics: telluric-mask checkboxes now show band wavelengths (H2O 9300 etc.);
  brightened illegible catalog colors (Ti II, Cr II, [Fe II]); removed the unused
  `savgol-polyorder-store`.

**Verified:** 18/18 unit tests pass; all key callbacks exercised end-to-end against
the running server via `/_dash-update-component` (redshift store with dotted paths,
rest-framing + label switch, group toggle on/off, clear-all, marker shapes/log
scale/telluric z-mapping); UI + rendered figure inspected via headless Chrome
screenshots. Note: headless Chrome without GPU can't render Scattergl (WebGL) —
fine in real browsers.

**Next steps (ideas, not started):**
- README screenshots are stale (input/*.png show the old UI) — retake.
- Possible features: template cross-correlation, continuum subtraction, line EW
  measurement from the UI (measure_ew exists in processing.py but is unused).
