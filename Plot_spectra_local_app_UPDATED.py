"""LEGACY LAUNCHER — the standalone app that used to live in this file was
superseded by the `spectra_plotter/` package (which has the redesigned UI and
the redshift/line-toggle fixes; see PROGRESS.md).

This file is kept only so old launch habits keep working: running it now
starts the canonical package app. The old standalone code is preserved in git
history (commit 984acdb and earlier).

Usage (all equivalent):
    python Plot_spectra_local_app_UPDATED.py
    python -m spectra_plotter
    spectra-plotter [--port PORT] [--folder DIR] [--window]
"""

from spectra_plotter.app import main

if __name__ == "__main__":
    main()
