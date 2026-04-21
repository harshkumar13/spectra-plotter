"""Tests for spectra_plotter.processing module."""

import numpy as np
import pytest

from spectra_plotter.processing import (
    clean_spectrum,
    apply_telluric_mask,
    cosmic_ray_clip,
    apply_savgol,
    apply_binning,
    apply_gaussian_smooth,
    compute_snr,
    measure_ew,
)


class TestCleanSpectrum:
    def test_removes_nans(self):
        wave = np.array([1.0, 2.0, np.nan, 4.0])
        flux = np.array([10.0, 20.0, 30.0, 40.0])
        w, f = clean_spectrum(wave, flux)
        assert len(w) == 3
        assert np.nan not in w

    def test_removes_negative_wavelengths(self):
        wave = np.array([-1.0, 0.0, 1.0, 2.0])
        flux = np.array([10.0, 20.0, 30.0, 40.0])
        w, f = clean_spectrum(wave, flux)
        assert w.min() > 0

    def test_sorts_by_wavelength(self):
        wave = np.array([3.0, 1.0, 2.0])
        flux = np.array([30.0, 10.0, 20.0])
        w, f = clean_spectrum(wave, flux)
        assert np.all(np.diff(w) > 0)
        np.testing.assert_array_equal(f, [10.0, 20.0, 30.0])

    def test_removes_duplicates(self):
        wave = np.array([1.0, 1.0, 2.0, 3.0])
        flux = np.array([10.0, 11.0, 20.0, 30.0])
        w, f = clean_spectrum(wave, flux)
        assert len(w) == 3

    def test_empty_input(self):
        w, f = clean_spectrum(np.array([]), np.array([]))
        assert len(w) == 0


class TestTelluricMask:
    def test_masks_band(self):
        wave = np.arange(7500, 7800, 1.0)
        flux = np.ones_like(wave)
        masked = apply_telluric_mask(wave, flux, ["O2_A"])
        in_band = (wave >= 7600) & (wave <= 7700)
        assert np.all(np.isnan(masked[in_band]))
        assert np.all(np.isfinite(masked[~in_band]))

    def test_no_bands(self):
        wave = np.arange(7500, 7800, 1.0)
        flux = np.ones_like(wave)
        result = apply_telluric_mask(wave, flux, [])
        np.testing.assert_array_equal(result, flux)


class TestCosmicRayClip:
    def test_removes_spike(self):
        rng = np.random.default_rng(42)
        flux = np.ones(100) + rng.normal(0, 0.01, 100)
        flux[50] = 1000.0
        clipped = cosmic_ray_clip(flux)
        assert clipped[50] < 100.0

    def test_preserves_smooth_data(self):
        flux = np.sin(np.linspace(0, 2 * np.pi, 100)) + 5.0
        clipped = cosmic_ray_clip(flux)
        np.testing.assert_allclose(clipped, flux, atol=0.5)


class TestSavgol:
    def test_smooths_noisy_data(self):
        rng = np.random.default_rng(42)
        flux = np.sin(np.linspace(0, 4 * np.pi, 200)) + rng.normal(0, 0.3, 200)
        smoothed = apply_savgol(flux, window=21, polyorder=3)
        assert np.nanstd(smoothed) < np.nanstd(flux)

    def test_small_window_returns_original(self):
        flux = np.ones(10)
        result = apply_savgol(flux, window=1, polyorder=1)
        np.testing.assert_array_equal(result, flux)


class TestBinning:
    def test_reduces_size(self):
        wave = np.arange(100.0)
        flux = np.ones(100)
        w, f = apply_binning(wave, flux, 5)
        assert len(w) == 20

    def test_bin_size_one_returns_original(self):
        wave = np.arange(10.0)
        flux = np.ones(10)
        w, f = apply_binning(wave, flux, 1)
        assert len(w) == 10


class TestGaussianSmooth:
    def test_smooths_data(self):
        rng = np.random.default_rng(42)
        flux = np.ones(200) + rng.normal(0, 0.5, 200)
        smoothed = apply_gaussian_smooth(flux, sigma=5)
        assert np.nanstd(smoothed) < np.nanstd(flux)


class TestSNR:
    def test_high_snr(self):
        flux = np.ones(200) * 100.0 + np.random.default_rng(42).normal(0, 0.1, 200)
        snr = compute_snr(flux)
        assert snr > 10

    def test_low_snr(self):
        flux = np.random.default_rng(42).normal(0, 10, 200)
        snr = compute_snr(flux)
        assert snr < 5


class TestMeasureEW:
    def test_emission_line(self):
        wave = np.linspace(6540, 6590, 200)
        continuum = np.ones_like(wave)
        line = continuum + 2.0 * np.exp(-0.5 * ((wave - 6563) / 2) ** 2)
        ew = measure_ew(wave, line, 6563.0, width=20.0)
        assert ew < 0  # emission = negative EW

    def test_absorption_line(self):
        wave = np.linspace(6540, 6590, 200)
        continuum = np.ones_like(wave)
        line = continuum - 0.5 * np.exp(-0.5 * ((wave - 6563) / 2) ** 2)
        ew = measure_ew(wave, line, 6563.0, width=20.0)
        assert ew > 0  # absorption = positive EW
