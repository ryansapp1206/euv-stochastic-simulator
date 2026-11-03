"""Sanity tests for simulator v2. Run: python sim/test_simulator_v2.py"""
from dataclasses import replace

import numpy as np

from simulator_v2 import (Params, PHOTONS_PER_NM2_PER_MJCM2, aerial_image, simulate,
                          summarize, dose_to_size)


def test_photon_constant():
    # 92 eV photons: 1 mJ/cm^2 = 1e-17 J/nm^2 -> ~0.680 photons/nm^2
    assert abs(PHOTONS_PER_NM2_PER_MJCM2 - 0.680) < 0.002


def test_open_frame_normalisation():
    # a mask that is all hole must give intensity 1 everywhere
    p = Params(cd_target=40, mask_cd=40, pitch_x=40, pitch_y=40)
    I = aerial_image(p)
    assert np.allclose(I, 1.0, atol=1e-6), (I.min(), I.max())


def test_image_symmetry_and_range():
    # a square array under symmetric quasar illumination images identically in x and y,
    # every hole looks the same, and a dark-field hole never exceeds the open-frame dose
    p = Params()
    I = aerial_image(p)
    assert I.min() >= 0 and I.max() < 1
    assert np.allclose(I, I.T, atol=1e-6)
    assert np.allclose(I, np.roll(I, p.pitch_x, axis=1), atol=1e-6)


def test_photon_poisson_statistics():
    # absorbed photons per column: mean and variance both = lambda
    rng = np.random.default_rng(0)
    lam = 30 * PHOTONS_PER_NM2_PER_MJCM2 * Params().absorbed_frac
    x = rng.poisson(lam, 200_000)
    assert abs(x.mean() - lam) / lam < 0.01 and abs(x.var() - lam) / lam < 0.02


def test_noise_off_gives_identical_holes():
    p = dose_to_size(Params())
    s = summarize(simulate(replace(p, noise=()), 20, np.random.default_rng(1)))
    assert s["sigma_cd"] < 1e-9, s


def test_each_noise_source_adds_variation():
    p = dose_to_size(Params())
    full = summarize(simulate(p, 200, np.random.default_rng(2)))["sigma_cd"]
    for src in ("photon", "pag", "quencher", "pg"):
        only = summarize(simulate(replace(p, noise=(src,)), 200, np.random.default_rng(2)))["sigma_cd"]
        assert 0 < only < full, (src, only, full)


def test_acid_noise_budget_mack2019():
    # Mack 2019: at an acid yield of 1 (low conversion), acid shot noise equals photon
    # shot noise, so total acid variance is twice the photon-only variance.
    rng = np.random.default_rng(3)
    mu = 50.0                                 # absorbed photons in a volume
    n = 400_000
    photons = rng.poisson(mu, n)
    acids_photon_only = photons * 1.0         # QE = 1, deterministic conversion
    acids_both = rng.poisson(photons * 1.0)   # QE = 1, Poisson conversion
    ratio = acids_both.var() / acids_photon_only.var()
    assert 1.9 < ratio < 2.1, ratio


def test_more_dose_less_variation():
    a = dose_to_size(Params(dose=20))
    b = dose_to_size(Params(dose=60))
    sa = summarize(simulate(a, 200, np.random.default_rng(4)))["sigma_cd"]
    sb = summarize(simulate(b, 200, np.random.default_rng(4)))["sigma_cd"]
    assert sb < sa, (sa, sb)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
