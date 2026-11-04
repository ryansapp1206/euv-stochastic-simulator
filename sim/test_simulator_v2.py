"""Sanity tests for simulator v2. Run: python sim/test_simulator_v2.py"""
from dataclasses import replace

import numpy as np

from simulator_v2 import (Params, PHOTONS_PER_NM2_PER_MJCM2, aerial_image, simulate,
                          summarize, dose_to_size, _geometry, _measure_components)


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


def _disc(clear, cy, cx, r):
    yy, xx = np.mgrid[:clear.shape[0], :clear.shape[1]] + 0.5
    clear |= (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r


def test_hole_measurement_ignores_specks_and_finds_bridges():
    # pass 2: a hole is its largest connected region; specks don't count; a real bridge merges
    p = Params(cd_target=10, pitch_x=20, pitch_y=20, cells_x=2, cells_y=2)
    _, cell, nh = _geometry(p)
    c = np.zeros((40, 40), bool)
    for cy, cx in ((10, 10), (10, 30), (30, 10), (30, 30)):
        _disc(c, cy, cx, 5)
    disc_area = c[:20, :20].sum()
    a = c.copy()
    a[9, 19] = a[9, 20] = True                 # two specks touching across the cell border
    area, merged, _ = _measure_components(a[None], cell, nh)
    assert not merged.any(), merged
    assert np.all(area == disc_area), area
    b = c.copy()
    b[9, 14:26] = True                         # a 1-pixel-wide bridge between holes 0 and 1
    area, merged, group = _measure_components(b[None], cell, nh)
    assert merged[0].tolist() == [True, True, False, False], merged
    assert group[0, 0] == group[0, 1] >= 0 and (group[0, 2:] == -1).all(), group   # one merge event
    d = c.copy()
    d[29, 34:40] = True; d[29, 0:7] = True     # bridge through the periodic tile edge (holes 3 and 2)
    area, merged, _ = _measure_components(d[None], cell, nh)
    assert merged[0].tolist() == [False, False, True, True], merged
    e = np.zeros((40, 40), bool)
    e[5:10, 5:10] = True; e[10:15, 10:15] = True   # corner-to-corner only: not a bridge
    area, merged, _ = _measure_components(e[None], cell, nh)
    assert not merged.any() and area[0, 0] == 25, (area, merged)


def test_mask_area_is_exact_for_any_size():
    # pass 2: anti-aliased mask, so a 17.17 nm or 31 nm mask hole has exactly that area
    for m, pitch in ((17.17, 32), (31.0, 52), (24.0, 50), (22.5, 50)):
        mask, _, nh = _geometry(Params(cd_target=m, mask_cd=m, pitch_x=pitch, pitch_y=pitch))
        assert abs(mask.sum() / nh - m * m) < 1e-6, (m, mask.sum() / nh)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)
