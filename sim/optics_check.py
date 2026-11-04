"""Final check: the aerial image computed a second, independent way.

simulator_v2.aerial_image rasterises the mask onto 1 nm pixels and uses FFTs.
Here the image is built from the exact Fourier series of ideal square mask
holes (each hole's transform is m^2 sinc(u m) sinc(v m) times a phase), summing
only the diffraction orders each source point lets through the pupil.  No
pixels, no FFT.  The two should agree to within the small pixel-averaging
effect (sinc of the pixel at the frequencies the lens passes, ~0.3%).

Usage: python optics_check.py  ->  ../results/optics_check.json
"""
import json
import os

import numpy as np

from simulator_v2 import Params, aerial_image, _hole_centres, _tile_shape, _source_points, WAVELENGTH_NM

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def analytic_image(p: Params) -> np.ndarray:
    H, W = _tile_shape(p)
    Lx, Ly = W * p.px, H * p.px
    m = p.mask_cd if p.mask_cd is not None else p.cd_target
    cx, cy = _hole_centres(p)
    cutoff = p.na / WAVELENGTH_NM
    nmax_x = int(np.ceil(2 * cutoff * Lx)) + 1
    nmax_y = int(np.ceil(2 * cutoff * Ly)) + 1
    n = np.arange(-nmax_x, nmax_x + 1)
    k = np.arange(-nmax_y, nmax_y + 1)
    u, v = np.meshgrid(n / Lx, k / Ly)                       # tile harmonics, cycles/nm
    coef = (m * m / (Lx * Ly)) * np.sinc(u * m) * np.sinc(v * m) \
        * np.exp(-2j * np.pi * (u[..., None] * cx + v[..., None] * cy)).sum(-1)
    y = (np.arange(H) + 0.5) * p.px
    x = (np.arange(W) + 0.5) * p.px
    ex = np.exp(2j * np.pi * np.outer(n / Lx, x))             # (orders_x, W)
    ey = np.exp(2j * np.pi * np.outer(k / Ly, y))             # (orders_y, H)
    pts, w = _source_points(p)
    img = np.zeros((H, W))
    for (sx, sy), wk in zip(pts, w):
        passed = ((u + sx * cutoff) ** 2 + (v + sy * cutoff) ** 2) <= cutoff ** 2
        c = np.where(passed, coef, 0)
        field = ey.T @ c @ ex                                  # (H, W)
        img += wk * np.abs(field) ** 2
    return (1 - p.flare) * img + p.flare


CASES = {
    "P32_CD16_quasar": Params(cd_target=16, pitch_x=32, pitch_y=32),
    "P40_CD20_quasar": Params(cd_target=20, pitch_x=40, pitch_y=40),
    "P48_CD24_mask26_quasar": Params(cd_target=24, mask_cd=26, pitch_x=48, pitch_y=48),
    "P64_CD30_quasar": Params(cd_target=30, pitch_x=64, pitch_y=64),
    "P50_CD23_mask24_small_annular": Params(cd_target=23, mask_cd=24, pitch_x=50, pitch_y=50,
                                            illumination="annular", sigma_in=0.292, sigma_out=0.643),
    "P52_mask31_NA025_quasar": Params(cd_target=26, mask_cd=31, pitch_x=52, pitch_y=52, na=0.25),
    "hex_46x40_hexapole": Params(cd_target=24, pitch_x=46, pitch_y=40, stagger=True, cells_x=3, cells_y=4,
                                 illumination="hexapole"),
    "P36_CD18_mask17.3_quasar": Params(cd_target=18, mask_cd=17.3, pitch_x=36, pitch_y=36),
}

if __name__ == "__main__":
    out = {}
    for name, p in CASES.items():
        a, b = aerial_image(p), analytic_image(p)
        peak = b.max()
        row = dict(peak_fft=float(a.max()), peak_analytic=float(peak),
                   max_abs_diff_rel_peak=float(np.abs(a - b).max() / peak),
                   mean_ratio_in_holes=float((a[b > 0.5 * peak] / b[b > 0.5 * peak]).mean()),
                   min_fft=float(a.min()), min_analytic=float(b.min()))
        out[name] = row
        print(f"{name:32s} peak {row['peak_fft']:.4f} vs {row['peak_analytic']:.4f}   "
              f"max |diff| {100 * row['max_abs_diff_rel_peak']:.2f}% of peak   "
              f"min {row['min_fft']:.4f} vs {row['min_analytic']:.4f}")
    json.dump(out, open(RESULTS + "/optics_check.json", "w"), indent=1)
