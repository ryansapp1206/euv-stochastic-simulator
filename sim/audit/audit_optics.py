"""Independent optics audit: analytic diffraction-order image vs simulator FFT image,
pole-clipping footprint in the dataset, and sensitivity of LCDU to the source shape."""
import sys
from dataclasses import replace
import numpy as np
import pandas as pd

import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import simulator_v2 as S  # noqa: E402

LAM, NA = 13.5, 0.33

# ---------------------------------------------------------------- A. analytic image
def analytic_image(p: S.Params, pixel_sinc: bool):
    """Image from the closed-form Fourier series of a square hole, no FFT, no pixel mask."""
    pitch, a = p.pitch_x, (p.mask_cd if p.mask_cd is not None else p.cd_target)
    pts, w = S._source_points(p)
    H, W = S._tile_shape(p)
    # image period is one pitch; evaluate over the tile at pixel centres
    yy, xx = (np.mgrid[:H, :W] + 0.5) * p.px
    img = np.zeros((H, W))
    mmax = 6
    ms = np.arange(-mmax, mmax + 1)
    cutoff = NA / LAM
    for (sx, sy), wk in zip(pts, w):
        field = np.zeros((H, W), dtype=complex)
        for m in ms:
            for n in ms:
                fx, fy = m / pitch, n / pitch
                if (fx / cutoff + sx) ** 2 + (fy / cutoff + sy) ** 2 > 1.0 + 1e-12:
                    continue
                c = (a / pitch) ** 2 * np.sinc(m * a / pitch) * np.sinc(n * a / pitch) * np.exp(-1j * np.pi * (m + n))
                if pixel_sinc:                                   # box-averaging by the 1 nm pixel
                    c *= np.sinc(fx * p.px) * np.sinc(fy * p.px)
                field += c * np.exp(2j * np.pi * (fx * xx + fy * yy))
        img += wk * np.abs(field) ** 2
    return img


print("=" * 96)
print("A. ANALYTIC DIFFRACTION-ORDER IMAGE  vs  SIMULATOR (FFT of anti-aliased mask)")
print("=" * 96)
print(f"{'case':<22}{'sim peak':>10}{'analytic':>10}{'(+pixel)':>10}{'max |diff|':>12}{'rel to peak':>13}")
for pitch, cd in ((64, 30.0), (48, 24.0), (40, 20.0), (36, 18.0), (32, 16.0)):
    p = S.Params(pitch_x=pitch, pitch_y=pitch, cd_target=cd, mask_cd=cd)
    sim = S.aerial_image(p)
    a0 = analytic_image(p, False)
    a1 = analytic_image(p, True)
    d = np.abs(sim - a1).max()
    print(f"P{pitch} CD{cd:<5.0f}{'':<10}{sim.max():10.4f}{a0.max():10.4f}{a1.max():10.4f}{d:12.2e}{d/sim.max():13.2%}")
print("  -> 'max |diff|' compares sim to the analytic image with the 1 nm pixel box-average included.")

# ---------------------------------------------------------------- B. clipping footprint
print()
print("=" * 96)
print("B. QUASAR POLE CLIPPING (pole centre limited to sigma <= 0.9 - pole radius = 0.8)")
print("=" * 96)
p_clip = LAM / (np.sqrt(2) * NA * 0.8)
print(f"   ideal pole position exceeds 0.8 for pitch < {p_clip:.2f} nm")
R = pd.read_csv(os.environ.get("EUV_ROWS", "data/full/rows.csv"))
nom = R[np.isclose(R.dose_factor, 1.0)]
res = R.drop_duplicates("resist_id")
print(f"   accepted resists with pitch < {p_clip:.1f} nm: {(res.pitch_nm < p_clip).sum()} of {len(res)} "
      f"({100*(res.pitch_nm < p_clip).mean():.1f}%)")
bins = [31.5, 36.2, 40, 48, 56, 64.5]
nom = nom.assign(pb=pd.cut(nom.pitch_nm, bins))
g = nom.groupby("pb", observed=True).agg(n=("lcdu_3s_nm", "size"), lcdu_med=("lcdu_3s_nm", "median"),
                                          nils_med=("nils", "median"), cd_ratio=("cd_ratio", "median"))
print(g.to_string())

# diagonal-order capture for sources across the whole pole at P32, to see the actual image loss
for pitch in (32, 34, 36):
    cd = 0.5 * pitch
    p = S.Params(pitch_x=pitch, pitch_y=pitch, cd_target=cd, mask_cd=cd)
    pts, w = S._source_points(p)
    cut = NA / LAM
    f1 = 1.0 / pitch / cut
    n_diag = 0
    for sx, sy in pts:
        n_diag += (np.hypot(-f1 + sx, -f1 + sy) <= 1.0)
    pc = replace(p, pole_sigma=0.10)
    im = S.aerial_image(p)
    print(f"   P{pitch}: {n_diag}/{len(pts)} source points pass the (1,1) order; image peak {im.max():.3f}, "
          f"min {im.min():.4f}")

# ---------------------------------------------------------------- C. source-shape sensitivity
print()
print("=" * 96)
print("C. SENSITIVITY OF LCDU TO SOURCE SHAPE (P40 / 20 nm / 60 mJ/cm2, quencher re-fitted each time)")
print("=" * 96)
base = S.Params(pitch_x=40, pitch_y=40, cd_target=20.0, dose=60.0)
print(f"{'illumination':<34}{'quencher':>10}{'NILS-ish peak':>15}{'LCDU 3s':>10}{'rel':>8}")
ref = None
cases = [("quasar pole r=0.05", dict(pole_sigma=0.05, pole_samples=5)),
         ("quasar pole r=0.10 (model)", dict(pole_sigma=0.10)),
         ("quasar pole r=0.15", dict(pole_sigma=0.15)),
         ("quasar pole r=0.20", dict(pole_sigma=0.20)),
         ("annular 0.4-0.85", dict(illumination="annular"))]
for name, kw in cases:
    q = replace(base, **kw)
    q = S.dose_to_size(q, "rho_q", n_trials=100, rng_seed=3)
    s = S.summarize(S.simulate(q, 2500, np.random.default_rng(11)))
    pk = S.aerial_image(q).max()
    if ref is None and "model" in name:
        ref = s["lcdu_3s"]
    print(f"{name:<34}{q.rho_q:10.4f}{pk:15.3f}{s['lcdu_3s']:10.3f}")
print(f"   (reference: model quasar r=0.10 gives {ref:.3f} nm)")
