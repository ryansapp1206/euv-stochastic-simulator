"""Prototype stochastic EUV contact-hole printing simulator (teaching-grade).

Pipeline per trial, on a 1 nm grid:
  mask (disk of diameter d) -> optical blur (Gaussian) -> aerial image I(x,y)
  absorbed photons  ~ Poisson(dose_photons_per_nm2 * absorption * I)
  acids             ~ Poisson(acid_yield * photons)
  quencher          ~ Poisson(q) per pixel
  net = blur(acids) - blur(quencher)        (secondary electrons + acid diffusion)
  pixel dissolves if net > threshold       (positive-tone development)
Outputs per trial: printed CD (equivalent-circle diameter) and pass/fail.
"""
import numpy as np
from scipy.ndimage import gaussian_filter

PHOTON_EV = 1239.84 / 13.5               # ~91.8 eV per EUV photon
J_PER_EV = 1.602e-19
# 1 mJ/cm^2 = 1e-17 J/nm^2  -> photons per nm^2 per (mJ/cm^2)
PHOTONS_PER_NM2_PER_MJCM2 = 1e-17 / (PHOTON_EV * J_PER_EV)

G = 64                                    # grid size (nm)
yy, xx = np.mgrid[:G, :G] - (G - 1) / 2
RR = np.hypot(xx, yy)


def aerial(d, optical_sigma):
    mask = (RR <= d / 2).astype(float)
    return gaussian_filter(mask, optical_sigma)


def mean_net(dose, d, optical_sigma=4.0, absorption=0.2, acid_yield=2.0,
             q=1.0, blur=4.0):
    """Noise-free expected net acid map (used to calibrate the threshold)."""
    I = aerial(d, optical_sigma)
    lam = dose * PHOTONS_PER_NM2_PER_MJCM2 * absorption * I
    return gaussian_filter(acid_yield * lam, blur) - q


def simulate(dose, d, n, threshold, optical_sigma=4.0, absorption=0.2,
             acid_yield=2.0, q=1.0, blur=4.0, rng=None):
    rng = rng or np.random.default_rng()
    I = aerial(d, optical_sigma)
    lam = dose * PHOTONS_PER_NM2_PER_MJCM2 * absorption * I
    photons = rng.poisson(np.broadcast_to(lam, (n, G, G)))
    acids = rng.poisson(acid_yield * photons)
    quench = rng.poisson(q, size=(n, G, G))
    net = gaussian_filter((acids - quench).astype(np.float32), (0, blur, blur))
    open_ = net > threshold
    area = open_[:, RR <= d].sum(axis=1)          # only count near the hole
    cd = 2 * np.sqrt(area / np.pi)
    fail = area < 0.25 * np.pi * (d / 2) ** 2     # hole effectively missing
    return cd, fail


def calibrate_threshold(ref_dose=30.0, ref_d=20.0, **kw):
    """Resist threshold that prints a ref_d hole at ref_d with ref_dose."""
    m = mean_net(ref_dose, ref_d, **kw)
    c = int(round((G - 1) / 2))
    edge = (G - 1) / 2 + ref_d / 2
    return float(np.interp(edge, np.arange(G), m[c]))


if __name__ == "__main__":
    print(f"photons/nm^2 per mJ/cm^2: {PHOTONS_PER_NM2_PER_MJCM2:.3f}")
    thr = calibrate_threshold()
    print(f"threshold: {thr:.3f}")
    rng = np.random.default_rng(0)
    for d in (16, 20, 24):
        for dose in (10, 20, 30, 40, 60):
            cd, fail = simulate(dose, d, 2000, thr, rng=rng)
            print(f"d={d} dose={dose}: meanCD={cd.mean():5.1f} sdCD={cd.std():4.2f} fail={fail.mean():.4f}")
