"""Pass 2 cross-check: a semi-analytical noise model versus the Monte Carlo.

The Monte Carlo draws random molecules and thresholds the result. This script
predicts the same hole-size variation WITHOUT any random numbers, using
first-order (linear) noise propagation, the method behind Gaussian random
field stochastic models (e.g. De Bisschop & Hansen; Gallatin/Naulleau).

For one hole, the printed area changes by
    dA = sum_j  r_j dtheta * dh(x_j) / (-dh/dr)(x_j)
along the mean resist contour (h = net acid, threshold h_th = ln 2 / K).
dA is a linear functional  sum_g U_g dh_g  of the net-acid noise field, so for
each independent Poisson source with mean lambda (per voxel) blurred by a
Gaussian G, its variance is  sum_y lambda_y (G * U)_y^2 / vox^2  (Campbell's
theorem).  Photon noise passes through electron blur, acid conversion
(gain QE * exp(-QE * E / PAG0)) and acid diffusion.  Protecting-group noise is
white per pixel: each pixel near the edge clears with probability
P(dep > rem) (Skellam), adding sum P(1 - P) to the area variance.
sigma_CD = sigma_A / (pi * CD / 2).

If the Monte Carlo is coded correctly, both should agree in the nominal
regime (where linear response holds) for every noise source separately.

Usage: python analytic_check.py             ->  ../results/analytic_check.json
       python analytic_check.py linearity   ->  ../results/analytic_linearity.json
"""
import os
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
import json
import sys
import time
from dataclasses import replace

import numpy as np
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.stats import norm, skellam

from simulator_v2 import (Params, PHOTONS_PER_NM2_PER_MJCM2, aerial_image, _hole_centres,
                          simulate, summarize, dose_to_size)


def analytic(p: Params, n_theta: int = 180):
    if p.layers != 1:
        raise ValueError("analytic model is for layers=1")
    I = aerial_image(p).astype(float)
    H, W = I.shape
    px = p.px
    vox = p.thickness * px * px
    Phi = p.dose * PHOTONS_PER_NM2_PER_MJCM2 * px * px * p.absorbed_frac * I   # photons / column
    pag0, q0, pg0 = p.rho_pag * vox, p.rho_q * vox, p.rho_pg * vox

    def g(x, s):
        return gaussian_filter(x, s / px, mode="wrap")

    E = g(Phi, p.sigma_e)
    c = np.exp(-p.qe * E / pag0)                       # unconverted PAG fraction
    A = pag0 * (1 - c)                                 # mean acids / voxel
    hbar = (g(A, p.sigma_a) - q0) / vox                # mean net acid / nm^3
    h_th = -np.log(1 - p.dev_threshold) / p.k_amp       # frac = 1 - exp(-K h) crosses the threshold
    gy, gx = np.gradient(hbar, px)                     # nm^-4 (interior; holes sit away from edges)

    cx, cy = _hole_centres(p)
    th = np.arange(n_theta) * 2 * np.pi / n_theta
    dth = 2 * np.pi / n_theta
    rgrid = np.arange(0, min(p.pitch_x, p.pitch_y) / 2, 0.01)
    out = []
    for x0, y0 in zip(cx, cy):
        # pixel-index coordinates: pixel centre (i + 0.5) * px
        ri = []
        for t in th:
            xs = (x0 + rgrid * np.cos(t)) / px - 0.5
            ys = (y0 + rgrid * np.sin(t)) / px - 0.5
            hv = map_coordinates(hbar, [ys, xs], order=1, mode="grid-wrap")
            below = np.nonzero(hv < h_th)[0]
            if len(below) == 0 or below[0] == 0:
                ri.append(np.nan); continue
            k = below[0]
            f = (hv[k - 1] - h_th) / (hv[k - 1] - hv[k])
            ri.append(rgrid[k - 1] + f * (rgrid[k] - rgrid[k - 1]))
        ri = np.array(ri)
        if np.isnan(ri).any():
            return None
        xe, ye = x0 + ri * np.cos(th), y0 + ri * np.sin(th)
        ie, je = ye / px - 0.5, xe / px - 0.5
        dhdx = map_coordinates(gx, [ie, je], order=1, mode="grid-wrap")
        dhdy = map_coordinates(gy, [ie, je], order=1, mode="grid-wrap")
        dhdr = dhdx * np.cos(th) + dhdy * np.sin(th)
        w = ri * dth / (-dhdr)                          # nm^2 of area per unit dh
        # spread the weights onto the grid (bilinear), giving the functional U
        U = np.zeros((H, W))
        i0, j0 = np.floor(ie).astype(int), np.floor(je).astype(int)
        fi, fj = ie - i0, je - j0
        for di, dj, ww in ((0, 0, (1 - fi) * (1 - fj)), (1, 0, fi * (1 - fj)),
                           (0, 1, (1 - fi) * fj), (1, 1, fi * fj)):
            np.add.at(U, ((i0 + di) % H, (j0 + dj) % W), w * ww)
        var = dict(
            photon=float(np.sum(Phi * g(p.qe * c * g(U, p.sigma_a), p.sigma_e) ** 2) / vox ** 2),
            pag=float(np.sum(A * g(U, p.sigma_a) ** 2) / vox ** 2),
            quencher=float(np.sum(q0 * g(U, p.sigma_q) ** 2) / vox ** 2),
        )
        area = 0.5 * np.sum(ri ** 2) * dth
        out.append((area, var))

    # protecting-group (white) noise and the PG-smeared mean area, per hole cell
    m = np.exp(-p.k_amp * np.maximum(hbar, 0))
    mu_dep, mu_rem = pg0 * (1 - m), pg0 * m
    if abs(p.dev_threshold - 0.5) < 1e-12:
        P = skellam.sf(0, np.maximum(mu_dep, 1e-12), np.maximum(mu_rem, 1e-12))  # P(dep - rem >= 1)
    else:   # normal approximation to P(dep / (dep + rem) > threshold)
        f = 1 - m
        P = norm.cdf((f - p.dev_threshold) / np.sqrt(np.maximum(f * (1 - f), 1e-12) / pg0))
    from simulator_v2 import _geometry
    _, cell, nh = _geometry(p)
    var_pg = np.array([np.sum((P * (1 - P))[cell == k]) for k in range(nh)]) * px ** 4
    area_pg = np.array([np.sum(P[cell == k]) for k in range(nh)]) * px * px

    res = {}
    areas = np.array([a for a, _ in out])
    cd = 2 * np.sqrt(areas / np.pi)
    res["mean_cd_contour"] = float(cd.mean())
    res["mean_cd_pixels_pg"] = float((2 * np.sqrt(area_pg / np.pi)).mean())
    sig = {}
    for src in ("photon", "pag", "quencher"):
        v = np.array([vv[src] for _, vv in out])
        sig[src] = float(np.sqrt(np.mean(v / (np.pi * cd / 2) ** 2)))
    sig["pg"] = float(np.sqrt(np.mean(var_pg / (np.pi * cd / 2) ** 2)))
    sig["all"] = float(np.sqrt(sum(s ** 2 for s in sig.values())))
    res["sigma_cd"] = sig
    return res


def mc(p, noise, n, seed):
    s = summarize(simulate(replace(p, noise=noise), n, np.random.default_rng(seed)))
    return dict(sigma_cd=s["sigma_cd"], mean_cd=s["mean_cd"], n_holes=s["n_holes"],
                p_missing=s["p_missing"], p_merged=s["p_merged"])


def cases():
    return {
        "C1_generic_P40_CD20_30mJ": dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=40, dose=30)),
        "C2_severi_H80V46_80mJ": dose_to_size(Params(cd_target=24, pitch_x=80, pitch_y=46, stagger=True,
                                                     cells_x=2, cells_y=4, illumination="annular", dose=80)),
        "C3_lowdose_P40_CD20_21mJ": dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=40, dose=21.1)),
        "C4_goethals_NA025_60nm": dose_to_size(Params(cd_target=26, mask_cd=31, pitch_x=52, pitch_y=52,
                                                      dose=31.7, na=0.25, thickness=60)),
        "C5_naulleau_P32": dose_to_size(Params(cd_target=17.17, pitch_x=32, pitch_y=32, dose=82, qe=3.0,
                                               rho_q=0.085, sigma_a=5.1, cells_x=4, cells_y=4), knob="k_amp"),
    }

def linearity():
    """Quencher noise has a short range (1.3 nm), where linear theory over-predicts. If that is
    the reason for the gap, the gap must vanish for longer quencher blur and not depend on pixel size."""
    p = dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=40, dose=30))
    rows = {}
    for lab, q in [("sigma_q_1.3", p), ("sigma_q_2.5", replace(p, sigma_q=2.5)),
                   ("sigma_q_4.0", replace(p, sigma_q=4.0)), ("sigma_q_1.3_px_0.5", replace(p, px=0.5))]:
        a = analytic(q)["sigma_cd"]["quencher"]
        m = mc(q, ("quencher",), 300 if q.px == 1 else 120, 11)["sigma_cd"]
        rows[lab] = dict(analytic=a, monte_carlo=m, ratio=m / a)
        print(lab, rows[lab], flush=True)
    json.dump(rows, open(RESULTS + "/analytic_linearity.json", "w"), indent=1)


if __name__ == "__main__" and sys.argv[1:] == ["linearity"]:
    linearity()
elif __name__ == "__main__":
    t0 = time.time()
    results = {}
    for name, p in cases().items():
        a = analytic(p)
        row = dict(analytic=a, monte_carlo={})
        n = 400 if p.pitch_x * p.pitch_y <= 2500 else 250
        for i, src in enumerate(("photon", "pag", "quencher", "pg")):
            row["monte_carlo"][src] = mc(p, (src,), n, 100 + i)
        row["monte_carlo"]["all"] = mc(p, ("photon", "pag", "quencher", "pg"), n, 200)
        results[name] = row
        print(name, json.dumps({"analytic": a["sigma_cd"], "mc": {k: round(v["sigma_cd"], 4)
                                for k, v in row["monte_carlo"].items()},
                                "cd": (round(a["mean_cd_contour"], 2), round(a["mean_cd_pixels_pg"], 2),
                                       round(row["monte_carlo"]["all"]["mean_cd"], 2))}), flush=True)
    results["_runtime_s"] = time.time() - t0
    json.dump(results, open(RESULTS + "/analytic_check.json", "w"), indent=1)
