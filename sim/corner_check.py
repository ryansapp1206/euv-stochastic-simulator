"""Final check: numerical convergence at the corners of the dataset domain.

Verification checked pixel size, tile size and source sampling at a few
calibration points. The dataset reaches further (32 nm pitch with 15 nm holes
and the sharpest blur, 64 nm pitch with the softest blur, high and low doses,
extreme development thresholds), so the same checks run here at those corners,
on the side cliffs as well as at nominal dose:

  px     1 nm vs 0.5 nm pixels (same resist; then refitted)
  tile   3x3 vs 5x5 holes per tile
  fit    resist fit with 120 vs 400 tiles per bisection step

Usage: python corner_check.py        ->  ../results/corner_check.json
       python corner_check.py sides  ->  ../results/corner_refit_sides.json (pixel size at matched hole size)
"""
import json
import os
import sys
import time
from dataclasses import replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from simulator_v2 import Params, simulate, summarize, dose_to_size  # noqa: E402
from make_dataset import to_params, FIXED  # noqa: E402

RESULTS = os.path.join(HERE, "..", "results")
BASE = dict(mask_bias_nm=1.0, qe=2.8, rho_pag=0.2, rho_pg=2.0, absorptivity_um=4.2, dev_threshold=0.5,
            thickness_nm=35.0)
CORNERS = {
    "P32_CD15_sharp_lowdose": dict(BASE, pitch_nm=32, cd_ratio=0.47, d0_mjcm2=20, sigma_e_nm=1.2, sigma_a_nm=4.0,
                                   thickness_nm=30.0),
    "P32_CD18.6_highdose_dev0.65": dict(BASE, pitch_nm=32, cd_ratio=0.58, d0_mjcm2=90, sigma_e_nm=2.5,
                                        sigma_a_nm=4.0, dev_threshold=0.65),
    "P64_CD26.9_soft_lowdose_dev0.35": dict(BASE, pitch_nm=64, cd_ratio=0.42, d0_mjcm2=20, sigma_e_nm=2.5,
                                            sigma_a_nm=6.0, dev_threshold=0.35, thickness_nm=50.0),
    "P40_CD20_mid": dict(BASE, pitch_nm=40, cd_ratio=0.50, d0_mjcm2=45, sigma_e_nm=1.7, sigma_a_nm=5.0),
}
FACTORS = (0.75, 1.0, 1.65)


def run(p, n, seed):
    s = summarize(simulate(p, n, np.random.default_rng(seed)))
    return {k: (round(v, 5) if isinstance(v, float) else v) for k, v in s.items()}


def refit_sides():
    """Pixel size at MATCHED hole size: refit the resist at 0.5 nm pixels, then compare the cliffs."""
    rows = {}
    for name in ("P32_CD15_sharp_lowdose", "P40_CD20_mid"):
        s = CORNERS[name]
        p1 = dose_to_size(to_params(s))
        p05 = dose_to_size(replace(to_params(s), px=0.5), n_trials=60)
        for f in (0.75, 0.85, 1.65):
            a = run(replace(p1, dose=p1.dose * f), 600, 11)
            b = run(replace(p05, dose=p05.dose * f), 300, 12)
            rows[f"{name}_f{f}"] = dict(px1=a, px05=b)
            print(name, f, "CD", a["mean_cd"], b["mean_cd"], "miss", a["p_missing"], b["p_missing"],
                  "merged", a["p_merged"], b["p_merged"], flush=True)
    json.dump(rows, open(os.path.join(RESULTS, "corner_refit_sides.json"), "w"), indent=1)


if __name__ == "__main__" and sys.argv[1:] == ["sides"]:
    refit_sides()
elif __name__ == "__main__":
    out = {}
    t0 = time.time()
    for name, s in CORNERS.items():
        p = dose_to_size(to_params(s))
        row = dict(rho_q=p.rho_q, cd_target=p.cd_target)
        for f in FACTORS:
            q = replace(p, dose=p.dose * f)
            row[f"px1_f{f}"] = run(q, 600, 1)
            row[f"px0.5_f{f}"] = run(replace(q, px=0.5), 300, 2)
            row[f"tile5_f{f}"] = run(replace(q, cells_x=5, cells_y=5), 220, 3)
        p05 = dose_to_size(replace(to_params(s), px=0.5), n_trials=60)
        row["px0.5_refit"] = dict(rho_q=p05.rho_q, **run(p05, 300, 4))
        p400 = dose_to_size(to_params(s), n_trials=400, rng_seed=7)
        row["fit400"] = dict(rho_q=p400.rho_q, **run(p400, 600, 1))
        out[name] = row
        print(name, json.dumps(row), flush=True)
    out["_runtime_s"] = time.time() - t0
    json.dump(out, open(os.path.join(RESULTS, "corner_check.json"), "w"), indent=1)

