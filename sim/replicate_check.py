"""Final check: are the per-row error bars honest?

Re-runs a random sample of dataset rows twice with fresh seeds and compares the
two runs: (a - b) / sqrt(se_a^2 + se_b^2) should have spread 1.0 if the stated
standard errors are right. Also compares failure counts with the binomial
expectation.

Usage: python replicate_check.py ../data/pilot [n_rows]  ->  <dataset>/replicate_check.json
"""
import json
import os
import sys
from dataclasses import replace

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from simulator_v2 import Params, simulate  # noqa: E402
from make_dataset import FIXED, label, seed_for  # noqa: E402


def params_from_row(r):
    return Params(cd_target=r.cd_target_nm, mask_cd=r.mask_cd_nm, pitch_x=r.pitch_nm, pitch_y=r.pitch_nm,
                  dose=r.dose_mjcm2, sigma_e=r.sigma_e_nm, sigma_a=r.sigma_a_nm, qe=r.qe, rho_pag=r.rho_pag,
                  rho_q=r.rho_q, rho_pg=r.rho_pg, thickness=r.thickness_nm, absorptivity=r.absorptivity_um / 1000,
                  dev_threshold=r.dev_threshold, **FIXED)


if __name__ == "__main__":
    path = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    meta = json.load(open(os.path.join(path, "metadata.json")))
    D = pd.read_csv(os.path.join(path, "rows.csv"))
    S = D.sample(n=min(n, len(D)), random_state=0)
    out = []
    for i, r in enumerate(S.itertuples()):
        p = params_from_row(r)
        # first: the stored row must regenerate exactly from its seed
        again = label(simulate(p, meta["trials"], np.random.default_rng(int(r.seed))))
        exact = np.isclose(again["mean_cd_nm"], r.mean_cd_nm, equal_nan=True) and again["n_missing"] == r.n_missing
        a = label(simulate(p, meta["trials"], np.random.default_rng(seed_for(424242, i, 1))))
        b = label(simulate(p, meta["trials"], np.random.default_rng(seed_for(424242, i, 2))))
        out.append(dict(resist_id=r.resist_id, dose_index=r.dose_index, regenerates=bool(exact),
                        **{f"{k}_a": a[k] for k in a}, **{f"{k}_b": b[k] for k in b}))
    R = pd.DataFrame(out)
    M = R[(R.n_ok_a >= 50) & (R.n_ok_b >= 50)]
    zm = (M.mean_cd_nm_a - M.mean_cd_nm_b) / np.hypot(M.mean_cd_se_a, M.mean_cd_se_b)
    zs = (M.sigma_cd_nm_a - M.sigma_cd_nm_b) / np.hypot(M.sigma_cd_se_a, M.sigma_cd_se_b)
    F = R[(R.n_missing_a + R.n_missing_b) > 0]
    pm = (F.n_missing_a + F.n_missing_b) / (F.n_holes_a + F.n_holes_b)
    zf = (F.n_missing_a - F.n_missing_b) / np.sqrt(2 * F.n_holes_a * pm * (1 - pm) + 1e-12)
    res = dict(rows=len(R), regenerate_exactly=int(R.regenerates.sum()), pairs_for_size=len(M),
               spread_z_mean_cd=float(zm.std()), spread_z_sigma=float(zs.std()),
               max_abs_z=float(np.abs(np.r_[zm, zs]).max()), pairs_with_missing=len(F),
               spread_z_missing=float(zf.std()) if len(F) > 1 else None)
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(path, "replicate_check.json"), "w"), indent=1)
