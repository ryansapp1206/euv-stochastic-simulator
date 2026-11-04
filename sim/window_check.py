"""Failure-free dose window, mirroring ASML 2021 (simulated): P40 dense holes, 0.33 NA.
ASML: a window with print failure rate < 1e-6 needs ~79 mJ/cm2 at 0.33 NA.
For a resist sized to print 20 nm holes at D0, sweep dose around D0, count missing and
merged holes, fit each cliff's tail (rates <= 0.2) as log10(rate) linear in dose, and find
where each reaches 1e-6 and where the two fitted cliffs cross (the bottom of the valley).
Pass 2: sweep extended to 2.2 x D0 because merging now needs a real bridge.
Usage: python window_check.py [layers]. Writes results/window_check_L<layers>.json."""
import os
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
import json
import sys
import time
from dataclasses import replace

import numpy as np

from simulator_v2 import Params, simulate, dose_to_size

layers = int(sys.argv[1]) if len(sys.argv) > 1 else 1
out = {"layers": layers}
t0 = time.time()
for D0 in (30, 50, 80):
    p0 = dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=40, dose=D0, layers=layers),
                      n_trials=80 if layers > 1 else 120)
    rows = []
    for f in (0.55, 0.65, 0.75, 0.85, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2):
        n = 1500 if layers == 1 else 300
        r = simulate(replace(p0, dose=D0 * f), n, np.random.default_rng(int(D0 * f * 10)))
        ok = ~(r["missing"] | r["merged"])
        rows.append(dict(dose=D0 * f, holes=int(r["cd"].size), missing=int(r["missing"].sum()),
                         merged=int(r["merged"].sum()),
                         mean_cd=float(r["cd"][ok].mean()) if ok.any() else float("nan"),
                         sigma_cd=float(r["cd"][ok].std()) if ok.any() else float("nan")))
        print(D0, rows[-1], flush=True)

    def tail(key, side):
        pts = [(x["dose"], x[key] / x["holes"]) for x in rows if 0 < x[key] and x[key] <= 0.2 * x["holes"]]
        pts = [q for q in pts if (q[0] < D0 if side == "low" else q[0] > D0)]
        if len(pts) < 2:
            return None
        d, pr = np.array(pts).T
        return np.polyfit(d, np.log10(pr), 1)              # slope, intercept

    fm, fg = tail("missing", "low"), tail("merged", "high")
    lo = float((-6 - fm[1]) / fm[0]) if fm is not None else None   # missing < 1e-6 above this dose
    hi = float((-6 - fg[1]) / fg[0]) if fg is not None else None   # merged < 1e-6 below this dose
    cross = None
    if fm is not None and fg is not None:
        dc = (fg[1] - fm[1]) / (fm[0] - fg[0])
        cross = dict(dose=float(dc), rate=float(10 ** (fm[0] * dc + fm[1])))
    out[str(D0)] = dict(rho_q=p0.rho_q, sweep=rows, missing_below_1e6_above=lo, merged_below_1e6_below=hi,
                        window_exists=(lo is not None and hi is not None and lo < hi), cliffs_cross=cross)
    print("D0", D0, "missing<1e-6 above", lo, "merged<1e-6 below", hi, "cross", cross, flush=True)
out["_runtime_s"] = time.time() - t0
json.dump(out, open(RESULTS + f"/window_check_L{layers}.json", "w"), indent=1, default=float)
print("done in %.0f s" % out["_runtime_s"])
