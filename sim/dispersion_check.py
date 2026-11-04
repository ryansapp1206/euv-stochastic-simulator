"""Final check: are failure counts binomial?

Models of failure rate will be trained with a binomial loss, which assumes the
holes in a row fail independently. Here one setting on each cliff is run 40
times with independent seeds; the variance of the failure count across runs is
compared with the binomial variance n p (1 - p). Ratio ~1 = binomial.
Merged holes come in pairs (one bridge flags two holes), so for merging the
check is repeated on merge events out of neighbour pairs.

Usage: python dispersion_check.py  ->  ../results/dispersion_check.json
"""
import json
import os
import sys
from dataclasses import replace

import numpy as np
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from simulator_v2 import Params, simulate, dose_to_size  # noqa: E402

RUNS, TRIALS = 40, 200
base = dose_to_size(Params(cd_target=18, pitch_x=36, pitch_y=36, dose=40))
cases = {"missing_side_f0.80": (replace(base, dose=32.0), "missing"),
         "merged_side_f1.60": (replace(base, dose=64.0), "merged"),
         "missing_side_f0.75_P48": (replace(dose_to_size(Params(cd_target=24, pitch_x=48, pitch_y=48, dose=40)),
                                            dose=30.0), "missing")}
out = {}
for name, (p, key) in cases.items():
    counts, per_tile, events = [], [], []
    for k in range(RUNS):
        r = simulate(p, TRIALS, np.random.default_rng(1000 + k))
        counts.append(int(r[key].sum()))
        g = r["merge_group"]
        events.append(int(np.unique(g[g >= 0]).size))
        per_tile.append(r[key].reshape(TRIALS, -1).sum(1))
    counts = np.array(counts)
    n = TRIALS * 9
    ph = counts.mean() / n
    ratio = counts.var(ddof=1) / (n * ph * (1 - ph))
    stat = (RUNS - 1) * ratio
    pval = float(2 * min(chi2.cdf(stat, RUNS - 1), chi2.sf(stat, RUNS - 1)))
    # within-tile clustering: variance of failures per tile vs binomial with 9 holes
    t = np.concatenate(per_tile)
    tile_ratio = t.var(ddof=1) / (9 * ph * (1 - ph))
    out[name] = dict(p=ph, mean_count=float(counts.mean()), variance_ratio=float(ratio), p_value=pval,
                     per_tile_variance_ratio=float(tile_ratio))
    if key == "merged":   # the independent unit: merge events out of neighbour pairs (2 per hole)
        ev = np.array(events)
        pe = ev.mean() / (2 * n)
        out[name]["events_variance_ratio"] = float(ev.var(ddof=1) / (2 * n * pe * (1 - pe)))
        out[name]["mean_events"] = float(ev.mean())
    print(name, json.dumps(out[name]), flush=True)
json.dump(out, open(os.path.join(HERE, "..", "results", "dispersion_check.json"), "w"), indent=1)
