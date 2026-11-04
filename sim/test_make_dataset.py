"""Tests for the dataset generator. Run: python -m pytest -q test_make_dataset.py (about a minute)."""
import csv
import os

import numpy as np

from make_dataset import DOMAIN, DOSE_STRATA, sample_domain, dose_factors, seed_for, to_params, main, process


def test_samples_stay_in_domain_and_are_latin_hypercube():
    S = sample_domain(50, seed=1)
    assert S == sample_domain(50, seed=1)                      # reproducible
    for k, (lo, hi, scale) in DOMAIN.items():
        v = np.array([s[k] for s in S])
        if k == "pitch_nm":
            lo, hi = round(lo), round(hi)
        assert v.min() >= lo - 1e-9 and v.max() <= hi + 1e-9, k
        if k != "pitch_nm":                                    # one sample per stratum
            u = (np.log(v / lo) / np.log(hi / lo)) if scale == "log" else (v - lo) / (hi - lo)
            assert sorted(np.floor(u * 50).astype(int)) == list(range(50)), k


def test_dose_ladder_covers_each_stratum_once():
    f = dose_factors(np.random.default_rng(0))
    assert len(f) == len(DOSE_STRATA) - 1 and f.count(1.0) == 1
    for x, lo, hi in zip(f, DOSE_STRATA[:-1], DOSE_STRATA[1:]):
        assert lo <= x <= hi


def test_seeds_are_unique():
    seeds = {seed_for(7, r, j) for r in range(300) for j in range(10)}
    assert len(seeds) == 3000


def test_params_mapping():
    s = sample_domain(1, seed=3)[0]
    p = to_params(s)
    assert abs(p.cd_target - s["cd_ratio"] * s["pitch_nm"]) < 1e-12
    assert abs(p.mask_cd - p.cd_target - s["mask_bias_nm"]) < 1e-12
    assert abs(p.absorptivity * 1000 - s["absorptivity_um"]) < 1e-12 and p.pitch_x == p.pitch_y == s["pitch_nm"]


CFG = dict(trials=12, fit_trials=12, seed=5, replicate_every=0, analytic=False)


def test_rows_are_reproducible():
    s = next(x for x in sample_domain(8, seed=2) if x["cd_ratio"] * x["pitch_nm"] >= 15)
    a, b = process((s, CFG)), process((s, CFG))
    assert a[1] == b[1] and len(a[1]) in (0, len(DOSE_STRATA) - 1)


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def test_resume_after_interruption(tmp_path):
    args = ["--n-resists", "3", "--trials", "12", "--fit-trials", "12", "--workers", "1",
            "--replicate-every", "0", "--no-analytic", "--seed", "11"]
    full = tmp_path / "full"
    main(["--out", str(full)] + args)
    # simulate a crash after the last resist's rows were written but before its record
    cut = tmp_path / "cut"
    os.makedirs(cut)
    for f in ("metadata.json", "rows.csv", "replicates.csv"):
        if (full / f).exists():
            (cut / f).write_bytes((full / f).read_bytes())
    recs = read(full / "resists.csv")
    with open(cut / "resists.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(recs[0]))
        w.writeheader()
        w.writerows(recs[:-1])
    main(["--out", str(cut)] + args)
    key = lambda r: (int(r["resist_id"]), int(r["dose_index"]))
    assert sorted(read(full / "rows.csv"), key=key) == sorted(read(cut / "rows.csv"), key=key)
