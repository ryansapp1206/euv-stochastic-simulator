"""Generate the ML dataset from simulator v2.

Design
------
* Each sample is one *resist + layout*: pitch, hole size, mask bias, resist
  parameters and a sizing dose D0, drawn by Latin hypercube over the domain
  where the simulator was verified (DOMAIN below, see results/DATASET.md).
* The resist's quencher loading is fitted so holes print at the target size at
  D0, exactly as in the calibration. Samples whose fitted quencher falls
  outside the published range (0.012-0.09 /nm^3), or whose fit does not
  reproduce the target size, are REJECTED and logged, not kept.
* Each accepted resist is then exposed at 8 doses spanning both failure
  cliffs: 0.60-1.90 x D0, one jittered dose per stratum plus D0 itself.
* Every row has its own independent random stream derived from
  (seed, resist_id, dose_index), so any row can be regenerated exactly and a
  run can be resumed.
* Rows carry counts and standard errors, so the models can weight by
  precision and the binomial failure counts are kept exact. Missing holes are
  independent per hole (binomial over n_holes). Merges are not: each bridge
  flags 2+ holes, so model n_merge_events out of n_neighbor_pairs instead.

Outputs (in --out):  rows.csv (one row per resist x dose), resists.csv (one row
per sampled resist, accepted or not, with QA numbers), replicates.csv (repeat
runs of some rows with other seeds, for checking the error bars),
metadata.json (simulator version, domain, settings).

Usage:
  python make_dataset.py --out ../data/pilot --n-resists 60 --workers 2
  python make_dataset.py --out ../data/full --n-resists 1500 --workers 8   # resumable
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from dataclasses import replace, asdict

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import scipy
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from simulator_v2 import (Params, PHOTONS_PER_NM2_PER_MJCM2, simulate, summarize, dose_to_size,  # noqa: E402
                          aerial_image)

# ---------------------------------------------------------------- domain
# (low, high, "lin"|"log"). Every range is one the simulator was checked over
# (VERIFICATION.md sensitivity study and comparisons).
DOMAIN = {
    "pitch_nm":        (32.0, 64.0, "lin"),   # square dense arrays; P30 was a poor match, so stop at 32
    "cd_ratio":        (0.42, 0.58, "lin"),   # target hole size / pitch
    "mask_bias_nm":    (0.0, 4.0, "lin"),     # mask hole minus target (Goethals used +5, PROLITH +1)
    "d0_mjcm2":        (20.0, 90.0, "log"),   # sizing dose
    "sigma_e_nm":      (1.2, 2.5, "lin"),     # electron blur
    "sigma_a_nm":      (4.0, 6.0, "lin"),     # acid diffusion
    "qe":              (2.0, 3.5, "lin"),     # acids per absorbed photon
    "rho_pag":         (0.12, 0.30, "lin"),   # PAG /nm^3
    "rho_pg":          (1.5, 3.0, "lin"),     # protecting groups /nm^3
    "thickness_nm":    (30.0, 50.0, "lin"),
    "absorptivity_um": (4.2, 5.3, "lin"),     # /um (CAR)
    "dev_threshold":   (0.35, 0.65, "lin"),   # deprotected fraction that dissolves
}
MIN_CD_NM = 15.0                    # smallest hole size checked in verification
RHO_Q_RANGE = (0.012, 0.090)        # published quencher loadings
FIT_TOL_NM = 0.5                    # |printed - target| allowed at D0
DOSE_STRATA = (0.60, 0.75, 0.85, 0.95, 1.05, 1.20, 1.40, 1.65, 1.90)
FIXED = dict(na=0.33, illumination="quasar", sigma_q=1.3, k_amp=60.0, cells_x=3, cells_y=3, px=1.0, layers=1)

ROW_FIELDS = ["resist_id", "dose_index", "seed", "pitch_nm", "cd_target_nm", "cd_ratio", "mask_cd_nm",
              "mask_bias_nm", "dose_mjcm2", "d0_mjcm2", "dose_factor", "sigma_e_nm", "sigma_a_nm", "sigma_q_nm",
              "qe", "rho_pag", "rho_q", "rho_pg", "k_amp", "thickness_nm", "absorptivity_um", "dev_threshold",
              "na", "absorbed_frac", "photons_abs_per_nm2", "image_peak", "nils",
              "n_holes", "n_ok", "n_missing", "n_merged", "n_merge_events", "n_neighbor_pairs",
              "p_missing", "p_merged",
              "mean_cd_nm", "mean_cd_se", "sigma_cd_nm", "sigma_cd_se", "lcdu_3s_nm"]
RESIST_FIELDS = ["resist_id", "accepted", "reason", "pitch_nm", "cd_target_nm", "cd_ratio", "mask_bias_nm",
                 "d0_mjcm2", "sigma_e_nm", "sigma_a_nm", "qe", "rho_pag", "rho_pg", "thickness_nm",
                 "absorptivity_um", "dev_threshold", "rho_q", "fit_cd_nm", "fit_error_nm",
                 "analytic_sigma_nm", "mc_sigma_nm", "mc_over_analytic", "seconds"]
REPLICATE_FIELDS = ["resist_id", "dose_index", "seed", "n_ok", "mean_cd_nm", "mean_cd_se", "sigma_cd_nm",
                    "sigma_cd_se", "n_missing", "n_merged", "n_merge_events", "n_holes"]


# ---------------------------------------------------------------- sampling
def sample_domain(n: int, seed: int) -> list[dict]:
    """n Latin hypercube samples over DOMAIN (before the validity filter)."""
    keys = list(DOMAIN)
    u = qmc.LatinHypercube(d=len(keys), seed=seed).random(n)
    out = []
    for i, row in enumerate(u):
        s = {}
        for k, x in zip(keys, row):
            lo, hi, scale = DOMAIN[k]
            s[k] = float(np.exp(np.log(lo) + x * (np.log(hi) - np.log(lo)))) if scale == "log" \
                else float(lo + x * (hi - lo))
        s["pitch_nm"] = float(round(s["pitch_nm"]))           # tile needs a whole number of pixels
        s["resist_id"] = i
        out.append(s)
    return out


def to_params(s: dict) -> Params:
    cd = s["cd_ratio"] * s["pitch_nm"]
    return Params(cd_target=cd, mask_cd=cd + s["mask_bias_nm"], pitch_x=s["pitch_nm"], pitch_y=s["pitch_nm"],
                  dose=s["d0_mjcm2"], sigma_e=s["sigma_e_nm"], sigma_a=s["sigma_a_nm"], qe=s["qe"],
                  rho_pag=s["rho_pag"], rho_pg=s["rho_pg"], thickness=s["thickness_nm"],
                  absorptivity=s["absorptivity_um"] / 1000.0, dev_threshold=s["dev_threshold"], **FIXED)


def dose_factors(rng: np.random.Generator) -> list[float]:
    """One jittered factor per dose stratum; the stratum holding 1.0 uses exactly 1.0."""
    f = []
    for lo, hi in zip(DOSE_STRATA[:-1], DOSE_STRATA[1:]):
        f.append(1.0 if lo < 1.0 < hi else float(rng.uniform(lo, hi)))
    return f


def seed_for(base: int, *keys: int) -> int:
    """63-bit seed unique to (base, keys): no collisions in practice, and fits a CSV int64 column."""
    return int(np.random.SeedSequence([base, *keys]).generate_state(1, dtype=np.uint64)[0] & ((1 << 63) - 1))


# ---------------------------------------------------------------- physics features
def image_features(p: Params) -> tuple[float, float]:
    """Peak intensity and NILS (target CD x image log-slope) at the target edge, along x through a hole centre."""
    from scipy.ndimage import map_coordinates
    I = aerial_image(p)
    xc, yc = p.pitch_x / 2, p.pitch_y / 2                     # first hole centre, nm
    xe = xc + p.cd_target / 2
    xs = np.array([xe - 0.25, xe + 0.25])
    vals = map_coordinates(I, [np.full(2, yc / p.px - 0.5), xs / p.px - 0.5], order=1, mode="grid-wrap")
    slope = abs(np.log(vals[1]) - np.log(vals[0])) / 0.5
    return float(I.max()), float(p.cd_target * slope)


def label(res: dict) -> dict:
    s = summarize(res)
    ok = ~(res["missing"] | res["merged"])
    n, n_ok = res["cd"].size, int(ok.sum())
    sd = s["sigma_cd"]
    g = res["merge_group"]
    n_events = int(np.unique(g[g >= 0]).size)            # one bridge joining 2+ holes = one event
    return dict(n_holes=n, n_ok=n_ok, n_missing=int(res["missing"].sum()), n_merged=int(res["merged"].sum()),
                n_merge_events=n_events, n_neighbor_pairs=2 * n,
                p_missing=s["p_missing"], p_merged=s["p_merged"], mean_cd_nm=s["mean_cd"],
                mean_cd_se=float(sd / np.sqrt(n_ok)) if n_ok > 1 else float("nan"),
                sigma_cd_nm=sd, sigma_cd_se=float(sd / np.sqrt(2 * (n_ok - 1))) if n_ok > 1 else float("nan"),
                lcdu_3s_nm=s["lcdu_3s"])


# ---------------------------------------------------------------- one resist
def process(args):
    """One sample -> (resist record, rows, replicate rows). Errors are logged as rejections, never lost."""
    s, cfg = args
    try:
        return _process(s, cfg)
    except Exception as e:                                   # keep a long run alive; the checker reports these
        rec = {k: s.get(k) for k in RESIST_FIELDS if k in s}
        rec.update(resist_id=s["resist_id"], accepted=False, reason=f"error: {type(e).__name__}: {e}"[:200])
        return rec, [], []


def _process(s, cfg):
    t0 = time.time()
    rid = s["resist_id"]
    rec = {k: s.get(k) for k in RESIST_FIELDS if k in s}
    rec.update(resist_id=rid, cd_target_nm=s["cd_ratio"] * s["pitch_nm"])
    rows, reps = [], []
    if rec["cd_target_nm"] < MIN_CD_NM:
        rec.update(accepted=False, reason=f"cd<{MIN_CD_NM}", seconds=time.time() - t0)
        return rec, rows, reps
    p0 = to_params(s)
    p = dose_to_size(p0, knob="rho_q", n_trials=cfg["fit_trials"], rng_seed=seed_for(cfg["seed"], rid, 999))
    rec["rho_q"] = p.rho_q
    if not (RHO_Q_RANGE[0] <= p.rho_q <= RHO_Q_RANGE[1]):
        rec.update(accepted=False, reason="rho_q_out_of_range", seconds=time.time() - t0)
        return rec, rows, reps
    rng = np.random.default_rng(seed_for(cfg["seed"], rid, 998))
    peak_nom, nils_nom = image_features(p)
    for j, f in enumerate(dose_factors(rng)):
        q = replace(p, dose=p.dose * f)
        sd = seed_for(cfg["seed"], rid, j)
        lab = label(simulate(q, cfg["trials"], np.random.default_rng(sd)))
        pk, nl = (peak_nom, nils_nom)          # the image does not depend on dose
        rows.append(dict(resist_id=rid, dose_index=j, seed=sd, pitch_nm=p.pitch_x, cd_target_nm=p.cd_target,
                         cd_ratio=s["cd_ratio"], mask_cd_nm=p.mask_cd, mask_bias_nm=s["mask_bias_nm"],
                         dose_mjcm2=q.dose, d0_mjcm2=p.dose, dose_factor=f, sigma_e_nm=p.sigma_e,
                         sigma_a_nm=p.sigma_a, sigma_q_nm=p.sigma_q, qe=p.qe, rho_pag=p.rho_pag, rho_q=p.rho_q,
                         rho_pg=p.rho_pg, k_amp=p.k_amp, thickness_nm=p.thickness,
                         absorptivity_um=p.absorptivity * 1000, dev_threshold=p.dev_threshold, na=p.na,
                         absorbed_frac=p.absorbed_frac,
                         photons_abs_per_nm2=q.dose * PHOTONS_PER_NM2_PER_MJCM2 * p.absorbed_frac,
                         image_peak=pk, nils=nl, **lab))
        if f == 1.0:
            rec["fit_cd_nm"] = lab["mean_cd_nm"]
            rec["fit_error_nm"] = lab["mean_cd_nm"] - p.cd_target
            rec["mc_sigma_nm"] = lab["sigma_cd_nm"]
        if cfg["replicate_every"] and rid % cfg["replicate_every"] == 0 and j in (0, 3, 7):
            sd2 = seed_for(cfg["seed"] + 1, rid, j)
            l2 = label(simulate(q, cfg["trials"], np.random.default_rng(sd2)))
            reps.append(dict(resist_id=rid, dose_index=j, seed=sd2, **{k: l2[k] for k in REPLICATE_FIELDS
                                                                      if k in l2}))
    if not np.isfinite(rec.get("fit_error_nm", np.nan)) or abs(rec["fit_error_nm"]) > FIT_TOL_NM:
        rec.update(accepted=False, reason="fit_check_failed", seconds=time.time() - t0)
        return rec, [], []
    if cfg["analytic"]:
        from analytic_check import analytic
        a = analytic(p)
        if a is not None:
            rec["analytic_sigma_nm"] = a["sigma_cd"]["all"]
            rec["mc_over_analytic"] = rec["mc_sigma_nm"] / a["sigma_cd"]["all"]
    rec.update(accepted=True, reason="", seconds=time.time() - t0)
    return rec, rows, reps


# ---------------------------------------------------------------- bookkeeping
def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=10).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def file_sha256(path: str) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def append(path: str, fields: list[str], rows: list[dict]):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        if new:
            w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-resists", type=int, default=60, help="Latin hypercube samples (before the filter)")
    ap.add_argument("--trials", type=int, default=200, help="tiles per row (x 9 holes)")
    ap.add_argument("--fit-trials", type=int, default=120, help="tiles per bisection step of the resist fit")
    ap.add_argument("--seed", type=int, default=20261116)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 0))
    ap.add_argument("--replicate-every", type=int, default=10, help="repeat 3 rows of every k-th resist (0=off)")
    ap.add_argument("--no-analytic", action="store_true", help="skip the semi-analytical cross-check")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    cfg = dict(trials=a.trials, fit_trials=a.fit_trials, seed=a.seed, replicate_every=a.replicate_every,
               analytic=not a.no_analytic)
    meta_path = os.path.join(a.out, "metadata.json")
    sim_path = os.path.join(HERE, "simulator_v2.py")
    meta = dict(simulator="simulator_v2.py", simulator_sha256=file_sha256(sim_path), git_commit=git_commit(),
                domain=DOMAIN, min_cd_nm=MIN_CD_NM, rho_q_range=RHO_Q_RANGE, fit_tol_nm=FIT_TOL_NM,
                dose_strata=DOSE_STRATA, fixed=FIXED, n_resists=a.n_resists, **cfg,
                python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
    if os.path.exists(meta_path):                     # resuming: settings must match
        old = json.load(open(meta_path))
        for k in ("simulator_sha256", "seed", "trials", "fit_trials", "n_resists", "domain"):
            if json.loads(json.dumps(old.get(k))) != json.loads(json.dumps(meta[k])):
                sys.exit(f"refusing to resume: {k} differs from the existing run in {a.out}")
    else:
        meta["started"] = time.strftime("%Y-%m-%d %H:%M:%S")
        json.dump(meta, open(meta_path, "w"), indent=1)

    done = set()
    rpath = os.path.join(a.out, "resists.csv")
    if os.path.exists(rpath):
        with open(rpath) as fh:
            done = {int(r["resist_id"]) for r in csv.DictReader(fh)}
        # drop rows of a resist that was interrupted after its rows were written
        for name, fields in (("rows.csv", ROW_FIELDS), ("replicates.csv", REPLICATE_FIELDS)):
            path = os.path.join(a.out, name)
            if os.path.exists(path):
                with open(path) as fh:
                    keep = [r for r in csv.DictReader(fh) if int(r["resist_id"]) in done]
                os.replace(path, path + ".bak")
                append(path, fields, keep)
                os.remove(path + ".bak")
    todo = [s for s in sample_domain(a.n_resists, a.seed) if s["resist_id"] not in done]
    print(f"{len(done)} resists done, {len(todo)} to go, {a.workers} workers", flush=True)
    t0 = time.time()
    jobs = [(s, cfg) for s in todo]
    if a.workers > 1:
        import multiprocessing as mp
        pool = mp.get_context("spawn").Pool(a.workers)
        it = pool.imap_unordered(process, jobs)
    else:
        pool, it = None, map(process, jobs)
    n_acc = 0
    for k, (rec, rows, reps) in enumerate(it, 1):
        # rows first, then the resist record: a resist counts as done only once its rows are on disk
        append(os.path.join(a.out, "rows.csv"), ROW_FIELDS, rows)
        append(os.path.join(a.out, "replicates.csv"), REPLICATE_FIELDS, reps)
        append(rpath, RESIST_FIELDS, [rec])
        n_acc += bool(rec["accepted"])
        if k % 10 == 0 or k == len(jobs):
            el = time.time() - t0
            print(f"{k}/{len(jobs)} resists ({n_acc} accepted), {el / 60:.1f} min, "
                  f"~{el / k * (len(jobs) - k) / 60:.0f} min left", flush=True)
    if pool:
        pool.close(); pool.join()
    meta = json.load(open(meta_path))
    meta["finished"] = time.strftime("%Y-%m-%d %H:%M:%S")
    json.dump(meta, open(meta_path, "w"), indent=1)


if __name__ == "__main__":
    main()
