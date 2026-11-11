"""New like-for-like comparisons against fully specified real CAR measurements.

Case A  Shin-Etsu patents: NXE:3300/3400, NA 0.33, quadrupole sigma 0.9/0.6, 23 nm holes,
        pitch 46 nm, mask +20% bias, 60 nm CAR film, PEB 60 s, 2.38% TMAH.
        Measured CDU (3s of 50 holes, one dose shot, CD-SEM): 3.0-3.4 nm at Eop 30-36 mJ/cm2.
Case B  imec/ASML NXE:3400B, regular (square) contact holes, 18 nm HP (P36), 33 mJ/cm2 typical
        for imec's P36 CAR; measured local CDU 4.6 nm (default) -> 3.6 nm (improved illumination).
Case C  imec/ASML NXE:3400, P40 (20 nm HP) dense contact holes, LCDU 3.0 nm (2.6 with pupil tuning).
Nothing in the simulator is modified: the real quadrupole is added by wrapping _source_points.
"""
import os
import sys
from dataclasses import replace
from multiprocessing import Pool

import numpy as np

import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import simulator_v2 as S  # noqa: E402

_orig = S._source_points


def _quad_points(p):
    """Quadrupole with poles at 45/135/225/315 deg, sigma in [in, out], opening angle in degrees.
    Encoded as illumination = 'quad:<in>:<out>:<opening>'."""
    _, s_in, s_out, opening = p.illumination.split(":")
    s_in, s_out, opening = float(s_in), float(s_out), float(opening)
    pts, w = [], []
    for centre in (45, 135, 225, 315):
        for r in np.linspace(s_in, s_out, 7):
            for a in np.linspace(centre - opening / 2, centre + opening / 2, 9):
                pts.append((r * np.cos(np.radians(a)), r * np.sin(np.radians(a))))
                w.append(r)                                # area weight ~ r
    w = np.array(w)
    return np.array(pts), w / w.sum()


def _patched(p):
    if p.illumination.startswith("quad:"):
        return _quad_points(p)
    return _orig(p)


S._source_points = _patched


def run(job):
    name, kw, fit_trials, n_trials = job
    p = S.Params(**kw)
    p = S.dose_to_size(p, "rho_q", n_trials=fit_trials, rng_seed=3)
    r = S.simulate(p, n_trials, np.random.default_rng(21))
    s = S.summarize(r)
    return name, kw, p.rho_q, s


if __name__ == "__main__":
    jobs = []
    # ---- Case A
    for a_um in (4.2, 5.3):
        for d0 in (30.0, 36.0):
            for ill in ("quad:0.6:0.9:30", "quad:0.6:0.9:90"):
                kw = dict(pitch_x=46, pitch_y=46, cd_target=23.0, mask_cd=23.0 * 1.2, dose=d0, thickness=60.0,
                          absorptivity=a_um / 1000.0, illumination=ill)
                jobs.append((f"A  abs={a_um} D0={d0:.0f} {ill.split(':',1)[1]}", kw, 100, 1500))
    # model's own idealised quasar, same resist, for reference
    for d0 in (30.0, 36.0):
        kw = dict(pitch_x=46, pitch_y=46, cd_target=23.0, mask_cd=23.0 * 1.2, dose=d0, thickness=60.0,
                  absorptivity=0.0042, illumination="quasar")
        jobs.append((f"A  abs=4.2 D0={d0:.0f} model-quasar(r=0.1)", kw, 100, 1500))
    # ---- Case B  (P36, 18 nm)
    for d0 in (33.0, 45.0, 60.0):
        for ill in ("quasar", "annular", "quad:0.6:0.9:30"):
            kw = dict(pitch_x=36, pitch_y=36, cd_target=18.0, dose=d0, illumination=ill)
            jobs.append((f"B  P36/18nm D0={d0:.0f} {ill}", kw, 100, 1500))
    # ---- Case C  (P40, 20 nm)
    for d0 in (40.0, 50.0, 60.0, 70.0):
        for ill in ("quasar", "quad:0.6:0.9:30"):
            kw = dict(pitch_x=40, pitch_y=40, cd_target=20.0, dose=d0, illumination=ill)
            jobs.append((f"C  P40/20nm D0={d0:.0f} {ill}", kw, 100, 1500))

    nproc = min(os.cpu_count() or 2, 8)
    print(f"{len(jobs)} runs on {nproc} processes", flush=True)
    with Pool(nproc) as pool:
        out = pool.map(run, jobs, chunksize=1)
    print(f"{'case':<46}{'quencher':>9}{'ok?':>5}{'LCDU3s':>8}{'meanCD':>8}{'miss':>8}{'merge':>8}")
    for name, kw, q, s in out:
        inr = "yes" if 0.012 <= q <= 0.09 else "NO"
        print(f"{name:<46}{q:9.4f}{inr:>5}{s['lcdu_3s']:8.2f}{s['mean_cd']:8.2f}{s['p_missing']:8.4f}{s['p_merged']:8.4f}", flush=True)
