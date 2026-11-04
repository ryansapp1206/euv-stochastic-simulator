"""Pass 2 cross-verification: new external comparisons and fresh-seed reruns.

P1  Severi et al. 2021 (SPIE 11609, Tables I, III, IV): stochastic-only LCDU of
    24 nm holes in three staggered arrays H80V46, H70V40, H66V38 at 80 mJ/cm^2,
    NXE:3400, custom hexapole, 50 nm CAR.  Measured (two runs):
    2.64/2.80, 2.29/2.48, 2.35/2.55 nm.  All three pitch pairs have the ratio
    sqrt(3), i.e. hexagonal arrays with centre-to-centre 46 / 40 / 38 nm
    (Fractilia 2023 describes the same imec layout this way).  Pass 1 had read
    H80V46 as rows 46 nm apart, a lattice with half the hole density.
P2  'Line-edge roughness performance targets for EUV lithography' (SPIE 2017, co-author
    C. Mack; lithoguru.com PDF), a PROLITH stochastic simulation (a different simulator):
    dense holes on 50 nm centres, mask 24 nm printing 23 nm at 72 mJ/cm^2,
    0.33 NA small annular 0.292-0.643, 30 nm CAR: LCDU 3.12 nm, no fails at
    mask CD >= 23.5 nm, fails rising below, MEEF > 10 below 23 nm.
P3  ASML hexagonal P40 CAR data (2018 2.7 nm at 66 mJ/cm^2, 2021 2.0 nm at
    53 mJ/cm^2) re-run with hexapole illumination instead of annular.
P4  Key pass-1 comparisons re-run with new random seeds and 2.5x more holes.

Usage: python pass2_checks.py P1 P2 ...  ->  ../results/pass2/<section>.json
"""
import os
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
import json
import time
from dataclasses import replace

import numpy as np

from simulator_v2 import Params, simulate, summarize, dose_to_size

import sys

out = {}


def run(p, n, seed):
    return summarize(simulate(p, n, np.random.default_rng(seed)))


def log(k, v):
    out[k] = v
    print(k, json.dumps(v, default=lambda x: round(float(x), 4)), flush=True)


def P1():
    # P1: Severi 2021, three hexagonal pitches ------------------------------------------
    severi = {"H80V46": (46, 40, (2.64, 2.80)), "H70V40": (40, 35, (2.29, 2.48)), "H66V38": (38, 33, (2.35, 2.55))}
    for tag, (px_, py_, meas) in severi.items():
        for ill in ("hexapole", "annular"):
            p = dose_to_size(Params(cd_target=24, pitch_x=px_, pitch_y=py_, stagger=True, cells_x=3, cells_y=4,
                                    illumination=ill, dose=80, thickness=50))
            log(f"P1_severi_{tag}_{ill}", dict(measured=meas, rho_q=p.rho_q, **run(p, 500, 31)))
    # pass-1 reading of H80V46 (rows 46 nm apart, 35 nm film, annular) for reference
    p = dose_to_size(Params(cd_target=24, pitch_x=80, pitch_y=46, stagger=True, cells_x=2, cells_y=4,
                            illumination="annular", dose=80))
    log("P1_severi_H80V46_pass1_geometry", dict(measured=severi["H80V46"][2], rho_q=p.rho_q, **run(p, 500, 31)))


def P2():
    # P2: Mack 2017 PROLITH cross-simulator check ---------------------------------------
    pm = dose_to_size(Params(cd_target=23, mask_cd=24, pitch_x=50, pitch_y=50, dose=72, illumination="annular",
                             sigma_in=0.292, sigma_out=0.643, thickness=30))
    rows = []
    for mask in (22.0, 22.5, 23.0, 23.5, 24.0, 24.5, 25.0):
        s = run(replace(pm, mask_cd=mask), 600, 41)
        rows.append(dict(mask_cd=mask, **s))
        print("   mask", mask, {k: round(v, 4) for k, v in s.items()}, flush=True)
    cds = [r["mean_cd"] for r in rows]
    meef = [(cds[i + 1] - cds[i]) / 0.5 for i in range(len(cds) - 1)]
    log("P2_mack2017_prolith", dict(published_lcdu=3.12, rho_q=pm.rho_q, sweep=rows, meef=meef))


def P3():
    # P3: ASML hexagonal P40 with hexapole -----------------------------------------------
    for yr, dose, lcdu in [(2018, 66, 2.7), (2021, 53, 2.0)]:
        p = dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=35, stagger=True, cells_x=3, cells_y=4,
                                illumination="hexapole", dose=dose, thickness=40))
        log(f"P3_asml_hexP40_{yr}_hexapole", dict(published_lcdu=lcdu, rho_q=p.rho_q, **run(p, 500, 51)))


def P4():
    # P4: fresh seeds, 2.5x more holes ----------------------------------------------------
    fresh = {
        "vanIngenSchenau2013": (Params(cd_target=24, pitch_x=48, pitch_y=48, dose=31, illumination="annular",
                                       sigma_in=0.1, sigma_out=0.4, thickness=40), 3.2, "rho_q"),
        "goethals2013_C31P52": (Params(cd_target=26, mask_cd=31, pitch_x=52, pitch_y=52, dose=31.7, na=0.25,
                                       thickness=60), 3.3, "rho_q"),
        "asml_P40": (Params(cd_target=20, pitch_x=40, pitch_y=40, dose=63, thickness=40), 2.4, "rho_q"),
        "asml_P36": (Params(cd_target=18, pitch_x=36, pitch_y=36, dose=66, thickness=40), 2.6, "rho_q"),
        "asml_P34": (Params(cd_target=17, pitch_x=34, pitch_y=34, dose=66, thickness=40), 2.7, "rho_q"),
        "asml_P32": (Params(cd_target=16, pitch_x=32, pitch_y=32, dose=72, thickness=40), 3.0, "rho_q"),
        "naulleau_model": (Params(cd_target=17.17, pitch_x=32, pitch_y=32, dose=82, qe=3.0, rho_q=0.085,
                                  sigma_a=5.1, cells_x=4, cells_y=4), 3.45, "k_amp"),
    }
    for i, (tag, (p, pub, knob)) in enumerate(fresh.items()):
        p = dose_to_size(p, knob=knob, rng_seed=900 + i)
        n = 1500 if p.cells_x * p.cells_y == 9 else 900
        log(f"P4_{tag}", dict(published_lcdu=pub, knob=knob, fitted=getattr(p, knob), **run(p, n, 7000 + i)))


def P1b():
    # Severi, pass-1 reading of the layout (rows V nm apart, holes H nm apart in a row), 50 nm film
    meas = {"H80V46": (80, 46, (2.64, 2.80)), "H70V40": (70, 40, (2.29, 2.48)), "H66V38": (66, 38, (2.35, 2.55))}
    for tag, (h, v, m) in meas.items():
        for ill in ("hexapole", "annular"):
            p = dose_to_size(Params(cd_target=24, pitch_x=h, pitch_y=v, stagger=True, cells_x=2, cells_y=4,
                                    illumination=ill, dose=80, thickness=50))
            log(f"P1b_severi_{tag}_rows_{ill}", dict(measured=m, rho_q=p.rho_q, **run(p, 500, 32)))


def P1c():
    # Severi, hexagonal reading, hole size scaled with the hole spacing (the size at the two
    # tighter pitches is not stated; 24 nm is only given for H80V46)
    for tag, (px_, py_, cc, m) in {"H70V40": (40, 35, 40.31, (2.29, 2.48)),
                                   "H66V38": (38, 33, 38.08, (2.35, 2.55))}.items():
        cd = round(24 * cc / 46.14, 1)
        p = dose_to_size(Params(cd_target=cd, pitch_x=px_, pitch_y=py_, stagger=True, cells_x=3, cells_y=4,
                                illumination="hexapole", dose=80, thickness=50))
        log(f"P1c_severi_{tag}_cd{cd}_hexapole", dict(measured=m, cd=cd, rho_q=p.rho_q, **run(p, 500, 33)))


def P5():
    # Goethals 2013 (NXE:3100, NA 0.25, quasar, 60 nm): four resists printing 26 nm holes (C31P52)
    # dose-to-size / 1-sigma LCDU: 18.1/1.5, 19.4/1.2, 26.3/1.1, 31.7/1.1 (totals incl. mask+metrology;
    # the slide's column labels are garbled, the dose-LCDU pairs are as printed)
    for dose, s1 in [(18.1, 1.5), (19.4, 1.2), (26.3, 1.1), (31.7, 1.1)]:
        p = dose_to_size(Params(cd_target=26, mask_cd=31, pitch_x=52, pitch_y=52, dose=dose, na=0.25,
                                thickness=60))
        log(f"P5_goethals_{dose}", dict(published_sigma=s1, rho_q=p.rho_q, **run(p, 500, 34)))


def P6():
    # How LCDU scales with dose depends on HOW the resist is made faster or slower (Maas et al.
    # 2018, ASML: dose^-0.5 when the development threshold changes, dose^-0.33 when QE changes).
    # Compare our quencher knob with a threshold-like knob (deprotection rate K; h_th = ln2/K).
    for knob, base in (("rho_q", Params(cd_target=20, pitch_x=40, pitch_y=40)),
                       ("k_amp", Params(cd_target=20, pitch_x=40, pitch_y=40, rho_q=0.0156))):
        rows = []
        for dose in (20, 30, 45, 67, 100):
            p = dose_to_size(replace(base, dose=dose), knob=knob)
            rows.append(dict(dose=dose, fitted=getattr(p, knob), **run(p, 400, 35)))
        slope = float(np.polyfit(np.log([r["dose"] for r in rows]), np.log([r["sigma_cd"] for r in rows]), 1)[0])
        log(f"P6_dose_scaling_{knob}", dict(loglog_slope=slope, sweep=rows))


def P2b():
    # Same Mack 2017 mask sweep with the film split into 5 depth layers (3D blur, light weakening
    # with depth, hole must open at the bottom). Does 3D development make the small-CD cliff steeper?
    pm = dose_to_size(Params(cd_target=23, mask_cd=24, pitch_x=50, pitch_y=50, dose=72, illumination="annular",
                             sigma_in=0.292, sigma_out=0.643, thickness=30, layers=5), n_trials=60)
    rows = []
    for mask in (22.0, 22.5, 23.0, 23.5, 24.0):
        s = run(replace(pm, mask_cd=mask), 300, 42)
        rows.append(dict(mask_cd=mask, **s))
        print("   L5 mask", mask, {k: round(v, 4) for k, v in s.items()}, flush=True)
    cds = [r["mean_cd"] for r in rows]
    log("P2b_mack2017_layers5", dict(published_lcdu=3.12, rho_q=pm.rho_q, sweep=rows,
                                     meef=[(cds[i + 1] - cds[i]) / 0.5 for i in range(len(cds) - 1)]))


if __name__ == "__main__":
    os.makedirs(os.path.join(RESULTS, "pass2"), exist_ok=True)
    for name in sys.argv[1:]:
        out.clear()
        t0 = time.time()
        globals()[name]()
        out["_runtime_s"] = time.time() - t0
        json.dump(out, open(os.path.join(RESULTS, "pass2", f"{name}.json"), "w"), indent=1, default=float)
        print(name, "done in %.0f s" % out["_runtime_s"], flush=True)
