"""Run simulator v2 against published calibration targets.

Each target uses the published geometry and dose. The only fitted quantity is
the resist's speed (quencher loading, or the deprotection rate for the model
cross-check), chosen so holes print at the published size, exactly as a
fab picks a resist formulation for a given dose. Hole-size variation, failure
rates and trends are then predictions, not fits.

Writes results/calibration.json and results/calibration_*.png.
"""
import os
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
import json
import time
from dataclasses import replace, asdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulator_v2 import Params, simulate, summarize, dose_to_size

RNG = 2026
N = 600          # trials per final evaluation (x holes per tile)
results = {}


def run(p, n=N, seed=RNG):
    return summarize(simulate(p, n, np.random.default_rng(seed)))


def log(name, **kw):
    results[name] = kw
    print(name, json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in kw.items()}), flush=True)


t0 = time.time()

# ---- Target 1: Severi et al. 2021 (imec / Fractilia) --------------------------------
# 24 nm holes, staggered array 80 x 46 nm pitch, 80 mJ/cm2; stochastic-only LCDU ~2.6 nm (3 sigma)
p1 = Params(cd_target=24, pitch_x=80, pitch_y=46, stagger=True, cells_x=2, cells_y=4,
            illumination="annular", dose=80)
p1 = dose_to_size(p1)
s1 = run(p1)
log("T1_severi2021", published_lcdu=2.6, rho_q=p1.rho_q, **s1)

# ---- Target 2: Wang et al. 2019 (PSI, EUV interference lithography) ------------------
# best organic CAR: HP24 (pitch 48) LCDU 2.3 nm at 16.4 mJ/cm2; HP20 (pitch 40) LCDU 2.2 nm at 21.1
for tag, pitch, dose, lcdu in [("T2a_wang2019_hp24", 48, 16.4, 2.3), ("T2b_wang2019_hp20", 40, 21.1, 2.2)]:
    p = Params(cd_target=pitch / 2, pitch_x=pitch, pitch_y=pitch, dose=dose)
    p = dose_to_size(p)
    log(tag, published_lcdu=lcdu, rho_q=p.rho_q, **run(p))

# ---- Target 3: Goethals et al. 2013 (imec, NXE:3100, NA 0.25, 60 nm film, quasar) ------
# C31P52 = 31 nm mask holes, 52 nm pitch, printing 26 nm holes: dose-to-size 31.7 mJ/cm2,
# 1-sigma LCDU 1.1 nm. C29P48 = 29 nm mask holes, 48 nm pitch, printing 24 nm: 33.5 mJ/cm2,
# 1-sigma 1.7 nm ("no process window"). LCDU totals include metrology and mask error.
for tag, cd, mask, pitch, dose, s1sig in [("T3a_goethals2013_c31p52", 26, 31, 52, 31.7, 1.1),
                                           ("T3b_goethals2013_c29p48", 24, 29, 48, 33.5, 1.7)]:
    p = Params(cd_target=cd, mask_cd=mask, pitch_x=pitch, pitch_y=pitch, dose=dose, na=0.25, thickness=60)
    p = dose_to_size(p)
    log(tag, published_lcdu=3 * s1sig, rho_q=p.rho_q, **run(p))

# ---- Target 4: Naulleau et al. model cross-check -------------------------------------
# 16 nm half-pitch contacts, NA 0.33 quadrupole pole 0.1, their CAR parameter set:
# mean CD 17.17 nm, sigma 1.15 nm (2.6M contacts). The paper does not state the dose for
# this result; 82 mJ/cm2 is the mean dose of its later dose-noise analysis (an assumption).
p4 = Params(cd_target=17.17, pitch_x=32, pitch_y=32, dose=82, qe=3.0, rho_q=0.085,
            sigma_e=1.7, sigma_a=5.1, sigma_q=1.3, cells_x=4, cells_y=4)
p4 = dose_to_size(p4, knob="k_amp")
s4 = run(p4)
log("T4_naulleau_model", published_sigma=1.15, published_lcdu=3.45, k_amp=p4.k_amp, **s4)

# ---- Target 5: Yen 2018 (ASML stochastic CAR *simulation*, 10,000 holes per dose) ------
# CD 15 -> 17.4 -> 19.5 -> 21.2 nm at 30 -> 33 -> 36 -> 39 mJ/cm2 (slope ~0.69 nm per mJ/cm2);
# ~10% of holes close at CD ~15 nm. Pitch, film and resist not given; we use 40 nm pitch,
# resist sized at 36 mJ/cm2.
p5 = dose_to_size(Params(cd_target=19.5, pitch_x=40, pitch_y=40, dose=36))
rows = []
for dose in [27, 30, 33, 36, 39]:
    s = run(replace(p5, dose=dose), n=400)
    rows.append(dict(dose=dose, **s))
slope = np.polyfit([r["dose"] for r in rows[1:]], [r["mean_cd"] for r in rows[1:]], 1)[0]
log("T5_yen2018", published_slope=0.69, slope=float(slope), rho_q=p5.rho_q, sweep=rows)

# ---- Target 6: double cliff (imec): dense array, wide dose sweep ---------------------
p6 = dose_to_size(Params(cd_target=18, pitch_x=36, pitch_y=36, dose=40))
cliff = []
for dose in [24, 28, 32, 36, 40, 46, 52, 60, 70, 80]:
    s = run(replace(p6, dose=dose), n=400)
    cliff.append(dict(dose=dose, **s))
log("T6_double_cliff", rho_q=p6.rho_q, sweep=cliff)

# ---- Target 7: size trend (Park 2021): LCDU grows as holes shrink --------------------
size_rows = []
for cd in [16, 18, 20, 24, 28]:
    p = dose_to_size(Params(cd_target=cd, pitch_x=2 * cd, pitch_y=2 * cd, dose=40), knob="dose")
    size_rows.append(dict(cd=cd, dose_to_size=p.dose, **run(p, n=400)))
log("T7_size_trend", sweep=size_rows)

# ---- Target 8: dose scaling (Mack 2019): sigma ~ 1/sqrt(dose) for fixed geometry -----
scale_rows = []
for dose in [20, 30, 45, 67, 100]:
    p = dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=40, dose=dose))
    scale_rows.append(dict(dose=dose, rho_q=p.rho_q, **run(p, n=400)))
d = np.log([r["dose"] for r in scale_rows]); s = np.log([r["sigma_cd"] for r in scale_rows])
log("T8_dose_scaling", loglog_slope=float(np.polyfit(d, s, 1)[0]), sweep=scale_rows)

results["_runtime_s"] = time.time() - t0
with open(RESULTS + "/calibration.json", "w") as f:
    json.dump(results, f, indent=1)
print("done in %.0f s" % results["_runtime_s"])
