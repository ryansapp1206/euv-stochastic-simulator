"""Extra real-world comparisons. Writes results/extra_checks.json.
F: ASML through-pitch CAR data at 0.33 NA (40 nm CAR film), LCDU vs pitch at stated doses.
   Conditions partly unknown (pattern type, hole CD); we assume square arrays, CD = pitch/2.
G: ASML CAR progress at P40 hexagonal CH: LCDU 2.7/2.2/2.1/2.0 nm at 66/65/61/53 mJ/cm2.
   Approximated as a staggered array 40 x 35 nm (near-hexagonal), annular illumination.
H: Wang 2019 interference lithography re-test with an interference image and a 25 nm film.
I: van Ingen Schenau 2013, NXE:3300: 24 nm HP regular CH, LCDU 3.2 nm (small annular, 31 mJ/cm2)."""
import os
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
import json, time
from dataclasses import replace
import numpy as np
from simulator_v2 import Params, simulate, summarize, dose_to_size

out = {}
def run(p, n=400, seed=21):
    return summarize(simulate(p, n, np.random.default_rng(seed)))
def log(k, v):
    out[k] = v; print(k, json.dumps(v, default=float), flush=True)
t0 = time.time()

for pitch, dose, lcdu in [(40, 63, 2.4), (36, 66, 2.6), (34, 66, 2.7), (32, 72, 3.0), (30, 73, 3.6)]:
    p = dose_to_size(Params(cd_target=pitch / 2, pitch_x=pitch, pitch_y=pitch, dose=dose, thickness=40))
    log(f"F_asml_P{pitch}", dict(published_lcdu=lcdu, rho_q=p.rho_q, **run(p)))

for yr, dose, lcdu in [(2018, 66, 2.7), (2021, 53, 2.0)]:
    p = dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=35, stagger=True, cells_x=3, cells_y=4,
                            illumination="annular", dose=dose, thickness=40))
    log(f"G_asml_hexP40_{yr}", dict(published_lcdu=lcdu, rho_q=p.rho_q, **run(p)))

for hp, dose, lcdu in [(24, 16.4, 2.3), (20, 21.1, 2.2)]:
    p = dose_to_size(Params(cd_target=hp, pitch_x=2 * hp, pitch_y=2 * hp, dose=dose,
                            illumination="interference", thickness=25))
    log(f"H_wang_IL_hp{hp}", dict(published_lcdu=lcdu, rho_q=p.rho_q, **run(p)))

p = dose_to_size(Params(cd_target=24, pitch_x=48, pitch_y=48, dose=31, illumination="annular",
                        sigma_in=0.1, sigma_out=0.4, thickness=40))
log("I_vanIngenSchenau2013_24hp_small_annular", dict(published_lcdu=3.2, rho_q=p.rho_q, **run(p)))

out["_runtime_s"] = time.time() - t0
json.dump(out, open(RESULTS + "/extra_checks.json", "w"), indent=1, default=float)
print("done in %.0f s" % out["_runtime_s"])
