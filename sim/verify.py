"""Verification runs: convergence, sensitivity, corrected Goethals targets, photon counts.
Writes results/verification.json."""
import json, time
from dataclasses import replace
import numpy as np
from simulator_v2 import Params, simulate, summarize, dose_to_size, aerial_image, PHOTONS_PER_NM2_PER_MJCM2

out = {}
def run(p, n=500, seed=11):
    return summarize(simulate(p, n, np.random.default_rng(seed)))
def log(k, v):
    out[k] = v; print(k, json.dumps(v, default=float), flush=True)
t0 = time.time()

# ---------------- A. numerical convergence (20 nm holes, 40 nm pitch, 30 mJ/cm2)
base = dose_to_size(Params())
ref = run(base)
log("A0_reference_px1", dict(rho_q=base.rho_q, **ref))
fine = replace(base, px=0.5)
log("A1_px0.5_same_resist", run(fine, n=300))
fine_fit = dose_to_size(replace(Params(), px=0.5), n_trials=80)
log("A2_px0.5_refit", dict(rho_q=fine_fit.rho_q, **run(fine_fit, n=300)))
log("A3_pole_samples_11", run(replace(base, pole_samples=11)))
log("A4_tile_5x5", run(replace(base, cells_x=5, cells_y=5), n=200))
log("A5_seed_repeat", [run(base, n=500, seed=s)["lcdu_3s"] for s in (1, 2, 3)])

# ---------------- B. sensitivity of the Severi 2021 match (refit quencher each time)
sev = Params(cd_target=24, pitch_x=80, pitch_y=46, stagger=True, cells_x=2, cells_y=4,
             illumination="annular", dose=80)
variants = {
    "baseline": {},
    "acid_blur_4.0": dict(sigma_a=4.0), "acid_blur_6.0": dict(sigma_a=6.0),
    "electron_blur_1.2": dict(sigma_e=1.2), "electron_blur_2.5": dict(sigma_e=2.5),
    "qe_2.0": dict(qe=2.0), "qe_3.5": dict(qe=3.5),
    "pag_0.12": dict(rho_pag=0.12), "pag_0.30": dict(rho_pag=0.30),
    "protect_1.5": dict(rho_pg=1.5), "protect_3.0": dict(rho_pg=3.0),
    "dev_threshold_0.3": dict(dev_threshold=0.3), "dev_threshold_0.7": dict(dev_threshold=0.7),
    "absorptivity_5.3": dict(absorptivity=0.0053), "thickness_30": dict(thickness=30.0),
    "annular_0.3_0.7": dict(sigma_in=0.3, sigma_out=0.7),
    "layers_5": dict(layers=5),
}
sens = {}
for name, kw in variants.items():
    p = dose_to_size(replace(sev, **kw))
    s = run(p, n=300)
    sens[name] = dict(lcdu=s["lcdu_3s"], rho_q=p.rho_q)
    print("B", name, round(s["lcdu_3s"], 3), round(p.rho_q, 4), flush=True)
out["B_sensitivity_severi"] = sens

# ---------------- C. corrected Goethals 2013 (NXE:3100, NA 0.25, 60 nm film, quasar)
for tag, cd, mask, pitch, dose, s1 in [("C1_goethals_26nm_C31P52", 26, 31, 52, 31.7, 1.1),
                                        ("C2_goethals_24nm_C29P48", 24, 29, 48, 33.5, 1.7)]:
    p = dose_to_size(Params(cd_target=cd, mask_cd=mask, pitch_x=pitch, pitch_y=pitch, dose=dose,
                            na=0.25, thickness=60))
    log(tag, dict(published_lcdu=3 * s1, rho_q=p.rho_q, **run(p)))

# ---------------- E. absorbed photons per hole vs Yen 2018 slide (961 photons, 27 nm cylinder?)
p5 = dose_to_size(Params(cd_target=19.5, pitch_x=40, pitch_y=40, dose=36))
I = aerial_image(replace(p5, dose=30))
H, W = I.shape
yy, xx = np.mgrid[:H, :W] + 0.5
r = np.hypot(xx - 20, yy - 20)
per_nm2 = 30 * PHOTONS_PER_NM2_PER_MJCM2 * p5.absorbed_frac * I
log("E_photons_30mJ", dict(cyl27=float(per_nm2[r <= 13.5].sum()), cyl15=float(per_nm2[r <= 7.5].sum()),
                           yen_slide=[961, 378]))

out["_runtime_s"] = time.time() - t0
json.dump(out, open("../results/verification.json", "w"), indent=1, default=float)
print("done in %.0f s" % out["_runtime_s"])
