"""Independent check of the simulator's physical constants and per-stage formulas.

Each quantity is re-derived here from first principles (CODATA constants and the
published model equations) and compared with what simulator_v2 actually computes.
Randomness is switched off (Params.noise = ()), so every stage has an exact
expected value and the comparison is deterministic.

Usage: python constants_check.py  ->  ../results/constants_check.json
Exit code 1 if any check fails.
"""
import json
import os
import sys
from dataclasses import replace

import numpy as np
from scipy.ndimage import label

import simulator_v2 as sv
from simulator_v2 import Params, simulate, aerial_image

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")

# CODATA 2018 exact values, typed in here rather than taken from the simulator
H_PLANCK = 6.62607015e-34      # J s
C_LIGHT = 299792458.0          # m / s
E_CHARGE = 1.602176634e-19     # C (J per eV)
LAMBDA_M = 13.5e-9             # EUV wavelength

checks = []


def check(name, expected, got, rtol=1e-4, atol=0.0, note=""):
    expected, got = np.asarray(expected, float), np.asarray(got, float)
    ok = bool(np.allclose(got, expected, rtol=rtol, atol=atol))
    err = float(np.max(np.abs(got - expected) / np.maximum(np.abs(expected), 1e-12)))
    checks.append(dict(check=name, ok=ok, expected=float(np.mean(expected)), got=float(np.mean(got)),
                       max_rel_error=err, note=note))


def fields(p, n=2):
    res = simulate(p, n, rng=np.random.default_rng(0), return_fields=True)
    return res, res["fields"]


# deterministic base case: no noise, no blur, one layer, 1 nm pixels
# (quencher and dose chosen so the noise-free holes open at about the 20 nm target)
base = Params(cd_target=20, pitch_x=40, pitch_y=40, dose=50.0, thickness=35.0, absorptivity=0.0042,
              qe=2.8, rho_pag=0.20, rho_q=0.06, rho_pg=2.0, sigma_e=0.0, sigma_a=0.0, sigma_q=0.0,
              k_amp=60.0, dev_threshold=0.5, layers=1, px=1.0, noise=())
t = base.thickness

# 1. photon energy E = h c / lambda
e_photon_ev = H_PLANCK * C_LIGHT / LAMBDA_M / E_CHARGE
check("photon energy at 13.5 nm (eV)", e_photon_ev, sv.PHOTON_EV, note="h c / lambda, CODATA")

# 2. photons per nm^2 per mJ/cm^2: 1 mJ/cm^2 = 1e-3 J / 1e14 nm^2
phi = 1e-3 / 1e14 / (e_photon_ev * E_CHARGE)
check("photons per nm^2 per mJ/cm^2", phi, sv.PHOTONS_PER_NM2_PER_MJCM2)

# 3. Beer-Lambert absorbed photons per pixel = dose * phi * I * (1 - exp(-alpha t))
_, f = fields(base)
expected_photons = base.dose * phi * f["image"] * (1 - np.exp(-base.absorptivity * t))
check("absorbed photons (Beer-Lambert)", expected_photons, f["photons"])

# 4. splitting the film into layers must not change the total absorbed
_, f4 = fields(replace(base, layers=4))
check("layer sum equals single-layer absorption", f["photons"], f4["photons"])

# 5. open-frame normalisation: a fully open mask images to exactly 1
open_img = aerial_image(replace(base, mask_cd=40.0))
check("open frame images to 1.0", np.ones_like(open_img), open_img, atol=1e-6)

# 6. acid generation, first-order exposure kinetics (Dill C):
#    acid = rho_PAG t (1 - exp(-QE n_abs / (rho_PAG t)))
pag0 = base.rho_pag * t
expected_acid = pag0 * (1 - np.exp(-base.qe * f["photons"] / pag0))
check("acid generation (Dill kinetics)", expected_acid, f["acid"])

# 7. low-dose limit: acids per absorbed photon -> QE
_, flo = fields(replace(base, dose=1e-3))
ratio = flo["acid"].sum() / flo["photons"].sum()
check("low-dose limit: acids per photon = QE", base.qe, ratio, rtol=1e-3)

# 8. high-dose saturation: every PAG converts, acid -> rho_PAG t where light is strong
_, fhi = fields(replace(base, dose=1e5))
bright = fhi["image"] > 0.5 * fhi["image"].max()
check("high-dose saturation: acid = rho_PAG t", np.full(bright.sum(), pag0), fhi["acid"][bright], rtol=1e-3)

# 9. quencher neutralisation: net acid per nm^3 = max(acid - rho_Q t, 0) / t
expected_h = np.maximum(f["acid"] - base.rho_q * t, 0) / t
check("quencher neutralisation (net acid)", expected_h, f["net_acid"], atol=1e-7)

# 10. deprotection: deprotected fraction = 1 - exp(-K h)
expected_frac = 1 - np.exp(-base.k_amp * f["net_acid"])
check("deprotection 1 - exp(-K h)", expected_frac, f["deprotected"], atol=1e-6)

# 11. development: a column dissolves where deprotected fraction > threshold
cleared_expected = f["deprotected"] > base.dev_threshold
agree = float((cleared_expected == f["clear"]).mean())
check("development threshold", 1.0, agree, rtol=0, atol=0, note="fraction of pixels that agree")

# 12. printed CD = equivalent-circle diameter of each hole's dissolved area
res, f = fields(base, n=1)
lab, n = label(f["clear"])
areas = np.sort(np.bincount(lab.ravel())[1:])
expected_cd = np.sort(2 * np.sqrt(areas / np.pi))
check("printed CD = 2 sqrt(area / pi)", expected_cd, np.sort(res["cd"][:n]),
      note=f"{n} holes in the tile")

n_ok = sum(c["ok"] for c in checks)
for c in checks:
    print(f"{'PASS' if c['ok'] else 'FAIL'}  {c['check']:44s} expected {c['expected']:.6g}  "
          f"got {c['got']:.6g}  max rel err {c['max_rel_error']:.1e}")
print(f"\n{n_ok} of {len(checks)} match")
os.makedirs(RESULTS, exist_ok=True)
json.dump(dict(passed=n_ok, total=len(checks), checks=checks), open(RESULTS + "/constants_check.json", "w"), indent=1)
sys.exit(0 if n_ok == len(checks) else 1)
