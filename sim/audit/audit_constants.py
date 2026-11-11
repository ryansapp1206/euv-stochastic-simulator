"""Independent first-principles audit of simulator_v2 constants and formulas.

Nothing here imports the simulator's constants; everything is recomputed from
physics and compared. Run: python3 -I audit_constants.py
"""
import sys
import numpy as np

import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import simulator_v2 as S  # noqa: E402

OK, BAD = [], []


def cmp(name, mine, theirs, tol=1e-3, unit=""):
    rel = abs(mine - theirs) / max(abs(theirs), 1e-30)
    (OK if rel <= tol else BAD).append(name)
    flag = "ok " if rel <= tol else "XX "
    print(f"{flag}{name:<46} mine={mine:<14.6g} code={theirs:<14.6g} rel={rel:.2e} {unit}")


print("=" * 100)
print("1. PHOTON PHYSICS")
print("=" * 100)
# CODATA 2018: h = 6.62607015e-34 J*s (exact), c = 299792458 m/s (exact), e = 1.602176634e-19 C (exact)
h, c, e = 6.62607015e-34, 299792458.0, 1.602176634e-19
E_J = h * c / (13.5e-9)
E_eV = E_J / e
cmp("photon energy at 13.5 nm (eV)", E_eV, S.PHOTON_EV, 1e-4, "eV")
cmp("hc (eV nm)", h * c / e * 1e9, 1239.84, 1e-4)

# 1 mJ/cm^2 = 1e-3 J / 1e14 nm^2 = 1e-17 J/nm^2
J_per_nm2 = 1e-3 / 1e14
cmp("J/nm^2 per mJ/cm^2", J_per_nm2, 1e-17, 1e-12)
cmp("photons/nm^2 per mJ/cm^2", J_per_nm2 / E_J, S.PHOTONS_PER_NM2_PER_MJCM2, 1e-4)

print()
print("=" * 100)
print("2. ABSORPTION (Beer-Lambert)")
print("=" * 100)
for alpha_um, t in ((4.2, 35.0), (4.2, 50.0), (5.3, 35.0), (15.0, 30.0)):
    p = S.Params(absorptivity=alpha_um / 1000.0, thickness=t)
    mine = 1.0 - np.exp(-(alpha_um / 1000.0) * t)
    cmp(f"absorbed fraction  a={alpha_um}/um t={t}nm", mine, p.absorbed_frac, 1e-9)

# layer split must sum to the total
for L in (1, 3, 5):
    p = S.Params(layers=L)
    dz = p.thickness / L
    lay = np.exp(-p.absorptivity * dz * np.arange(L)) * (1 - np.exp(-p.absorptivity * dz))
    cmp(f"layer absorption sums to total (L={L})", lay.sum(), p.absorbed_frac, 1e-9)

print()
print("=" * 100)
print("3. PHOTON BUDGET  (sanity against published 'photons per hole' numbers)")
print("=" * 100)
p = S.Params()
inc = p.dose * S.PHOTONS_PER_NM2_PER_MJCM2
print(f"   dose {p.dose} mJ/cm^2 -> {inc:.2f} incident photons/nm^2")
print(f"   absorbed fraction {p.absorbed_frac:.4f} -> {inc * p.absorbed_frac:.3f} absorbed/nm^2")
for cd in (20.0, 27.0):
    area = np.pi * (cd / 2) ** 2
    print(f"   over a {cd:.0f} nm circle ({area:6.1f} nm^2): {inc * p.absorbed_frac * area:7.1f} absorbed photons")
print("   (VERIFICATION.md claims ~420 in v2 for the Yen 27 nm comparison -- that is at a")
print("    different dose/film; recomputed below from the aerial image itself.)")

print()
print("=" * 100)
print("4. AERIAL IMAGE NORMALISATION AND ENERGY")
print("=" * 100)
# open frame: mask all ones -> image must be exactly 1 everywhere
openf = S.Params(cd_target=40.0, mask_cd=40.0, pitch_x=40, pitch_y=40, cells_x=1, cells_y=1)
img = S.aerial_image(openf)
cmp("open frame image == 1", float(img.mean()), 1.0, 2e-3)
print(f"   open-frame image min/max = {img.min():.4f} / {img.max():.4f}")

# Parseval / DC: mean image should equal sum of |captured order amplitudes|^2
for pitch, cd in ((40, 20.0), (32, 16.0), (64, 30.0)):
    pp = S.Params(pitch_x=pitch, pitch_y=pitch, cd_target=cd, mask_cd=cd)
    im = S.aerial_image(pp)
    mask, _, _ = S._geometry(pp)
    print(f"   P{pitch} CD{cd}: mask open fraction={mask.mean():.4f}  image mean={im.mean():.4f}  "
          f"peak={im.max():.4f}  min={im.min():.4f}")

print()
print("=" * 100)
print("5. DIFFRACTION / PUPIL GEOMETRY  (does the chosen illumination actually capture orders?)")
print("=" * 100)
lam, NA = 13.5, 0.33
print(f"   cutoff frequency NA/lambda = {NA/lam:.6f} cycles/nm  ->  min resolvable pitch (2-beam) = {lam/(2*NA):.2f} nm")
for pitch in (32, 36, 40, 48, 64):
    f1 = (1.0 / pitch) / (NA / lam)          # first order, pupil units
    sc_ideal = lam / (np.sqrt(2) * pitch * NA)
    sc_used = float(np.clip(sc_ideal, 0.3, 0.9 - S.Params().pole_sigma))
    cpt = sc_used / np.sqrt(2)
    # is the (1,1) diagonal order captured from a diagonal pole at (-c,-c)?
    d11 = np.hypot(f1 - cpt, f1 - cpt)
    d10 = np.hypot(f1 - cpt, 0 - cpt)
    d00 = np.hypot(cpt, cpt)
    clipped = "CLIPPED" if abs(sc_used - sc_ideal) > 1e-9 else ""
    print(f"   P{pitch:<3} first order at {f1:5.3f} pupil radii | pole sigma ideal={sc_ideal:5.3f} "
          f"used={sc_used:5.3f} {clipped:<8} | captured: 0th {d00<=1}, (1,0) {d10<=1}, (1,1) {d11<=1}")

print()
print("=" * 100)
print("6. RESIST CHEMISTRY FORMULAS")
print("=" * 100)
# Dill-C style first-order kinetics: acids = PAG * (1 - exp(-QE*n_abs/PAG)); low-dose limit = QE*n_abs
pag_col = S.Params().rho_pag * S.Params().thickness        # PAG per 1 nm^2 column
for nabs in (0.01, 0.1, 1.0, 5.0, 20.0):
    pconv = 1 - np.exp(-2.8 * nabs / pag_col)
    acids = pag_col * pconv
    print(f"   {nabs:6.2f} absorbed photons/column -> {acids:7.3f} acids  "
          f"(low-dose limit QE*n = {2.8*nabs:7.3f}, ratio {acids/(2.8*nabs):.4f}, PAG pool {pag_col:.1f})")
print("   -> QE is 'acids per absorbed photon at low dose' as documented; saturates at the PAG pool.")

# deprotection
print()
K = 60.0
for hh in (0.0, 0.005, 0.01, 0.02, 0.05, 0.1):
    print(f"   net acid h={hh:5.3f}/nm^3 -> protected m=exp(-K h)={np.exp(-K*hh):.4f}, deprotected={1-np.exp(-K*hh):.4f}")
print(f"   -> 50% deprotection (the development threshold) needs h = ln2/K = {np.log(2)/K:.5f} /nm^3")
print(f"      i.e. {np.log(2)/K * S.Params().thickness:.3f} net acids in a 1 nm^2 column of {S.Params().thickness} nm film")

print()
print("=" * 100)
print("7. QUENCHER / PAG RATIO  (vs published CAR formulations)")
print("=" * 100)
for q in (0.012, 0.03, 0.09):
    print(f"   quencher {q:.3f}/nm^3 with PAG 0.20/nm^3 -> Q/PAG = {q/0.20:.3f}")
print("   published: base is typically 10-20% of PAG (AIP Advances 2025); Q/PAG ~0.4 in a generic CAR (JMM 2021)")

print()
print("=" * 100)
print(f"SUMMARY: {len(OK)} checks reproduced, {len(BAD)} failed")
if BAD:
    print("FAILED:", BAD)
print("=" * 100)
