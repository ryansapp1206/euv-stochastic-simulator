"""Audit of what the size labels actually mean near the failure cliffs."""
import sys
from dataclasses import replace
import numpy as np
import pandas as pd

import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import simulator_v2 as S  # noqa: E402

R = pd.read_csv(os.environ.get("EUV_ROWS", "data/full/rows.csv"))
R["surv"] = R.n_ok / R.n_holes
R["fail_frac"] = 1 - R.surv

print("=" * 96)
print("1. HOW MUCH OF THE DATASET HAS SURVIVOR-CONDITIONED SIZE LABELS?")
print("=" * 96)
print(f"   rows: {len(R)}")
for t in (0.0, 0.01, 0.05, 0.20, 0.50):
    m = R.fail_frac > t
    print(f"   rows where > {100*t:4.0f}% of holes are missing/merged (excluded from sigma): {m.sum():6d} ({100*m.mean():5.1f}%)")
print(f"   rows with < 100 surviving holes (sigma unusable): {(R.n_ok < 100).sum()}")

print()
print("=" * 96)
print("2. SIZE-LABEL BIAS: sigma over survivors vs sigma over ALL holes, one resist swept in dose")
print("=" * 96)
base = S.Params(pitch_x=40, pitch_y=40, cd_target=20.0, dose=60.0)
base = S.dose_to_size(base, "rho_q", n_trials=100, rng_seed=3)
print(f"   P40 / 20 nm / D0 = 60 mJ/cm2, fitted quencher {base.rho_q:.4f}")
print(f"{'dose/D0':>8}{'missing':>9}{'merged':>8}{'mean CD (surv)':>16}{'sigma surv':>12}{'sigma all':>11}{'bias':>8}")
for f in (0.70, 0.78, 0.85, 0.92, 1.0, 1.4, 1.7, 1.9):
    q = replace(base, dose=base.dose * f)
    r = S.simulate(q, 1500, np.random.default_rng(5))
    ok = ~(r["missing"] | r["merged"])
    s_ok = r["cd"][ok].std() if ok.sum() > 5 else np.nan
    s_all = r["cd"].std()
    mcd = r["cd"][ok].mean() if ok.sum() > 5 else np.nan
    print(f"{f:8.2f}{r['missing'].mean():9.3f}{r['merged'].mean():8.3f}{mcd:16.2f}{s_ok:12.3f}{s_all:11.3f}"
          f"{(s_ok / s_all - 1) if s_all > 0 else np.nan:8.0%}")

print()
print("=" * 96)
print("3. DOSE COVERAGE vs REAL SCANNER USE")
print("=" * 96)
print(f"   D0 (sizing dose) range in accepted resists: {R.drop_duplicates('resist_id').d0_mjcm2.min():.1f} - "
      f"{R.drop_duplicates('resist_id').d0_mjcm2.max():.1f} mJ/cm2, median {R.drop_duplicates('resist_id').d0_mjcm2.median():.1f}")
res = R.drop_duplicates("resist_id")
for lo, hi in ((20, 30), (30, 45), (45, 60), (60, 90)):
    print(f"   D0 in [{lo},{hi}): {((res.d0_mjcm2 >= lo) & (res.d0_mjcm2 < hi)).sum():5d} resists "
          f"({100*((res.d0_mjcm2 >= lo) & (res.d0_mjcm2 < hi)).mean():4.1f}%)")
print(f"   absolute dose range across rows: {R.dose_mjcm2.min():.1f} - {R.dose_mjcm2.max():.1f} mJ/cm2")
print(f"   rows with dose > 100 mJ/cm2: {(R.dose_mjcm2 > 100).sum()} ({100*(R.dose_mjcm2>100).mean():.1f}%)")

print()
print("=" * 96)
print("4. NOMINAL-DOSE LCDU DISTRIBUTION vs PUBLISHED TYPICAL CAR (2-5 nm) AND BEST (2.0-2.7 nm)")
print("=" * 96)
nom = R[np.isclose(R.dose_factor, 1.0)]
print(f"   nominal-dose LCDU quantiles (nm): " + ", ".join(f"{q:.0%}={nom.lcdu_3s_nm.quantile(q):.2f}" for q in (0.05, 0.25, 0.5, 0.75, 0.95)))
print(f"   fraction below 2.7 nm (best published CAR region): {(nom.lcdu_3s_nm < 2.7).mean():.1%}")
print(f"   fraction below 2.0 nm: {(nom.lcdu_3s_nm < 2.0).mean():.2%}")
print(f"   fraction above 5.0 nm: {(nom.lcdu_3s_nm > 5.0).mean():.2%}")
