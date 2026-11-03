# Dose-to-size LCDU (3 sigma of CD) for isolated holes, prototype model, to compare with published numbers
import numpy as np
from simulator import simulate, calibrate_threshold
rng = np.random.default_rng(7)
for d, dose in [(24, 80), (24, 33), (20, 21), (24, 16.4), (20, 30)]:
    thr = calibrate_threshold(ref_dose=dose, ref_d=d)
    cd, fail = simulate(dose, d, 4000, thr, rng=rng)
    print(f"hole {d} nm @ {dose} mJ/cm2: mean CD {cd.mean():.1f}, LCDU(3s) {3*cd.std():.2f} nm, fail {fail.mean():.4f}")
