import numpy as np, csv
from simulator import simulate, calibrate_threshold
thr = calibrate_threshold(); rng = np.random.default_rng(1)
N = 2000
with open('results/sweep.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['size_nm', 'dose', 'fail_rate', 'fails', 'mean_cd', 'sd_cd'])
    for d in (16, 18, 20, 24):
        for dose in range(12, 37, 2):
            cd, fail = simulate(dose, d, N, thr, rng=rng)
            w.writerow([d, dose, fail.mean(), int(fail.sum()), round(cd.mean(), 2), round(cd.std(), 3)])
            print(d, dose, fail.mean(), round(cd.mean(),1), round(cd.std(),2))
