# euv-stochastic-simulator

CS 427 Project 3 (Hood College) — Ryan Sapp

A Monte Carlo simulator of how contact holes print under EUV lithography, plus machine-learning surrogate models (linear regression, random forest, neural network) that predict failure rate and hole-size variation from process settings.

**Simulator v2** follows the structure of published stochastic resist models: partially coherent NA 0.33 optics, photon shot noise, PAG, quencher and protecting-group counting, electron blur, acid and quencher diffusion, and threshold development. Its parameters come from published values. After two independent verification passes, hole-size variation for typical resist processes lands within about ±30% of published numbers (−14% to +29% across eight comparisons, including +2% against the PROLITH simulator), and its Monte Carlo agrees with a separate semi-analytical noise model within 1–3% for photon and acid noise. Failure rates are order-of-magnitude estimates, probably optimistic on missing holes. See [`results/VERIFICATION.md`](results/VERIFICATION.md) and [`results/CALIBRATION.md`](results/CALIBRATION.md).

## Status

- [x] Prototype simulator (v1) and failure-cliff sweep
- [x] Simulator v2 with literature parameters
- [x] Calibration against published data (`results/CALIBRATION.md`)
- [x] Verification, two passes: sources, real-world data, another simulator, a semi-analytical check (`results/VERIFICATION.md`)
- [x] Unit tests (`sim/test_simulator_v2.py`)
- [x] Dataset generator with built-in quality checks; pilot passed (`results/DATASET.md`)
- [ ] Full dataset generation (Latin hypercube over process settings)
- [ ] Model comparison: linear regression, random forest, neural network
- [ ] Report and presentation

## Layout

```
sim/
  simulator_v2.py        v2 model (use this)
  calibrate.py           runs all calibration targets -> results/calibration.json
  verify.py              convergence, sensitivity, corrected Goethals -> results/verification.json
  extra_checks.py        more real-world comparisons -> results/extra_checks.json
  window_check.py        failure-free dose window at 40 nm pitch -> results/window_check_L1.json
  analytic_check.py      semi-analytical noise model vs Monte Carlo -> results/analytic_check.json
  pass2_checks.py        pass 2 comparisons (Severi pitches, PROLITH, ...) -> results/pass2/*.json
  compare_passes.py      pass 1 vs pass 2 table (pass 1 outputs in results/pass1/)
  make_dataset.py        ML dataset generator (resumable, multi-core) -> data/<name>/
  check_dataset.py       dataset quality checks -> data/<name>/qa_report.md
  replicate_check.py     re-runs dataset rows to test regeneration and error bars
  optics_check.py, corner_check.py, dispersion_check.py   final pre-dataset checks
  test_make_dataset.py   generator tests
  make_figures.py        calibration figures -> results/*.png
  make_verification_figure.py   parity plot -> results/verification_parity.png
  test_simulator_v2.py   sanity tests
  simulator.py           v1 prototype (kept for reference)
  sweep.py, lcdu_check.py, stages.py   v1 scripts
results/      figures, calibration report, small CSV/JSON
notebooks/    exploration, calibration, modeling
report/       written report
data/         large generated datasets (not committed)
```

## Run

```bash
pip install -r requirements.txt
python sim/test_simulator_v2.py   # ~35 s
python sim/simulator_v2.py        # one quick dose-to-size example
python sim/calibrate.py           # ~5 min
python sim/make_figures.py

# dataset (see results/DATASET.md)
cd sim && python make_dataset.py --out ../data/full --n-resists 1500 --workers 8
python check_dataset.py ../data/full
```

Scripts save to `results/` wherever they are run from.

## Quick use

```python
import numpy as np
from simulator_v2 import Params, simulate, summarize, dose_to_size

p = Params(cd_target=20, pitch_x=40, pitch_y=40, dose=30)   # 3x3 hole array per trial
p = dose_to_size(p)                  # fit quencher so holes print at 20 nm at this dose
s = summarize(simulate(p, 400, np.random.default_rng(0)))
print(s)   # mean_cd, sigma_cd, lcdu_3s, p_missing, p_merged
```
