# euv-stochastic-simulator

CS 427 Project 3 (Hood College) — Ryan Sapp

A Monte Carlo simulator of how contact holes print under EUV lithography, plus machine-learning surrogate models (linear regression, random forest, neural network) that predict failure rate and hole-size variation from process settings.

**Simulator v2** follows the structure of published stochastic resist models: partially coherent NA 0.33 optics, photon shot noise, PAG, quencher and protecting-group counting, electron blur, acid and quencher diffusion, and threshold development. Its parameters come from published values. After two independent verification passes, hole-size variation for typical resist processes lands within about ±30% of published numbers (−14% to +29% across eight comparisons, including +2% against the PROLITH simulator), and its Monte Carlo agrees with a separate semi-analytical noise model within 1–3% for photon and acid noise. Failure rates are order-of-magnitude estimates, probably optimistic on missing holes. ASML's best 2018/2021 resists beat it by 46–110%, but those are non-CAR resists, so not like-for-like. See [`results/VERIFICATION.md`](results/VERIFICATION.md) and [`results/CALIBRATION.md`](results/CALIBRATION.md).

## Status

All Phase 2 work is complete.

- [x] Prototype simulator (v1) and failure-cliff sweep
- [x] Simulator v2 with literature parameters
- [x] Calibration against published data (`results/CALIBRATION.md`)
- [x] Verification, two passes plus a Nov 25 source re-check: sources, real-world data, another simulator, a semi-analytical check (`results/VERIFICATION.md`)
- [x] Unit tests: 16 (`sim/test_simulator_v2.py`, `sim/test_make_dataset.py`)
- [x] Dataset generator with built-in quality checks; pilot passed (`results/DATASET.md`)
- [x] Full dataset: 3,000 resists sampled by Latin hypercube, 2,520 accepted, 20,160 rows; QA passed with no warnings (`results/qa_full.md`)
- [x] Held-out test set: 400 resists from a different seed, 337 accepted, 2,696 rows
- [x] Model comparison: linear regression, random forest and neural network on three targets, with grouped cross-validation, held-out scores and a leakage check (`notebooks/ml_pipeline_ryan.ipynb`, `results/ml/`)
- [x] Written report
- [ ] Class presentation (Phase 3)

## Results

R² on the 337 held-out resists, which no model saw during training or selection (`results/ml/heldout_scores.csv`):

| Target | Linear | Random forest | Neural net |
| --- | --- | --- | --- |
| Failure rate (log10) | 0.503 | 0.974 | **0.989** |
| Mean hole size | 0.981 | 0.996 | **0.9997** |
| Size variation | 0.719 | 0.945 | **0.982** |

Held-out scores match the grouped cross-validation scores (`cv_scores.csv`), and an ordinary random split changes R² by at most 0.004 (`leakage_comparison.csv`).

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
notebooks/
  ml_pipeline_ryan.ipynb feature selection, grouped CV, held-out test, leakage check -> results/ml/
results/      figures, calibration and verification reports, QA reports, small CSV/JSON
  ml/                    cv_scores.csv, heldout_scores.csv, leakage_comparison.csv
data/
  pilot/                 50-resist pilot dataset (committed)
  full/, test/           full and held-out datasets (not committed; regenerate below)
```

## Run

```bash
pip install -r requirements.txt
python sim/test_simulator_v2.py   # ~35 s
python sim/simulator_v2.py        # one quick dose-to-size example
python sim/calibrate.py           # ~5 min
python sim/make_figures.py

# datasets (see results/DATASET.md); about 1 h 50 min for the full set on 8 cores
cd sim
python make_dataset.py --out ../data/full --n-resists 3000 --workers 8
python make_dataset.py --out ../data/test --n-resists 400 --workers 8 --seed 999
python check_dataset.py ../data/full

# models: open notebooks/ml_pipeline_ryan.ipynb and run all cells
# (reads ../data/full/rows.csv and ../data/test/rows.csv)
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
