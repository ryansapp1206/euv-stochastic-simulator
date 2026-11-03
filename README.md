# euv-stochastic-simulator

CS 427 Project 3 (Hood College) — Ryan Sapp

A Monte Carlo simulator of how contact holes print under EUV lithography, plus machine-learning surrogate models (linear regression, random forest, neural network) that predict failure rate and hole-size variation from process settings.

**Simulator v2** follows the structure of published stochastic resist models: partially coherent NA 0.33 optics, photon shot noise, PAG, quencher and protecting-group counting, electron blur, acid and quencher diffusion, and threshold development. Its parameters come from published values. On both like-for-like references, hole-size variation lands within 8% of the published numbers. See [`results/CALIBRATION.md`](results/CALIBRATION.md).

## Status

- [x] Prototype simulator (v1) and failure-cliff sweep
- [x] Simulator v2 with literature parameters
- [x] Calibration against published data (`results/CALIBRATION.md`)
- [x] Unit tests (`sim/test_simulator_v2.py`)
- [ ] Dataset generation (Latin hypercube over process settings)
- [ ] Model comparison: linear regression, random forest, neural network
- [ ] Report and presentation

## Layout

```
sim/
  simulator_v2.py        v2 model (use this)
  calibrate.py           runs all calibration targets -> results/calibration.json
  make_figures.py        calibration figures -> results/*.png
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
```

Run commands from the repository root.

## Quick use

```python
import numpy as np
from simulator_v2 import Params, simulate, summarize, dose_to_size

p = Params(cd_target=20, pitch_x=40, pitch_y=40, dose=30)   # 3x3 hole array per trial
p = dose_to_size(p)                  # fit quencher so holes print at 20 nm at this dose
s = summarize(simulate(p, 400, np.random.default_rng(0)))
print(s)   # mean_cd, sigma_cd, lcdu_3s, p_missing, p_merged
```
