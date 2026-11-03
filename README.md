# euv-stochastic-simulator

CS 427 Project 3 (Hood College) — Ryan Sapp

A simplified Monte Carlo simulator of how contact holes print under EUV lithography, plus machine-learning surrogate models (linear regression, random forest, neural network) that predict failure rate and hole-size variation from process settings.

The simulator is a teaching-grade model: it is calibrated against published measurements, not against any specific fab.

## Status

- [x] Prototype simulator (photon shot noise, acid generation, quencher, blur, threshold)
- [x] Dose sweep showing the failure cliff (`results/sweep.csv`)
- [x] Calibration gap identified: prototype hole-size variation is 2–4x lower than published values
- [ ] Physics upgrades 1–6 (see project plan)
- [ ] Calibration checkpoint (Nov 14)
- [ ] Dataset generation, model comparison, report

## Layout

```
sim/          simulator and analysis scripts
  simulator.py   core Monte Carlo model
  sweep.py       failure rate vs dose for several hole sizes
  lcdu_check.py  hole-size variation (3 sigma) vs published values
  stages.py      figure of one hole at each simulation stage
results/      generated figures and small CSVs
notebooks/    exploration, calibration, modeling
report/       written report
data/         large generated datasets (not committed)
```

## Run

```bash
pip install -r requirements.txt
python sim/simulator.py     # quick sanity check
python sim/sweep.py         # writes results/sweep.csv (~45 s)
python sim/lcdu_check.py    # prints hole-size variation at published settings
python sim/stages.py        # writes results/hole_stages.png (~1 min)
```

Run commands from the repository root.
