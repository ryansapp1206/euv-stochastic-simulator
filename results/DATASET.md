# ML dataset: specification and final checks

*Nov 4, 2026. Generator: `sim/make_dataset.py`. Quality checks: `sim/check_dataset.py`. Pilot QA report: `data/pilot/qa_report.md`. The 400-row pilot is committed in `data/pilot/` so the ML code can be built on it; the full dataset is not committed. Every check below can be re-run from `sim/`.*

**Status: ready to generate.** A 60-sample pilot passed every quality check, with no warnings. The final checks found three problems, all fixed: merged holes come in pairs, so merges are now also counted as events; 32-bit row seeds could collide in a full run, so they are now 63-bit; and one bad sample could crash a long run, so errors are now logged as rejections.

## What one row is

Each sample is one **resist and layout**: a square array of contact holes at a given pitch, hole size and mask bias, with a resist (blur, acid yield, PAG, protecting groups, film, absorption, development threshold) and a sizing dose D0. Samples are drawn by Latin hypercube over the domain below. For each sample:

1. The resist's quencher loading is **fitted** so holes print at the target size at D0, as a resist formulator would.
2. The sample is **rejected and logged** if the fitted quencher falls outside the published range (0.012–0.09 /nm³), or if a fresh check run misses the target size by more than 0.5 nm.
3. Each accepted resist is exposed at **8 doses from 0.60× to 1.90× D0**: one jittered dose in each band (0.60–0.75, 0.75–0.85, 0.85–0.95, D0 exactly, 1.05–1.20, 1.20–1.40, 1.40–1.65, 1.65–1.90). That spans both failure cliffs.
4. Each row simulates 200 tiles of 3 × 3 holes, so **1,800 holes per row**. Each row has its own random stream derived from (seed, resist, dose). Any row regenerates exactly, and an interrupted run resumes where it stopped.

## Domain (only where the simulator was verified)

| Input | Range | Why this range |
| --- | --- | --- |
| Pitch | 32–64 nm, square arrays | ASML P32–P40, Goethals P52, PROLITH P50 checked; P30 was a poor match |
| Hole size / pitch | 0.42–0.58 (and hole ≥ 15 nm) | Size trend checked from 16 to 28 nm |
| Mask bias | 0–4 nm | Goethals used +5 nm, PROLITH +1 nm |
| Sizing dose D0 | 20–90 mJ/cm² (log-uniform) | Calibration and comparisons span 16–100 mJ/cm² |
| Electron blur σ | 1.2–2.5 nm | Sensitivity study range (Bhattarai 2017) |
| Acid diffusion σ | 4–6 nm | Sensitivity study range (Gronheid 2011) |
| Acids per photon | 2.0–3.5 | Sensitivity study range |
| PAG | 0.12–0.30 /nm³ | Sensitivity study range |
| Protecting groups | 1.5–3.0 /nm³ | Sensitivity study range |
| Film | 30–50 nm | Fallica 2015; Severi 50 nm |
| Absorption | 4.2–5.3 /µm | Fallica 2015 (CAR) |
| Development threshold | 0.35–0.65 | Sensitivity study tested 0.3–0.7 |
| Fixed | NA 0.33, quasar illumination, quencher blur 1.3 nm, K = 60 nm³, 1 nm pixels, one depth layer | Calibrated settings |

In the pilot, 50 of 60 samples were accepted. All 10 rejections needed more quencher than any published resist has (0.092–0.143 /nm³), and all had D0 above 50 mJ/cm². Very slow resists at high dose therefore appear less often than the Latin hypercube alone would give.

## Columns (`rows.csv`)

| Column | Unit | Meaning |
| --- | --- | --- |
| resist_id, dose_index | | Sample and dose step. **Use resist_id to split train and test.** |
| seed | | Random stream of this row (regenerates it exactly) |
| pitch_nm, cd_target_nm, cd_ratio, mask_cd_nm, mask_bias_nm | nm | Layout |
| dose_mjcm2, d0_mjcm2, dose_factor | mJ/cm² | Dose of this row, sizing dose, and their ratio |
| sigma_e_nm, sigma_a_nm, sigma_q_nm, qe, rho_pag, rho_q, rho_pg, k_amp, thickness_nm, absorptivity_um, dev_threshold, na | various | Resist and scanner (rho_q is the fitted quencher) |
| absorbed_frac, photons_abs_per_nm2 | –, /nm² | Share of light the film absorbs; absorbed photons per nm² at this dose |
| image_peak, nils | – | Peak image intensity (open frame = 1) and normalized image log-slope at the target edge: physics features for the linear model |
| n_holes, n_ok, n_missing, n_merged | count | Holes simulated, printed fine, missing (< 25% of target area), merged |
| n_merge_events, n_neighbor_pairs | count | Bridges (one joins 2 or more holes) and the neighbour pairs they could form (2 per hole) |
| p_missing, p_merged | – | n_missing / n_holes, n_merged / n_holes |
| mean_cd_nm, mean_cd_se | nm | Mean printed size of the holes that printed, and its standard error |
| sigma_cd_nm, sigma_cd_se, lcdu_3s_nm | nm | Hole-to-hole standard deviation, its standard error, and 3σ (LCDU) |

`resists.csv` has one row per sample, accepted or not, with the reason for each rejection, the fit check, the semi-analytical noise prediction and runtime. `replicates.csv` repeats three rows of every 10th resist with other seeds. `metadata.json` records the simulator's git commit and file hash, the domain and all settings.

## Final checks before generating

| Check | Result |
| --- | --- |
| Optics, computed a second way (exact Fourier series of the mask, no pixels or FFT; 8 layouts incl. NA 0.25, annular, hexapole, sub-pixel masks) | Images agree within 0.35% of peak intensity |
| Pixel size 0.5 vs 1 nm at four domain corners | The same resist prints 0.3–1.1 nm smaller at 0.5 nm. After re-fitting, as the dataset does for every resist, LCDU agrees within −4% to 0%, and both failure cliffs agree within noise (missing 0.31 vs 0.33, 0.10 vs 0.10, 0.06 vs 0.06; merged 0.209 vs 0.209) |
| Tile 3 × 3 vs 5 × 5 holes at the corners | Same size within 0.2 nm; failures within noise |
| Resist fit with 120 vs 400 tiles per step | Quencher within 0.7%, size within 0.05 nm |
| Rows regenerate from their stored seed | 80 of 80 exactly |
| Error bars, 80 rows re-run twice | Differences scatter by 0.98× (mean size), 0.96× (σ) and 0.92× (missing counts) of the stated errors; 1.0 is perfect |
| Are failure counts binomial? (40 repeat runs per setting) | Missing holes: yes (variance 1.10× and 1.30× binomial, consistent with 1; per tile 0.98–0.99×). **Merged holes: no, 1.84×**, because one bridge flags two holes. Merge events per neighbour pair: yes (0.91×). So the dataset now records merge events |
| Simulation vs the semi-analytical noise model, all 50 pilot resists | Ratio 0.91–1.01, median 0.96 (verification found 0.94–0.99) |
| Pilot QA (`check_dataset.py`), 400 rows | All PASS: bookkeeping, no duplicates or reused seeds, no missing values, everything inside the domain, fit error median 0.03 nm (max 0.11), mean size rises with dose in all 347 dose steps, missing falls and merging rises with dose, nominal LCDU 2.5–4.9 nm (median 3.4), no resist failing at D0, LCDU rank-correlates with 1/(NILS × √photons) at 0.83 |

**Cliff coverage in the pilot:** missing holes appear in 98%, 66% and 30% of rows in the three lowest dose bands. Merging appears in 6%, 16%, 28% and 44% of rows in the four highest bands, almost all at pitches of 40 nm or less (31% of those rows, against 2% above 48 nm). Within 1.9× dose, wider pitches simply don't merge.

## How precise the labels are

| Label | Typical error per row | Spread across the dataset | Best R² a model can honestly reach |
| --- | --- | --- | --- |
| Mean printed size | 0.026 nm | 7.6 nm (s.d.) | ≈ 1.000 |
| Hole-to-hole σ | 1.7% | 33% (coefficient of variation) | ≈ 0.997 (0.99 for D0 rows only) |
| Failure rates | Binomial over 1,800 holes (missing) or 3,600 pairs (merges) | 0 to ~1 | Use a binomial loss, not R² |

A model that beats those R² values is fitting simulation noise.

## How to use it for the models

- **Split by resist_id** (scikit-learn `GroupKFold` or `GroupShuffleSplit`). The 8 rows of one resist share everything but dose, so a row-level split leaks.
- **Features:** the layout, dose and resist columns, plus `photons_abs_per_nm2` and `nils` for the linear baseline. Never use `seed`, `resist_id` or any count, rate, mean or σ column as a feature.
- **Size targets:** `mean_cd_nm` and `sigma_cd_nm`, weighted by 1/SE². For σ, use rows with at least 100 printed holes, because near a cliff σ describes only the survivors.
- **Failure targets:** missing = binomial, `n_missing` out of `n_holes`. Merging = binomial, `n_merge_events` out of `n_neighbor_pairs`. Fit log-odds with a binomial loss (logistic regression on counts); don't take the log of zero rates.
- **Don't extrapolate** beyond the domain table, and say so when predicting rare-event rates below 1 in 1,800 per row.
- **Absolute doses** run from 12 to 171 mJ/cm² (0.6× the lowest D0 to 1.9× the highest), wider than the 16–100 mJ/cm² the comparisons covered. Rows outside that band are the same physics pushed further, and are almost all failure-cliff rows. Filter on `dose_mjcm2` if a model should stay inside the checked band.

## Known biases (from VERIFICATION.md)

- LCDU is within about ±30% of published typical CAR processes; the best modern resists are 46–110% better.
- Missing-hole rates near the small-size cliff are probably too low, by up to ~10× (two other simulators close holes more readily).
- Square arrays, NA 0.33 and quasar illumination only. No focus, mask 3D effects or flare. Simple threshold development. The 1 nm pixel is part of the model's definition, and the pixel check above shows it has no effect once each resist is fitted.

## Generating the full dataset

```bash
cd sim
python make_dataset.py --out ../data/full --n-resists 1500 --workers 8
python check_dataset.py ../data/full
python replicate_check.py ../data/full 80
```

- **Runtime:** about 34 s per sample on one core (42 s accepted, 12 s rejected). 1,500 samples give ~1,250 resists and ~10,000 rows: ~14 hours on one core, ~2 hours with 8 workers. With 1,000 samples (~6,600 rows) the time drops by a third.
- **Resuming:** re-run the same command. It refuses to resume if the simulator or the settings changed.
- Set `--workers` to the number of CPU cores, and allow about 0.5 GB of memory per worker (peak at 64 nm pitch). Tested here on Linux; it uses only standard multiprocessing, so it should also run on Windows, macOS and Colab. Colab's free tier has 2 cores, so 1,500 samples take ~7 hours; put `--out` on Google Drive so a run can resume across sessions.
