# Simulator v2 verification

*Nov 3, 2026. Every number below is reproducible: `python sim/verify.py`, `python sim/extra_checks.py`, `python sim/window_check.py 1`, then `python sim/make_verification_figure.py`.*

## Bottom line

For hole-to-hole size variation (LCDU), simulator v2 behaves like a **typical chemically amplified resist (CAR) on a 0.33 NA EUV scanner**:

- Against three measured references with known conditions, it lands at −17% to +8%.
- Against ASML industry data across pitches 40–32 nm, it lands at +14% to +27%.
- It correctly predicts that variation grows as pitch shrinks.

The **best modern resists** (ASML's 2021 data) do 30–50% better than the model, so the model is conservative for state-of-the-art processes.

**Failure rates** (missing and merged holes) are right in shape and roughly right in level: the two failure cliffs meet around one in 10⁵–10⁶ at 40 nm pitch, consistent with imec measurements. They are uncertain by about an order of magnitude, and one ASML simulation predicts about 10× more closing at small hole sizes.

**The project is safe to describe as:** "a literature-parameterized stochastic simulator that reproduces measured hole-size variation for typical CAR processes within about 25% (−17% to +27% across eight comparisons), with order-of-magnitude failure-rate estimates." It should not be described as matching any specific fab or resist.

![Parity plot](verification_parity.png)

## Corrections to the Nov 2 calibration

| Item | What the Nov 2 report said | What checking the source showed | Effect |
| --- | --- | --- | --- |
| Goethals 2013 | 31 nm and 29 nm holes | C31P52 / C29P48 are 31 / 29 nm **mask** holes printing **26 / 24 nm** | Re-run: 2.75 nm against 3.3 nm measured total (−17%) |
| Naulleau model cross-check | At ~82 mJ/cm² | The paper doesn't state the dose for that result; 82 comes from a different section | Comparison kept, labeled as an assumed dose |
| Yen 2018 (ASML) | Called "data" | ASML stochastic **simulations** (10,000 holes per dose); pitch, film and resist not stated | Relabeled as a simulation |
| Yen 2018 photon count | His holes absorb ~5× more photons | Misread. Columns are most likely cylinder diameters: 961 photons in 27 nm vs 420 in v2 | **About 2.2×**, not 5× |
| Wang 2019 (interference lithography) | Peak intensity ~4× the stated dose | That was my reasoning, not the paper's. PSI uses double-exposure four-beam interference, ~25 nm film, and a dose calibrated against a reference resist | Still not like-for-like, for the stated reasons. My attempt to simulate an interference image merged most holes, so it is not informative. |

## Source checks

| Claim used | Source | How verified | Status |
| --- | --- | --- | --- |
| 24 nm holes, staggered 80 × 46 nm, 80 mJ/cm², after develop, stochastic LCDU 2.6 nm | Severi et al. 2021 | Abstract + authors' PDF (via search) | Confirmed |
| CAR parameter set (QE 3, PAG 0.2, quencher 0.085, protecting groups 2 /nm³, 4.2 /µm, 35 nm, ranges 12/4/3 nm, threshold development) | Naulleau, Bhattarai, Neureuther | Full text read | Confirmed; "range" undefined in the paper |
| CD 15 / 17.4 / 19.5 / 21.2 nm at 30 / 33 / 36 / 39 mJ/cm²; ~10% closing at 15 nm | Yen 2018 slides | Full slides read | Confirmed; simulations |
| Dose-to-size 31.7 / 33.5 mJ/cm², 1σ LCDU 1.1 / 1.7 nm, NA 0.25, 60 nm film, quasar | Goethals et al. 2013 | Full slides read | Confirmed; corrected CDs |
| 0.33 NA, 40 nm CAR: LCDU 2.4 / 2.6 / 2.7 / 3.0 / 3.6 nm at P40–P30, 63–73 mJ/cm² | ASML 2024 investor presentation (SEC exhibit 99.4) | Full text read | Confirmed numbers; measured vs simulated and pattern type not stated |
| P40 CAR LCDU 2.7 / 2.2 / 2.1 / 2.0 nm at 66 / 65 / 61 / 53 mJ/cm² (2018–2021) | ASML, EUVL Workshop 2021 | Slides read | Confirmed numbers; conditions partly unknown |
| P40 dense holes at 0.33 NA need 79 mJ/cm² for <10⁻⁶ defect rate | ASML, EUVL Workshop 2021 | Slides read | Confirmed; simulated |
| 24 nm HP regular holes: LCDU 3.2 nm (small annular, 31 mJ/cm² resist) | van Ingen Schenau 2013 (ASML) | Slides via search | Confirmed numbers; sigma values assumed |
| Cliffs "touch at around 1 ppm" at 40 nm pitch | imec via SemiWiki 2019 | Article text | Confirmed (dose not stated) |
| Absorptivity 4.2–5.3 /µm (CAR), 15–19 /µm (metal oxide); films 25–45 nm | Fallica 2015 | Full poster read | Confirmed |
| ~2 acids per absorbed photon; electron blur 1.8–2 nm | Bhattarai 2017 | Full text read | Confirmed |
| Total blur ~4 nm at lowest bake; acid diffusion target ~5 nm | Gronheid 2011 | Abstract | Confirmed |
| PAG 0.04–0.2 /nm³; base 10–20% of PAG | AIP Advances 2025 | Abstract | Confirmed |
| Quencher/PAG 0.4 in a generic CAR | Contribution of EUV resist counting statistics (JMM 2021) | Abstract | Confirmed |
| 92 eV photons; 0.680 photons/nm² per mJ/cm² | Physics (E = hc/λ) | Computed; unit test | Confirmed |

## Real-world LCDU comparisons

| Reference | Type | Published | v2 | Difference |
| --- | --- | --- | --- | --- |
| Severi 2021, 24 nm holes, 80 mJ/cm² | Measured, stochastic only | 2.6 nm | 2.80 nm | +8% |
| van Ingen Schenau 2013, 24 nm HP, 31 mJ/cm² resist | Measured | 3.2 nm | 2.95 nm | −8% |
| Goethals 2013, 26 nm holes, NA 0.25 | Measured total | 3.3 nm | 2.75 nm | −17% (about +6% if ~80% of the total is stochastic, as in Severi) |
| ASML 2024, P40 / P36 / P34 / P32 | Industry, conditions partly unknown | 2.4 / 2.6 / 2.7 / 3.0 nm | 2.93 / 3.03 / 3.07 / 3.80 nm | +22% / +16% / +14% / +27% |
| ASML 2024, P30 | Industry, near resolution limit | 3.6 nm | 5.48 nm | +52%; v2's square-array setup is past its limit (26% merged) |
| ASML best CAR 2018 / 2021, P40 | Industry; hexagonal pattern approximated | 2.7 / 2.0 nm | 3.63 / 3.83 nm | +35% / +92% |
| Naulleau model, 16 nm HP | Published model | 3.45 nm | 3.55 nm | +3% |
| Goethals 2013, 24 nm (no process window) | Measured total | 5.1 nm | 2.76 nm | Not comparable |
| Wang 2019 (interference lithography) | Measured, different tool | 2.2–2.3 nm | 4.5 nm | Not comparable |

**Reading this table:** the model is most accurate for typical processes and is conservative (noisier) than the best recent resists and optimized illumination. The ASML "best resist" points also used an approximate hexagonal geometry with plain annular illumination, which handicaps the model.

## Numerical checks

| Check | Result |
| --- | --- |
| Pixel size 0.5 nm vs 1 nm | LCDU 3.76–3.89 vs 3.94 nm (1–5% lower); 1 nm is adequate |
| Light-source sampling 11 vs 7 points per pole | Identical, as expected: for an in-focus periodic pattern only which diffraction orders pass matters |
| Tile 5 × 5 vs 3 × 3 holes | 3.93 vs 3.94 nm |
| Repeat with three random seeds | 3.82 / 3.90 / 3.98 nm (±2%) |
| 8 unit tests | All pass |

## Sensitivity of the Severi match

Each uncertain parameter was moved across its published range, with the resist re-fitted every time. The match stays between **2.64 and 3.20 nm** against 2.6 published.

| Change | LCDU (nm) |
| --- | --- |
| Baseline | 2.82 |
| Acid blur 4.0 / 6.0 nm | 2.68 / 2.92 |
| Electron blur 1.2 / 2.5 nm | 2.74 / 2.88 |
| Acids per photon 2.0 / 3.5 | 2.86 / 2.76 |
| PAG 0.12 / 0.30 /nm³ | 3.20 / 2.64 |
| Protecting groups 1.5 / 3.0 /nm³ | 2.81 / 2.79 |
| Development threshold 0.3 / 0.7 | 2.85 / 2.78 |
| Absorptivity 5.3 /µm | 2.70 |
| Film 30 nm | 3.00 |
| Annular 0.3–0.7 | 2.72 |
| Five depth layers | 2.96 |

## Failure-rate checks

- **Dose windows at 40 nm pitch** (resist sized for 20 nm holes at 30, 50 and 80 mJ/cm²; 13,500 holes per point). Missing holes fall from ~1–5% at 75% of the sizing dose to ~10⁻³–10⁻⁴ at 85%. Merged holes appear at 120–160% of the sizing dose. Extending both cliffs, they meet around 10⁻⁵–10⁻⁶. That matches imec's measured "cliffs touch at around 1 ppm" at 40 nm pitch, and is consistent with ASML's simulation that a <10⁻⁶ window needs high dose (79 mJ/cm²). The windows widen as the sizing dose rises, as ASML reports.
- **At nominal dose,** 108,000 holes at 30 mJ/cm² gave zero failures, so the rate there is below 3 × 10⁻⁵ (rule of three).
- **Against Yen 2018 (ASML simulation),** v2 closes about 1% of holes at a 15 nm mean CD; ASML reports ~10%, with a resist absorbing ~2× more photons. Pitch and resist aren't stated, but this suggests v2 may under-predict closing at small sizes by up to ~10×.
- **Conclusion:** treat failure probabilities as order-of-magnitude estimates. Present the *shape* of the cliffs and how they move with dose, pitch and size as reliable trends, and absolute rates as uncertain within a factor of ~10.

## What would make it more realistic

Not needed for the project, but worth listing as future work:

- A full development-rate model and 3D dissolution paths (likely raises closing rates at small CD).
- Mask bias, attenuated phase-shift masks, and optimized pupils (hexapole, dipole), which would close part of the gap to the best modern processes.
- Calibration against a single resist's full dataset (dose, focus, CD, LCDU and defect counts together), which only a fab or resist supplier has.
