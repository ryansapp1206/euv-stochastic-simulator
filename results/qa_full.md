# Dataset QA: fullset

**Overall: PASS.** 20160 rows from 2520 accepted resists (3000 sampled). Simulator git 111048d, 200 tiles per row.

| Check | Result | Detail |
| --- | --- | --- |
| simulator version | PASS | dataset made with simulator sha256 52c642bf5b26, git 111048d; current file matches (recorded on a CRLF checkout, before the hash was newline-normalised) |
| run finished | PASS | 3000 of 3000 sampled resists processed |
| acceptance | PASS | 2520/3000 accepted (84%); rejected: {'rho_q_out_of_range': 431, 'cd<15.0': 49} |
| schema | PASS | missing columns: none |
| rows per resist | PASS | 20160 rows; every accepted resist has 8 rows: True; rows from rejected resists: 0; accepted resists without rows: 0 |
| duplicates | PASS | duplicate (resist, dose) rows: 0; reused random seeds: 0 |
| count bookkeeping | PASS | 0 rows where ok + missing + merged != holes; 0 rows where merge events don't fit merged holes (each event joins 2 or more holes) |
| missing values | PASS | NaN features: 0; NaN labels where holes survived: 0; rows with no surviving hole (labels NaN by design): 52 |
| inside the verified domain | PASS | out of range: none; quencher in (0.012, 0.09): True; dose factors in range: True |
| domain coverage | PASS | accepted resists span {'pitch_nm': '32-64', 'cd_ratio': '0.42-0.58', 'mask_bias_nm': '0.00128-4', 'd0_mjcm2': '20-89.9', 'sigma_e_nm': '1.2-2.5', 'sigma_a_nm': '4-6', 'qe': '2-3.5', 'rho_pag': '0.12-0.3', 'rho_pg': '1.5-3', 'thickness_nm': '30-50', 'absorptivity_um': '4.2-5.3', 'dev_threshold': '0.35-0.65'} |
| resist fit | PASS | |printed - target| at D0: median 0.028 nm, max 0.196 nm (limit 0.5) |
| more dose, bigger holes | PASS | 0 of 17473 dose steps where the mean size drops by more than 3 standard errors |
| more dose, fewer missing / more merged | PASS | missing-hole rate rising with dose: 0; merged rate falling with dose: 0 |
| nominal LCDU plausible | PASS | at D0: LCDU 1.75-6.27 nm, median 3.36 (published typical CAR 2-5 nm) |
| failures at nominal dose | PASS | 6 of 2520 resists lose >1% of holes even at D0 (no process window; kept, flagged) |
| LCDU follows image quality and dose | PASS | Spearman rank correlation of nominal LCDU with 1/(NILS x sqrt(photons)): 0.82 (p=0.0e+00) |
| rejected resists | INFO | 431 needed quencher 0.005-0.165 /nm^3 (published max 0.09); their D0 20-90 mJ/cm2 |
| cliff coverage | INFO | 0.60-0.75x: missing in 95%, merging in 0% of rows; 0.75-0.85x: missing in 64%, merging in 0% of rows; 0.85-0.95x: missing in 28%, merging in 1% of rows; 1.00-1.00x: missing in 3%, merging in 2% of rows; 1.05-1.20x: missing in 0%, merging in 5% of rows; 1.20-1.40x: missing in 0%, merging in 13% of rows; 1.40-1.65x: missing in 0%, merging in 27% of rows; 1.65-1.90x: missing in 0%, merging in 44% of rows |
| where merging happens | INFO | rows with merging: 28% at pitch <= 40 nm, 3% above 48 nm |
| label precision | INFO | mean size error 0.025 nm vs spread 8.0 nm; sigma error 1.7% vs spread 33% -> best honest R2 ~0.998 (0.99 on D0 rows only) |
| Monte Carlo vs semi-analytical model | PASS | ratio median 0.960, range 0.886-1.063; 0 of 2520 outside 0.80-1.10 (verification found 0.94-0.99 at nominal) |
| error bars match repeat runs | PASS | 722 repeated rows: spread of (difference / stated error) is 1.03 for mean size and 1.05 for sigma (1.0 = error bars exactly right); largest |z| 3.1 |
