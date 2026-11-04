# Dataset QA: pilot

**Overall: PASS.** 400 rows from 50 accepted resists (60 sampled). Simulator git 0433188, 200 tiles per row.

| Check | Result | Detail |
| --- | --- | --- |
| simulator version | PASS | dataset made with simulator sha256 6991abacdfee, git 0433188; current file matches |
| run finished | PASS | 60 of 60 sampled resists processed |
| acceptance | PASS | 50/60 accepted (83%); rejected: {'rho_q_out_of_range': 10} |
| schema | PASS | missing columns: none |
| rows per resist | PASS | 400 rows; every accepted resist has 8 rows: True; rows from rejected resists: 0; accepted resists without rows: 0 |
| duplicates | PASS | duplicate (resist, dose) rows: 0; reused random seeds: 0 |
| count bookkeeping | PASS | 0 rows where ok + missing + merged != holes; 0 rows where merge events don't fit merged holes (each event joins 2 or more holes) |
| missing values | PASS | NaN features: 0; NaN labels where holes survived: 0; rows with no surviving hole (labels NaN by design): 0 |
| inside the verified domain | PASS | out of range: none; quencher in (0.012, 0.09): True; dose factors in range: True |
| domain coverage | PASS | accepted resists span {'pitch_nm': '32-63', 'cd_ratio': '0.422-0.579', 'mask_bias_nm': '0.124-3.84', 'd0_mjcm2': '20.5-89.5', 'sigma_e_nm': '1.23-2.48', 'sigma_a_nm': '4.01-6', 'qe': '2.02-3.48', 'rho_pag': '0.121-0.299', 'rho_pg': '1.5-2.99', 'thickness_nm': '30.1-49.8', 'absorptivity_um': '4.21-5.3', 'dev_threshold': '0.35-0.641'} |
| resist fit | PASS | |printed - target| at D0: median 0.029 nm, max 0.108 nm (limit 0.5) |
| more dose, bigger holes | PASS | 0 of 347 dose steps where the mean size drops by more than 3 standard errors |
| more dose, fewer missing / more merged | PASS | missing-hole rate rising with dose: 0; merged rate falling with dose: 0 |
| nominal LCDU plausible | PASS | at D0: LCDU 2.53-4.88 nm, median 3.38 (published typical CAR 2-5 nm) |
| failures at nominal dose | PASS | 0 of 50 resists lose >1% of holes even at D0 (no process window; kept, flagged) |
| LCDU follows image quality and dose | PASS | Spearman rank correlation of nominal LCDU with 1/(NILS x sqrt(photons)): 0.83 (p=1.1e-13) |
| rejected resists | INFO | 10 needed quencher 0.092-0.143 /nm^3 (published max 0.09); their D0 50-85 mJ/cm2 |
| cliff coverage | INFO | 0.60-0.74x: missing in 98%, merging in 0% of rows; 0.75-0.85x: missing in 66%, merging in 0% of rows; 0.85-0.95x: missing in 30%, merging in 0% of rows; 1.00-1.00x: missing in 2%, merging in 2% of rows; 1.05-1.20x: missing in 0%, merging in 6% of rows; 1.20-1.40x: missing in 0%, merging in 16% of rows; 1.40-1.65x: missing in 0%, merging in 28% of rows; 1.65-1.90x: missing in 0%, merging in 44% of rows |
| where merging happens | INFO | rows with merging: 31% at pitch <= 40 nm, 2% above 48 nm |
| label precision | INFO | mean size error 0.026 nm vs spread 7.6 nm; sigma error 1.7% vs spread 33% -> best honest R2 ~0.997 (0.99 on D0 rows only) |
| Monte Carlo vs semi-analytical model | PASS | ratio median 0.960, range 0.910-1.006; 0 of 50 outside 0.80-1.10 (verification found 0.94-0.99 at nominal) |
| error bars match repeat runs | PASS | 15 repeated rows: spread of (difference / stated error) is 1.10 for mean size and 0.71 for sigma (1.0 = error bars exactly right); largest |z| 2.3 |
