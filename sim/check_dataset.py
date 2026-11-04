"""Quality checks for a dataset made by make_dataset.py.

Every check prints PASS / WARN / FAIL with its numbers, and the full report is
written to <dataset>/qa_report.md and qa_report.json. A dataset is fit to use
when nothing FAILs; each WARN says what to keep in mind.

Usage: python check_dataset.py ../data/pilot
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_dataset import (DOMAIN, MIN_CD_NM, RHO_Q_RANGE, FIT_TOL_NM, DOSE_STRATA, ROW_FIELDS,  # noqa: E402
                          RESIST_FIELDS, file_sha256)

results = []


def check(name, status, detail):
    results.append(dict(check=name, status=status, detail=detail))
    print(f"[{status}] {name}: {detail}", flush=True)


def main(path):
    meta = json.load(open(os.path.join(path, "metadata.json")))
    R = pd.read_csv(os.path.join(path, "resists.csv"))
    D = pd.read_csv(os.path.join(path, "rows.csv"))
    rep_path = os.path.join(path, "replicates.csv")
    P = pd.read_csv(rep_path) if os.path.exists(rep_path) and os.path.getsize(rep_path) > 0 else pd.DataFrame()

    # ---- provenance
    sim_path = os.path.join(HERE, "simulator_v2.py")
    sha = file_sha256(sim_path)                                  # newline-normalised (LF)
    body = open(sim_path, "rb").read().replace(b"\r\n", b"\n")
    sha_crlf = hashlib.sha256(body.replace(b"\n", b"\r\n")).hexdigest()
    rec = meta["simulator_sha256"]
    # Datasets made before file_sha256 normalised newlines recorded the raw digest of
    # whatever git checked out, so a run on Windows recorded the CRLF digest. Accept it.
    same = rec in (sha, sha_crlf)
    note = "matches" if rec == sha else \
        "matches (recorded on a CRLF checkout, before the hash was newline-normalised)" if same else "DIFFERS"
    check("simulator version", "PASS" if same else "WARN",
          f"dataset made with simulator sha256 {rec[:12]}, git {meta['git_commit'][:7]}; current file {note}")
    check("run finished", "PASS" if "finished" in meta else "WARN",
          f"{len(R)} of {meta['n_resists']} sampled resists processed")

    # ---- accounting
    acc = R[R.accepted.astype(str) == "True"]
    reasons = R[R.accepted.astype(str) != "True"].reason.value_counts().to_dict()
    check("acceptance", "PASS" if len(acc) / max(len(R), 1) >= 0.5 else "WARN",
          f"{len(acc)}/{len(R)} accepted ({100 * len(acc) / max(len(R), 1):.0f}%); rejected: {reasons}")
    missing_cols = [c for c in ROW_FIELDS if c not in D.columns] + [c for c in RESIST_FIELDS if c not in R.columns]
    check("schema", "PASS" if not missing_cols else "FAIL", f"missing columns: {missing_cols or 'none'}")
    per = D.groupby("resist_id").size()
    n_dose = len(DOSE_STRATA) - 1
    bad_count = per[per != n_dose]
    orphans = set(D.resist_id) - set(acc.resist_id)
    lost = set(acc.resist_id) - set(D.resist_id)
    check("rows per resist", "PASS" if bad_count.empty and not orphans and not lost else "FAIL",
          f"{len(D)} rows; every accepted resist has {n_dose} rows: {bad_count.empty}; rows from rejected "
          f"resists: {len(orphans)}; accepted resists without rows: {len(lost)}")
    dup = D.duplicated(["resist_id", "dose_index"]).sum()
    dup_seed = D.seed.duplicated().sum()
    check("duplicates", "PASS" if dup == 0 and dup_seed == 0 else "FAIL",
          f"duplicate (resist, dose) rows: {dup}; reused random seeds: {dup_seed}")
    tot = (D.n_ok + D.n_missing + D.n_merged != D.n_holes).sum() + (np.abs(D.p_missing - D.n_missing / D.n_holes)
                                                                     > 1e-12).sum()
    ev_bad = int(((D.n_merge_events > D.n_merged / 2 + 1e-9) | ((D.n_merged > 0) != (D.n_merge_events > 0))).sum())
    check("count bookkeeping", "PASS" if tot == 0 and ev_bad == 0 else "FAIL",
          f"{tot} rows where ok + missing + merged != holes; {ev_bad} rows where merge events don't fit merged holes "
          f"(each event joins 2 or more holes)")

    # ---- missing values
    feats = ["pitch_nm", "cd_target_nm", "mask_cd_nm", "dose_mjcm2", "d0_mjcm2", "dose_factor", "sigma_e_nm",
             "sigma_a_nm", "qe", "rho_pag", "rho_q", "rho_pg", "thickness_nm", "absorptivity_um", "dev_threshold",
             "photons_abs_per_nm2", "image_peak", "nils"]
    nan_f = int(D[feats].isna().sum().sum())
    lab_nan = D[["mean_cd_nm", "sigma_cd_nm"]].isna().any(axis=1)
    bad_lab = int((lab_nan & (D.n_ok > 1)).sum())
    check("missing values", "PASS" if nan_f == 0 and bad_lab == 0 else "FAIL",
          f"NaN features: {nan_f}; NaN labels where holes survived: {bad_lab}; rows with no surviving hole "
          f"(labels NaN by design): {int((D.n_ok <= 1).sum())}")

    # ---- domain
    out = []
    for k, (lo, hi, _) in DOMAIN.items():
        col = acc[k]
        if k == "pitch_nm":
            lo, hi = round(lo), round(hi)
        if (col < lo - 1e-9).any() or (col > hi + 1e-9).any():
            out.append(k)
    q_ok = acc.rho_q.between(*RHO_Q_RANGE).all() and (acc.cd_target_nm >= MIN_CD_NM - 1e-9).all()
    f_ok = D.dose_factor.between(DOSE_STRATA[0], DOSE_STRATA[-1]).all()
    check("inside the verified domain", "PASS" if not out and q_ok and f_ok else "FAIL",
          f"out of range: {out or 'none'}; quencher in {RHO_Q_RANGE}: {q_ok}; dose factors in range: {f_ok}")
    cover = {k: f"{acc[k].min():.3g}-{acc[k].max():.3g}" for k in DOMAIN}
    check("domain coverage", "PASS", f"accepted resists span {cover}")

    # ---- fit
    fe = acc.fit_error_nm.abs()
    check("resist fit", "PASS" if fe.max() <= FIT_TOL_NM else "FAIL",
          f"|printed - target| at D0: median {fe.median():.3f} nm, max {fe.max():.3f} nm (limit {FIT_TOL_NM})")

    # ---- physics trends within each resist
    inc_viol = miss_viol = merg_viol = pairs = 0
    for _, g in D.sort_values("dose_factor").groupby("resist_id"):
        g = g.reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            if a.n_ok >= 50 and b.n_ok >= 50:
                pairs += 1
                z = (b.mean_cd_nm - a.mean_cd_nm) / np.hypot(a.mean_cd_se, b.mean_cd_se)
                inc_viol += z < -3
            for key in ("n_missing", "n_merged"):
                pa, pb = a[key] / a.n_holes, b[key] / b.n_holes
                se = np.sqrt(max(pa * (1 - pa) / a.n_holes + pb * (1 - pb) / b.n_holes, 1e-12))
                z = (pb - pa) / se
                if key == "n_missing" and z > 3 and b[key] >= 5:
                    miss_viol += 1
                if key == "n_merged" and z < -3 and a[key] >= 5:
                    merg_viol += 1
    check("more dose, bigger holes", "PASS" if inc_viol == 0 else "FAIL",
          f"{inc_viol} of {pairs} dose steps where the mean size drops by more than 3 standard errors")
    check("more dose, fewer missing / more merged", "PASS" if miss_viol + merg_viol == 0 else "FAIL",
          f"missing-hole rate rising with dose: {miss_viol}; merged rate falling with dose: {merg_viol}")

    nom = D[np.isclose(D.dose_factor, 1.0)]
    lc = nom.lcdu_3s_nm
    check("nominal LCDU plausible", "PASS" if lc.between(1.5, 8).all() else "WARN",
          f"at D0: LCDU {lc.min():.2f}-{lc.max():.2f} nm, median {lc.median():.2f} "
          f"(published typical CAR 2-5 nm)")
    nf = nom[(nom.p_missing > 0.01) | (nom.p_merged > 0.01)]
    check("failures at nominal dose", "PASS" if len(nf) <= 0.05 * len(nom) else "WARN",
          f"{len(nf)} of {len(nom)} resists lose >1% of holes even at D0 (no process window; kept, flagged)")

    # literature trend: LCDU ~ 1 / (NILS * sqrt(absorbed photons)) (ILS correlator, shot-noise scaling)
    x = 1 / (nom.nils * np.sqrt(nom.photons_abs_per_nm2))
    rho, pv = spearmanr(x, lc)
    check("LCDU follows image quality and dose", "PASS" if rho > 0.5 else "WARN",
          f"Spearman rank correlation of nominal LCDU with 1/(NILS x sqrt(photons)): {rho:.2f} (p={pv:.1e})")

    # ---- informational: rejections, cliff coverage, label precision (do not affect the verdict)
    rj = R[R.reason.astype(str) == "rho_q_out_of_range"]
    if len(rj):
        check("rejected resists", "INFO",
              f"{len(rj)} needed quencher {rj.rho_q.min():.3f}-{rj.rho_q.max():.3f} /nm^3 (published max "
              f"{RHO_Q_RANGE[1]}); their D0 {rj.d0_mjcm2.min():.0f}-{rj.d0_mjcm2.max():.0f} mJ/cm2")
    cov = D.groupby("dose_index").agg(lo=("dose_factor", "min"), hi=("dose_factor", "max"),
                                      miss=("n_missing", lambda x: (x > 0).mean()),
                                      merg=("n_merge_events", lambda x: (x > 0).mean()))
    check("cliff coverage", "INFO", "; ".join(f"{r.lo:.2f}-{r.hi:.2f}x: missing in {100 * r.miss:.0f}%, merging in "
                                              f"{100 * r.merg:.0f}% of rows" for r in cov.itertuples()))
    small, wide = D[D.pitch_nm <= 40], D[D.pitch_nm > 48]
    check("where merging happens", "INFO",
          f"rows with merging: {100 * (small.n_merge_events > 0).mean():.0f}% at pitch <= 40 nm, "
          f"{100 * (wide.n_merge_events > 0).mean():.0f}% above 48 nm")
    ok = D[D.n_ok >= 50]
    rs, cv = (ok.sigma_cd_se / ok.sigma_cd_nm).median(), ok.sigma_cd_nm.std() / ok.sigma_cd_nm.mean()
    rs0 = (nom.sigma_cd_se / nom.sigma_cd_nm).median()
    cv0 = nom.sigma_cd_nm.std() / nom.sigma_cd_nm.mean()
    check("label precision", "INFO",
          f"mean size error {ok.mean_cd_se.median():.3f} nm vs spread {ok.mean_cd_nm.std():.1f} nm; sigma error "
          f"{100 * rs:.1f}% vs spread {100 * cv:.0f}% -> best honest R2 ~{1 - rs ** 2 / cv ** 2:.3f} "
          f"({1 - rs0 ** 2 / cv0 ** 2:.2f} on D0 rows only)")

    # ---- semi-analytical cross-check
    if "mc_over_analytic" in acc and acc.mc_over_analytic.notna().any():
        r = acc.mc_over_analytic.dropna()
        n_out = int((~r.between(0.80, 1.10)).sum())
        check("Monte Carlo vs semi-analytical model", "PASS" if n_out <= 0.05 * len(r) else "WARN",
              f"ratio median {r.median():.3f}, range {r.min():.3f}-{r.max():.3f}; {n_out} of {len(r)} outside "
              f"0.80-1.10 (verification found 0.94-0.99 at nominal)")

    # ---- honest error bars: replicate rows with fresh seeds
    if not P.empty:
        M = P.merge(D, on=["resist_id", "dose_index"], suffixes=("_r", ""))
        M = M[(M.n_ok >= 50) & (M.n_ok_r >= 50)]
        zm = (M.mean_cd_nm_r - M.mean_cd_nm) / np.hypot(M.mean_cd_se_r, M.mean_cd_se)
        zs = (M.sigma_cd_nm_r - M.sigma_cd_nm) / np.hypot(M.sigma_cd_se_r, M.sigma_cd_se)
        ok = 0.7 <= zm.std() <= 1.4 and 0.7 <= zs.std() <= 1.4 and (np.abs(np.r_[zm, zs]) < 4.5).all()
        check("error bars match repeat runs", "PASS" if ok else "WARN",
              f"{len(M)} repeated rows: spread of (difference / stated error) is {zm.std():.2f} for mean size "
              f"and {zs.std():.2f} for sigma (1.0 = error bars exactly right); largest |z| "
              f"{np.abs(np.r_[zm, zs]).max():.1f}")

    # ---- report
    worst = "FAIL" if any(r["status"] == "FAIL" for r in results) else \
        "WARN" if any(r["status"] == "WARN" for r in results) else "PASS"
    lines = [f"# Dataset QA: {os.path.basename(os.path.abspath(path))}", "",
             f"**Overall: {worst}.** {len(D)} rows from {len(acc)} accepted resists "
             f"({len(R)} sampled). Simulator git {meta['git_commit'][:7]}, {meta['trials']} tiles per row.", "",
             "| Check | Result | Detail |", "| --- | --- | --- |"]
    lines += [f"| {r['check']} | {r['status']} | {r['detail']} |" for r in results]
    open(os.path.join(path, "qa_report.md"), "w").write("\n".join(lines) + "\n")
    json.dump(dict(overall=worst, checks=results), open(os.path.join(path, "qa_report.json"), "w"), indent=1,
              default=str)
    print("overall:", worst)
    return worst


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "data", "pilot"))
