"""Print pass 1 vs pass 2 LCDU for every comparison in calibration.json and extra_checks.json.

Pass 1 results are kept in results/pass1/.  Usage: python compare_passes.py
"""
import json
import os

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def load(name, sub=""):
    return json.load(open(os.path.join(R, sub, name)))


rows = []
for fname in ("calibration.json", "extra_checks.json"):
    a, b = load(fname, "pass1"), load(fname)
    for k, x in a.items():
        if isinstance(x, dict) and "lcdu_3s" in x:
            y = b[k]
            pub = x.get("published_lcdu")
            rows.append((k, pub, x["lcdu_3s"], y["lcdu_3s"], x["p_merged"], y["p_merged"]))

print("| Comparison | Published | Pass 1 | Pass 2 | Pass 1 vs published | Pass 2 vs published | Merged (pass 1 -> 2) |")
print("| --- | --- | --- | --- | --- | --- | --- |")
for k, pub, l1, l2, m1, m2 in rows:
    d1 = f"{(l1 / pub - 1) * 100:+.0f}%" if pub else ""
    d2 = f"{(l2 / pub - 1) * 100:+.0f}%" if pub else ""
    print(f"| {k} | {pub} | {l1:.2f} | {l2:.2f} | {d1} | {d2} | {m1:.4f} -> {m2:.4f} |")
