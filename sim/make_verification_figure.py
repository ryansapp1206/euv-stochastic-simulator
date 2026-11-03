"""Parity plot: published vs simulated hole-size variation (LCDU, 3 sigma) for every
real-world comparison. Reads results/calibration.json, verification.json, extra_checks.json."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

cal = json.load(open("results/calibration.json"))
ver = json.load(open("results/verification.json"))
ext = json.load(open("results/extra_checks.json"))

BLUE, AMBER, GREY = "#2563eb", "#d97706", "#6b7280"
pts = [
    # (label, published, simulated, group)
    ("Severi 2021 (imec)", 2.6, cal["T1_severi2021"]["lcdu_3s"], "measured"),
    ("van Ingen Schenau 2013", 3.2, ext["I_vanIngenSchenau2013_24hp_small_annular"]["lcdu_3s"], "measured"),
    ("Goethals 2013, 26 nm*", 3.3, cal["T3a_goethals2013_c31p52"]["lcdu_3s"], "measured"),
    ("ASML P40", 2.4, ext["F_asml_P40"]["lcdu_3s"], "industry"),
    ("ASML P36", 2.6, ext["F_asml_P36"]["lcdu_3s"], "industry"),
    ("ASML P34", 2.7, ext["F_asml_P34"]["lcdu_3s"], "industry"),
    ("ASML P32", 3.0, ext["F_asml_P32"]["lcdu_3s"], "industry"),
    ("ASML best resist 2018", 2.7, ext["G_asml_hexP40_2018"]["lcdu_3s"], "industry"),
    ("ASML best resist 2021", 2.0, ext["G_asml_hexP40_2021"]["lcdu_3s"], "industry"),
    ("Naulleau model", 3.45, cal["T4_naulleau_model"]["lcdu_3s"], "model"),
]
style = {"measured": (BLUE, "o", "Measured, conditions known"),
         "industry": (AMBER, "s", "Industry data, conditions partly unknown"),
         "model": (GREY, "^", "Published model (cross-check)")}

fig, ax = plt.subplots(figsize=(7.2, 6.2))
x = np.linspace(1.5, 4.5, 50)
ax.fill_between(x, 0.7 * x, 1.3 * x, color="#e5e7eb", alpha=0.6, lw=0, label="±30%")
ax.fill_between(x, 0.8 * x, 1.2 * x, color="#d1d5db", alpha=0.8, lw=0, label="±20%")
ax.plot(x, x, color="black", lw=0.8)
seen = set()
offsets = {"ASML P36": (6, -10), "ASML P34": (6, 2), "ASML P40": (-58, 6), "Severi 2021 (imec)": (6, -12),
           "Naulleau model": (6, -2), "van Ingen Schenau 2013": (6, -4), "Goethals 2013, 26 nm*": (6, -8),
           "ASML P32": (6, -2), "ASML best resist 2018": (6, -2), "ASML best resist 2021": (6, -2)}
for lab, pub, sim, g in pts:
    c, m, name = style[g]
    ax.plot(pub, sim, m, color=c, ms=8, label=None if g in seen else name)
    seen.add(g)
    ax.annotate(f"{lab} ({(sim - pub) / pub:+.0%})", (pub, sim), xytext=offsets.get(lab, (6, 0)),
                textcoords="offset points", fontsize=8.5, color="#111")
ax.set_xlim(1.5, 4.2); ax.set_ylim(1.5, 4.2)
ax.set_xlabel("Published LCDU, 3σ (nm)")
ax.set_ylabel("Simulator v2 LCDU, 3σ (nm)")
ax.set_title("Typical processes land within −17% to +27%; the best modern resists beat the model",
             loc="left", fontsize=10.5)
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.text(0.01, 0.01, "*Published total includes metrology and mask error. Points on the diagonal match exactly.",
         fontsize=7.5, color="#555")
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig("results/verification_parity.png", dpi=150)
print("written")
