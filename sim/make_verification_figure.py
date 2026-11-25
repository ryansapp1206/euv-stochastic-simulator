"""Parity plot: published vs simulated hole-size variation (LCDU, 3 sigma) for every
real-world comparison, after the pass 2 fixes.  Reads results/calibration.json,
extra_checks.json and results/pass2/*.json.  Writes results/verification_parity.png."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
cal = json.load(open(RESULTS + "/calibration.json"))
ext = json.load(open(RESULTS + "/extra_checks.json"))
p2 = {}
for f in ("P1", "P1b", "P2", "P3", "P4"):
    p2.update(json.load(open(f"{RESULTS}/pass2/{f}.json")))

# reference palette, first three categorical slots (validated all-pairs, light mode)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, SURFACE = "#0b0b0b", "#52514e", "#fcfcfb"

mean = lambda m: sum(m) / len(m)
pts = [  # (label, published, simulated, group)
    ("van Ingen Schenau 2013", 3.2, p2["P4_vanIngenSchenau2013"]["lcdu_3s"], "measured"),
    ("Goethals 2013, 26 nm*", 3.3, p2["P4_goethals2013_C31P52"]["lcdu_3s"], "measured"),
    ("ASML P40", 2.4, p2["P4_asml_P40"]["lcdu_3s"], "industry"),
    ("ASML P36", 2.6, p2["P4_asml_P36"]["lcdu_3s"], "industry"),
    ("ASML P34", 2.7, p2["P4_asml_P34"]["lcdu_3s"], "industry"),
    ("ASML P32", 3.0, p2["P4_asml_P32"]["lcdu_3s"], "industry"),
    # ASML/imec P40 hexagonal pillars/CH: the source slides label these resists non-CAR,
    # so they are not like-for-like with v2's chemically amplified resist model
    ("ASML non-CAR 2018", 2.7, p2["P3_asml_hexP40_2018_hexapole"]["lcdu_3s"], "industry"),
    ("ASML non-CAR 2021", 2.0, p2["P3_asml_hexP40_2021_hexapole"]["lcdu_3s"], "industry"),
    ("Naulleau model", 3.45, p2["P4_naulleau_model"]["lcdu_3s"], "model"),
    ("PROLITH 2017 study", 3.12, next(r for r in p2["P2_mack2017_prolith"]["sweep"] if r["mask_cd"] == 24.0)["lcdu_3s"],
     "model"),
]
severi = [(tag, mean(p2[f"P1_severi_{tag}_hexapole"]["measured"]),
           p2[f"P1b_severi_{tag}_rows_hexapole"]["lcdu_3s"], p2[f"P1_severi_{tag}_hexapole"]["lcdu_3s"])
          for tag in ("H80V46", "H70V40", "H66V38")]

style = {"measured": (BLUE, "o", "Measured, conditions known"),
         "industry": (ORANGE, "s", "Industry data, conditions partly unknown"),
         "model": (AQUA, "^", "Other models (cross-check)")}

fig, ax = plt.subplots(figsize=(7.6, 6.6), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
x = np.linspace(1.5, 6, 50)
ax.fill_between(x, 0.7 * x, 1.3 * x, color="#ecebe8", lw=0, label="±30%")
ax.fill_between(x, 0.8 * x, 1.2 * x, color="#dcdbd7", lw=0, label="±20%")
ax.plot(x, x, color=INK2, lw=0.8)

seen = set()
offsets = {"ASML P36": (7, -9), "ASML P34": (7, 3), "ASML P40": (-60, 4), "Naulleau model": (7, -3),
           "van Ingen Schenau 2013": (7, -3), "Goethals 2013, 26 nm*": (7, -9), "ASML P32": (7, -3),
           "ASML non-CAR 2018": (7, -3), "ASML non-CAR 2021": (7, -3), "PROLITH 2017 study": (7, 2)}
for lab, pub, sim, g in pts:
    c, m, name = style[g]
    ax.plot(pub, sim, m, color=c, ms=8, mec=SURFACE, mew=1.5, label=None if g in seen else name, zorder=3)
    seen.add(g)
    right = lab in ("ASML non-CAR 2018", "ASML P40")
    ax.annotate(f"{lab} ({(sim - pub) / pub:+.0%})", (pub, sim),
                xytext=(-9, -3) if lab == "ASML non-CAR 2018" else (-12, 7) if lab == "ASML P40"
                else offsets.get(lab, (7, 0)),
                textcoords="offset points", fontsize=8, color=INK, ha="right" if right else "left")

# Severi 2021: two readings of the staggered layout, joined per pitch
for i, (tag, pub, rows, hexa) in enumerate(severi):
    ax.plot([pub, pub], [rows, hexa], color=BLUE, lw=1.2, zorder=2)
    ax.plot(pub, rows, "D", color=BLUE, ms=7, mec=SURFACE, mew=1.5, zorder=3,
            label="Severi 2021, rows reading (pass 1)" if i == 0 else None)
    ax.plot(pub, hexa, "D", mfc=SURFACE, mec=BLUE, mew=1.8, ms=7, zorder=3,
            label="Severi 2021, hexagonal reading" if i == 0 else None)
    if tag == "H70V40":
        ax.annotate(f"Severi {tag}", (pub, rows), xytext=(-8, -4), textcoords="offset points",
                    fontsize=8, color=INK, ha="right")
    else:
        ax.annotate(f"Severi {tag}", (pub, hexa), xytext=(8, -3), textcoords="offset points", fontsize=8, color=INK)

ax.set_xlim(1.5, 4.0); ax.set_ylim(1.5, 5.6)
ax.set_xlabel("Published LCDU, 3σ (nm)", color=INK)
ax.set_ylabel("Simulator v2 LCDU, 3σ (nm)", color=INK)
# headline range is computed from the eight typical-process points, not typed in
typical = [(s - p) / p for lab, p, s, g in pts if "non-CAR" not in lab]
ax.set_title(f"Typical processes: {min(typical):+.0%} to {max(typical):+.0%}. ".replace("-", "−") + "ASML's non-CAR resists and the tightest\n"
             "hexagonal reading do not match, and the Severi layout is ambiguous.",
             loc="left", fontsize=10.5, color=INK)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.grid(color="#e7e6e2", lw=0.6); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(INK2)
ax.tick_params(colors=INK2)
sev_note = ", ".join(f"{t} {(r - p) / p:+.0%} / {(h - p) / p:+.0%}" for t, p, r, h in severi)
fig.text(0.01, 0.035, "Severi 2021, rows reading / hexagonal reading: " + sev_note + ".", fontsize=7.2, color=INK2)
fig.text(0.01, 0.01, "*Published total includes metrology and mask error. Pass 2: fresh seeds, 13,500+ holes per point. "
         "ASML P30 (+79%) is off this chart.", fontsize=7.2, color=INK2)
fig.tight_layout(rect=[0, 0.055, 1, 1])
fig.savefig(RESULTS + "/verification_parity.png", dpi=150, facecolor=SURFACE)
print("written")
