"""Figures for the calibration report. Run after calibrate.py (reads results/calibration.json)."""
import os
RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simulator_v2 import Params, simulate, dose_to_size

r = json.load(open(RESULTS + "/calibration.json"))
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
BLUE, GREY, RED = "#2563eb", "#9ca3af", "#dc2626"

# ---- 1. hole-size variation: published vs simulated --------------------------------
rows = [
    ("Severi 2021\n24 nm, 80 mJ/cm²", r["T1_severi2021"], True),
    ("Naulleau model\n17 nm (dose assumed)", r["T4_naulleau_model"], True),
    ("Goethals 2013, 26 nm\n(total, incl. metrology)", r["T3a_goethals2013_c31p52"], False),
    ("Goethals 2013, 24 nm\n(no process window)", r["T3b_goethals2013_c29p48"], False),
    ("Wang 2019 HP24\n(interference litho)", r["T2a_wang2019_hp24"], False),
    ("Wang 2019 HP20\n(interference litho)", r["T2b_wang2019_hp20"], False),
]
fig, ax = plt.subplots(figsize=(8, 4.6))
for i, (lab, d, comparable) in enumerate(rows):
    pub, sim = d["published_lcdu"], d["lcdu_3s"]
    ax.fill_betweenx([i - 0.3, i + 0.3], pub * 0.8, pub * 1.2, color=GREY, alpha=0.25, lw=0)
    ax.plot(pub, i, "o", color="black", ms=7, label="published" if i == 0 else None)
    ax.plot(sim, i, "D", color=BLUE if comparable else GREY, ms=7,
            label="simulator v2" if i == 0 else None)
    ax.text(max(pub, sim) * 1.25 + 0.2, i, f"{(sim - pub) / pub:+.0%}", va="center",
            color=BLUE if comparable else GREY)
ax.set_yticks(range(len(rows)), [x[0] for x in rows])
ax.invert_yaxis()
ax.set_xlabel("Hole-to-hole size variation, LCDU 3σ (nm)")
ax.set_xlim(0, 7.5)
ax.set_title("Simulator v2 matches both like-for-like references within 8%", loc="left", fontsize=11)
ax.legend(loc="upper right", frameon=False)
fig.text(0.01, 0.01, "Grey band: ±20% of published. Grey diamonds: not like-for-like comparisons (see CALIBRATION.md).",
         fontsize=8, color="#555")
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(RESULTS + "/calibration_lcdu.png", dpi=150)

# ---- 2. CD vs dose (Yen 2018) ----------------------------------------------------------
y = r["T5_yen2018"]["sweep"]
fig, ax = plt.subplots(figsize=(5.2, 3.8))
ax.plot([30, 33, 36, 39], [15.0, 17.4, 19.5, 21.2], "o-", color="black", label="Yen 2018 (ASML simulation)")
ax.plot([s["dose"] for s in y], [s["mean_cd"] for s in y], "D--", color=BLUE, label="simulator v2")
ax.set_xlabel("Dose (mJ/cm²)"); ax.set_ylabel("Mean hole size (nm)")
ax.set_title("Printed size vs dose matches ASML's simulation within 5%", loc="left", fontsize=11)
ax.legend(frameon=False)
fig.tight_layout(); fig.savefig(RESULTS + "/calibration_cd_vs_dose.png", dpi=150)

# ---- 3. double cliff ----------------------------------------------------------------
c = r["T6_double_cliff"]["sweep"]
dose = np.array([s["dose"] for s in c])
n = c[0]["n_holes"]
floor = 1.0 / n
fig, ax = plt.subplots(figsize=(6, 3.9))
ax.semilogy(dose, np.maximum([s["p_missing"] for s in c], floor * 0.5), "o-", color=RED, label="missing holes")
ax.semilogy(dose, np.maximum([s["p_merged"] for s in c], floor * 0.5), "s-", color=BLUE, label="merged holes")
ax.axhline(floor, color=GREY, ls=":", lw=1)
ax.text(dose[-1], floor * 1.3, f"1 in {n:,} (measurement floor; points below = none observed)", ha="right", fontsize=8, color="#555")
ax.set_xlabel("Dose (mJ/cm²)"); ax.set_ylabel("Share of holes that fail")
ax.set_title("Dense 18 nm holes, 36 nm pitch: a failure cliff on each side", loc="left", fontsize=11)
ax.legend(frameon=False, loc="upper center")
fig.tight_layout(); fig.savefig(RESULTS + "/calibration_double_cliff.png", dpi=150)

# ---- 4. one tile at each stage ------------------------------------------------------
p = dose_to_size(Params(cd_target=20, pitch_x=40, pitch_y=40, dose=30))
f = simulate(p, 1, np.random.default_rng(4), return_fields=True)["fields"]
panels = [("Aerial image\n(share of open-frame dose)", f["image"], "magma"),
          ("Absorbed photons\n(per nm² column)", f["photons"], "magma"),
          ("Acid generated\n(per column)", f["acid"], "viridis"),
          ("Deprotected fraction\nafter bake", f["deprotected"], "viridis"),
          ("Printed holes\n(white = dissolved)", f["clear"], "gray")]
fig, axs = plt.subplots(1, 5, figsize=(15, 3.6))
for a, (t, img, cm) in zip(axs, panels):
    im = a.imshow(img, cmap=cm, interpolation="nearest")
    a.set_title(t, fontsize=10); a.set_xticks([]); a.set_yticks([])
fig.suptitle("Simulator v2: 3×3 array of 20 nm holes at 40 nm pitch, 30 mJ/cm² (one random trial)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.9]); fig.savefig(RESULTS + "/v2_stages.png", dpi=130)
print("figures written")
