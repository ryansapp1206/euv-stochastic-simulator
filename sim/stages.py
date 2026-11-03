import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from simulator import aerial, PHOTONS_PER_NM2_PER_MJCM2, RR, G, calibrate_threshold

d, dose, absorption, yld, q, blur = 20, 18, 0.2, 2.0, 1.0, 4.0
thr = calibrate_threshold()
I = aerial(d, 4.0)
lam = dose * PHOTONS_PER_NM2_PER_MJCM2 * absorption * I

def trial(rng):
    ph = rng.poisson(lam)
    ac = rng.poisson(yld * ph)
    qu = rng.poisson(q, size=(G, G))
    net = gaussian_filter((ac - qu).astype(float), blur)
    hole = net > thr
    area = hole[RR <= d].sum()
    return ph, net, hole, area

rng = np.random.default_rng(11)
ok = bad = None
while ok is None or bad is None:
    r = trial(rng)
    if r[3] <= 5 and bad is None: bad = r
    if r[3] >= 0.6 * np.pi * (d / 2) ** 2 and ok is None: ok = r

c = slice(12, 52)
mask = (RR <= d / 2)
fig, axes = plt.subplots(2, 5, figsize=(15, 7.4))
vmax_ph = max(ok[0][c, c].max(), bad[0][c, c].max())
vmin_net = min(ok[1][c, c].min(), bad[1][c, c].min()); vmax_net = max(ok[1][c, c].max(), bad[1][c, c].max())
cols = ['1. Mask pattern\n(20 nm hole)', '2. Light reaching resist\n(blurred by optics)', '3. Absorbed photons\n(random counts per nm²)', '4. Acid after bake\n(blurred, minus quencher)', '5. Printed result\n(dissolves above threshold)']
for row, (ph, net, hole, area) in enumerate([ok, bad]):
    failed = row == 1
    inside = ph[RR <= d / 2].sum()
    panels = [mask[c, c], I[c, c], ph[c, c], net[c, c], hole[c, c]]
    cmaps = ['gray', 'magma', 'magma', 'viridis', 'gray']
    for j, (p, cm) in enumerate(zip(panels, cmaps)):
        ax = axes[row, j]
        kw = {2: dict(vmin=0, vmax=vmax_ph), 3: dict(vmin=vmin_net, vmax=vmax_net)}.get(j, {})
        ax.imshow(p, cmap=cm, interpolation='nearest', **kw)
        if j == 3:
            ax.contour(net[c, c], levels=[thr], colors='white', linewidths=1.2)
        ax.set_xticks([]); ax.set_yticks([])
        if row == 0: ax.set_title(cols[j], fontsize=11)
    label = ('Trial A: hole opens' if not failed else 'Trial B: SAME settings,\nhole stays closed')
    axes[row, 0].set_ylabel(label, fontsize=12, fontweight='bold')
    axes[row, 2].text(1, 2, f'{inside} photons absorbed\ninside the hole outline', color='white', fontsize=9, va='top')
    cd = 2 * np.sqrt(area / np.pi)
    axes[row, 4].text(1, 2, f'printed size {cd:.0f} nm' if area else 'hole never opened', color='black', fontsize=9, va='top', bbox=dict(facecolor='white', edgecolor='none', pad=2))
fig.suptitle('Same dose (18 mJ/cm²), same 20 nm hole, different random luck', fontsize=14, fontweight='bold')
fig.text(0.5, 0.01, 'Prototype simulator · 1 pixel = 1 nm · white line in panel 4 = development threshold', ha='center', fontsize=9, color='gray')
plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=2.0)
plt.savefig('results/hole_stages.png', dpi=130)
print('ok', ok[0][RR <= d/2].sum(), ok[3], 'bad', bad[0][RR <= d/2].sum(), bad[3])
