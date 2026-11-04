"""Stochastic EUV contact-hole simulator, version 2.

A 2D "column" model: each 1 nm x 1 nm pixel stands for a column through the
resist film. It follows the structure of published stochastic resist models
(Naulleau/Gallatin MPPM; Mack/Biafore PROLITH stochastic model):

  1. Aerial image      partially coherent imaging (Abbe sum over source points)
                       of a periodic array of mask holes, NA 0.33, 13.5 nm,
                       normalised so an open frame = 1 (scanner dose units).
  2. Photons           absorbed photons per column ~ Poisson
                       (dose * I * absorbed fraction of a film of thickness t).
  3. Electron blur     secondary electrons spread the exposure (Gaussian).
  4. Acid              PAG molecules per column ~ Poisson(rho_PAG * t); each
                       converts with p = 1 - exp(-QE * n_abs / (rho_PAG * t))
                       (first-order exposure kinetics, Dill C).  By Poisson
                       thinning, acids per column ~ Poisson(rho_PAG * t * p).
  5. Quencher          quencher per column ~ Poisson(rho_Q * t).
  6. Bake (PEB)        acid and quencher diffuse (Gaussian blurs) and
                       neutralise instantly; net acid h (per nm^3) drives
                       deprotection: protected fraction m = exp(-K * h).
  7. Protecting groups per column ~ Poisson(rho_PG * t); deprotected and
                       remaining counts are independent Poisson (thinning).
  8. Development       a column dissolves when its deprotected fraction
                       exceeds the development threshold (0.5).

Per hole it reports the printed CD (equivalent-circle diameter), whether the
hole is missing (cleared area < 25% of target area) and whether it merged with
a neighbour.  By default (hole_metric="component") a hole is the largest
connected dissolved region in its cell (4-connected, periodic tile), its area
is the size of that region, and it is merged when that region is also the
largest region of another hole.  Stray dissolved specks are ignored.  The
older rule (hole_metric="cell": all dissolved pixels in the cell count, and
any dissolved pixels touching across a cell border count as a merge) is kept
for comparison.

Every noise source can be switched off for tests and noise budgets.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from functools import lru_cache

import numpy as np
from scipy.ndimage import gaussian_filter, label as ndi_label
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

# ---------------------------------------------------------------- constants
WAVELENGTH_NM = 13.5
PHOTON_EV = 1239.84 / WAVELENGTH_NM          # ~91.8 eV
J_PER_EV = 1.602176634e-19
# 1 mJ/cm^2 = 1e-17 J/nm^2
PHOTONS_PER_NM2_PER_MJCM2 = 1e-17 / (PHOTON_EV * J_PER_EV)   # ~0.680


@dataclass(frozen=True)
class Params:
    # ---- pattern
    cd_target: float = 20.0       # nm, target printed hole diameter
    mask_cd: float | None = None  # nm, square mask hole size (None = cd_target)
    pitch_x: int = 40             # nm
    pitch_y: int = 40             # nm
    stagger: bool = False         # offset alternate rows by pitch_x / 2
    cells_x: int = 3              # holes per tile in x
    cells_y: int = 3              # holes per tile (rows) in y
    # ---- optics
    na: float = 0.33
    illumination: str = "quasar"  # "quasar" (square arrays), "hexapole" (staggered), "annular", "interference"
    pole_sigma: float = 0.10      # quasar pole radius
    sigma_in: float = 0.40        # annular
    sigma_out: float = 0.85
    pole_samples: int = 7         # source points across each pole (convergence tests)
    flare: float = 0.0            # fraction of open-frame dose added uniformly
    # ---- exposure
    dose: float = 30.0            # mJ/cm^2 (scanner / open-frame dose)
    # ---- resist (defaults: chemically amplified resist, Naulleau et al.)
    thickness: float = 35.0       # nm
    absorptivity: float = 0.0042  # 1/nm  (4.2 / um)
    qe: float = 2.8               # acids per absorbed photon at low dose
    rho_pag: float = 0.20         # 1/nm^3
    rho_q: float = 0.09           # 1/nm^3
    rho_pg: float = 2.0           # 1/nm^3 protecting groups
    sigma_e: float = 1.7          # nm, electron blur
    sigma_a: float = 5.0          # nm, acid diffusion blur
    sigma_q: float = 1.3          # nm, quencher diffusion blur
    k_amp: float = 60.0           # nm^3, deprotection rate x bake time
    dev_threshold: float = 0.5    # deprotected fraction needed to dissolve
    layers: int = 1               # layers through the film (1 = 2D column model)
    px: float = 1.0               # nm per pixel (for convergence tests; pitches must be multiples)
    # ---- failure definitions
    missing_area_frac: float = 0.25
    hole_metric: str = "component"  # "component" (default) or "cell" (old rule)
    # ---- noise switches
    noise: tuple = ("photon", "pag", "quencher", "pg")

    @property
    def absorbed_frac(self) -> float:
        return 1.0 - np.exp(-self.absorptivity * self.thickness)


# ---------------------------------------------------------------- geometry
def _hole_centres(p: Params):
    cx, cy = [], []
    for j in range(p.cells_y):
        off = (p.pitch_x / 2) if (p.stagger and j % 2) else 0.0
        for i in range(p.cells_x):
            cx.append(i * p.pitch_x + p.pitch_x / 2 + off)
            cy.append(j * p.pitch_y + p.pitch_y / 2)
    return np.array(cx), np.array(cy)


def _tile_shape(p: Params):
    if p.stagger and p.cells_y % 2:
        raise ValueError("staggered arrays need an even number of rows")
    H, W = p.cells_y * p.pitch_y / p.px, p.cells_x * p.pitch_x / p.px
    if abs(H - round(H)) > 1e-9 or abs(W - round(W)) > 1e-9:
        raise ValueError("pitch must be a multiple of px")
    return int(round(H)), int(round(W))   # pixels


@lru_cache(maxsize=64)
def _geometry(p: Params):
    """Mask, cell map (nearest hole, periodic) for a parameter set."""
    H, W = _tile_shape(p)
    yy, xx = (np.mgrid[:H, :W] + 0.5) * p.px
    Hn, Wn = H * p.px, W * p.px
    cx, cy = _hole_centres(p)
    m = p.mask_cd if p.mask_cd is not None else p.cd_target
    mask = np.zeros((H, W))
    best = np.full((H, W), np.inf)
    cell = np.zeros((H, W), dtype=np.int32)
    for k, (x0, y0) in enumerate(zip(cx, cy)):
        dx = (xx - x0 + Wn / 2) % Wn - Wn / 2
        dy = (yy - y0 + Hn / 2) % Hn - Hn / 2
        # exact fraction of each pixel covered by the square mask hole (anti-aliased),
        # so mask sizes are not rounded to whole pixel pairs
        cov_x = np.clip(m / 2 - np.abs(dx) + p.px / 2, 0.0, p.px) / p.px
        cov_y = np.clip(m / 2 - np.abs(dy) + p.px / 2, 0.0, p.px) / p.px
        mask = np.maximum(mask, cov_x * cov_y)
        d2 = dx * dx + dy * dy
        closer = d2 < best
        best[closer] = d2[closer]
        cell[closer] = k
    return mask, cell, len(cx)


# ---------------------------------------------------------------- optics
def _source_points(p: Params):
    """Source points (sigma_x, sigma_y) and weights."""
    pts = []
    if p.illumination == "quasar":
        # Centre the four 45-degree poles so the 0th and first diagonal orders
        # sit symmetrically in the pupil (standard choice for square arrays).
        sc = WAVELENGTH_NM / (np.sqrt(2) * min(p.pitch_x, p.pitch_y) * p.na)
        sc = float(np.clip(sc, 0.3, 0.9 - p.pole_sigma))
        r = np.linspace(-p.pole_sigma, p.pole_sigma, p.pole_samples)
        for sx_c, sy_c in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
            c = sc / np.sqrt(2)
            for ax in r:
                for ay in r:
                    if ax * ax + ay * ay <= p.pole_sigma ** 2 + 1e-12:
                        pts.append((sx_c * c + ax, sy_c * c + ay))
    elif p.illumination == "hexapole":
        # Staggered (near-hexagonal) arrays: six poles at half of each first
        # diffraction order g, so every pole images the 0th order and the
        # opposite first order symmetrically (the hexagonal analogue of the
        # quasar rule).  Lattice a1 = (px, 0), a2 = (px/2, py).
        if not p.stagger:
            raise ValueError("hexapole illumination is for staggered arrays")
        b1 = np.array([1.0 / p.pitch_x, -0.5 / p.pitch_y])
        b2 = np.array([0.0, 1.0 / p.pitch_y])
        r = np.linspace(-p.pole_sigma, p.pole_sigma, p.pole_samples)
        for g in (b1, -b1, b2, -b2, b1 + b2, -(b1 + b2)):
            c = WAVELENGTH_NM * g / (2 * p.na)
            rad = np.hypot(*c)
            c = c * float(np.clip(rad, 0.3, 0.9 - p.pole_sigma)) / rad
            for ax in r:
                for ay in r:
                    if ax * ax + ay * ay <= p.pole_sigma ** 2 + 1e-12:
                        pts.append((c[0] + ax, c[1] + ay))
    elif p.illumination == "interference":
        return None, None             # handled analytically in aerial_image
    elif p.illumination == "annular":
        for rad in np.linspace(p.sigma_in, p.sigma_out, 5):
            n = max(8, int(40 * rad))
            for t in np.linspace(0, 2 * np.pi, n, endpoint=False):
                pts.append((rad * np.cos(t), rad * np.sin(t)))
    else:
        raise ValueError(p.illumination)
    pts = np.array(pts)
    return pts, np.full(len(pts), 1.0 / len(pts))


@lru_cache(maxsize=64)
def aerial_image(p: Params) -> np.ndarray:
    """Intensity relative to an open frame (1.0 = full scanner dose)."""
    mask, _, _ = _geometry(p)
    H, W = mask.shape
    if p.illumination == "interference":
        # Double-exposure four-beam interference (PSI-style EUV-IL): two crossed
        # two-beam line gratings, each a sinusoid; holes sit at the shared peaks.
        # Scaled so the average intensity equals the stated dose (1.0).
        yy, xx = (np.mgrid[:H, :W] + 0.5) * p.px
        gx = 1 + np.cos(2 * np.pi * (xx - p.pitch_x / 2) / p.pitch_x)
        gy = 1 + np.cos(2 * np.pi * (yy - p.pitch_y / 2) / p.pitch_y)
        return 0.5 * (gx + gy)
    M = np.fft.fft2(mask)
    fy = np.fft.fftfreq(H, d=p.px)[:, None]  # cycles per nm
    fx = np.fft.fftfreq(W, d=p.px)[None, :]
    cutoff = p.na / WAVELENGTH_NM
    pts, w = _source_points(p)
    img = np.zeros((H, W))
    for (sx, sy), wk in zip(pts, w):
        pupil = ((fx + sx * cutoff) ** 2 + (fy + sy * cutoff) ** 2) <= cutoff ** 2
        field = np.fft.ifft2(M * pupil)
        img += wk * np.abs(field) ** 2
    return (1 - p.flare) * img + p.flare


# ---------------------------------------------------------------- resist
def _draw(rng, lam, on):
    return rng.poisson(lam).astype(np.float32) if on else lam.astype(np.float32)


_ST4 = np.zeros((3, 3, 3), dtype=bool)
_ST4[1] = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]   # 4-connected in-plane, never across trials


def _measure_components(clear: np.ndarray, cell: np.ndarray, n_holes: int):
    """Per-hole area (pixels) and merged flag from connected dissolved regions.

    clear: (b, H, W) bool.  The tile is periodic, so regions that touch across
    the tile edges are joined.  Each hole's body is the region with the most
    pixels inside the hole's cell; area = full size of that region.  A hole is
    merged when its body is also the body of another hole in the same tile.
    Also returns group: the body id for merged holes, else -1. Holes sharing a
    group id were joined by one merge event.
    """
    b = clear.shape[0]
    lab, n = ndi_label(clear, structure=_ST4)
    if n == 0:
        return (np.zeros((b, n_holes)), np.zeros((b, n_holes), dtype=bool),
                np.full((b, n_holes), -1, dtype=np.int64))
    # join regions across the periodic edges
    pa, pb = [], []
    for e0, e1 in ((lab[:, :, -1], lab[:, :, 0]), (lab[:, -1, :], lab[:, 0, :])):
        both = (e0 > 0) & (e1 > 0)
        pa.append(e0[both]); pb.append(e1[both])
    pa, pb = np.concatenate(pa), np.concatenate(pb)
    if len(pa):
        g = coo_matrix((np.ones(len(pa), dtype=np.int8), (pa, pb)), shape=(n + 1, n + 1))
        _, root = connected_components(g, directed=False)
    else:
        root = np.arange(n + 1)
    r = root[lab[clear]]                                   # region id per dissolved pixel
    t = np.broadcast_to(np.arange(b)[:, None, None], clear.shape)[clear]
    k = np.broadcast_to(cell[None], clear.shape)[clear]
    size = np.bincount(r)
    key = r.astype(np.int64) * n_holes + k
    uk, cnt = np.unique(key, return_counts=True)
    ur, uc = uk // n_holes, uk % n_holes
    trial_of = np.zeros(size.size, dtype=np.int64)
    trial_of[r] = t
    ut = trial_of[ur]
    o = np.lexsort((-cnt, uc, ut))                         # biggest share first per (trial, hole)
    ut, uc, ur = ut[o], uc[o], ur[o]
    first = np.ones(len(o), dtype=bool)
    first[1:] = (ut[1:] != ut[:-1]) | (uc[1:] != uc[:-1])
    body = np.full((b, n_holes), -1, dtype=np.int64)
    body[ut[first], uc[first]] = ur[first]
    area = np.where(body >= 0, size[np.maximum(body, 0)], 0).astype(float)
    same = (body[:, :, None] == body[:, None, :]) & ~np.eye(n_holes, dtype=bool)[None]
    merged = same.any(axis=2) & (body >= 0)
    group = np.where(merged, body, -1)                    # same id = holes joined by one merge event
    return area, merged, group


def simulate(p: Params, n_trials: int, rng=None, batch: int | None = None, return_fields: bool = False):
    """Run n_trials tiles. Returns dict of per-hole arrays (n_trials * holes).

    With p.layers > 1 the film is split into layers through its thickness:
    light weakens with depth (Beer-Lambert), every layer gets its own random
    molecules, blur and diffusion act in 3D, and a column only opens if it
    clears in every layer (the developer must reach the bottom).
    """
    rng = rng if rng is not None else np.random.default_rng()
    I = aerial_image(p).astype(np.float32)
    _, cell, n_holes = _geometry(p)
    H, W = I.shape
    L = p.layers
    dz = p.thickness / L
    nz = set(p.noise)
    # share of incident photons absorbed in each layer (sums to absorbed_frac)
    layer_abs = (np.exp(-p.absorptivity * dz * np.arange(L))
                 * (1 - np.exp(-p.absorptivity * dz))).astype(np.float32)
    a_px = p.px * p.px                                          # pixel area, nm^2
    mean_abs = p.dose * PHOTONS_PER_NM2_PER_MJCM2 * a_px * layer_abs[:, None, None] * I[None]
    vox = dz * a_px                                              # voxel volume, nm^3
    pag0, q0, pg0 = p.rho_pag * vox, p.rho_q * vox, p.rho_pg * vox  # per voxel
    target_area = np.pi * (p.cd_target / 2) ** 2
    cell_flat = cell.ravel()
    modes = ("wrap", "reflect", "wrap", "wrap")    # periodic in x, y; film surfaces reflect

    def blur(x, s):
        return gaussian_filter(x, (0, s / dz, s / p.px, s / p.px), mode=modes)

    batch = batch or max(20, 250 // L)
    cds, missing, merged, groups = [], [], [], []
    fields = None
    done = 0
    while done < n_trials:
        b = min(batch, n_trials - done)
        lam = np.broadcast_to(mean_abs, (b, L, H, W))
        photons = _draw(rng, lam, "photon" in nz)
        exposure = blur(photons, p.sigma_e)
        p_conv = 1.0 - np.exp(-p.qe * exposure / pag0)
        acid = _draw(rng, pag0 * p_conv, "pag" in nz)
        quench = _draw(rng, np.full((b, L, H, W), q0, np.float32), "quencher" in nz)
        a_b = blur(acid, p.sigma_a)
        q_b = blur(quench, p.sigma_q)
        h = np.maximum(a_b - q_b, 0.0) / vox                # net acid per nm^3
        m = np.exp(-p.k_amp * h)                             # protected fraction
        if "pg" in nz:
            dep = rng.poisson(pg0 * (1 - m)).astype(np.float32)
            rem = rng.poisson(pg0 * m).astype(np.float32)
            frac = dep / np.maximum(dep + rem, 1.0)
        else:
            frac = 1 - m
        clear = (frac > p.dev_threshold).all(axis=1)         # open through every layer
        photons, acid, h, frac = photons.sum(1), acid.sum(1), h.mean(1), frac.min(1)

        if p.hole_metric == "component":
            area, mg, grp = _measure_components(clear, cell, n_holes)
            # make merge-event ids unique across batches
            groups.append(np.where(grp >= 0, grp + (done + 1) * (1 << 40), -1))
        elif p.hole_metric == "cell":
            # old rule: every dissolved pixel in the cell counts ...
            flat = clear.reshape(b, -1)
            area = np.zeros((b, n_holes))
            for k in range(n_holes):
                area[:, k] = flat[:, cell_flat == k].sum(axis=1)
            # ... and dissolved pixels touching across a cell border = merged
            mg = np.zeros((b, n_holes), dtype=bool)
            for axis in (1, 2):
                nb = np.roll(clear, -1, axis=axis)
                cnb = np.roll(cell, -1, axis=axis - 1)
                bridge = clear & nb & (cell != cnb)[None]
                if bridge.any():
                    idx_b, yy_, xx_ = np.nonzero(bridge)
                    mg[idx_b, cell[yy_, xx_]] = True
                    mg[idx_b, cnb[yy_, xx_]] = True
        else:
            raise ValueError(f"unknown hole_metric {p.hole_metric!r}")
        area = area * a_px                                   # nm^2
        cds.append(2 * np.sqrt(area / np.pi))
        missing.append(area < p.missing_area_frac * target_area)
        merged.append(mg)
        if p.hole_metric != "component":
            groups.append(np.full((b, n_holes), -1, dtype=np.int64))
        if return_fields and fields is None:
            fields = dict(image=I, photons=photons[0], acid=acid[0], net_acid=h[0],
                          deprotected=frac[0], clear=clear[0], cell=cell)
        done += b

    out = dict(cd=np.concatenate(cds).ravel(), missing=np.concatenate(missing).ravel(),
               merged=np.concatenate(merged).ravel(), merge_group=np.concatenate(groups).ravel())
    if return_fields:
        out["fields"] = fields
    return out


def summarize(res: dict) -> dict:
    ok = ~(res["missing"] | res["merged"])
    cd = res["cd"][ok]
    return dict(n_holes=res["cd"].size,
                mean_cd=float(cd.mean()) if cd.size else float("nan"),
                sigma_cd=float(cd.std()) if cd.size else float("nan"),
                lcdu_3s=float(3 * cd.std()) if cd.size else float("nan"),
                p_missing=float(res["missing"].mean()),
                p_merged=float(res["merged"].mean()))


# knob: (lower bound, upper bound, +1 if a larger value prints a larger hole)
KNOBS = {
    "rho_q": (0.002, 0.30, -1),   # quencher loading: how formulators set resist speed
    "k_amp": (0.5, 5000.0, +1),   # deprotection rate x bake time
    "dose": (2.0, 400.0, +1),     # classic dose-to-size
}


def dose_to_size(p: Params, knob: str = "rho_q", n_trials: int = 120, rng_seed: int = 0,
                 iters: int = 14) -> Params:
    """Return params whose mean printed CD equals cd_target, by tuning one knob.

    knob="rho_q" (default): keep the published dose and tune the resist's
    quencher loading, i.e. pick the resist formulation that prints to size at
    that dose. knob="dose": keep the resist and tune the dose.
    Bisection in log space; common random numbers keep it monotone.
    """
    lo_v, hi_v, sign = KNOBS[knob]
    lo, hi = np.log(lo_v), np.log(hi_v)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        q = replace(p, **{knob: float(np.exp(mid))})
        s = summarize(simulate(q, n_trials, np.random.default_rng(rng_seed)))
        if np.isfinite(s["mean_cd"]):
            cd = s["mean_cd"]
        else:   # every hole failed: decide which way
            cd = np.inf if s["p_merged"] > s["p_missing"] else 0.0
        too_small = cd < p.cd_target
        if too_small == (sign > 0):
            lo = mid
        else:
            hi = mid
    return replace(p, **{knob: float(np.exp(0.5 * (lo + hi)))})


if __name__ == "__main__":
    p = dose_to_size(Params())
    s = summarize(simulate(p, 400, np.random.default_rng(1)))
    print(f"quencher tuned to {p.rho_q:.4f} /nm^3 for {p.cd_target} nm holes at {p.dose} mJ/cm^2")
    print({k: round(v, 4) if isinstance(v, float) else v for k, v in s.items()})
