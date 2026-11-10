"""C0.4 (posterior part): the C5 quantities on a sub-box at half the grid step (PROTOCOL section 2).

Sub-box fixed before any posterior was seen: Omega_m 0.25..0.55, w0 -1.4..0.0, wa -8..1 (inside both the
DES5Y box and the Pantheon+ wa-variant box). Pantheon+ and DES5Y under the paper's (1+z_HD). On this
sub-box the posterior (truncated to it) is computed from the main grid nodes (step 0.005, 0.02, 0.05, read
from raw/grids) and from a half-step grid (0.0025, 0.01, 0.025, computed here with the same code path as
checks/grids.py); means, sds and the Planck-point levels (definition as checks/posteriors.py) are compared.
A level change > 0.01 is reported as numerical uncertainty. Output: checks/out/grid_convergence.out
"""
import os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn

HERE = os.path.dirname(os.path.abspath(__file__))
GDIR = os.path.join(sn.RAW, "grids")
OM = np.round(np.linspace(0.01, 0.81, 161), 6)
W0 = np.round(np.linspace(-2.0, 1.0, 151), 6)
WA = np.round(np.linspace(-10.0, 2.0, 241), 6)
SUB = [(0.25, 0.55), (-1.4, 0.0), (-8.0, 1.0)]
PLANCK = (sn.PLANCK_OM, -1.0, 0.0)
PROJ = {"Om-w0": (0, 1), "Om-wa": (0, 2), "w0-wa": (1, 2)}


def bilinear(p2, ax1, ax2, x, y):
    i = int(np.clip(np.searchsorted(ax1, x) - 1, 0, len(ax1) - 2)); j = int(np.clip(np.searchsorted(ax2, y) - 1, 0, len(ax2) - 2))
    tx = (x - ax1[i]) / (ax1[i + 1] - ax1[i]); ty = (y - ax2[j]) / (ax2[j + 1] - ax2[j])
    return p2[i, j] * (1 - tx) * (1 - ty) + p2[i + 1, j] * tx * (1 - ty) + p2[i, j + 1] * (1 - tx) * ty + p2[i + 1, j + 1] * tx * ty


def summary(chi2, ax):
    post = np.exp(-0.5 * (chi2 - chi2.min())); post /= post.sum()
    out = {}
    for k, nm in enumerate(("Omega_m", "w0", "wa")):
        p = post.sum(axis=tuple(i for i in range(3) if i != k)); m = (ax[k] * p).sum()
        out[nm] = (m, np.sqrt(((ax[k] - m) ** 2 * p).sum()))
    for proj, (a, b) in PROJ.items():
        p2 = post.sum(axis=3 - a - b)
        d = bilinear(p2, ax[a], ax[b], PLANCK[a], PLANCK[b])
        out[proj] = float(p2[p2 > d].sum())
    return out


# coarse: main-grid nodes in the sub-box
sel = [(a >= lo - 1e-9) & (a <= hi + 1e-9) for a, (lo, hi) in zip((OM, W0, WA), SUB)]
coarse_ax = [a[s] for a, s in zip((OM, W0, WA), sel)]
# fine: half step over the same sub-box
fine_ax = [np.round(np.linspace(lo, hi, 2 * (len(c) - 1) + 1), 6) for (lo, hi), c in zip(SUB, coarse_ax)]
SAMPLES = {"pantheon": sn.pantheon(), "des5y": sn.des5y()}
fine = {k: np.empty(tuple(len(a) for a in fine_ax)) for k in SAMPLES}

sn.ZGRID = np.linspace(0.0, 2.4, 1001)   # as checks/grids.py
a3 = (1 + sn.ZGRID) ** 3
dz = sn.ZGRID[1] - sn.ZGRID[0]
W0g, WAg = np.meshgrid(fine_ax[1], fine_ax[2], indexing="ij")
pairs = np.c_[W0g.ravel(), WAg.ravel()]
t0 = time.time()
CH = 2048
for k0 in range(0, len(pairs), CH):
    w0, wa = pairs[k0:k0 + CH, 0], pairs[k0:k0 + CH, 1]
    F = (1 + sn.ZGRID[None, :]) ** (3 * (1 + w0 + wa))[:, None] * np.exp(-3 * wa[:, None] * sn.ZGRID[None, :] / (1 + sn.ZGRID[None, :]))
    idx = np.unravel_index(np.arange(k0, k0 + len(w0)), (len(fine_ax[1]), len(fine_ax[2])))
    for i, om in enumerate(fine_ax[0]):
        inv = 1.0 / np.sqrt(om * a3[None, :] + (1 - om) * F)
        val = np.zeros_like(inv)
        val[:, 1:] = np.cumsum(0.5 * (inv[:, 1:] + inv[:, :-1]) * dz, axis=1)
        der = np.gradient(inv, dz, axis=1, edge_order=2)
        val[:, 1:] -= dz ** 2 / 12.0 * (der[:, 1:] - der[:, :1])
        for name, s in SAMPLES.items():
            dc = sn.interp_rows((val, inv), s.zhd)
            mth = 5 * np.log10((1 + s.zhd)[None, :] * dc)
            A, B = s.chi2_parts(mth)
            fine[name][i][idx] = A - B ** 2 / s.E
    print(f"{k0 + len(w0)}/{len(pairs)} pairs, {(time.time() - t0) / 60:.1f} min", flush=True)

lines = [f"sub-box {SUB}; coarse {[len(a) for a in coarse_ax]} nodes, half-step {[len(a) for a in fine_ax]} nodes; (1+zHD); {(time.time() - t0) / 60:.1f} min"]
for name in SAMPLES:
    c = np.load(os.path.join(GDIR, f"{name}_hd.npy"))[np.ix_(*sel)]
    # the coarse nodes are every other fine node: direct check of the two code paths
    lines.append(f"{name}: max |coarse - fine| on shared nodes = {np.abs(c - fine[name][::2, ::2, ::2]).max():.2e}")
    sc, sf = summary(c, coarse_ax), summary(fine[name], fine_ax)
    for key in ("Omega_m", "w0", "wa"):
        lines.append(f"  {key:8s} mean +- sd: step {sc[key][0]:.4f} +- {sc[key][1]:.4f}; half step {sf[key][0]:.4f} +- {sf[key][1]:.4f}; change {sf[key][0] - sc[key][0]:+.4f}, {sf[key][1] - sc[key][1]:+.4f}")
    for key in PROJ:
        flag = "  > 0.01: numerical uncertainty" if abs(sf[key] - sc[key]) > 0.01 else ""
        lines.append(f"  Planck-point level {key} (truncated to the sub-box): step {sc[key]:.4f}; half step {sf[key]:.4f}; change {sf[key] - sc[key]:+.4f}{flag}")
print("\n".join(lines))
open(os.path.join(HERE, "out", "grid_convergence.out"), "w").write("\n".join(lines) + "\n")
