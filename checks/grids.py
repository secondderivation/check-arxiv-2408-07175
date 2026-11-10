"""Likelihood grids for C5, C7, C8, M1c, M2 (PROTOCOL section 1, Addendum A1: both D_L conventions).

Grid (conventions.json posteriors): Omega_m 0.01..0.81 step 0.005 (161), w0 -2..1 step 0.02 (151),
wa -10..2 step 0.05 (241). For every point and every SN data vector: chi2 profiled over the offset
(the marginalised -2lnL differs by a constant). BAO: -2 ln L marginalised over alpha for both
alpha priors, DESI DR1 and DR2. Saved to raw/grids/<name>.npy (float64; outside Git, hashes in
checks/GRIDS.sha256). Spot check against sn.Sample.chi2 at random points: checks/out/grids.out.

Usage: python3 checks/grids.py
"""
import hashlib, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn

HERE = os.path.dirname(os.path.abspath(__file__))
GDIR = os.path.join(sn.RAW, "grids")
os.makedirs(GDIR, exist_ok=True)
OM = np.round(np.linspace(0.01, 0.81, 161), 6)
W0 = np.round(np.linspace(-2.0, 1.0, 151), 6)
WA = np.round(np.linspace(-10.0, 2.0, 241), 6)
SHAPE = (len(OM), len(W0), len(WA))

# data vectors per sample (the shifted DES5Y vectors share the DES5Y covariance and redshifts)
P, D, V = sn.pantheon(), sn.des5y(), sn.dovekie()
SAMPLES = {"pantheon": (P, {"pantheon": P.m}),
           "des5y": (D, {"des5y": D.m, "des5y_shift_survey": sn.des5y(shift=0.04).m, "des5y_shift_z": sn.des5y(shift=0.04, shift_rule="z").m}),
           "dovekie": (V, {"dovekie": V.m})}
BAOS = {"dr1": sn.BAO("dr1"), "dr2": sn.BAO("dr2")}
ALPHA = {k: sn.alpha_nodes(k) for k in ("primary", "variant")}

out = {}
for sname, (s, vecs) in SAMPLES.items():
    for v in vecs:
        for zc in ("hd", "hel"):
            out[f"{v}_{zc}"] = np.empty(SHAPE)
for b in BAOS:
    for pr in ALPHA:
        out[f"bao_{b}_{pr}"] = np.empty(SHAPE)

FINE = sn.ZGRID
sn.ZGRID = np.linspace(0.0, 2.4, 1001)   # coarser table for the grid (Hermite error O(h^4) ~ 1e-9; the spot check below uses the fine table)
a3 = (1 + sn.ZGRID) ** 3
W0g, WAg = np.meshgrid(W0, WA, indexing="ij")
pairs = np.c_[W0g.ravel(), WAg.ravel()]
CH = 2048
t0 = time.time()
for k0 in range(0, int(os.environ.get('SD_GRID_LIMIT', len(pairs))), CH):
    w0, wa = pairs[k0:k0 + CH, 0], pairs[k0:k0 + CH, 1]
    F = (1 + sn.ZGRID[None, :]) ** (3 * (1 + w0 + wa))[:, None] * np.exp(-3 * wa[:, None] * sn.ZGRID[None, :] / (1 + sn.ZGRID[None, :]))
    idx = np.unravel_index(np.arange(k0, k0 + len(w0)), (len(W0), len(WA)))
    for i, om in enumerate(OM):
        inv = 1.0 / np.sqrt(om * a3[None, :] + (1 - om) * F)
        dz = sn.ZGRID[1] - sn.ZGRID[0]
        val = np.zeros_like(inv)
        val[:, 1:] = np.cumsum(0.5 * (inv[:, 1:] + inv[:, :-1]) * dz, axis=1)
        der = np.gradient(inv, dz, axis=1, edge_order=2)
        val[:, 1:] -= dz ** 2 / 12.0 * (der[:, 1:] - der[:, :1])
        tab = (val, inv)
        for sname, (s, vecs) in SAMPLES.items():
            dc = sn.interp_rows(tab, s.zhd)
            for zc in ("hd", "hel"):
                fac = (1 + s.zhd) if zc == "hd" else (1 + s.zhel)
                mth = 5 * np.log10(fac[None, :] * dc)
                X = mth @ s.icov
                q = np.einsum("ij,ij->i", X, mth)
                be = mth @ s.e
                for v, m in vecs.items():
                    A = m @ s.icov @ m - 2 * X @ m + q
                    B = s.e @ m - be
                    out[f"{v}_{zc}"][i][idx] = A - B ** 2 / s.E
        for b, bao in BAOS.items():
            dm = sn.interp_rows(tab, bao.z)
            dh = 1.0 / sn.efunc(bao.z[None, :], om, w0[:, None], wa[:, None])
            dv = (bao.z[None, :] * dm ** 2 * dh) ** (1 / 3)
            S = np.empty_like(dm)
            for j, kind in enumerate(bao.kind):
                S[:, j] = {"DM_over_rs": dm, "DH_over_rs": dh, "DV_over_rs": dv}[kind][:, j]
            for pr in ALPHA:
                out[f"bao_{b}_{pr}"][i][idx] = sn.bao_m2lnL(bao, S, pr)
    el = time.time() - t0
    print(f"{k0 + len(w0)}/{len(pairs)} pairs, {el / 60:.1f} min", flush=True)

sn.ZGRID = FINE
lines = []
rng = np.random.default_rng(2408)
for name in ("pantheon_hd", "des5y_hel", "des5y_shift_survey_hd", "dovekie_hd"):
    worst = 0
    for _ in range(20):
        i, j, k = rng.integers(SHAPE[0]), rng.integers(SHAPE[1]), rng.integers(SHAPE[2])
        base = name.rsplit("_", 1)
        s = {"pantheon": P, "des5y": D, "dovekie": V, "des5y_shift_survey": sn.des5y(shift=0.04)}[base[0]]
        worst = max(worst, abs(s.chi2(OM[i], W0[j], WA[k], zconv=base[1])[0][0] - out[name][i, j, k]))
        assert np.isfinite(out[name]).all()
    lines.append(f"spot check {name}: max |grid - sn.Sample.chi2| over 20 random points = {worst:.2e}")
bao = BAOS["dr1"]
worst = 0
for _ in range(20):
    i, j, k = rng.integers(SHAPE[0]), rng.integers(SHAPE[1]), rng.integers(SHAPE[2])
    worst = max(worst, abs(sn.bao_m2lnL(bao, bao.shapes(OM[i], W0[j], WA[k]), "primary")[0] - out["bao_dr1_primary"][i, j, k]))
lines.append(f"spot check bao_dr1_primary: max |grid - sn.bao_m2lnL| = {worst:.2e}")
sums = []
for name, arr in out.items():
    p = os.path.join(GDIR, name + ".npy")
    np.save(p, arr)
    sums.append(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {name}.npy")
open(os.path.join(HERE if not sn.REPLAY else sn.OUT, "GRIDS.sha256"), "w").write("\n".join(sums) + "\n")
lines.append(f"grid {SHAPE}, {len(out)} arrays")
print("\n".join(lines))
print(f"run time {(time.time() - t0) / 60:.1f} min")   # stdout only, so grids.out is the same on every run
open(os.path.join(sn.OUT, "grids.out"), "w").write("\n".join(lines) + "\n")
