"""C1-C4, C7 (best fits) and M1b: counts, chi2 at the Planck point, LCDM and w0waCDM best fits,
Omega_m with errors, Delta chi2 and sigma conversions. Output: checks/out/fits.out, fits.json (G14).

Usage: python3 checks/fits.py
"""
import json, os, sys
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import chi2 as chi2dist
from scipy.special import erfinv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn

OUT = sn.OUT
lines, res = [], {}
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); lines.append(s)

def nsig(dchi2, k):
    p = chi2dist.sf(dchi2, k)
    return float(np.sqrt(2) * erfinv(1 - p))

STARTS = [(0.3, -1, 0), (0.35, -0.8, -1), (0.45, -0.6, -5), (0.25, -1.1, 1)]
BOX = {"Pantheon+": [(0.2, 0.6), (-1.4, 0.4), (-4, -0.4)], "DES5Y": [(0.01, 0.81), (-2, 1), (-10, 1)]}

def chi2f(s, **kw):
    return lambda x: float(s.chi2(*x, **kw)[0][0]) if 0 < x[0] < 1 else 1e9

def fit_lcdm(s):
    r = minimize_scalar(lambda om: float(s.chi2(om)[0][0]), bounds=(0.01, 0.99), method="bounded", options={"xatol": 1e-7})
    return r.x, r.fun

def fit_w0wa(s, box=None):
    """Five starts (four fixed + lowest point of a coarse scan); Nelder-Mead then Powell."""
    om, w0, wa = np.meshgrid(np.linspace(0.05, 0.75, 15), np.linspace(-1.8, 0.8, 14), np.linspace(-9.5, 1.5, 12), indexing="ij")
    pts = np.c_[om.ravel(), w0.ravel(), wa.ravel()]
    if box is not None:
        pts = pts[np.all([(pts[:, i] >= box[i][0]) & (pts[:, i] <= box[i][1]) for i in range(3)], axis=0)]
    c = np.concatenate([s.chi2(*pts[i:i + 400].T)[0] for i in range(0, len(pts), 400)])
    starts = STARTS + [tuple(pts[np.argmin(c)])]
    f = chi2f(s)
    if box is not None:
        g = f
        f = lambda x: g(x) if all(box[i][0] <= x[i] <= box[i][1] for i in range(3)) else 1e9
    mins = []
    for x0 in starts:
        r = minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-7, "maxiter": 20000, "maxfev": 20000})
        r = minimize(f, r.x, method="Powell", options={"xtol": 1e-6, "ftol": 1e-9}) if box is None else r
        mins.append((r.fun, tuple(r.x)))
    mins.sort()
    return mins

def om_errors(s):
    """Posterior mean/sd (flat Omega_m prior [0.01, 0.99], offset marginalised) and profile half-width."""
    g = np.linspace(0.01, 0.99, 9801)
    c = np.concatenate([s.chi2(g[i:i + 2000])[0] for i in range(0, len(g), 2000)])
    p = np.exp(-0.5 * (c - c.min())); p /= p.sum()
    mean = (g * p).sum(); sd = np.sqrt(((g - mean) ** 2 * p).sum())
    inside = g[c <= c.min() + 1.0]
    return mean, sd, (inside.max() - inside.min()) / 2, np.median(np.repeat(g, 1))  # last unused

def report(s, label, printed=None, zconv="hel"):
    say(f"\n== {label}: N = {s.n}")
    s.zconv = zconv
    cP, offP = s.chi2(sn.PLANCK_OM)
    cP = float(cP[0])
    say(f"  chi2 at the Planck point (Omega_m = {sn.PLANCK_OM}): {cP:.3f}   (offset M' = {float(offP[0]):.4f})")
    om, cL = fit_lcdm(s)
    mean, sd, hw, _ = om_errors(s)
    say(f"  LCDM best fit: Omega_m = {om:.4f}, chi2 = {cL:.3f}; Delta chi2 (Planck - LCDM) = {cP - cL:.3f}")
    say(f"  LCDM Omega_m: posterior mean {mean:.4f} +- {sd:.4f} (sd); profile half-width (Delta chi2 = 1) {hw:.4f}")
    mins = fit_w0wa(s)
    cW, xW = mins[0]
    say(f"  w0waCDM best fit: Omega_m = {xW[0]:.4f}, w0 = {xW[1]:.4f}, wa = {xW[2]:.4f}, chi2 = {cW:.3f}")
    say("    the five minima:", "; ".join(f"{c:.3f} at ({x[0]:.3f}, {x[1]:.3f}, {x[2]:.3f})" for c, x in mins))
    out = {"N": s.n, "chi2_planck": cP, "lcdm": {"Omega_m": om, "chi2": cL, "mean": mean, "sd": sd, "profile_hw": hw},
           "w0wa": {"Omega_m": xW[0], "w0": xW[1], "wa": xW[2], "chi2": cW, "minima": [[c, list(x)] for c, x in mins]}}
    d1, d3, d2 = cP - cL, cP - cW, cL - cW
    say(f"  Delta chi2: Planck - w0wa = {d3:.3f}; LCDM - w0wa = {d2:.3f}")
    say(f"  sigma: Planck vs LCDM fit k=1 {nsig(d1, 1):.2f}; Planck vs w0wa k=3 {nsig(d3, 3):.2f} (k=2 {nsig(d3, 2):.2f}, k=1 {nsig(d3, 1):.2f}, sqrt {np.sqrt(max(d3, 0)):.2f}); "
        f"LCDM vs w0wa k=2 {nsig(d2, 2):.2f} (k=1 {nsig(d2, 1):.2f}, sqrt {np.sqrt(max(d2, 0)):.2f})")
    gs = abs(mean - sn.PLANCK_OM) / np.hypot(sd, 0.0081)
    say(f"  LCDM Omega_m vs Planck: |mean - 0.3135|/sqrt(sd^2 + 0.0081^2) = {gs:.2f}; best fit: {abs(om - sn.PLANCK_OM) / np.hypot(sd, 0.0081):.2f}; /sd only {abs(mean - sn.PLANCK_OM) / sd:.2f}")
    out["dchi2"] = {"planck_lcdm": d1, "planck_w0wa": d3, "lcdm_w0wa": d2}
    out["sigma"] = {"planck_lcdm_k1": nsig(d1, 1), "planck_w0wa_k3": nsig(d3, 3), "planck_w0wa_k2": nsig(d3, 2), "planck_w0wa_k1": nsig(d3, 1),
                    "lcdm_w0wa_k2": nsig(d2, 2), "lcdm_w0wa_k1": nsig(d2, 1), "om_gauss_mean": gs}
    if label in BOX:
        mb = fit_w0wa(s, BOX[label])[0]
        say(f"  w0waCDM best fit inside the printed prior box: chi2 = {mb[0]:.3f} at ({mb[1][0]:.4f}, {mb[1][1]:.4f}, {mb[1][2]:.4f})")
        out["w0wa_in_box"] = {"chi2": mb[0], "x": list(mb[1])}
    # conventions variants at the best fits
    for zc in ({"hel": "hd", "hd": "hel"}[zconv],):
        say(f"  variant (1+z{zc.upper()}): Planck {s.chi2(sn.PLANCK_OM, zconv=zc)[0][0]:.3f}, at LCDM best fit {s.chi2(om, zconv=zc)[0][0]:.3f}, at w0wa best fit {s.chi2(*xW, zconv=zc)[0][0]:.3f}")
    res[label] = out
    return out

# C1: counts
import pandas as pd
P_all, D_all = sn.pantheon(cut=False), sn.des5y(cut=False)
say("== C1 counts")
say(f"  Pantheon+ entries: {P_all.n}; IDSURVEY 10 (DES Y3): {(P_all.meta.IDSURVEY == 10).sum()}; unique CIDs: {P_all.meta.CID.nunique()}")
dm = D_all.meta
say(f"  DES5Y SN: {D_all.n}; IDSURVEY 10: {(dm.IDSURVEY == 10).sum()} with zHD {dm.zHD[dm.IDSURVEY == 10].min()}..{dm.zHD[dm.IDSURVEY == 10].max()}; others: {(dm.IDSURVEY != 10).sum()} with zHD max {dm.zHD[dm.IDSURVEY != 10].max()}")
say("  DES5Y low-z by IDSURVEY:", dict(dm.IDSURVEY[dm.IDSURVEY != 10].value_counts().sort_index()))
say(f"  DES5Y MUERR_FINAL > 1: {(dm.MUERR_FINAL > 1).sum()}; zHD outside [0.02, 1.2]: {((dm.zHD < 0.02) | (dm.zHD > 1.2)).sum()}")
say(f"  after cuts: Pantheon+ {sn.pantheon().n} (printed 1417), DES5Y {sn.des5y().n} (printed 1754)")
res["C1"] = {"pantheon_all": P_all.n, "pantheon_des": int((P_all.meta.IDSURVEY == 10).sum()), "des_all": D_all.n,
             "des_des": int((dm.IDSURVEY == 10).sum()), "des_lowz": int((dm.IDSURVEY != 10).sum()),
             "pantheon_cut": sn.pantheon().n, "des_cut": sn.des5y().n}

PRINTED = {  # paper_values.json figures: legend best fits (Figs 1 and 4)
    "Pantheon+": ((0.333,), (0.327, -1.01, 0.18), 1258.4, 1258.4),
    "DES5Y": ((0.349,), (0.474, -0.54, -6.56), 1639.6, 1632.2),
    "DES5Y shifted (IDSURVEY != 10)": ((0.310,), (0.426, -1.068, -2.63), 1635.7, 1633.4),
}
SAMPLES = [("Pantheon+", lambda: sn.pantheon()), ("DES5Y", lambda: sn.des5y()),
           ("DES5Y shifted (IDSURVEY != 10)", lambda: sn.des5y(shift=0.04)),
           ("DES5Y shifted (zHD < 0.1)", lambda: sn.des5y(shift=0.04, shift_rule="z")),
           ("M1 DES-Dovekie", lambda: sn.dovekie())]
say("\n== New check (Addendum A1): chi2 at the printed best-fit parameters, both conventions")
for label, (pl, pw, cl, cw) in PRINTED.items():
    s = dict(SAMPLES)[label]()
    for zc in ("hel", "hd"):
        a, b = s.chi2(*pl, zconv=zc)[0][0], s.chi2(*pw, zconv=zc)[0][0]
        say(f"  {label:32s} (1+z{zc.upper()}): LCDM {pl} chi2 {a:.3f} (printed {cl}); w0wa {pw} chi2 {b:.3f} (printed {cw})")
        res.setdefault("printed_params", {}).setdefault(label, {})[zc] = {"lcdm": a, "w0wa": b}
for zc in ("hd", "hel"):
    say(f"\n######## convention (1+z{zc.upper()}) " + ("(paper's, Addendum A1)" if zc == "hd" else "(Pantheon+/DES/Cobaya)"))
    for label, mk in SAMPLES:
        res_key = f"{label} [{zc}]"
        out = report(mk(), label, zconv=zc)
        res[res_key] = res.pop(label)

json.dump(res, open(os.path.join(OUT, "fits.json"), "w"), indent=1, default=float)
open(os.path.join(OUT, "fits.out"), "w").write("\n".join(lines) + "\n")
