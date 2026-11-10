"""C5, C7 (posterior part), C8, M1c, M2: grid posteriors from the saved likelihood grids (checks/grids.py).

PROTOCOL section 1 and 3, conventions.json posteriors/priors/bao, Addendum A1 (both D_L conventions;
outcomes judged under the paper's (1+z_HD)). Posterior = exp(-(chi2_SN + m2lnL_BAO)/2) on the grid nodes
inside the prior box (box edges included). 2D credible level of a point = posterior mass of the cells
whose 2D marginal density exceeds the bilinearly interpolated density at the point; a point outside
the prior box has level 1 by definition (reported as 'outside prior').
Output: checks/out/posteriors.out, posteriors.json. Usage: python3 checks/posteriors.py
"""
import json, os, sys
import numpy as np
from scipy.special import erfinv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sn.OUT
GDIR = os.path.join(sn.RAW, "grids")
OM = np.round(np.linspace(0.01, 0.81, 161), 6)
W0 = np.round(np.linspace(-2.0, 1.0, 151), 6)
WA = np.round(np.linspace(-10.0, 2.0, 241), 6)
AX = (OM, W0, WA)
NAMES = ("Omega_m", "w0", "wa")
BOX = {"pantheon_printed": [(0.2, 0.6), (-1.4, 0.4), (-4.0, -0.4)],
       "pantheon_wa_variant": [(0.2, 0.6), (-1.4, 0.4), (-10.0, 2.0)],
       "des": [(0.01, 0.81), (-2.0, 1.0), (-10.0, 1.0)]}
PLANCK = (sn.PLANCK_OM, -1.0, 0.0)
PROJ = {"Om-w0": (0, 1), "Om-wa": (0, 2), "w0-wa": (1, 2)}
lines, res = [], {}
def say(s=""):
    print(s, flush=True); lines.append(s)

_cache = {}
def grid(name):
    if name not in _cache:
        _cache[name] = np.load(os.path.join(GDIR, name + ".npy"))
    return _cache[name]

def mask(box):
    m = [((a >= lo - 1e-9) & (a <= hi + 1e-9)) for a, (lo, hi) in zip(AX, box)]
    return m[0][:, None, None] & m[1][None, :, None] & m[2][None, None, :]

def nsig2(level):
    return float(np.sqrt(2) * erfinv(min(level, 1 - 1e-16)))

def in_box(pt, box):
    return all(lo - 1e-9 <= x <= hi + 1e-9 for x, (lo, hi) in zip(pt, box))

def bilinear(p2, ax1, ax2, x, y):
    i = int(np.clip(np.searchsorted(ax1, x) - 1, 0, len(ax1) - 2)); j = int(np.clip(np.searchsorted(ax2, y) - 1, 0, len(ax2) - 2))
    tx = (x - ax1[i]) / (ax1[i + 1] - ax1[i]); ty = (y - ax2[j]) / (ax2[j + 1] - ax2[j])
    return (p2[i, j] * (1 - tx) * (1 - ty) + p2[i + 1, j] * tx * (1 - ty) + p2[i, j + 1] * (1 - tx) * ty + p2[i + 1, j + 1] * tx * ty)

def level(post, pt, box, proj):
    a, b = PROJ[proj]
    if not (in_box([pt[a]], [box[a]]) and in_box([pt[b]], [box[b]])):
        return 1.0, "outside prior"
    other = 3 - a - b
    p2 = post.sum(axis=other)
    d = bilinear(p2, AX[a], AX[b], pt[a], pt[b])
    return float(p2[p2 > d].sum()), ""

def stats1d(post, k):
    p = post.sum(axis=tuple(i for i in range(3) if i != k)); x = AX[k]
    c = np.cumsum(p) - 0.5 * p          # mid-node cumulative mass
    q = lambda f: float(np.interp(f, c, x))
    mean = float((x * p).sum()); sd = float(np.sqrt(((x - mean) ** 2 * p).sum()))
    return {"mean": mean, "sd": sd, "median": q(0.5), "lo68": q(0.15865), "hi68": q(0.84135)}

def posterior(terms, box):
    m2 = sum(grid(t) for t in terms)
    msk = mask(box)
    m2 = np.where(msk, m2, np.inf)
    mn = m2.min()
    post = np.exp(-0.5 * (m2 - mn)); post /= post.sum()
    ijk = np.unravel_index(np.argmin(m2), m2.shape)
    return post, mn, tuple(float(AX[i][ijk[i]]) for i in range(3))

def run(label, terms, box_name, extra_points=()):
    box = BOX[box_name]
    post, mn, at = posterior(terms, box)
    r = {"terms": list(terms), "box": box_name, "grid_min": float(mn), "grid_min_at": at}
    say(f"\n== {label}  [{' + '.join(terms)}; prior {box_name}]")
    say(f"   grid minimum of -2lnL {mn:.3f} at (Om, w0, wa) = ({at[0]:.3f}, {at[1]:.2f}, {at[2]:.2f})")
    for k in range(3):
        s = stats1d(post, k); r[NAMES[k]] = s
        say(f"   {NAMES[k]:8s} mean {s['mean']:.4f} +- {s['sd']:.4f}; median {s['median']:.4f}; 68% [{s['lo68']:.4f}, {s['hi68']:.4f}]")
    pw0 = float(post[:, W0 > -1.0 + 1e-9, :].sum() + 0.5 * post[:, np.isclose(W0, -1.0), :].sum())
    r["P(w0>-1)"] = pw0
    say(f"   P(w0 > -1) = {pw0:.4f}")
    r["levels"] = {}
    for proj in PROJ:
        lv, note = level(post, PLANCK, box, proj)
        r["levels"][proj] = {"level": lv, "nsigma": nsig2(lv), "note": note}
        say(f"   Planck point (0.3135, -1, 0) in {proj}: level {lv:.4f} ({nsig2(lv):.2f} sigma) {note}")
    res[label] = r
    return r

CONV = ("hd", "hel")
for zc in CONV:
    tag = "paper's (1+zHD)" if zc == "hd" else "(1+zHEL), Pantheon+/DES/Cobaya"
    say(f"\n################ convention {tag}")
    say("\n# C5: SN-only posteriors")
    run(f"C5 Pantheon+ printed prior [{zc}]", [f"pantheon_{zc}"], "pantheon_printed")
    run(f"C5 Pantheon+ wa variant [{zc}]", [f"pantheon_{zc}"], "pantheon_wa_variant")
    run(f"C5 DES5Y [{zc}]", [f"des5y_{zc}"], "des")
    say("\n# C7: DES5Y with 0.04 mag subtracted from the low-z sub-samples")
    run(f"C7 DES5Y shifted, IDSURVEY != 10 [{zc}]", [f"des5y_shift_survey_{zc}"], "des")
    run(f"C7 DES5Y shifted, zHD < 0.1 [{zc}]", [f"des5y_shift_z_{zc}"], "des")
    say("\n# C8: DESI DR1 BAO (alpha prior primary; variant beside), DES5Y prior box")
    for pr in ("primary", "variant"):
        if zc == "hd":
            run(f"C8 DESI DR1 alone, alpha {pr}", [f"bao_dr1_{pr}"], "des")
        run(f"C8 DESI DR1 + Pantheon+, alpha {pr} [{zc}]", [f"bao_dr1_{pr}", f"pantheon_{zc}"], "des")
        run(f"C8 DESI DR1 + DES5Y, alpha {pr} [{zc}]", [f"bao_dr1_{pr}", f"des5y_{zc}"], "des")
    say("\n# M1c: DES-Dovekie SN-only (data trace)")
    run(f"M1c DES-Dovekie [{zc}]", [f"dovekie_{zc}"], "des")
    say("\n# M2: DESI DR2 BAO (data trace), alpha primary")
    if zc == "hd":
        run("M2 DESI DR2 alone, alpha primary", ["bao_dr2_primary"], "des")
    for s in ("pantheon", "des5y", "dovekie"):
        run(f"M2 DESI DR2 + {s}, alpha primary [{zc}]", ["bao_dr2_primary", f"{s}_{zc}"], "des")
    if zc == "hd":
        for s in ("pantheon", "des5y", "dovekie"):
            run(f"M2 DESI DR2 + {s}, alpha variant [{zc}]", ["bao_dr2_variant", f"{s}_{zc}"], "des")

# C0.4 (part): grid LCDM minimum vs the minimiser (fits.json)
fits = json.load(open(os.path.join(OUT, "fits.json")))
say("\n# C0.4: grid minimum on the LCDM line (w0 = -1, wa = 0) vs the Brent minimiser")
j0, k0 = int(np.argmin(abs(W0 + 1))), int(np.argmin(abs(WA)))
for name, key in (("pantheon", "Pantheon+"), ("des5y", "DES5Y"), ("des5y_shift_survey", "DES5Y shifted (IDSURVEY != 10)"), ("dovekie", "M1 DES-Dovekie")):
    for zc in CONV:
        g = grid(f"{name}_{zc}")[:, j0, k0]
        f = fits[f"{key} [{zc}]"]
        say(f"   {name}_{zc}: grid LCDM min {g.min():.3f} at Om {OM[np.argmin(g)]:.3f}; minimiser {f['lcdm']['chi2']:.3f} at {f['lcdm']['Omega_m']:.4f}; diff {g.min() - f['lcdm']['chi2']:.3f}; "
            f"w0wa grid min (DES box) {np.where(mask(BOX['des']), grid(f'{name}_{zc}'), np.inf).min():.3f} vs minimiser {f['w0wa']['chi2']:.3f}")

json.dump(res, open(os.path.join(OUT, "posteriors.json"), "w"), indent=1)
open(os.path.join(OUT, "posteriors.out"), "w").write("\n".join(lines) + "\n")
