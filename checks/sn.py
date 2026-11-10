"""SN and BAO likelihoods for case 2408-07175, written from the paper and conventions.json (G13).

Data come from data/fetch_data.py (RAW, default papers/2408-07175/raw). Magnitudes:
Pantheon+ m_b_corr; DES5Y MU - 19.33 (Eq. 2); DES-Dovekie MU - 19.33. The offset M' is profiled
(chi2) or marginalised analytically (posterior); the two differ by a constant.
"""
import os
import sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.environ.get("SD_RAW", os.path.join(HERE, "..", "raw"))
if "SD_RAW" not in os.environ and not os.path.isdir(RAW):   # published package: fetch_data.py writes data/raw/
    RAW = os.path.join(HERE, "..", "data", "raw")
C_KMS = 299792.458
CALIB = 19.33                       # Eq. 2, l.240
PLANCK_OM = 0.3135                  # l.259
# Outputs: checks/out/, or with --replay (reproduce.py in the published package) a separate folder beside the
# work copy, so a replay never rewrites the published files it is compared against (VALIDATION.md).
REPLAY = "--replay" in sys.argv
OUT = os.path.join(HERE, "..", "..", "replay-out") if REPLAY else os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
LOWZ_SURVEYS = {150: "FOUND", 63: "CFA3S", 64: "CFA3K", 65: "CFA4P2", 66: "CFA4P3", 5: "CSP"}
ZGRID = np.linspace(0.0, 2.4, 4001)


# ---------------------------------------------------------------- background
def efunc(z, om, w0=-1.0, wa=0.0, omr=0.0):
    """E(z) for flat w0waCDM (Eq. 1); om, w0, wa may be arrays broadcast against z."""
    a1 = 1.0 + z
    fde = a1 ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * z / a1)
    return np.sqrt(om * a1 ** 3 + omr * a1 ** 4 + (1.0 - om - omr) * fde)


def dc_table(om, w0, wa, omr=0.0):
    """Dimensionless comoving distance H0 D_C / c on ZGRID for parameter arrays of shape (P,),
    with its derivative 1/E: returns (D, D'), both (P, len(ZGRID)). Cumulative trapezoid with the
    end correction (h^2/12)(f'(a) - f'(b)), f' by finite differences, so the error is O(h^4)."""
    om, w0, wa = (np.atleast_1d(x)[:, None] for x in (om, w0, wa))
    inv = 1.0 / efunc(ZGRID[None, :], om, w0, wa, omr)
    dz = ZGRID[1] - ZGRID[0]
    out = np.zeros_like(inv)
    out[:, 1:] = np.cumsum(0.5 * (inv[:, 1:] + inv[:, :-1]) * dz, axis=1)
    der = np.gradient(inv, dz, axis=1, edge_order=2)
    out[:, 1:] -= dz ** 2 / 12.0 * (der[:, 1:] - der[:, :1])
    return out, inv


def interp_rows(table, z):
    """Cubic Hermite interpolation of each row of a dc_table (values and derivatives) at redshifts z
    (error O(h^4), checked against quad in C0.4)."""
    val, der = table
    dz = ZGRID[1] - ZGRID[0]
    i = np.clip((z / dz).astype(int), 0, len(ZGRID) - 2)
    t = (z - ZGRID[i]) / dz
    h00, h10, h01, h11 = 2 * t**3 - 3 * t**2 + 1, t**3 - 2 * t**2 + t, -2 * t**3 + 3 * t**2, t**3 - t**2
    return val[:, i] * h00 + der[:, i] * dz * h10 + val[:, i + 1] * h01 + der[:, i + 1] * dz * h11


# ---------------------------------------------------------------- SN samples
class Sample:
    """m (data), zhd, zhel, cov (or icov), with the chi2 profiled / marginalised over the offset."""

    def __init__(self, name, m, zhd, zhel, cov=None, icov=None, meta=None):
        self.name, self.m, self.zhd, self.zhel, self.meta = name, np.asarray(m, float), np.asarray(zhd, float), np.asarray(zhel, float), meta
        self.icov = np.linalg.inv(cov) if icov is None else icov
        self.icov = 0.5 * (self.icov + self.icov.T)
        self.e = self.icov.sum(axis=0)          # C^-1 1
        self.E = self.e.sum()                   # 1' C^-1 1
        self.n = len(self.m)

    def mth(self, om, w0=-1.0, wa=0.0, omr=0.0, zconv="hel"):
        dc = interp_rows(dc_table(om, w0, wa, omr), self.zhd)
        fac = (1.0 + self.zhel) if zconv == "hel" else (1.0 + self.zhd)
        return 5.0 * np.log10(fac[None, :] * dc)

    def chi2_parts(self, mth):
        """mth: (P, N). Returns A = r'C^-1 r, B = 1'C^-1 r for r = m - mth."""
        r = self.m[None, :] - mth
        Cr = r @ self.icov
        return np.einsum("ij,ij->i", Cr, r), r @ self.e

    zconv = "hel"   # default D_L convention of the instance: "hel" (1+z_HEL) or "hd" (1+z_HD)

    def chi2(self, om, w0=-1.0, wa=0.0, omr=0.0, zconv=None):
        """chi2 profiled over the offset M'; also returns the best offset (m - mth = M')."""
        zconv = zconv or self.zconv
        A, B = self.chi2_parts(self.mth(om, w0, wa, omr, zconv))
        return A - B ** 2 / self.E, B / self.E

    def chi2_marg(self, om, w0=-1.0, wa=0.0):
        """-2 ln L with M' marginalised (flat prior): profiled chi2 + ln(E/2pi) (constant)."""
        return self.chi2(om, w0, wa)[0] + np.log(self.E / (2 * np.pi))


def pantheon(cut=True):
    d = pd.read_csv(os.path.join(RAW, "PantheonPlus/Pantheon+SH0ES.dat"), sep=r"\s+")
    with open(os.path.join(RAW, "PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov")) as f:
        n = int(f.readline())
        cov = np.loadtxt(f).reshape(n, n)
    keep = ((d.zHD >= 0.02) & (d.zHD <= 1.2)).to_numpy() if cut else np.ones(len(d), bool)
    dd = d[keep].reset_index(drop=True)
    return Sample("Pantheon+", dd.m_b_corr, dd.zHD, dd.zHEL, cov=cov[np.ix_(keep, keep)], meta=dd)


def des5y(cut=True, shift=0.0, shift_rule="survey"):
    d = pd.read_csv(os.path.join(RAW, "DESY5/DES-SN5YR_HD.csv"))
    with open(os.path.join(RAW, "DESY5/covsys_000.txt")) as f:
        n = int(f.readline())
        cov = np.loadtxt(f).reshape(n, n)
    cov = cov + np.diag(d.MUERR_FINAL.to_numpy() ** 2)
    keep = np.ones(len(d), bool)
    if cut:
        keep = ((d.zHD >= 0.02) & (d.zHD <= 1.2) & (d.MUERR_FINAL <= 1.0)).to_numpy()
    dd = d[keep].reset_index(drop=True)
    m = dd.MU.to_numpy() - CALIB
    if shift:
        low = (dd.IDSURVEY != 10).to_numpy() if shift_rule == "survey" else (dd.zHD < 0.1).to_numpy()
        m = m - shift * low
    return Sample("DES5Y" + (f" shifted {shift}" if shift else ""), m, dd.zHD, dd.zHEL, cov=cov[np.ix_(keep, keep)], meta=dd)


def dovekie(cut=True):
    p = os.path.join(RAW, "desrepo/DES-Dovekie_HD.csv")
    rows = [l.split()[1:] for l in open(p) if l.startswith("SN:")]
    cols = [l.split()[1:] for l in open(p) if l.startswith("VARNAMES:")][0]
    d = pd.DataFrame(rows, columns=cols)
    for c in cols[1:]:
        d[c] = pd.to_numeric(d[c])
    z = np.load(os.path.join(RAW, "desrepo/STAT+SYS.npz"))
    n = int(z["nsn"][0])
    icov = np.zeros((n, n))
    icov[np.triu_indices(n)] = z["cov"]
    icov = icov + np.triu(icov, 1).T
    keep = ((d.zHD >= 0.02) & (d.zHD <= 1.2)).to_numpy() if cut else np.ones(len(d), bool)
    if not keep.all():  # sub-sample: invert to the covariance, take the block, invert back
        cov = np.linalg.inv(icov)[np.ix_(keep, keep)]
        icov = np.linalg.inv(cov)
    dd = d[keep].reset_index(drop=True)
    return Sample("DES-Dovekie", dd.MU.to_numpy() - CALIB, dd.zHD, dd.zHEL, icov=icov, meta=dd)


# ---------------------------------------------------------------- BAO
class BAO:
    """DESI Gaussian BAO likelihood; observables = alpha * dimensionless distance, alpha = c/(H0 r_d)."""

    def __init__(self, release="dr1"):
        base = os.path.join(RAW, "bao_data", "desi_2024_gaussian_bao_ALL_GCcomb" if release == "dr1" else "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb")
        rows = [l.split() for l in open(base + "_mean.txt") if not l.startswith("#") and l.strip()]
        self.z = np.array([float(r[0]) for r in rows])
        self.val = np.array([float(r[1]) for r in rows])
        self.kind = [r[2] for r in rows]
        self.icov = np.linalg.inv(np.loadtxt(base + "_cov.txt"))
        self.release = release

    def shapes(self, om, w0, wa, omr=0.0):
        """Dimensionless observables (alpha = 1) for parameter arrays: (P, n)."""
        tab = dc_table(om, w0, wa, omr)
        dm = interp_rows(tab, self.z)
        om_, w0_, wa_ = (np.atleast_1d(x)[:, None] for x in (om, w0, wa))
        dh = 1.0 / efunc(self.z[None, :], om_, w0_, wa_, omr)
        dv = (self.z[None, :] * dm ** 2 * dh) ** (1.0 / 3.0)
        out = np.empty_like(dm)
        for j, k in enumerate(self.kind):
            out[:, j] = {"DM_over_rs": dm, "DH_over_rs": dh, "DV_over_rs": dv}[k][:, j]
        return out

    def chi2_alpha(self, shapes, alphas):
        """chi2 for every parameter row and alpha node: (P, K)."""
        S = shapes
        a = S @ self.icov @ self.val            # s'C^-1 d
        b = np.einsum("ij,jk,ik->i", S, self.icov, S)  # s'C^-1 s
        c = self.val @ self.icov @ self.val
        al = np.asarray(alphas)[None, :]
        return c - 2 * al * a[:, None] + al ** 2 * b[:, None]


R_D, SIG_RD = 147.31, 0.31          # l.650
H0_P, SIG_H0 = 67.44, 0.58          # l.303


def alpha_nodes(prior="primary", n=241):
    """ln-alpha quadrature nodes and log-weights (per unit ln alpha) for the two alpha priors."""
    if prior == "primary":   # H0 flat on [20, 100], r_d Gaussian: alpha density ~ alpha^-2, i.e. alpha^-1 per ln alpha
        lo, hi = np.log(C_KMS / (100 * R_D)), np.log(C_KMS / (20 * R_D))
        la = np.linspace(lo, hi, n)
        lw = -la
    else:                    # plus the Planck H0 prior: alpha log-normal
        mu = np.log(C_KMS / (H0_P * R_D))
        s = np.hypot(SIG_H0 / H0_P, SIG_RD / R_D)
        la = np.linspace(mu - 6 * s, mu + 6 * s, n)
        lw = -0.5 * ((la - mu) / s) ** 2
    lw = lw + np.log(la[1] - la[0])
    return np.exp(la), lw


def bao_m2lnL(bao, shapes, prior="primary"):
    """-2 ln L_BAO marginalised over alpha (up to a constant)."""
    al, lw = alpha_nodes(prior)
    x = -0.5 * bao.chi2_alpha(shapes, al) + lw[None, :]
    mx = x.max(axis=1, keepdims=True)
    return -2 * (mx[:, 0] + np.log(np.exp(x - mx).sum(axis=1)))


def bao_chi2_min(bao, shapes):
    """chi2 minimised over alpha (no prior), for reporting."""
    a = shapes @ bao.icov @ bao.val
    b = np.einsum("ij,jk,ik->i", shapes, bao.icov, shapes)
    c = bao.val @ bao.icov @ bao.val
    return c - a ** 2 / b, a / b
