"""PROTOCOL section 2: controls C0.1-C0.4. Output: checks/out/c0_controls.out (G14)."""
import hashlib, io, os, sys, urllib.request
import numpy as np
import pandas as pd
from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
os.makedirs(OUT, exist_ok=True)
lines = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); lines.append(s)

POINTS = [(0.3135, -1.0, 0.0), (0.35, -0.8, -1.0), (0.45, -0.6, -5.0)]

# C0.1 data identity: Cobaya DES5Y file vs the DES release (tag 1.3, commit e3493cb3)
META_URL = "https://raw.githubusercontent.com/des-science/DES-SN5YR/e3493cb3b9505fc1f3f392d887364b50dae21439/4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv"
META_SHA = "80e081d13145a7ee66730b59826ed83e6414087fd1e26a5a206a196f650f02c5"
mp = os.path.join(sn.RAW, "desrepo/DES-SN5YR_HD+MetaData.csv")
if not os.path.exists(mp):
    urllib.request.urlretrieve(META_URL, mp)
assert hashlib.sha256(open(mp, "rb").read()).hexdigest() == META_SHA
meta = pd.read_csv(mp)
cob = pd.read_csv(os.path.join(sn.RAW, "DESY5/DES-SN5YR_HD.csv"))
j = cob.merge(meta, on="CID", suffixes=("", "_rel"))
say("C0.1 DES5Y Cobaya file vs DES release 1.3:", len(cob), "rows,", len(j), "matched by CID")
for c in ["MU", "MUERR_FINAL", "zHD", "zHEL", "IDSURVEY"]:
    say(f"  max |{c} - release| = {np.abs(j[c] - j[c + '_rel']).max():.2e}")
say("  same order as the release:", bool((cob.CID.astype(str).values == meta.CID.astype(str).values).all()))

# C0.2 radiation and redshift convention
P, D = sn.pantheon(), sn.des5y()
OMR = 4.18e-5 / 0.6744 ** 2
say("C0.2 Omega_r =", f"{OMR:.3e}")
for s in (P, D):
    for pt in POINTS:
        c0 = s.chi2(*pt)[0][0]
        cr = s.chi2(*pt, omr=OMR)[0][0]
        cz = s.chi2(*pt, zconv="hd")[0][0]
        say(f"  {s.name:10s} {pt}: chi2 {c0:.3f}; with radiation {cr - c0:+.4f}; (1+zHD) instead of (1+zHEL) {cz - c0:+.4f}")

# C0.3 Cobaya 3.6.2 likelihoods with our background (minimal Theory)
from cobaya.theory import Theory
from cobaya.model import get_model

class Bg(Theory):
    params = {"om": None, "w0": None, "wa": None, "H0": None, "rdrag": None}
    def get_requirements(self):
        return {}
    def must_provide(self, **req):
        return {}
    def get_can_provide(self):
        return ["angular_diameter_distance", "comoving_radial_distance", "Hubble", "rdrag"]
    def calculate(self, state, want_derived=True, **p):
        state["p"] = dict(p)
    def _dc(self, z):
        p = self.current_state["p"]
        z = np.atleast_1d(z)
        tab = sn.dc_table(p["om"], p["w0"], p["wa"])
        return sn.C_KMS / p["H0"] * sn.interp_rows(tab, z)[0]
    def get_angular_diameter_distance(self, z):
        return self._dc(z) / (1 + np.atleast_1d(z))
    def get_comoving_radial_distance(self, z):
        return self._dc(z)
    def get_Hubble(self, z, units="km/s/Mpc"):
        p = self.current_state["p"]
        h = p["H0"] * sn.efunc(np.atleast_1d(z), p["om"], p["w0"], p["wa"])
        return h / sn.C_KMS if units == "1/Mpc" else h
    def get_param(self, name):
        return self.current_state["p"][name]

pk = os.path.join(sn.RAW, "_cobaya")
os.makedirs(os.path.join(pk, "data"), exist_ok=True)
for name, src in (("sn_data", sn.RAW), ("bao_data", os.path.join(sn.RAW, "bao_data"))):
    link = os.path.join(pk, "data", name)
    if not os.path.exists(link):
        os.symlink(os.path.abspath(src), link)
    open(os.path.join(src, "version.dat"), "w").write({"sn_data": "v1.8", "bao_data": "v2.6"}[name])
like = {"sn.pantheonplus": None, "sn.desy5": None, "bao.desi_2024_bao_all": None}
info = {"likelihood": like, "theory": {"bg": Bg},
        "params": {k: {"prior": {"min": -20, "max": 200}} for k in ("om", "w0", "wa", "H0", "rdrag")}, "packages_path": pk}
model = get_model(info)
Pfull = sn.pantheon(cut=False)
keep = (Pfull.meta.zHD > 0.01).to_numpy()
with open(os.path.join(sn.RAW, "PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov")) as f:
    n = int(f.readline()); covP = np.loadtxt(f).reshape(n, n)
P01 = sn.Sample("P+ z>0.01", Pfull.m[keep], Pfull.zhd[keep], Pfull.zhel[keep], cov=covP[np.ix_(keep, keep)])
Dfull = sn.des5y(cut=False)
bao = sn.BAO("dr1")
worst = 0.0
say("C0.3 our chi2 minus Cobaya 3.6.2 (-2 loglike), same samples as Cobaya cuts them:")
for pt in POINTS:
    lp = model.loglikes({"om": pt[0], "w0": pt[1], "wa": pt[2], "H0": 67.44, "rdrag": 147.31}, as_dict=True)[0]
    ours_p = P01.chi2(*pt)[0][0]                     # Cobaya projects the offset out: profiled chi2
    ours_d = Dfull.chi2(*pt)[0][0]                   # Cobaya's DESy5 inherits the Pantheon+ projection: profiled too
    alpha = sn.C_KMS / (67.44 * 147.31)
    ours_b = bao.chi2_alpha(bao.shapes(*pt), [alpha])[0, 0]
    dp, dd, db = ours_p + 2 * lp["sn.pantheonplus"], ours_d + 2 * lp["sn.desy5"], ours_b + 2 * lp["bao.desi_2024_bao_all"]
    worst = max(worst, abs(dp), abs(dd), abs(db))
    say(f"  {pt}: Pantheon+ {dp:+.5f}  DES5Y {dd:+.5f}  DESI DR1 {db:+.5f}")
say("  C0.3", "pass" if worst < 0.01 else "FAIL", f"(max |difference| {worst:.2e}, tolerance 0.01)")

# C0.4 numerics: distance integral vs quad
worst = 0.0
for pt in POINTS + [(0.25, -1.1, 1.0), (0.8, 0.9, -9.9)]:
    tab = sn.dc_table(*pt)
    for z in (0.02, 0.3, 1.2, 2.33):
        q = quad(lambda x: 1 / sn.efunc(x, *pt), 0, z, epsabs=0, epsrel=1e-12)[0]
        worst = max(worst, abs(sn.interp_rows(tab, np.array([z]))[0, 0] / q - 1))
say("C0.4 distance table vs quad: max relative error", f"{worst:.1e}", "pass" if worst < 1e-6 else "FAIL")

open(os.path.join(OUT, "c0_controls.out"), "w").write("\n".join(lines) + "\n")
