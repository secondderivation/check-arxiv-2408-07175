"""C6 and M1a: magnitudes of SN common to Pantheon+ and DES5Y (and DES-Dovekie).
Delta m = m_b_corr(Pantheon+) - (MU - 19.33), conventions.json common_sn. Output: checks/out/crossmatch.out, .json
"""
import json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn

OUT = sn.OUT
lines, res = [], {}
def say(*a):
    s = " ".join(str(x) for x in a); print(s); lines.append(s)

P = sn.pantheon(cut=False).meta.copy()
P["CID"] = P.CID.astype(str)
NAMES = {10: "DES5Y", **sn.LOWZ_SURVEYS}
ORDER = ["DES5Y", "FOUND", "CFA3S", "CFA3K", "CFA4P2", "CFA4P3", "CSP"]
TABLE1 = json.load(open(os.path.join(sn.HERE, "..", "paper_values.json")))["table1"]


def match(des, rule):
    des = des.copy(); des["CID"] = des.CID.astype(str)
    if rule == "same_survey":
        j = des.merge(P[["CID", "IDSURVEY", "m_b_corr", "zHD"]], on=["CID", "IDSURVEY"], suffixes=("", "_P"))
    else:  # any survey, Pantheon+ duplicates averaged
        pm = P.groupby("CID").agg(m_b_corr=("m_b_corr", "mean"), zHD_P=("zHD", "mean"), nP=("m_b_corr", "size")).reset_index()
        j = des.merge(pm, on="CID")
    if "zHD_P" not in j:
        j = j.rename(columns={"zHD_P": "zHD_P"})
    j["dm"] = j.m_b_corr - (j.MU - sn.CALIB)
    j["group"] = j.IDSURVEY.map(NAMES)
    return j


def table(j, label, compare=True):
    say(f"\n== {label}: {len(j)} matched SN")
    say(f"  {'sample':10s} {'N':>4s} {'z_med(DES)':>10s} {'z_med(P+)':>9s} {'mean':>8s} {'error':>7s} {'disp':>6s}   printed (N, z_med, mean +- error)")
    out = {}
    groups = [(g, j[j.group == g]) for g in ORDER] + [("All low z", j[j.group != "DES5Y"])]
    for g, x in groups:
        if not len(x):
            continue
        n = len(x); mean = x.dm.mean(); sd = x.dm.std(ddof=1) if n > 1 else float("nan")
        err = sd / np.sqrt(n)
        zc = "zHD_P" if "zHD_P" in x else "zHD_P"
        row = {"N": n, "z_median": float(x.zHD.median()), "z_median_P": float(x[zc].median()), "mean": mean, "error": err, "dispersion": sd}
        out[g] = row
        pr = TABLE1.get(g)
        ps = f"({pr['N']}, {pr['z_median']}, {pr['mean']} +- {pr['error']})" if (pr and compare) else ""
        say(f"  {g:10s} {n:4d} {row['z_median']:10.3f} {row['z_median_P']:9.3f} {mean:8.4f} {err:7.4f} {sd:6.3f}   {ps}")
    diff = out["All low z"]["mean"] - out["DES5Y"]["mean"]
    derr = np.hypot(out["All low z"]["error"], out["DES5Y"]["error"])
    say(f"  differential (all low z) - (DES): {diff:.4f} +- {derr:.4f}  ({abs(diff) / derr:.1f} scatter errors from zero; no cross-catalogue covariance); printed: -0.0482 - (-0.0122) = -0.0360")
    out["differential"] = {"value": diff, "error": derr}
    return out


des = sn.des5y(cut=False).meta
say("C6 inputs: Pantheon+ entries", len(P), "; DES Y3 rows (IDSURVEY 10)", int((P.IDSURVEY == 10).sum()),
    "; Foundation rows (150)", int((P.IDSURVEY == 150).sum()), "; DES5Y SN", len(des))
fz = P[P.IDSURVEY == 150].zHD
say("  Foundation rows (printed '150 SN from the Foundation Survey', l.389): all", len(fz), "; 0.02 <= zHD <= 1.2:", int(((fz >= 0.02) & (fz <= 1.2)).sum()),
    "; zHD > 0.01:", int((fz > 0.01).sum()), "; unique CIDs:", P[P.IDSURVEY == 150].CID.nunique())
fc = P[P.IDSURVEY == 150].zCMB   # post_hoc diagnostic: the same cut on zCMB
say("  post_hoc diagnostic: Foundation rows with 0.02 <= zCMB <= 1.2:", int(((fc >= 0.02) & (fc <= 1.2)).sum()),
    "; whole Pantheon+ with that zCMB cut:", int(((P.zCMB >= 0.02) & (P.zCMB <= 1.2)).sum()), "(the paper's 1417 needs zHD)")
dupe = P[P.IDSURVEY.isin([10, 150, 5, 63, 64, 65, 66])].groupby(["CID", "IDSURVEY"]).size()
say("  Pantheon+ (CID, IDSURVEY) pairs with more than one entry in these surveys:", int((dupe > 1).sum()))
for rule in ("same_survey", "any_survey"):
    j = match(des, rule)
    res[f"DES5Y {rule}"] = table(j, f"DES5Y, rule {rule}")
    if rule == "same_survey":
        j.to_csv(os.path.join(OUT, "crossmatch_des5y.csv"), index=False, float_format="%.5f",
                 columns=["CID", "IDSURVEY", "group", "zHD", "zHD_P", "MU", "MUERR_FINAL", "m_b_corr", "dm"])
        lowz = des[des.IDSURVEY != 10]
        miss = sorted(set(lowz.CID.astype(str)) - set(j.CID))
        say("  DES5Y low-z SN not matched (same survey):", len(miss), miss)
        for c in miss:
            say(f"    {c}: in Pantheon+ as IDSURVEY {sorted(P[P.CID == c].IDSURVEY.tolist())}")

# weighting variant: inverse-variance with the catalogue diagonal errors (reported only)
j = match(des, "same_survey")
pe = sn.pantheon(cut=False).meta.copy(); pe["CID"] = pe.CID.astype(str)
j = j.merge(pe[["CID", "IDSURVEY", "m_b_corr_err_DIAG"]], on=["CID", "IDSURVEY"])
w = 1 / (j.MUERR_FINAL ** 2 + j.m_b_corr_err_DIAG ** 2)
for g, x in [("DES5Y", j.group == "DES5Y"), ("All low z", j.group != "DES5Y")]:
    say(f"  variant, inverse-variance weights (diagonal errors added): {g} mean {np.average(j.dm[x], weights=w[x]):.4f}")

# M1a: DES-Dovekie (data trace, update role)
dov = sn.dovekie(cut=False).meta
res["Dovekie same_survey"] = table(match(dov, "same_survey"), "M1a DES-Dovekie, rule same_survey", compare=False)
res["Dovekie any_survey"] = table(match(dov, "any_survey"), "M1a DES-Dovekie, rule any_survey", compare=False)

json.dump(res, open(os.path.join(OUT, "crossmatch.json"), "w"), indent=1, default=float)
open(os.path.join(OUT, "crossmatch.out"), "w").write("\n".join(lines) + "\n")
