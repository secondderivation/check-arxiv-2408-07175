"""post_hoc diagnostic (G5): which redshift column gives the printed Table 1 z_median? Output: checks/out/c6_zmedian.out"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sn
OUT = sn.OUT
j = pd.read_csv(os.path.join(OUT, "crossmatch_des5y.csv")); j["CID"] = j.CID.astype(str)
d = pd.read_csv(os.path.join(sn.RAW, "DESY5/DES-SN5YR_HD.csv")); d["CID"] = d.CID.astype(str)
p = sn.pantheon(cut=False).meta; p["CID"] = p.CID.astype(str)
j = j.merge(d[["CID", "zCMB", "zHEL"]], on="CID").merge(p[["CID", "IDSURVEY", "zCMB", "zHEL"]], on=["CID", "IDSURVEY"], suffixes=("", "_P"))
printed = {"DES5Y": 0.329, "FOUND": 0.043, "CFA3S": 0.037, "CFA3K": 0.033, "CFA4P2": 0.038, "CFA4P3": 0.033, "CSP": 0.038}
lines = ["median redshift of the matched SN per sample, by column (DES5Y file: zHD zCMB zHEL; Pantheon+: _P), printed Table 1 value last"]
for g, pv in printed.items():
    x = j[j.group == g]
    lines.append(f"{g:7s} N={len(x):3d} " + " ".join(f"{c}={x[c].median():.4f}" for c in ["zHD", "zCMB", "zHEL", "zHD_P", "zCMB_P", "zHEL_P"])
                 + f"  lower-median zHD={np.sort(x.zHD)[(len(x) - 1) // 2]:.4f}  printed {pv}")
x = j[j.group != "DES5Y"]
lines.append("All low z " + " ".join(f"{c}={x[c].median():.4f}" for c in ["zHD", "zCMB", "zHEL"]) + "  printed 0.038")
print("\n".join(lines)); open(os.path.join(OUT, "c6_zmedian.out"), "w").write("\n".join(lines) + "\n")
