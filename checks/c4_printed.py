"""C4 from the printed numbers alone (conventions.json sigma_statements). Output: checks/out/c4_printed.out"""
import os
import numpy as np
from scipy.stats import chi2
from scipy.special import erfinv

import sys
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "replay-out") if "--replay" in sys.argv else os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
OUT = os.path.join(OUT, "c4_printed.out")
nsig = lambda d, k: float(np.sqrt(2) * erfinv(1 - chi2.sf(d, k)))
lines = [
    f"'~2.2 sigma' (l.280): DES5Y Omega_m 0.354 +- 0.017 vs Planck 0.3135 +- 0.0081: {(0.354 - 0.3135) / np.hypot(0.017, 0.0081):.2f} (sd only {(0.354 - 0.3135) / 0.017:.2f})",
    f"'roughly 3 sigma' (l.289): Delta chi2 12.5 (Planck - w0wa), k=3 {nsig(12.5, 3):.2f}, k=2 {nsig(12.5, 2):.2f}, k=1 {nsig(12.5, 1):.2f}",
    f"'~2 sigma' (l.290): Delta chi2 7.4 (LCDM - w0wa), k=2 {nsig(7.4, 2):.2f}, k=1 {nsig(7.4, 1):.2f}",
    f"C7 'no significant preference': Delta chi2 2.3 (legends) k=2 {nsig(2.3, 2):.2f}; 2.5 (text) k=2 {nsig(2.5, 2):.2f}",
]
print("\n".join(lines)); open(OUT, "w").write("\n".join(lines) + "\n")
