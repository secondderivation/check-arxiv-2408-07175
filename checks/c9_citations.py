"""C9: numbers the paper attributes to other papers, checked against the cited papers' TeX.

Sources (arXiv e-prints, outside Git, hashes below) are read from SD_SRC (default ~/work/src,
layout <arxiv id>/src-v<N>/<file>.tex). Each check asserts that the quoted string is on the given
line of the cited source and of the paper, then prints both lines. Output: checks/out/c9_citations.out
"""
import hashlib, os

SRC = os.environ.get("SD_SRC", os.path.expanduser("~/work/src"))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "c9_citations.out")
FILES = {  # file -> sha256 of the e-print file used
    "paper": ("2408.07175/src-v3/SNpaper.tex", "b52743a408aae9a954e372200d0086aff14612beaf23103ab05249b9d0a1c583"),
    "planck18": ("1807.06209/src-v4/ms.tex", "fe0833d2d0fda17427a709c9559dc7b3a988b955f47bb7de7a02c900ffdc21cf"),
    "desi24": ("2404.03002/src-v3/JournalReview_adds_v2_kp7_paper1.tex", "88d7ca6ea395ecc55b0e4264bd3be6dd26bdb0c209a45ee3582560163b67168c"),
    "v25": ("2501.06664/src-v1/main.tex", "45c7a8c684180ae098b0c8c6e4ca04a4204c455030a0204e442ddd1708994e0c"),
    "v25app": ("2501.06664/src-v1/appendix.tex", "07b70ee3798b50b5cdb046de488bbd37dd2d338b39651dc266e813c3186f6919"),
    "eg21": ("1910.00483/src-v3/LCDM_compare_12.5.tex", "3e2d27c966cbeb6293fc61a786ea916d8d37a35e584462ee3f7bdcc73dd30e83"),
    "brout22": ("2202.04077/src-v2/main.tex", "5bd398947f9ec0355facec937e5865c80154455c2a29331cfd5355f6f8d06408"),
}
lines = []
def say(s=""):
    print(s); lines.append(s)

TEXT = {}
for k, (p, sha) in FILES.items():
    raw = open(os.path.join(SRC, p), "rb").read()
    got = hashlib.sha256(raw).hexdigest()
    assert got == sha, (p, got)
    TEXT[k] = raw.decode("utf-8", "replace").splitlines()
    say(f"source {k:9s} {p}  sha256 {sha[:16]}...")

def at(src, line, needle):
    s = TEXT[src][line - 1]
    assert needle in s, (src, line, needle)
    i = s.index(needle)
    return s[max(0, i - 60): i + len(needle) + 60].strip()

# (label, paper line, paper string, cited source, cited line, cited string, verdict)
CHECKS = [
    ("Planck 2018 w (Planck+lensing+SNe+BAO)", 133, r"$w=-1.028 \pm 0.031$", "planck18", 4258, r"w_0 = -1.028\pm0.031", "match"),
    ("DESI24 w, DESI+CMB+Pantheon+", 139, r"$w = - 0.997 \pm", "desi24", 670, r"DESI+CMB+Panth. & $0.3095\pm 0.0069$ & $67.74\pm 0.71$ &---& $-0.997\pm 0.025$", "match"),
    ("DESI24 2.5/3.5/3.9 sigma (Pantheon+, Union3, DES5Y)", 167, r"$\sim 2.5\sigma$, $3.5\sigma$ and $3.9 \sigma$", "desi24", 907, r"$2.5\sigma$, $3.5\sigma$ and $3.9\sigma$", "match"),
    ("DESI24 2.2/3.1/3.5 sigma (DESI+SDSS combined BAO, appendix)", 168, r"$\sim 2.2 \sigma$, $3.1\sigma$ and $3.5\sigma$", "desi24", 1284, r"$2.2\sigma$, $3.1\sigma$ and $3.5\sigma$", "match"),
    ("V25 total offset -0.042", 471, r"$\Delta \mu_{\rm offset} = -0.042$", "v25", 378, r"$-$0.042", "match (V25 defines it for Foundation vs DES common SN, l.370; the paper says 'low redshift SN')"),
    ("V25 intrinsic scatter 0.008", 477, r"$-0.008$ mag.", "v25", 380, r"\textbf{0.008} &\textbf{$-$0.034}", "match (sign convention: V25 lists the change on reverting)"),
    ("V25 host mass 0.010", 478, r"$-0.010$  mag.", "v25", 381, r"\textbf{0.010} &\textbf{$-$0.024}", "match"),
    ("V25 selection 0.016, remaining -0.008", 475, r"$-0.016$ mag.", "v25", 383, r"\textbf{0.016} & \textbf{$-$0.008}", "match"),
    ("V25 appendix BBC update ~0.01 low vs high z", 479, r"$\sim 0.01$ mag", "v25app", 71, r"${\sim} {+0.005}$~mag at low-$z$", "match (+0.005 low z, -0.005 high z, l.71-72)"),
    ("Planck EG21 TTTEEE Omega_m (the 'Planck point')", 259, r"$\Omega_m = 0.3135 \pm 0.0081", "eg21", 15, r"$0.3135 \pm 0.0081$", "match"),
    ("Planck EG21 TTTEEE H0 (prior in Fig. 1 fits)", 304, r"$67.44 \pm 0.58", "eg21", 13, r"$67.44\pm 0.58$", "match"),
    ("Brout22 Pantheon+ LCDM Omega_m ('in agreement')", 266, r"$\Omega_m = 0.333 \pm 0.018$", "brout22", 218, r"\def\fLCDMOMPan{$0.334\pm0.018$}", "agrees (0.333 vs 0.334)"),
    ("DESI24 DESI+CMB w0wa Omega_m", 657, r"$\Omega_m = 0.0344 ^{+0.032}_{-0.027}$", "desi24", 677, r"DESI+CMB & $0.344^{+0.032}_{-0.027}$", "TYPO: paper prints 0.0344, DESI24 gives 0.344 (S0)"),
]
say()
for lab, pl, ps, src, cl, cs, verdict in CHECKS:
    say(f"== {lab}: {verdict}")
    say(f"   paper l.{pl}: ...{at('paper', pl, ps)}...")
    say(f"   {src} l.{cl}: ...{at(src, cl, cs)}...")
# the paper's sum over V25's Table 1 and its 'remaining -0.026'
say()
parts = {"selection": -0.016, "scatter": -0.008, "host mass": -0.010, "unexplained": -0.008}
say(f"sum of the paper's V25 parts {parts}: {sum(parts.values()):.3f} (paper and V25: -0.042); -0.042 - (-0.016) = {-0.042 + 0.016:.3f} (paper l.477: -0.026)")
assert abs(sum(parts.values()) + 0.042) < 1e-9
open(OUT, "w").write("\n".join(lines) + "\n")
