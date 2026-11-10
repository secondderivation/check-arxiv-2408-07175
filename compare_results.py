"""Release replay checks for case 2408.07175, distinct from the frozen scientific claim criteria (see VALIDATION.md)."""
import math
import re
from fractions import Fraction
from pathlib import Path

NUM = re.compile(r'(?<![\w.])[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?')


def text_compare(old, new, atol=1e-9, rtol=1e-9):
    if NUM.sub('#', old) != NUM.sub('#', new):
        return ['Output structure or verdict changed']
    a, b = NUM.findall(old), NUM.findall(new)
    if len(a) != len(b):
        return ['Numeric replay tolerance exceeded']
    # Compare the printed decimal values exactly, including the allowed boundary.
    # Binary floats can turn a difference of 0.001 into 0.001000000000000334.
    absolute, relative = Fraction(str(atol)), Fraction(str(rtol))
    for x, y in zip(a, b):
        if not (math.isfinite(float(x)) and math.isfinite(float(y))):
            return ['Numeric replay tolerance exceeded']
        x, y = Fraction(x), Fraction(y)
        if abs(x - y) > max(absolute, relative * max(abs(x), abs(y))):
            return ['Numeric replay tolerance exceeded']
    return []


# A level of 1.0000 converts to a sigma set by floating-point rounding (8.29 here); any value >= 5 is the same verdict.
HIGH_SIGMA = re.compile(r'\b(\d+\.\d+) sigma')


def _cap_sigma(s):
    return HIGH_SIGMA.sub(lambda m: '>=5 sigma' if float(m.group(1)) >= 5 else m.group(0), s)


def posterior_compare(old, new):
    """Grid posteriors from regenerated grids: printed values within 1e-3 (one unit in the last printed place of a
    level; PROTOCOL C0.4 calls 0.01 numerical uncertainty), words and verdicts unchanged."""
    return text_compare(_cap_sigma(old), _cap_sigma(new), atol=1e-3, rtol=0)


SPOT = re.compile(r'^spot check .*= ([0-9.eE+-]+)$')
SHAPE = re.compile(r'^grid \((\d+, \d+, \d+)\), (\d+) arrays')


def grids_compare(old, new):
    """grids.py prints progress lines and its run time; only the spot checks and the grid shape count.
    A spot check passes below 1e-4 in chi2 (the grids are re-computed, not copied)."""
    pick = lambda s, rx: [rx.match(x) for x in s.splitlines() if rx.match(x)]
    so, sn = pick(old, SPOT), pick(new, SPOT)
    go, gn = pick(old, SHAPE), pick(new, SHAPE)
    if len(sn) != len(so) or len(gn) != 1:
        return ['Grid replay output missing or incomplete']
    if gn[0].groups() != go[0].groups():
        return ['Grid shape or array count changed']
    if any(float(m.group(1)) >= 1e-4 for m in sn):
        return ['Grid spot check above 1e-4']
    return []


REF = Path(__file__).resolve().parent  # the published package, never the work copy the checks rewrite


def compare_run(run, work, stdout):
    """Compares against the published output in the package root; `work` is unused because the checks rewrite
    their own outputs in the work copy, which would make the comparison a self-comparison."""
    old = (REF / 'checks' / run['out']).read_text()
    if run['out'].endswith('grids.out'):
        return grids_compare(old, stdout)
    if run['out'].endswith('posteriors.out'):
        return posterior_compare(old, stdout)
    if run['out'].endswith('fits.out'):
        return text_compare(old, stdout, atol=1e-3, rtol=0)  # minimiser results, printed to 3-4 decimals
    return text_compare(old, stdout)
