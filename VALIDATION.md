# Release validation contract

Dated 15 October 2026 during package preparation. These are engineering rules for replaying the package. They are not part of the preregistration and do not replace the scientific tolerances in `PROTOCOL.md`.

## Three separate outcomes

1. **Scientific reproduction.** This uses the frozen claim criteria in `PROTOCOL.md` with addenda A1 and A2. The final adjudications after review are the claim outcomes in the README, and a replay PASS changes none of them.
2. **Numerical replay.** This re-runs what can be re-run from the package and compares the result with the saved outputs, under the rules below.
3. **Distribution integrity.** Every distributed file is checked against `MANIFEST.txt` (SHA-256), and downloaded data must match their pinned hashes. This verifies consistency, not scientific correctness.

## What is replayed, and what is not

Each check runs with `--replay`, which sends everything it writes (outputs, `fits.json` for the posteriors, `GRIDS.sha256`) to `_runs/<time>/replay-out/` instead of `checks/out/`, so the work copy stays identical to the manifest. Each run's stdout is compared with the published output in the package root (`compare_results.py`). `test_validation.py` tests this end to end.

- **Deterministic checks** (`c4_printed`, `crossmatch`, `c6_zmedian`). The text structure and verdict words must match exactly, with every number within 1e-9.
- **Best fits** (`fits`). Minimisations from fixed starts. The words must match exactly, with numbers within 1e-3. Values are printed to 3–4 decimals, so this allows a different last digit across platforms and nothing more.
- **Likelihood grids** (`grids`, slow, about 70 minutes on 4 cores).
  - The 14 grids (161 × 151 × 241 points, 0.66 GB) are rebuilt, not downloaded.
  - Only the spot checks and the grid shape are compared. Every spot check against a direct χ² evaluation must stay below 1e-4. The progress lines and run time go to stdout only and are ignored.
  - `checks/GRIDS.sha256` holds the hashes of our grids. A rebuild on another platform may differ in the last bits, so a hash mismatch alone is not a failure. On our machine the rebuilt grids were bit-identical.
- **Posteriors** (`posteriors`, slow, runs on the rebuilt grids). The words and verdicts must match, with every printed value within 1e-3; PROTOCOL C0.4 counts 0.01 as numerical uncertainty. A level of 1.0000 converts to a σ set by floating-point rounding, so any σ of 5 or more counts as the same value.
- **Not replayed**:
  - `c0_controls` needs Cobaya and the DES metadata file.
  - `c9_citations` needs the TeX sources of the cited papers, which the package does not redistribute. Set `SD_SRC` to a folder holding the arXiv e-prints 2408.07175v3, 1807.06209v4, 2404.03002v3, 2501.06664v1, 1910.00483v3 and 2202.04077v2, unpacked as `<id>/src-v<n>/`, to run it; the script checks each source's SHA-256.
  - `grid_convergence` takes a further 24 minutes.

  Their saved outputs are the record.

The suite is complete only with `python3 reproduce.py --all`. No CI step publishes, creates a release or records owner approval.

## Engineering correction — 25 October 2026

The verifier found that binary floating-point subtraction could reject an allowed printed difference of exactly 0.001 (for example, -7.995 versus -7.996). Numeric text comparisons now use exact rational values of the printed decimals and the same absolute/relative tolerances. Non-finite numeric values remain rejected. Regression tests accept the exact boundary and reject values above it, including through the best-fit replay path. No scientific criterion, replay tolerance, reference output or input data changed.

The README now distinguishes best-fit-point agreement from Planck-point agreement, corrects the count of affected Table 1 entries, and clarifies the error convention, figure-legend subtraction, CMB dataset label and alternate-prior gate. The manual citation-replay source list now includes the required Pantheon+ e-print 2202.04077v2. These are release-verification corrections, not changes to the frozen scientific protocol.
