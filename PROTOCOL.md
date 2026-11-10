# PROTOCOL: arXiv:2408.07175v3 (G. Efstathiou, "Evolving Dark Energy or Supernovae Systematics?", MNRAS 538 (2025) 875)

Public edition prepared 25 October 2026; internal role descriptions omitted. The original is retained privately. Scientific criteria, numerical thresholds and dated amendments are preserved.

Originally committed 2026-10-15 before any computation on the paper's data (G2). Changes go in a dated addendum at the end, also before the run they affect. Anything not planned here is labelled `post_hoc` in the report.

Inputs: fixed paper version and source hashes, `claims/C1–C9-r1.json`, `conventions.json`, `paper_values.json` (printed numbers, parsed from the TeX by `checks/v1_transcribe.py`; figure legends read from renderings), `datasets.json` and `data/fetch_data.py` (pinned files and hashes).


Information available before the original computation: the paper's TeX (all of it), renderings of its figures (Figs 1, 2, 3, 4, 5, 7; legends of Figs 1 and 4 transcribed), the arXiv abstract page, the INSPIRE publication record. From the data, only column headers were read and rows counted to plan the cross-match: the DES5Y file has 1,635 rows with IDSURVEY 10 and 194 others (118, 31, 19, 15, 8, 3 in surveys 150, 64, 65, 63, 5, 66); 190 of the 194 low-z CIDs and 145 DES CIDs occur in the Pantheon+ file by name; the Pantheon+ file has 203 rows with IDSURVEY 10. No magnitude, fit, χ² or posterior has been computed. The preparation did not include reading V25 (2501.06664), DESI24 (2404.03002) or the Dovekie papers beyond their titles. Prior knowledge included that DES later recalibrated its sample (Dovekie) and that the reanalysis lowered DES's evidence for evolving dark energy; these are not used for any claim outcome.

## 1. Implementation (fixed now)

**Code.** Own numpy/scipy code written from the paper's equations and conventions.json (`checks/sn.py`: data loading, cuts, covariances, distances, χ² profiled and marginalised over the offset; `checks/bao.py`: DESI Gaussian BAO likelihood). No author code is public or used (G13). Distances: D_C by cumulative trapezoid in z on 4,000 points to z = 2.4 (relative error < 1e-6, checked in C0.4), interpolated to each SN and BAO redshift.

**Best fits.** As conventions.json `best_fits`: LCDM by bounded scalar minimisation (Brent); w0waCDM by Nelder–Mead then Powell from the lowest grid point and four fixed starts ((0.3, −1, 0), (0.35, −0.8, −1), (0.45, −0.6, −5), (0.25, −1.1, 1)); the lowest minimum is the best fit, the spread of the five minima is reported.

**Posteriors.** Exact grid evaluation (conventions.json `posteriors`), no sampler: 161 × 151 × 241 points, likelihood with the offset marginalised analytically. Each sample's likelihood grid is computed once and saved (`checks/out/grids/*.npz`, kept outside Git, SHA-256 in `checks/GRIDS.sha256`); priors are masks on it. Grid convergence: C0.4.

**Statistics.** 1D: marginal mean, median, standard deviation, 68 % equal-tail interval. 2D: credible level of a point (conventions.json). σ conversions as conventions.json `sigma_statements`.

## 2. Controls (run first)

- **C0.1 Data identity.** DES-SN5YR MU, MUERR_FINAL and zHD in the Cobaya file vs the DES release (des-science/DES-SN5YR tag 1.3, `DES-SN5YR_HD+MetaData.csv`, commit e3493cb3): equal for all 1,829 SN.
- **C0.2 Radiation and redshift convention.** At the Planck point and at each best fit: χ² with and without radiation (Ω_r h² = 4.18e-5, h = 0.6744); D_L with (1+z_HEL) vs (1+z_HD). Reported; |Δχ²| > 0.3 from radiation goes in an addendum before the claim runs.
- **C0.3 Independent likelihood check.** Our χ² vs Cobaya 3.6.2's `sn.pantheonplus`, `sn.desy5` and `bao.desi_2024_bao_all` at three fixed points (Planck point; (0.35, −0.8, −1); (0.45, −0.6, −5)), fed by a minimal Cobaya Theory with our background, full samples as Cobaya cuts them. Agreement within 0.01 in χ² after accounting for the known differences (Cobaya marginalises the offset; Cobaya's Pantheon+ cut z > 0.01): pass/fail.
- **C0.4 Numerics.** Distance integral vs scipy quad at 5 points (< 1e-6 relative); grid best χ² minus minimiser best χ² (< 0.05 expected for LCDM, reported for w0wa); posterior means and the C5 levels on a grid with half the step in a sub-box around the peak (changes reported; > 0.01 in a level goes in the report as numerical uncertainty).

## 3. Checks, tolerances and outcomes

Tolerances are for `held:reproduced`. A miss is first checked against our transcription and the conventions in conventions.json (z convention, δm cut, offset treatment) before it is called a discrepancy (G5).

| Claim | Check | Tolerance |
|---|---|---|
| C1 | Row counts after the cuts; DES5Y composition by IDSURVEY and redshift range; Pantheon+ DES Y3 rows | Exact. A count miss is investigated (alternative reading of δm, z column) and reported with the reading that reproduces it, if any |
| C2 | χ² at the Planck point; LCDM best fit and χ²; Δχ² (Planck − LCDM fit); Ω_m with error three ways | **Fidelity gate**, both samples: Planck χ² within ± 1.0 of 1259.6 and 1644.68; LCDM best-fit Ω_m within ± 0.003 of the Fig. 1 legends (0.333, 0.349); LCDM χ² within ± 1.0 of 1258.4 and 1639.6; Δχ² within ± 0.3 of 1.2 and 5.1. Reported, not gating: text Ω_m (0.333 ± 0.018, 0.354 ± 0.017) against our posterior mean ± sd and profile half-width (± 0.005 on the value, ± 0.003 on the error) |
| C3 | w0waCDM best fits | χ² within ± 1.0 of 1258.4 and 1632.2; Δχ² within ± 0.5 of 12.5 (vs Planck) and 7.4 (vs DES5Y LCDM); Pantheon+ 'no improvement': Δχ²(LCDM fit − w0wa fit) < 1.0. Best-fit parameters reported against the legends (w0, wa are poorly constrained by SN alone: differences along a flat direction with χ² within 0.3 are not a miss) |
| C4 | σ statements, from the printed numbers and from ours | 2.2σ: Gaussian difference (conventions) within ± 0.25. 'Roughly 3σ' and 'about 3σ': k = 3 conversion of Δχ² within [2.5, 3.5]. '~2σ': k = 2 conversion within [1.5, 2.5]. If only another k fits, the statement holds with an S0 note on the unstated convention |
| C5 | Grid posteriors, printed priors (and the pre-declared Pantheon+ wa variant) | Pantheon+: credible level of (Ω_m, w0, wa) = (0.3135, −1, 0) ≤ 0.6827 in all three 2D projections. DES5Y: marginal median Ω_m ≥ 0.40, P(w0 > −1) ≥ 0.84, median wa ≤ −2. Levels of the Planck point for DES5Y reported (related to C4) |
| C6 | Cross-match (conventions `common_sn`) | Counts exact (203, 145, 150, 118 and the Table 1 N); each mean within ± 0.003 of the printed value; errors within ± 20 % of the printed value; dispersion 0.067 within ± 0.005; z_median within ± 0.002. Headline: (all-low-z mean) − (DES mean) within [−0.05, −0.03] ('~0.04 mag') |
| C7 | Shifted DES5Y (conventions `lowz_shift`): best fits and grid posterior | LCDM Ω_m within ± 0.003 of 0.310 and χ² within ± 1.0 of 1635.7; w0wa χ² within ± 1.0 of 1633.4; Δχ²(LCDM − w0wa) within ± 0.5 of the legends' 2.3 (text: 2.5); 'not significant': k = 2 conversion < 2σ; 'near perfect agreement': |Ω_m − 0.3135| < 0.25 × the SN-only Ω_m error. Posterior: level of (w0, wa) = (−1, 0) ≤ 0.9545 in the w0–wa projection ('consistent'; ≤ 0.6827 reported); 'skewed to high Ω_m': marginal median Ω_m > 0.3135 |
| C8 | Joint DESI DR1 BAO + SN grid posteriors (conventions `bao`), DES5Y prior box | BAO alone: level of (w0, wa) = (−1, 0) between 0.6827 and 0.9876 ('mild', 1σ–2.5σ). DESI + Pantheon+: level of the Planck point ≤ 0.80 in all three projections ('touches the 1σ contours'). DESI + DES5Y: level of (−1, 0) in the w0–wa projection converts to [2.5, 3.5]σ ('~3σ'); joint marginal median Ω_m < 0.36 while the DES5Y-alone median is ≥ 0.40 ('high Ω_m disfavoured by DESI BAO'). Judged on the primary α prior; the variant is reported |
| C9 | Quoted numbers against the cited TeX; V25 sums | Exact after half-up rounding to the printed digits; sums exact. A mismatch is an S0 note unless the quoted number carries a conclusion of this paper (then S1) |

**Severity, fixed now (G18, 05 §3).** C1/C2 misses: transcription and convention checks first, then Track R questions (exact selection, χ² definition); not a finding unless a conclusion changes. C3/C4: a Δχ² off by more than 1.0 or a σ statement off by more than 0.5σ in the abstract's direction is S1; one that moves the abstract's "about 3σ" below 2.5σ under every k is an S2 candidate (full lane). C5: a Pantheon+ level > 0.6827 in one projection is S1 (wording), > 0.9545 in any is an S2 candidate on "consistent with Planck LCDM". Printed priors that contradict the plotted contours: S0. C6: a Table 1 error that disagrees with the scatter it summarises: S0 if the text or figure gives the right value; if our all-low-z minus DES differential is outside [−0.05, −0.03] or its error makes it consistent with zero at 2σ, S2 candidate on the abstract's 0.04 mag (full lane). C7: if the shifted DES5Y still prefers w0waCDM over LCDM by more than 2σ (k = 2) or its LCDM Ω_m differs from 0.3135 by more than 2 × the SN error, S2 candidate on the abstract's "can bring the DES5Y sample into good agreement" (full lane, impact run = this check). C8: a level > 0.9545 for DESI + Pantheon+ in any projection is S1 (appendix wording); a DESI + DES5Y level below 2σ is S1. C9: S0–S1.

## 4. Data re-tests (only after the C2 fidelity gate passes)

Data traces, not verification; they never change a claim outcome (G17).
- **M1 DES-Dovekie** (update role; des-science file with CIDs and its STAT+SYS inverse covariance, full sample as distributed, no cuts beyond 0.02 ≤ z_HD ≤ 1.2): (a) the C6 cross-match against Pantheon+ (same rule, MU − 19.33) and the differential (all low z − DES); (b) SN-only LCDM and w0waCDM best fits, χ² at the Planck point, Δχ² and σ (C2–C4 recipe); (c) grid posterior, level of the Planck point (C5 recipe, DES5Y prior box).
- **M2 DESI DR2 BAO** (update role): the C8 recipe with DR2 alone and joined with Pantheon+, DES5Y and DES-Dovekie: level of (w0, wa) = (−1, 0) and of the Planck point.
Metric: the numbers above with the same conventions; no claim outcome changes.

## 5. Literature (G21), after the runs

INSPIRE (record 2818064, citing list saved on the Mac on 2026-10-15, 180 records): papers that re-examine the DES5Y/Pantheon+ low-z offset or redo the comparison (V25, the Dovekie papers 2506.05471 and 2511.07517, and others found by title); errata; the v1 → v2 → v3 changes (PDFs kept). Results to `literature.md` with queries and dates.

## 6. Stop rules

Claim cap 10, 9 used. Every number in the report is printed by a script in `checks/` into `checks/out/` (G14). Effort is measured from the session log at the end of each stage (G28). Compute cap for the original run: 3 h of wall clock from the first likelihood grid for all grids, best fits and re-tests; what has not finished by then is `not_run`.

## Addendum A1 (2026-10-15 09:35 UTC, after the first fits, before any grid or posterior): redshift convention of D_L

**What has run.** Controls C0.1–C0.4 (`checks/out/c0_controls.out`): all pass; C0.2 shows that the D_L convention (1+z_HEL) vs (1+z_HD) moves χ² by up to 0.97, radiation by at most 0.05. Then `checks/fits.py` computed, in one script, C1 counts, the C2 gate, the C3 and C7 best fits, the C4 conversions and the M1 Dovekie fits, all under the primary (1+z_HEL) convention, with χ² under (1+z_HD) evaluated at the same best fits. The checks had therefore produced the C3, C4, C7 and M1 numbers under (1+z_HEL) before this addendum. (A smoke test of `checks/sn.py` printed the two Planck-point χ² before the controls ran; no other number.)

**Gate result.** Under (1+z_HEL) every gate number passes except one: DES5Y Δχ²(Planck − LCDM) = 5.49 against 5.1 ± 0.3. Under (1+z_HD), at our (1+z_HEL) best fits, the four gate χ² are 1259.52, 1258.45, 1644.54, 1639.63 against the printed 1259.6, 1258.4, 1644.68, 1639.6, i.e. all within 0.15.

**Decision (the §3 diagnosis, planned).** The paper's own convention is taken to be D_L = (1+z_HD) D_C(z_HD), subject to a direct test run next: χ² at the printed best-fit parameters of Figs 1 and 4 under both conventions (whichever reproduces the printed χ² is the paper's). Then:
- Every claim is reported under both conventions. Claim outcomes are judged under the paper's convention (G26: the paper's own inputs first); the (1+z_HEL) results, the convention of the Pantheon+ and DES releases and of Cobaya, are the "best current inputs" and are reported beside them. A claim that holds under one and not the other gets a note.
- The C2 gate is re-evaluated under the paper's convention with re-minimised fits.
- Grids for C5, C7 and C8 are computed under both conventions (the cap of §6 is unchanged).
- New check, not in §3 (planned now, before it runs): χ² at the printed best-fit parameters of Figs 1 and 4. Our w0waCDM minima lie below the printed χ² under (1+z_HEL); if the printed best fits are not minima under either convention, that is an S0 note (the best-fit method is unstated) unless it moves a Δχ² beyond the C3/C7 tolerances, which is S1.

## Addendum A2 (2026-10-15 15:00 UTC, after the scientific runs): figure numbers

The original figure numbering was one too low from the offset figure onward. In TeX caption order (`SNpaper.tex`):

| Fig. | File | Content |
|---|---|---|
| 1 | `Panfit`, `DESfit` | residuals and best fits |
| 2 | `PanDES*` | SN-only contours |
| 3 | `PanvDES1/2` | common-SN magnitude differences |
| 4 | `pgPanDES_offset` | mean offsets by survey |
| 5 | `ODESfit` | shifted DES5Y residuals; the 2.3 legend |
| 6 | `PanCDES*` | shifted DES5Y contours |
| 7 | `bias` | bias corrections |
| 8 | `DESI+*` | DESI DR1 BAO with SN (placed before `\appendix`) |

The reference to "Figs 1 and 4" in A1 means Figs 1 and 5. The frozen locations of C7 ("Fig. 4 legends, Fig. 5") and C8 ("Fig. 7") are corrected in new revision files `claims/C7-r2.json` and `claims/C8-r2.json`. Only the location changes; the statement, method and criteria are those of r1. No number and no outcome changes: every check read its values from `paper_values.json` by key, and its key `fig4_DES5Y_offset` holds the Fig. 5 legend values.
