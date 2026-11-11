# Evolving Dark Energy or Supernovae Systematics?

[2408.07175v3](https://arxiv.org/abs/2408.07175v3) · Review 0.2.0-rc3

**Reproduced in scope.** The supernova fits and common-sample magnitude offset reproduce under the inferred distance convention and stated prior variants.

**Limit:** The offset does not establish a full-sample calibration error; correlated cross-catalogue systematic uncertainty remains unmeasured.

## Checked scope

- Sample selection and counts for Pantheon+ and DES5Y (C1)
- χ² of the Planck ΛCDM cosmology and the ΛCDM fits to both catalogues (C2)
- The w0waCDM fits and their χ² improvements (C3)
- The significance statements in the text and abstract (C4)
- The supernova-only posteriors and where Planck ΛCDM lies in them (C5)
- Table 1: the magnitude offsets of supernovae common to both catalogues, and the low-redshift minus DES difference (C6)
- DES5Y with 0.04 mag subtracted from its low-redshift supernovae: fits and posterior (C7)
- Appendix A: DESI DR1 BAO alone and combined with each supernova catalogue (C8)
- Numbers quoted from Planck 2018, DESI 2024 and the DES analysis of the offset (Vincenzi et al. 2025) (C9)

Not checked: Whether the offset seen in the common supernovae implies a calibration error in the full DES5Y sample. That needs the light curves, bias-correction simulations and cross-catalogue systematic covariance, which the public distance tables do not carry. The paper's 'can bring … into good agreement' is conditional, and the shift experiment supports it; The DES analysis's own decomposition of the offset (selection, host-mass and intrinsic-scatter terms), beyond checking that the quoted numbers match their source; The statement that binned and unbinned likelihoods give nearly identical results (l.257), and the attribution of the remaining asymmetry of the shifted posterior to supernovae at z \> 1 (l.425–426); The paper's own chains and plotting settings, which are not public. Our posteriors are grids over the prior box; for the Pantheon+ posterior we used a w\_a range that contains the plotted contours (note N2); The ± 0.008 mag on the offset is the scatter of the common supernovae only, without correlated systematic errors between the catalogues, so it is not a calibrated significance

## Review record

| Claim | Method | Outcome / adequacy | Evidence |
|---|---|---|---|
| C1: C1 in fixed reviewed version — Counts and DES5Y composition exact (1417 and 1754 after the cuts; A and B independently) | source\_check | held:reproduced / pass | [checks/out/fits.out](checks/out/fits.out) |
| C2: C2 in fixed reviewed version — Fidelity gate passed under D\_L = (1+zHD) D\_C(zHD), the paper's unstated convention; (1+zHEL) misses the DES5Y Delta chi2 (5.49 vs 5.1) (N1) | source\_check | held:reproduced / pass | [checks/out/fits.out](checks/out/fits.out) |
| C3: C3 in fixed reviewed version — w0waCDM chi2 and improvements within tolerance; printed best-fit points 0.14/0.15 above the minima along the flat direction (N3) | source\_check | held:reproduced / pass | [checks/out/fits.out](checks/out/fits.out) |
| C4: C4 in fixed reviewed version — ~2.2 sigma 2.15 (printed numbers) / 1.98; roughly 3 sigma 2.76 (k=3); ~2 sigma 2.29 (k=2) | source\_check | held:reproduced / pass | [checks/out/fits.out](checks/out/fits.out) |
| C5: C5 in fixed reviewed version — Planck point inside 1 sigma for Pantheon+ with a wa range containing the plotted contours (0.429/0.383/0.070); the printed range \[-4,-0.4\] cannot have been used (N2); DES5Y thresholds pass | source\_check | held:discrepancy / pass | [checks/out/posteriors.out](checks/out/posteriors.out) |
| C6: C6 in fixed reviewed version — All counts and means reproduce, differential -0.0360 +- 0.0080 (scatter only); Table 1 DES/FOUND errors x10 typo, three z\_median, '150 Foundation' (N4-N6) | source\_check | held:discrepancy / pass | [checks/out/crossmatch.out](checks/out/crossmatch.out) |
| C7: C7 in fixed reviewed version — Shifted DES5Y: LCDM Omega\_m 0.310, no significant preference (1.19 sigma); printed w0wa point not a minimum, Delta chi2 2.89 vs legend 2.3 (N7) | source\_check | held:discrepancy / pass | [checks/out/fits.out](checks/out/fits.out) |
| C8: C8 in fixed reviewed version — Appendix statements hold under the flat-H0 prior; unstated H0 prior moves the DESI+Pantheon+ w0-wa level 0.74 -\> 0.84 (N10; protocol rule S0) | source\_check | held:discrepancy / pass | [checks/out/posteriors.out](checks/out/posteriors.out) |
| C9: C9 in fixed reviewed version — 12 of 13 quoted numbers match; 0.0344 for 0.344 (N8); V25's weighted Foundation-DES -0.042 compared with the unweighted all-low-z offset (N9) | source\_check | held:discrepancy / pass | [checks/out/c9_citations.out](checks/out/c9_citations.out) |

[Protocol and provenance](PROTOCOL.md). This review record is retrospective.

## Review provenance

Prepared with AI assistance. Historical Claude first pass and Codex independent review/critique remain distinct; this assembly reads saved evidence and is neither blind B nor new P2. No scientific code was run in this assembly. Before the replay command below, extract legacy-original.zip into a new directory and read its README.md and requirements.txt for installation, dependency, data and scope instructions.

Tools: Claude (historical A); Codex (historical B/C and current informed assembly). Human involvement: Owner selected the preparation scope and approved the predecessor release published on 10 November 2026; the owner approved publication of this editorial revision on 11 November 2026.

## Replay

~~~sh
mkdir replay-work && python3 -m zipfile -e legacy-original.zip replay-work && cd replay-work && python3 reproduce.py --all
~~~

See original source-payload README and recorded P2; not rerun during assembly.. Unpacks the unchanged original candidate into a new work directory and runs its documented scoped replay. Package validation checked the artifact and inherited the earlier replay after verifying the relevant files were unchanged; it did not rerun the scientific calculation.
