---
title: 'ggc-kinetics: exact enzyme progress-curve solvers, plate-level detection statistics and a plate-reader parser'
tags:
  - Python
  - enzyme kinetics
  - progress curves
  - Michaelis-Menten
  - Lambert W
  - plate reader
authors:
  - name: Bruce McGinley
    orcid: 0009-0001-9263-6174
    affiliation: 1
affiliations:
  - name: University of Massachusetts Amherst
    index: 1
date: 20 November 2026
bibliography: paper.bib
---

# Summary

`ggc-kinetics` is a small Python package (NumPy and SciPy only) for analysing whole enzyme progress curves rather
than initial rates. It provides solvers that need no general-purpose ODE integration for six common rate laws: the
integrated Michaelis-Menten equation through the Lambert W function, with product inhibition or first-order enzyme
decay; a Strang-split extension with non-enzymatic substrate loss; a coupled-assay lag with a decaying primary rate;
slow-onset inhibition by time warping; tight binding (Morrison) by a Newton solve; and substrate inhibition with a
per-well pipetting prior. Two detection statistics come with it, a plate-level saturation-curvature score and a
melt-transition score for thermal-shift curves. A parser reads text exports from common plate readers and turns them
into a JSON and Markdown report. A command-line tool, `ggc-fit`, wraps the parser, and an optional Streamlit app
shows the parsed traces. Version 0.2.0 adds an opt-in noise-scale guard for the substrate-inhibition test.

Every number the package or its documentation reports carries an evidence tag: MEASURED (published data),
MODELED (simulation, or a model fitted to data, written "MODELED fit of MEASURED data") or ASSUMED (a chosen value
such as a threshold). The tags are there so that thresholds calibrated on simulated plates are never mistaken for
validated real-data performance.

# Statement of need

Initial-rate analysis throws away most of a progress curve and is biased when substrate depletion, product
inhibition or enzyme decay bend the curve early [@johnson2013century]. Fitting the full curve avoids this, but it
usually means numerical integration inside a nonlinear least-squares loop, with step-size and stiffness choices
that users rarely check. Closed forms exist for the simplest case [@schnell1997closed], and interactive tools such
as ICEKAT [@olp2020icekat] make initial-rate fitting easy, but there is no small, dependency-light library that
(i) gives exact or split-step solutions for the common extensions, (ii) states the numerical error of each one
against a reference integrator, and (iii) says plainly which of its decision thresholds have been tested on real
data. `ggc-kinetics` is aimed at enzymologists and assay developers who script their analysis in Python and want a
fast, checkable forward model and a parser that fails with an error code, not a silent wrong number.

# Key capabilities and architecture

The package has four modules.

* `solvers`: `kp_kd_exact` and `fit_kp_vs_kd` use the Lambert W function [@corless1996lambertw] for the integrated
  rate law with product inhibition (Kp) or enzyme decay (kd), and choose between them by BIC.
  `unified_kp_strang` adds non-enzymatic hydrolysis by Strang splitting [@strang1968splitting].
  `coupled_lag_decay` and `slow_binding_warped` are closed forms. `morrison_exact` solves the tight-binding
  quadratic [@morrison1969tight] along the curve by Newton iteration. `k57_map_fit` fits a substrate-inhibition term
  with a prior on per-well pipetting error, with an analytic Jacobian. Its optional `dispersion_guard` rescales the
  readout SD by a one-sided Birge ratio from the null fit; the default leaves 0.1.0 results bit-identical.
* `detectors`: `z_rao_orth_robust` is a robust score test for saturation curvature across the S0 levels of one
  plate, made orthogonal to first-order decay and a lag term. `z_dsf_asym_boltz` scores a sigmoid melt transition
  over a quadratic baseline.
* `parser`: a frozen copy of the archive's plate parser, pinned by SHA-256 in the tests. It reads CSV, TSV and text
  exports in long or wide form, decimal commas and SoftMax Pro blocks, skips instrument preambles, and returns
  errors as codes (`EMPTY`, `NO_WELL_TABLE`, `TIME_NOT_MONOTONIC`).
* `cli`: `ggc-fit FILE --evidence MEASURED|MODELED --out DIR` writes `FILE.report.json` and `FILE.report.md`, and
  exits 0 on success and 2 on a parse error.

All functions are copied without changes to their mathematics from a 66-cycle research archive. A parity script
checks the release against the archive: 20 of 20 checks pass, most with a difference of exactly 0.

# Verification and real-data benchmarks

**Numerical verification (MEASURED numerics on MODELED inputs).** The Morrison solver agrees with fourth-order
Runge-Kutta to 2.5e-9 with residual 6.2e-15; the coupled-lag form to below 1e-11; slow binding to below 1e-8. The
substrate-inhibition Jacobian matches finite differences to 5.3e-10 (gate 1e-6), and the solver residual is at most
4.5e-14 over 600 simulated plates (gate 1e-8). The 18 unit tests run in about 2 s, and two runs give identical logs
apart from timing.

**Real data (MODELED fits of MEASURED data, no re-tuning).** On the 1913 invertase data
[@michaelis2011translation], the frozen solver gives Km 16.52 mM (+0.11% against the reference) and prefers product
inhibition over enzyme decay by dBIC -19.29. The published paraoxonase-1 comparison [@golicnik2020pon1], in which
a Lambert-W fit stays within 0.69% of an ODE fit and DynaFit deviates by up to 17.2%, is pinned in the tests as
MEASURED reference values (published numbers, not refitted here). In a pre-registered
release benchmark on ICEKAT supplementary files [@olp2020icekat] and horseradish-peroxidase source data
[@kont2020lpmo], 3 of 5 gates passed. The parse, the peroxidase Km comparison (10 of 12 series within 10% of the
authors' values) and the runtime gates passed. Two gates failed: an exact-identity check between two copies of one
file (one cell missing, others equal to 2.2e-16), and a residual gate of 1e-8, which reached 3.6e-7 at a point
where only 2e-9 of the substrate remains. The second failure is a conditioning effect at deep depletion; the gate was
not changed after the fact.

On one real laccase plate [@haussler2024darus], which had already been inspected and is therefore not a fresh
test, the plate-level curvature score separated the three enzyme rows (3.48 to 3.82) from the three no-enzyme rows
(0.32 to 0.50). The controls stayed below the 4.60 threshold, but the enzyme rows did not reach it, so the
pre-registered active-row gate failed (2 of 3 gates). On the same plate the noise-adaptive K59 variant (tau 0.056,
threshold 1.768, both ASSUMED) did separate all three enzyme rows (2.26 to 2.61) from the three controls (0.045 to
0.071); the verdict stays 2 of 3, and K59 still needs a fresh real plate. A second laccase deposit [@ngubane2022darus] gives replicate
Km values too scattered (0.60 dex) to call one Km.

**Noise scale of the substrate-inhibition rule (MODELED).** The rule treats the readout SD as known. Two
pre-registered studies on fresh simulated plates from an independent mass-action model, with every threshold
frozen, tested what happens when it is not. In Q64 a noise floor twice the nominal one raised the false-call rate
to a Wilson 90% upper bound of 0.099, against a gate of 0.060. Q65 (7 arms of 450 plates) confirmed the failure:
the nominal SD gave 0.034 at the nominal floor, 0.076 at 2x and 0.131 at 3x. Two fixes passed every arm without
re-tuning. V1 passes an SD estimated from blank reads (w90 at most 0.037, recall 0.991 and 0.982), but its 2%
proportional noise term is ASSUMED and matches the generator, which favours it. V2, the Birge guard now shipped
as `dispersion_guard=True`, needs no blanks; its worst arm, the 2x floor, sat close to the gate (w90 0.0555
against 0.060). Neither fix has been tested on a real plate.

**Limits.** The detector thresholds (4.60 and 1.768 for the curvature score, 5.72 for the melt score) and the
substrate-inhibition decision rule are ASSUMED values calibrated on MODELED plates. The substrate-inhibition rule is
valid in simulation only for peak absorbance up to 2.0 and pipetting CV up to 5%. The Morrison solver needs the
enzyme concentration to within 10%. None of these has yet passed a test on a fresh real plate: two pre-registered searches screened 12 public
datasets for a fresh plate with raw multi-level progress curves and matched blanks, and none qualified, so K59
remains unvalidated on fresh real data. The parser does not
read `.xlsx` files. These limits are listed in the README and close only when a fresh, openly licensed real plate
with several substrate levels and blank rows passes the pre-registered gates.

# Acknowledgements

We thank the authors of the openly licensed datasets used in the benchmarks. The package builds on NumPy
[@harris2020numpy] and SciPy [@virtanen2020scipy].

# References
