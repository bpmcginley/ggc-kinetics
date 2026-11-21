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

# State of the field

General kinetic-modelling programs such as DynaFit [@kuzmic1996dynafit] and KinTek Explorer fit arbitrary
mechanisms by numerical integration; they are powerful but are separate applications, not Python libraries, and
leave the choice and checking of the integrator to the user. ICEKAT [@olp2020icekat] is a web tool for initial
rates, not whole curves. Closed forms exist in the literature for single cases [@schnell1997closed;
@golicnik2020pon1] but are scattered across papers. General Python fitting tools (SciPy's `curve_fit` or ODE
solvers [@virtanen2020scipy]) supply the machinery but no enzyme-specific rate laws, plate parsing or detection
rules.

*Build versus contribute.* The solvers depend only on NumPy and SciPy and use `scipy.special.lambertw` rather than
re-implementing it. A new package was still needed because none of the tools above combines (i) exact forward
models for these extensions, (ii) a stated error bound for each against a reference integrator, and (iii) frozen,
evidence-tagged decision thresholds; adding (iii) to a general fitting library would impose a policy that does not
belong there. The parser is self-contained because plate-reader exports vary by vendor and need explicit error
codes.

# Software design

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

Three design choices shape the package. First, **forward models without a general ODE integrator**: closed forms
or split-step solutions are fast inside a least-squares loop and have no step-size setting to get wrong; the cost is
that each rate law is a separate function, and each one is checked against fourth-order Runge-Kutta. Second,
**frozen decision rules**: thresholds are fixed constants copied from the archive, and new behaviour is opt-in
(for example `dispersion_guard`), so results from an older version stay reproducible. Third, **errors as data**:
the parser returns an error code inside the report instead of raising, so a batch run over many plates never stops
on one bad file and never returns a silent wrong number.

# Verification and real-data benchmarks

**Numerical verification (MEASURED numerics on MODELED inputs).** The Morrison solver agrees with fourth-order
Runge-Kutta to 2.5e-9; the coupled-lag form to below 1e-11; slow binding to below 1e-8. The substrate-inhibition
Jacobian matches finite differences to 5.3e-10, and its solver residual is at most 4.5e-14 over 600 simulated
plates. The 18 unit tests run on Python 3.10 and 3.12 in CI.

**Real data (MODELED fits of MEASURED data, no re-tuning).** On the 1913 invertase data
[@michaelis2011translation] the frozen solver gives Km 16.52 mM (+0.11% against the reference) and prefers product
inhibition over enzyme decay by dBIC -19.29. Published paraoxonase-1 values [@golicnik2020pon1] are pinned in the
tests. A pre-registered benchmark on ICEKAT files [@olp2020icekat] and peroxidase data [@kont2020lpmo] passed 3 of
5 gates; the two failures (one missing cell, and a residual of 3.6e-7 at deep substrate depletion) are reported, not
re-gated. On one real, already-inspected laccase plate [@haussler2024darus] the curvature score separated enzyme
rows (3.48 to 3.82) from no-enzyme rows (0.32 to 0.50) but did not reach the 4.60 threshold (2 of 3 gates). A
second deposit [@ngubane2022darus] was too scattered to call one Km.

**Limits.** All detector thresholds and the substrate-inhibition rule are ASSUMED values calibrated on MODELED
plates; simulated studies (Q64, Q65) showed the rule needs a correct readout SD, which `dispersion_guard` or an
estimated SD provides. None of these has passed a test on a fresh real plate: two pre-registered searches of 12
public datasets found no plate with raw multi-level curves and matched blanks. The parser does not read `.xlsx`.

# Research impact statement

`ggc-kinetics` was extracted from the author's own GGC-inf research archive (66 analysis cycles), where these
solvers and detectors were used for every progress-curve fit and simulated-plate study, and it is the analysis code
behind the real-data benchmarks described above. It is released with tests, CI, a command-line tool and an archived
DOI so that others can reuse it. There are, as yet, no external users or publications that cite it.
**TODO (author): before submitting, add concrete evidence of research use (a preprint, poster, lab or course using
the package), as JOSS requires.**

# AI usage disclosure

**TODO (author): complete and confirm this statement; JOSS treats an incomplete disclosure as an ethical breach.**
Claude Opus 5.5 (Anthropic), through Claude Code, was used during the v0.2.0 release for packaging and CI set-up,
for adjusting floating-point tolerances in the tests, for copy-editing the README and changelog, and for drafting
the State of the field, Software design, Research impact and AI usage sections of this paper.
**TODO (author): list any other tools and versions used in the GGC-inf research archive (code generation,
refactoring, test scaffolding, review) and where.** The author reviewed, edited and validated all AI-assisted
output and made the core design decisions. **TODO (author): confirm this sentence is accurate.**

# Acknowledgements

We thank the authors of the openly licensed datasets used in the benchmarks. The package builds on NumPy
[@harris2020numpy] and SciPy [@virtanen2020scipy]. **TODO (author): state funding sources, or that the work
received no specific funding, and any conflicts of interest.**

# References
