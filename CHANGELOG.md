# Changelog

All notable changes are listed here. Versions follow [semantic versioning](https://semver.org/).

## Unreleased

- Added `CONTRIBUTING.md` (how to report problems, project rules, support expectations) and this changelog.
- Added the author ORCID to `CITATION.cff` and `paper/paper.md`.
- Paper: added State of the field, Software design, Research impact statement and AI usage disclosure sections.
- Tests: fitted iksi and km in the guard-off comparison now use relative 1e-6 (an intermittent CI failure on
  Python 3.10 showed a relative 1.6e-8 run-to-run difference). The called flags are still compared exactly.
- CI: added a workflow that builds the paper PDF with the Open Journals draft action.

## 0.2.0 (2026-11-20)

DOI: [10.5281/zenodo.24409789](https://doi.org/10.5281/zenodo.24409789). On PyPI as `ggc-kinetics==0.2.0`.

- `k57_map_fit`: new opt-in `dispersion_guard=True`. It refits each plate at `sig * birge` when the null fit's
  Birge ratio is above 1. The result carries `"birge"` (all 1.0 when the guard is off). The default output is
  unchanged from 0.1.0.
- README: new "K57 noise scale" section (Q65 study, MODELED) on passing an estimated `sig`.
- Tests: 18 unit tests (4 new for the guard). Floating-point comparisons for the guard use rtol 1e-10 so the suite
  passes on Linux as well as Windows; test fixtures are stored byte-exact.
- Streamlit demo app checked by hand in a browser.
- CI on Python 3.10 and 3.12; PyPI upload by Trusted Publishing from GitHub Actions.

## 0.1.0

- First packaged release of the solvers, detectors, parser and `ggc-fit` CLI, copied without changes to the
  mathematics from the GGC-inf research archive. Not published to PyPI.
