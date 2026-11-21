# Contributing to ggc-kinetics

Thanks for your interest. Bug reports, real-data test cases and pull requests are all welcome.

## Reporting a problem

Open an issue at https://github.com/bpmcginley/ggc-kinetics/issues. Please include:

- the `ggc-kinetics` version (`python -c "import ggc_kinetics; print(ggc_kinetics.__version__)"`), Python version and OS;
- the smallest input that shows the problem (a plate export, or the arrays passed to a solver);
- what you expected and what you got. For a parse problem, attach the `FILE.report.json` that `ggc-fit` wrote.

Real plates are the most useful thing you can send. If you have an openly licensed plate with raw progress curves at
several substrate levels and matched blank wells, please open an issue: the detector thresholds have not yet been
tested on such a plate (see "Not claimed" in the README).

## Making a change

```
git clone https://github.com/bpmcginley/ggc-kinetics
cd ggc-kinetics
pip install -e ".[app]"
python -m unittest discover -s tests -v
python examples/quickstart.py
```

1. Open an issue first for anything larger than a small fix, so the approach can be agreed before you write code.
2. Work on a branch and open a pull request against `main`. CI runs the tests on Python 3.10 and 3.12; it must be
   green before merging.
3. Add or update a test for every behaviour change.

## Project rules

- **Evidence tags.** Every number in the docs, tests or reports carries a tag: MEASURED (from real data), MODELED
  (from simulation, or a model fitted to data) or ASSUMED (a chosen value). New numbers need one.
- **Frozen thresholds.** Detector thresholds and decision rules are not re-tuned to make a dataset pass. A change to
  a threshold needs a pre-registered study with the gates written down before the data is seen.
- **Default output is stable.** Changes must not alter the default output of an existing function. New behaviour goes
  behind an opt-in argument (as `dispersion_guard` did in 0.2.0) and is listed in `CHANGELOG.md`.
- **Parser errors are codes.** The parser returns an error code, never a silent wrong number or an uncaught
  exception.

## Support

ggc-kinetics is maintained by one person (Bruce McGinley) on a best-effort basis. Questions, bugs and feature
requests go to the issue tracker; the aim is a first reply within two weeks. There is no paid support. Releases
follow semantic versioning, and every release is archived on Zenodo.

## Code of conduct

Be respectful and constructive. Harassment or personal attacks are not accepted in issues, pull requests or
discussions, and the maintainer may remove such content.
