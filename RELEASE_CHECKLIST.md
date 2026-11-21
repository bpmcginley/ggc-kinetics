# Release checklist: ggc-kinetics v0.2.0

Everything up to the upload commands has been run and checked (Release Phase 6, v0.2.0-rc1). The steps below are manual.

## 0. The only two fields you must confirm

| Placeholder | Where | Replace with |
|---|---|---|
| `GITHUB_HANDLE` | `pyproject.toml` (`[project.urls]`), `CITATION.cff` (`repository-code`) | your GitHub user or org name |
| ORCID `0000-0000-0000-0000` | `CITATION.cff`, `paper/paper.md` | your real ORCID iD |

Find every occurrence (PowerShell, from this folder):

```powershell
Select-String -Path pyproject.toml, CITATION.cff, paper\paper.md -Pattern "GITHUB_HANDLE", "0000-0000-0000-0000"
```

Also check that the author name and affiliation are right in all three files.

## 1. GitHub

The repository ships `.github/workflows/ci.yml` (listed in `MANIFEST.in`). On every push and pull request it runs
on ubuntu-latest with Python 3.10 and 3.12: `pip install -e .`, `python -m unittest discover -s tests`, a `ggc-fit`
CLI smoke run on `tests/data/plate_healthy_biotek_wide.csv`, and `examples/quickstart.py`. Check that the first
run after the push is green before going on.

```powershell
git init
git add .
git commit -m "ggc-kinetics v0.2.0"
git branch -M main
git remote add origin https://github.com/bpmcginley/ggc-kinetics.git
git push -u origin main
git tag v0.2.0
git push origin v0.2.0
```

## 2. Build and check

```powershell
Remove-Item -Recurse -Force dist -ErrorAction SilentlyContinue
python -m pip install --upgrade build twine
python -m build
python -m twine check --strict dist/*
```

`dist/` must hold exactly two files: `ggc_kinetics-0.2.0-py3-none-any.whl` and `ggc_kinetics-0.2.0.tar.gz`.

## 3. TestPyPI, then a fresh-venv install

v0.2.0 was uploaded with PyPI Trusted Publishing instead of `twine upload`: `.github/workflows/publish.yml`
(manual `workflow_dispatch`, input `target` = `testpypi` or `pypi`, GitHub environments of the same names). Register
a pending trusted publisher on test.pypi.org and pypi.org for owner `bpmcginley`, repo `ggc-kinetics`, workflow
`publish.yml`, and the matching environment. The `twine upload` lines below are the token-based fallback.

```powershell
python -m twine upload --repository testpypi dist/*
```

```powershell
python -m venv $env:TEMP\ggc_tp; & $env:TEMP\ggc_tp\Scripts\python.exe -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ ggc-kinetics
```

```powershell
& $env:TEMP\ggc_tp\Scripts\ggc-fit.exe tests\data\plate_healthy_biotek_wide.csv --evidence MODELED --scale 0.05 --offset 0.040 --out $env:TEMP\ggc_tp_out
```

That should exit 0 and print P1 KILL (expected; see README). Then upload to PyPI:

```powershell
python -m twine upload dist/*
```

## 4. Browser check of the Streamlit app

```powershell
pip install ".[app]"
streamlit run app/streamlit_app.py
```

Upload `tests/data/plate_healthy_biotek_wide.csv`, set evidence MODELED, scale 0.05 and offset 0.040, and check that
the report and traces appear. Only after this works, remove the sentence "It has not been used in a browser." from
`README.md`.

## 5. Zenodo DOI (before JOSS)

- Log in to https://zenodo.org with GitHub and switch on the `ggc-kinetics` repository under GitHub integration.
- Make a GitHub release from tag `v0.2.0`. Zenodo archives it and mints a DOI.
- Add the DOI to `CITATION.cff` (`doi:`) and the README, and keep it for the JOSS submission form.

## 6. Paper

- Submit `paper/paper.md` and `paper/paper.bib` to JOSS (https://joss.theoj.org/papers/new); it needs the public
  GitHub repository from step 1.
- Post the preprint to bioRxiv.
