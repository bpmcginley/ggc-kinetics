# ggc-kinetics 0.2.0

Exact enzyme progress-curve solvers that need no ODE integration, two detection statistics, and a plate-reader
parser. Every function is copied without changes to the math from the GGC-inf research archive (66 cycles) and
checked against the archive by `release_parity_check.py`: 20/20 checks, 0.0 difference for most.

**Read this first.** Every number in this project carries a tag: MEASURED (from real data), MODELED (from
simulation) or ASSUMED (a chosen value). The solvers are exact numerics. The **thresholds and validity limits come
from MODELED studies only**. They have not yet been confirmed on a real plate. On real published data so far
(see `benchmarks/REAL_DATA_BENCHMARKS.md`), 3 of 5 pre-registered gates passed.

## Install and run

```
pip install .            # numpy, scipy
pip install ".[app]"     # adds streamlit, matplotlib for the demo app
ggc-fit plate.csv --evidence MEASURED --out reports/
python -m unittest discover -s tests -v              # from the repo root, package installed; < 1 s
python examples/quickstart.py                         # 1913 invertase fit + one simulated plate; < 1 s
```

`tests/data/plate_healthy_biotek_wide.csv` is a MODELED healthy plate that gives `P1: KILL` (replicate |ln Km ratio| 0.305 > 0.277): it is one draw from the false-reject tail of `P1` v2, whose healthy-plate power is 0.86 on 300 fresh MODELED plates, so about 14% of healthy plates fail.
It is kept on purpose, not swapped for a passing draw, to show that the replicate check fails loudly.

`ggc-fit FILE [--layout L.json] [--mode substrate|product] [--evidence MEASURED|MODELED] [--out DIR]` writes
`FILE.report.json` and `FILE.report.md`. It exits 0 on OK and 2 on a parse error (for example `NO_WELL_TABLE`).

## What is inside

| Function | What it does | Limits and evidence |
|---|---|---|
| `kp_kd_exact`, `fit_kp_vs_kd` | Lambert-W progress curve with product inhibition (Kp) and enzyme decay (kd); BIC choice between them | 1913 invertase: Km 16.52 mM (+0.11%), Vmax/Km 0.0461 (+0.32%), dBIC −19.29 (MODELED fit of MEASURED data) |
| `unified_kp_strang` | Kp + kd + non-enzymatic hydrolysis by Strang splitting | exact without hydrolysis; < 0.005 vs RK4 with it (MODELED) |
| `coupled_lag_decay` | coupled-assay lag with a decaying primary rate | < 1e-11 vs RK4 |
| `slow_binding_warped` | slow-onset inhibition by time warping | < 1e-8 vs 2-ODE RK4 |
| `morrison_exact` | tight binding, E0 comparable to Km | residual ≤ 1e-10. **E0 must be known to ±10%** (K66, MODELED) |
| `k57_map_fit` | substrate-inhibition test with a per-well pipetting prior | called = dBIC > −2.5 and 1/Ksi > 0.045 (ASSUMED). **Valid for A_peak ≤ 2.00, pipetting CV ≤ 5%, and `sig` close to the plate's real readout SD** (all MODELED; see below). Opt-in `dispersion_guard=True` (v0.2.0). Jacobian FD error ≤ 1e-6, solver residual ≤ 1e-8 |
| `z_rao_orth_robust` | saturation-curvature score for a whole plate | K59: tau 0.056, threshold 1.768; baseline threshold 4.60 (ASSUMED). Needs several S0 levels on one plate: one well gives 0 |
| `z_dsf_asym_boltz` | melt-transition score for DSF curves | call at ≥ 5.72 (K60, ASSUMED); fixed Tm and width |
| `parse_file`, `process_bytes` | CSV / TSV / text exports, long or wide, decimal commas, SoftMax blocks; Gen5 and Tecan headers skipped | no .xlsx. Errors come back as codes, not exceptions |

## K57 noise scale (new in 0.2.0)

`k57_map_fit` treats the readout SD you pass as `sig` as known. Its false-alarm bound (Wilson 90% upper ≤ 0.060,
MODELED) holds only if `sig` is close to the plate's real readout SD. In the pre-registered Q65 study (7 arms × 450
fresh MODELED plates, thresholds frozen), the nominal `sig` gave false-call w90 0.034 at the nominal noise floor,
0.076 at a 2x floor and 0.131 at a 3x floor (all MODELED). Do not use the default `sig` without checking it.

Two ways to fix the noise scale, neither of which re-tunes a threshold:

1. **Pass an estimated `sig` (preferred when blanks exist).** `sig` is a 1-D array with one SD per well (per S0
   level). Estimate it from ≥ 3 noise-only blank or replicate reads per S0 level, for example
   `sig = 0.02 * first_read + floor_hat`, and call `k57_map_fit(Y, s0n, t_grid, sig)`. In Q65 (variant V1) this
   gave false-call w90 0.020-0.037 across arms and recall 0.991 / 0.982 (MODELED). **Caveat:** V1 used the
   generator's own noise form (a 2% proportional term plus a floor). The 2% term is ASSUMED, which favours V1; on a
   real plate the proportional part must be estimated too.
2. **`dispersion_guard=True` (fallback with no blanks).** The null fit's chi-square sets
   `birge = max(1, sqrt(chi2_null / dof_null))`, with `dof_null = max(W*M - 4 - W, 1)`, and each plate with
   `birge > 1` is refitted at `sig * birge`. The result carries `"birge"` (all 1.0 when the guard is off). In Q65
   (variant V2) this passed every arm, but with a thin margin at the 2x floor: w90 0.0555 against the 0.060 gate
   (MODELED). The guard is one-sided: it never shrinks `sig`. The default `dispersion_guard=False` gives output
   bit-identical to 0.1.0.

Neither fix has been checked on a real plate.

## Not claimed

- The detectors and K57 rule have not yet been shown to work on a fresh real plate. On one real, **non-fresh**
  plate (darus-3867 SLAC / ABTS, Release Phase 2) Z_Rao scored enzyme rows 3.48-3.82 and no-enzyme rows 0.32-0.50
  (MODELED fit of MEASURED data): the controls stayed below 4.60, but the enzyme rows did not reach it (G1 FAIL,
  2/3 gates). The K59 statistic put the enzyme rows at 2.26-2.61 and the controls at 0.05-0.07 against 1.768.
- **K59 has not been validated on a fresh real plate with matched blanks.** K59 (tau 0.056, threshold 1.768,
  ASSUMED) is calibrated on MODELED plates and checked only on the non-fresh darus-3867 plate above. Two
  pre-registered searches screened 12 public candidates and none qualified, so no fresh plate was run
  (`benchmarks/REAL_DATA_BENCHMARKS.md`). Public deposits mostly give initial rates or processed values, not raw
  multi-S0 progress curves with blanks.
- `ggc-fit` reports the parse and the frozen contract tests; it does not fit every model listed above.
- The Streamlit app (`app/streamlit_app.py`) passed a headless `AppTest` smoke run (streamlit 1.67.0) on one
  MODELED plate, with the upload widget stubbed because AppTest cannot drive it. It was also checked once by hand in
  Chrome (20 November 2026): upload of the healthy plate, MODELED, scale 0.05, offset 0.040; title, report and trace
  plot rendered.

MIT licence. Data sources: `ATTRIBUTION.md`.
