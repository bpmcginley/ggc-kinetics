# GGC-inf real-data benchmarks (release v0.1)

Every number here comes from public, openly licensed, published data, or from a fit to that data with frozen
operators (no re-tuning). Tags (L7, as in PLAYBOOK.md): **MEASURED** = the published data themselves;
**MODELED fit of MEASURED data** = any number a frozen model computes from those data (Km, Vmax, dBIC, z, IC50);
**ASSUMED** = a choice made by us (gates, thresholds, units). No number here comes from simulated data.
Source logs are named in the last column.

| Dataset (licence) | Enzyme / assay | What was tested | Result | Verdict | Log |
|---|---|---|---|---|---|
| PMC3381512, Michaelis & Menten 1913 (public domain) | Invertase, sucrose | G2_Kp vs G1_kd fit, frozen Lambert-W solver | Km 16.52 mM (+0.11% vs reference), Eff 0.0461 (+0.32%), dBIC -19.29 (MODELED fit of MEASURED data) | PASS (suite R1) | sprint15 suite (R1/R2) |
| PMC6896510 (CC BY) | PON1 | Integrated MM fit vs GraphPad / DynaFit | max error 0.69% vs GraphPad LW; DynaFit 17.2% (MEASURED: published values) | PASS (suite R2) | sprint15 suite (R1/R2) |
| DaRUS darus-3867 (CC BY 4.0) | SLAC laccase, ABTS | Two wavelengths; SoftMax plate through plate_parser | Km 27.733 uM (340 nm), 42.756 uM (420 nm); SoftMax Km 34.28 uM; replicate spread 0.248-0.281 dex (MODELED fit of MEASURED data) | parser PASS (R15/R16); P1 KILL on that plate | cycle 40 logs |
| DaRUS darus-2467 (CC BY 4.0) | Trametes laccase, ABTS | Per-replicate W = 4 Km | 381.9 / 271.6 / 1076.5 uM; spread 0.598 dex (MODELED fit of MEASURED data) | K86 (spread too large to call one Km) | cycle51_run1, cycle52_run1/2 |
| ICEKAT Additional file 3 (CC BY 4.0) | not stated (sample file; paper's case study is SIRT1) | Global Schnell-Mendoza fit, 8 traces x 76 reads | Km 7.82 uM, Vmax 0.349 uM/s, RMSE 0.0029 (MODELED fit of MEASURED data); W = 8 initial-rate Km 9.64 uM = 1.23x higher, as predicted | **G3 FAIL**: max \|F\| 3.6e-7 > 1e-8 (see note 1) | release_supp_bench_run1/2 |
| ICEKAT GitHub test.csv vs file 3 | same data | Exact identity | max \|diff\| 2.2e-16; file 3 lacks one cell (500 uM at 600 s) | **G2 FAIL** on the strict "= 0" gate (note 2) | same |
| ICEKAT Additional file 4 (CC BY 4.0) | not stated | W = 8 rates, 4PL IC50 | IC50 0.853 (header units), Hill 0.18 (MODELED fit of MEASURED data; poor shape, no comparator) | report only | same |
| ICEKAT Additional file 5 (CC BY 4.0) | not stated (HTS sample) | plate_parser wide parse; per-well slopes | 96 wells, 121 reads, duplicate D2 flagged "first kept", 0 non-numeric; 6 wells \|z\| > 3 (B8, B9, D8-D11) | parse PASS; no control map, so no Z' | same |
| PMC7666214 Source Data (CC BY 4.0) | Horseradish peroxidase, ABTS + H2O2 | Plain MM per series vs the authors' own apparent Km | 10/12 within 10% (9 match to < 0.01%); misses 0.5 mM C +20.1%, 1 mM A -18.4%; implied [HRP] = Vmax/kcat = 2.25 nM in 9/12 series | **G4 PASS** (predicted >= 11/12: missed) | same |

Release benchmark survival: **3/5** pre-registered gates (G1 parse, G4 HRP, G5 runtime 0.1 s). Determinism: run1 and
run2 logs differ only in the tag and wall-time lines.

## Notes

1. **G3 |F| failure is a conditioning effect, not a solver error** (post-hoc diagnosis, does not change the verdict).
   The largest |F| is at S/S0 = 2e-9, deep in depletion, where F = Km ln(S0/S) + (S0 - S) - Vmax t is
   ill-conditioned: an S error of 2.3e-15 uM gives |F| 3.6e-7. Over points with S/S0 > 1e-3, max |F| is 8.0e-13.
   A residual gate for real data should be stated relative to Vmax t, or only on points above a depletion floor.
   That would be a new rule for a later cycle; it was not applied here.
2. **G2 failure is one floating-point unit.** The two files print the same numbers with different digits, so one
   cell parses 2.2e-16 apart. The pre-registered gate demanded exact equality and so fails.
3. **HRP misses.** Dropping the 200 uM H2O2 point reproduces the authors' 1 mM series-A Km to 0.1%, so the authors
   probably left it out of that fit. No single dropped point explains 0.5 mM C or 2 mM C.
4. **Z_Rao_Orth_Robust on file 5 is 0 by construction.** With one well per plate the curvature direction is
   collinear with the first-order direction and is projected out. On file 3, treated as one plate of 8 S0 levels,
   z = 15.33 (above the 4.60 baseline). That is a post-hoc report with no prediction registered beforehand.
5. **Not fitted:** PMC7666214 Suppl Figure 10 (ascorbate lag, which is stoichiometric and outside every frozen model),
   and every xlsx file through plate_parser (it reads text only; the xlsx was read with a stdlib reader).
6. **Screen.** See `logs/release_real_screen.txt`. Both new datasets are IN BOUNDS. The ICEKAT CSVs were fetched
   before the written screen was finished, and the log says so.

## Release Phase 2: Z_Rao_Orth_Robust on a real plate, enzyme rows vs no-enzyme rows

Script `ggc_inf_kit/benchmark_real_zrao_darus3867.py` (pre-registered in its docstring before run1), logs
`release_zrao_darus3867_run1/2.txt` (identical apart from the tag and wall-time lines). Data: darus-3867 SLAC / ABTS
plate, 340 nm, 9 S0 levels (5-200 uM) per row, 11 reads over 15 min (MEASURED). **L9: this plate is NOT FRESH**; its
Km was inspected in Cycles 39/40, though Z_Rao had never been computed on it. Frozen detector, zero re-tuning, a = 0.

| Row | Type | z_base (thr 4.60) | z_K59 (thr 1.768) | Mean conversion |
|---|---|---|---|---|
| A | enzyme | 3.480 | 2.256 | 14.6% |
| B | enzyme | 3.824 | 2.610 | 14.1% |
| C | enzyme | 3.678 | 2.481 | 13.9% |
| D | no enzyme | 0.496 | 0.070 | 0.4% |
| E | no enzyme | 0.502 | 0.071 | 0.2% |
| F | no enzyme | 0.315 | 0.045 | 0.4% |

z values are MODELED fits of MEASURED data; conversion is MEASURED; thresholds are ASSUMED. Release function vs
archive `plate_parser.z_rao_orthogonal_robust`: max |diff| 4.4e-16.

Gates: **G1 FAIL** (all enzyme rows clear K59 1.768 but none clears the base 4.60), **G2 PASS** (all control rows
< 4.60, and also < 1.768), **G3 PASS** (0.02 s). Survival **2/3**.

Post-hoc reading (no verdict change): both statistics rank every enzyme row above every control row by a wide margin,
but at 14% conversion and 11 reads the base statistic falls short of its MODELED-calibrated 4.60. The K59 rule would
have separated this plate; that observation was not the pre-registered gate and is one non-fresh plate. Report-only
variant with nominal S0 / 50 as the S0 input: enzyme z_base 2.71-3.16, controls -2.93 to +3.01, so the first-read S0
convention matters on real data.

## Release Phase 4: fresh real-plate K59 benchmark (pre-registered): NO_QUALIFYING_FRESH_PLATE_FOUND

Pre-registration: `benchmarks/benchmark_fresh_real_k59.py`, docstring frozen and hashed before any search
(sha256 c871c5cd...c7128715). Eligibility: CC BY / CC0; textbook-safe enzyme; raw progress curves at >= 4 S0 levels
on one plate; matched no-enzyme or blank wells on the same plate; never used in GGC-inf; script-readable export.
Screening used metadata, README text and licences only. No dataset was downloaded.

| # | Query | Candidate | Licence | Enzyme / assay | Decision and reason |
|---|---|---|---|---|---|
| C1 | Q1 | zenodo 7632173 | CC BY 4.0 | not stated (mostly simulated) | REJECT: spectrometer export, not one plate; no matched blanks described; 135 MB zip needed to check |
| C2 | Q1 | zenodo 8250616 | CC BY 4.0 | chitinase (mAMCase) | REJECT: pH series at one substrate level; no blanks described |
| C3 | Q2 | zenodo 5852254 | CC BY 4.0 | beta-galactosidase | REJECT: paper PDF only, no raw data |
| C4 | Q1 | zenodo 21627132 | CC BY 4.0 | LinB dehalogenase | REJECT: one S0 per plate; buffer control, not no-enzyme; not a listed textbook enzyme |
| C5 | Q4 | EnzymeML laccase/ABTS hits | CC BY 4.0 | laccase / ABTS | REJECT: only hits are darus-3867 and darus-2467 (excluded, already used) |
| C6 | Q5 | PMC OA hits (renz, Agilent note, FAIRDOMHub 923) | mixed | beta-galactosidase / ONPG | REJECT: simulated curves or no deposited data |

Screening stopped at the pre-registered limit of 6, so Q6 was not run. Verdict: **NO_QUALIFYING_FRESH_PLATE_FOUND**.
Gates G1-G3 were NOT RUN. Nothing simulated was substituted and no excluded dataset was reused. The K59 rule
(tau 0.056, thr 1.768, ASSUMED) is still validated on MODELED plates plus one non-fresh real plate (darus-3867), and
the depletion-fit check was not run because no plate qualified. Logs: `logs/phase4_fresh_k59_run{1,2}.txt`
(identical except for the tag line).

## Release Phase 5: fresh real-plate K59 search, pass 2 (tiered blanks, pre-registered): NO_QUALIFYING_FRESH_PLATE_FOUND_PASS2

Pre-registration: `benchmarks/benchmark_fresh_real_k59_pass2.py`, hashed before any Phase 5 search
(sha256 65cb02d2...050ba12ac). Frozen: K59 tau 0.056, thr 1.768; base 4.60; gates G1-G3 as in Phase 4 (all ASSUMED).
New here: blanks may be Tier 1 SAME_PLATE or Tier 2 SAME_STUDY_SEPARATE_FILE (same study and instrument, separate
file or tab; logged as Tier 2). The six Phase 4 candidates and every earlier exclusion were not re-screened.
Screening used metadata and article pages only. No dataset was downloaded.

| # | Query | Candidate | Licence | Enzyme / assay | Blank tier | Decision and reason |
|---|---|---|---|---|---|---|
| C1 | R1 | DaRUS darus-2469 | CC BY 4.0 | glucose oxidase | NONE | REJECT: enzyme not in scope; deactivation series at one substrate level |
| C2 | R1 | DaRUS darus-2470 | CC BY 4.0 | AHAS carboligation | NONE | REJECT: enzyme not in scope; pH series |
| C3 | R2 | PMC4908277 (Data in Brief 2016) | CC BY 4.0 | beta-gal / ONPG, whole cells | NONE | REJECT: deposited data are Miller units, not raw curves; one ONPG level |
| C4 | R4 | zenodo 17298898 | not confirmed | lipase / p-NP palmitate | NONE | REJECT: hypothetical (simulated) table; one time course |
| C5 | R4 | PMC3963902 (LpEst1 esterase) | CC BY | esterase / p-NP esters | NONE | REJECT: acyl-chain series, not an S0 ladder; no raw curves found |
| C6 | R5 | EnzymeML GitHub examples | not confirmed | glucose assay, carboligation | NONE | REJECT: no in-scope .omex absorbance dataset; only in-scope-adjacent hit is darus-3867 (excluded) |

R3 (Dryad) returned no dataset hit. R6 was not run because the limit of 6 was reached. Verdict:
**NO_QUALIFYING_FRESH_PLATE_FOUND_PASS2**. Gates G1-G3 NOT RUN, nothing simulated was substituted, and no excluded
dataset was reused. Widening the blank rule to Tier 2 did not help: no candidate had raw multi-S0 progress curves
for an in-scope enzyme. Logs: `logs/phase5_fresh_k59_pass2_run{1,2}.txt` (identical except for the tag line).
