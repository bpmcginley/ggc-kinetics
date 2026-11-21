# ggc-kinetics v0.2.0 release notes (DRAFT, Phase 5)

**Status (Phase 6):** merged into README.md, paper/paper.md and RELEASE_CHECKLIST.md for v0.2.0. V2 now ships as
`k57_map_fit(..., dispersion_guard=True)` (default False). The text below is the Phase 5 draft, kept unchanged.

Tags: MEASURED, MODELED, "MODELED fit of MEASURED data", ASSUMED (L7). This is a draft. Nothing is published, and
`src/ggc_kinetics/solvers.py` is unchanged from v0.1.0. The released v0.1.0 wheel and sdist are untouched.

## 1. New validity condition for the K57 substrate-inhibition rule (from Q64)

The frozen K57 MAP rule (thr -2.5, iksi 0.045, sp 0.05, all ASSUMED) treats the readout SD as known:
SIG = 0.015 x max(S0, 0.5) (ASSUMED). If the real noise floor is higher, the false-alarm rate rises:

| Noise floor | False calls (Row 0) | Wilson 90% upper | Source |
|---|---|---|---|
| 0.0075 (1x) | 9/450 | 0.034 | Q65 arm N0 (MODELED) |
| 0.015 (2x) | 34/450 (Q64), 25/450 (Q65) | 0.099 (Q64), 0.076 (Q65) | MODELED |
| 0.0225 (3x) | 47/450 | 0.131 | Q65 arm C3 (MODELED) |

README text to add next to "CV <= 5%" and "A_peak <= 2": *The K57 false-alarm bound (<= 6%, MODELED) holds only if
the `sig` you pass is close to your plate's real readout SD. Do not use the default SIG without checking it.
Estimate SD from blank or replicate wells (see 2).*

## 2. Q65 outcome: two noise-scale fixes pass, with no threshold re-tuned (MODELED)

Q65 was pre-registered and hashed before running (file sha256 f1680c43...41d35a). It used 7 arms x 450 fresh
plates from an independent mass-action ODE (fresh `q65_*` seeds). The rule's thresholds were frozen, and only the
`sig` argument changed.

| Variant | N0 | C2 | C3 | A | B | G2 (all w90 <= 0.060) | P1 recall | P2 recall | Q65 |
|---|---|---|---|---|---|---|---|---|---|
| Row 0 (nominal SIG) | 9 (0.034) | 25 (0.076) | 47 (0.131) | 0 (0.006) | 14 (0.048) | FAIL (C2, C3) | 0.980 | 0.964 | FAIL |
| V1 blank-estimated sig | 4 (0.020) | 10 (0.037) | 9 (0.034) | 0 (0.006) | 9 (0.034) | PASS | 0.991 | 0.982 | PASS |
| V2 Birge guard | 7 (0.028) | 17 (0.056) | 10 (0.037) | 0 (0.006) | 11 (0.040) | PASS | 0.980 | 0.951 | PASS |

Cells are false calls /450 (w90), all MODELED. G1 numerics passed (MEASURED): Jacobian FD 1.18e-11, |F| 1.17e-13,
tolerance halving 8.25e-09. G4 passed: every arm invocation took 11-25 s (MEASURED).

- **V1** (recommended if blanks exist): estimate the floor from >= 3 noise-only blank reads per S0 level, and pass
  sig[w] = 0.02 x first read + floor_hat. It had the best specificity and the best recall.
- **V2** (fallback with no blanks): fit the null at nominal SIG, then birge = max(1, sqrt(chi2_null / 104)), and refit
  at SIG x birge. It passes, but arm C2 is close to the limit (w90 0.0555 vs 0.060).
- **Caveat:** these are MODELED plates from one generator family. The noise model in V1 matches the generator's
  (2% proportional + floor), which favours V1. Neither fix has been checked on a real plate.
- **Not yet done:** solvers.py has not been changed. Per the rules, adding V1/V2 as an option waits for an audit of Q65.

## 3. K59 fresh real-plate search: still no qualifying plate

Pass 1 (Phase 4) and pass 2 (Phase 5, with Tier 2 separate-file blanks allowed) together screened 12 candidates.
None qualified, and no plate was run. Verdicts: NO_QUALIFYING_FRESH_PLATE_FOUND and
NO_QUALIFYING_FRESH_PLATE_FOUND_PASS2. README and paper must say plainly: *K59 (tau 0.056, thr 1.768, ASSUMED) is
calibrated on MODELED plates and checked on one non-fresh real plate (darus-3867). It has not been validated on a
fresh real plate with matched blanks.* The common gap: public deposits give initial rates or processed values,
not raw multi-S0 progress curves with blanks. A purpose-made open plate (wet lab, outside this project's scope) may be
the only route.

## 4. Paper changes

- Methods: add the noise-validity condition and the Q65 design, including the V1/V2 formulas as pre-registered.
- Results: add the Q65 table. Report the Row 0 failure at 2x and 3x floors as a real limitation, not a footnote.
- Limitations: no fresh real-plate validation for K59 or for V1/V2. All Q64/Q65 rates are MODELED.
