# GGC-inf Release Phase 4 (v0.2.0-dev): summary

Tags: MEASURED, MODELED, "MODELED fit of MEASURED data", ASSUMED (L7). Every script ran twice; run1 and run2 logs
differ only in the tag and wall-time lines.

## 1. Fresh real-plate K59 benchmark: NO_QUALIFYING_FRESH_PLATE_FOUND

- Pre-registration `benchmark_fresh_real_k59.py` was hashed before any search (sha256 c871c5cd...c7128715).
- 6 candidates were screened from metadata only (table in `REAL_DATA_BENCHMARKS.md` and `fresh_k59_screening.json`).
  All 6 were rejected. Main reasons: no matched no-enzyme wells on the same plate, fewer than 4 S0 levels on one plate,
  no raw data, simulated data, or already used (darus-3867, darus-2467).
- Gates G1-G3 were NOT RUN. Nothing was downloaded, nothing simulated was substituted, nothing excluded was reused.
- So K59 (tau 0.056, threshold 1.768, ASSUMED) still rests on MODELED plates plus one non-fresh real plate.

## 2. Q64 confounder check of the frozen K57 MAP rule: FAIL (G2, arm C)

Independent mass-action ODE (no code shared with `ggc_kinetics.solvers`); `k57_map_fit` frozen (thr -2.5, iksi 0.045).
450 null plates per arm (3 seed blocks x 150), 1,350 total. Pre-registration sha256 a9e1da62...9b133352.

| Arm | Mechanism (MODELED, rates ASSUMED) | False calls | Rate | Wilson 90% upper | G2 (<= 0.060) |
|---|---|---|---|---|---|
| A | product inhibition, Kp = 3 Km | 3/450 | 0.0067 | 0.0166 | PASS |
| B | ES inactivated at 2 kd, free E at kd | 14/450 | 0.0311 | 0.0476 | PASS |
| C | MM + decay, noise floor 0.015 | 34/450 | 0.0756 | 0.0987 | FAIL |
| pooled | | 51/1350 | 0.0378 | 0.0473 | PASS |

- G1 numerics (MEASURED): Jacobian vs central FD max rel err 1.03e-11 (<= 1e-6); invariants max |F| 1.07e-13
  (<= 1e-8); tolerance halving max |d readout| 5.21e-09 (<= 1e-6). PASS.
- G3 (MEASURED): arm invocations 5.07-5.93 s each (< 30 s). PASS.
- L21: the number that moved is arm C's false-alarm rate, from 10/450 (Cycle 66, floor 0.0075) to 34/450 (floor
  0.015). The mechanism is noise mis-specification: the rule takes its readout SD as known
  (SIG = 0.015 x max(S0, 0.5)), and a noise floor 2x higher than that inflates dBIC on null plates. No threshold was
  re-tuned.
- My pre-registered predictions were wrong in two places. I expected B to be the most likely to fail and C to pass;
  the opposite happened.

## 3. Regression

- Archive suite `sprint15_cycle30_regression_suite.py`: 42/42 PASSED, twice (only the wall_time line differs).
- Unit tests: 14/14 OK, twice, run against a fresh venv with the rc1 wheel installed (sha256 b88157b2...ed9c5f).

## Caveats for the v0.2.0 README and paper

- The K59 detector has still not been tested on a fresh real plate with matched blanks. Say so plainly.
- The K57 rule's false-alarm guarantee needs the user's readout SD to be close to the assumed SIG. If the real
  noise floor is about 2x higher, the false-alarm rate rose to 7.6% (w90 0.099) in MODELED data. Advise estimating SD
  from blank or replicate wells, and state this as a validity condition next to CV <= 5% and A_peak <= 2.
- Product inhibition (Kp = 3 Km) and faster ES inactivation (2 kd) did not raise false alarms above the bound
  (MODELED).
