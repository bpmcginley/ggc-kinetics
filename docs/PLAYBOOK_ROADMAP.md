# ggc-kinetics: Playbook roadmap (v0.3 → v1.0)

Drafted 2026-11-20. Runs under `INNOVATION_PLAYBOOK_7.1.md` (Full + Recursive mode).
Nothing here is a result yet. Every number below is a **target** (ASSUMED) until a cycle measures it.

---

## The bet in one paragraph

Today a scientist runs a plate, exports it, picks "initial rate" windows by hand, fits Michaelis-Menten in a
spreadsheet or GraphPad, and only learns afterwards that the plate could not answer the question (Km outside
the substrate range, too few points, an artifact). ggc-kinetics already has fast closed-form solvers. The bet:
use that speed to **simulate thousands of virtual copies of the user's own assay** in seconds, so the tool can
(1) tell you *before* you pipette whether the plate design can answer your question, (2) fit the full progress
curves instead of hand-picked initial rates, (3) attach a calibrated false-call rate to every flag *for your
specific plate design*, and (4) tell you which plate to run next. Nobody does (1), (3) or (4) as a routine
step today. That is the "we didn't have access to this before" part. Speed is the enabler, not the pitch.

**North-star demo:** drop a raw plate file in a browser tab → in under 10 s get parameters with honest
intervals, flagged wells with their false-call rates, a "this design was / wasn't adequate" verdict, and a
recommended next plate. Nothing to install.

---

## 0. Project card

```
PROJECT:              ggc-kinetics, plate-to-decision enzyme kinetics with design-before-data
WHAT COUNTS AS A WIN: (a) on simulated plates with an independent truth model, Km/kcat interval coverage
                       within 90 +/- 3 % and median interval width >= 2x narrower than an initial-rate
                       MM fit on the same wells; (b) a recommended design that reaches the same Km
                       precision with >= 30 % fewer wells than a standard 2-fold series; (c) one real
                       public plate re-analysed where our result agrees with the published value within
                       its reported interval; (d) a timed study: >= 10x less wall time than the
                       user's own manual workflow (n >= 5 people)
REAL-WORLD CHECK:     public raw progress curves (DaRUS 10.18419/darus-3867, darus-2467; EnzymeML /
                       STRENDA DB; Zenodo/Dryad), old plate exports from UMass labs, timed-study volunteers
BUDGET PER CYCLE:     one to two weeks of evenings, $0, laptop CPU
FAST-LOOP LIMIT:      30 s per test (closed-form solvers make 10^4 virtual plates fit inside this)
TIER-1 ENGINE:        vectorised Lambert-W / closed-form virtual-plate generator (NumPy), MAP-Elites
                       over assay designs
VERIFIER COMMAND:     python -m unittest discover -s tests + run-twice log diff + fresh-seed calibration
BOUNDS:               no wet lab; no patient or human-subject data; no proprietary instrument files
                       without a licence; never claim real-data validation from simulation
EXCLUDED SOURCES:     any dataset used to tune a threshold may not be used to test it; ../ggc_inf_kit/
                       is read-only history and is not modified
CONTROLS:             null plates (no enzyme, blank-only, pure noise), plates with a planted wrong model,
                       plates whose Km sits outside the substrate range (must be called "not identifiable")
BORDERLINE CASES:     Km at the top substrate level; inhibitor Ki ~ [E]; 2 replicates per level
FP CAP:               6 % (Wilson upper 90 %)
SEED RULE:            crc32("ggck_c<cycle>_<set>")
ID PREFIXES:          GGCK-K#, GGCK-Q#, GGCK-V#, GGCK-P#   (archive IDs such as K57/Q65 stay as they are)
PROJECT TAG:          GGCK
```

---

## G0. Inventory (what is free right now)

| Kind | Items |
|---|---|
| Code | closed-form MM / tight-binding / inhibition solvers (`solvers.py`, 398 lines); detectors (`detectors.py`); 944-line multi-format plate parser; `k57_map_fit` with `dispersion_guard`; CLI; Streamlit app; 18 tests; CI; Trusted Publishing |
| Data | 3 synthetic plate fixtures (BioTek wide, Tecan transposed, long TSV); 2 public DaRUS sets already screened (K59 benchmarks) |
| Closed forms | Lambert-W progress curve (Schnell–Mendoza), Morrison tight-binding, Strang splitting for coupled cases |
| Known gaps | detectors 3/5 and 2/3 real-data gates; no external users; no timing data; no design tools |
| Open questions | Q64 (ODE confounders), Q65 follow-ups (noise scale) |

Thin spot, recorded as a finding: **only 2 public real-data sets screened**. Real data is the binding constraint,
so every cycle ends with a real-world falsifier step (G7) even when it is blocked.

---

## G1. Twelve goals

| # | GOAL (falsifiable) | FALSIFIER | COST | Kind |
|---|---|---|---|---|
| 1 | Initial-rate MM fits bias Km by > 10 % under realistic depletion and lag, and a global full-curve fit removes the bias | 2,000 virtual plates from an **independent** mass-action ODE truth | free | reported-number definition |
| 2 | A design score computed *before* the run predicts the realised Km interval width within +/- 20 % | design → simulate → fit → compare predicted vs realised | free | order of steps |
| 3 | A MAP-Elites design search finds plates that reach the same Km precision with >= 30 % fewer wells than a 2-fold series | compare best-found vs standard designs on fresh seeds | free | the object (design) |
| 4 | Per-plate calibration (simulate nulls on the user's own design) keeps every detector's FP <= 6 %, where the fixed global threshold exceeds it on >= 1 realistic design | sweep designs, count null calls | free | admission criterion |
| 5 | "Not identifiable" calls (Km outside range, flat likelihood) catch >= 95 % of planted unidentifiable plates with <= 5 % false alarms | planted controls | free | admission criterion |
| 6 | A next-plate recommender (expected information gain) separates competitive vs mixed inhibition in fewer plates than a standard matrix | simulated sequential campaigns | bounded | order of steps |
| 7 | Readout artifacts (inner filter, saturation, evaporation, temperature lag) are detectable from curve shape with FP <= 6 % | planted artifacts vs clean plates | free | measurement |
| 8 | Model selection by dBIC with a design-calibrated threshold picks the true mechanism >= 90 % of the time | planted mechanisms | free | criterion |
| 9 | The full pipeline runs in a browser (Pyodide) on a 96-well plate in < 10 s | build and time it | bounded | delivery |
| 10 | The parser auto-detects layout on >= 80 % of public plate exports with zero configuration | collect >= 20 public exports | bounded | input |
| 11 | EnzymeML / STRENDA export validates against the official schema | schema validator | free | output |
| 12 | Users finish a plate analysis >= 10x faster than with their own manual workflow | timed study, n >= 5 | bounded | real world |

## G2. Divergence filter (keep >= +3)

Scores: Vocabulary, Resources, Predicate (x2), Definition, Asymmetry, Negative result.

| # | V | R | P x2 | D | A | N | Total | Fate |
|---|---|---|---|---|---|---|---|---|
| 1 | +1 | +1 | +2 | +1 | +1 | +1 | **+7** | keep: everyone reports initial rates, almost nobody questions the window choice |
| 2 | +1 | +1 | +2 | +1 | +1 | +1 | **+7** | keep: design adequacy *before* data is a step nobody runs |
| 4 | +1 | +1 | +2 | 0 | +1 | +1 | **+6** | keep: thresholds calibrated to your plate, not a global constant |
| 5 | +1 | +1 | +2 | 0 | +1 | +1 | +6 | merge into 2 (same machinery) |
| 3 | 0 | +1 | 0 | 0 | +1 | +1 | +3 | run as Tier-1 engine inside cycle 2 |
| 8 | 0 | +1 | +2 | 0 | 0 | 0 | +3 | merge into 4 |
| 7 | 0 | +1 | 0 | 0 | +1 | +1 | +3 | red-team archive feed (stress cases) |
| 6 | 0 | 0 | 0 | 0 | +1 | +1 | +2 | deferred to cycle 5; needs 2 and 4 first |
| 11 | +1 | +1 | 0 | 0 | 0 | 0 | +2 | product chore, do alongside |
| 9 | -1 | 0 | 0 | 0 | 0 | 0 | -1 | product wedge, not a research goal: do it, but not as a cycle |
| 10 | 0 | -1 | 0 | 0 | 0 | +1 | 0 | blocked on data collection; run in background |
| 12 | 0 | -1 | 0 | 0 | 0 | +1 | 0 | the north-star check; needs volunteers |

The filter measures openness, not worth. 9 and 12 are where the "jaws drop" lives, so they stay as product
work and real-world checks, just not as research cycles.

## G3. Contradictions (survivors)

- **Goal 1:** precision vs bias. A longer rate window narrows the interval but adds depletion bias. Both are scored
  by the same measure, Km error, so only a structural change can work: fit the whole curve with the closed form.
- **Goal 2:** wells spent vs answer quality. The same measure scores both, which points to a structural change:
  move the check to *before* the run.
- **Goal 4:** sensitivity vs false calls. The same threshold scores both. The structural change is to make the
  threshold a function of the design, calibrated by simulation, not a constant.

---

## The engine (3-tier mapping, §3.5)

Build this once in cycle 0. Every later cycle reuses it.

| Tier | What it is here | Gate |
|---|---|---|
| **Tier 1** (10^6–10^8 evals) | Vectorised closed-form virtual-plate generator. MAP-Elites over designs (substrate levels, time grid, [E], replicates). Axes: **wells used** x **identifiability margin**. | L0 MODELED. It may nominate designs, never declare a win. |
| **Tier 2** (10^2–10^4) | Detector and design-rule code mutated as diffs, with failing test IDs attached. Islands: sensitivity-first, FP-first, speed-first. | Unit tests plus bit-exact parity with the slow reference, plus calibration on fresh seeds (L1) |
| **Tier 3** (5–20 per campaign) | **Independent truth:** a mass-action ODE (`solve_ivp`) that never uses the Lambert-W code, so the fit cannot just agree with itself. **Red-team archive:** inner filter, evaporation, temperature lag, photobleaching, product inhibition, pipetting CV, plate edge effects, noise learned from blank wells. | Fresh seeds plus the red-team archive plus one real plate (L2) |

**Run-twice rule.** Every simulation script writes a log. A second run must match it after the run tag and the
wall time are dropped.

---

## Cycle plan

Each cycle ends in a tagged release (v0.3, v0.4, …). That also builds the iterative public history JOSS asks for.

| Cycle | Title | Target release | Ships to users |
|---|---|---|---|
| **C0** | Virtual-assay engine + independent ODE truth + red-team library v1 | v0.3 (Dec 2026) | `ggc_kinetics.sim`: `simulate_plate(design, mechanism, noise, artifacts, seed)` |
| **C1** | Initial rates vs full-curve fit (goal 1) | v0.4 (Jan 2027) | `global_fit()` opt-in; report shows both, with intervals |
| **C2** | Design before data (goals 2, 3, 5) | v0.5 (Feb 2027) | `ggc-design`: "will this plate answer your question?" plus a designer |
| **C3** | Per-plate calibrated calls (goals 4, 8) | v0.6 (Mar 2027) | every flag carries "FP x % on your design, n sims" |
| **C4** | Browser wedge (goal 9) plus EnzymeML export (goal 11) | v0.7 (Apr 2027) | static Pyodide page: drop a file, get a report |
| **C5** | Real-world falsifiers: public-data re-analysis plus timed study (goals 10, 12, 6) | v1.0 (May 2027) | paper results; JOSS resubmission |

### C0 pre-registration (fill the date when it starts)

```
CYCLE C0 — Virtual-assay engine and independent truth
Pre-registered on <date>, before the first run.
FRESH?          yes
PROMOTION TIER: L1

CLAIM:          the ODE truth model and the closed-form solver agree on noise-free MM curves to
                relative 1e-6, and disagree as expected (> 1 %) once a planted non-MM mechanism is on
DATA:           simulated only; MODELED; screened against BOUNDS: yes
GATES (ASSUMED):
  G1 max rel. diff, closed form vs ODE, 200 noise-free MM designs <= 1e-6   prediction: SURVIVED
  G2 10^4 virtual 96-well plates generated in <= 30 s on the laptop         prediction: SURVIVED
  G3 every red-team artifact shifts at least one curve feature > 3 sigma    prediction: SURVIVED
NEGATIVE CONTROL: 200 blank plates; the generator must yield zero-slope curves within noise
CALIBRATION:    seed set ggck_c0_cal; on failure print GGCK_C0_CAL_FAIL
FRESH TEST:     seed set ggck_c0_fresh, never used for tuning
CLOSURE CONDITIONS: the ODE path imports nothing from solvers.py (checked by a test)
SELF-KILL:      fresh draws x3; stiff regimes ([S] >> Km, [E] ~ [S]); wrong noise model
RUNTIME LIMIT:  30 s per gate
SEEDS:          crc32("ggck_c0_cal"), crc32("ggck_c0_fresh")
RERUN:          run twice, diff logs
```

### C1 pre-registration (skeleton)

```
CLAIM:  on ODE-truth plates with realistic depletion (10–30 % substrate used) and a 0–60 s lag, a manual-style
        initial-rate MM fit has median |Km error| > 10 %, while global_fit has <= 3 % with 90 % interval
        coverage of 87–93 %
GATES:  G1 initial-rate median |dKm| > 10 %   G2 global median |dKm| <= 3 %   G3 coverage 87–93 %
        G4 the global fit adds <= 1 s per plate
NEGATIVE CONTROL: no-enzyme plates; global_fit must report "no activity", never a Km
BORDERLINE: Km at the top substrate level; it must be reported as a lower bound, not a number
SELF-KILL: noise learned from blanks; red-team archive v1; perturb depletion and lag ranges
KILL MEANS: if the initial-rate bias is <= 10 % on realistic settings, the headline is dead. Record it in
        KILLS and pivot C1 to "intervals only".
```

C2 and C3 get the same block, written before their first run, never after.

---

## G7. Real-world falsifiers (none needs a wet lab)

1. **Public raw curves.** Do a pre-registered search of EnzymeML / STRENDA DB, DaRUS, Zenodo and Dryad for raw
   progress curves with a published Km. Re-analyse each one blind, then compare. Start in C0 and run it in the
   background.
2. **Lab exports.** Email 5–10 UMass enzymology or biochemistry labs and ask for old raw plate exports they can
   share openly. A draft email is to be written; nothing is sent without Bruce's OK.
3. **Timed study.** Ask 5 or more students or lab members to analyse the same plate their usual way, then with
   the tool, and time both. This turns the 15–60 minute versus seconds idea into a MEASURED number.
4. **GraphPad / DynaFit parity.** Run the same public curves through the free or trial versions where licensing
   allows, and report agreement.

Until these run, every speed or accuracy claim stays MODELED, in the README as well as here.

---

## Meta step (between cycles)

- Append to KILLS / LEVERS / QUESTIONS / PRIMITIVES. These live in `docs/registries/` and are created in C0.
- Promote verified motifs to PRIMITIVES. Candidates: "closed-form global fit", "design-calibrated threshold".
- Add every stress case that broke a champion to the red-team archive.
- Pin each newly locked rule as a behaviour test, with exact counts on one fixed seed.
- Write a CHANGELOG entry plus a playbook changelog note, then end with the Gemini update prompt
  (see the GGC rule).
- **Stop rule:** two cycles in a row with no new kind of surviving goal.

## What would make this "industry standard"

These are the claims each cycle is trying to earn. None of them is earned yet.

1. "It told me my plate couldn't answer my question *before* I ran it." (C2)
2. "Its Km didn't depend on how I picked the initial-rate window, because it doesn't ask." (C1)
3. "Every flag told me how often it's wrong on a plate like mine." (C3)
4. "I dropped the file in a browser and had a STRENDA-ready report in seconds." (C4)
5. "It agreed with the published values on real data, and it saved our lab hours." (C5)

If C1's bias claim dies, the remaining ones (2, 3, 4) still stand on their own. That is why the cycles are
ordered this way.
