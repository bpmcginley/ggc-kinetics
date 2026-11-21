"""
GGC-inf Release Phase 4 (v0.2.0-dev), item 1: PRE-REGISTERED fresh real-plate K59 benchmark.
PRE-REGISTRATION (this block was written and frozen BEFORE any web search or download in Phase 4; not edited after).
Every number carries an L7 tag (MEASURED / MODELED / "MODELED fit of MEASURED data" / ASSUMED).

HYPOTHESIS: on a fresh real plate (never downloaded or inspected in GGC-inf), the frozen K59 statistic
  ggc_kinetics.detectors.z_rao_orth_robust(..., tau=0.056) scores active-enzyme units above 1.768 and matched
  no-enzyme / blank units below 1.768 (both values ASSUMED, calibrated on MODELED plates only).
FROZEN PARAMETERS (zero re-tuning): K59 tau 0.056, threshold 1.768; row-0 baseline = base Z_Rao (tau None),
  threshold 4.60 (all ASSUMED). a_plate = 0 (no inner filter, ASSUMED). Normalisation = parser substrate mode:
  Y = (A - offset) / scale, scale = absorbance per 50 uM from the dataset's own standard curve or extinction
  coefficient (taken from the dataset or its paper, never fitted to the kinetic traces), t in minutes, first read
  = t 0, S0 input = first-read Y floored at 1e-4 (parser conventions, as in Release Phase 2).
  One "unit" (plate in the detector's sense) = the wells of one replicate series across all S0 levels.
DEPLETION-FIT CHECK (only if multi-S0 progress curves are present): Schnell-Mendoza closed form
  S = Km W((S0/Km) exp((S0 - Vmax t)/Km)) fitted per active unit (MODELED fit of MEASURED data); report raw
  max |F| and relative max |F|/(Vmax t) over points with S/S0 > 1e-3, both logged; no gate is moved.
GATES (ASSUMED):
  G1 Schema/Parse: ggc_kinetics.parser.parse_file reads the file with no custom hacks (status OK), or fails loudly
     with a documented error code (then G2 is NOT RUN and the verdict is the error code).
  G2 K59 separation: >= 80% of active units have z_K59 > 1.768 AND <= 10% of matched no-enzyme / blank units have
     z_K59 > 1.768. (Base z vs 4.60 reported alongside as row 0, not gated.)
  G3 Runtime: <= 5.0 s per plate (MEASURED).
ELIGIBILITY (all required, screened from metadata / README / licence BEFORE any download):
  licence CC BY or CC0; safe textbook enzyme and assay (beta-galactosidase, alkaline phosphatase, laccase,
  peroxidase, lipase, esterase or similar colorimetric hydrolase / oxidase), within SAFETY_BOUNDS.md, not
  Lagerman, not PMC5149054; raw time-resolved progress curves (not only initial rates or fitted parameters) at
  >= 4 S0 levels on one plate; matched no-enzyme or blank wells on the same plate; never used in GGC-inf
  (excluded: darus-3867, darus-2467, ICEKAT PMC7222511 files, HRP PMC7666214 source data, Michaelis-Menten 1913,
  PON1 PMC6896510); a plain-text / CSV / XLSX-converted-to-CSV export that a script can read.
SEARCH QUERIES (frozen, run in this order; screening stops at 6 candidates):
  Q1 zenodo.org: "enzyme kinetics plate reader raw data progress curve blank"   (CC BY / CC0 filter)
  Q2 zenodo.org: "beta-galactosidase ONPG kinetics plate reader"
  Q3 zenodo.org: "alkaline phosphatase pNPP kinetics raw data"
  Q4 EnzymeML / STRENDA DB: "EnzymeML progress curve absorbance laccase OR peroxidase OR esterase"
  Q5 PMC open access: "progress curves" "no enzyme control" "source data" pNPP OR ONPG OR ABTS
  Q6 DaRUS / Dataverse: "enzyme kinetics absorbance time course" (excluding darus-3867 and darus-2467)
STOP RULE: if no candidate qualifies after screening up to 6, verdict NO_QUALIFYING_FRESH_PLATE_FOUND. No simulated
  data is substituted and no excluded dataset is reused. The screened table goes into the log and
  benchmarks/REAL_DATA_BENCHMARKS.md.
No randomness (no seed). Run twice (run1, run2) and diff; only tag and wall-time lines may differ.

Run:  python -I -X utf8 benchmarks/benchmark_fresh_real_k59.py TAG SCREENING_JSON [PLATE_FILE LAYOUT_JSON]
  SCREENING_JSON: the screened-candidate table (written during screening). LAYOUT_JSON (only if a plate qualifies,
  written from the dataset's own metadata before any z is computed): {"scale":..,"offset":..,"S0_uM":{well:..},
  "active_units":[[wells..],..],"blank_units":[[wells..],..]}.
Writes ../ggc_inf_kit/logs/phase4_fresh_k59_TAG.txt.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.dirname(HERE)
LOGS = os.path.join(os.path.dirname(REL), "ggc_inf_kit", "logs")
try:
    import ggc_kinetics  # noqa: F401
except ImportError:
    sys.path.insert(0, os.path.join(REL, "src"))
from ggc_kinetics import detectors, parser  # noqa: E402

THR_K59, THR_BASE = 1.768, 4.60


def run_plate(A, plate, layout):
    import numpy as np
    t0 = time.perf_counter()
    rep, parsed = parser.parse_file(plate, evidence="MEASURED", mode="substrate",
                                    scale=layout["scale"], offset=layout["offset"])
    if rep.get("status") != "OK" or parsed is None:
        code = rep.get("error", {}).get("code", "UNKNOWN")
        A(f"  G1 parse: FAIL loudly with code {code}; G2 NOT RUN")
        return f"PARSE_ERROR_{code}", [("G1", False), ("G2", None), ("G3", time.perf_counter() - t0 <= 5.0)]
    A("  G1 parse: status OK (parser substrate mode, no custom hacks) ... PASS")
    from ggc_kinetics._plate_parser import well_index
    tmin = parsed["times_s"] / 60.0
    Yall = (parsed["A"] - layout["offset"]) / layout["scale"]

    def score(units):
        out = []
        for u in units:
            Y = np.stack([Yall[well_index(w)] for w in u])[None]
            s0 = np.maximum(Y[:, :, 0], 1e-4)
            z0 = np.zeros_like(s0)
            out.append((float(detectors.z_rao_orth_robust(Y, s0, tmin, z0)[0]),
                        float(detectors.z_rao_orth_robust(Y, s0, tmin, z0, tau=detectors.Z_RAO_K59_TAU)[0])))
        return out
    act, blk = score(layout["active_units"]), score(layout["blank_units"])
    for lab, sc in (("active", act), ("blank", blk)):
        for i, (zb, zk) in enumerate(sc):
            A(f"  {lab} unit {i + 1}: z_base {zb:.3f} (vs 4.60), z_K59 {zk:.3f} (vs 1.768) (MODELED fit of MEASURED data)")
    fa = sum(zk > THR_K59 for _, zk in act) / max(len(act), 1)
    fb = sum(zk > THR_K59 for _, zk in blk) / max(len(blk), 1)
    g2 = fa >= 0.80 and fb <= 0.10
    dt = time.perf_counter() - t0
    A(f"  G2 active above 1.768: {fa:.3f} (need >= 0.80); blank above 1.768: {fb:.3f} (need <= 0.10) ... "
      f"{'PASS' if g2 else 'FAIL'}")
    return "RUN", [("G1", True), ("G2", g2), ("G3", dt <= 5.0)]


def main():
    tag, scr = sys.argv[1], sys.argv[2]
    t0 = time.perf_counter()
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE PHASE 4: pre-registered fresh real-plate K59 benchmark | {tag}")
    A("=" * 78)
    A("  frozen: K59 tau 0.056 thr 1.768; base thr 4.60; a = 0 (ASSUMED). Gates G1 parse, G2 >= 80% active and "
      "<= 10% blank above 1.768, G3 <= 5.0 s/plate (ASSUMED)")
    table = json.load(open(scr, encoding="utf-8"))
    A(f"  screened candidates: {len(table['candidates'])} (limit 6). Queries Q1-Q6 as pre-registered.")
    for i, c in enumerate(table["candidates"], 1):
        A(f"  C{i} [{c['query']}] {c['id']}")
        A(f"     url {c['url']}")
        A(f"     licence {c['licence']}; enzyme/assay {c['enzyme']}; safety {c['safety']}")
        A(f"     decision {c['decision']}: {c['reason']}")
    q = [c for c in table["candidates"] if c["decision"] == "ACCEPT"]
    if len(sys.argv) > 4 and q:
        layout = json.load(open(sys.argv[4], encoding="utf-8"))
        verdict, gates = run_plate(A, sys.argv[3], layout)
    else:
        verdict, gates = "NO_QUALIFYING_FRESH_PLATE_FOUND", [("G1", None), ("G2", None), ("G3", None)]
        A("  STOP: no candidate met every eligibility rule; no download, no plate run, no simulated substitute,")
        A("        no reuse of darus-3867 / darus-2467 / ICEKAT / HRP.")
    A("  gates: " + ", ".join(f"{k} {'NOT RUN' if v is None else ('PASS' if v else 'FAIL')}" for k, v in gates))
    A(f"  VERDICT: {verdict}")
    A(f"  wall time {time.perf_counter() - t0:.2f} s (MEASURED)")
    os.makedirs(LOGS, exist_ok=True)
    open(os.path.join(LOGS, f"phase4_fresh_k59_{tag}.txt"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
