"""
GGC-inf Release Phase 5 (v0.2.0-dev), item 2: PRE-REGISTERED fresh real-plate K59 search, pass 2 (tiered blanks).
PRE-REGISTRATION (this block was written and hashed BEFORE any Phase 5 web search or download; not edited after).
Every number carries an L7 tag (MEASURED / MODELED / "MODELED fit of MEASURED data" / ASSUMED).

HYPOTHESIS (unchanged from Phase 4): on a fresh real plate the frozen K59 statistic
  ggc_kinetics.detectors.z_rao_orth_robust(..., tau=0.056) scores active-enzyme units above 1.768 and matched
  no-enzyme / blank units below 1.768 (ASSUMED values, calibrated on MODELED plates only).
FROZEN (zero re-tuning): K59 tau 0.056, threshold 1.768; row-0 base Z_Rao threshold 4.60; a_plate 0; parser
  substrate-mode normalisation exactly as in benchmark_fresh_real_k59.py (Phase 4). Gates G1-G3 unchanged:
  G1 parse status OK or a loud documented error code; G2 >= 80% active units z_K59 > 1.768 AND <= 10% blank units
  z_K59 > 1.768; G3 <= 5.0 s per plate (MEASURED).
WHAT CHANGES vs Phase 4: only the blank-eligibility rule is widened into two tiers, each logged:
  Tier 1 SAME_PLATE: matched no-enzyme / blank wells on the same plate as the active wells.
  Tier 2 SAME_STUDY_SEPARATE_FILE: blanks from the same study and the same instrument, in a separate file or sheet
    tab. A Tier 2 result is reported as Tier 2 and never relabelled as Tier 1.
ELIGIBILITY (all required, from metadata / README / licence BEFORE any download): CC BY or CC0; textbook-safe
  enzyme and assay within SAFETY_BOUNDS.md (beta-galactosidase / ONPG, alkaline phosphatase / pNPP, peroxidase /
  ABTS, lipase, esterase, urease, catalase); not Lagerman; not PMC5149054; raw time-resolved progress curves at
  >= 4 S0 levels; Tier 1 or Tier 2 blanks; never used in GGC-inf; a text / CSV / XLSX export a script can read.
EXCLUDED (never re-screened): darus-3867, darus-2467, ICEKAT PMC7222511 files, HRP PMC7666214 source data,
  Michaelis-Menten 1913, PON1 PMC6896510, and the six Phase 4 candidates (zenodo 7632173, 8250616, 5852254,
  21627132; the Phase 4 EnzymeML laccase hits; renz BMC Bioinf. simulated curves, Agilent ONPG note, FAIRDOMHub 923).
SEARCH QUERIES (frozen, run in this order; screening stops at 6 NEW candidates):
  R1 DaRUS / Dataverse: "enzyme kinetics absorbance time course" (Phase 4 Q6, never run)
  R2 Figshare: "beta-galactosidase ONPG kinetics raw data" OR "alkaline phosphatase pNPP kinetics plate reader"
  R3 Dryad: "enzyme kinetics plate reader progress curves"
  R4 Figshare / Zenodo: "lipase OR esterase p-nitrophenyl kinetics raw absorbance"
  R5 GitHub / EnzymeML: "EnzymeML omex" absorbance progress curves (urease, catalase, peroxidase, esterase)
  R6 GitHub paper supplements: "enzyme kinetics" "plate reader" raw data csv blank
STOP RULE: if no candidate qualifies, verdict NO_QUALIFYING_FRESH_PLATE_FOUND_PASS2. No simulated substitute, no
  reuse of excluded data. Table in benchmarks/fresh_k59_pass2_screening.json and REAL_DATA_BENCHMARKS.md.
  If one qualifies: download into its own new empty directory, write LAYOUT_JSON from the dataset's own metadata
  before any z is computed, run this script twice under python -I, diff, append the verdict (with its tier).
No randomness (no seed). Run twice (run1, run2) and diff; only tag and wall-time lines may differ.

Run:  python -I -X utf8 benchmarks/benchmark_fresh_real_k59_pass2.py TAG SCREENING_JSON [PLATE_FILE LAYOUT_JSON]
  Each candidate in SCREENING_JSON carries "tier" (SAME_PLATE / SAME_STUDY_SEPARATE_FILE / NONE).
Writes ../ggc_inf_kit/logs/phase5_fresh_k59_pass2_TAG.txt.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import benchmark_fresh_real_k59 as p4  # noqa: E402  (Phase 4 run_plate reused verbatim: same gates, same parser)

LOGS = p4.LOGS


def main():
    tag, scr = sys.argv[1], sys.argv[2]
    t0 = time.perf_counter()
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE PHASE 5: pre-registered fresh real-plate K59 search, pass 2 (tiered blanks) | {tag}")
    A("=" * 78)
    A("  frozen: K59 tau 0.056 thr 1.768; base thr 4.60; a = 0 (ASSUMED). Gates G1 parse, G2 >= 80% active and "
      "<= 10% blank above 1.768, G3 <= 5.0 s/plate (ASSUMED). Blank tiers: SAME_PLATE, SAME_STUDY_SEPARATE_FILE")
    table = json.load(open(scr, encoding="utf-8"))
    A(f"  screened NEW candidates: {len(table['candidates'])} (limit 6). Queries R1-R6 as pre-registered.")
    for i, c in enumerate(table["candidates"], 1):
        A(f"  C{i} [{c['query']}] {c['id']}")
        A(f"     url {c['url']}")
        A(f"     licence {c['licence']}; enzyme/assay {c['enzyme']}; safety {c['safety']}; blank tier {c['tier']}")
        A(f"     decision {c['decision']}: {c['reason']}")
    q = [c for c in table["candidates"] if c["decision"] == "ACCEPT"]
    if len(sys.argv) > 4 and q:
        A(f"  running accepted candidate {q[0]['id']} (blank tier {q[0]['tier']})")
        layout = json.load(open(sys.argv[4], encoding="utf-8"))
        verdict, gates = p4.run_plate(A, sys.argv[3], layout)
        verdict = f"{verdict}_TIER_{q[0]['tier']}"
    else:
        verdict, gates = "NO_QUALIFYING_FRESH_PLATE_FOUND_PASS2", [("G1", None), ("G2", None), ("G3", None)]
        A("  STOP: no candidate met every eligibility rule (either tier); no download, no plate run, no simulated")
        A("        substitute, no reuse of excluded datasets.")
    A("  gates: " + ", ".join(f"{k} {'NOT RUN' if v is None else ('PASS' if v else 'FAIL')}" for k, v in gates))
    A(f"  VERDICT: {verdict}")
    A(f"  wall time {time.perf_counter() - t0:.2f} s (MEASURED)")
    os.makedirs(LOGS, exist_ok=True)
    open(os.path.join(LOGS, f"phase5_fresh_k59_pass2_{tag}.txt"), "w", encoding="utf-8",
         newline="\n").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
