"""
GGC-inf Release Phase 6: build tests/data/k57_guard_ref.json BEFORE solvers.py is edited.
Run with a venv that has ONLY the released v0.1.0 wheel installed, so the stored outputs are v0.1.0 outputs:
  python -I -X utf8 benchmarks/make_k57_guard_ref.py
Contents (all inputs MODELED, crc32 seeds; outputs are MEASURED numerics of the frozen rule):
  v010: k57_map_fit outputs (dbic, iksi, km, called) of v0.1.0 on the unit-test 11-well design, tags
        pkg_k57_guard_null (6 null plates) and pkg_k57_guard_pos (4 Ksi = 5 Km plates). Plates are regenerated in the
        test with the same generator (tests.gen_k57), so only outputs are stored.
  q65:  3 high-noise null plates from the Q65 generator, arm C3 (floor 0.0225 = 3x), tag q65_C3_unittest_ref (fresh),
        9-well Q65 design. Stores Y, the Q65 birge_null reference, Row 0 dbic/called, and the Q65 V2 result
        (k57_map_fit at SIG x birge, v0.1.0). Plates picked: the first 2 Row 0 false calls and the first non-call.
"""
import importlib.util
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.dirname(HERE)
import ggc_kinetics as gk  # noqa: E402

assert gk.__version__ == "0.1.0" and "site-packages" in gk.__file__, (gk.__version__, gk.__file__)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


tests = load("t_ggc", os.path.join(REL, "tests", "test_ggc_kinetics.py"))
q65 = load("q65", os.path.join(HERE, "q65_k57_dispersion_guard.py"))


def out(r):
    return {k: [float(x) for x in np.asarray(r[k], float)] for k in ("dbic", "iksi", "km")} | \
        {"called": [bool(x) for x in r["called"]]}


ref = {"version_used": gk.__version__, "v010": {}}
for tag, B, mode in (("pkg_k57_guard_null", 6, 0), ("pkg_k57_guard_pos", 4, 1)):
    Y = tests.gen_k57(tag, B, mode)[0]
    s0n = np.broadcast_to(tests.S0_NOM[None, :], (B, 11)).copy()
    ref["v010"][tag] = out(gk.k57_map_fit(Y, s0n, tests.T_GRID, tests.SIG))

Y, _, _ = q65.gen("C3", "q65_C3_unittest_ref", n=40)
nom = np.broadcast_to(q65.S0_NOM[None, :], Y.shape[:2]).copy()
r0 = gk.k57_map_fit(Y, nom, q65.T_GRID, q65.SIG)
c = np.asarray(r0["called"], bool)
pick = [int(i) for i in np.flatnonzero(c)[:2]] + [int(np.flatnonzero(~c)[0])]
Yp, nomp = Y[pick], nom[pick]
bg = q65.birge_null(Yp, nomp)
v2 = [out(gk.k57_map_fit(Yp[i:i + 1], nomp[i:i + 1], q65.T_GRID, q65.SIG * bg[i])) for i in range(len(pick))]
ref["q65"] = {"tag": "q65_C3_unittest_ref", "picked": pick, "s0_nom": [float(x) for x in q65.S0_NOM],
              "Y": Yp.tolist(), "birge": [float(x) for x in bg],
              "row0": out(gk.k57_map_fit(Yp, nomp, q65.T_GRID, q65.SIG)),
              "v2": {k: [v[k][0] for v in v2] for k in ("dbic", "iksi", "km", "called")}}
p = os.path.join(REL, "tests", "data", "k57_guard_ref.json")
with open(p, "w", encoding="utf-8", newline="\n") as fh:
    json.dump(ref, fh, indent=1)
print("wrote", os.path.relpath(p, REL), "picked", pick, "birge", [round(x, 4) for x in bg],
      "row0 called", ref["q65"]["row0"]["called"], "v2 called", ref["q65"]["v2"]["called"])
