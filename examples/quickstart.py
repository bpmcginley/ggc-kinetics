"""ggc-kinetics quickstart (< 1 s). Run from the repo root: python examples/quickstart.py

1. Fit the 1913 Michaelis & Menten invertase data (MEASURED data; PMC3381512 Table 1) with the frozen Lambert-W
   solver: product inhibition (Kp) vs enzyme decay (kd), chosen by BIC. The fit is a MODELED fit of MEASURED data.
2. Parse a MODELED (simulated) plate export from tests/data and print the contract verdicts. P1 = KILL on this
   plate is the frozen archive result (stored Cycle 34 report: replicate Km |ln ratio| 0.305 > 0.277).
   The plate is a draw from P1 v2's false-reject tail (healthy power 0.86 on 300 MODELED plates, ~14% fail),
   kept on purpose instead of a passing draw, to show the replicate check failing loudly.
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
try:
    import ggc_kinetics as gk                              # installed package first
except ImportError:                                        # else a source checkout, not installed
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import ggc_kinetics as gk
from ggc_kinetics import parser                            # noqa: E402

# (S0 mM sucrose, t min, fraction converted) per series. Michaelis & Menten 1913, Table 1.
MM1913 = [
    (333.0, [7.0, 14.0, 26.0, 49.0, 75.0, 117.0, 1052.0], [0.0164, 0.0316, 0.0528, 0.0923, 0.1404, 0.2137, 0.9834]),
    (166.7, [8.0, 16.0, 28.0, 52.0, 82.0, 103.0], [0.0350, 0.0636, 0.1080, 0.1980, 0.3000, 0.3780]),
    (83.0, [49.5, 90.0, 125.0, 151.0, 208.0], [0.352, 0.575, 0.690, 0.766, 0.900]),
    (41.6, [10.25, 30.75, 61.75, 90.75, 112.70, 132.70, 154.70], [0.1147, 0.3722, 0.615, 0.747, 0.850, 0.925, 0.940]),
    (20.8, [17.0, 27.0, 38.0, 62.0, 95.0, 1372.0], [0.331, 0.452, 0.611, 0.736, 0.860, 0.990]),
]


def main():
    t0 = time.perf_counter()
    S0 = np.concatenate([np.full(len(t), s) for s, t, _ in MM1913])
    t = np.concatenate([t for _, t, _ in MM1913])
    f = np.concatenate([f for _, _, f in MM1913])
    r = gk.fit_kp_vs_kd(S0, t, f)
    print(f"invertase 1913 (MODELED fit of MEASURED data): Km {r['Km']:.2f} mM, Vmax/Km {r['eff']:.4f}, "
          f"dBIC(Kp - kd) {r['dbic']:.2f}")
    plate = os.path.join(ROOT, "tests", "data", "plate_healthy_biotek_wide.csv")
    # The simulated plate was written with 0.05 absorbance per 50 uM and a 0.040 background (ASSUMED inputs, Cycle 33).
    rep, _ = parser.parse_file(plate, evidence="MODELED", mode="substrate", scale=0.05, offset=0.040)
    print(f"plate {os.path.basename(plate)} (MODELED plate): status {rep['status']}")
    for name, test in sorted(rep.get("tests", {}).items()):
        print(f"  {name}: {test.get('verdict')}")
    print("  (P1 KILL is expected here: this MODELED plate is a kept false-reject tail draw, ~14% of healthy "
          "MODELED plates; it shows the replicate check failing loudly)")
    print(f"wall time {time.perf_counter() - t0:.2f} s")
    return 0 if rep["status"] == "OK" else 1


if __name__ == "__main__":
    sys.exit(main())
