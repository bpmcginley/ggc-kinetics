"""
GGC-inf Release Phase 4 (v0.2.0-dev), item 2 (Q64): frozen K57 MAP rule vs three mass-action CONFOUNDER mechanisms.
PRE-REGISTERED: this docstring was written and hashed before the first run (L8). Simulation only (MODELED); no
downloads, no wet lab. Every number carries an L7 tag.

QUESTION (L21): the release rule ggc_kinetics.k57_map_fit (frozen: dBIC thr -2.5, iksi 0.045, prior sp 0.05) calls
  substrate inhibition. Cycle 66 (Q61) showed 10/450 false calls (w90 0.0368, MODELED) on plain MM + decay data from
  an independent ODE. The movable number here is the false-alarm rate when the null data come from mechanisms that
  also bend progress curves but contain NO substrate inhibition. A pass says the rule does not mistake these for
  substrate inhibition; a fail names the mechanism and the number that moved. No threshold is re-tuned on a fail.
GENERATOR (independent: zero code shared with ggc_kinetics.solvers; own RHS, Jacobian, IFE, noise, Wilson bound).
  Safe textbook mechanisms only (no organism, no synthesis). Species S, E, ES, P, EP, Ed. All rates ASSUMED.
  E + S <-> ES -> E + P: k_on 100, k_off 99, k_cat 1 (Km = (k_off + k_cat)/k_on = 1.000), E0 0.02.
  Per plate: rule-time scale lam = vm/(k_cat E0), vm ~ U(0.18, 0.45), multiplies every binding/catalytic rate;
  kd ~ U(0, 0.025) per rule time unit (not scaled). Pipetting S0 = S0_NOM x clip(N(1, 0.05), 0.7, 1.3).
  Readout = unconverted substrate S + ES through inner filter a = 1/16, plus N(0, sd), sd = 0.02 x signal + floor.
  Ladder A1.000: geomspace(0.2, 16, 9) (9 wells); 13 reads on 0..10 (ASSUMED, Cycle 66 design).
  ARM A product inhibition: E + P <-> EP, kp_on 100, kp_off 300 (Kp = 3.000 = 3 Km), E, ES, EP all decay at kd
        (EP -> Ed + P). Floor 0.0075.
  ARM B differential inactivation: free E -> Ed at kd, ES -> Ed + S at 2 kd. No EP. Floor 0.0075.
  ARM C MM + uniform decay at kd (E, ES), readout noise floor 0.015 (2x Cycle 66). No EP.
  Solver: scipy solve_ivp BDF, rtol 1e-8, atol 1e-12, block-diagonal sparse (bsr) Jacobian written by hand.
RULE INPUT (frozen): k57_map_fit(Y, nominal S0, t grid, sig) with sig = 0.015 x max(S0_NOM, 0.5) (the rule's own
  ASSUMED SIG, re-derived as constants here). Called iff dBIC > -2.5 AND iksi > 0.045.
SIZE: per arm 3 seed blocks x 150 null plates = 450; 1,350 plates total. Seeds zlib.crc32 of tags
  q64_{A,B,C}_null_d{0,1,2}, q64_{A,B,C}_tolcheck, q64_{A,B,C}_jaccheck. Never used before.
  Reading of the prompt (disclosed): "3 ladders x 150" is read as 3 seed blocks x 150 plates on the A1.000 design;
  "W = 8 of 12 reads" is unclear, so the Cycle 66 design (9 wells, 13 reads) is kept.
SPLIT: one invocation per arm (3 blocks in 3 worker processes; writes JSON), then one "pool" invocation that reads
  the three JSONs and applies the gates. Each invocation < 30 s.
GATES (ASSUMED rules):
  G1 numerics per arm: generator Jacobian vs central FD max rel err <= 1e-6; invariants S + ES + P + EP = S0 and
     E + ES + EP + Ed = E0 held to max |F| <= 1e-8 on every plate; tolerance halving (rtol, atol / 2) on 20 plates
     max |d readout| <= 1e-6 (MEASURED numerics).
  G2 false alarms: per arm Wilson 90% one-sided upper bound (z 1.6448536) w90 <= 0.060, AND pooled (1,350) w90 <= 0.060
     (MODELED).
  G3 each invocation wall time < 30 s (MEASURED).
MY PREDICTIONS (written before run1): A (product inhibition): high-S0 wells make more product and slow more; the decay
  term absorbs most of it; false calls ~ 10-25/450, w90 0.03-0.07 (borderline). B (ES inactivated 2x faster): high-S0
  wells hold more ES and lose activity faster, which is the closest mimic of substrate inhibition; w90 0.05-0.12,
  most likely to FAIL. C (noise floor 2x the rule's SIG at low wells): inflates chi2 at the low wells, roughly
  symmetric between models; w90 0.03-0.06. Pooled w90 ~ 0.05-0.07.
Run:  python -I -X utf8 benchmarks/q64_ode_confounders.py TAG A   (then B, C, then pool)
Writes ../ggc_inf_kit/logs/phase4_q64_{ARM|pool}_TAG.{txt,json}. Run twice (run1, run2) and diff; only tag and wall-time
lines may differ.
"""
import json
import math
import os
import sys
import time
import zlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.dirname(HERE)
LOGS = os.path.join(os.path.dirname(REL), "ggc_inf_kit", "logs")
try:
    import ggc_kinetics  # noqa: F401
except ImportError:
    sys.path.insert(0, os.path.join(REL, "src"))
from ggc_kinetics.solvers import k57_map_fit, K57_THR, K57_IKSI  # noqa: E402  (the frozen rule under test)

K_ON, K_OFF, K_CAT, E0 = 100.0, 99.0, 1.0, 0.02
KP_ON, KP_OFF = 100.0, 300.0
S0_NOM = np.geomspace(0.2, 16.0, 9)
T_GRID = np.linspace(0.0, 10.0, 13)
SIG = 0.015 * np.maximum(S0_NOM, 0.5)
A_IFE = 1.0 / 16.0
NB, NBLK = 150, 3
RTOL, ATOL = 1e-8, 1e-12
ZW = 1.6448536
ARMS = {"A": dict(pi=1.0, fes=1.0, floor=0.0075, name="product inhibition E + P <-> EP, Kp = 3 Km"),
        "B": dict(pi=0.0, fes=2.0, floor=0.0075, name="ES inactivated at 2 kd, free E at kd"),
        "C": dict(pi=0.0, fes=1.0, floor=0.015, name="MM + decay, noise floor 0.015")}


def seed(tag):
    return zlib.crc32(tag.encode("utf-8")) & 0xFFFFFFFF


def wilson_up(k, n, z=ZW):
    p = k / n
    den = 1.0 + z * z / n
    return (p + z * z / (2 * n) + z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / den


def ife(S, a):
    k = a * math.log(10.0)
    return -np.expm1(-k * S) / k


def rhs_jac(y, lam, kE, kES, kEP, pi):
    """y (M, 6) = S, E, ES, P, EP, Ed. Returns f (M, 6), J (M, 6, 6)."""
    S, E, ES, P, EP = y[:, 0], y[:, 1], y[:, 2], y[:, 3], y[:, 4]
    a, b, c = lam * K_ON, lam * K_OFF, lam * K_CAT
    pn, pf = lam * KP_ON * pi, lam * KP_OFF * pi
    f = np.empty_like(y)
    f[:, 0] = -a * E * S + b * ES + kES * ES
    f[:, 1] = -a * E * S + (b + c) * ES - kE * E - pn * E * P + pf * EP
    f[:, 2] = a * E * S - (b + c + kES) * ES
    f[:, 3] = c * ES - pn * E * P + pf * EP + kEP * EP
    f[:, 4] = pn * E * P - (pf + kEP) * EP
    f[:, 5] = kE * E + kES * ES + kEP * EP
    J = np.zeros((y.shape[0], 6, 6))
    J[:, 0, 0], J[:, 0, 1], J[:, 0, 2] = -a * E, -a * S, b + kES
    J[:, 1, 0], J[:, 1, 1], J[:, 1, 2], J[:, 1, 3], J[:, 1, 4] = -a * E, -a * S - kE - pn * P, b + c, -pn * E, pf
    J[:, 2, 0], J[:, 2, 1], J[:, 2, 2] = a * E, a * S, -(b + c + kES)
    J[:, 3, 1], J[:, 3, 2], J[:, 3, 3], J[:, 3, 4] = -pn * P, c, -pn * E, pf + kEP
    J[:, 4, 1], J[:, 4, 3], J[:, 4, 4] = pn * P, pn * E, -(pf + kEP)
    J[:, 5, 1], J[:, 5, 2], J[:, 5, 4] = kE, kES, kEP
    return f, J


def rates(arm, kd, W):
    p = ARMS[arm]
    kdm = np.repeat(kd, W)
    return kdm, p["fes"] * kdm, kdm * (1.0 if p["pi"] else 0.0), np.full(kdm.size, p["pi"])


def integrate(arm, S0, lam, kd, rtol=RTOL, atol=ATOL):
    from scipy.integrate import solve_ivp
    from scipy.sparse import bsr_matrix
    B, W = S0.shape
    M = B * W
    lm = np.repeat(lam, W)
    kE, kES, kEP, pi = rates(arm, kd, W)
    y0 = np.zeros((M, 6))
    y0[:, 0], y0[:, 1] = S0.ravel(), E0
    idx, ptr = np.arange(M), np.arange(M + 1)

    def fun(_, z):
        return rhs_jac(z.reshape(M, 6), lm, kE, kES, kEP, pi)[0].ravel()

    def jac(_, z):
        return bsr_matrix((rhs_jac(z.reshape(M, 6), lm, kE, kES, kEP, pi)[1], idx, ptr), shape=(6 * M, 6 * M))
    sol = solve_ivp(fun, (T_GRID[0], T_GRID[-1]), y0.ravel(), method="BDF", t_eval=T_GRID, jac=jac,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    Y = sol.y.reshape(M, 6, T_GRID.size)
    inv1 = np.max(np.abs(Y[:, 0] + Y[:, 2] + Y[:, 3] + Y[:, 4] - S0.ravel()[:, None]))
    inv2 = np.max(np.abs(Y[:, 1] + Y[:, 2] + Y[:, 4] + Y[:, 5] - E0))
    return (Y[:, 0] + Y[:, 2]).reshape(B, W, T_GRID.size), float(max(inv1, inv2))


def gen(arm, tag, n=NB, rtol=RTOL, atol=ATOL):
    W = S0_NOM.size
    rng = np.random.default_rng(seed(tag))
    vm = rng.uniform(0.18, 0.45, size=n)
    kd = rng.uniform(0.0, 0.025, size=n)
    S0 = S0_NOM[None, :] * np.clip(rng.normal(1.0, 0.05, size=(n, W)), 0.7, 1.3)
    tot, F = integrate(arm, S0, vm / (K_CAT * E0), kd, rtol, atol)
    clean = ife(tot, A_IFE)
    Y = clean + rng.normal(0.0, 1.0, size=clean.shape) * (0.02 * clean + ARMS[arm]["floor"])
    return Y, clean, F


def block(spec):
    arm, d = spec
    t0 = time.perf_counter()
    Y, _, F = gen(arm, f"q64_{arm}_null_d{d}")
    nom = np.broadcast_to(S0_NOM[None, :], Y.shape[:2]).copy()
    r = k57_map_fit(Y, nom, T_GRID, SIG)
    return d, int(np.sum(r["called"])), int(Y.shape[0]), F, time.perf_counter() - t0


def jac_check(arm):
    rng = np.random.default_rng(seed(f"q64_{arm}_jaccheck"))
    M = 8
    y = np.abs(rng.normal(0.0, 1.0, size=(M, 6))) * np.array([5.0, 0.01, 0.01, 2.0, 0.005, 0.005])
    lam, kd = rng.uniform(9.0, 23.0, M), rng.uniform(0.0, 0.025, M)
    kE, kES, kEP, pi = rates(arm, kd, 1)
    _, J = rhs_jac(y, lam, kE, kES, kEP, pi)
    worst = 0.0
    for j in range(6):
        h = 1e-7 * np.maximum(1.0, np.abs(y[:, j]))
        yp, ym = y.copy(), y.copy()
        yp[:, j] += h
        ym[:, j] -= h
        fd = (rhs_jac(yp, lam, kE, kES, kEP, pi)[0] - rhs_jac(ym, lam, kE, kES, kEP, pi)[0]) / (2 * h[:, None])
        worst = max(worst, float(np.max(np.abs(fd - J[:, :, j])) / max(float(np.max(np.abs(J))), 1e-30)))
    return worst


def write(name, L, out):
    os.makedirs(LOGS, exist_ok=True)
    with open(os.path.join(LOGS, f"{name}.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(out, fh, indent=1)
    with open(os.path.join(LOGS, f"{name}.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")
    print("\n".join(L))


def run_arm(tag, arm):
    from concurrent.futures import ProcessPoolExecutor
    t0 = time.perf_counter()
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE PHASE 4 (Q64): frozen K57 MAP rule vs confounder arm {arm} | {tag}")
    A("=" * 78)
    A(f"  arm {arm}: {ARMS[arm]['name']} (MODELED; rates ASSUMED)")
    A(f"  Km = (k_off + k_cat)/k_on = {(K_OFF + K_CAT) / K_ON:.3f}; Kp = kp_off/kp_on = {KP_OFF / KP_ON:.3f}; E0 {E0}; "
      f"A_peak {S0_NOM[-1] * A_IFE:.3f} (ASSUMED)")
    A("  ladder A1.000 (ASSUMED): " + ", ".join(f"{x:.4f}" for x in S0_NOM) + f"; reads {T_GRID.size} on 0..10")
    sdg = 0.02 * ife(S0_NOM, A_IFE) + ARMS[arm]["floor"]
    A("  generator sd / rule SIG at t 0 per well: " + ", ".join(f"{g / s:.2f}" for g, s in zip(sdg, SIG))
      + " (MODELED vs ASSUMED)")
    jc = jac_check(arm)
    A(f"  G1 generator Jacobian vs central FD: max rel err {jc:.2e} (limit 1e-6, MEASURED numerics)")
    _, ca, Fa = gen(arm, f"q64_{arm}_tolcheck", n=20)
    _, cb, Fb = gen(arm, f"q64_{arm}_tolcheck", n=20, rtol=RTOL / 2, atol=ATOL / 2)
    dY = float(np.max(np.abs(ca - cb)))
    A(f"  G1 tolerance halving (20 plates): max |d readout| {dY:.2e} (limit 1e-6, MEASURED numerics)")
    with ProcessPoolExecutor(max_workers=NBLK) as ex:
        res = sorted(ex.map(block, [(arm, d) for d in range(NBLK)]))
    Fmax = max([r[3] for r in res] + [Fa, Fb])
    A(f"  G1 invariants on every plate: max |F| {Fmax:.2e} (limit 1e-8, MEASURED numerics)")
    for d, k, n, _, _ in res:
        A(f"  block d{d}: false calls {k}/{n} (MODELED)")
    k = sum(r[1] for r in res)
    n = sum(r[2] for r in res)
    w = wilson_up(k, n)
    A(f"  arm {arm} false calls {k}/{n} = {k / n:.4f}, Wilson 90% upper {w:.4f} (limit 0.060) (MODELED)")
    g1 = jc <= 1e-6 and dY <= 1e-6 and Fmax <= 1e-8
    dt = time.perf_counter() - t0
    A(f"  wall time {dt:.2f} s (MEASURED; {NBLK} worker processes)")
    out = {"arm": arm, "tag": tag, "jac": jc, "dY": dY, "Fmax": Fmax, "k": k, "n": n, "w90": w, "g1": bool(g1),
           "wall": dt}
    write(f"phase4_q64_{arm}_{tag}", L, out)


def run_pool(tag):
    t0 = time.perf_counter()
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE PHASE 4 (Q64): pooled gates over arms A, B, C | {tag}")
    A("=" * 78)
    A(f"  rule frozen: k57_map_fit thr {K57_THR}, iksi {K57_IKSI} (ASSUMED). Baseline Cycle 66 MM + decay: 10/450, "
      "w90 0.0368 (MODELED)")
    o = {}
    for arm in ARMS:
        with open(os.path.join(LOGS, f"phase4_q64_{arm}_{tag}.json"), encoding="utf-8") as fh:
            o[arm] = json.load(fh)
    g1 = all(o[a]["g1"] for a in o)
    A(f"  G1 worst over arms: Jacobian FD {max(o[a]['jac'] for a in o):.2e} (<= 1e-6), |F| "
      f"{max(o[a]['Fmax'] for a in o):.2e} (<= 1e-8), tolerance halving {max(o[a]['dY'] for a in o):.2e} (<= 1e-6) "
      f"... {'PASS' if g1 else 'FAIL'}")
    fails = []
    for a in o:
        r = o[a]
        ok = r["w90"] <= 0.060
        A(f"  arm {a} ({ARMS[a]['name']}): {r['k']}/{r['n']} = {r['k'] / r['n']:.4f}, w90 {r['w90']:.4f} ... "
          f"{'PASS' if ok else 'FAIL'} (MODELED)")
        if not ok:
            fails.append(a)
    K = sum(o[a]["k"] for a in o)
    N = sum(o[a]["n"] for a in o)
    wp = wilson_up(K, N)
    A(f"  pooled: {K}/{N} = {K / N:.4f}, w90 {wp:.4f} ... {'PASS' if wp <= 0.060 else 'FAIL'} (MODELED)")
    g2 = not fails and wp <= 0.060
    g3 = all(o[a]["wall"] < 30.0 for a in o)
    A(f"  G3 arm invocations < 30 s: {g3} (MEASURED; per-arm wall times are in the arm logs)")
    gates = [("G1", g1), ("G2", g2), ("G3", g3)]
    A("  gates: " + ", ".join(f"{g} {'PASS' if v else 'FAIL'}" for g, v in gates))
    if g1 and g2 and g3:
        A("  VERDICT: Q64 PASS -- no confounder arm pushes the frozen rule's false-alarm bound above 0.060 (MODELED)")
    else:
        for a in fails:
            A(f"  L21: mechanism '{ARMS[a]['name']}' moved the false-alarm rate from 10/450 (Cycle 66) to "
              f"{o[a]['k']}/{o[a]['n']} (w90 {o[a]['w90']:.4f}); no threshold re-tuned")
        A("  VERDICT: Q64 FAIL -- " + ", ".join(g for g, v in gates if not v))
    A(f"  wall time {time.perf_counter() - t0:.2f} s (MEASURED)")
    write(f"phase4_q64_pool_{tag}", L, {"K": K, "N": N, "w90": wp, "gates": dict(gates)})


if __name__ == "__main__":
    if sys.argv[2] == "pool":
        run_pool(sys.argv[1])
    else:
        run_arm(sys.argv[1], sys.argv[2])
