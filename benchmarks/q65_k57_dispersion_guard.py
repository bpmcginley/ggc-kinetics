"""
GGC-inf Release Phase 5 (v0.2.0-dev), item 1 (Q65): zero-re-tuning noise-scale fix for the K57 MAP rule.
PRE-REGISTERED: this docstring was written and hashed before the first run (L8). Simulation only (MODELED); no
downloads, no wet lab. Every number carries an L7 tag.

QUESTION (L21): Q64 found the frozen rule's false-alarm bound rises to w90 0.0987 (MODELED) when the readout noise
  floor is 2x the rule's assumed SIG. The movable number is the per-arm false-alarm bound under noise
  mis-specification, with the decision rule itself unchanged. Two ways to set the noise scale are compared with the
  frozen baseline. FROZEN (ASSUMED, never changed): K57_THR -2.5, K57_IKSI 0.045, MAP prior sp 0.05, n_iter 24.
  src/ggc_kinetics/solvers.py is NOT modified; its functions are only called.
VARIANTS (all call the frozen ggc_kinetics.solvers.k57_map_fit; only the sig argument differs):
  Row 0  sig = SIG = 0.015 x max(S0_NOM, 0.5) (the rule's own ASSUMED SIG), one batched call per plate set.
  V1     blank-estimated sig, per plate: 3 independent noise-only blank reads per S0 level (27 per plate, drawn as
         N(0, floor) from a separate seed stream q65_{arm}_blank_d{d}; blank signal is 0, so blanks carry the floor
         only). floor_hat = sqrt(sum x^2 / 27) (known zero mean: blanks are offset-subtracted).
         sig_V1[w] = 0.02 x max(Y[w, read 0], 0) + floor_hat (2% proportional term on the well's own first read).
  V2     one-sided Birge dispersion guard, per plate: null fit = the FIRST null stage of k57_map_fit, reproduced with
         the same initialisation (IFE grid, initial-rate start) and the frozen solvers._k57_lm (free = [1,1,1,1,0],
         n_iter 24) at sig = SIG. chi2_null = LM cost - sum (U / sp)^2 (data part only).
         dof_null = n_obs - n_params_null = W x M - (4 kinetic + W per-well u) = 117 - 13 = 104 (u counted as free:
         conservative). birge = max(1, sqrt(chi2_null / dof_null)). sig_eff = SIG x birge, and the plate is
         re-fitted with k57_map_fit(..., sig_eff) (null AND SI fits re-evaluated at sig_eff, prior sp 0.05 unchanged).
         Plates with birge = 1 keep the Row 0 result (identical inputs; batch vs single-plate calls were checked
         equal, max |d dBIC| 0.0, MEASURED). No new solver or Jacobian is introduced by V2.
GENERATOR (independent: zero code shared with ggc_kinetics.solvers; own RHS, Jacobian, IFE, noise, Wilson bound).
  7 species S, E, ES, ESS, P, EP, Ed: the union of the Q64 (EP) and Cycle 66 (ESS) generators (disclosed; the prompt
  says "6-species"). Safe textbook mechanisms only. All rates ASSUMED.
  E + S <-> ES -> E + P: k_on 100, k_off 99, k_cat 1 (Km 1.000), E0 0.02. ES + S <-> ESS: ksi_on 20, ksi_off 100
  (Ksi 5.000; positive arms only). E + P <-> EP: kp_on 100, kp_off 300 (Kp 3.000; arm A only). Decay: E at kd,
  ES at fes x kd, ESS at kd (releases 2 S), EP at kd (releases P). Bound substrate is released on decay.
  Per plate lam = vm/(k_cat E0), vm ~ U(0.18, 0.45); kd ~ U(0, 0.025). S0 = S0_NOM x clip(N(1, 0.05), 0.7, 1.3).
  Readout = S + ES + 2 ESS through IFE a = 1/16 plus N(0, 0.02 x signal + floor). Ladder A1.000 = geomspace(0.2, 16, 9),
  13 reads on 0..10 (ASSUMED). BDF rtol 1e-8, atol 1e-12, hand-written bsr Jacobian.
ARMS (3 seed blocks x 150 = 450 plates each; seeds zlib.crc32 of q65_{arm}_null_d{d} / q65_{arm}_pos_d{d}, fresh):
  N0 MM + decay, floor 0.0075 | C2 floor 0.015 | C3 floor 0.0225 | A product inhibition Kp = 3 Km, floor 0.0075 |
  B ES inactivated at 2 kd, floor 0.0075 | P1 substrate inhibition Ksi 5, floor 0.0075 | P2 Ksi 5, floor 0.015.
SPLIT: one invocation per arm; each block is split in 4 slices over 12 worker processes (each worker regenerates its
  block deterministically and fits only its slice). Then one "pool" invocation applies the gates. Each < 30 s.
GATES (ASSUMED rules), applied to each variant separately:
  G1 numerics per arm: generator Jacobian vs central FD max rel err <= 1e-6; invariants S + ES + 2 ESS + P + EP = S0
     and E + ES + ESS + EP + Ed = E0 to max |F| <= 1e-8 on every plate; tolerance halving (20 plates) max |d readout|
     <= 1e-6 (MEASURED numerics).
  G2 specificity: the variant passes iff EVERY null arm N0, C2, C3, A, B has Wilson 90% one-sided upper bound
     (z 1.6448536) w90 <= 0.060 (MODELED).
  G3 power: P1 recall >= 0.90 AND P1 recall >= Row 0 P1 recall - 0.05; P2 recall reported (MODELED).
  G4 every invocation < 30 s (MEASURED).
  Q65 verdict: PASS for a variant iff G1-G4 all pass for it; report Row 0 against the same gates as a reference.
MY PREDICTIONS (before run1): Row 0 fails G2 on C2 (about 30-40/450, as Q64) and C3 (worse, w90 > 0.15). V1 matches the
  generator's noise model closely, so it passes G2 on every arm; P1 recall drops a little (0.88-0.95), so G3 is
  borderline. V2 inflates sig uniformly; the generator's sd/SIG ratio is not uniform across wells (about 2 at low
  wells, 0.55 at the top well), so V2 should pass C2 and may be borderline on C3; P1 recall 0.85-0.95.
Run:  python -I -X utf8 benchmarks/q65_k57_dispersion_guard.py TAG ARM   (ARM in N0 C2 C3 A B P1 P2, then pool)
Writes ../ggc_inf_kit/logs/phase5_q65_{ARM|pool}_TAG.{txt,json}. Run twice (run1, run2) and diff; only tag and
wall-time lines may differ.
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
from ggc_kinetics import solvers  # noqa: E402  (frozen rule under test; only called)

K_ON, K_OFF, K_CAT, E0 = 100.0, 99.0, 1.0, 0.02
KSI_ON, KSI_OFF = 20.0, 100.0
KP_ON, KP_OFF = 100.0, 300.0
S0_NOM = np.geomspace(0.2, 16.0, 9)
T_GRID = np.linspace(0.0, 10.0, 13)
SIG = 0.015 * np.maximum(S0_NOM, 0.5)
A_IFE = 1.0 / 16.0
NB, NBLK, NSL = 150, 3, 4
RTOL, ATOL = 1e-8, 1e-12
ZW = 1.6448536
N_BLANK = 3
DOF_NULL = S0_NOM.size * T_GRID.size - (4 + S0_NOM.size)
ARMS = {"N0": dict(pi=0, si=0, fes=1.0, floor=0.0075, name="MM + decay, floor 0.0075"),
        "C2": dict(pi=0, si=0, fes=1.0, floor=0.015, name="MM + decay, floor 0.015 (2x)"),
        "C3": dict(pi=0, si=0, fes=1.0, floor=0.0225, name="MM + decay, floor 0.0225 (3x)"),
        "A": dict(pi=1, si=0, fes=1.0, floor=0.0075, name="product inhibition Kp = 3 Km"),
        "B": dict(pi=0, si=0, fes=2.0, floor=0.0075, name="ES inactivated at 2 kd"),
        "P1": dict(pi=0, si=1, fes=1.0, floor=0.0075, name="substrate inhibition Ksi 5, floor 0.0075"),
        "P2": dict(pi=0, si=1, fes=1.0, floor=0.015, name="substrate inhibition Ksi 5, floor 0.015")}
NULLS, POS = ("N0", "C2", "C3", "A", "B"), ("P1", "P2")
VARS = ("row0", "v1", "v2")


def seed(tag):
    return zlib.crc32(tag.encode("utf-8")) & 0xFFFFFFFF


def wilson_up(k, n, z=ZW):
    p = k / n
    den = 1.0 + z * z / n
    return (p + z * z / (2 * n) + z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / den


def ife(S, a):
    k = a * math.log(10.0)
    return -np.expm1(-k * S) / k


def rhs_jac(y, lam, kE, kES, kESS, kEP, pi, si):
    """y (M, 7) = S, E, ES, ESS, P, EP, Ed. Returns f (M, 7), J (M, 7, 7)."""
    S, E, ES, ESS, P, EP = (y[:, i] for i in range(6))
    a, b, c = lam * K_ON, lam * K_OFF, lam * K_CAT
    d, e = lam * KSI_ON * si, lam * KSI_OFF * si
    pn, pf = lam * KP_ON * pi, lam * KP_OFF * pi
    f = np.empty_like(y)
    f[:, 0] = -a * E * S + b * ES + kES * ES - d * ES * S + e * ESS + 2 * kESS * ESS
    f[:, 1] = -a * E * S + (b + c) * ES - kE * E - pn * E * P + pf * EP
    f[:, 2] = a * E * S - (b + c + kES) * ES - d * ES * S + e * ESS
    f[:, 3] = d * ES * S - (e + kESS) * ESS
    f[:, 4] = c * ES - pn * E * P + pf * EP + kEP * EP
    f[:, 5] = pn * E * P - (pf + kEP) * EP
    f[:, 6] = kE * E + kES * ES + kESS * ESS + kEP * EP
    J = np.zeros((y.shape[0], 7, 7))
    J[:, 0, 0], J[:, 0, 1], J[:, 0, 2], J[:, 0, 3] = -a * E - d * ES, -a * S, b + kES - d * S, e + 2 * kESS
    J[:, 1, 0], J[:, 1, 1], J[:, 1, 2], J[:, 1, 4], J[:, 1, 5] = -a * E, -a * S - kE - pn * P, b + c, -pn * E, pf
    J[:, 2, 0], J[:, 2, 1], J[:, 2, 2], J[:, 2, 3] = a * E - d * ES, a * S, -(b + c + kES) - d * S, e
    J[:, 3, 0], J[:, 3, 2], J[:, 3, 3] = d * ES, d * S, -(e + kESS)
    J[:, 4, 1], J[:, 4, 2], J[:, 4, 4], J[:, 4, 5] = -pn * P, c, -pn * E, pf + kEP
    J[:, 5, 1], J[:, 5, 4], J[:, 5, 5] = pn * P, pn * E, -(pf + kEP)
    J[:, 6, 1], J[:, 6, 2], J[:, 6, 3], J[:, 6, 5] = kE, kES, kESS, kEP
    return f, J


def rates(arm, kd, W):
    p = ARMS[arm]
    kdm = np.repeat(kd, W)
    n = kdm.size
    return kdm, p["fes"] * kdm, kdm, kdm * p["pi"], np.full(n, float(p["pi"])), np.full(n, float(p["si"]))


def integrate(arm, S0, lam, kd, rtol=RTOL, atol=ATOL):
    from scipy.integrate import solve_ivp
    from scipy.sparse import bsr_matrix
    B, W = S0.shape
    M = B * W
    lm = np.repeat(lam, W)
    rp = rates(arm, kd, W)
    y0 = np.zeros((M, 7))
    y0[:, 0], y0[:, 1] = S0.ravel(), E0
    idx, ptr = np.arange(M), np.arange(M + 1)

    def fun(_, z):
        return rhs_jac(z.reshape(M, 7), lm, *rp)[0].ravel()

    def jac(_, z):
        return bsr_matrix((rhs_jac(z.reshape(M, 7), lm, *rp)[1], idx, ptr), shape=(7 * M, 7 * M))
    sol = solve_ivp(fun, (T_GRID[0], T_GRID[-1]), y0.ravel(), method="BDF", t_eval=T_GRID, jac=jac,
                    rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    Y = sol.y.reshape(M, 7, T_GRID.size)
    inv1 = np.max(np.abs(Y[:, 0] + Y[:, 2] + 2 * Y[:, 3] + Y[:, 4] + Y[:, 5] - S0.ravel()[:, None]))
    inv2 = np.max(np.abs(Y[:, 1] + Y[:, 2] + Y[:, 3] + Y[:, 5] + Y[:, 6] - E0))
    return (Y[:, 0] + Y[:, 2] + 2 * Y[:, 3]).reshape(B, W, T_GRID.size), float(max(inv1, inv2))


def block_tag(arm, d):
    return f"q65_{arm}_{'pos' if arm in POS else 'null'}_d{d}"


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


def birge_null(Y, s0n):
    """First null stage of k57_map_fit (same initialisation), frozen solvers._k57_lm at sig = SIG."""
    B = Y.shape[0]
    ag = np.linspace(0.0, 0.10, 61)
    chi = np.sum(((Y[:, :, 0][:, :, None] - solvers.ife_forward(s0n[:, :, None], ag[None, None, :]))
                  / SIG[None, :, None]) ** 2, axis=1)
    a0 = ag[np.argmin(chi, axis=1)]
    Sr = solvers.ife_inverse(Y, a0[:, None, None])
    v = np.maximum((Sr[:, :, 0] - Sr[:, :, 2]) / (T_GRID[2] - T_GRID[0]), 1e-4)
    P0 = np.stack([np.log(np.clip(np.max(v, axis=1) * 1.25, 1e-3, 10.0)), np.zeros(B),
                   np.full(B, 0.015), a0, np.zeros(B)], axis=1)
    _, c0, U0 = solvers._k57_lm(Y, s0n, P0, [1, 1, 1, 1, 0], T_GRID, SIG, solvers.K57_SP, 24)
    chi2 = c0 - np.sum((U0 / solvers.K57_SP) ** 2, axis=1)
    return np.maximum(1.0, np.sqrt(chi2 / DOF_NULL))


def work(spec):
    arm, d, s = spec
    Y, _, F = gen(arm, block_tag(arm, d))
    rngb = np.random.default_rng(seed(f"q65_{arm}_blank_d{d}"))
    blanks = rngb.normal(0.0, ARMS[arm]["floor"], size=(NB, S0_NOM.size * N_BLANK))
    sl = np.array_split(np.arange(NB), NSL)[s]
    Y, blanks = Y[sl], blanks[sl]
    nom = np.broadcast_to(S0_NOM[None, :], Y.shape[:2]).copy()
    r0 = solvers.k57_map_fit(Y, nom, T_GRID, SIG)
    c0 = np.asarray(r0["called"], bool)
    fh = np.sqrt(np.mean(blanks ** 2, axis=1))
    bg = birge_null(Y, nom)
    c1, c2 = np.zeros(len(sl), bool), c0.copy()
    for i in range(len(sl)):
        s1 = 0.02 * np.maximum(Y[i, :, 0], 0.0) + fh[i]
        c1[i] = bool(solvers.k57_map_fit(Y[i:i + 1], nom[i:i + 1], T_GRID, s1)["called"][0])
        if bg[i] > 1.0:
            c2[i] = bool(solvers.k57_map_fit(Y[i:i + 1], nom[i:i + 1], T_GRID, SIG * bg[i])["called"][0])
    return d, s, F, c0.tolist(), c1.tolist(), c2.tolist(), fh.tolist(), bg.tolist()


def jac_check(arm):
    rng = np.random.default_rng(seed(f"q65_{arm}_jaccheck"))
    M = 8
    y = np.abs(rng.normal(0.0, 1.0, size=(M, 7))) * np.array([5.0, 0.01, 0.01, 0.005, 2.0, 0.005, 0.005])
    lam, kd = rng.uniform(9.0, 23.0, M), rng.uniform(0.0, 0.025, M)
    rp = rates(arm, kd, 1)
    _, J = rhs_jac(y, lam, *rp)
    worst = 0.0
    for j in range(7):
        h = 1e-7 * np.maximum(1.0, np.abs(y[:, j]))
        yp, ym = y.copy(), y.copy()
        yp[:, j] += h
        ym[:, j] -= h
        fd = (rhs_jac(yp, lam, *rp)[0] - rhs_jac(ym, lam, *rp)[0]) / (2 * h[:, None])
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
    A(f"GGC-inf RELEASE PHASE 5 (Q65): K57 noise-scale variants, arm {arm} | {tag}")
    A("=" * 78)
    A(f"  arm {arm}: {ARMS[arm]['name']} (MODELED; rates ASSUMED); {'positive' if arm in POS else 'null'} arm")
    A(f"  Km {(K_OFF + K_CAT) / K_ON:.3f}, Ksi {KSI_OFF / KSI_ON:.3f}, Kp {KP_OFF / KP_ON:.3f}, E0 {E0}; frozen rule thr "
      f"{solvers.K57_THR}, iksi {solvers.K57_IKSI}, sp {solvers.K57_SP} (ASSUMED); dof_null {DOF_NULL}")
    sdg = 0.02 * ife(S0_NOM, A_IFE) + ARMS[arm]["floor"]
    A("  generator sd / SIG at t 0 per well: " + ", ".join(f"{g / s:.2f}" for g, s in zip(sdg, SIG))
      + " (MODELED vs ASSUMED)")
    jc = jac_check(arm)
    A(f"  G1 generator Jacobian vs central FD: max rel err {jc:.2e} (limit 1e-6, MEASURED numerics)")
    _, ca, Fa = gen(arm, f"q65_{arm}_tolcheck", n=20)
    _, cb, Fb = gen(arm, f"q65_{arm}_tolcheck", n=20, rtol=RTOL / 2, atol=ATOL / 2)
    dY = float(np.max(np.abs(ca - cb)))
    A(f"  G1 tolerance halving (20 plates): max |d readout| {dY:.2e} (limit 1e-6, MEASURED numerics)")
    specs = [(arm, d, s) for d in range(NBLK) for s in range(NSL)]
    with ProcessPoolExecutor(max_workers=len(specs)) as ex:
        res = sorted(ex.map(work, specs), key=lambda r: (r[0], r[1]))
    Fmax = max([r[2] for r in res] + [Fa, Fb])
    A(f"  G1 invariants on every plate: max |F| {Fmax:.2e} (limit 1e-8, MEASURED numerics)")
    out = {"arm": arm, "tag": tag, "jac": jc, "dY": dY, "Fmax": Fmax,
           "g1": bool(jc <= 1e-6 and dY <= 1e-6 and Fmax <= 1e-8)}
    fh = np.concatenate([r[6] for r in res])
    bg = np.concatenate([r[7] for r in res])
    A(f"  V1 floor_hat median {np.median(fh):.5f} (true floor {ARMS[arm]['floor']}); V2 birge median {np.median(bg):.3f}, "
      f"range {bg.min():.3f}-{bg.max():.3f}, birge > 1 on {int(np.sum(bg > 1))}/{bg.size} plates "
      "(MODELED fit of MODELED data)")
    for j, v in enumerate(VARS):
        per = {}
        for d in range(NBLK):
            per[d] = sum(sum(r[3 + j]) for r in res if r[0] == d)
        k = sum(per.values())
        n = NB * NBLK
        w = wilson_up(k, n)
        out[v] = {"k": k, "n": n, "rate": k / n, "w90": w}
        lab = "recall" if arm in POS else "false calls"
        A(f"  {v}: {lab} {k}/{n} = {k / n:.4f} (blocks " + ", ".join(f"d{d} {per[d]}" for d in range(NBLK))
          + (f"), w90 {w:.4f}" if arm in NULLS else ")") + " (MODELED)")
    dt = time.perf_counter() - t0
    out["wall"] = dt
    A(f"  wall time {dt:.2f} s (MEASURED; {len(specs)} worker processes)")
    write(f"phase5_q65_{arm}_{tag}", L, out)


def run_pool(tag):
    t0 = time.perf_counter()
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE PHASE 5 (Q65): gates for Row 0, V1 (blank sig), V2 (Birge guard) | {tag}")
    A("=" * 78)
    o = {}
    for arm in ARMS:
        with open(os.path.join(LOGS, f"phase5_q65_{arm}_{tag}.json"), encoding="utf-8") as fh:
            o[arm] = json.load(fh)
    g1 = all(o[a]["g1"] for a in o)
    g4 = all(o[a]["wall"] < 30.0 for a in o)
    A(f"  G1 worst over arms: Jacobian FD {max(o[a]['jac'] for a in o):.2e} (<= 1e-6), |F| "
      f"{max(o[a]['Fmax'] for a in o):.2e} (<= 1e-8), tolerance halving {max(o[a]['dY'] for a in o):.2e} (<= 1e-6) "
      f"... {'PASS' if g1 else 'FAIL'} (MEASURED numerics)")
    A(f"  G4 every arm invocation < 30 s: {g4} (MEASURED; wall times in the arm logs)")
    rec0 = o["P1"]["row0"]["rate"]
    verdicts = {}
    for v in VARS:
        A(f"  -- {v} --")
        bad = []
        for a in NULLS:
            r = o[a][v]
            ok = r["w90"] <= 0.060
            if not ok:
                bad.append(a)
            A(f"    {a}: {r['k']}/{r['n']}, w90 {r['w90']:.4f} ... {'PASS' if ok else 'FAIL'} (MODELED)")
        g2 = not bad
        p1, p2 = o["P1"][v]["rate"], o["P2"][v]["rate"]
        g3 = p1 >= 0.90 and p1 >= rec0 - 0.05
        A(f"    G2 specificity ... {'PASS' if g2 else 'FAIL ' + ', '.join(bad)}")
        A(f"    G3 P1 recall {p1:.4f} (need >= 0.90 and >= Row 0 {rec0:.4f} - 0.05) ... {'PASS' if g3 else 'FAIL'}; "
          f"P2 recall {p2:.4f} (reported) (MODELED)")
        ok = g1 and g2 and g3 and g4
        verdicts[v] = {"G1": g1, "G2": g2, "G3": g3, "G4": g4, "pass": ok, "g2_fail_arms": bad}
        A(f"    {v}: G1 {'PASS' if g1 else 'FAIL'}, G2 {'PASS' if g2 else 'FAIL'}, G3 {'PASS' if g3 else 'FAIL'}, "
          f"G4 {'PASS' if g4 else 'FAIL'} -> {'PASS' if ok else 'FAIL'}")
    win = [v for v in ("v1", "v2") if verdicts[v]["pass"]]
    if win:
        A(f"  VERDICT: Q65 PASS for {', '.join(win)} (MODELED); frozen thresholds unchanged; solvers.py untouched "
          "pending audit")
    else:
        A("  VERDICT: Q65 FAIL -- neither V1 nor V2 passes all gates (MODELED); no threshold re-tuned")
    A(f"  wall time {time.perf_counter() - t0:.2f} s (MEASURED)")
    write(f"phase5_q65_pool_{tag}", L, {"verdicts": verdicts, "winners": win})


if __name__ == "__main__":
    if sys.argv[2] == "pool":
        run_pool(sys.argv[1])
    else:
        run_arm(sys.argv[1], sys.argv[2])
