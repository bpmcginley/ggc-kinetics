"""ggc-kinetics v0.2.0 regression tests (target: < 5 s). Run: PYTHONPATH=src python -m unittest -v tests.test_ggc_kinetics

Evidence tags (L7):
  invertase (1913, PMC3381512 Table 1, 31 points)  MEASURED data, fit pinned to the archive suite R1
  PON1 (PMC6896510 Table 1)                         MEASURED published errors, pinned to the suite R2
  RK4 parity, Morrison, K57 FD and |F|               MODELED inputs, MEASURED numerics
  parser fixtures in tests/data                      MODELED plates (Cycle 33), stored reports from Cycle 34
"""
import hashlib
import io
import json
import math
import os
import tempfile
import unittest
import zlib
from contextlib import redirect_stdout

import numpy as np

import ggc_kinetics as gk
from ggc_kinetics import cli, parser

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
PP_SHA256 = "7c4485128b7a3c6348035398ae531f4869625113f40ebb39712537d59c791c48"

MM1913 = [
    (333.0, [7.0, 14.0, 26.0, 49.0, 75.0, 117.0, 1052.0], [0.0164, 0.0316, 0.0528, 0.0923, 0.1404, 0.2137, 0.9834]),
    (166.7, [8.0, 16.0, 28.0, 52.0, 82.0, 103.0], [0.0350, 0.0636, 0.1080, 0.1980, 0.3000, 0.3780]),
    (83.0, [49.5, 90.0, 125.0, 151.0, 208.0], [0.352, 0.575, 0.690, 0.766, 0.900]),
    (41.6, [10.25, 30.75, 61.75, 90.75, 112.70, 132.70, 154.70], [0.1147, 0.3722, 0.615, 0.747, 0.850, 0.925, 0.940]),
    (20.8, [17.0, 27.0, 38.0, 62.0, 95.0, 1372.0], [0.331, 0.452, 0.611, 0.736, 0.860, 0.990]),
]

# K57 test design (Cycle 37 grid): 11 nominal S0 levels, 13 reads over 10 min, readout SD 1.5% of max(S0, 0.5).
S0_NOM = np.sort(np.concatenate([np.geomspace(0.125, 16.0, 9), [5.5, 11.0]]))
T_GRID = np.linspace(0.0, 10.0, 13)
SIG = 0.015 * np.maximum(S0_NOM, 0.5)


def rng_for(tag):
    return np.random.default_rng(zlib.crc32(tag.encode("utf-8")) & 0xFFFFFFFF)


def gen_k57(tag, B, mode, cv=0.05):
    """Copy of the Cycle 35/37 generator (MODELED): mode 0 null, 1 Ksi = 5 Km."""
    rng = rng_for(tag)
    Km = 10.0 ** rng.uniform(-0.2, 0.35, size=B)
    Vm = rng.uniform(0.18, 0.45, size=B)
    kd = rng.uniform(0.0, 0.025, size=B)
    S0 = S0_NOM[None, :] * np.clip(rng.normal(1.0, cv, size=(B, 11)), 0.7, 1.3)
    ik = np.zeros(B) if mode == 0 else 1.0 / (5.0 * Km)
    S, _ = gk.ksi_solve_conv(S0, T_GRID, np.log(Vm), np.log(Km), kd, ik)
    Y = gk.solvers.ife_forward(S, 1.0 / 16.0) + rng.normal(0.0, SIG[None, :, None], size=S.shape)
    return Y, S0, Km, Vm, kd, ik


class TestSolvers(unittest.TestCase):
    def test_invertase_1913(self):
        S0 = np.concatenate([np.full(len(t), s) for s, t, _ in MM1913])
        t = np.concatenate([t for _, t, _ in MM1913])
        f = np.concatenate([f for _, _, f in MM1913])
        r = gk.fit_kp_vs_kd(S0, t, f)
        self.assertAlmostEqual(r["Km"], 16.52, delta=0.005)
        self.assertAlmostEqual(100 * (r["Km"] - 16.5) / 16.5, 0.11, delta=0.005)
        self.assertAlmostEqual(r["eff"], 0.0461, delta=0.00005)
        self.assertAlmostEqual(100 * (r["eff"] - 0.0460) / 0.0460, 0.32, delta=0.005)
        self.assertAlmostEqual(r["dbic"], -19.29, delta=0.005)

    def test_pon1_published_errors(self):
        gp = np.array([-0.56, -0.58, -0.46, -0.69, -0.53, -0.57])     # GraphPad LW vs ODE, %
        df = np.array([-11.2, -15.3, 16.5, 5.1, -17.2, 1.7])           # DynaFit LW vs ODE, %
        self.assertLessEqual(float(np.max(np.abs(gp))), 0.69)
        self.assertGreater(float(np.max(np.abs(df))), 15.0)

    def test_coupled_lag_rk4_parity(self):
        v0, kd, k2 = 12.0, 0.04, 0.35
        tg = np.linspace(0.0, 20.0, 41)
        dt, I, P, out = 0.002, 0.0, 0.0, [0.0]
        for step in range(1, 10001):
            t0 = (step - 1) * dt

            def f(tt, ii):
                return v0 * math.exp(-kd * tt) - k2 * ii, k2 * ii
            a1, b1 = f(t0, I)
            a2, b2 = f(t0 + dt / 2, I + dt / 2 * a1)
            a3, b3 = f(t0 + dt / 2, I + dt / 2 * a2)
            a4, b4 = f(t0 + dt, I + dt * a3)
            I += dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            P += dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
            if step % 250 == 0:
                out.append(P)
        self.assertLessEqual(float(np.max(np.abs(gk.coupled_lag_decay(tg, v0, kd, k2) - np.array(out)))), 1e-11)

    def test_slow_binding_rk4_parity(self):
        S0, Km, Vm, ko, ess = 200.0, 100.0, 18.0, 0.25, 0.2
        dt, S, E, out = 0.002, S0, 1.0, [S0]
        for step in range(1, 10001):
            def f(s, e):
                return -Vm * e * s / (Km + s), -ko * (e - ess)
            a1, b1 = f(S, E)
            a2, b2 = f(S + dt / 2 * a1, E + dt / 2 * b1)
            a3, b3 = f(S + dt / 2 * a2, E + dt / 2 * b2)
            a4, b4 = f(S + dt * a3, E + dt * b3)
            S += dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            E += dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
            if step % 250 == 0:
                out.append(S)
        tg = np.linspace(0.0, 20.0, 41)
        self.assertLessEqual(float(np.max(np.abs(gk.slow_binding_warped(tg, S0, Km, Vm, ko, ess) - np.array(out)))),
                             1e-8)

    def test_morrison_residual(self):
        S0 = np.array([[0.25, 1.0, 4.0, 16.0], [0.5, 2.0, 8.0, 0.1]])
        S, res = gk.morrison_exact(S0, np.linspace(0, 10, 13), np.array([0.3, 1.2]), np.array([0.05, 0.5]),
                                   np.array([1.0, 0.3]))
        self.assertLessEqual(res, 1e-10)
        self.assertTrue(np.all(S <= S0[:, :, None] + 1e-12))

    def test_strang_exact_without_hydrolysis(self):
        S0 = np.array([0.5, 2.0, 6.0])
        tg = np.linspace(0, 30, 31)
        P = gk.unified_kp_strang(S0, tg, 0.6, 1.5, 0.0, 0.0, 1e-12)
        ref = np.array([s - s * gk.kp_kd_exact(s, tg, 1.5, 0.6, 1e12, 0.0) for s in S0])
        self.assertLessEqual(float(np.max(np.abs((S0[:, None] - P) - ref))), 1e-9)


class TestK57(unittest.TestCase):
    def test_fd_jacobian_and_residual(self):
        Y, _, Km, _, _, _ = gen_k57("pkg_k57_fdcheck", 4, 1)
        rng = rng_for("pkg_k57_fdcheck_point")
        s0n = np.broadcast_to(S0_NOM[None, :], (4, 11)).copy()
        P = np.stack([np.full(4, math.log(0.3)), np.log(Km), np.full(4, 0.01), np.full(4, 1 / 16), 1 / (5 * Km)], 1)
        U = rng.normal(0.0, 0.04, size=(4, 11))
        errs = gk.k57_fd_check(Y, s0n, P, U, T_GRID, SIG)
        self.assertEqual(len(errs), 16)
        self.assertLessEqual(max(errs), 1e-6)
        worst = 0.0
        for mode in (0, 1):
            _, S0, Kmr, Vm, kd, ik = gen_k57(f"pkg_k57_resid_m{mode}", 150, mode)
            worst = max(worst, gk.ksi_residual(S0, T_GRID, np.log(Vm), np.log(Kmr), kd, ik))
        self.assertLessEqual(worst, 1e-8)

    def test_map_fit_runs_and_calls(self):
        Y, _, _, _, _, _ = gen_k57("pkg_k57_fit_pos", 8, 1)
        s0n = np.broadcast_to(S0_NOM[None, :], (8, 11)).copy()
        r = gk.k57_map_fit(Y, s0n, T_GRID, SIG)
        self.assertEqual(r["called"].shape, (8,))
        self.assertTrue(np.all(np.isfinite(r["dbic"])))
        self.assertGreaterEqual(int(np.sum(r["called"])), 6)   # Ksi = 5 Km positives (MODELED)


class TestK57DispersionGuard(unittest.TestCase):
    """v0.2.0 opt-in Birge guard (Q65 V2). Reference: tests/data/k57_guard_ref.json, built from the released v0.1.0
    wheel by benchmarks/make_k57_guard_ref.py before solvers.py changed. All plates MODELED, crc32 seeds."""

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(DATA, "k57_guard_ref.json"), encoding="utf-8") as fh:
            cls.ref = json.load(fh)

    # Float tolerance: the reference was made on Windows; Linux numpy/scipy differ around the 14th digit (seen in CI,
    # relative ~5e-12). Numbers are compared to rtol 1e-10; the called flags stay exact.
    RTOL = 1e-10
    # Fitted iksi / km sit on a flat objective; the optimizer stops at slightly different points across platforms
    # and even between CI runs (absolute up to ~3e-8, relative ~2e-8 seen). Relative 1e-6 plus absolute 1e-8.
    ATOL_FIT = 1e-8
    RTOL_FIT = 1e-6

    def test_guard_off_matches_v010(self):
        for tag, B, mode in (("pkg_k57_guard_null", 6, 0), ("pkg_k57_guard_pos", 4, 1)):
            Y = gen_k57(tag, B, mode)[0]
            s0n = np.broadcast_to(S0_NOM[None, :], (B, 11)).copy()
            for r in (gk.k57_map_fit(Y, s0n, T_GRID, SIG), gk.k57_map_fit(Y, s0n, T_GRID, SIG, dispersion_guard=False)):
                for k in ("dbic", "iksi", "km"):
                    atol, rtol = (0.0, self.RTOL) if k == "dbic" else (self.ATOL_FIT, self.RTOL_FIT)
                    np.testing.assert_allclose(r[k], self.ref["v010"][tag][k], rtol=rtol, atol=atol, err_msg=f"{tag} {k}")
                self.assertEqual([bool(x) for x in r["called"]], self.ref["v010"][tag]["called"])
                self.assertTrue(np.array_equal(r["birge"], np.ones(B)))

    def test_guard_on_matches_q65_birge_reference(self):
        q = self.ref["q65"]
        Y = np.array(q["Y"])
        s0 = np.array(q["s0_nom"])
        nom = np.broadcast_to(s0[None, :], Y.shape[:2]).copy()
        sig = 0.015 * np.maximum(s0, 0.5)                       # Q65 nominal SIG (ASSUMED)
        r = gk.k57_map_fit(Y, nom, np.linspace(0.0, 10.0, 13), sig, dispersion_guard=True)
        np.testing.assert_allclose(r["birge"], q["birge"], rtol=self.RTOL, atol=0)
        self.assertTrue(np.all(r["birge"] > 1.0))              # 3x noise floor (MODELED)
        np.testing.assert_allclose(r["dbic"], q["v2"]["dbic"], rtol=self.RTOL, atol=1e-9)
        self.assertEqual([bool(x) for x in r["called"]], q["v2"]["called"])
        self.assertTrue(any(q["row0"]["called"]))              # Row 0 false calls exist on these plates ...
        self.assertFalse(any(bool(x) for x in r["called"]))    # ... and the guard suppresses them
        for i, c in enumerate(q["row0"]["called"]):
            if c:
                self.assertLess(r["dbic"][i], q["row0"]["dbic"][i])

    def test_low_noise_plate_unchanged(self):
        rng = rng_for("pkg_k57_guard_lownoise")
        B = 4
        Km = 10.0 ** rng.uniform(-0.2, 0.35, size=B)
        S0 = S0_NOM[None, :] * np.clip(rng.normal(1.0, 0.05, size=(B, 11)), 0.7, 1.3)
        S, _ = gk.ksi_solve_conv(S0, T_GRID, np.log(rng.uniform(0.18, 0.45, size=B)), np.log(Km),
                                 rng.uniform(0.0, 0.025, size=B), np.zeros(B))
        Y = gk.solvers.ife_forward(S, 1.0 / 16.0) + rng.normal(0.0, 0.5 * SIG[None, :, None], size=S.shape)
        s0n = np.broadcast_to(S0_NOM[None, :], (B, 11)).copy()
        r0 = gk.k57_map_fit(Y, s0n, T_GRID, SIG)
        r1 = gk.k57_map_fit(Y, s0n, T_GRID, SIG, dispersion_guard=True)
        self.assertTrue(np.array_equal(r1["birge"], np.ones(B)))
        for k in ("dbic", "iksi", "km", "called"):
            self.assertTrue(np.array_equal(r0[k], r1[k]), k)

    def test_custom_empirical_sig_accepted(self):
        Y = gen_k57("pkg_k57_guard_empsig", 3, 0)[0]
        s0n = np.broadcast_to(S0_NOM[None, :], (3, 11)).copy()
        sig_emp = 0.02 * np.abs(Y[:, :, 0]).mean(axis=0) + 0.0075   # V1-style blank estimate (ASSUMED form)
        for guard in (False, True):
            r = gk.k57_map_fit(Y, s0n, T_GRID, sig_emp, dispersion_guard=guard)
            for k in ("dbic", "iksi", "km", "called", "birge"):
                self.assertEqual(r[k].shape, (3,), k)
            self.assertTrue(np.all(np.isfinite(r["dbic"])) and np.all(r["birge"] >= 1.0))


class TestDetectors(unittest.TestCase):
    def test_dsf_null_and_signal(self):
        T = np.linspace(25, 95, 71)
        rng = rng_for("pkg_dsf")
        base = 0.02 * T + 1e-4 * T ** 2 + 0.01 * rng.standard_t(4, 71)
        self.assertLess(gk.z_dsf_asym_boltz(T, base), 5.72)
        self.assertGreaterEqual(gk.z_dsf_asym_boltz(T, base + 0.3 / (1 + np.exp(-(T - 55) / 2.5))), 5.72)

    def test_zrao_single_well_is_degenerate(self):
        Y = np.linspace(1.0, 0.6, 13)[None, None, :]
        z = gk.z_rao_orth_robust(Y, np.array([[1.0]]), T_GRID, np.zeros((1, 1)))
        self.assertAlmostEqual(float(z[0]), 0.0, places=6)


class TestParser(unittest.TestCase):
    KW = dict(scale=0.05, offset=0.040, mode="substrate", evidence="MODELED")

    def test_vendored_parser_hash(self):
        p = os.path.join(os.path.dirname(gk.__file__), "_plate_parser.py")
        with open(p, "rb") as fh:
            self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), PP_SHA256)

    def test_reports_byte_identical(self):
        # Gen5 (BioTek) and Tecan i-control preambles above the well table are skipped; long TSV also read.
        for stem, fname in (("plate_healthy", "plate_healthy_biotek_wide.csv"),
                            ("plate_dead", "plate_dead_tecan_transposed.csv"),
                            ("plate_coldC6", "plate_coldC6_long_tsv.tsv")):
            rep, _ = parser.parse_file(os.path.join(DATA, fname), **self.KW)
            with open(os.path.join(DATA, stem + ".report.json"), encoding="utf-8") as fh:
                self.assertEqual(parser.dumps(rep), fh.read(), fname)

    def test_error_codes(self):
        for fname, code in (("bad_empty.csv", "EMPTY"), ("bad_no_wells.csv", "NO_WELL_TABLE"),
                            ("bad_time_order.csv", "TIME_NOT_MONOTONIC")):
            rep, _ = parser.parse_file(os.path.join(DATA, fname), **self.KW)
            self.assertEqual(rep["status"], "ERROR")
            self.assertEqual(rep["error"]["code"], code)

    def test_cli_writes_reports(self):
        with tempfile.TemporaryDirectory() as d:
            buf = io.StringIO()
            with redirect_stdout(buf):
                rc = cli.main([os.path.join(DATA, "plate_healthy_biotek_wide.csv"), "--evidence", "MODELED",
                               "--scale", "0.05", "--offset", "0.040", "--out", d])
            self.assertEqual(rc, 0)
            with open(os.path.join(d, "plate_healthy_biotek_wide.report.json"), encoding="utf-8") as fh:
                self.assertEqual(json.load(fh)["status"], "OK")
            self.assertTrue(os.path.exists(os.path.join(d, "plate_healthy_biotek_wide.report.md")))
            with redirect_stdout(io.StringIO()):
                self.assertEqual(cli.main([os.path.join(DATA, "bad_no_wells.csv"), "--out", d]), 2)


if __name__ == "__main__":
    unittest.main()
