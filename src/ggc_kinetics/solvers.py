"""Frozen zero-ODE kinetic solvers from the GGC-inf archive (Cycles 10-66), packaged without changes to the math.

Every function here is a copy of an archived, regression-pinned operator, and its source is named. Evidence tags
(L7): the solvers are exact or converged numerics. The statements about when they are valid (E0 known to +-10%,
A_peak <= 2.00, pipetting CV <= 5%) come from MODELED studies and are not yet confirmed on real plates.
"""
import math

import numpy as np
from scipy.optimize import least_squares
from scipy.special import lambertw

LN10 = math.log(10.0)

__all__ = ["kp_kd_exact", "kp_warped_step", "unified_kp_strang", "coupled_lag_decay", "slow_binding_warped",
           "morrison_G", "morrison_ES", "morrison_exact", "ife_forward", "ife_inverse", "ksi_solve_conv",
           "ksi_residual", "k57_evaluate", "k57_map_fit", "k57_called", "k57_fd_check", "fit_kp_vs_kd",
           "K57_THR", "K57_IKSI", "K57_SP"]


# ----------------------------------------------------------------------------- Lambert-W family (Cycles 10, 20, 26)
def kp_kd_exact(S0, t, Km, Vmax, Kp, kd):
    """Fraction converted (S0 - S)/S0 for dS/dt = -Vmax e^{-kd t} S / (Km (1 + (S0 - S)/Kp) + S).

    Three regimes (Kp > Km, Kp = Km, Kp < Km) on the principal Lambert-W branch; Kp > 1e8 is plain
    Michaelis-Menten (the Schnell-Mendoza solution). Verbatim from the regression suite (Cycles 10 and 30).
    """
    tau = (1.0 - np.exp(-kd * t)) / kd if kd > 1e-9 else t
    if Kp > 1e8:
        u = (S0 / Km) * np.exp((S0 - Vmax * tau) / Km)
        S = Km * np.real(lambertw(u))
        return np.clip((S0 - S) / S0, 0.0, 1.0)
    alpha = 1.0 - Km / Kp
    beta = Km * (1.0 + S0 / Kp)
    if abs(alpha) < 1e-6:
        S = S0 * np.exp(-Vmax * tau / beta)
    elif alpha > 0.0:
        c = alpha / beta
        u = c * S0 * np.exp(c * S0 - Vmax * tau / beta)
        S = np.real(lambertw(u, k=0)) / c
    else:
        c = alpha / beta
        u = c * S0 * np.exp(c * S0 - Vmax * tau / beta)
        u = np.maximum(u, -math.exp(-1.0) + 1e-15)
        S = np.real(lambertw(u, k=0)) / c
    return np.clip((S0 - S) / S0, 0.0, 1.0)


def fit_kp_vs_kd(S0, t, frac, start_kp=(16.0, 0.75, 35.0), start_kd=(20.0, 0.75, 0.001)):
    """Fit product inhibition (Km, Vmax, Kp) and enzyme decay (Km, Vmax, kd) to fraction-converted data and
    compare them by BIC (dbic < 0 favours product inhibition). The suite's R1 test on the 1913 invertase data."""
    S0, t, frac = (np.asarray(a, float) for a in (S0, t, frac))
    n = len(frac)
    r_kp = least_squares(lambda p: kp_kd_exact(S0, t, *np.exp(p[:2]), math.exp(p[2]), 0.0) - frac,
                         np.log(start_kp), method="lm")
    r_kd = least_squares(lambda p: kp_kd_exact(S0, t, *np.exp(p[:2]), 1e12, math.exp(p[2])) - frac,
                         np.log(start_kd), method="lm")
    Km, Vmax, Kp = np.exp(r_kp.x)
    rss_kp, rss_kd = float(np.sum(r_kp.fun ** 2)), float(np.sum(r_kd.fun ** 2))
    bic_kp = n * math.log(rss_kp / n) + 3 * math.log(n)
    bic_kd = n * math.log(rss_kd / n) + 3 * math.log(n)
    return {"Km": float(Km), "Vmax": float(Vmax), "Kp": float(Kp), "eff": float(Vmax / Km),
            "dbic": float(bic_kp - bic_kd), "kd_fit": [float(v) for v in np.exp(r_kd.x)]}


def kp_warped_step(S_in, S0_tot, dtau, Vmax, Km, inv_Kp):
    """Exact step over warped time dtau of -dS/dtau = Vmax S / (Km (1 + (S0_tot - S) inv_Kp) + S).
    Verbatim from tier0_kp_kd_cycle10.py."""
    if Vmax <= 0.0 or dtau <= 0.0:
        return S_in.copy()
    alpha = 1.0 - Km * inv_Kp
    if abs(alpha) < 1e-6:
        denom = Km * (1.0 + S0_tot * inv_Kp)
        return S_in * np.exp(-(Vmax * dtau) / np.maximum(denom, 1e-12))
    elif alpha > 0.0:
        K_eff = Km * (1.0 + S0_tot * inv_Kp) / alpha
        V_eff = Vmax / alpha
        u = S_in / K_eff
        arg = np.clip(u - (V_eff * dtau) / K_eff, -500.0, 500.0)
        z = u * np.exp(arg)
        log_z = arg + np.log(np.maximum(u, 1e-30))
        w_asym = log_z - np.log(np.maximum(log_z, 1.0))
        w_val = np.where(arg > 40.0, w_asym, np.real(lambertw(np.maximum(z, 0.0))))
        return K_eff * w_val
    else:
        beta = -alpha
        K_neg = Km * (1.0 + S0_tot * inv_Kp) / beta
        V_neg = Vmax / beta
        u = S_in / K_neg
        arg = np.clip(-(u + (V_neg * dtau) / K_neg), -500.0, 0.0)
        z = -u * np.exp(arg)
        w_val = np.real(lambertw(np.clip(z, -math.exp(-1.0) + 1e-15, 0.0)))
        return -K_neg * w_val


def unified_kp_strang(S0, t_grid, Vmax, Km, kd, knon, inv_Kp):
    """Product formed P(t) = S0 - S(t), shape (W, M), for
    dS/dt = -Vmax e^{-kd t} S / (Km (1 + (S0 - S) inv_Kp) + S) - knon S,
    by symmetric (2nd-order) Strang splitting of the exact enzyme step and exact first-order hydrolysis.
    Exact when knon = 0; Cycle 10 measured < 1e-5 vs RK4 for Kp + kd and < 0.005 with knon (MODELED).
    Verbatim from tier0_kp_kd_cycle10.unified_kp_strang_predict."""
    S0 = np.asarray(S0, dtype=float)
    t = np.asarray(t_grid, dtype=float)
    W, M = S0.shape[0], t.shape[0]
    S_curr = S0.copy()
    P_out = np.empty((W, M), dtype=float)
    t_prev = 0.0
    for j in range(M):
        tj = t[j]
        dt = tj - t_prev
        if dt > 0.0:
            eh = math.exp(-0.5 * knon * dt) if knon > 0.0 else 1.0
            S_mid = S_curr * eh
            dtau = (math.exp(-kd * t_prev) - math.exp(-kd * tj)) / kd if kd > 1e-9 else dt
            S_mid2 = kp_warped_step(S_mid, S0, dtau, Vmax, Km, inv_Kp)
            S_curr = S_mid2 * eh
        P_out[:, j] = S0 - S_curr
        t_prev = tj
    return P_out


def coupled_lag_decay(t, v0, kd, k2):
    """Coupled-assay product P(t) when the primary rate decays, v(t) = v0 e^{-kd t}, and the intermediate is
    turned over at k2: dI/dt = v - k2 I, dP/dt = k2 I. Exact closed form (Cycle 20; suite R3, |exact - RK4| < 1e-11).
    Requires kd > 0 and k2 != kd."""
    t = np.asarray(t, dtype=float)
    return (v0 / kd) * (1.0 - np.exp(-kd * t)) - (v0 / (k2 - kd)) * (np.exp(-kd * t) - np.exp(-k2 * t))


def slow_binding_warped(t, S0, Km, Vmax, k_obs, e_ss):
    """S(t) for slow-onset inhibition: the active-enzyme fraction relaxes 1 -> e_ss at rate k_obs, so
    tau = e_ss t + (1 - e_ss)/k_obs (1 - e^{-k_obs t}) and S = Km W((S0/Km) exp((S0 - Vmax tau)/Km)).
    Cycle 26; suite R4 (|warped LW - 2-ODE RK4| < 1e-8)."""
    t = np.asarray(t, dtype=float)
    tau = e_ss * t + ((1.0 - e_ss) / k_obs) * (1.0 - np.exp(-k_obs * t))
    u = (S0 / Km) * np.exp((S0 - Vmax * tau) / Km)
    return Km * np.real(lambertw(u))


# ----------------------------------------------------------------------------- Morrison tight binding (Cycle 32)
def morrison_G(S, E0, Km):
    """Antiderivative G with kcat t = G(S0) - G(S) for the quadratic (Morrison) rate law. Verbatim, Cycle 32."""
    p = Km - E0
    q = Km + E0
    R = np.sqrt(S * S + 2 * p * S + q * q)
    return (2 * q * np.log(S) + S + R + p * np.log(R + S + p) - q * np.log(q * q + p * S + q * R)) / (2 * E0)


def morrison_ES(S, E0, Km):
    """[ES] from the quadratic, using the stable small root 2 E0 S / (b + sqrt(b^2 - 4 E0 S))."""
    b = E0 + S + Km
    return 2.0 * E0 * S / (b + np.sqrt(np.maximum(b * b - 4 * E0 * S, 0.0)))


def morrison_exact(S0, t_grid, kcat, E0, Km):
    """Exact S(t) when [E0] is comparable to Km (no free-substrate approximation): -dS/dt = kcat [ES](S).

    Shapes: S0 (B, W); t_grid (M,); kcat, E0, Km (B,). Returns (S (B, W, M), max implicit-equation residual).
    Newton in u = ln S; the archive pins residual <= 1e-10 and |S - RK4|/S0 <= 1e-8 (suite R9).

    VALIDITY (K66, MODELED): the closure needs the active-site E0 known to better than about +-10%. With E0
    misspecified by x0.9 the Km error was 0.082 dex, and by x1.1 it was 0.134 dex (at E0 = 1.5 Km). Measure E0
    by an independent active-site titration, or fit it as a parameter.
    """
    S03 = S0[:, :, None]
    kc, e0, km = kcat[:, None, None], E0[:, None, None], Km[:, None, None]
    t3 = t_grid[None, None, :]
    G0 = morrison_G(S03, e0, km)
    v0 = kc * morrison_ES(S03, e0, km)
    u = np.log(S03) - (v0 / S03) * t3
    u = np.maximum(u, np.log(S03) - 60.0)
    lnS0 = np.log(S03)
    for _ in range(40):
        S = np.exp(u)
        Hf = G0 - morrison_G(S, e0, km) - kc * t3
        du = Hf * morrison_ES(S, e0, km) / S
        u = np.clip(u + du, lnS0 - 60.0, lnS0)
        if np.max(np.abs(du)) < 1e-13:
            break
    S = np.exp(u)
    resid = np.max(np.abs(G0 - morrison_G(S, e0, km) - kc * t3) * kc / np.maximum(kc, 1e-300))
    return S, resid


# ----------------------------------------------------------------------------- K57 MAP rule (Cycles 37, 43, 46)
K57_THR, K57_IKSI, K57_SP = -2.5, 0.045, 0.05
K57_LO = np.array([-7.0, -3.5, 0.0, 0.0, 0.0])
K57_HI = np.array([3.0, 4.0, 0.30, 0.12, 2.0])
K57_U_LIM = 0.35


def ife_forward(S, a):
    """Inner-filter readout Y = (1 - e^{-a ln10 S}) / (a ln10); Y = S when a = 0."""
    k = np.asarray(a) * LN10
    ks = np.maximum(k, 1e-12)
    return np.where(k > 1e-12, -np.expm1(-ks * S) / ks, S)


def ife_inverse(Y, a):
    k = np.asarray(a) * LN10
    ks = np.maximum(k, 1e-12)
    arg = np.clip(1.0 - ks * Y, 1e-12, None)
    return np.where(k > 1e-12, -np.log(arg) / ks, Y)


def ksi_solve_conv(S0, t_grid, log_vm, log_km, kd, iksi):
    """Substrate-inhibition progress curve with decay: Km ln(S0/S) + (S0 - S) + iksi (S0^2 - S^2)/2 = Vmax tau,
    tau = (1 - e^{-kd t})/kd. Newton in x = ln S (8 steps; Cycle 43). Returns S (B, W, M) and
    dS/d[ln Vm, ln Km, kd, iksi] (B, W, M, 4). Verbatim from sprint28_cycle43_k57_jufix.ksi_solve_conv."""
    Vm = np.exp(log_vm)[:, None, None]
    Km = np.exp(log_km)[:, None, None]
    kd3 = kd[:, None, None]
    ik = iksi[:, None, None]
    t3 = t_grid[None, None, :]
    S03 = S0[:, :, None]
    kt = kd3 * t3
    big = kd3 > 1e-6
    kds = np.maximum(kd3, 1e-12)
    tau = np.where(big, -np.expm1(-kt) / kds, t3 * (1.0 - 0.5 * kt))
    dtau = np.where(big, (kt * np.exp(-kt) + np.expm1(-kt)) / (kds * kds), -0.5 * t3 * t3)
    v0 = Vm * S03 / (Km + S03 + ik * S03 * S03)
    S = np.minimum(np.maximum(S03 * np.exp(-(v0 / np.maximum(S03, 1e-9)) * tau), 1e-9), S03)
    lS03 = np.log(S03)
    x = np.log(S)
    for _ in range(8):
        S = np.exp(x)
        F = Km * (lS03 - x) + (S03 - S) + 0.5 * ik * (S03 * S03 - S * S) - Vm * tau
        x = np.minimum(x + F / (Km + S + ik * S * S), lS03)
    S = np.exp(x)
    inv = 1.0 / (Km / S + 1.0 + ik * S)
    J = np.stack([(-Vm * tau) * inv,
                  (Km * np.log(np.maximum(S03 / S, 1e-30))) * inv,
                  (-Vm * dtau) * inv,
                  (0.5 * (S03 * S03 - S * S)) * inv], axis=-1)
    return S, J


def ksi_residual(S0, t_grid, log_vm, log_km, kd, iksi):
    """max |F| of the implicit equation at the ksi_solve_conv solution (exported for audits; gate <= 1e-8)."""
    S, _ = ksi_solve_conv(S0, t_grid, log_vm, log_km, kd, iksi)
    t3 = t_grid[None, None, :]
    kd3 = kd[:, None, None]
    kt = kd3 * t3
    tau = np.where(kd3 > 1e-6, -np.expm1(-kt) / np.maximum(kd3, 1e-12), t3 * (1.0 - 0.5 * kt))
    K3, V3 = np.exp(log_km)[:, None, None], np.exp(log_vm)[:, None, None]
    I3, S03 = iksi[:, None, None], S0[:, :, None]
    F = K3 * np.log(S03 / S) + (S03 - S) + 0.5 * I3 * (S03 * S03 - S * S) - V3 * tau
    return float(np.max(np.abs(F)))


def k57_evaluate(Y, s0n, P, U, t_grid, sig, sp=K57_SP):
    """Residuals, kinetic Jacobian Jg (B, W, M, 5), per-well log-S0 Jacobian Ju (B, W, M; with the Cycle 43
    Ju x S0 fix), and MAP cost with the per-well log-normal S0 prior of width sp.
    P columns: ln Vmax, ln Km, kd, a (IFE absorbance per norm unit), iksi = 1/Ksi. U: per-well ln(S0/S0_nominal)."""
    S0 = s0n * np.exp(U)
    S, JS = ksi_solve_conv(S0, t_grid, P[:, 0], P[:, 1], P[:, 2], P[:, 4])
    k = (P[:, 3] * LN10)[:, None, None]
    ks = np.maximum(k, 1e-12)
    x = ks * S
    Yp = np.where(k > 1e-12, -np.expm1(-x) / ks, S)
    g = np.exp(-k * S)
    dYdk = np.where(x > 1e-4, (x * np.exp(-x) + np.expm1(-x)) / (ks * ks), -0.5 * S * S + k * S ** 3 / 3.0)
    JY = JS * g[..., None]
    sg = sig[None, :, None]
    Jg = np.concatenate([JY[..., :3], (dYdk * LN10)[..., None], JY[..., 3:4]], axis=-1) / sg[..., None]
    Km = np.exp(P[:, 1])[:, None, None]
    ik = P[:, 4][:, None, None]
    S03 = S0[:, :, None]
    Ju = g * (Km + S03 + ik * S03 * S03) / (Km / S + 1.0 + ik * S) / S03 / sg
    Ju = Ju * S03                                        # Cycle 43: dY/du = dY/dS0 * S0
    r = (Yp - Y) / sg
    cost = np.sum(r * r, axis=(1, 2)) + np.sum((U / sp) ** 2, axis=1)
    return r, Jg, Ju, cost


def _k57_lm(Y, s0n, P0, free, t_grid, sig, sp, n_iter, U0=None):
    """Levenberg-Marquardt with the per-well u eliminated by a Schur complement (Cycle 37 lm, 'local' prior)."""
    B, W = s0n.shape
    P = P0.copy()
    U = np.zeros((B, W)) if U0 is None else U0.copy()
    fm = np.array(free, dtype=float)
    r, Jg, Ju, c = k57_evaluate(Y, s0n, P, U, t_grid, sig, sp)
    lam = np.full(B, 1e-2)
    for _ in range(n_iter):
        Jgm = Jg * fm
        A = np.einsum('bwmi,bwmj->bij', Jgm, Jgm)
        gg = np.einsum('bwmi,bwm->bi', Jgm, r)
        dA = np.maximum(np.diagonal(A, axis1=1, axis2=2), 1e-8)
        Ad = A + np.eye(5)[None] * (lam[:, None] * dA + 1e-9)[:, None, :] + np.diag(1.0 - fm)[None]
        Bm = np.einsum('bwmi,bwm->biw', Jgm, Ju)
        D = np.sum(Ju * Ju, axis=2) + 1.0 / sp ** 2
        gu = np.sum(Ju * r, axis=2) + U / sp ** 2
        D = D * (1.0 + lam[:, None])
        Sch = Ad - np.einsum('biw,bjw->bij', Bm / D[:, None, :], Bm)
        rhs = -gg + np.einsum('biw,bw->bi', Bm, gu / D)
        dg = np.linalg.solve(Sch, rhs[:, :, None])[:, :, 0] * fm
        du = (-gu - np.einsum('biw,bi->bw', Bm, dg)) / D
        Pc = np.clip(P + dg, K57_LO, K57_HI)
        Uc = np.clip(U + du, -K57_U_LIM, K57_U_LIM)
        rc, Jgc, Juc, cc = k57_evaluate(Y, s0n, Pc, Uc, t_grid, sig, sp)
        acc = cc < c
        P = np.where(acc[:, None], Pc, P)
        U = np.where(acc[:, None], Uc, U)
        r = np.where(acc[:, None, None], rc, r)
        Jg = np.where(acc[:, None, None, None], Jgc, Jg)
        Ju = np.where(acc[:, None, None], Juc, Ju)
        c = np.where(acc, cc, c)
        lam = np.where(acc, np.maximum(lam * 0.3, 1e-7), np.minimum(lam * 5.0, 1e7))
    return P, c, U


def k57_map_fit(Y, s0n, t_grid, sig, sp=K57_SP, n_iter=24, dispersion_guard=False):
    """Frozen K57 MAP substrate-inhibition test (Cycle 46 rule), generalized from the 11-well x 13-read design to
    any (B, W, M) plate. Y in normalized units (1 norm unit = the assay's S0 scale), s0n (B, W) nominal S0,
    t_grid (M,) starting at 0, sig (W,) readout SD. Returns dict with dbic, iksi, km, called.

    Rule: called = dbic > -2.5 and iksi > 0.045 (ASSUMED thresholds, calibrated on 1,200 MODELED null plates).
    VALIDITY (MODELED only, not yet confirmed on a real plate): pipetting CV <= 5% (the 5% log-normal S0 prior),
    IFE peak absorbance A_peak <= 2.00 (Cycle 62, R38), and a 12-13-read grid. The initial-rate start uses
    reads 0 and 2, so M >= 3. Code: sprint22_cycle37_k57_s0prior.fit with the Cycle 43 fixes (Cycle 46 install).
    NOISE SCALE (Q64/Q65, MODELED): the false-alarm bound needs sig close to the plate's real readout SD. A noise
    floor 2x or 3x above sig raised false calls to w90 0.076-0.131. Prefer sig estimated from blank or replicate
    wells; otherwise set dispersion_guard=True.

    dispersion_guard (opt-in, default False = v0.1.0 behaviour, bit-identical): one-sided Birge guard (Q65 V2).
    From the first null stage, chi2_null = cost - sum((U0 / sp)^2), dof_null = max(W*M - 4 - W, 1), and
    birge = max(1, sqrt(chi2_null / dof_null)). Each plate with birge > 1 is refitted at sig * birge (thresholds and
    sp unchanged). The returned dict always carries "birge" (all 1.0 when the guard is off).
    """
    Y = np.asarray(Y, float)
    s0n = np.asarray(s0n, float)
    t_grid = np.asarray(t_grid, float)
    sig = np.asarray(sig, float)
    B = Y.shape[0]
    ag = np.linspace(0.0, 0.10, 61)
    chi = np.sum(((Y[:, :, 0][:, :, None] - ife_forward(s0n[:, :, None], ag[None, None, :]))
                  / sig[None, :, None]) ** 2, axis=1)
    a0 = ag[np.argmin(chi, axis=1)]
    Sr = ife_inverse(Y, a0[:, None, None])
    v = np.maximum((Sr[:, :, 0] - Sr[:, :, 2]) / (t_grid[2] - t_grid[0]), 1e-4)
    P0 = np.stack([np.log(np.clip(np.max(v, axis=1) * 1.25, 1e-3, 10.0)), np.zeros(B),
                   np.full(B, 0.015), a0, np.zeros(B)], axis=1)
    f0, f1 = [1, 1, 1, 1, 0], [1, 1, 1, 1, 1]
    Ph0, c0, U0 = _k57_lm(Y, s0n, P0, f0, t_grid, sig, sp, n_iter)
    c0_first = c0
    P1s = Ph0.copy()
    P1s[:, 4] = 0.05
    Ph1, c1, U1 = _k57_lm(Y, s0n, P1s, f1, t_grid, sig, sp, n_iter, U0)
    P0b = Ph1.copy()
    P0b[:, 4] = 0.0
    Ph0b, c0b, _ = _k57_lm(Y, s0n, P0b, f0, t_grid, sig, sp, n_iter, U1)
    Ph0 = np.where((c0b < c0)[:, None], Ph0b, Ph0)
    c0 = np.minimum(c0, c0b)
    better1 = c1 <= c0
    c1 = np.minimum(c1, c0)
    n = Y.shape[1] * Y.shape[2]
    dbic = c0 - c1 - math.log(n)
    iksi = np.where(better1, Ph1[:, 4], 0.0)
    out = {"dbic": dbic, "iksi": iksi, "km": np.exp(np.where(better1, Ph1[:, 1], Ph0[:, 1])),
           "called": k57_called(dbic, iksi), "birge": np.ones(B)}
    if dispersion_guard:
        dof = max(Y.shape[1] * Y.shape[2] - 4 - Y.shape[1], 1)
        chi2 = c0_first - np.sum((U0 / sp) ** 2, axis=1)
        birge = np.maximum(1.0, np.sqrt(np.maximum(chi2, 0.0) / dof))
        out["birge"] = birge
        for i in np.flatnonzero(birge > 1.0):
            ri = k57_map_fit(Y[i:i + 1], s0n[i:i + 1], t_grid, sig * birge[i], sp, n_iter)
            for k in ("dbic", "iksi", "km", "called"):
                out[k][i] = ri[k][0]
    return out


def k57_called(dbic, iksi, thr=K57_THR, ikm=K57_IKSI):
    return (np.asarray(dbic) > thr) & (np.asarray(iksi) > ikm)


def k57_fd_check(Y, s0n, P, U, t_grid, sig, sp=K57_SP, h=1e-6):
    """Max relative central-difference error of the exported analytic Jacobian, per column (5 kinetic + W u)."""
    def res(P_, U_):
        return k57_evaluate(Y, s0n, P_, U_, t_grid, sig, sp)[0]
    _, Jg, Ju, _ = k57_evaluate(Y, s0n, P, U, t_grid, sig, sp)
    errs = []
    for i in range(5):
        hi = h * max(1.0, float(np.max(np.abs(P[:, i]))))
        Pp, Pm = P.copy(), P.copy()
        Pp[:, i] += hi
        Pm[:, i] -= hi
        f = (res(Pp, U) - res(Pm, U)) / (2 * hi)
        errs.append(float(np.max(np.abs(f - Jg[..., i])) / np.max(np.abs(f))))
    for w in range(s0n.shape[1]):
        Up, Um = U.copy(), U.copy()
        Up[:, w] += h
        Um[:, w] -= h
        f = (res(P, Up) - res(P, Um)) / (2 * h)
        an = np.zeros_like(f)
        an[:, w, :] = Ju[:, w, :]
        errs.append(float(np.max(np.abs(f - an)) / np.max(np.abs(f))))
    return errs
