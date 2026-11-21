"""
plate_parser.py  (GGC-inf Cycle 33; SoftMax reader and --layout added in Cycle 40)  --  96-well kinetic plate export -> Cycle 29 contract -> P1-P5 verdicts.

Usage
  python plate_parser.py PLATE.csv [--scale A_PER_NORM] [--offset A0] [--mode substrate|product]
                         [--evidence MEASURED|MODELED] [--out-dir DIR]
Writes PLATE.report.json and PLATE.report.md. On a file it cannot read it writes a JSON error
report and exits with code 2 (no traceback).

What it reads (auto-detected, no vendor name needed)
  * wide        one row per read: a time column plus well columns A1..H12 (BioTek/Gen5-like)
  * transposed  one row per well: a "Time" row, then rows A1..H12 (Tecan/Magellan-like)
  * long        one row per (well, time, value) with a header naming Well / Time / value
  Encodings utf-8, utf-8 BOM, utf-16, cp1252. Delimiters comma, semicolon, tab. Decimal commas.
  Times as h:mm:ss, mm:ss, or numbers with a unit in the header ([s], (min), h).
  * softmax_kinetic  SoftMax Pro plate-format kinetic text export (##BLOCKS= header, one 8 x 12 block per
                read per wavelength; Cycle 40). --wavelength-index picks the wavelength (0 = first listed).
  --layout MAP.json  maps contract columns 1-12 to source wells and may declare its own 8-row S0 ladder
                (Cycle 40). A ladder other than the contract ladder is reported as design_match
                CUSTOM_LADDER: verdicts are computed, but their thresholds were only calibrated on the contract
                ladder (ASSUMED transfer).
  Well names A1 / A01 / a 1. Overflow and other non-numeric cells become NaN and are counted;
  interior gaps are filled by linear interpolation in time and reported.

The Cycle 29 contract (ASSUMED, from sprint14_cycle29_wetlab_blueprint.py)
  rows A..H = S0 ladder [12.5, 25, 50, 100, 200, 350, 550, 800] uM, normalised by 50 uM
  C1 cols 1-2 baseline | C2 col 3 buffer blank ladder | C3 cols 4-5 BSA matrix (2x hydrolysis)
  C4 cols 6-7 inert dye IFE | C5 cols 8-9 product inhibition (parsed, not scored: the blueprint has no
  test for it) | C6 col 10 thermal decay | C7a col 11 heat-killed null | C7b col 12 dye + drift null
  The readout is converted to normalised substrate units: Y = (A - offset) / scale (mode substrate),
  or Y = S0_nominal - (A - A(t0)) / scale (mode product). scale and offset are ASSUMED inputs.
  The first read is taken as t = 0 (ASSUMED); its original time is reported as t_offset_s.

Single-plate rules (ASSUMED thresholds; their power is MODELED in sprint18_cycle33_plate_parser.py)
  P1 baseline   PASS iff C1 is active (Z_Rao_Orth_Robust >= 4.60), fitted Km lies in
                [S0_min/2, 2*S0_max], Vmax is off its fit bound, and the two C1 replicate columns
                agree on Km and on Vmax/Km: |ln ratio| <= 0.277 (rule v2, default since Cycle 34;
                rule v1 = ln 1.15 is kept as --p1-rule v1).
  P2 IFE        PASS iff IFE-calibrated C4 Km is within 20% of C1 Km.
  P3 thermal    PASS iff kd(C6) - kd(C1) > 0.03 /min.
  P4 matrix     PASS iff the blank-only hydrolysis correction leaves C3 Km within 15% of C1 Km.
                KILL means the hard-blank correction is falsified on this plate (the designed
                2x-hydrolysis C3 is expected to KILL it).
  P5 nulls      PASS iff z(C1) >= 4.60 and z(C7a) < 4.60 and z(C7b) < 4.60.
  P2-P4 are INCONCLUSIVE when P1 does not PASS (no valid baseline). Any test whose wells are
  missing or more than 25% unreadable is INCONCLUSIVE.

Numerics are copied verbatim from sprint2/3/4 of the kit (ife_forward, ife_inverse, Halley
Lambert-W, mm_kd_eval_and_jac, batched_lm_readout_space, z_rao_orthogonal_robust,
conditional_ife_self_cal with a per-plate trigger mask added) so this file needs only numpy.
"""
import csv
import hashlib
import json
import math
import os
import re
import sys

import numpy as np

LN10 = math.log(10.0)
S0_UM = np.array([12.5, 25.0, 50.0, 100.0, 200.0, 350.0, 550.0, 800.0])
NORM_UM = 50.0
S0_NORM = S0_UM / NORM_UM
ROWS = "ABCDEFGH"
Z_THR = 4.60
# P1 replicate-agreement threshold on |ln ratio| (ASSUMED). v1 = ln 1.15 (Cycle 33, uncalibrated, K67);
# v2 = 0.277 (post-hoc in Cycle 33, pre-registered and tested on fresh seeds in Cycle 34).
P1_REP_THR = {"v1": math.log(1.15), "v2": 0.277}
P1_RULE = "v2"   # switched in Cycle 34 after G1-G5 passed on fresh c34_test_* seeds (MODELED)
BLOCKS = {"C1": [0, 1], "C2": [2], "C3": [3, 4], "C4": [5, 6], "C5": [7, 8], "C6": [9], "C7a": [10], "C7b": [11]}
TEST_BLOCKS = {"P1": ["C1", "C2"], "P2": ["C1", "C2", "C4"], "P3": ["C1", "C2", "C6"],
               "P4": ["C1", "C2", "C3"], "P5": ["C1", "C7a", "C7b"]}
SCHEMA = "ggc.plate_report/1"


# ============================================================================
# Numerics (verbatim copies from the kit; see module docstring)
# ============================================================================
def ife_forward(S, a_eff):
    kappa = a_eff * LN10
    ks = kappa * S
    return np.where(ks < 1e-6, S * (1.0 - 0.5 * ks), (1.0 - np.exp(-ks)) / np.maximum(kappa, 1e-12))


def ife_inverse(y, a_eff):
    kappa = a_eff * LN10
    ceil = 1.0 / np.maximum(kappa, 1e-12)
    y_clip = np.clip(y, 0.0, 0.98 * ceil)
    ky = kappa * y_clip
    return np.where(kappa < 1e-7, y_clip, -np.log(np.maximum(1.0 - ky, 1e-4)) / np.maximum(kappa, 1e-12))


def lambert_w0_halley_log(log_u):
    w = np.where(log_u > 1.5, log_u - np.log(np.maximum(log_u, 1.0)),
                 np.where(log_u > -3.0, np.log1p(np.exp(np.clip(log_u, -3.0, 1.5))), np.exp(np.clip(log_u, -40.0, -3.0))))
    w = np.maximum(w, 1e-12)
    for _ in range(4):
        f = w + np.log(w) - log_u
        wp1 = w + 1.0
        step = (w * f) / np.maximum(wp1 + 0.5 * f / wp1, 1e-8)
        w = np.maximum(w - step, 1e-12)
    return w


def mm_kd_eval_and_jac(S0, t_grid, log_vm, log_km, kd):
    Vm = np.exp(log_vm)[:, None, None]
    Km = np.exp(log_km)[:, None, None]
    kd3 = kd[:, None, None]
    t3 = t_grid[None, None, :]
    S03 = S0[:, :, None]
    kt = kd3 * t3
    tau = np.where(kd3 > 1e-6, (1.0 - np.exp(-kt)) / np.maximum(kd3, 1e-12), t3 * (1.0 - 0.5 * kt))
    dtau_dkd = np.where(kd3 > 1e-6, (kt * np.exp(-kt) - (1.0 - np.exp(-kt))) / np.maximum(kd3 * kd3, 1e-12), -0.5 * t3 * t3)
    u = S03 / Km
    log_z = np.log(np.maximum(u, 1e-30)) + u - (Vm * tau) / Km
    w = lambert_w0_halley_log(log_z)
    S = Km * w
    dw_dlogz = w / (1.0 + w)
    dS_dlvm = Km * dw_dlogz * (-(Vm * tau) / Km)
    dlogz_dlkm = -1.0 - u + (Vm * tau) / Km
    dS_dlkm = Km * (w + dw_dlogz * dlogz_dlkm)
    dS_dkd = Km * dw_dlogz * (-(Vm / Km) * dtau_dkd)
    return S, np.stack([dS_dlvm, dS_dlkm, dS_dkd], axis=-1)


def conditional_ife_self_cal(Y_obs, S0_nom, path_mult, chi2_thr=6.63):
    """Kit version, but returns the per-plate trigger mask instead of its mean."""
    B = Y_obs.shape[0]
    y0 = Y_obs[:, :, 0]
    sig0 = 0.015 * np.maximum(S0_nom, 0.5)[None, :]
    a_grid = np.linspace(0.0, 0.10, 61)[None, None, :]
    a_eff = a_grid * path_mult[None, :, None]
    y0_pred = ife_forward(S0_nom[None, :, None], a_eff)
    chi2_all = np.sum(((y0[:, :, None] - y0_pred) / sig0[:, :, None]) ** 2, axis=1)
    chi2_0 = chi2_all[:, 0]
    best_idx = np.argmin(chi2_all, axis=1)
    chi2_best = chi2_all[np.arange(B), best_idx]
    a_best = a_grid[0, 0, best_idx]
    apply_mask = (chi2_0 - chi2_best) > chi2_thr
    a_used = np.where(apply_mask, a_best, 0.0)
    S_corr = ife_inverse(Y_obs, a_used[:, None, None] * path_mult[None, :, None])
    return S_corr, a_used, apply_mask


def batched_lm_readout_space(Y_obs, S0_init, t_grid, a_plate, n_iter=16):
    B, W, M = Y_obs.shape
    w_scale = np.maximum(S0_init, 0.5)[:, :, None]
    kappa = (a_plate * LN10)[:, :, None]
    S_rough = ife_inverse(Y_obs, a_plate[:, :, None])
    v_est = np.maximum((S_rough[:, :, 0] - S_rough[:, :, 2]) / (t_grid[2] - t_grid[0]), 1e-4)
    log_vm = np.log(np.clip(np.max(v_est, axis=1) * 1.25, 1e-3, 10.0))
    log_km = np.zeros(B)
    kd = np.full(B, 0.015)
    lam = np.full(B, 1e-2)
    S_pred, J_S = mm_kd_eval_and_jac(S0_init, t_grid, log_vm, log_km, kd)
    Y_pred = ife_forward(S_pred, a_plate[:, :, None])
    J_Y = J_S * np.exp(-kappa * S_pred)[:, :, :, None]
    res = (Y_pred - Y_obs) / w_scale
    cost = np.sum(res * res, axis=(1, 2))
    for _ in range(n_iter):
        Jw = J_Y / w_scale[:, :, :, None]
        H = np.einsum('bwmi,bwmj->bij', Jw, Jw)
        g = np.einsum('bwmi,bwm->bi', Jw, res)
        diag = np.maximum(np.diagonal(H, axis1=1, axis2=2), 1e-4)
        H_d = H.copy()
        H_d[:, 0, 0] += lam * diag[:, 0]
        H_d[:, 1, 1] += lam * diag[:, 1]
        H_d[:, 2, 2] += lam * diag[:, 2] + 1e-3
        try:
            delta = -np.linalg.solve(H_d, g[:, :, None])[:, :, 0]
        except np.linalg.LinAlgError:
            delta = np.zeros((B, 3))
        cand_lvm = np.clip(log_vm + delta[:, 0], -7.0, 3.0)
        cand_lkm = np.clip(log_km + delta[:, 1], -3.5, 4.0)
        cand_kd = np.clip(kd + delta[:, 2], 0.0, 0.30)
        S_c, J_Sc = mm_kd_eval_and_jac(S0_init, t_grid, cand_lvm, cand_lkm, cand_kd)
        Y_c = ife_forward(S_c, a_plate[:, :, None])
        res_c = (Y_c - Y_obs) / w_scale
        cost_c = np.sum(res_c * res_c, axis=(1, 2))
        acc = cost_c < cost
        log_vm = np.where(acc, cand_lvm, log_vm)
        log_km = np.where(acc, cand_lkm, log_km)
        kd = np.where(acc, cand_kd, kd)
        J_Y = np.where(acc[:, None, None, None], J_Sc * np.exp(-kappa * S_c)[:, :, :, None], J_Y)
        res = np.where(acc[:, None, None], res_c, res)
        cost = np.where(acc, cost_c, cost)
        lam = np.where(acc, np.maximum(lam * 0.3, 1e-6), np.minimum(lam * 5.0, 1e5))
    return np.exp(log_vm), np.exp(log_km), kd, cost


def fit_c3_cauchy(Y_obs, S0_init, t_grid, gamma, n_iter=24):
    """Cycle 35 (Q32): C3 matrix-hydrolysis term. Y_obs is C3 AFTER the hard-blank correction (a = 0).
    Model Y = S_MM+kd(t) * exp(-delta t), delta >= 0 with a half-Cauchy(0, gamma) prior (gamma ASSUMED,
    calibrated in sprint20_cycle35). sigma^2 is taken from the delta = 0 fit residual. Minimises
    0.5 * cost / sigma^2 + log(1 + (delta / gamma)^2) by damped Gauss-Newton (IRLS weight for the prior).
    Returns Vmax, Km, kd, delta (all (B,))."""
    B, W, M = Y_obs.shape
    zero = np.zeros((B, W))
    vm0, km0, kd0, c0 = batched_lm_readout_space(Y_obs, S0_init, t_grid, zero)
    s2 = np.maximum(c0 / max(W * M - 3, 1), 1e-10)
    w = np.maximum(S0_init, 0.5)[:, :, None]
    t3 = t_grid[None, None, :]
    lo = np.array([-7.0, -3.5, 0.0, 0.0])
    hi = np.array([3.0, 4.0, 0.30, 0.10])

    def ev(P):
        S, JS = mm_kd_eval_and_jac(S0_init, t_grid, P[:, 0], P[:, 1], P[:, 2])
        e = np.exp(-P[:, 3][:, None, None] * t3)
        Yp = S * e
        J = np.concatenate([JS * e[..., None], (-t3 * Yp)[..., None]], axis=-1) / w[..., None]
        r = (Yp - Y_obs) / w
        d = P[:, 3]
        F = 0.5 * np.sum(r * r, axis=(1, 2)) / s2 + np.log1p((d / gamma) ** 2)
        return r, J, F
    P = np.stack([np.log(vm0), np.log(km0), kd0, np.zeros(B)], axis=1)
    r, J, F = ev(P)
    lam = np.full(B, 1e-2)
    for _ in range(n_iter):
        d = P[:, 3]
        H = np.einsum('bwmi,bwmj->bij', J, J) / s2[:, None, None]
        g = np.einsum('bwmi,bwm->bi', J, r) / s2[:, None]
        H[:, 3, 3] += 2.0 / (gamma * gamma + d * d)
        g[:, 3] += 2.0 * d / (gamma * gamma + d * d)
        dg = np.maximum(np.diagonal(H, axis1=1, axis2=2), 1e-8)
        Hd = H + np.eye(4)[None] * (lam[:, None] * dg + 1e-9)[:, None, :]
        Pc = np.clip(P - np.linalg.solve(Hd, g[:, :, None])[:, :, 0], lo, hi)
        rc, Jc, Fc = ev(Pc)
        acc = Fc < F
        P = np.where(acc[:, None], Pc, P)
        r = np.where(acc[:, None, None], rc, r)
        J = np.where(acc[:, None, None, None], Jc, J)
        F = np.where(acc, Fc, F)
        lam = np.where(acc, np.maximum(lam * 0.3, 1e-7), np.minimum(lam * 5.0, 1e7))
    return np.exp(P[:, 0]), np.exp(P[:, 1]), P[:, 2], P[:, 3]


def c3_hard_blank(Yc, t_min):
    """C3 block after the hard-blank (C2) correction, and its S0_init, exactly as score_plates builds them."""
    c2, c3 = Yc[:, :, 2, :], 0.5 * (Yc[:, :, 3, :] + Yc[:, :, 4, :])
    frac = np.clip(c2[:, -3:, -1] / np.maximum(c2[:, -3:, 0], 1e-4), 0.1, 0.999)
    kn = np.median(-np.log(frac) / t_min[-1], axis=1)
    Hc = np.exp(kn[:, None, None] * t_min[None, None, :])
    return c3 * Hc, np.maximum(c3[:, :, 0], 1e-4)


def z_rao_orthogonal_robust(Y_obs, S0_init, t_grid, a_plate):
    B, W, M = Y_obs.shape
    w_scale = np.maximum(S0_init, 0.5)[:, :, None]
    kappa = (a_plate * LN10)[:, :, None]
    k_grid = np.linspace(-0.015, 0.18, 40)[None, None, None, :]
    S_fo = S0_init[:, :, None, None] * np.exp(-k_grid * t_grid[None, None, :, None])
    Y_fo = ife_forward(S_fo, a_plate[:, :, None, None])
    res_all = np.clip((Y_obs[:, :, :, None] - Y_fo) / w_scale[:, :, :, None], -0.06, 0.06)
    best_idx = np.argmin(np.sum(res_all * res_all, axis=(1, 2)), axis=1)
    k0_hat = k_grid[0, 0, 0, best_idx]
    S_n = S0_init[:, :, None] * np.exp(-k0_hat[:, None, None] * t_grid[None, None, :])
    Y_n = ife_forward(S_n, a_plate[:, :, None])
    dy_dS = np.exp(-kappa * S_n)
    r_raw = (Y_obs - Y_n) / w_scale
    med_r = np.median(r_raw, axis=(1, 2), keepdims=True)
    sig_rob = np.maximum(1.4826 * np.median(np.abs(r_raw - med_r), axis=(1, 2), keepdims=True), 0.008)
    r_clip = np.clip(r_raw, -3.0 * sig_rob, 3.0 * sig_rob)
    u0 = (-dy_dS * S_n * t_grid[None, None, :]) / w_scale
    u1 = (-dy_dS * (1.0 - np.exp(-0.25 * t_grid[None, None, :]))) / w_scale
    e0 = u0 / np.sqrt(np.sum(u0 * u0, axis=(1, 2), keepdims=True) + 1e-12)
    u1_perp = u1 - np.sum(u1 * e0, axis=(1, 2), keepdims=True) * e0
    e1 = u1_perp / np.sqrt(np.sum(u1_perp * u1_perp, axis=(1, 2), keepdims=True) + 1e-12)
    r_orth = r_clip - np.sum(r_clip * e0, axis=(1, 2), keepdims=True) * e0 \
                    - np.sum(r_clip * e1, axis=(1, 2), keepdims=True) * e1
    z_max = np.full(B, -1e9)
    sig_flat = sig_rob[:, 0, 0]
    for s_ref in (0.5, 2.0, 6.0):
        sat_w = S0_init[:, :, None] / (s_ref + S0_init[:, :, None])
        psi = (dy_dS * S_n * t_grid[None, None, :] * sat_w) / w_scale
        psi_orth = psi - np.sum(psi * e0, axis=(1, 2), keepdims=True) * e0 \
                       - np.sum(psi * e1, axis=(1, 2), keepdims=True) * e1
        psi_norm = np.sqrt(np.sum(psi_orth * psi_orth, axis=(1, 2)) + 1e-12)
        z_max = np.maximum(z_max, np.sum(r_orth * psi_orth, axis=(1, 2)) / (sig_flat * psi_norm))
    return z_max


# ============================================================================
# Parsing
# ============================================================================
class PlateParseError(Exception):
    def __init__(self, code, message):
        super().__init__(message)
        self.code = code
        self.message = message


WELL_RE = re.compile(r'^([A-Ha-h])\s*0*([1-9]|1[0-2])$')
NUM_RE = re.compile(r'^[+-]?(\d+([.]\d*)?|[.]\d+)([eE][+-]?\d+)?$')
DECOMMA_RE = re.compile(r'^[+-]?\d+,\d+$')
DELIMS = (("tab", "\t"), ("semicolon", ";"), ("comma", ","))


def _clean(cell):
    return cell.strip().strip('"').strip().replace('−', '-').replace('\xa0', '')


def well_id(cell):
    m = WELL_RE.match(_clean(cell))
    return f"{m.group(1).upper()}{int(m.group(2))}" if m else None


def well_index(wid):
    return ROWS.index(wid[0]), int(wid[1:]) - 1


def decode_bytes(raw):
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        return raw.decode('utf-16'), 'utf-16'
    if raw[:3] == b'\xef\xbb\xbf':
        return raw[3:].decode('utf-8', errors='replace'), 'utf-8-sig'
    try:
        return raw.decode('utf-8'), 'utf-8'
    except UnicodeDecodeError:
        return raw.decode('cp1252', errors='replace'), 'cp1252'


def split_rows(text, delim):
    lines = text.replace('\r\n', '\n').replace('\r', '\n').split('\n')
    return [list(r) for r in csv.reader(lines, delimiter=delim)]


def sniff(text):
    best = None
    for name, d in DELIMS:
        rows = split_rows(text, d)
        n = sum(1 for r in rows for c in r if well_id(c))
        if best is None or n > best[0]:
            best = (n, name, d, rows)
    if best[0] < 8:
        raise PlateParseError("NO_WELL_TABLE", "no table with well names (A1..H12) was found")
    return best[1], best[2], best[3]


def time_factor_s(label):
    """Seconds per unit implied by a header label; None if no unit is stated."""
    s = label.lower()
    if re.search(r'\bms\b|\[ms\]|\(ms\)', s):
        return 0.001
    if re.search(r'\[s\]|\(s\)|\bsec|\bseconds?\b|\bs\b', s):
        return 1.0
    if re.search(r'\bmin|\[m\]|\(m\)', s):
        return 60.0
    if re.search(r'\[h\]|\(h\)|\bhours?\b|\bhr\b', s):
        return 3600.0
    return None


class _Ctx:
    def __init__(self, decimal_comma):
        self.decimal_comma = decimal_comma
        self.flag_tokens = {}
        self.warnings = []

    def num(self, cell):
        """float, or NaN for a non-numeric token (counted), or None for an empty cell."""
        s = _clean(cell)
        if s == "":
            return None
        if self.decimal_comma and DECOMMA_RE.match(s):
            s = s.replace(',', '.')
        if NUM_RE.match(s):
            return float(s)
        tok = s.upper()[:16]
        self.flag_tokens[tok] = self.flag_tokens.get(tok, 0) + 1
        return float('nan')

    def time_s(self, cell, factor):
        s = _clean(cell)
        if s == "":
            return None
        if ':' in s:
            parts = s.split(':')
            try:
                vals = [float(p.replace(',', '.')) for p in parts]
            except ValueError:
                return None
            if len(vals) == 3:
                return vals[0] * 3600.0 + vals[1] * 60.0 + vals[2]
            if len(vals) == 2:
                return vals[0] * 60.0 + vals[1]
            return None
        if self.decimal_comma and DECOMMA_RE.match(s):
            s = s.replace(',', '.')
        if not NUM_RE.match(s):
            return None
        return float(s) * factor


def _numeric_time_factor(ctx, label, samples):
    f = time_factor_s(label)
    if f is not None:
        return f
    vals = [ctx.time_s(c, 1.0) for c in samples]
    vals = [v for v in vals if v is not None]
    if any(':' in _clean(c) for c in samples):
        return 1.0
    if vals and max(vals) > 100.0:
        ctx.warnings.append(f"time header '{label}' has no unit; values up to {max(vals):g} read as seconds (ASSUMED)")
        return 1.0
    ctx.warnings.append(f"time header '{label}' has no unit; read as minutes (ASSUMED)")
    return 60.0


def _parse_wide(rows, hdr, ctx):
    head = rows[hdr]
    wcols = [(j, well_id(c)) for j, c in enumerate(head) if well_id(c)]
    tcol = next((j for j, c in enumerate(head) if 'time' in _clean(c).lower()), 0)
    data = []
    k = hdr + 1
    while k < len(rows) and any(_clean(c) for c in rows[k]):
        data.append(rows[k])
        k += 1
    if not data:
        raise PlateParseError("NO_READS", "well header found but no data rows follow it")
    factor = _numeric_time_factor(ctx, head[tcol] if tcol < len(head) else "", [r[tcol] for r in data if tcol < len(r)])
    times, table = [], []
    for r in data:
        t = ctx.time_s(r[tcol], factor) if tcol < len(r) else None
        if t is None:
            break
        times.append(t)
        table.append(r)
    extra = sum(1 for r in rows[k:] if sum(1 for c in r if well_id(c)) >= 8)
    if extra:
        ctx.warnings.append(f"{extra} further well table(s) after the first were ignored")
    wells = {}
    for j, wid in wcols:
        if wid in wells:
            ctx.warnings.append(f"duplicate column for well {wid}; first kept")
            continue
        wells[wid] = [ctx.num(r[j]) if j < len(r) else None for r in table]
    return np.array(times), wells


def _parse_transposed(rows, ctx):
    wrows = [r for r in rows if r and well_id(r[0])]
    first = next(i for i, r in enumerate(rows) if r and well_id(r[0]))
    trow = None
    for r in rows[:first][::-1]:
        if r and 'time' in _clean(r[0]).lower():
            trow = r
            break
    if trow is None:
        raise PlateParseError("NO_TIME", "per-well rows found but no 'Time' row above them")
    factor = _numeric_time_factor(ctx, trow[0], trow[1:])
    idx, times = [], []
    for j in range(1, len(trow)):
        t = ctx.time_s(trow[j], factor)
        if t is not None:
            idx.append(j)
            times.append(t)
    wells = {}
    for r in wrows:
        wid = well_id(r[0])
        if wid in wells:
            ctx.warnings.append(f"duplicate row for well {wid}; first kept")
            continue
        wells[wid] = [ctx.num(r[j]) if j < len(r) else None for j in idx]
    return np.array(times), wells


def _parse_long(rows, hdr, wc, ctx):
    head = [_clean(c).lower() for c in rows[hdr]]
    tc = next((j for j, c in enumerate(head) if 'time' in c), None)
    if tc is None:
        raise PlateParseError("NO_TIME", "long-format header has no Time column")
    skip = {wc, tc}
    cand = [j for j, c in enumerate(head) if j not in skip and not re.search(r'temp|content|sample|type|group|name', c)]
    pref = [j for j in cand if re.search(r'abs|od|value|signal|read|raw|\d{3}', head[j])]
    if not (pref or cand):
        raise PlateParseError("NO_VALUE", "long-format header has no value column")
    vc = (pref or cand)[0]
    body = [r for r in rows[hdr + 1:] if len(r) > max(wc, tc, vc) and well_id(r[wc])]
    factor = _numeric_time_factor(ctx, rows[hdr][tc], [r[tc] for r in body])
    trip = []
    for r in body:
        t = ctx.time_s(r[tc], factor)
        if t is None:
            ctx.warnings.append("a row with an unreadable time was skipped")
            continue
        trip.append((well_id(r[wc]), t, r[vc]))
    times = np.unique(np.array([t for _, t, _ in trip]))
    pos = {t: i for i, t in enumerate(times.tolist())}
    wells, seen = {}, set()
    for wid, t, v in trip:
        if (wid, t) in seen:
            ctx.warnings.append(f"duplicate read {wid} at {t:g} s; first kept")
            continue
        seen.add((wid, t))
        wells.setdefault(wid, [None] * len(times))[pos[t]] = ctx.num(v)
    return times, wells


def _parse_softmax(text, wl_index):
    """SoftMax Pro plate-format kinetic export (Cycle 40). Fails loud on anything it does not recognise."""
    lines = text.splitlines()
    hdr = next((ln for ln in lines if ln.startswith("Plate:")), None)
    if hdr is None:
        raise PlateParseError("SOFTMAX_NO_PLATE", "##BLOCKS header found but no 'Plate:' section")
    h = hdr.split("\t")
    if len(h) < 16 or h[4] != "Kinetic" or h[3] != "PlateFormat":
        raise PlateParseError("SOFTMAX_NOT_KINETIC", "only PlateFormat Kinetic SoftMax sections are supported")
    try:
        n_reads, n_wl = int(h[8]), int(h[14])
    except ValueError:
        raise PlateParseError("SOFTMAX_BAD_HEADER", "read count or wavelength count is not an integer")
    wls = h[15].split()
    if n_wl < 1 or len(wls) != n_wl:
        raise PlateParseError("SOFTMAX_BAD_HEADER", f"declared {n_wl} wavelengths but listed {wls}")
    if not 0 <= wl_index < n_wl:
        raise PlateParseError("SOFTMAX_WAVELENGTH", f"wavelength index {wl_index} not in 0..{n_wl - 1} ({wls})")
    start = lines.index(hdr) + 2
    times, blocks, cur = [], [], None
    for ln in lines[start:]:
        if ln.startswith("~End"):
            break
        f = ln.split("\t")
        if f[0] and f[0][0].isdigit() and f[0].count(":") == 2:
            try:
                hh, mm, ss = (int(x) for x in f[0].split(":"))
            except ValueError:
                raise PlateParseError("SOFTMAX_BAD_TIME", f"cannot read time {f[0]!r}")
            times.append(float(3600 * hh + 60 * mm + ss))
            cur = []
            blocks.append(cur)
        elif ln == "":
            continue
        elif cur is not None and len(cur) == 8 and len(blocks) == n_reads:
            break  # the declared reads are complete; what follows is a reduced-data section
        elif f[0] != "" or cur is None or len(cur) >= 8:
            raise PlateParseError("SOFTMAX_BAD_BLOCK", f"unexpected line in kinetic block: {ln[:40]!r}")
        if len(f) < 2 + 13 * n_wl - 1:
            raise PlateParseError("SOFTMAX_BAD_BLOCK", f"a plate row has {len(f)} fields, need {2 + 13 * n_wl - 1}")
        c0 = 2 + 13 * wl_index
        cur.append([_clean(x) for x in f[c0:c0 + 12]])
    if len(blocks) != n_reads:
        raise PlateParseError("SOFTMAX_READ_COUNT", f"header declares {n_reads} reads, found {len(blocks)}")
    if any(len(b) != 8 for b in blocks):
        raise PlateParseError("SOFTMAX_BAD_BLOCK", "a read block does not have 8 plate rows")
    wells = {}
    for r in range(8):
        for c in range(12):
            vals = []
            for b in blocks:
                x = b[r][c]
                if x == "":
                    vals.append(None)
                elif NUM_RE.match(x):
                    vals.append(float(x))
                else:
                    raise PlateParseError("SOFTMAX_BAD_VALUE", f"non-numeric cell {x!r} in {ROWS[r]}{c + 1}")
            if any(v is not None for v in vals):
                wells[f"{ROWS[r]}{c + 1}"] = vals
    if not wells:
        raise PlateParseError("SOFTMAX_EMPTY", "no numeric wells in the selected wavelength")
    return np.array(times), wells, wls[wl_index]


def load_layout(raw):
    """--layout JSON: {"row_S0_uM": [8 floats] (optional), "columns": {"1": [8 source wells], ...}}."""
    try:
        spec = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError):
        raise PlateParseError("LAYOUT_BAD_JSON", "the layout file is not UTF-8 JSON")
    if not isinstance(spec, dict) or not isinstance(spec.get("columns"), dict) or not spec["columns"]:
        raise PlateParseError("LAYOUT_BAD_SPEC", "the layout needs a non-empty 'columns' object")
    lad = spec.get("row_S0_uM", S0_UM.tolist())
    if (not isinstance(lad, list) or len(lad) != 8 or not all(isinstance(x, (int, float)) and x > 0 for x in lad)
            or any(b <= a for a, b in zip(lad, lad[1:]))):
        raise PlateParseError("LAYOUT_BAD_LADDER", "row_S0_uM must be 8 increasing positive numbers")
    cols = {}
    for k, v in spec["columns"].items():
        if not (isinstance(k, str) and k.isdigit() and 1 <= int(k) <= 12):
            raise PlateParseError("LAYOUT_BAD_SPEC", f"contract column {k!r} is not 1..12")
        if not isinstance(v, list) or len(v) != 8 or not all(isinstance(w, str) and well_id(w) for w in v):
            raise PlateParseError("LAYOUT_BAD_SPEC", f"column {k} needs 8 source well names")
        cols[int(k) - 1] = [well_id(w) for w in v]
    return {"row_S0_uM": [float(x) for x in lad], "columns": cols}


def apply_layout(parsed, lay):
    M = parsed["A"].shape[2]
    A = np.full((8, 12, M), np.nan)
    present = np.zeros((8, 12), dtype=bool)
    for j, src in lay["columns"].items():
        for r, w in enumerate(src):
            rr, cc = well_index(w)
            A[r, j] = parsed["A"][rr, cc]
            present[r, j] = parsed["present"][rr, cc]
    out = dict(parsed)
    out["A"], out["present"] = A, present
    return out


def parse_bytes(raw, wavelength_index=0):
    """Parse a plate export. Returns dict with times_s (M,), A (8,12,M) with NaN gaps, meta."""
    if not raw or not raw.strip():
        raise PlateParseError("EMPTY", "the file is empty")
    text, enc = decode_bytes(raw)
    if text.lstrip("\ufeff \r\n").startswith("##BLOCKS="):
        times, wells, wl = _parse_softmax(text, wavelength_index)
        return _assemble(times, wells, _Ctx(False), enc, "tab", False, "softmax_kinetic", {"wavelength": wl})
    dname, delim, rows = sniff(text)
    cells = [_clean(c) for r in rows for c in r]
    decimal_comma = delim != ',' and any(DECOMMA_RE.match(c) for c in cells)
    ctx = _Ctx(decimal_comma)
    layout = None
    for i, r in enumerate(rows):
        low = [_clean(c).lower() for c in r]
        wc = next((j for j, c in enumerate(low) if c.startswith('well')), None)
        if wc is not None and any('time' in c for c in low) and sum(1 for c in r if well_id(c)) == 0:
            layout = "long"
            times, wells = _parse_long(rows, i, wc, ctx)
            break
        if sum(1 for c in r if well_id(c)) >= 8:
            layout = "wide"
            times, wells = _parse_wide(rows, i, ctx)
            break
    if layout is None:
        if sum(1 for r in rows if r and well_id(r[0])) >= 8:
            layout = "transposed"
            times, wells = _parse_transposed(rows, ctx)
        else:
            raise PlateParseError("NO_WELL_TABLE", "well names found, but not in a wide, transposed or long table")
    return _assemble(times, wells, ctx, enc, dname, decimal_comma, layout, {})


def _assemble(times, wells, ctx, enc, dname, decimal_comma, layout, extra):
    M = len(times)
    if M < 4:
        raise PlateParseError("TOO_FEW_READS", f"only {M} time points were read (need >= 4)")
    if np.any(np.diff(times) <= 0):
        raise PlateParseError("TIME_NOT_MONOTONIC", "read times are not strictly increasing")
    A = np.full((8, 12, M), np.nan)
    present = np.zeros((8, 12), dtype=bool)
    n_empty = 0
    for wid, vals in wells.items():
        r, c = well_index(wid)
        present[r, c] = True
        for m, v in enumerate(vals):
            if v is None:
                n_empty += 1
            else:
                A[r, c, m] = v
    if times[0] != 0.0:
        ctx.warnings.append(f"first read at {times[0]:g} s taken as t = 0 (ASSUMED)")
    meta = {"encoding": enc, "delimiter": dname, "decimal_comma": bool(decimal_comma), "layout": layout,
            "n_reads": int(M), "flag_tokens": dict(sorted(ctx.flag_tokens.items())), "empty_cells": int(n_empty),
            "warnings": ctx.warnings, "t_offset_s": float(times[0]), **extra}
    return {"times_s": times - times[0], "A": A, "present": present, "meta": meta}


def fill_gaps(A, present, max_bad_frac=0.25):
    """Linear interpolation in time inside each well. Returns filled array, #filled cells, bad-well mask."""
    Y = A.copy()
    M = A.shape[2]
    bad = ~present.copy()
    n_fill = 0
    for r in range(8):
        for c in range(12):
            if not present[r, c]:
                Y[r, c, :] = 0.0
                continue
            ok = np.isfinite(A[r, c])
            nbad = int(M - ok.sum())
            if nbad == 0:
                continue
            if nbad > max_bad_frac * M or ok.sum() < 2:
                bad[r, c] = True
                Y[r, c, :] = np.where(ok, A[r, c], 0.0)
                continue
            idx = np.arange(M)
            Y[r, c] = np.interp(idx, idx[ok], A[r, c, ok])
            n_fill += nbad
    return Y, n_fill, bad


# ============================================================================
# Contract scoring (batched over B plates; B = 1 for a real file)
# ============================================================================
def score_plates(Yc, t_min, s0_norm=None):
    """Yc (B, 8, 12, M) normalised readout; t_min (M,). Returns dict of (B,) arrays."""
    s0_norm = S0_NORM if s0_norm is None else s0_norm
    B = Yc.shape[0]
    col = lambda j: Yc[:, :, j, :]
    c1a, c1b = col(0), col(1)
    c1 = 0.5 * (c1a + c1b)
    c2, c3, c4 = col(2), 0.5 * (col(3) + col(4)), 0.5 * (col(5) + col(6))
    c6, c7a, c7b = col(9), col(10), col(11)
    frac = np.clip(c2[:, -3:, -1] / np.maximum(c2[:, -3:, 0], 1e-4), 0.1, 0.999)
    kn = np.median(-np.log(frac) / t_min[-1], axis=1)
    Hc = np.exp(kn[:, None, None] * t_min[None, None, :])
    pm = np.ones(8)
    s0 = lambda Y: np.maximum(Y[:, :, 0], 1e-4)
    S_c4, a4, trig4 = conditional_ife_self_cal(c4, s0_norm, pm)
    c4h = ife_forward(ife_inverse(c4, a4[:, None, None]) * Hc, a4[:, None, None])
    zero = np.zeros((B, 8))
    groups = [(c1 * Hc, s0(c1), zero), (c1a * Hc, s0(c1a), zero), (c1b * Hc, s0(c1b), zero),
              (c3 * Hc, s0(c3), zero), (c6 * Hc, s0(c6), zero), (c4, s0(c4), zero),
              (c4h, s0(S_c4), a4[:, None] * pm[None, :])]
    vm, km, kd, _ = batched_lm_readout_space(np.concatenate([g[0] for g in groups]),
                                             np.concatenate([g[1] for g in groups]), t_min,
                                             np.concatenate([g[2] for g in groups]))
    sl = lambda x, i: x[i * B:(i + 1) * B]
    S7b, a7b, _ = conditional_ife_self_cal(c7b, s0_norm, pm)
    z = z_rao_orthogonal_robust(np.concatenate([c1, c7a, c7b]),
                                np.concatenate([s0(c1), s0(c7a), s0(S7b)]), t_min,
                                np.concatenate([zero, zero, a7b[:, None] * pm[None, :]]))
    out = {"kn_blank": kn, "vm1": sl(vm, 0), "km1": sl(km, 0), "kd1": sl(kd, 0),
           "vm1a": sl(vm, 1), "km1a": sl(km, 1), "vm1b": sl(vm, 2), "km1b": sl(km, 2),
           "km3": sl(km, 3), "kd6": sl(kd, 4), "km4_uncal": sl(km, 5), "km4": sl(km, 6),
           "a4": a4, "ife_trig4": trig4.astype(float), "a7b": a7b,
           "z1": sl(z, 0), "z7a": sl(z, 1), "z7b": sl(z, 2)}
    return out


def verdicts(s, p1_rule=None, s0_norm=None):
    """Return dict test -> int array (1 PASS, 0 KILL, -1 INCONCLUSIVE) plus the derived values."""
    s0_norm = S0_NORM if s0_norm is None else s0_norm
    thr = P1_REP_THR[p1_rule or P1_RULE]
    km_lo, km_hi = s0_norm[0] / 2.0, 2.0 * s0_norm[-1]
    lr = lambda a, b: np.abs(np.log(np.maximum(a, 1e-12) / np.maximum(b, 1e-12)))
    eff1a, eff1b = s["vm1a"] / s["km1a"], s["vm1b"] / s["km1b"]
    active = s["z1"] >= Z_THR
    in_range = (s["km1"] >= km_lo) & (s["km1"] <= km_hi) & (np.log(s["vm1"]) < 3.0 - 1e-6)
    rep_km, rep_eff = lr(s["km1a"], s["km1b"]), lr(eff1a, eff1b)
    p1 = active & in_range & (rep_km <= thr) & (rep_eff <= thr)
    p2 = np.where(p1, (lr(s["km4"], s["km1"]) <= math.log(1.20)).astype(int), -1)
    p3 = np.where(p1, ((s["kd6"] - s["kd1"]) > 0.03).astype(int), -1)
    p4 = np.where(p1, (lr(s["km3"], s["km1"]) <= math.log(1.15)).astype(int), -1)
    p5 = (active & (s["z7a"] < Z_THR) & (s["z7b"] < Z_THR)).astype(int)
    derived = {"active_c1": active, "km_in_range": in_range, "rep_km_lnratio": rep_km, "rep_eff_lnratio": rep_eff,
               "p2_lnratio": lr(s["km4"], s["km1"]), "p3_dkd": s["kd6"] - s["kd1"],
               "p4_bias": s["km3"] / s["km1"] - 1.0}
    return {"P1": p1.astype(int), "P2": p2, "P3": p3, "P4": p4, "P5": p5}, derived


# ============================================================================
# Report
# ============================================================================
VNAME = {1: "PASS", 0: "KILL", -1: "INCONCLUSIVE"}


def _r(x):
    if isinstance(x, dict):
        return {k: _r(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_r(v) for v in x]
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        return None if not math.isfinite(x) else float(f"{x:.6g}")
    return x


def to_normalised(A, times_s, scale, offset, mode, s0_norm=None):
    s0_norm = S0_NORM if s0_norm is None else s0_norm
    if mode == "substrate":
        return (A - offset) / scale
    return s0_norm[:, None, None] - (A - A[:, :, :1]) / scale


def build_report(parsed, name, sha, scale=1.0, offset=0.0, mode="substrate", evidence="MEASURED", p1_rule=None,
                 layout=None):
    p1_rule = p1_rule or P1_RULE
    s0_um = np.array(layout["row_S0_uM"]) if layout else S0_UM
    s0n = s0_um / NORM_UM
    if layout:
        parsed = apply_layout(parsed, layout)
    p1_thr = P1_REP_THR[p1_rule]
    Afill, n_fill, bad = fill_gaps(parsed["A"], parsed["present"])
    t_min = parsed["times_s"] / 60.0
    Y = to_normalised(Afill, parsed["times_s"], scale, offset, mode, s0n)
    s = score_plates(Y[None], t_min, s0n)
    v, d = verdicts(s, p1_rule, s0n)
    bad_wells = [f"{ROWS[r]}{c + 1}" for r in range(8) for c in range(12) if bad[r, c]]
    missing = [f"{ROWS[r]}{c + 1}" for r in range(8) for c in range(12) if not parsed["present"][r, c]]
    block_bad = {b: any(bad[:, j].any() for j in cols) for b, cols in BLOCKS.items()}
    one = lambda a: a[0]
    tests = {}
    rules = {
        "P1": (f"C1 active (z >= 4.60), Km in [S0min/2, 2*S0max], Vmax off bound, replicate Km and Vmax/Km "
               f"|ln ratio| <= {p1_thr:.4g} (rule {p1_rule})"),
        "P2": "IFE-calibrated C4 Km within 20% of C1 Km",
        "P3": "kd(C6) - kd(C1) > 0.03 /min",
        "P4": "C3 Km within 15% of C1 Km after blank-only hydrolysis correction (KILL = hard blank falsified)",
        "P5": "z(C1) >= 4.60 and z(C7a) < 4.60 and z(C7b) < 4.60",
    }
    values = {
        "P1": {"z_c1": one(s["z1"]), "Km_uM": one(s["km1"]) * NORM_UM, "Vmax_uM_per_min": one(s["vm1"]) * NORM_UM,
               "kd_per_min": one(s["kd1"]), "Km_rep_uM": [one(s["km1a"]) * NORM_UM, one(s["km1b"]) * NORM_UM],
               "rep_Km_lnratio": one(d["rep_km_lnratio"]), "rep_Eff_lnratio": one(d["rep_eff_lnratio"]),
               "kn_blank_per_min": one(s["kn_blank"])},
        "P2": {"Km_c4_uncal_uM": one(s["km4_uncal"]) * NORM_UM, "Km_c4_cal_uM": one(s["km4"]) * NORM_UM,
               "a_hat": one(s["a4"]), "ife_triggered": bool(one(s["ife_trig4"])), "lnratio_vs_C1": one(d["p2_lnratio"])},
        "P3": {"kd_c1": one(s["kd1"]), "kd_c6": one(s["kd6"]), "dkd": one(d["p3_dkd"])},
        "P4": {"Km_c3_uM": one(s["km3"]) * NORM_UM, "bias_vs_C1": one(d["p4_bias"])},
        "P5": {"z_c1": one(s["z1"]), "z_c7a": one(s["z7a"]), "z_c7b": one(s["z7b"]), "a_hat_c7b": one(s["a7b"])},
    }
    for p in ("P1", "P2", "P3", "P4", "P5"):
        code = int(v[p][0])
        reason = None
        badb = [b for b in TEST_BLOCKS[p] if block_bad[b]]
        if badb:
            code, reason = -1, "unreadable or missing wells in " + ", ".join(badb)
        elif code == -1:
            reason = "P1 did not pass, so there is no valid C1 baseline"
        tests[p] = {"verdict": VNAME[code], "rule": rules[p], "values": values[p], "evidence": evidence,
                    "reason": reason}
    if block_bad["C1"]:
        for p in ("P2", "P3", "P4"):
            tests[p]["verdict"], tests[p]["reason"] = "INCONCLUSIVE", "unreadable or missing wells in C1"
    rep = {
        "schema": SCHEMA,
        "status": "OK",
        "input": {"file": name, "sha256": sha, **{k: parsed["meta"][k] for k in ("encoding", "delimiter", "decimal_comma", "layout")},
                  **({"wavelength_nm": parsed["meta"]["wavelength"]} if "wavelength" in parsed["meta"] else {})},
        "assumptions": {"tag": "ASSUMED", "mode": mode, "scale_A_per_norm": scale, "offset_A": offset,
                        "norm_unit_uM": NORM_UM, "row_S0_uM": s0_um.tolist(),
                        "columns": {b: [c + 1 for c in cols] for b, cols in BLOCKS.items()},
                        "t0": "first read", "t_offset_s": parsed["meta"]["t_offset_s"], "z_threshold": Z_THR,
                        "p1_rule": p1_rule, "p1_rep_lnratio_max": p1_thr, "C5": "parsed, not scored"},
        "qc": {"n_reads": parsed["meta"]["n_reads"], "t_end_min": float(t_min[-1]),
               "t_offset_s": parsed["meta"]["t_offset_s"],
               "non_numeric_tokens": parsed["meta"]["flag_tokens"], "empty_cells": parsed["meta"]["empty_cells"],
               "interpolated_cells": int(n_fill), "missing_wells": missing,
               "bad_wells": [w for w in bad_wells if w not in missing], "warnings": parsed["meta"]["warnings"]},
        "tests": tests,
        **({"design_match": layout["design_match"]} if layout else {}),
        "summary": {k: [p for p in tests if tests[p]["verdict"] == k] for k in ("PASS", "KILL", "INCONCLUSIVE")},
    }
    return _r(rep)


def report_markdown(rep):
    L = []
    if rep.get("status") != "OK":
        return f"# Plate report: ERROR\n\n`{rep['error']['code']}`: {rep['error']['message']}\n"
    i, q = rep["input"], rep["qc"]
    L.append(f"# Plate report: {i['file']}\n")
    L.append(f"SHA-256 `{i['sha256'][:16]}…` · layout **{i['layout']}** · delimiter {i['delimiter']} · "
             f"encoding {i['encoding']} · decimal comma {i['decimal_comma']}\n")
    L.append("| Test | Verdict | Evidence | Key values | Rule |")
    L.append("|---|---|---|---|---|")
    for p, tst in rep["tests"].items():
        vals = ", ".join(f"{k}={v}" for k, v in tst["values"].items())
        why = f" ({tst['reason']})" if tst["reason"] else ""
        L.append(f"| {p} | **{tst['verdict']}**{why} | {tst['evidence']} | {vals} | {tst['rule']} |")
    L.append("")
    L.append(f"**QC:** {q['n_reads']} reads to {q['t_end_min']:.4g} min; non-numeric cells {q['non_numeric_tokens'] or 'none'}; "
             f"empty cells {q['empty_cells']}; interpolated {q['interpolated_cells']}; "
             f"missing wells {q['missing_wells'] or 'none'}; unusable wells {q['bad_wells'] or 'none'}.")
    for w in q["warnings"]:
        L.append(f"- warning: {w}")
    a = rep["assumptions"]
    L.append("")
    L.append(f"**ASSUMED:** mode {a['mode']}, scale {a['scale_A_per_norm']} A per {a['norm_unit_uM']:g} uM, offset {a['offset_A']} A, "
             f"rows S0 = {a['row_S0_uM']} uM, t0 = first read (offset {a['t_offset_s']:g} s subtracted), "
             f"z threshold {a['z_threshold']}, P1 rule {a['p1_rule']} (replicate |ln ratio| <= {a['p1_rep_lnratio_max']}). "
             f"C5 is parsed but not scored.")
    return "\n".join(L) + "\n"


def process_bytes(raw, name, wavelength_index=0, layout_raw=None, layout_name=None, **kw):
    sha = hashlib.sha256(raw).hexdigest()
    try:
        parsed = parse_bytes(raw, wavelength_index)
        if layout_raw is not None:
            lay = load_layout(layout_raw)
            used = sorted(lay["columns"])
            absent = [b for b, cols in BLOCKS.items() if not all(j in lay["columns"] for j in cols)]
            same = bool(np.allclose(lay["row_S0_uM"], S0_UM))
            lay["design_match"] = {
                "layout_file": layout_name, "layout_sha256": hashlib.sha256(layout_raw).hexdigest(),
                "status": "CONTRACT_LADDER" if same else "CUSTOM_LADDER",
                "note": None if same else "thresholds were calibrated on the contract ladder only (ASSUMED transfer)",
                "row_S0_uM": lay["row_S0_uM"], "contract_row_S0_uM": S0_UM.tolist(),
                "mapped_columns": [j + 1 for j in used], "absent_blocks": absent}
            kw["layout"] = lay
    except PlateParseError as e:
        return {"schema": SCHEMA, "status": "ERROR", "input": {"file": name, "sha256": sha},
                "error": {"code": e.code, "message": e.message}}, None
    return build_report(parsed, name, sha, **kw), parsed


def dumps(rep):
    return json.dumps(rep, indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def main(argv):
    import argparse
    ap = argparse.ArgumentParser(description="96-well kinetic export -> Cycle 29 P1-P5 report")
    ap.add_argument("path")
    ap.add_argument("--scale", type=float, default=1.0, help="absorbance per 50 uM substrate (ASSUMED)")
    ap.add_argument("--offset", type=float, default=0.0, help="background absorbance (ASSUMED)")
    ap.add_argument("--mode", choices=("substrate", "product"), default="substrate")
    ap.add_argument("--evidence", choices=("MEASURED", "MODELED"), default="MEASURED",
                    help="MEASURED for a real plate, MODELED for a simulated one")
    ap.add_argument("--p1-rule", choices=sorted(P1_REP_THR), default=None,
                    help="P1 replicate threshold rule (default: the module's P1_RULE)")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--wavelength-index", type=int, default=0, help="SoftMax exports: which wavelength (0 = first)")
    ap.add_argument("--layout", default=None, help="JSON map of contract columns to source wells (Cycle 40)")
    a = ap.parse_args(argv)
    lay_raw = None
    try:
        with open(a.path, "rb") as fh:
            raw = fh.read()
        if a.layout:
            with open(a.layout, "rb") as fh:
                lay_raw = fh.read()
    except OSError as e:
        print(json.dumps({"status": "ERROR", "error": {"code": "IO", "message": str(e)}}))
        return 2
    name = os.path.basename(a.path)
    rep, _ = process_bytes(raw, name, wavelength_index=a.wavelength_index, layout_raw=lay_raw,
                           layout_name=os.path.basename(a.layout) if a.layout else None,
                           scale=a.scale, offset=a.offset, mode=a.mode, evidence=a.evidence, p1_rule=a.p1_rule)
    out = a.out_dir or os.path.dirname(os.path.abspath(a.path))
    os.makedirs(out, exist_ok=True)
    stem = os.path.join(out, os.path.splitext(name)[0] + ".report")
    with open(stem + ".json", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(dumps(rep))
    with open(stem + ".md", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(report_markdown(rep))
    if rep["status"] != "OK":
        print(f"ERROR {rep['error']['code']}: {rep['error']['message']}")
        return 2
    print(" ".join(f"{p}={t['verdict']}" for p, t in rep["tests"].items()))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
