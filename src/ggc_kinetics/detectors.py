"""Frozen detection statistics from the GGC-inf archive.

z_rao_orth_robust : saturation-curvature score for a progress-curve plate (Cycles 19-36). Needs several S0 levels
                    on one plate: with one well per plate the curvature direction is projected out and z = 0.
z_dsf_asym_boltz  : sigmoid-transition score for a thermal-shift (DSF) melt curve (Cycles 27-28, K60).
Thresholds are ASSUMED values calibrated on MODELED nulls; they are not yet confirmed on real plates.
"""
import math

import numpy as np

from .solvers import ife_forward

LN10 = math.log(10.0)
Z_RAO_BASELINE = 4.60          # Cycle 19 clean-trained threshold, no noise inflation
Z_RAO_K59_TAU = 0.056          # Cycle 36 (K59): noise-floor inflation tau
Z_RAO_K59_THR = 1.768          # Cycle 36 (K59): threshold for the tau-inflated statistic
Z_DSF_THR = 5.72               # Cycle 28 (K60): mixed Gaussian + t4 calibrated threshold

__all__ = ["z_rao_orth_robust", "z_rao_call", "z_dsf_asym_boltz", "Z_RAO_BASELINE", "Z_RAO_K59_TAU",
           "Z_RAO_K59_THR", "Z_DSF_THR"]


def z_rao_orth_robust(Y, S0, t_grid, a_plate, tau=None):
    """Z_Rao_Orth_Robust per plate. Y (B, W, M) normalized readout, S0 (B, W) normalized substrate, t_grid (M,)
    in minutes, a_plate (B, W) IFE absorbance per norm unit (0 = no inner filter).

    tau=None returns the base statistic (compare with 4.60). tau=0.056 returns the K59 noise-inflated statistic
    num / sqrt(sigma^2 + tau^2) (compare with 1.768). Same math as plate_parser.z_rao_orthogonal_robust and
    sprint21_cycle36_k59_m1anchor.z_rao.
    """
    Y = np.asarray(Y, float)
    S0 = np.asarray(S0, float)
    t_grid = np.asarray(t_grid, float)
    a_plate = np.asarray(a_plate, float)
    B = Y.shape[0]
    w = np.maximum(S0, 0.5)[:, :, None]
    kap = (a_plate * LN10)[:, :, None]
    kg = np.linspace(-0.015, 0.18, 40)[None, None, None, :]
    Yfo = ife_forward(S0[:, :, None, None] * np.exp(-kg * t_grid[None, None, :, None]), a_plate[:, :, None, None])
    ra = np.clip((Y[:, :, :, None] - Yfo) / w[:, :, :, None], -0.06, 0.06)
    k0 = kg[0, 0, 0, np.argmin(np.sum(ra * ra, axis=(1, 2)), axis=1)]
    Sn = S0[:, :, None] * np.exp(-k0[:, None, None] * t_grid[None, None, :])
    Yn = ife_forward(Sn, a_plate[:, :, None])
    dy = np.exp(-kap * Sn)
    r = (Y - Yn) / w
    med = np.median(r, axis=(1, 2), keepdims=True)
    sig = np.maximum(1.4826 * np.median(np.abs(r - med), axis=(1, 2), keepdims=True), 0.008)
    rc = np.clip(r, -3.0 * sig, 3.0 * sig)

    def nrm(v):
        return v / np.sqrt(np.sum(v * v, axis=(1, 2), keepdims=True) + 1e-12)

    def pr(v, e):
        return v - np.sum(v * e, axis=(1, 2), keepdims=True) * e
    e0 = nrm((-dy * Sn * t_grid) / w)
    e1 = nrm(pr((-dy * (1.0 - np.exp(-0.25 * t_grid))) / w, e0))
    ro = pr(pr(rc, e0), e1)
    num = np.full(B, -1e9)
    for s_ref in (0.5, 2.0, 6.0):
        psi = (dy * Sn * t_grid * (S0[:, :, None] / (s_ref + S0[:, :, None]))) / w
        po = pr(pr(psi, e0), e1)
        num = np.maximum(num, np.sum(ro * po, axis=(1, 2)) / np.sqrt(np.sum(po * po, axis=(1, 2)) + 1e-12))
    s = sig[:, 0, 0]
    return num / s if tau is None else num / np.sqrt(s * s + tau * tau)


def z_rao_call(Y, S0, t_grid, a_plate, rule="k59"):
    """Boolean call per plate: rule 'k59' (tau 0.056, thr 1.768) or 'baseline' (thr 4.60)."""
    if rule == "k59":
        return z_rao_orth_robust(Y, S0, t_grid, a_plate, tau=Z_RAO_K59_TAU) >= Z_RAO_K59_THR
    if rule == "baseline":
        return z_rao_orth_robust(Y, S0, t_grid, a_plate) >= Z_RAO_BASELINE
    raise ValueError("rule must be 'k59' or 'baseline'")


def z_dsf_asym_boltz(T, y, Tm=55.0, width=2.5, clip=3.5):
    """Robust score for a Boltzmann melt transition at (Tm, width) on top of a quadratic baseline.
    The template sigmoid is QR-orthogonalized against [1, T, T^2]; residuals are clipped at 3.5 MAD.
    Call a transition at z >= 5.72 (K60). This is a fixed-template test: scan Tm yourself if it is unknown,
    and note that a scan needs its own multiplicity-corrected threshold (not calibrated here)."""
    T = np.asarray(T, float)
    y = np.asarray(y, float)
    sig = 1.0 / (1.0 + np.exp(-(T - Tm) / width))
    X = np.column_stack([np.ones_like(T), T, T ** 2])
    Q, _ = np.linalg.qr(X)
    so = sig - Q @ (Q.T @ sig)
    so /= np.linalg.norm(so)
    r = y - Q @ (Q.T @ y)
    mad = np.median(np.abs(r - np.median(r))) * 1.4826 + 1e-9
    return float(np.dot(so, np.clip(r, -clip * mad, clip * mad)) / mad)
