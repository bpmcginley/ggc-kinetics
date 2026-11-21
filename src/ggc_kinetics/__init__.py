"""ggc-kinetics: frozen zero-ODE enzyme-kinetics solvers, detectors and a plate-reader parser from GGC-inf."""
from .solvers import (coupled_lag_decay, fit_kp_vs_kd, k57_called, k57_fd_check, k57_map_fit, kp_kd_exact,  # noqa: F401
                      ksi_residual, ksi_solve_conv, morrison_exact, slow_binding_warped, unified_kp_strang)
from .detectors import z_dsf_asym_boltz, z_rao_call, z_rao_orth_robust  # noqa: F401
from .parser import parse_file, process_bytes  # noqa: F401

__version__ = "0.2.0"
