"""
GGC-inf Release Phase 6 / 6b (v0.2.0): clean-venv wheel smoke test. Copy of ggc_inf_kit/release_wheel_smoke.py
(Phase 2); changes: package root and log paths, version 0.2.0 required in S1, and S4 expects 18 tests (14 + 4 guard).
Usage: python -I -X utf8 benchmarks/release_wheel_smoke_v020.py TAG VENV_PYTHON [PREFIX]   (PREFIX default phase6; Phase 6b uses phase6b)
  VENV_PYTHON is the python.exe of a fresh venv that has ONLY the built wheel installed (plus numpy / scipy).
Checks (all ASSUMED pass/fail rules; the plates in tests/data are MODELED):
  S1 installed package is imported from the venv's site-packages (not from src/).
  S2 ggc-fit (the installed console script) on tests/data/plate_healthy_biotek_wide.csv with --evidence MODELED
     --scale 0.05 --offset 0.040 exits 0 and its JSON report equals the stored Cycle 34 report byte for byte.
  S3 ggc-fit on tests/data/bad_no_wells.csv exits 2.
  S4 python -I -X utf8 -m unittest discover -s tests passes 18/18.
  S5 python -I -X utf8 examples/quickstart.py exits 0.
Writes ../ggc_inf_kit/logs/PREFIX_wheel_smoke_TAG.txt. Only the tag and wall-time lines may differ between run1 and run2.
"""
import os
import re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.dirname(HERE)
LOGS = os.path.join(os.path.dirname(REL), "ggc_inf_kit", "logs")
DATA = os.path.join(REL, "tests", "data")


def run(cmd, cwd=REL):
    env = dict(os.environ, PYTHONIOENCODING="utf8")
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, encoding="utf-8", env=env)
    return p.returncode, p.stdout, p.stderr


def main():
    tag, py = sys.argv[1], os.path.abspath(sys.argv[2])
    prefix = sys.argv[3] if len(sys.argv) > 3 else "phase6"
    t0 = time.perf_counter()
    scripts = os.path.dirname(py)
    exe = os.path.join(scripts, "ggc-fit.exe" if os.name == "nt" else "ggc-fit")
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE {prefix.upper()} (v0.2.0): clean-venv wheel smoke test | {tag}")
    A("=" * 78)
    ok = []
    rc, out, err = run([py, "-I", "-c", "import ggc_kinetics as g, os; print(g.__version__); "
                        "print(os.path.relpath(g.__file__, os.path.dirname(os.path.dirname(os.__file__))))"])
    lines = out.split()
    s1 = rc == 0 and lines[:1] == ["0.2.0"] and "site-packages" in out and "src" not in out.replace("\\", "/").split("/")
    A(f"  S1 import: version {lines[0] if lines else '?'}, from {lines[1].replace(os.sep, '/') if len(lines) > 1 else '?'}"
      f" ... {'PASS' if s1 else 'FAIL'}")
    ok.append(s1)
    with tempfile.TemporaryDirectory() as d:
        rc, out, err = run([exe, os.path.join(DATA, "plate_healthy_biotek_wide.csv"), "--evidence", "MODELED",
                            "--scale", "0.05", "--offset", "0.040", "--out", d])
        rep = os.path.join(d, "plate_healthy_biotek_wide.report.json")
        got = open(rep, encoding="utf-8").read() if os.path.exists(rep) else ""
        ref = open(os.path.join(DATA, "plate_healthy.report.json"), encoding="utf-8").read()
        s2 = rc == 0 and got == ref
        A(f"  S2 ggc-fit healthy plate: exit {rc}, report == stored Cycle 34 report: {got == ref} "
          f"({len(got)} bytes) ... {'PASS' if s2 else 'FAIL'}")
        A(f"     stdout: {out.strip().splitlines()[-1] if out.strip() else ''}")
        ok.append(s2)
        rc, out, err = run([exe, os.path.join(DATA, "bad_no_wells.csv"), "--out", d])
        s3 = rc == 2
        A(f"  S3 ggc-fit bad_no_wells.csv: exit {rc} (want 2) ... {'PASS' if s3 else 'FAIL'}")
        ok.append(s3)
    rc, out, err = run([py, "-I", "-X", "utf8", "-m", "unittest", "discover", "-s", "tests", "-v"])
    m = re.search(r"Ran (\d+) tests? in [\d.]+s", err)
    n = int(m.group(1)) if m else 0
    res_ok = "\nOK" in err
    s4 = rc == 0 and n == 18 and res_ok
    A(f"  S4 unittest (installed package): ran {n}, result {'OK' if res_ok else 'NOT OK'} ... "
      f"{'PASS' if s4 else 'FAIL'}")
    for ln in err.splitlines():
        if ln.endswith("... ok") or ln.endswith("... FAIL") or ln.endswith("... ERROR"):
            A(f"     {ln}")
    ok.append(s4)
    rc, out, err = run([py, "-I", "-X", "utf8", os.path.join("examples", "quickstart.py")])
    s5 = rc == 0
    A(f"  S5 examples/quickstart.py: exit {rc} ... {'PASS' if s5 else 'FAIL'}")
    for ln in out.splitlines():
        if not ln.startswith("wall time"):
            A(f"     {ln}")
    ok.append(s5)
    A(f"  SMOKE {sum(ok)}/{len(ok)}")
    A(f"  wall time {time.perf_counter() - t0:.2f} s (MEASURED)")
    os.makedirs(LOGS, exist_ok=True)
    open(os.path.join(LOGS, f"{prefix}_wheel_smoke_{tag}.txt"), "w", encoding="utf-8", newline="\n").write("\n".join(L) + "\n")
    print("\n".join(L))
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
