"""
GGC-inf Release Phase 6b (v0.2.0-rc2): headless Streamlit smoke test. Copy of ggc_inf_kit/release_streamlit_smoke.py
(Phase 2); changes: package root and log paths, and T2 also requires the rendered st.title to contain v0.2.0.
Usage: python -I -X utf8 benchmarks/release_streamlit_smoke_v020.py TAG    (run with the smoke venv's python, streamlit installed)

streamlit.testing.v1.AppTest has no API for st.file_uploader, so the harness script below replaces
st.file_uploader with a stub that returns the bytes of tests/data/plate_healthy_biotek_wide.csv (MODELED plate)
for the plate widget and None for the layout widget, then runs the unchanged app file. The app's own number inputs
are set to the plate's calibration (scale 0.05, offset 0.040, ASSUMED) and the evidence box to MODELED.
Checks (ASSUMED rules):
  T1 the app runs with zero exceptions and no st.error.
  T2 the rendered title contains v0.2.0 and the rendered Markdown contains the report (mentions P1 and P5).
  T3 the P1..P5 verdicts shown equal the stored Cycle 34 report (P1 KILL, P2-P4 INCONCLUSIVE, P5 PASS).
  T4 a matplotlib figure was sent to st.pyplot (traces plotted).
Writes ../ggc_inf_kit/logs/phase6b_streamlit_smoke_TAG.txt. Only the tag and wall-time lines may differ between run1 and run2.
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REL = os.path.dirname(HERE)
LOGS = os.path.join(os.path.dirname(REL), "ggc_inf_kit", "logs")
APP = os.path.join(REL, "app", "streamlit_app.py")
PLATE = os.path.join(REL, "tests", "data", "plate_healthy_biotek_wide.csv")
REF = os.path.join(REL, "tests", "data", "plate_healthy.report.json")

HARNESS = '''
import io, runpy
import streamlit as st
import matplotlib
matplotlib.use("Agg")
_raw = open(PLATE_PATH, "rb").read()
class _Up(io.BytesIO):
    name = "plate_healthy_biotek_wide.csv"
    def getvalue(self):
        return _raw
def _fu(label, *a, **k):
    return _Up(_raw) if label.startswith("Plate export") else None
st.file_uploader = _fu
_figs = []
_orig = st.pyplot
def _pp(fig=None, *a, **k):
    _figs.append(fig)
    st.session_state["_n_figs"] = len(_figs)
    return _orig(fig, *a, **k)
st.pyplot = _pp
runpy.run_path(APP_PATH, run_name="__main__")
'''


def main():
    tag = sys.argv[1]
    t0 = time.perf_counter()
    L = []
    A = L.append
    A("=" * 78)
    A(f"GGC-inf RELEASE PHASE 6b (v0.2.0-rc2): headless Streamlit smoke test | {tag}")
    A("=" * 78)
    try:
        import streamlit
        from streamlit.testing.v1 import AppTest
    except Exception as e:  # noqa: BLE001
        A(f"  streamlit not importable ({type(e).__name__}); smoke test NOT RUN, README keeps the compile-only note")
        out = os.path.join(LOGS, f"phase6b_streamlit_smoke_{tag}.txt")
        open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
        print("\n".join(L))
        return 1
    import matplotlib
    A(f"  streamlit {streamlit.__version__}, matplotlib {matplotlib.__version__}; upload stubbed (AppTest has no "
      f"file_uploader API)")
    src = HARNESS.replace("PLATE_PATH", repr(PLATE)).replace("APP_PATH", repr(APP))
    at = AppTest.from_string(src, default_timeout=60)
    at.run()
    at.selectbox[0].set_value("MODELED")
    at.number_input[0].set_value(0.05)
    at.number_input[1].set_value(0.040)
    at.run()
    exc = list(at.exception)
    errs = [e.value for e in at.error]
    t1 = not exc and not errs
    A(f"  T1 exceptions {len(exc)}, st.error {len(errs)} ... {'PASS' if t1 else 'FAIL'}")
    for e in exc:
        A(f"     exception: {e.value}")
    for e in errs:
        A(f"     error: {e}")
    md = "\n".join(m.value for m in at.markdown)
    title = at.title[0].value if len(at.title) else ""
    t2 = "v0.2.0" in title and len(md) > 0 and "P1" in md and "P5" in md
    A(f"  T2 title {title!r} contains v0.2.0: {'v0.2.0' in title}; report markdown: {len(md)} chars, "
      f"mentions P1 and P5: {'P1' in md and 'P5' in md} ... "
      f"{'PASS' if t2 else 'FAIL'}")
    ref = json.load(open(REF, encoding="utf-8"))
    want = {k: v["verdict"] for k, v in ref["tests"].items()}
    raw = "\n".join(c.value for c in at.code)
    got = {}
    try:
        rep = json.loads(raw)
        got = {k: v["verdict"] for k, v in rep["tests"].items()}
        ev = rep["tests"]["P1"]["evidence"]
    except Exception:  # noqa: BLE001
        ev = "?"
    t3 = got == want and ev == "MODELED"
    A(f"  T3 verdicts shown {got} (evidence {ev}); stored {want} ... {'PASS' if t3 else 'FAIL'}")
    nf = at.session_state["_n_figs"] if "_n_figs" in at.session_state else 0
    t4 = nf >= 1
    A(f"  T4 figures sent to st.pyplot: {nf} ... {'PASS' if t4 else 'FAIL'}")
    ok = [t1, t2, t3, t4]
    A(f"  STREAMLIT SMOKE {sum(ok)}/{len(ok)}")
    A(f"  wall time {time.perf_counter() - t0:.2f} s (MEASURED)")
    out = os.path.join(LOGS, f"phase6b_streamlit_smoke_{tag}.txt")
    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
