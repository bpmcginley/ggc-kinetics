"""ggc-kinetics demo app. Run: pip install "ggc-kinetics[app]"; streamlit run app/streamlit_app.py

Upload a plate-reader export, see the contract report, and plot the parsed traces. Every number is shown with the
evidence tag you choose (MEASURED for a real plate). No fitting knobs are exposed: the rules are frozen.
"""
import json

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from ggc_kinetics import __version__
from ggc_kinetics.parser import dumps, process_bytes, report_markdown

st.set_page_config(page_title="ggc-kinetics", layout="wide")
st.title(f"ggc-kinetics v{__version__} plate report")
st.caption("Frozen GGC-inf rules. Thresholds are ASSUMED values calibrated on MODELED plates.")

up = st.file_uploader("Plate export (CSV / TSV / TXT)", type=["csv", "tsv", "txt"])
lay = st.file_uploader("Optional layout JSON", type=["json"])
c1, c2, c3, c4 = st.columns(4)
evidence = c1.selectbox("Evidence", ["MEASURED", "MODELED"])
mode = c2.selectbox("Readout", ["substrate", "product"])
scale = c3.number_input("Absorbance per 50 uM (ASSUMED)", value=1.0, format="%.4f")
offset = c4.number_input("Background absorbance (ASSUMED)", value=0.0, format="%.4f")

if up is not None:
    rep, parsed = process_bytes(up.getvalue(), up.name, layout_raw=lay.getvalue() if lay else None,
                                layout_name=lay.name if lay else None, scale=scale, offset=offset, mode=mode,
                                evidence=evidence)
    if rep["status"] != "OK":
        st.error(f"{rep['error']['code']}: {rep['error']['message']}")
    else:
        st.markdown(report_markdown(rep))
        if parsed is not None:
            fig, ax = plt.subplots(figsize=(7, 4))
            t = np.asarray(parsed["times_s"]) / 60.0
            A, present = np.asarray(parsed["A"]), np.asarray(parsed["present"])
            for r in range(A.shape[0]):
                for c in range(A.shape[1]):
                    if present[r, c]:
                        ax.plot(t, A[r, c], lw=0.7)
            ax.set_xlabel("time (min)")
            ax.set_ylabel("absorbance")
            st.pyplot(fig)
    st.download_button("Download report JSON", dumps(rep), file_name="report.json")
    with st.expander("Raw report"):
        st.code(json.dumps(rep, indent=1)[:20000])
