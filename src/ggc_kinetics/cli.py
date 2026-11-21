"""ggc-fit: parse a plate-reader export and write <name>.report.json and <name>.report.md.

Usage: ggc-fit FILE [--layout LAYOUT.json] [--mode substrate|product] [--evidence MEASURED|MODELED] [--out DIR]
       [--scale A_PER_50uM] [--offset A] [--wavelength-index N]
Exit codes: 0 report written with status OK; 2 parse error or unreadable file (an error report is still written
when the file could be read).
"""
import argparse
import json
import os
import sys

from .parser import dumps, parse_file, report_markdown


def main(argv=None):
    ap = argparse.ArgumentParser(prog="ggc-fit", description="96-well kinetic export -> GGC contract report")
    ap.add_argument("path")
    ap.add_argument("--layout", default=None, help="JSON map of contract columns to source wells")
    ap.add_argument("--mode", choices=("substrate", "product"), default="substrate")
    ap.add_argument("--evidence", choices=("MEASURED", "MODELED"), default="MEASURED",
                    help="MEASURED for a real plate, MODELED for a simulated one")
    ap.add_argument("--out", default=None, help="output directory (default: next to the input)")
    ap.add_argument("--scale", type=float, default=1.0, help="absorbance per 50 uM substrate (ASSUMED)")
    ap.add_argument("--offset", type=float, default=0.0, help="background absorbance (ASSUMED)")
    ap.add_argument("--wavelength-index", type=int, default=0, help="SoftMax exports: which wavelength")
    a = ap.parse_args(argv)
    try:
        rep, _ = parse_file(a.path, layout=a.layout, evidence=a.evidence, mode=a.mode, scale=a.scale,
                            offset=a.offset, wavelength_index=a.wavelength_index)
    except OSError as e:
        print(json.dumps({"status": "ERROR", "error": {"code": "IO", "message": str(e)}}))
        return 2
    out = a.out or os.path.dirname(os.path.abspath(a.path))
    os.makedirs(out, exist_ok=True)
    stem = os.path.join(out, os.path.splitext(os.path.basename(a.path))[0] + ".report")
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
    sys.exit(main())
