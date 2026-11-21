"""Plate-reader export parser: a thin wrapper over the frozen GGC-inf plate_parser.py (Cycles 33-40).

_plate_parser.py is a byte-for-byte copy of the archive's plate_parser.py (its sha256 is pinned in the tests).
It reads CSV / TSV / text exports, long or wide, with comma or decimal-comma numbers, SoftMax Pro text blocks,
and instrument preambles (Gen5, Tecan and similar): the wide reader starts at the first row that holds at least
8 well names (A1..H12), so header lines above the table are skipped. It does not read .xlsx files.
Errors come back as a report with status ERROR and a code (EMPTY, NO_WELL_TABLE, TIME_NOT_MONOTONIC, NO_READS...).
"""
import os

from . import _plate_parser as _pp

process_bytes = _pp.process_bytes
dumps = _pp.dumps
report_markdown = _pp.report_markdown
z_rao_orthogonal_robust = _pp.z_rao_orthogonal_robust
PlateParseError = _pp.PlateParseError

__all__ = ["parse_file", "process_bytes", "dumps", "report_markdown", "PlateParseError"]


def parse_file(path, layout=None, evidence="MEASURED", mode="substrate", scale=1.0, offset=0.0,
               wavelength_index=0, p1_rule=None):
    """Parse one export and build the Cycle 29 contract report. Returns (report dict, parsed or None).
    evidence: 'MEASURED' for a real plate, 'MODELED' for a simulated one (the tag goes into every number)."""
    with open(path, "rb") as fh:
        raw = fh.read()
    lay_raw = None
    if layout:
        with open(layout, "rb") as fh:
            lay_raw = fh.read()
    return process_bytes(raw, os.path.basename(path), wavelength_index=wavelength_index, layout_raw=lay_raw,
                         layout_name=os.path.basename(layout) if layout else None, scale=scale, offset=offset,
                         mode=mode, evidence=evidence, p1_rule=p1_rule)
