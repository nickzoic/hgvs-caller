"""Comparison with VarSim's extraction, on the fixture tests/data/varsim_extract_results.json.

VarSim (https://github.com/liu-sun/VarSim, release 2.0.1, MIT) was run on our test cases by
tools/make_varsim_fixture.py; the fixture records its version and commit. These tests do not
need VarSim installed.

* Both tools' descriptions give the variant sequence.
* Where the descriptions differ, the difference is one of the known kinds: VarSim describes
  duplications as insertions and inversions as delins (it never writes dup or inv), and it describes
  several changes as one delins covering them (it finds only the first and last differing bases).
  Any other difference is a disagreement to investigate.
"""

import json
import pathlib
import re

import pytest
from hgvs_apply import apply_nucleotide, norm

from hgvs_extractor import extract_variant_string

FIXTURE = json.loads((pathlib.Path(__file__).parent / "data" / "varsim_extract_results.json").read_text())
ENTRIES = FIXTURE["entries"]
_PARAMS = [pytest.param(e, id=f"{i}-{e['source'][:12]}") for i, e in enumerate(ENTRIES)]


def _strip(description: str) -> str:
    """'g.3G>T' -> '3G>T' (the fixture stores VarSim descriptions with the g. prefix)"""
    return description[2:] if description.startswith("g.") else description


def _span(description: str) -> tuple[int, int]:
    """(first, last) reference position a description covers"""
    numbers = [int(n) for n in re.findall(r"\d+", description)]
    return numbers[0], numbers[-1] if len(numbers) > 1 else numbers[0]


def explain(ours: str, varsim: str, ref: str) -> str:
    """the known kind of difference between our description and VarSim's, or '' if unexplained"""
    if ours == varsim:
        return "same"
    ours_edits = ours.strip("[]").split(";")
    if len(ours_edits) > 1 and ";" not in varsim:  # several changes, one delins (VarSim)
        return "several changes, VarSim one delins"
    if m := re.fullmatch(r"(\d+)(?:_(\d+))?dup", ours):  # duplication: VarSim inserts the same bases
        a, b = int(m.group(1)), int(m.group(2) or m.group(1))
        m_ins = re.fullmatch(r"(\d+)_(\d+)ins([ACGT]+)", varsim)
        if m_ins and m_ins.group(3) == ref[a - 1 : b]:
            return "duplication, VarSim insertion"
    if m := re.fullmatch(r"(\d+)_(\d+)ins(\d+)_(\d+)", ours):  # insertion of a copy of the reference
        s, e = int(m.group(3)), int(m.group(4))
        m_ins = re.fullmatch(r"(\d+)_(\d+)ins([ACGT]+)", varsim)
        if m_ins and m_ins.group(1) == m.group(1) and m_ins.group(3) == ref[s - 1 : e]:
            return "insertion of a reference copy, VarSim literal bases"
    m_end = re.fullmatch(r"(\d+)_(\d+)ins[ACGT]+", varsim)  # insertion at an end, VarSim position outside
    m_ours = re.fullmatch(r"(\d+)delins[ACGT]+", ours)
    if m_end and m_ours:
        n = len(ref)
        first, last = int(m_end.group(1)), int(m_end.group(2))
        if (first == 0 and int(m_ours.group(1)) == 1) or (last == n + 1 and int(m_ours.group(1)) == n):
            return "insertion at a reference end, VarSim position outside the reference"
    if m := re.fullmatch(r"(\d+)_(\d+)inv", ours):  # inversion: VarSim a delins of the same range
        if re.fullmatch(rf"{m.group(1)}_{m.group(2)}delins[ACGT]+", varsim):
            return "inversion, VarSim delins"
    return ""


@pytest.mark.parametrize("entry", _PARAMS)
def test_varsim_description_gives_the_variant(entry):
    """VarSim's description, applied to the reference, gives the observed sequence"""
    varsim = entry["varsim"]
    if varsim == "=":
        assert norm(entry["ref"]) == norm(entry["alt"])
    else:
        assert apply_nucleotide("g." + _strip(varsim), entry["ref"]) == norm(entry["alt"]), varsim


@pytest.mark.parametrize("entry", _PARAMS)
def test_our_description_agrees_with_varsim_or_is_a_known_difference(entry):
    """our description is the same as VarSim's, or differs in a known way (see explain)"""
    ours = extract_variant_string("g.", entry["ref"], entry["alt"])[2:]
    varsim = _strip(entry["varsim"]) if entry["varsim"] != "=" else "="
    if ours == "=" or varsim == "=":
        assert ours == varsim
        return
    reason = explain(ours, varsim, norm(entry["ref"]))
    assert reason, f"unexplained difference: ours {ours}, VarSim {varsim}"


def test_fixture_metadata():
    """the fixture records the VarSim version it was made with"""
    assert FIXTURE["tool"] == "VarSim extract_with_positions"
    assert FIXTURE["varsim_version"] == "2.0.1"
    assert len(ENTRIES) > 900
