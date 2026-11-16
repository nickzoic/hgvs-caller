"""Test cases and results from Mutalyzer's description-extractor.

https://github.com/mutalyzer/description-extractor (commit 155b79d), whose
``describe_dna(reference, sample)`` gives HGVS descriptions without a prefix.

* EXAMPLES: description-extractor's tests/test_describe.py.  Its tests/test_extractor.py checks the
  same sequence pairs at the level of its internal alignment, and its protein tests are
  disabled (and take protein sequences), so they are not included.
* data/description_extractor_results.json: 500 random sequence pairs (from
  tools/compare_description_extractor.py --seed 1) for which description-extractor and variant_extractor give
  the same description.  The tool reports where they differ: mostly where description-extractor joins
  separate changes in one delins, which HGVS describes separately.
"""

import json
import pathlib

import pytest
from hgvs_apply import apply_nucleotide

from variant_extractor import extract_variant_string

REFERENCE = "ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT"  # used by test5 to test18

# (test name in test_describe.py, reference, sample, description-extractor's description)
EXAMPLES = [
    (
        "test1",
        "ATGATGATCAGATACAGTGTGATACAGGTAGTTAGACAA",
        "ATGATTTGATCAGATACATGTGATACCGGTAGTTAGGACAA",
        "[5_6insTT;17del;26A>C;35dup]",
    ),
    (
        "test2",
        "TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTACTGTG",
        "TAAGCACCAGGAGTCCATGAAGAAGCTGGATCCTCCCATGGAATCCCCTACTCTACTGTG",
        "[26A>C;30C>A;35G>C]",
    ),
    (
        "test3",
        "TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA",
        "TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGGAATCCCCTACTCTA",
        "[26_29inv;30C>G]",
    ),
    (
        "test4",
        "TAAGCACCAGGAGTCCATGAAGAAGATGGCTCCTGCCATGGAATCCCCTACTCTA",
        "TAAGCACCAGGAGTCCATGAAGAAGCCATGTCCTGCCATGAATCCCCTACTCTA",
        "[26_29inv;30C>G;41del]",
    ),
    ("test5", REFERENCE, "ACGTCGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "="),
    ("test6", REFERENCE, "ACGTCGGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7A>G"),
    ("test7", REFERENCE, "ACGTCGTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7del"),
    ("test8", REFERENCE, "ACGTCGTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7_8del"),
    ("test9", REFERENCE, "ACGTCGCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "6_7insC"),
    ("test10", REFERENCE, "ACGTCGCCATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "6_7insCC"),
    ("test11", REFERENCE, "ACGTCGAATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7dup"),
    ("test12", REFERENCE, "ACGTCGAGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "6_7dup"),
    ("test13", REFERENCE, "ACGTCGACGATTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "5_7dup"),
    ("test14", REFERENCE, "ACGTCGCGAATCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7_11inv"),
    ("test15", REFERENCE, "ACGTCGCCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7delinsCC"),
    ("test16", REFERENCE, "ACGTCGATTCGCTAGCTTCGTTTTGATAGATAGAGATATAGAGAT", "21_23delinsTTTT"),
    ("test17", REFERENCE, "ACGTCTCTTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "6_7inv"),
    ("test18", REFERENCE, "ACGTCGTCTCGCTAGCTTCGGGGGATAGATAGAGATATAGAGAT", "7_8delinsTC"),
]


@pytest.mark.parametrize("name,ref,alt,expected", EXAMPLES, ids=[e[0] for e in EXAMPLES])
def test_description_extractor_example(name, ref, alt, expected):  # pylint: disable=unused-argument
    """same description as description-extractor"""
    assert extract_variant_string("g.", ref, alt) == "g." + expected


@pytest.mark.parametrize("name,ref,alt,expected", EXAMPLES, ids=[e[0] for e in EXAMPLES])
def test_description_extractor_example_reproduces_sample(name, ref, alt, expected):  # pylint: disable=unused-argument
    """sanity check of the examples: description-extractor's description gives the sample"""
    assert apply_nucleotide("g." + expected, ref) == alt


RESULTS = json.loads((pathlib.Path(__file__).parent / "data" / "description_extractor_results.json").read_text())


@pytest.mark.parametrize("ref,alt,expected", RESULTS)
def test_description_extractor_result(ref, alt, expected):
    """same description as description-extractor gave"""
    assert extract_variant_string("g.", ref, alt) == "g." + expected
