"""Tests of the public API and of behaviour not covered by the specification examples."""

import random

import pytest
from hgvs_apply import apply_nucleotide

import variant_extractor
import variant_extractor.extractor
from variant_extractor import (
    TooManyVariationsException,
    extract_variant_string,
    extract_variants,
    reverse_complement,
    translate,
)


@pytest.mark.parametrize(
    "prefix,ref,alt,options,expected",
    [
        # no change
        ("g.", "GATTACA", "GATTACA", {}, "g.="),
        ("p.", "ATGGATTAA", "ATGGATTAA", {}, "p.="),
        # input: any case, RNA or DNA; output follows the prefix
        ("g.", "gattaca", "gatgaca", {}, "g.4T>G"),
        ("g.", "GAUUACA", "GAUGACA", {}, "g.4T>G"),
        ("r.", "GATTACA", "GATGACA", {}, "r.4u>g"),
        # a prefix may include the reference sequence
        ("NM_004006.2:c.", "AGAAGTAGAGG", "AGAAGTTAGAGG", {"offset": 14}, "NM_004006.2:c.20dup"),
        # c. numbering has no position 0
        ("c.", "CCCATG", "CCCGTG", {"offset": -4}, "c.1A>G"),
        ("c.", "CCCATG", "CCGATG", {"offset": -4}, "c.-1C>G"),
        # several changes
        ("g.", "GATTACA", "GTTTAGA", {}, "g.[2A>T;6C>G]"),
        ("g.", "GATTACA", "GTGTAACA", {"minus_strand": True}, "g.[3dup;5_6delinsCA]"),
        # two changes one nucleotide apart: delins only when they affect one codon
        ("c.", "ATGGCTTAA", "ATGACATAA", {}, "c.4_6delinsACA"),
        ("c.", "ATGGCTTAA", "ATGACATAA", {"merge_one_nt": "never"}, "c.[4G>A;6T>A]"),
        ("g.", "ATGGCTTAA", "ATGACATAA", {}, "g.[4G>A;6T>A]"),
        ("g.", "ATGGCTTAA", "ATGACATAA", {"merge_one_nt": "always"}, "g.4_6delinsACA"),
        # duplication at the 3' end; insertions at either end have only one flank
        ("g.", "GATTACA", "GATTACACA", {}, "g.6_7dup"),
        ("g.", "GATTACA", "CGATTACA", {}, "g.1delinsCG"),
        ("g.", "GATTACA", "GATTACAG", {}, "g.7delinsAG"),
        # inversion
        ("g.", "GATTACA", "GTAATCA", {}, "g.2_5inv"),
        # inserted copies (inverted copies too) are referenced from 10 nucleotides on
        (
            "g.",
            "AAAAAAAAAACGTACGTACGTTTTTTTTTT",
            "AAAAAAAAAACGTACGTACGTTTTTCGTACGTACGTTTTT",
            {},
            "g.25_26ins11_20",
        ),
        (
            "g.",
            "AAAAAAAAAACGTACGTACGTTTTTTTTTT",
            "AAAAAAAAAACGTACGTACGTTTTTCGTACGTACGTTTTT",
            {"min_copy_length": 20},
            "g.25_26insCGTACGTACG",
        ),
        ("g.", "AAAAACGATTGCAACTTTTTCCCCC", "AAAAACGATTGCAACTTTTTGTTGCAATCGCCCCC", {}, "g.20_21ins6_15inv"),
    ],
)
def test_nucleotide(prefix, ref, alt, options, expected):
    """nucleotide descriptions"""
    assert extract_variant_string(prefix, ref, alt, **options) == expected


@pytest.mark.parametrize(
    "ref,alt,options,expected",
    [
        ("ATGGTTGGTTCAGGTTAA", "ATGGCTGGTTCAGCTTAA", {}, "p.[Val2Ala;Gly5Ala]"),
        # offset counts nucleotides; translation starts at the next codon boundary
        ("ATGGTTGGTTCAGGTTAA", "ATGGCTGGTTCAGCTTAA", {"offset": 3}, "p.[Val3Ala;Gly6Ala]"),
        ("CATGGTTGGTTCATAA", "CATGGTTGCTTCATAA", {"offset": 1}, "p.Val4Leu"),
        # synonymous
        ("ATGGTTGGTTCATAA", "ATGGTTGGCTCATAA", {}, "p.Gly3="),
        ("ATGGTTGGTTCATAA", "ATGGTTGGTTCATAG", {}, "p.Ter5="),
        # in-frame deletion, duplication, deletion-insertion
        ("ATGGTTGGTTCATAA", "ATGGTTTCATAA", {}, "p.Gly3del"),
        ("ATGGTTGGTTCATAA", "ATGGTTGGTGGTTCATAA", {}, "p.Gly3dup"),
        # new stop codons
        ("ATGGTTGGTTCATAA", "ATGGTTTGATCATAA", {}, "p.Gly3Ter"),
        ("ATGGTTGGTTCATAA", "ATGGTTGGTTAGTAA", {}, "p.Ser4Ter"),
        ("ATGGTTGGTTCATAA", "ATGGTTCTGTAGTCATAA", {}, "p.Gly3delinsLeuTer"),
        # an insertion containing a stop is an insertion (HGVS protein insertion examples)
        ("ATGAAAATGCATTGGTAA", "ATGAAAATGGGGTAGCATTGGTAA", {}, "p.Met3_His4insGlyTer"),
        # (c.6_7insTCATCATGAGCTCCC, as in the spec's NP_004371.2:p.(Pro46_Asn47insSerSerTer))
        ("ATGCCTAACGGCTTTTAA", "ATGCCTTCATCATGAGCTCCCAACGGCTTTTAA", {}, "p.Pro2_Asn3insSerSerTer"),
        # stop codon changed: extension (with other changes: in an allele)
        ("ATGGTTGGTTCATAA", "ATGGTTGGTTCACAAGGTTGAAAA", {}, "p.Ter5GlnextTer2"),
        (
            "ATGGCTCCAAATCACGTGTACACTCAAACCATTTAA",
            "ATGGCTCCAAATCACATGTACACTCAAACCATTCAA",
            {},
            "p.[Val6Met;Ter12GlnextTer?]",
        ),
        # stop codon unchanged: residues added before it are an insertion or duplication
        ("ATGGTTGGTTCATAA", "ATGGTTGGTTCAGGTTAA", {}, "p.Ser4_Ter5insGly"),
        ("ATGGTTGGTTCATAA", "ATGGTTGGTTCATCATAA", {}, "p.Ser4dup"),
        # a new stop before the lost reference stop
        ("ATGCAGGTGTTGCGCTGCTTGTAA", "ATGCATGTGTTGCGCTGCTAGCAA", {}, "p.[Gln2His;Leu7Ter]"),
        # start codon lost
        ("ATGGTTGGTTCATAA", "ATAGTTGGTTCATAA", {}, "p.Met1?"),
        # frameshift, long and short form; no stop codon in the given sequence: Ter?
        ("ATGGTTGGTTCAAAACCC", "ATGGTTGTTCAAAACCC", {}, "p.Gly3ValfsTer?"),
        ("ATGGTTGGTTCAAAACCC", "ATGGTTGTTCAAAACCC", {"fs_short": True}, "p.Gly3fs"),
        ("ATGGTTGGTTCATAAGG", "ATGGTTGGTTTCATAAGG", {}, "p.Ser4PhefsTer?"),
        # effect beyond the end of the given sequence
        ("ATGGTTGGTTCATAA", "ATGGTTGGTTCATA", {}, "p.?"),
    ],
)
def test_protein(ref, alt, options, expected):
    """protein descriptions"""
    assert extract_variant_string("p.", ref, alt, **options) == expected


def test_extract_variants_returns_descriptions_without_prefix():
    assert extract_variants("g.", "GATTACA", "GTTTAGA") == ["2A>T", "6C>G"]
    assert not extract_variants("g.", "GATTACA", "GATTACA")


def test_max_mutations():
    assert extract_variant_string("g.", "GATTACA", "GTTTAGA", max_mutations=2) == "g.[2A>T;6C>G]"
    with pytest.raises(TooManyVariationsException):
        extract_variant_string("g.", "GATTACA", "GTTTAGA", max_mutations=1)
    assert issubclass(TooManyVariationsException, ValueError)


@pytest.mark.parametrize(
    "args,options",
    [
        (("g.", "GATXACA", "GATTACA"), {}),
        (("g.", "GATTACA", "GA-TACA"), {}),
        (("g.", "", "GATTACA"), {}),
        (("x.", "GATTACA", "GATGACA"), {}),
        (("g", "GATTACA", "GATGACA"), {}),
        (("g.", "GATTACA", "GATGACA"), {"merge_one_nt": "sometimes"}),
    ],
)
def test_invalid_input(args, options):
    with pytest.raises(ValueError):
        extract_variant_string(*args, **options)


def test_helpers():
    assert reverse_complement("AACGTN") == "NACGTT"
    assert translate("ATGGATTAAGGG") == "MD*"
    assert translate("ATGGATTAAGGG", to_stop=False) == "MD*G"


def test_version_matches_metadata():
    importlib_metadata = pytest.importorskip("importlib.metadata")
    try:
        version = importlib_metadata.version("variant_extractor")
    except importlib_metadata.PackageNotFoundError:
        pytest.skip("package not installed")
    assert variant_extractor.__version__ == version


@pytest.mark.parametrize(
    "prefix,ref,alt,options,expected",
    [
        # an insertion at each end of the reference, the first merged with the change after it
        ("g.", "AAA", "CACAAC", {}, "g.[1delinsCAC;3delinsAC]"),
        ("g.", "CA", "ACCAC", {}, "g.1_2delinsACCAC"),
        # a copy reference needs a range of at least two
        ("g.", "CT", "A", {"min_copy_length": 1}, "g.1_2delinsA"),
        # the variant protein continues after the end of the reference translation
        ("p.", "AAA", "AAAAAA", {}, "p.Lys1dup"),
        ("p.", "AAAN", "AAANTAANTT", {"offset": 68}, "p.?"),
        # no complete codon in the reference
        ("p.", "AAA", "AAAAAA", {"offset": 1}, "p.?"),
        # stops inserted after the stop codon do not change the protein
        ("p.", "AAAACAATAA", "AAAAAAACAATAAAAATAA", {"offset": 134}, "p.Lys46dup"),
        (
            "p.",
            "ATGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATAA",
            "ATGATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTATAATAATAA",
            {},
            "p.Phe2Ile",
        ),
    ],
)
def test_regressions_found_by_fuzzing(prefix, ref, alt, options, expected):
    """cases found by the Hypothesis fuzz tests (tests/test_fuzz.py)"""
    assert extract_variant_string(prefix, ref, alt, **options) == expected


@pytest.mark.parametrize("prefix", ["g.", "n.", "m.", "o.", "p."])
def test_negative_offset_needs_c_or_r_coordinates(prefix):
    """found by fuzzing: g.0dup; only c. and r. coordinates go below 1"""
    with pytest.raises(ValueError):
        extract_variant_string(prefix, "ATGAAA", "ATGAAAA", offset=-1)


@pytest.mark.parametrize("region,spacing", [(80, 6), (260, 6)])
def test_large_clusters_of_changes(region, spacing):
    """a densely changed region: searched with fewer inversion candidates (above 50 x 50),
    or not refined at all (above 200), and still described correctly"""
    rng = random.Random(region)
    ref = "".join(rng.choices("ACGT", k=region + 40))
    alt = list(ref)
    for i in range(20, 20 + region, spacing):  # a change every few bases: one cluster
        alt[i] = "ACGT"[("ACGT".index(alt[i]) + 1) % 4]
    alt = "".join(alt)
    description = extract_variant_string("g.", ref, alt)
    assert apply_nucleotide(description, ref) == alt


def test_long_sequence_with_scattered_changes():
    """found by benchmarking: a long changed region is aligned in a band, so changes far
    apart are described separately (not as one delins of everything between them)"""
    rng = random.Random(10)
    ref = "".join(rng.choices("ACGT", k=6000))
    alt = ref[:1000] + "T" + ref[1001:3000] + ref[3002:5000] + "GGA" + ref[5000:]
    if ref[1000] == "T":
        alt = alt[:1000] + "G" + alt[1001:]
    description = extract_variant_string("g.", ref, alt)
    assert apply_nucleotide(description, ref) == alt
    assert description.count(";") == 2, description


@pytest.mark.parametrize("aligner", ["_align_edits", "_align"])
def test_too_different_to_align(aligner, monkeypatch):
    """when the alignment band would be too large, the region is one change"""
    # pylint: disable=protected-access
    monkeypatch.setattr(variant_extractor.extractor, "_edit_distance", lambda ref, alt: 10_000)
    rng = random.Random(11)
    ref, alt = "".join(rng.choices("ACGT", k=1000)), "".join(rng.choices("ACGT", k=1000))
    with pytest.raises(variant_extractor.extractor._RegionTooDifferent):
        getattr(variant_extractor.extractor, aligner)(ref, alt)
    assert extract_variant_string("g.", ref, alt).startswith("g.1_1000delins")
