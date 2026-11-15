"""Test variant calling against examples from the HGVS Nomenclature specification.

Examples are taken from HGVS Nomenclature 21.1.5 (https://hgvs-nomenclature.org/stable/,
text CC0), from the "Examples" sections and the Notes / Q&A passages, selecting only
those where the text gives both the reference and the variant sequence. Each example
was checked against the spec: applying the spec's description to the reference
reproduces the variant sequence.

Adaptations, since the caller works on two nucleotide sequences:

* only the part after the accession (``NM_004006.2:``) is compared;
* ``offset`` is set so that the first base of ``ref`` gets the coordinate the spec gives it;
* protein examples given as amino acids are back-translated with one fixed codon per
  amino acid, so no synonymous changes are introduced;
* p.Glu5ValfsTer5 and p.Glu5fs are the same change in the long and short frameshift
  forms, so the short form is requested with ``fs_short=True``.
"""

import re

import pytest

from hgvs_apply import AA1, apply_nucleotide, apply_protein, norm

from hgvs_caller import find_variant_string

SPEC = "https://hgvs-nomenclature.org/stable/recommendations/"

# (expected, kind, ref, alt, offset, source page, options for find_variant_string)
#   kind: "dna" / "rna"   -- ref and alt are nucleotide sequences
#         "protein"       -- ref and alt are amino acid sequences (3-letter codes)
#         "dna-protein"   -- protein description, ref and alt given as coding DNA
#   offset: passed to find_variant_string (position of first base of ref, minus 1;
#           for protein calling it is a nucleotide offset)
EXAMPLES = [
    # --- deletion
    ("g.33344591del", "dna", "TGTGATTCT", "TGTGTTCT", 33344586, "DNA/deletion/", {}),
    ("c.5697del", "dna", "ATTGAAAAAAAATTAG", "ATTGAAAAAAATTAG", 5685, "DNA/deletion/", {}),
    ("g.32343183del", "dna", "CTAATTTTTTTTCAAT", "CTAATTTTTTTCAAT", 32343171, "DNA/deletion/", {}),
    ("g.33344590_33344592del", "dna", "GTGTGATTCTG", "GTGTTCTG", 33344585, "DNA/deletion/", {}),
    ("g.5del", "dna", "TGGGGCAT", "TGGGCAT", 0, "general/", {}),
    # --- duplication
    ("c.20dup", "dna", "AGAAGTAGAGG", "AGAAGTTAGAGG", 14, "DNA/duplication/", {}),
    ("c.5697dup", "dna", "ATTGAAAAAAAATTAG", "ATTGAAAAAAAAATTAG", 5685, "DNA/duplication/", {}),
    ("g.32343183dup", "dna", "CTAATTTTTTTTCAAT", "CTAATTTTTTTTTCAAT", 32343171, "DNA/duplication/", {}),
    ("c.20_23dup", "dna", "AGAAGTAGAGG", "AGAAGTAGATAGAGG", 14, "DNA/duplication/", {}),
    # --- substitution
    ("c.2077G>A", "dna", "CATGACA", "CATAACA", 2073, "DNA/substitution/", {}),
    ("r.2077g>a", "rna", "caugaca", "cauaaca", 2073, "RNA/delins/", {}),
    # --- insertion
    (
        "g.17_18ins5_16",
        "dna",
        "ATCGATCGATCGATCGAGGGTCCC",
        "ATCGATCGATCGATCGAATCGATCGATCGGGGTCCC",
        0,
        "DNA/insertion/",
        {},
    ),
    (
        "r.17_18ins5_16",
        "rna",
        "aucgaucgaucgaucgaggguccc",
        "aucgaucgaucgaucgaaucgaucgaucgggguccc",
        0,
        "RNA/insertion/",
        {},
    ),
    # --- deletion-insertion
    ("g.32386323delinsGA", "dna", "CAGCTCTTT", "CAGCGACTTT", 32386318, "DNA/delins/", {}),
    ("c.6775_6777delinsC", "dna", "GGAAGAGTTGC", "GGAACTTGC", 6770, "DNA/delins/", {}),
    ("g.12_13delinsTG", "dna", "GAAGCCAG", "GAATGCAG", 8, "DNA/substitution/", {}),
    ("g.4_5delinsTG", "dna", "TGTGCCA", "TGTTGCA", 0, "DNA/delins/", {}),
    ("c.2077delinsATA", "dna", "CATGACA", "CATATAACA", 2073, "DNA/delins/", {}),
    ("r.775delinsga", "rna", "aggcucauu", "aggcgacauu", 770, "RNA/delins/", {}),
    ("r.775_777delinsc", "rna", "aggcucauu", "aggccuu", 770, "RNA/delins/", {}),
    ("r.5_6delinsug", "rna", "augugcca", "auguugca", 0, "RNA/delins/", {}),
    ("r.2077delinsaua", "rna", "caugaca", "cauauaaca", 2073, "RNA/delins/", {}),
    # --- inversion
    ("g.32361330_32361333inv", "dna", "CATCAGCCT", "CACTGACCT", 32361327, "DNA/inversion/", {}),
    ("c.5657_5660inv", "dna", "AGGCTGATG", "AGGTCAGTG", 5653, "DNA/inversion/", {}),
    ("r.177_180inv", "rna", "aggcugauu", "aggucaguu", 173, "RNA/inversion/", {}),
    # --- protein
    ("p.Trp4del", "protein", "MetLeuTrpTrpGlu", "MetLeuTrpGlu", 0, "protein/deletion/", {}),
    ("p.Ser5del", "protein", "MetTrpSerSerSerHisAsp", "MetTrpSerSerHisAsp", 0, "protein/deletion/", {}),
    ("p.Trp4dup", "protein", "MetLeuTrpTrpGlu", "MetLeuTrpTrpTrpGlu", 0, "protein/duplication/", {}),
    ("p.Ser5dup", "protein", "MetTrpSerSerSerHisAsp", "MetTrpSerSerSerSerHisAsp", 0, "protein/duplication/", {}),
    (
        "p.His4_Gln5insAla",
        "protein",
        "MetLysGlyHisGlnGlnCys",
        "MetLysGlyHisAlaGlnGlnCys",
        0,
        "protein/insertion/",
        {},
    ),
    (
        "p.Lys2_Gly3insGlnSerLys",
        "protein",
        "MetLysGlyHisGlnGlnCys",
        "MetLysGlnSerLysGlyHisGlnGlnCys",
        0,
        "protein/insertion/",
        {},
    ),
    (
        "p.His7_Gln8insGly4_Ser6",
        "protein",
        "MetArgThrGlySerSerHisGlnTrpPhe",
        "MetArgThrGlySerSerHisGlySerSerGlnTrpPhe",
        0,
        "protein/insertion/",
        {},
    ),
    # codon only; placed so that it is residue 188
    ("p.Cys188=", "dna-protein", "TGC", "TGT", 187 * 3, "protein/substitution/", {}),
    # coding DNA starting at the ATG (c.1)
    (
        "p.Glu5ValfsTer5",
        "dna-protein",
        "ATGGATGCATACGAGATGAGG",
        "ATGGATGCATACGTGCATACGAGATGAGG",
        0,
        "protein/frameshift/",
        {},
    ),
    (
        "p.Glu5fs",
        "dna-protein",
        "ATGGATGCATACGAGATGAGG",
        "ATGGATGCATACGTGCATACGAGATGAGG",
        0,
        "protein/frameshift/",
        {"fs_short": True},
    ),
]

CODON_FOR = {
    "Ala": "GCT", "Arg": "CGT", "Asn": "AAT", "Asp": "GAT", "Cys": "TGT", "Gln": "CAA", "Glu": "GAA",
    "Gly": "GGT", "His": "CAT", "Ile": "ATT", "Leu": "CTG", "Lys": "AAA", "Met": "ATG", "Phe": "TTT",
    "Pro": "CCT", "Ser": "TCT", "Thr": "ACT", "Trp": "TGG", "Tyr": "TAT", "Val": "GTT",
}  # fmt: skip


def back_translate(protein: str) -> str:
    """amino acid sequence in 3-letter codes -> DNA, one fixed codon per amino acid"""
    return "".join(CODON_FOR[aa] for aa in re.findall(r"[A-Z][a-z]{2}", protein))


def call(expected: str, kind: str, ref: str, alt: str, offset: int, options: dict) -> str:
    """call the variant with find_variant_string"""
    if kind == "protein":
        ref, alt = back_translate(ref), back_translate(alt)
    return find_variant_string(expected[:2], ref, alt, offset=offset, **options)


def _params(examples):
    return [
        pytest.param(expected, kind, ref, alt, offset, SPEC + page, options, id=expected)
        for expected, kind, ref, alt, offset, page, options in examples
    ]


@pytest.mark.parametrize("expected,kind,ref,alt,offset,source,options", _params(EXAMPLES))
def test_hgvs_spec_example(expected, kind, ref, alt, offset, source, options):
    """the caller gives the same description as the HGVS specification"""
    assert call(expected, kind, ref, alt, offset, options) == expected, f"see {source}"


# Independent of the exact description: applying the called variant to the reference
# must give the variant sequence.


@pytest.mark.parametrize(
    "expected,kind,ref,alt,offset,source,options", _params(e for e in EXAMPLES if e[1] in ("dna", "rna"))
)
def test_called_variant_reproduces_sequence(expected, kind, ref, alt, offset, source, options):
    """applying the called variant to the reference gives the variant sequence"""
    called = call(expected, kind, ref, alt, offset, options)
    assert apply_nucleotide(called, ref, offset) == norm(alt), f"called {called}, see {source}"


def test_apply_helper_with_spec_descriptions():
    """sanity check of the test helpers: the spec's own descriptions reproduce the variant sequences"""
    for expected, kind, ref, alt, offset, _, _ in EXAMPLES:
        if kind in ("dna", "rna"):
            assert apply_nucleotide(expected, ref, offset) == norm(alt), expected


def test_apply_protein_helper_with_spec_descriptions():
    """sanity check of apply_protein with the spec's in-frame protein examples"""

    def one_letter(protein: str) -> str:
        return "".join(AA1[aa] for aa in re.findall(r"[A-Z][a-z]{2}", protein))

    for expected, kind, ref, alt, _, _, _ in EXAMPLES:
        if kind == "protein":
            assert apply_protein(expected, one_letter(ref)) == one_letter(alt), expected
