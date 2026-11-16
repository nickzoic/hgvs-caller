"""Fuzz tests with Hypothesis.

Hypothesis generates inputs, including unrelated sequence pairs, repetitive sequences,
ambiguous bases (N) and arbitrary text, and shrinks any failure to a minimal example.

* robustness: any input gives a well-formed description or a ValueError, nothing else
* grammar: every description parses as HGVS
* correctness: descriptions reproduce the variant sequence / protein
* metamorphic properties: case, T/U, strand, offset and option invariants

Run more examples with ``pytest tests/test_fuzz.py --hypothesis-profile=fuzz``.
"""

import re

from hgvs_apply import apply_nucleotide, check_protein, revcomp
from hypothesis import assume, given
from hypothesis import strategies as st

import hgvs_caller.caller as caller_module
from hgvs_caller import TooManyVariationsException, find_variant_string, find_variants

# --------------------------------------------------------------------------- grammar

_POSITIVE = r"[1-9]\d*"
_POS = {"cr": r"-?[1-9]\d*", "gmno": _POSITIVE}
_AA = r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|Ser|Thr|Trp|Tyr|Val|Xaa)"
_RES = rf"(?:{_AA}|Ter){_POSITIVE}"


def _nucleotide_edit(kind: str) -> str:
    pos = _POS["cr"] if kind in "cr" else _POS["gmno"]
    base = "[acgun]" if kind == "r" else "[ACGTN]"
    span = rf"{pos}(?:_{pos})?"
    inserted = rf"(?:{base}+|{pos}_{pos}(?:inv)?)"
    return rf"(?:{pos}{base}>{base}|{span}(?:del|dup|inv)|{pos}_{pos}ins{inserted}|{span}delins{inserted})"


def _allele(edit: str) -> str:
    return rf"(?:=|{edit}|\[{edit}(?:;{edit})+\])"


NUCLEOTIDE_GRAMMAR = {kind: re.compile(_allele(_nucleotide_edit(kind))) for kind in "gcrnmo"}

_P_EDIT = (
    rf"(?:{_RES}(?:_{_RES})?="  # synonymous
    rf"|{_RES}(?:{_AA}|Ter)"  # substitution, nonsense
    rf"|{_RES}(?:_{_RES})?(?:del|dup)"
    rf"|{_RES}_{_RES}ins(?:{_AA}+(?:Ter)?|Ter|{_RES}_{_RES})"
    rf"|{_RES}(?:_{_RES})?delins(?:{_AA}+(?:Ter)?|Ter|{_RES}_{_RES})"
    rf"|Ter{_POSITIVE}{_AA}extTer(?:{_POSITIVE}|\?))"
)
PROTEIN_GRAMMAR = re.compile(rf"(?:\?|Met1\?|{_RES}{_AA}fsTer(?:{_POSITIVE}|\?)|{_RES}fs|{_allele(_P_EDIT)})")


def assert_well_formed(description: str, prefix: str) -> None:
    """the description is the prefix followed by a description in HGVS grammar"""
    assert description.startswith(prefix), description
    kind, body = prefix[-2], description[len(prefix) :]
    grammar = PROTEIN_GRAMMAR if kind == "p" else NUCLEOTIDE_GRAMMAR[kind]
    assert grammar.fullmatch(body), description


def check_nucleotide_positions(description: str, prefix: str, ref_length: int, offset: int) -> None:
    """ranges are ordered, insertions are between adjacent positions, changes are 5' to 3'
    and do not overlap, and everything is within the reference"""
    kind, body = prefix[-2], description[len(prefix) :]
    if body == "=":
        return
    skip_zero = kind in "cr" and offset + 1 <= 0

    def index(pos: str) -> int:
        p = int(pos)
        return p - offset - 1 - (1 if skip_zero and p > 0 else 0)

    previous_end = -1
    for edit in body.strip("[]").split(";"):
        m = re.match(r"(-?\d+)(?:_(-?\d+))?(>|.*?)(del|dup|inv|ins|delins|>)?", edit)
        assert m
        positions = [int(p) for p in re.findall(r"-?\d+", re.split(r"ins", edit)[0])]
        first, last = index(str(positions[0])), index(str(positions[-1]))
        assert 0 <= first <= last < ref_length, description
        assert len(positions) == 1 or first < last, description
        if re.fullmatch(r"-?\d+_-?\d+ins.*", edit):  # between two positions: may follow a change
            assert last == first + 1 and first >= previous_end, description
            previous_end = first
        else:
            assert first > previous_end, description
            previous_end = last


# --------------------------------------------------------------------------- strategies

PREFIXES = ["g.", "c.", "r.", "n.", "m.", "o.", "NM_004006.2:c.", "NC_000023.11:g."]
dna = st.text(alphabet="ACGTACGTACGTN", min_size=1, max_size=50)
repetitive_dna = st.text(alphabet="AT", min_size=1, max_size=40)
any_dna = st.one_of(dna, repetitive_dna)


@st.composite
def mutated(draw, ref: str) -> str:
    """ref with 1-4 random substitutions, deletions, insertions, duplications or inversions"""
    seq = ref
    for _ in range(draw(st.integers(1, 4))):
        i = draw(st.integers(0, len(seq)))
        n = draw(st.integers(1, 6))
        op = draw(st.sampled_from(["sub", "del", "ins", "dup", "inv"]))
        if op == "sub" and i < len(seq):
            seq = seq[:i] + draw(st.sampled_from("ACGTN")) + seq[i + 1 :]
        elif op == "del":
            seq = seq[:i] + seq[i + n :]
        elif op == "ins":
            seq = seq[:i] + draw(st.text(alphabet="ACGT", min_size=1, max_size=n)) + seq[i:]
        elif op == "dup":
            seq = seq[:i] + seq[max(0, i - n) : i] + seq[i:]
        else:
            seq = seq[:i] + revcomp(seq[i : i + n]) + seq[i + n :]
    return seq or "A"


@st.composite
def sequence_pairs(draw):
    """(ref, alt): alt is a mutated copy of ref, or an unrelated sequence"""
    ref = draw(any_dna)
    alt = draw(st.one_of(mutated(ref), any_dna))
    return ref, alt


def offsets(kind: str):
    return st.integers(-200, 10**7) if kind in "cr" else st.integers(0, 10**7)


SENSE = [codon for codon, aa in caller_module.CODON_TABLE.items() if aa != "*"]


@st.composite
def coding_pairs(draw):
    """(ref, alt): an open reading frame and a mutated copy"""
    ref = "ATG" + "".join(draw(st.lists(st.sampled_from(SENSE), min_size=2, max_size=15))) + "TAA"
    return ref, draw(mutated(ref))


# --------------------------------------------------------------------------- robustness


@given(
    prefix=st.one_of(st.sampled_from(PREFIXES + ["p."]), st.text(max_size=6)),
    ref=st.one_of(any_dna, st.text(max_size=20)),
    alt=st.one_of(any_dna, st.text(max_size=20)),
    offset=st.integers(-(10**6), 10**6),
    fs_short=st.booleans(),
    merge=st.sampled_from(["codon", "always", "never", "", "x"]),
    min_copy=st.one_of(st.none(), st.integers(1, 20)),
)
def test_any_input_gives_description_or_value_error(prefix, ref, alt, offset, fs_short, merge, min_copy):
    try:
        description = find_variant_string(
            prefix, ref, alt, offset=offset, fs_short=fs_short, merge_one_nt=merge, min_copy_length=min_copy
        )
    except ValueError:
        return
    assert_well_formed(description, prefix)


# --------------------------------------------------------------------------- nucleotides


@given(data=st.data(), pair=sequence_pairs(), prefix=st.sampled_from(PREFIXES))
def test_nucleotide_descriptions(data, pair, prefix):
    """well-formed, positions consistent, and reproduces the variant, with every option"""
    ref, alt = pair
    kind = prefix[-2]
    offset = data.draw(offsets(kind), label="offset")
    merge = data.draw(st.sampled_from(["codon", "always", "never"]), label="merge_one_nt")
    min_copy = data.draw(st.one_of(st.none(), st.integers(2, 12)), label="min_copy_length")
    description = find_variant_string(prefix, ref, alt, offset=offset, merge_one_nt=merge, min_copy_length=min_copy)
    assert_well_formed(description, prefix)
    check_nucleotide_positions(description, prefix, len(ref), offset)
    assert apply_nucleotide(description[len(prefix) - 2 :], ref, offset) == alt


@given(ref=any_dna, data=st.data())
def test_single_deletion_and_duplication(ref, data):
    """a deletion or tandem duplication is described as such, at its most 3' position"""
    i = data.draw(st.integers(0, len(ref) - 1), label="start")
    n = data.draw(st.integers(1, len(ref) - i), label="length")
    for alt, op in ((ref[:i] + ref[i + n :], "del"), (ref[: i + n] + ref[i : i + n] + ref[i + n :], "dup")):
        if not alt:
            continue
        description = find_variant_string("g.", ref, alt)
        m = re.fullmatch(rf"g\.(\d+)(?:_(\d+))?{op}", description)
        assert m, (ref, alt, description)
        start, end = int(m.group(1)) - 1, int(m.group(2) or m.group(1)) - 1
        assert end - start + 1 == n, (ref, alt, description)
        assert end + 1 >= len(ref) or ref[start] != ref[end + 1], (ref, alt, description)


# --------------------------------------------------------------------------- metamorphic


@given(pair=sequence_pairs(), prefix=st.sampled_from(["g.", "c.", "r."]))
def test_case_and_rna_input_do_not_matter(pair, prefix):
    ref, alt = pair
    expected = find_variant_string(prefix, ref, alt)
    assert find_variant_string(prefix, ref.lower(), alt.lower()) == expected
    assert find_variant_string(prefix, ref.replace("T", "U"), alt.replace("T", "U").lower()) == expected


@given(pair=sequence_pairs())
def test_minus_strand_is_reverse_complement(pair):
    ref, alt = pair
    expected = find_variant_string("g.", revcomp(ref), revcomp(alt))
    assert find_variant_string("g.", ref, alt, minus_strand=True) == expected


@given(pair=sequence_pairs(), offset=st.integers(0, 10**7))
def test_offset_only_shifts_coordinates(pair, offset):
    ref, alt = pair
    unshifted = find_variant_string("g.", ref, alt)
    shifted = re.sub(r"\d+", lambda m: str(int(m.group()) + offset), unshifted)
    assert find_variant_string("g.", ref, alt, offset=offset) == shifted


@given(ref=any_dna, prefix=st.sampled_from(PREFIXES + ["p."]))
def test_no_change(ref, prefix):
    assert find_variant_string(prefix, ref, ref) == prefix + "="


@given(pair=sequence_pairs(), limit=st.integers(0, 5))
def test_max_mutations_agrees_with_find_variants(pair, limit):
    ref, alt = pair
    count = len(find_variants("g.", ref, alt))
    try:
        find_variant_string("g.", ref, alt, max_mutations=limit)
        assert count <= limit
    except TooManyVariationsException:
        assert count > limit


# --------------------------------------------------------------------------- protein


@given(pair=coding_pairs(), fs_short=st.booleans())
def test_protein_descriptions(pair, fs_short):
    """well-formed and consistent with the translated sequences"""
    ref, alt = pair
    description = find_variant_string("p.", ref, alt, fs_short=fs_short)
    assert_well_formed(description, "p.")
    if not fs_short:
        check_protein(ref, alt, description)


@given(pair=sequence_pairs(), offset=st.integers(0, 10**6))
def test_protein_arbitrary_sequences(pair, offset):
    """not open reading frames, any frame: still well-formed"""
    ref, alt = pair
    assume(len(ref) >= 3)
    assert_well_formed(find_variant_string("p.", ref, alt, offset=offset), "p.")
