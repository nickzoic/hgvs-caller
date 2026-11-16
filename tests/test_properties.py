"""Randomized property tests (seeded, so they are reproducible).

For random reference sequences and random changes, check that:

* applying the called description to the reference gives the variant sequence;
* deletions and duplications are at their most 3' position (3' rule);
* protein descriptions reproduce the variant protein, frameshifts give the right first
  changed residue and stop position, and a new stop codon is always described as Ter;
* the banded alignment gives exactly the same result as a full alignment.

Set HGVS_CALLER_PROPERTY_CASES to change the number of random cases per test.
"""

import os
import random
import re

import pytest
from hgvs_apply import apply_nucleotide, check_protein, revcomp

import hgvs_caller.caller as caller_module
from hgvs_caller import find_variant_string

CASES = int(os.environ.get("HGVS_CALLER_PROPERTY_CASES", "3000"))
SENSE_CODONS = [codon for codon, aa in caller_module.CODON_TABLE.items() if aa != "*"]


def mutate(rng: random.Random, seq: str) -> str:
    """apply 1-3 random substitutions, deletions, insertions, duplications, inversions or delins"""
    s = list(seq)
    for _ in range(rng.randint(1, 3)):
        i = rng.randrange(len(s))
        op = rng.choice("sdiDvx")
        if op == "s":
            s[i] = rng.choice("ACGT")
        elif op == "d":
            del s[i : i + rng.randint(1, 4)]
        elif op == "i":
            s[i:i] = rng.choices("ACGT", k=rng.randint(1, 5))
        elif op == "D":
            length = rng.randint(1, 4)
            s[i:i] = s[max(0, i - length) : i]
        elif op == "v":
            length = rng.randint(2, 6)
            s[i : i + length] = list(revcomp("".join(s[i : i + length])))
        else:
            s[i : i + rng.randint(1, 3)] = rng.choices("ACGT", k=rng.randint(1, 4))
        if not s:
            s = list(rng.choice("ACGT"))
    return "".join(s)


def nucleotide_cases(seed: int):
    rng = random.Random(seed)
    for n in range(CASES):
        alphabet = "ACGT" if n % 2 else "AT"  # every other case highly repetitive
        ref = "".join(rng.choices(alphabet, k=rng.randint(8, 40)))
        alt = mutate(rng, ref)
        if alt == ref:
            continue
        kind = rng.choice("gcr")
        offset = 0 if kind == "g" else rng.choice([0, 99, -30])
        merge = rng.choice(["codon", "always", "never"])
        yield kind, ref, alt, offset, merge


@pytest.mark.parametrize("seed", [1, 2])
def test_nucleotide_round_trip(seed):
    """applying the called description to the reference gives the variant"""
    for kind, ref, alt, offset, merge in nucleotide_cases(seed):
        called = find_variant_string(kind + ".", ref, alt, offset=offset, merge_one_nt=merge)
        assert apply_nucleotide(called, ref, offset) == alt, (ref, alt, offset, called)


@pytest.mark.parametrize("seed", [3])
def test_nucleotide_3prime_rule(seed):
    """a single deletion or duplication cannot be moved further 3'"""
    for _, ref, alt, _, _ in nucleotide_cases(seed):
        called = find_variant_string("g.", ref, alt)
        if m := re.fullmatch(r"g\.(\d+)(?:_(\d+))?(del|dup)", called):
            start, end = int(m.group(1)) - 1, int(m.group(2) or m.group(1)) - 1
            assert end + 1 >= len(ref) or ref[start] != ref[end + 1], (ref, alt, called)


def protein_cases(seed: int):
    rng = random.Random(seed)
    for _ in range(CASES):
        ref = "ATG" + "".join(rng.choice(SENSE_CODONS) for _ in range(rng.randint(4, 15))) + "TAA"
        s = list(ref)
        i = rng.randrange(3, len(ref) - 3)
        kind = rng.choice(["inframe", "frameshift", "substitution"])
        if kind == "inframe":
            length = rng.choice([3, 6])
            if rng.random() < 0.5:
                del s[i : i + length]
            elif rng.random() < 0.5:
                s[i:i] = rng.choices("ACGT", k=length)
            else:
                s[i:i] = s[max(3, i - length) : i]
        elif kind == "frameshift":
            if rng.random() < 0.5:
                del s[i : i + rng.choice([1, 2, 4])]
            else:
                s[i:i] = rng.choices("ACGT", k=rng.choice([1, 2, 4]))
        else:
            for _ in range(rng.randint(1, 3)):
                s[rng.randrange(3, len(s))] = rng.choice("ACGT")
        yield ref, "".join(s)


@pytest.mark.parametrize("seed", [4, 5])
def test_protein_consequence(seed):
    """protein descriptions are consistent with the translated sequences"""
    for ref, alt in protein_cases(seed):
        check_protein(ref, alt, find_variant_string("p.", ref, alt))


def test_protein_check_catches_unlabelled_new_stop():
    """the check above rejects a description that hides a new stop codon"""
    with pytest.raises(AssertionError):
        check_protein("ATGGTTGGTTCATAA", "ATGGTTTGATCATAA", "p.Gly3_Ser4del")


def full_alignment(ref: str, alt: str) -> list[str]:
    """unbanded version of caller._align (fewest edits, then fewest separate changes), for comparison"""
    n, m = len(ref), len(alt)
    k_edit, inf = n + m + 2, (n + m + 1) * (n + m + 2)
    same = [[inf] * (m + 1) for _ in range(n + 1)]
    edit = [[inf] * (m + 1) for _ in range(n + 1)]
    same[0][0] = 0
    for j in range(1, m + 1):
        edit[0][j] = j * k_edit + 1
    for i in range(1, n + 1):
        edit[i][0] = i * k_edit + 1
        for j in range(1, m + 1):
            if ref[i - 1] == alt[j - 1]:
                same[i][j] = min(same[i - 1][j - 1], edit[i - 1][j - 1])
                best = inf
            else:
                best = min(same[i - 1][j - 1] + k_edit + 1, edit[i - 1][j - 1] + k_edit)
            best = min(best, same[i - 1][j] + k_edit + 1, edit[i - 1][j] + k_edit)
            edit[i][j] = min(best, same[i][j - 1] + k_edit + 1, edit[i][j - 1] + k_edit)
    ops: list[str] = []
    i, j = n, m
    in_edit = edit[n][m] < same[n][m]
    while i or j:
        if not in_edit:
            ops.append("=")
            here = same[i][j]
            i, j = i - 1, j - 1
            in_edit = same[i][j] != here
            continue
        here = edit[i][j]
        for di, dj, op in ((1, 1, "X"), (1, 0, "D"), (0, 1, "I")):
            if i < di or j < dj or (op == "X" and ref[i - 1] == alt[j - 1]):
                continue
            if edit[i - di][j - dj] + k_edit == here:
                in_edit = True
            elif same[i - di][j - dj] + k_edit + 1 == here:
                in_edit = False
            else:
                continue
            ops.append(op)
            i, j = i - di, j - dj
            break
    return ops


@pytest.mark.parametrize("use_rapidfuzz", [True, False])
def test_banded_alignment_equals_full_alignment(use_rapidfuzz, monkeypatch):
    """the banded alignment gives the same alignment, with or without rapidfuzz"""
    # pylint: disable=protected-access
    if use_rapidfuzz and caller_module._edit_distance is None:
        pytest.skip("rapidfuzz not installed")
    if not use_rapidfuzz:
        monkeypatch.setattr(caller_module, "_edit_distance", None)
    rng = random.Random(6)
    for n in range(CASES):
        ref = "".join(rng.choices("ACGT" if n % 2 else "AT", k=rng.randint(1, 60)))
        alt = mutate(rng, ref)
        assert caller_module._align(ref, alt) == full_alignment(ref, alt), (ref, alt)


def test_round_trip_without_rapidfuzz(monkeypatch):
    """results do not depend on whether rapidfuzz is installed"""
    monkeypatch.setattr(caller_module, "_edit_distance", None)
    for kind, ref, alt, offset, merge in nucleotide_cases(7):
        called = find_variant_string(kind + ".", ref, alt, offset=offset, merge_one_nt=merge)
        assert apply_nucleotide(called, ref, offset) == alt, (ref, alt, offset, called)
