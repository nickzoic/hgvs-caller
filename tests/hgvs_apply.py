"""Apply HGVS descriptions to a reference sequence.

Written independently of hgvs_caller, so tests can check that a called description
reproduces the variant sequence rather than comparing against the caller's own logic.
"""

import re

COMPLEMENT = str.maketrans("ACGT", "TGCA")

AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}  # fmt: skip


def norm(seq: str) -> str:
    """RNA or DNA -> uppercase DNA"""
    return seq.upper().replace("U", "T")


def revcomp(seq: str) -> str:
    """reverse complement of uppercase DNA"""
    return seq[::-1].translate(COMPLEMENT)


def split_allele(description: str) -> tuple[str, list[str]]:
    """'g.[3dup;5del]' -> ('g.', ['3dup', '5del'])"""
    prefix, body = description[:2], description[2:]
    return prefix, body[1:-1].split(";") if body.startswith("[") else [body]


def apply_nucleotide(description: str, ref: str, offset: int = 0) -> str:
    """Apply a nucleotide description ("g.5del", "c.[3dup;5del]", ...) to ref, whose first
    base has coordinate offset + 1 (c. and r. numbering skips position 0); returns DNA."""
    prefix, edits = split_allele(description)
    skip_zero = prefix in ("c.", "r.") and offset + 1 <= 0
    original = seq = norm(ref)

    def idx(pos: str) -> int:
        p = int(pos)
        if skip_zero and p > 0:
            p -= 1
        return p - offset - 1

    def first_position(edit: str) -> int:
        m = re.match(r"-?\d+", edit)
        assert m, f"cannot apply {edit}"
        return idx(m.group())

    # apply 3' to 5' so positions stay valid
    for edit in sorted(edits, key=first_position, reverse=True):
        if m := re.fullmatch(r"(-?\d+)([ACGTUacgtu])>([ACGTUacgtu])", edit):
            i = idx(m.group(1))
            assert seq[i] == norm(m.group(2)), f"{edit}: reference base is {seq[i]}"
            seq = seq[:i] + norm(m.group(3)) + seq[i + 1 :]
            continue
        m = re.fullmatch(r"(-?\d+)(?:_(-?\d+))?(delins|del|dup|inv|ins)(.*)", edit)
        assert m, f"cannot apply {edit}"
        start = idx(m.group(1))
        end = idx(m.group(2)) if m.group(2) else start
        op, arg = m.group(3), m.group(4)
        if ref_range := re.fullmatch(r"(-?\d+)_(-?\d+)(inv)?", arg):  # inserted copy of reference
            arg = original[idx(ref_range.group(1)) : idx(ref_range.group(2)) + 1]
            if ref_range.group(3):
                arg = revcomp(arg)
        arg = norm(arg)
        if op == "del":
            seq = seq[:start] + seq[end + 1 :]
        elif op == "dup":
            seq = seq[: end + 1] + seq[start : end + 1] + seq[end + 1 :]
        elif op == "inv":
            seq = seq[:start] + revcomp(seq[start : end + 1]) + seq[end + 1 :]
        elif op == "delins":
            seq = seq[:start] + arg + seq[end + 1 :]
        else:  # ins between start and end (= start + 1)
            seq = seq[: start + 1] + arg + seq[start + 1 :]
    return seq


def apply_protein(description: str, protein: str) -> str:
    """Apply an in-frame protein description ("p.Trp4del", "p.[Gly3Ala;Ser5del]", ...) to a
    one-letter protein sequence starting at residue 1.  An extension (Ter5GlnextTer2) ends
    the result with "?": the rest of the extended protein is not part of the description."""
    _, edits = split_allele(description)

    def idx(residue: str) -> int:
        return int(re.search(r"\d+", residue).group()) - 1  # type: ignore[union-attr]

    def residues(text: str) -> str:
        if m := re.fullmatch(r"([A-Z][a-z]{2}\d+)_([A-Z][a-z]{2}\d+)", text):  # copy of reference
            return protein[idx(m.group(1)) : idx(m.group(2)) + 1]
        return "".join(AA1[t] for t in re.findall(r"[A-Z][a-z]{2}", text))

    seq = list(protein)
    for edit in sorted(edits, key=idx, reverse=True):
        if m := re.fullmatch(r"Ter(\d+)([A-Z][a-z]{2})extTer(\d+|\?)", edit):
            # stop codon changed: the protein continues; "?" marks the unknown continuation
            start = idx(edit)
            assert protein[start] == "*", f"{edit}: reference residue is {protein[start]}"
            seq[start:] = [AA1[m.group(2)], "?"]
            continue
        m = re.fullmatch(r"([A-Z][a-z]{2}\d+)(?:_([A-Z][a-z]{2}\d+))?(delins|del|dup|ins|=)?(.*)", edit)
        assert m, f"cannot apply {edit}"
        start = idx(m.group(1))
        end = idx(m.group(2)) if m.group(2) else start
        assert AA1[m.group(1)[:3]] == protein[start], f"{edit}: reference residue is {protein[start]}"
        op, arg = m.group(3), m.group(4)
        if op == "del":
            seq[start : end + 1] = []
        elif op == "dup":
            seq[end + 1 : end + 1] = seq[start : end + 1]
        elif op == "ins":
            seq[start + 1 : start + 1] = residues(arg)
        elif op == "delins":
            seq[start : end + 1] = residues(arg)
        elif op is None:  # substitution
            seq[start] = AA1[arg]
    return "".join(seq)
