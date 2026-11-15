"""Call HGVS variant descriptions from a reference and a variant nucleotide sequence.

>>> find_variant_string("g.", "TGTGATTCT", "TGTGTTCT", offset=33344586)
'g.33344591del'
>>> find_variant_string("c.", "AGAAGTAGAGG", "AGAAGTAGATAGAGG", offset=14)
'c.20_23dup'
>>> find_variant_string("g.", "GATTACA", "GTTTAGA")
'g.[2A>T;6C>G]'
>>> find_variant_string("r.", "aggcucauu", "aggcgacauu", offset=770)
'r.775delinsga'
>>> find_variant_string("p.", "ATGGATGCATACGAGATGAGG", "ATGGATGCATACGTGCATACGAGATGAGG")
'p.Glu5ValfsTer5'

The signature is compatible with countess_variant_caller.find_variant_string.
Follows HGVS Nomenclature 21.1 (https://hgvs-nomenclature.org/stable/), in particular:

* the 3' rule: deletions, duplications and insertions are placed as far 3' as possible;
* an insertion that is a copy of the sequence directly 5' of it is a duplication;
* changes of two or more consecutive residues are a deletion-insertion (delins);
* changes separated by one or more residues are described individually, except that two
  changes one nucleotide apart affecting one amino acid (c. and r. coordinates) are a delins;
* a delins whose inserted sequence is the reverse complement of the deleted one is an inversion.

Not supported: intronic positions (c.88+1), repeat notation (c.89_118AGC[13]), and
descriptions of uncertain positions.
"""

from typing import Callable, Optional

_edit_distance: Optional[Callable[[str, str], int]]
try:  # optional: exact edit distance in C, used only to size the alignment band
    from rapidfuzz.distance.Levenshtein import distance as _edit_distance
except ImportError:  # pragma: no cover
    _edit_distance = None

__all__ = ["find_variants", "find_variant_string", "TooManyVariationsException", "reverse_complement", "translate"]

# Inserted sequences at least this long are described by reference to a copy elsewhere in
# the reference sequence when one exists (e.g. g.17_18ins5_16, p.His7_Gln8insGly4_Ser6).
MIN_COPY_LENGTH_NT = 10
MIN_COPY_LENGTH_AA = 3

# Above this many cells the middle region is not aligned but described as one change.
MAX_ALIGNMENT_CELLS = 4_000_000

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu", "G": "Gly",
    "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser",
    "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter", "X": "Xaa",
}  # fmt: skip

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}  # fmt: skip # pylint: disable=line-too-long

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class TooManyVariationsException(ValueError):
    """Raised when the number of variations found exceeds max_mutations"""


def reverse_complement(seq: str) -> str:
    """reverse complement of an (uppercase) DNA sequence"""
    return seq.translate(_COMPLEMENT)[::-1]


def translate(dna: str, to_stop: bool = True) -> str:
    """translate DNA to one-letter amino acids, by default stopping after the first stop codon"""
    protein = []
    for i in range(0, len(dna) - 2, 3):
        aa = CODON_TABLE.get(dna[i : i + 3], "X")
        protein.append(aa)
        if aa == "*" and to_stop:
            break
    return "".join(protein)


# --------------------------------------------------------------------------- alignment


def _align(ref: str, alt: str) -> list[str]:
    """Levenshtein alignment as a list of operations in reverse order: "=" match,
    "X" substitution, "D" deletion from ref, "I" insertion from alt.  On ties a
    match/substitution is preferred, so changed residues stay together rather than
    becoming an insertion and a deletion.

    Only a band of diagonals is computed: those between 0 and len(alt) - len(ref), widened
    by `width` on each side.  A path leaving the band needs more than 2 * width insertions
    and deletions, so when the distance found is at most `width` the result is the same as
    the full alignment; otherwise the band is doubled and the alignment repeated."""
    n, m = len(ref), len(alt)
    inf = n + m + 1
    width = max(1, _edit_distance(ref, alt)) if _edit_distance is not None else 2
    while True:
        low = min(0, m - n) - width  # column j of row i is stored at index j - i - low
        cols = max(0, m - n) + width - low + 1
        dist = [[inf] * cols for _ in range(n + 1)]
        row = dist[0]
        for j in range(0, min(m, cols + low - 1) + 1):
            row[j - low] = j
        for i in range(1, n + 1):
            prev, row, r_char = dist[i - 1], dist[i], ref[i - 1]
            for k in range(max(0, i + low) - i - low, min(m, i + low + cols - 1) - i - low + 1):
                j = i + low + k
                if j == 0:
                    row[k] = i
                    continue
                best = prev[k] + (r_char != alt[j - 1])  # diagonal: same k in previous row
                if k + 1 < cols and prev[k + 1] + 1 < best:  # from above: (i-1, j)
                    best = prev[k + 1] + 1
                if k and row[k - 1] + 1 < best:  # from the left: (i, j-1)
                    best = row[k - 1] + 1
                row[k] = best
        if dist[n][m - n - low] <= width or width >= n + m:
            break
        width *= 2

    def get(i: int, j: int) -> int:
        k = j - i - low
        return dist[i][k] if 0 <= k < cols else inf

    ops = []
    i, j = n, m
    while i or j:
        here = get(i, j)
        if i and j and here == get(i - 1, j - 1) + (ref[i - 1] != alt[j - 1]):
            ops.append("=" if ref[i - 1] == alt[j - 1] else "X")
            i, j = i - 1, j - 1
        elif i and here == get(i - 1, j) + 1:
            ops.append("D")
            i -= 1
        else:
            ops.append("I")
            j -= 1
    return ops


def _blocks(ref: str, alt: str, inversions: bool) -> list[list[int]]:
    """Changed regions between ref and alt as [ref_start, ref_end, alt_start, alt_end]
    (python-style half-open ranges), in order.

    The common prefix is removed first, then the common suffix, which on its own places
    a single change at its most 3' position.  Only the remaining middle is aligned."""
    prefix = 0
    limit = min(len(ref), len(alt))
    while prefix < limit and ref[prefix] == alt[prefix]:
        prefix += 1
    suffix = 0
    limit -= prefix
    while suffix < limit and ref[-1 - suffix] == alt[-1 - suffix]:
        suffix += 1
    r_mid = ref[prefix : len(ref) - suffix]
    a_mid = alt[prefix : len(alt) - suffix]
    single = [[prefix, prefix + len(r_mid), prefix, prefix + len(a_mid)]]

    if not r_mid and not a_mid:
        return []
    if not r_mid or not a_mid or len(r_mid) * len(a_mid) > MAX_ALIGNMENT_CELLS:
        return single
    if inversions and len(r_mid) == len(a_mid) > 1 and a_mid == reverse_complement(r_mid):
        return single

    ops = _align(r_mid, a_mid)
    blocks: list[list[int]] = []
    i, j = prefix, prefix
    for op in reversed(ops):
        if op == "=":
            i, j = i + 1, j + 1
            continue
        if not blocks or blocks[-1][1] != i or blocks[-1][3] != j:
            blocks.append([i, i, j, j])
        if op in "XD":
            i += 1
        if op in "XI":
            j += 1
        blocks[-1][1], blocks[-1][3] = i, j
    return blocks


def _shift_3prime(ref: str, alt: str, block: list[int], limit: int) -> None:
    """Move a pure deletion or insertion 3' while the result is the same; `limit` is the
    ref index it must stay at least one residue before (start of the next change)."""
    r_start, r_end, a_start, a_end = block
    if r_end > r_start and a_end == a_start:  # deletion
        while r_end < limit - 1 and ref[r_start] == ref[r_end]:
            r_start, r_end, a_start, a_end = r_start + 1, r_end + 1, a_start + 1, a_end + 1
    elif r_end == r_start and a_end > a_start:  # insertion
        while r_start < limit - 1 and alt[a_start] == ref[r_start]:
            r_start, r_end, a_start, a_end = r_start + 1, r_end + 1, a_start + 1, a_end + 1
    block[:] = [r_start, r_end, a_start, a_end]


def _normalize(ref: str, alt: str, blocks: list[list[int]]) -> list[list[int]]:
    """3'-shift each change, then give insertions at either end of ref a flanking residue
    (they cannot be described as an insertion between two residues): the block is widened
    to include the first / last residue and merged with any change it then touches.
    Insertions at the 3' end that are duplications are left alone."""
    limit = len(ref) + 1
    for block in reversed(blocks):
        _shift_3prime(ref, alt, block, limit)
        limit = block[0]
    for k, block in enumerate(blocks):
        if not block:  # already merged into its neighbour
            continue
        r_start, r_end, a_start, a_end = block
        if r_start != r_end:
            continue
        n = a_end - a_start
        if r_start == 0:
            block[:] = [0, 1, a_start, a_end + 1]
        elif r_start == len(ref) and not (r_start >= n and ref[r_start - n :] == alt[a_start:a_end]):
            block[:] = [r_start - 1, r_start, a_start - 1, a_end]
        else:
            continue
        for other in (k - 1, k + 1):
            if 0 <= other < len(blocks) and (blocks[other][1] == block[0] or blocks[other][0] == block[1]):
                first, second = sorted((block, blocks[other]))
                block[:] = [first[0], second[1], first[2], second[3]]
                blocks[other] = []
    return [block for block in blocks if block]


# --------------------------------------------------------------------------- nucleotides


class _Nucleotide:
    """describes changes between two nucleotide sequences"""

    def __init__(self, kind: str, ref: str, alt: str, offset: int, merge_one_nt: str, min_copy: int):
        self.kind, self.ref, self.alt = kind, ref, alt
        self.offset, self.merge_one_nt, self.min_copy = offset, merge_one_nt, min_copy

    def coord(self, index: int) -> int:
        """HGVS coordinate of 0-based ref index; c. and r. numbering has no position 0"""
        start = self.offset + 1
        pos = start + index
        if self.kind in "cr" and start <= 0 <= pos:
            pos += 1
        return pos

    def fmt(self, seq: str) -> str:
        """sequence as written in the description (RNA: lowercase with u)"""
        return seq.lower().replace("t", "u") if self.kind == "r" else seq

    def span(self, start: int, end: int) -> str:
        """coordinates of ref[start:end]"""
        if end - start == 1:
            return str(self.coord(start))
        return f"{self.coord(start)}_{self.coord(end - 1)}"

    def inserted(self, seq: str) -> str:
        """an inserted sequence, by reference to a copy in ref if long enough"""
        if len(seq) >= self.min_copy:
            idx = self.ref.rfind(seq)
            if idx >= 0:
                return self.span(idx, idx + len(seq))
            idx = self.ref.rfind(reverse_complement(seq))
            if idx >= 0:
                return self.span(idx, idx + len(seq)) + "inv"
        return self.fmt(seq)

    def mergeable(self, left: list[int], right: list[int]) -> bool:
        """two changes one nucleotide apart: one delins if they affect one amino acid"""
        if right[0] - left[1] != 1 or self.merge_one_nt == "never":
            return False
        if self.merge_one_nt == "always":
            return True
        if self.kind not in "cr":
            return False
        first, last = self.coord(left[0]), self.coord(max(right[1] - 1, right[0]))
        return first > 0 and last > 0 and (first - 1) // 3 == (last - 1) // 3

    def describe(self, block: list[int], prev_end: int = 0) -> str:
        """description of one change; prev_end is where the previous change ends in ref
        (a duplicated stretch must not overlap it)"""
        r_start, r_end, a_start, a_end = block
        ref, deleted, inserted = self.ref, self.ref[r_start:r_end], self.alt[a_start:a_end]
        if not deleted:  # insertion or duplication between r_start - 1 and r_start
            n = len(inserted)
            if r_start - n >= prev_end and ref[r_start - n : r_start] == inserted:
                return f"{self.span(r_start - n, r_start)}dup"
            return f"{self.coord(r_start - 1)}_{self.coord(r_start)}ins{self.inserted(inserted)}"
        if not inserted:
            return f"{self.span(r_start, r_end)}del"
        if len(deleted) == len(inserted) == 1:
            return f"{self.coord(r_start)}{self.fmt(deleted)}>{self.fmt(inserted)}"
        if len(deleted) == len(inserted) > 1 and inserted == reverse_complement(deleted):
            return f"{self.span(r_start, r_end)}inv"
        return f"{self.span(r_start, r_end)}delins{self.inserted(inserted)}"

    def variants(self) -> list[str]:
        """descriptions of all changes, 5' to 3'"""
        blocks: list[list[int]] = []
        for block in _blocks(self.ref, self.alt, inversions=True):
            if blocks and self.mergeable(blocks[-1], block):
                blocks[-1][1], blocks[-1][3] = block[1], block[3]
            else:
                blocks.append(block)
        blocks = _normalize(self.ref, self.alt, blocks)
        return [self.describe(block, blocks[k - 1][1] if k else 0) for k, block in enumerate(blocks)]


# --------------------------------------------------------------------------- protein


class _Protein:  # pylint: disable=too-many-instance-attributes
    """describes the protein consequence of changes between two coding sequences"""

    def __init__(self, ref_dna: str, alt_dna: str, offset: int, fs_short: bool, min_copy: int):
        frame = (3 - offset) % 3
        self.ref_dna, self.alt_dna = ref_dna[frame:], alt_dna[frame:]
        self.aa_offset = (offset + 2) // 3
        self.ref, self.alt = translate(self.ref_dna), translate(self.alt_dna)
        self.frameshift = (len(alt_dna) - len(ref_dna)) % 3 != 0
        self.fs_short, self.min_copy = fs_short, min_copy

    def res(self, index: int) -> str:
        """residue at 0-based index of the reference protein, e.g. 'Gly3'"""
        return f"{AA3[self.ref[index]]}{index + self.aa_offset + 1}"

    def span(self, start: int, end: int) -> str:
        """residues ref[start:end], e.g. 'Gly4' or 'Gly4_Ser6'"""
        return self.res(start) if end - start == 1 else f"{self.res(start)}_{self.res(end - 1)}"

    @staticmethod
    def aa(seq: str) -> str:
        """one-letter to three-letter amino acid codes"""
        return "".join(AA3[x] for x in seq)

    def inserted(self, seq: str) -> str:
        """inserted residues, by reference to a copy in ref if long enough"""
        if len(seq) >= self.min_copy and "*" not in seq:
            idx = self.ref.rfind(seq)
            if idx >= 0:
                return self.span(idx, idx + len(seq))
        return self.aa(seq)

    def synonymous(self) -> list[str]:
        """changed codons when the protein is unchanged, e.g. ['Cys188=']"""
        changed = [
            i
            for i in range(len(self.ref))
            if self.ref_dna[3 * i : 3 * i + 3] != self.alt_dna[3 * i : 3 * i + 3] and len(self.alt_dna) >= 3 * i + 3
        ]
        runs: list[list[int]] = []
        for i in changed:
            if runs and runs[-1][1] == i:
                runs[-1][1] = i + 1
            else:
                runs.append([i, i + 1])
        return [f"{self.span(start, end)}=" for start, end in runs]

    def extension(self, i: int) -> str:
        """the reference stop codon at index i is changed"""
        new_stop = self.alt.find("*", i + 1)
        tail = f"Ter{new_stop - i}" if new_stop >= 0 else "Ter?"
        return f"{self.res(i)}{AA3[self.alt[i]]}ext{tail}"

    def describe(self, block: list[int], prev_end: int = 0) -> tuple[str, bool]:
        """description of one in-frame change, and whether translation stops there"""
        r_start, r_end, a_start, a_end = block
        ref, deleted, inserted = self.ref, self.ref[r_start:r_end], self.alt[a_start:a_end]
        if deleted == "*" and r_start == len(ref) - 1 and inserted:  # stop codon changed: extension
            tail = f"Ter{inserted.index('*')}" if "*" in inserted else "Ter?"
            return f"{self.res(r_start)}{AA3[inserted[0]]}ext{tail}", True
        stops = "*" in inserted
        if stops:  # amino acids after a new stop are not listed
            inserted = inserted[: inserted.index("*") + 1]
            if inserted == "*" and deleted:  # nonsense
                return f"{self.res(r_start)}Ter", True
        if not deleted:
            n = len(inserted)
            if r_start - n >= prev_end and ref[r_start - n : r_start] == inserted:
                return f"{self.span(r_start - n, r_start)}dup", stops
            return f"{self.res(r_start - 1)}_{self.res(r_start)}ins{self.inserted(inserted)}", stops
        if not inserted:
            return f"{self.span(r_start, r_end)}del", stops
        if len(deleted) == len(inserted) == 1:
            return f"{self.res(r_start)}{AA3[inserted]}", stops
        return f"{self.span(r_start, r_end)}delins{self.inserted(inserted)}", stops

    def variants(self) -> list[str]:
        """descriptions of the protein changes, N- to C-terminal"""
        ref, alt = self.ref, self.alt
        if ref == alt:
            return self.synonymous()
        i = 0
        while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
            i += 1
        # The change lies beyond the end of one of the given sequences (a translation ran
        # off its end without a stop codon): the effect on the protein is unknown, "p.?".
        if i >= len(alt) or (i >= len(ref) and self.frameshift):
            return ["?"]
        if i == 0 and self.aa_offset == 0 and ref[0] == "M":
            return ["Met1?"]

        # Is the reference stop codon kept (at its position moved by an in-frame change)?
        # Only a change of the stop codon itself is an extension; residues added before an
        # unchanged stop codon are an insertion or duplication.
        through = translate(self.alt_dna, to_stop=False)
        stop_index = -1
        if ref.endswith("*") and not self.frameshift:
            stop_index = len(ref) - 1 + (len(self.alt_dna) - len(self.ref_dna)) // 3
        stop_kept = 0 <= stop_index < len(through) and through[stop_index] == "*"
        if ref[i] == "*" and not stop_kept:
            return [self.extension(i)]
        if alt[i] == "*":  # nonsense: the first changed residue is a stop
            return [f"{self.res(i)}Ter"]
        if self.frameshift:
            if self.fs_short:
                return [f"{self.res(i)}fs"]
            new_stop = alt.find("*", i)
            tail = f"Ter{new_stop - i + 1}" if new_stop >= 0 else "Ter?"
            return [f"{self.res(i)}{AA3[alt[i]]}fs{tail}"]

        # In frame: compare with the variant translated through any new stop codon, up to
        # the reference stop codon, so the alignment follows the DNA change (an inserted stop
        # is an insertion, a replaced codon a substitution).  Residues after a new stop are
        # dropped when describing.
        if stop_kept:
            through = through[: stop_index + 1]
        elif ref.endswith("*") and alt.endswith("*") and len(alt) - 1 < stop_index:
            ref, through = ref[:-1], alt  # stop codon lost, but a new stop comes first
        else:
            through = alt  # stop codon lost (an extension) or no stop codon in ref
        blocks = _normalize(ref, through, _blocks(ref, through, inversions=False))
        self.alt = through
        stop_removed = ref.endswith("*") and any(
            b[0] <= len(ref) - 1 < b[1] and b[3] - b[2] < b[1] - b[0] for b in blocks
        )
        if stop_removed and "*" not in through:  # and the given sequence ends: effect unknown
            return ["?"]
        result = []
        for k, block in enumerate(blocks):
            description, stops = self.describe(block, blocks[k - 1][1] if k else 0)
            result.append(description)
            if stops:
                break
        return result


# --------------------------------------------------------------------------- API


def _clean(seq: str, what: str) -> str:
    seq = seq.strip().upper().replace("U", "T")
    if not seq or seq.strip("ACGTN"):
        raise ValueError(f"Invalid {what} sequence")
    return seq


def find_variants(
    prefix: str,
    ref_seq: str,
    var_seq: str,
    offset: int = 0,
    minus_strand: bool = False,
    *,
    fs_short: bool = False,
    merge_one_nt: str = "codon",
    min_copy_length: Optional[int] = None,
) -> list[str]:
    """HGVS descriptions (without prefix) of the changes from ref_seq to var_seq.

    prefix:          coordinate type, the last two characters are used ("g.", "c.", "r.", "p.", ...)
    ref_seq/var_seq: DNA or RNA sequences (any case, U or T)
    offset:          coordinate of the first base of ref_seq, minus 1 (c./r. skip position 0);
                     for "p." the coding DNA starts at the first codon boundary after it
    minus_strand:    reverse-complement both sequences first
    fs_short:        describe frameshifts in the short form (p.Glu5fs)
    merge_one_nt:    two changes one nucleotide apart become one delins when they affect one
                     amino acid ("codon", current HGVS rule), "always" (proposed rule), "never"
    min_copy_length: shortest inserted sequence described as a copy (default 10 nt / 3 aa)
    """
    if merge_one_nt not in ("codon", "always", "never"):
        raise ValueError("merge_one_nt must be 'codon', 'always' or 'never'")
    ref, alt = _clean(ref_seq, "reference"), _clean(var_seq, "variant")
    if minus_strand:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    kind = prefix[-2:-1] if prefix.endswith(".") else ""
    if kind == "p":
        return _Protein(ref, alt, offset, fs_short, min_copy_length or MIN_COPY_LENGTH_AA).variants()
    if kind not in ("g", "c", "r", "n", "m", "o"):
        raise ValueError(f"Unsupported prefix {prefix!r}")
    return _Nucleotide(kind, ref, alt, offset, merge_one_nt, min_copy_length or MIN_COPY_LENGTH_NT).variants()


def find_variant_string(
    prefix: str,
    ref_seq: str,
    var_seq: str,
    max_mutations: Optional[int] = None,
    offset: int = 0,
    minus_strand: bool = False,
    **options,
) -> str:
    """HGVS description of the change from ref_seq to var_seq, e.g. 'c.20_23dup'.

    Several separate changes are returned as an allele, e.g. 'g.[2A>T;6C>G]', no change as
    'g.='.  For "p." a change whose effect lies beyond the end of the given sequence gives
    'p.?' (effect unknown).  Keyword options are passed to find_variants.

    >>> find_variant_string("c.", "CATGACA", "CATATAACA", offset=2073)
    'c.2077delinsATA'
    >>> find_variant_string("g.", "CATCAGCCT", "CACTGACCT", offset=32361327)
    'g.32361330_32361333inv'
    >>> find_variant_string("g.", "ATCGATCGATCGATCGAGGGTCCC", "ATCGATCGATCGATCGAATCGATCGATCGGGGTCCC")
    'g.17_18ins5_16'
    >>> find_variant_string("c.", "CAGGGT", "CATGGT", offset=-26)
    'c.-23G>T'
    >>> find_variant_string("p.", "ATGGTTGGTTCA", "ATGGTTTAGACA")
    'p.Gly3Ter'
    >>> find_variant_string("p.", "TGC", "TGT", offset=561)
    'p.Cys188='
    >>> find_variant_string("g.", "GATTACA", "GATTACA")
    'g.='
    """
    variants = find_variants(prefix, ref_seq, var_seq, offset, minus_strand, **options)
    if not variants:
        return prefix + "="
    if max_mutations is not None and len(variants) > max_mutations:
        raise TooManyVariationsException(f"Too many variations ({len(variants)}) in {var_seq}")
    if len(variants) == 1:
        return prefix + variants[0]
    return prefix + "[" + ";".join(variants) + "]"
