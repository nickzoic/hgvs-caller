"""Extract HGVS variant descriptions from a reference and a variant nucleotide sequence.

>>> extract_variant_string("g.", "TGTGATTCT", "TGTGTTCT", offset=33344586)
'g.33344591del'
>>> extract_variant_string("c.", "AGAAGTAGAGG", "AGAAGTAGATAGAGG", offset=14)
'c.20_23dup'
>>> extract_variant_string("g.", "GATTACA", "GTTTAGA")
'g.[2A>T;6C>G]'
>>> extract_variant_string("r.", "aggcucauu", "aggcgacauu", offset=770)
'r.775delinsga'
>>> extract_variant_string("p.", "ATGGATGCATACGAGATGAGG", "ATGGATGCATACGTGCATACGAGATGAGG")
'p.Glu5ValfsTer5'

Apart from its name, the signature of extract_variant_string is that of
countess_variant_caller.find_variant_string.
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

__all__ = [
    "extract_variants",
    "extract_variant_string",
    "TooManyVariationsException",
    "reverse_complement",
    "translate",
]

# Inserted sequences at least this long are described by reference to a copy elsewhere in
# the reference sequence when one exists (e.g. g.17_18ins5_16, p.His7_Gln8insGly4_Ser6).
MIN_COPY_LENGTH_NT = 10
MIN_COPY_LENGTH_AA = 3

# Clusters of changes this close together are aligned again, keeping each change in one
# piece, and searched for inversions; up to this long; inversions from this length.
CLUSTER_GAP = 8
MAX_CLUSTER_LENGTH = 200
MIN_INVERSION_LENGTH = 2

# If aligning the changed region needs a band of more cells than this (the sequences are
# too different), it is described as one change instead.
MAX_ALIGNMENT_CELLS = 4_000_000


class _RegionTooDifferent(Exception):
    """the alignment band would exceed MAX_ALIGNMENT_CELLS"""


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


def _align_edits(ref: str, alt: str) -> list[str]:
    """Levenshtein alignment as a list of operations in reverse order: "=" match,
    "X" substitution, "D" deletion from ref, "I" insertion from alt: the fewest edits.
    On ties a match/substitution is preferred, so changed residues stay together rather
    than becoming an insertion and a deletion.  Used for the whole changed region; clusters
    of nearby changes are then aligned again with _align (see _refine).

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
        if (n + 1) * cols > MAX_ALIGNMENT_CELLS:
            raise _RegionTooDifferent
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


def _align(ref: str, alt: str) -> list[str]:  # pylint: disable=too-many-statements
    """Alignment of ref and alt as a list of operations in reverse order: "=" match,
    "X" substitution, "D" deletion from ref, "I" insertion from alt.

    Minimizes edits plus changes (runs of edits), an affine gap cost where starting a change
    costs one more than continuing it, then edits: so that e.g. a deletion is not split in
    two around a matching base, and a stretch of many small changes may be one larger
    change.  Remaining ties prefer a match or substitution, so changed residues stay
    together rather than becoming an insertion and a deletion.  The cost of a path is
    (edits + changes) * k + edits; the two states record whether the last operation was a
    match ("=") or an edit.  Used for clusters of nearby changes (see _refine).

    Only a band of diagonals is computed: those between 0 and len(alt) - len(ref), widened
    by `width` on each side.  A path leaving the band needs more than 2 * width insertions
    and deletions, so when the edits plus changes found are at most 2 * width + 2 the result
    is the same as the full alignment; otherwise the band is doubled and the alignment
    repeated."""
    n, m = len(ref), len(alt)
    k = n + m + 2  # more than any number of edits
    extend, start = k + 1, 2 * k + 1  # the cost of an edit continuing or starting a change
    inf = 3 * (n + m + 1) * k
    width = max(1, _edit_distance(ref, alt)) if _edit_distance is not None else 2
    while True:
        low = min(0, m - n) - width  # column j of row i is stored at index j - i - low
        cols = max(0, m - n) + width - low + 1
        if (n + 1) * cols > MAX_ALIGNMENT_CELLS:
            raise _RegionTooDifferent
        same = [[inf] * cols for _ in range(n + 1)]  # last operation a match (or the start)
        edit = [[inf] * cols for _ in range(n + 1)]  # last operation an edit
        same[0][-low] = 0
        for j in range(1, min(m, cols + low - 1) + 1):
            edit[0][j - low] = start + (j - 1) * extend
        for i in range(1, n + 1):
            s_prev, e_prev, s_row, e_row, r_char = same[i - 1], edit[i - 1], same[i], edit[i], ref[i - 1]
            for c in range(max(0, i + low) - i - low, min(m, i + low + cols - 1) - i - low + 1):
                j = i + low + c
                if j == 0:
                    e_row[c] = start + (i - 1) * extend
                    continue
                if r_char == alt[j - 1]:
                    s_row[c] = min(s_prev[c], e_prev[c])  # same c in the previous row
                    best = inf
                else:
                    best = min(s_prev[c] + start, e_prev[c] + extend)
                if c + 1 < cols:  # from above: (i-1, j), a deletion
                    best = min(best, s_prev[c + 1] + start, e_prev[c + 1] + extend)
                if c:  # from the left: (i, j-1), an insertion
                    best = min(best, s_row[c - 1] + start, e_row[c - 1] + extend)
                e_row[c] = best
        if min(same[n][m - n - low], edit[n][m - n - low]) // k <= 2 * width + 2 or width >= n + m:
            break
        width *= 2

    def get(table: list[list[int]], i: int, j: int) -> int:
        c = j - i - low
        return table[i][c] if 0 <= c < cols else inf

    ops = []
    i, j = n, m
    in_edit = get(edit, n, m) < get(same, n, m)
    while i or j:
        if not in_edit:  # a match
            ops.append("=")
            here = get(same, i, j)
            i, j = i - 1, j - 1
            in_edit = get(same, i, j) != here
            continue
        # an edit: a substitution if it gives this cost, else a deletion, else an insertion;
        # the previous state continues the change if that gives this cost, else starts it
        here = get(edit, i, j)

        def reaches(di: int, dj: int) -> bool:
            return here in (get(edit, i - di, j - dj) + extend, get(same, i - di, j - dj) + start)

        if i and j and ref[i - 1] != alt[j - 1] and reaches(1, 1):
            di, dj, op = 1, 1, "X"
        elif i and reaches(1, 0):
            di, dj, op = 1, 0, "D"
        else:
            di, dj, op = 0, 1, "I"
        in_edit = get(edit, i - di, j - dj) + extend == here
        ops.append(op)
        i, j = i - di, j - dj
    return ops


def _blocks(ref: str, alt: str, inversions: bool, align=_align_edits) -> list[list[int]]:
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
    if not r_mid or not a_mid:
        return single
    if inversions and len(r_mid) == len(a_mid) > 1 and a_mid == reverse_complement(r_mid):
        return single

    try:
        ops = align(r_mid, a_mid)
    except _RegionTooDifferent:
        return single
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


def _is_inversion(ref: str, alt: str, block: list[int]) -> bool:
    """the change replaces two or more residues by their reverse complement"""
    r_start, r_end, a_start, a_end = block
    return r_end - r_start == a_end - a_start > 1 and alt[a_start:a_end] == reverse_complement(ref[r_start:r_end])


def _cost(ref: str, alt: str, blocks: list[list[int]]) -> tuple[int, int]:
    """Edits plus changes (each change also costs one, as an affine gap penalty does), then
    edits.  Edits: two for an inversion (one for each breakpoint), otherwise the most of the
    deleted and inserted residues (a change contains no matching residues)."""
    edits = sum(2 if _is_inversion(ref, alt, b) else max(b[1] - b[0], b[3] - b[2]) for b in blocks)
    return edits + len(blocks), edits


def _offset_blocks(ref: str, alt: str, r_start: int, a_start: int) -> list[list[int]]:
    """_blocks of ref and alt aligned with _align, with positions moved by r_start and a_start"""
    return [[b[0] + r_start, b[1] + r_start, b[2] + a_start, b[3] + a_start] for b in _blocks(ref, alt, True, _align)]


def _refine_cluster(ref: str, alt: str, cluster: list[list[int]]) -> list[list[int]]:
    """The changes of a cluster, aligned again with _align so that no change is split in two;
    or an inversion inside the cluster with the changes either side of it, whichever costs
    less by _cost (e.g. [26_29inv;30C>G] rather than [25_26insCC;29_30del])."""
    r0, r1, a0, a1 = cluster[0][0], cluster[-1][1], cluster[0][2], cluster[-1][3]
    n, m = r1 - r0, a1 - a0
    if max(n, m) > MAX_CLUSTER_LENGTH:
        return cluster
    if len(cluster) > 1:  # the same number of edits, with each change in one piece
        cluster = _offset_blocks(ref[r0:r1], alt[a0:a1], r0, a0)
    if min(n, m) < MIN_INVERSION_LENGTH:
        return cluster
    inverted, region = reverse_complement(ref[r0:r1]), alt[a0:a1]
    # runs[p][q]: length of the common substring of inverted and region ending before p, q
    runs = [[0] * (m + 1) for _ in range(n + 1)]
    for p in range(1, n + 1):
        for q in range(1, m + 1):
            if inverted[p - 1] == region[q - 1]:
                runs[p][q] = runs[p - 1][q - 1] + 1
    # the longest common substrings (maximal: not extending further), and in small regions
    # also shorter ones, nearest first to the diagonals of the changes around them
    longest = sorted(
        (
            (runs[p][q], p, q)
            for p in range(1, n + 1)
            for q in range(1, m + 1)
            if runs[p][q] >= MIN_INVERSION_LENGTH and (p == n or q == m or runs[p + 1][q + 1] == 0)
        ),
        reverse=True,
    )[:20]
    nearest: list[tuple[int, int, int]] = []
    if n * m <= 2500:

        def distance(length: int, p: int, q: int) -> int:
            diagonal = (q - length) - (n - p)  # alt start minus ref start, within the region
            return min(abs(diagonal), abs(diagonal - (m - n)))

        nearest = sorted(
            (
                (length, p, q)
                for p in range(1, n + 1)
                for q in range(1, m + 1)
                for length in range(MIN_INVERSION_LENGTH, runs[p][q] + 1)
            ),
            key=lambda c: (distance(*c), -c[0]),
        )[:20]
    candidates = list(dict.fromkeys(longest + nearest))
    best, best_cost = cluster, _cost(ref, alt, cluster)
    for length, p, q in candidates:
        r_inv, a_inv = r0 + n - p, a0 + q - length  # inverted[p - length:p] is ref[r_inv:r_inv + length]
        # residues at the ends that inverting leaves unchanged are not part of the inversion
        # (a stretch that is its own reverse complement is not changed by inverting it)
        while length >= 2 and ref[r_inv] == alt[a_inv] and ref[r_inv + length - 1] == alt[a_inv + length - 1]:
            r_inv, a_inv, length = r_inv + 1, a_inv + 1, length - 2
        if length < MIN_INVERSION_LENGTH:
            continue
        blocks = (
            _offset_blocks(ref[r0:r_inv], alt[a0:a_inv], r0, a0)
            + [[r_inv, r_inv + length, a_inv, a_inv + length]]
            + _offset_blocks(ref[r_inv + length : r1], alt[a_inv + length : a1], r_inv + length, a_inv + length)
        )
        cost = _cost(ref, alt, blocks)
        if cost < best_cost:
            best, best_cost = blocks, cost
    return best


def _refine(ref: str, alt: str, blocks: list[list[int]]) -> list[list[int]]:
    """the changes, with each cluster of nearby changes aligned again so that no change is
    split in two, and with an inversion where that needs fewer edits"""
    result: list[list[int]] = []
    k = 0
    while k < len(blocks):
        j = k
        while (
            j + 1 < len(blocks)
            and blocks[j + 1][0] - blocks[j][1] <= CLUSTER_GAP
            and blocks[j + 1][2] - blocks[j][3] <= CLUSTER_GAP
        ):
            j += 1
        result += _refine_cluster(ref, alt, blocks[k : j + 1])
        k = j + 1
    return result


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
        # merge with the nearest remaining changes on either side while they touch or overlap
        merged = True
        while merged:
            merged = False
            before = next((j for j in range(k - 1, -1, -1) if blocks[j]), None)
            after = next((j for j in range(k + 1, len(blocks)) if blocks[j]), None)
            for other in (before, after):
                if other is not None and blocks[other][1] >= block[0] and blocks[other][0] <= block[1]:
                    first, second = sorted((block, blocks[other]))
                    block[:] = [first[0], max(first[1], second[1]), first[2], max(first[3], second[3])]
                    blocks[other] = []
                    merged = True
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
        if len(seq) >= max(2, self.min_copy):  # a copy is a range of at least two
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
        for block in _refine(self.ref, self.alt, _blocks(self.ref, self.alt, inversions=True)):
            # an inversion stays separate, even next to another change (as in the spec's
            # g.[776788_93191545inv;93191546T>C])
            if (
                blocks
                and self.mergeable(blocks[-1], block)
                and not _is_inversion(self.ref, self.alt, blocks[-1])
                and not _is_inversion(self.ref, self.alt, block)
            ):
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
        if len(seq) >= max(2, self.min_copy) and "*" not in seq:  # a copy is a range of at least two
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

    def describe(self, block: list[int], alt: str, prev_end: int = 0) -> tuple[str, bool]:
        """description of one in-frame change to alt, and whether translation stops there"""
        r_start, r_end, a_start, a_end = block
        ref, deleted, inserted = self.ref, self.ref[r_start:r_end], alt[a_start:a_end]
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
        if not ref:  # no complete codon in the reference: nothing to compare with
            return [] if self.ref_dna == self.alt_dna else ["?"]
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
        if i < len(ref) and ref[i] == "*" and not stop_kept:
            return [self.extension(i)]
        if alt[i] == "*":  # nonsense: the first changed residue is a stop
            # (a stop after the end of the reference translation: effect unknown)
            return [f"{self.res(i)}Ter"] if i < len(ref) else ["?"]
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
            # Two readings of the variant protein: through any new stop up to the moved
            # reference stop (so an inserted stop is an insertion), or up to its first stop
            # (residues and stops inserted after the stop codon are not translated).
            readings = [through[: stop_index + 1], alt]
        elif ref.endswith("*") and alt.endswith("*") and len(alt) - 1 < stop_index:
            ref, readings = ref[:-1], [alt]  # stop codon lost, but a new stop comes first
        else:
            readings = [alt]  # stop codon lost (an extension) or no stop codon in ref
        # the reading with the shortest description; on a tie, the first
        return min((self.describe_reading(ref, reading) for reading in readings), key=lambda d: d[0])[1]

    def describe_reading(self, ref: str, reading: str) -> tuple[int, list[str]]:
        """descriptions of the in-frame changes from ref to reading, and how many residues
        they list (residues after a new stop are not listed)"""
        blocks = _normalize(ref, reading, _blocks(ref, reading, inversions=False))
        stop_removed = ref.endswith("*") and any(
            b[0] <= len(ref) - 1 < b[1] and b[3] - b[2] < b[1] - b[0] for b in blocks
        )
        if stop_removed and "*" not in reading:  # and the given sequence ends: effect unknown
            return 0, ["?"]
        result, listed = [], 0
        for k, block in enumerate(blocks):
            description, stops = self.describe(block, reading, blocks[k - 1][1] if k else 0)
            result.append(description)
            inserted = reading[block[2] : block[3]]
            listed += block[1] - block[0] + len(inserted[: inserted.index("*") + 1] if "*" in inserted else inserted)
            if stops:
                break
        return listed, result


# --------------------------------------------------------------------------- API


def _clean(seq: str, what: str) -> str:
    seq = seq.strip().upper().replace("U", "T")
    if not seq or seq.strip("ACGTN"):
        raise ValueError(f"Invalid {what} sequence")
    return seq


def extract_variants(
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
    offset:          coordinate of the first base of ref_seq, minus 1 (c./r. skip position 0,
                     other coordinates start at 1, so offset must not be negative);
                     for "p." the coding DNA starts at the first codon boundary after it
    minus_strand:    reverse-complement both sequences first
    fs_short:        describe frameshifts in the short form (p.Glu5fs)
    merge_one_nt:    two changes one nucleotide apart become one delins when they affect one
                     amino acid ("codon", current HGVS rule), "always" (proposed rule), "never"
    min_copy_length: shortest inserted sequence described as a copy (default 10 nt / 3 aa;
                     never fewer than 2)
    """
    if merge_one_nt not in ("codon", "always", "never"):
        raise ValueError("merge_one_nt must be 'codon', 'always' or 'never'")
    ref, alt = _clean(ref_seq, "reference"), _clean(var_seq, "variant")
    if minus_strand:
        ref, alt = reverse_complement(ref), reverse_complement(alt)
    kind = prefix[-2:-1] if prefix.endswith(".") else ""
    if offset < 0 and kind not in ("c", "r"):
        raise ValueError(f"Negative offset: {prefix!r} coordinates start at 1")
    if kind == "p":
        return _Protein(ref, alt, offset, fs_short, min_copy_length or MIN_COPY_LENGTH_AA).variants()
    if kind not in ("g", "c", "r", "n", "m", "o"):
        raise ValueError(f"Unsupported prefix {prefix!r}")
    return _Nucleotide(kind, ref, alt, offset, merge_one_nt, min_copy_length or MIN_COPY_LENGTH_NT).variants()


def extract_variant_string(
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
    'p.?' (effect unknown).  Keyword options are passed to extract_variants.

    >>> extract_variant_string("c.", "CATGACA", "CATATAACA", offset=2073)
    'c.2077delinsATA'
    >>> extract_variant_string("g.", "CATCAGCCT", "CACTGACCT", offset=32361327)
    'g.32361330_32361333inv'
    >>> extract_variant_string("g.", "ATCGATCGATCGATCGAGGGTCCC", "ATCGATCGATCGATCGAATCGATCGATCGGGGTCCC")
    'g.17_18ins5_16'
    >>> extract_variant_string("c.", "CAGGGT", "CATGGT", offset=-26)
    'c.-23G>T'
    >>> extract_variant_string("p.", "ATGGTTGGTTCA", "ATGGTTTAGACA")
    'p.Gly3Ter'
    >>> extract_variant_string("p.", "TGC", "TGT", offset=561)
    'p.Cys188='
    >>> extract_variant_string("g.", "GATTACA", "GATTACA")
    'g.='
    """
    variants = extract_variants(prefix, ref_seq, var_seq, offset, minus_strand, **options)
    if not variants:
        return prefix + "="
    if max_mutations is not None and len(variants) > max_mutations:
        raise TooManyVariationsException(f"Too many variations ({len(variants)}) in {var_seq}")
    if len(variants) == 1:
        return prefix + variants[0]
    return prefix + "[" + ";".join(variants) + "]"
