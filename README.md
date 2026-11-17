# HGVS Extractor

[![tests](https://github.com/nickzoic/hgvs_extractor/actions/workflows/tests.yml/badge.svg)](https://github.com/nickzoic/hgvs_extractor/actions/workflows/tests.yml)

Extract [HGVS](https://hgvs-nomenclature.org/stable/) variant descriptions from two sequences:
give it a reference and a variant nucleotide sequence, get back a description such as
`c.20_23dup`, `r.775delinsga` or `p.Glu5ValfsTer5`.

```python
>>> from hgvs_extractor import extract_variant_string
>>> extract_variant_string("c.", "AGAAGTAGAGG", "AGAAGTAGATAGAGG", offset=14)
'c.20_23dup'

```

- Genomic, coding, non-coding, mitochondrial, circular DNA and RNA descriptions
  (`g.`, `c.`, `n.`, `m.`, `o.`, `r.`) and the predicted protein consequence (`p.`).
- Follows HGVS Nomenclature 21.1, and produces exactly the descriptions given by the
  specification for all of its examples that show both the reference and the variant sequence.
- No required dependencies; Python 3.10 or later.

## Installation

```sh
pip install hgvs_extractor           # pure Python
pip install "hgvs_extractor[fast]"   # also installs rapidfuzz: faster for reads with several changes
```

From a checkout of https://github.com/nickzoic/hgvs_extractor: `pip install .`
(or `pip install -e ".[dev]"` for development).

## Usage

```python
>>> from hgvs_extractor import extract_variant_string, extract_variants

```

`offset` places the reference on the coordinate system: the first base of `ref_seq` gets
coordinate `offset + 1`. In `c.` and `r.` numbering there is no position 0, so a reference
starting three bases before the start codon (`c.-3`) has `offset=-4`. Other coordinates start
at 1, so for them a negative `offset` raises `ValueError`.

```python
>>> extract_variant_string("g.", "TGTGATTCT", "TGTGTTCT", offset=33344586)
'g.33344591del'
>>> extract_variant_string("c.", "CCCATG", "CCGATG", offset=-4)
'c.-1C>G'
>>> extract_variant_string("NM_004006.2:c.", "AGAAGTAGAGG", "AGAAGTTAGAGG", offset=14)
'NM_004006.2:c.20dup'

```

The prefix may include the reference sequence; its last two characters select the coordinate
type. Input may be upper or lower case, with `T` or `U`; `r.` descriptions are written in
lower case with `u`:

```python
>>> extract_variant_string("r.", "aggcucauu", "aggcgacauu", offset=770)
'r.775delinsga'

```

Several changes are returned as an allele, no change as `=`:

```python
>>> extract_variant_string("g.", "GATTACA", "GTTTAGA")
'g.[2A>T;6C>G]'
>>> extract_variant_string("g.", "GATTACA", "GATTACA")
'g.='
>>> extract_variants("g.", "GATTACA", "GTTTAGA")
['2A>T', '6C>G']

```

### Protein consequences

With a `p.` prefix the sequences are coding DNA: the reading frame starts at the first codon
boundary at or after the start of `ref_seq` (with `offset=0`, `ref_seq` starts with codon 1).

```python
>>> extract_variant_string("p.", "ATGGATGCATACGAGATGAGG", "ATGGATGCATACGTGCATACGAGATGAGG")
'p.Glu5ValfsTer5'
>>> extract_variant_string("p.", "ATGGATGCATACGAGATGAGG", "ATGGATGCATACGTGCATACGAGATGAGG", fs_short=True)
'p.Glu5fs'
>>> extract_variant_string("p.", "ATGGTTGGTTCATAA", "ATGGTTTGATCATAA")
'p.Gly3Ter'
>>> extract_variant_string("p.", "ATGAAAATGCATTGGTAA", "ATGAAAATGGGGTAGCATTGGTAA")
'p.Met3_His4insGlyTer'
>>> extract_variant_string("p.", "ATGGTTGGTTCATAA", "ATGGTTGGTTCACAAGGTTGAAAA")
'p.Ter5GlnextTer2'
>>> extract_variant_string("p.", "TGC", "TGT", offset=561)
'p.Cys188='

```

When the consequence lies beyond the end of the given sequence (for example a stop codon
deleted at the end of the read), the result is `p.?`, "effect unknown".

### API

`extract_variant_string(prefix, ref_seq, var_seq, max_mutations=None, offset=0, minus_strand=False, **options)`
returns the description as a string. Apart from its name, the signature is that of
[countess-variant-caller](https://github.com/CountESS-Project/countess-variant-caller)'s
`find_variant_string`.

- `max_mutations`: raise `TooManyVariationsException` (a `ValueError`) if more changes are found.
- `minus_strand`: reverse-complement both sequences first.
- `fs_short`: describe frameshifts in the short form, `p.Glu5fs` instead of `p.Glu5ValfsTer5`.
- `merge_one_nt`: when two changes one nucleotide apart become one delins.
  `"codon"` (default) follows the current HGVS rule: only when they affect one amino acid,
  which is known for `c.` and `r.` coordinates. `"always"` follows the rule the HGVS SVD-WG has
  proposed (changes less than two nucleotides apart are a delins); `"never"` never merges.
- `min_copy_length`: inserted sequences at least this long are described by reference to a copy
  in the reference sequence, e.g. `g.17_18ins5_16` or `p.His7_Gln8insGly4_Ser6`.
  Defaults: 10 nucleotides, 3 amino acids.

`extract_variants(...)` takes the same arguments (without `max_mutations`) and returns the list of
descriptions without the prefix. Invalid sequences or prefixes raise `ValueError`.
`reverse_complement(dna)` and `translate(dna, to_stop=True)` are exported as well.

## How descriptions are chosen

- **3' rule**: deletions, duplications and insertions are placed as far 3' as possible,
  for DNA, RNA and protein.
- An insertion that copies the sequence directly 5' of it is a **duplication**.
- Changes of two or more consecutive residues are one **deletion-insertion**; a delins whose
  inserted sequence is the reverse complement of the deleted sequence is an **inversion**.
- An **inversion** inside a changed region is found too, and stays a separate change even
  next to another one, as in the spec's `g.[776788_93191545inv;93191546T>C]`:

  ```python
  >>> extract_variant_string("g.", "AAGAAGATGGCTCCTG", "AAGAAGCCATGTCCTG")
  'g.[7_10inv;11C>G]'

  ```
- Changes separated by one or more residues are described **individually**, as an allele,
  except as set by `merge_one_nt`.
- Insertions before the first or after the last base of the reference have only one flanking
  residue, so they are described as a delins of that residue (or as a duplication).
- **Protein**: descriptions compare the protein sequences. A stop as the first changed residue
  is a nonsense variant (`Gly3Ter`); amino acids after a new stop are not listed; an insertion
  containing a stop is an insertion (`insGlyTer`), not a delins of the rest of the protein.
  Where an in-frame change could be read either way, the reading with the shorter description
  is used (stops inserted after the stop codon change nothing). A change **of the stop codon**
  is an extension (`Ter5GlnextTer2`); residues added before an unchanged stop codon are an
  insertion or duplication (`Ser4_Ter5insGly`). Changing `Met1` gives `p.Met1?`.

Not supported: intronic positions (`c.88+1G>T`), repeat notation (`c.89_118AGC[13]`),
uncertain positions, and descriptions on more than one reference sequence.

### Method

The common prefix and then the common suffix of the two sequences are removed, which on its
own puts a single change at its most 3' position. The remaining middle is aligned with the
fewest edits, preferring substitutions on ties so that adjacent changes stay together, in a
band around the diagonal that is widened until the result is provably the same as a full
alignment. Unchanged stretches between changes are skipped: once the band's values settle
along a matching diagonal, every following row is the same while the match lasts, so those
rows are not computed (the alignment is the same as computing them).

Clusters of nearby changes (up to 8 bases apart) are then aligned again with an affine gap
cost, where starting a change costs one more than continuing it: edits plus changes, then
edits. So a deletion is not split around a matching base, and a stretch of many small changes
can become one larger change. Each cluster is also searched for an inversion (counted as two
edits, one per breakpoint), which is used where it costs less by the same measure, e.g.
`[26_29inv;30C>G]` rather than `[25_26insCC;29_30del]`. The changes are then merged per the
rules above, shifted 3', and classified. If `rapidfuzz` is installed it is used to compute the
edit distance that sizes the band; the results are the same with or without it.

Typical cost per read on a 300 bp amplicon:

| changes in the read     | `g.`, with rapidfuzz | `g.`, pure Python | `p.`     |
|-------------------------|----------------------|-------------------|----------|
| one                     | 14 µs                | 19 µs             | 9-11 µs  |
| three, far apart        | 131 µs               | 201 µs            | 17-18 µs |

A 10 kb sequence with five changes far apart takes about 0.7 ms with rapidfuzz, 1 ms without.

## Tests

```sh
pip install -e ".[test]"
pytest                                    # all tests
pytest --cov                              # with a coverage report (pytest-cov)
pytest tests/test_fuzz.py --hypothesis-profile=fuzz    # a longer fuzz run
```

Tests run from a checkout without installing too, and on GitHub Actions on every push:
Python 3.10 to 3.14 with and without rapidfuzz, each requiring 100% line and branch
coverage; a longer fuzz run; linting; and a test of the built package. They include:

- `tests/test_hgvs_spec_examples.py`: every example on the HGVS Nomenclature site
  (https://hgvs-nomenclature.org/stable/, version 21.1.5) whose text gives both the reference
  and the variant sequence, from the Examples sections and from the Notes and Q&A. Each is
  checked against the spec's description and by applying the description to the reference.
- `tests/test_api.py`: the API, options and protein edge cases.
- `tests/test_properties.py`: seeded randomized tests: the described change reproduces the
  variant sequence (or protein), the 3' rule holds, a new stop codon as the first changed
  residue is described as nonsense, and the banded alignment equals a full alignment. Set
  `HGVS_EXTRACTOR_PROPERTY_CASES` (default 3000) for more cases.
- `tests/test_fuzz.py`: fuzz tests with [Hypothesis](https://hypothesis.readthedocs.io/),
  which generates inputs and shrinks any failure to a minimal example: any input gives a
  well-formed HGVS description or a `ValueError`; every description matches the HGVS grammar
  and reproduces the variant; unrelated sequence pairs as well as similar ones; and case, T/U,
  strand and offset changes only change the result as they should. Cases they found are kept
  as regression tests in `tests/test_api.py`.
- `tests/test_description_extractor.py`: the test cases of Mutalyzer's
  [description-extractor](https://github.com/mutalyzer/description-extractor), and 500 of its
  results on random sequence pairs (see below); all give the same description.
- the doctests in the package and in this README.

Two examples from the specification are not included: the exon/intron junction examples
(their sequences are not contiguous), and `c.-23A>T` from the protein Extension Q&A, whose
sequence shows a `G>T` change (`CAGGGT` to `CATGGT`), i.e. `c.-23G>T`.

### Comparison with Mutalyzer's description-extractor

`tools/compare_description_extractor.py` runs both on random sequence pairs (it needs
description-extractor built: Python 3.10, SWIG from PyPI, biopython and
`git+https://github.com/mutalyzer/crossmapper-legacy.git`). On 20,000 pairs:

| share | |
|-------|-|
| 74.0% | the same description |
| 17.9% | description-extractor joins separate changes in one delins (`32_41delinsTCCCAT` for `[32_33del;36_38del;43dup]`); HGVS describes changes separated by one or more nucleotides individually |
| 7.0%  | both correct, different alignment: description-extractor minimizes the length of the description, hgvs_extractor edits plus changes (see Method) |
| 0.8%  | description-extractor's description does not give the variant sequence (e.g. an inversion one base out) |
| 0.3%  | description-extractor uses a position outside the reference (`0_1insAATAC`) |

Every description from hgvs_extractor gave the variant sequence. Where both are correct and
description-extractor does not join separate changes, scored by hgvs_extractor's measure (edits plus changes,
then edits): hgvs_extractor's description costs less for 1574 pairs, the same for 235, and more
for 13, all with several inversions close together.

## Benchmark

[BENCHMARK.md](BENCHMARK.md) compares hgvs_extractor with countess-variant-caller, Mutalyzer 2's
description-extractor, Mutalyzer 3's description extractor and mutalyzer-algebra on 7063 cases
(`tools/benchmark.py`). In brief, on 5000 random short sequence pairs:

| tool | correct, valid HGVS | correct, positions outside the reference | wrong | time per pair |
|---|---|---|---|---|
| hgvs_extractor | 100% | 0% | 0% | 164 µs |
| countess-variant-caller | 92.7% | 4.3% | 2.9% | 7 µs |
| Mutalyzer 2 description-extractor | 95.9% | 3.3% | 0.8% | 38 µs |
| Mutalyzer 3 description_extractor | 96.7% | 3.3% | 0% | 125 µs |
| mutalyzer-algebra | 95.7% | 4.3% | 0% | 112 µs |

hgvs_extractor also gave the HGVS specification's description for the most of its examples
(24 of 25, against 19 for Mutalyzer). countess-variant-caller is the fastest; hgvs_extractor
is about as fast as Mutalyzer 2's C++ (faster for one change, and on 10 kb sequences).

## Origin

This code was written almost entirely by [Claude](https://claude.ai/) and I'm not pretending
otherwise.  The prompts I gave it were (approximately):

> extract all the examples from the HGVS specification site at https://hgvs-nomenclature.org/stable/
> which show a sequence variation and its HGVS description and build a python list of tuples
> which can be used to test an implementation.
>
> include only examples where the before and after sequences can be found in the example text.
> output tests as a pytest compatible file.
>
> write a python function which takes two nucleotide sequences and returns an HGVS string representing
> the variant.  The function should pass all the tests above.
>
> make this into a python package with a README.md, pyproject.toml and make the tests run under pytest.

## Licence

CC0 1.0 Universal (public domain dedication): see [LICENSE](LICENSE).
The HGVS Nomenclature text the test examples are taken from is also published under CC0.
