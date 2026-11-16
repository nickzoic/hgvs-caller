# HGVS Caller

[![tests](https://github.com/nickzoic/hgvs-caller/actions/workflows/tests.yml/badge.svg)](https://github.com/nickzoic/hgvs-caller/actions/workflows/tests.yml)

Call [HGVS](https://hgvs-nomenclature.org/stable/) variant descriptions from two sequences:
give it a reference and a variant nucleotide sequence, get back a description such as
`c.20_23dup`, `r.775delinsga` or `p.Glu5ValfsTer5`.

```python
>>> from hgvs_caller import find_variant_string
>>> find_variant_string("c.", "AGAAGTAGAGG", "AGAAGTAGATAGAGG", offset=14)
'c.20_23dup'

```

- Genomic, coding, non-coding, mitochondrial, circular DNA and RNA descriptions
  (`g.`, `c.`, `n.`, `m.`, `o.`, `r.`) and the predicted protein consequence (`p.`).
- Follows HGVS Nomenclature 21.1, and produces exactly the descriptions given by the
  specification for all of its examples that show both the reference and the variant sequence.
- No required dependencies; Python 3.10 or later.

## Installation

```sh
pip install hgvs_caller           # pure Python
pip install "hgvs_caller[fast]"   # also installs rapidfuzz: faster for reads with several changes
```

From a checkout of https://github.com/nickzoic/hgvs-caller: `pip install .`
(or `pip install -e ".[dev]"` for development).

## Usage

```python
>>> from hgvs_caller import find_variant_string, find_variants

```

`offset` places the reference on the coordinate system: the first base of `ref_seq` gets
coordinate `offset + 1`. In `c.` and `r.` numbering there is no position 0, so a reference
starting three bases before the start codon (`c.-3`) has `offset=-4`. Other coordinates start
at 1, so for them a negative `offset` raises `ValueError`.

```python
>>> find_variant_string("g.", "TGTGATTCT", "TGTGTTCT", offset=33344586)
'g.33344591del'
>>> find_variant_string("c.", "CCCATG", "CCGATG", offset=-4)
'c.-1C>G'
>>> find_variant_string("NM_004006.2:c.", "AGAAGTAGAGG", "AGAAGTTAGAGG", offset=14)
'NM_004006.2:c.20dup'

```

The prefix may include the reference sequence; its last two characters select the coordinate
type. Input may be upper or lower case, with `T` or `U`; `r.` descriptions are written in
lower case with `u`:

```python
>>> find_variant_string("r.", "aggcucauu", "aggcgacauu", offset=770)
'r.775delinsga'

```

Several changes are returned as an allele, no change as `=`:

```python
>>> find_variant_string("g.", "GATTACA", "GTTTAGA")
'g.[2A>T;6C>G]'
>>> find_variant_string("g.", "GATTACA", "GATTACA")
'g.='
>>> find_variants("g.", "GATTACA", "GTTTAGA")
['2A>T', '6C>G']

```

### Protein consequences

With a `p.` prefix the sequences are coding DNA: the reading frame starts at the first codon
boundary at or after the start of `ref_seq` (with `offset=0`, `ref_seq` starts with codon 1).

```python
>>> find_variant_string("p.", "ATGGATGCATACGAGATGAGG", "ATGGATGCATACGTGCATACGAGATGAGG")
'p.Glu5ValfsTer5'
>>> find_variant_string("p.", "ATGGATGCATACGAGATGAGG", "ATGGATGCATACGTGCATACGAGATGAGG", fs_short=True)
'p.Glu5fs'
>>> find_variant_string("p.", "ATGGTTGGTTCATAA", "ATGGTTTGATCATAA")
'p.Gly3Ter'
>>> find_variant_string("p.", "ATGAAAATGCATTGGTAA", "ATGAAAATGGGGTAGCATTGGTAA")
'p.Met3_His4insGlyTer'
>>> find_variant_string("p.", "ATGGTTGGTTCATAA", "ATGGTTGGTTCACAAGGTTGAAAA")
'p.Ter5GlnextTer2'
>>> find_variant_string("p.", "TGC", "TGT", offset=561)
'p.Cys188='

```

When the consequence lies beyond the end of the given sequence (for example a stop codon
deleted at the end of the read), the result is `p.?`, "effect unknown".

### API

`find_variant_string(prefix, ref_seq, var_seq, max_mutations=None, offset=0, minus_strand=False, **options)`
returns the description as a string. The signature is compatible with
[countess-variant-caller](https://github.com/CountESS-Project/countess-variant-caller)'s function
of the same name.

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

`find_variants(...)` takes the same arguments (without `max_mutations`) and returns the list of
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
  >>> find_variant_string("g.", "AAGAAGATGGCTCCTG", "AAGAAGCCATGTCCTG")
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
fewest edits and, among those, the fewest separate changes (so a deletion is not split around
a matching base), preferring substitutions on remaining ties so that adjacent changes stay
together. The alignment is computed in a band around the diagonal, widened until the result
is provably the same as a full alignment. Clusters of nearby changes are searched for an
inversion that describes them more briefly. The changes are then merged per the rules above,
shifted 3', and classified. If `rapidfuzz` is installed it is used to compute the edit
distance that sizes the band; the results are the same with or without it.

Typical cost per read on a 300 bp amplicon:

| changes in the read     | `g.`, with rapidfuzz | `g.`, pure Python | `p.`     |
|-------------------------|----------------------|-------------------|----------|
| one                     | 18 µs                | 18 µs             | 9-12 µs  |
| three, far apart        | 430 µs               | 870 µs            | 14-21 µs |

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
  `HGVS_CALLER_PROPERTY_CASES` (default 3000) for more cases.
- `tests/test_fuzz.py`: fuzz tests with [Hypothesis](https://hypothesis.readthedocs.io/),
  which generates inputs and shrinks any failure to a minimal example: any input gives a
  well-formed HGVS description or a `ValueError`; every description matches the HGVS grammar
  and reproduces the variant; unrelated sequence pairs as well as similar ones; and case, T/U,
  strand and offset changes only change the result as they should. Cases they found are kept
  as regression tests in `tests/test_api.py`.
- `tests/test_description_extractor.py`: the test cases of the Mutalyzer
  [description extractor](https://github.com/mutalyzer/description-extractor), and 500 of its
  results on random sequence pairs (see below); all give the same description.
- the doctests in the package and in this README.

Two examples from the specification are not included: the exon/intron junction examples
(their sequences are not contiguous), and `c.-23A>T` from the protein Extension Q&A, whose
sequence shows a `G>T` change (`CAGGGT` to `CATGGT`), i.e. `c.-23G>T`.

### Comparison with the Mutalyzer description extractor

`tools/compare_description_extractor.py` runs both on random sequence pairs (it needs the
extractor built: Python 3.10, SWIG from PyPI, biopython and
`git+https://github.com/mutalyzer/crossmapper-legacy.git`). On 20,000 pairs:

| share | |
|-------|-|
| 73.4% | the same description |
| 18.0% | the extractor joins separate changes in one delins (`32_41delinsTCCCAT` for `[32_33del;36_38del;43dup]`); HGVS describes changes separated by one or more nucleotides individually |
| 7.5%  | both correct, different alignment: the extractor minimizes the length of the description, which prefers fewer, larger changes (`[5_14del;18_19insACC]`); hgvs_caller minimizes the number of edited bases |
| 0.8%  | the extractor's description does not give the variant sequence (e.g. an inversion one base out) |
| 0.3%  | the extractor uses a position outside the reference (`0_1insAATAC`) |

Every description from hgvs_caller gave the variant sequence.

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
