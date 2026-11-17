# Benchmark: variant extraction tools

Tools that extract an HGVS variant description from a reference and an observed sequence,
compared on the same inputs and checked the same way. (VariantValidator and the biocommons
`hgvs` package validate, map and normalize existing descriptions; they do not extract one from
two sequences, so they are not included.)

| tool | version | method | language |
|---|---|---|---|
| hgvs_extractor | 0.1.0 | fewest-edits banded alignment, clusters realigned with an affine gap cost, inversion search | Python (optional rapidfuzz) |
| [countess-variant-caller](https://github.com/CountESS-Project/countess-variant-caller) | 0.1.1 | rapidfuzz Levenshtein opcodes | Python, rapidfuzz (C++) |
| Mutalyzer 2 [description-extractor](https://github.com/mutalyzer/description-extractor) | git 155b79d | shortest description by weights (Vis et al., Bioinformatics 2015) | C++ with a Python wrapper |
| Mutalyzer 3 `mutalyzer.description_extractor` | mutalyzer 3.1.1 | description-extractor 3.0.0 and Mutalyzer 3's description model | C++ and Python |
| [mutalyzer-algebra](https://github.com/mutalyzer/algebra) `extract_sequence` | 1.5.2 | graph of all minimal alignments (Santcroos et al., 2025) | Python |

## Method

Datasets (`tools/benchmark.py`, all deterministic):

* **HGVS specification examples**: the 25 DNA and RNA examples of HGVS Nomenclature 21.1.5 that
  show both sequences (from `tests/test_hgvs_spec_examples.py`), as DNA and numbered from 1.
* **description-extractor test cases**: its 18 test cases; "as expected" is its own results.
* **random short pairs**: 5000 references of 5 to 60 bases (every third AT-rich) with 1 to 4
  random substitutions, deletions, insertions, duplications or inversions.
* **300 bp, 1 change** and **300 bp, 3 changes**: 1000 reads each of one random 300 bp
  reference with substitutions or 1-3 base deletions or insertions, far apart.
* **10 kb, 5 changes**: 20 sequences of a random 10 kb reference with 5 such changes.

Every description is applied to the reference by one checker (`apply_description`, which
understands all the tools' notations: repeats such as `insATCG[3]`, inserted copies such as
`ins5_16inv`, composite insertions such as `ins[2_13;GGTC]`):

* **correct**: it gives the observed sequence;
* **valid**: its positions are inside the reference (no position 0, none after the last), and
  its changes are in order and do not overlap;
* **as expected**: it is the description given by the specification (or by
  description-extractor, for its test cases).

All tools are called with `g.` coordinates. Time is the mean per call over one run, in one
Linux x86_64 environment with CPython 3.12 (the Mutalyzer 2 code needs Python 3.10).

To reproduce, install each tool (see the docstring of `tools/benchmark.py`), then:

```sh
PYTHONPATH=src:tests:tools python tools/benchmark.py run "hgvs_extractor" --out results/ve.json
PYTHONPATH=src:tests:tools python tools/benchmark.py report results/*.json
```

## Summary

* **Correctness.** hgvs_extractor gave a correct, valid description in every one of the 7063
  cases. Mutalyzer 3 and mutalyzer-algebra were always correct, but 3.3% and 4.3% of their
  descriptions of the random pairs use positions outside the reference: an insertion before the first base
  (`0_1insC`) or after the last (`7_8insG` for a 7-base reference), where HGVS has no flanking
  position (hgvs_extractor writes a delins of the first or last base). Mutalyzer 2 also gave
  41 wrong descriptions (0.8%; fixed in the description-extractor 3 used by Mutalyzer 3).
  countess-variant-caller gave 147 wrong descriptions of the random pairs (2.9%), and 11 of the
  1000 reads with three changes, e.g. `[128_129insT;128_129dup;...]` for `128_129insTGA`.
* **The HGVS specification.** hgvs_extractor gave the specification's description for 24 of
  the 25 examples; the other, `r.775_777delinsc`, needs `r.` coordinates (here all tools get `g.`;
  with `r.` it gives the specification's description). Mutalyzer 2 and 3 gave 19,
  mutalyzer-algebra 14 (it prefers repeat notation, e.g. `2_7AT[4]` for `6_7dup`) and
  countess-variant-caller 13.
* **Description size.** Mutalyzer 2 and 3 use the fewest changes (1.70 per random pair, against
  1.99): they minimize the length of the description, joining changes separated by unchanged
  bases into one delins, which HGVS describes separately (see "Comparison with Mutalyzer's
  description-extractor" in the README).
* **Speed.** countess-variant-caller is the fastest (5-12 µs per read), with its alignment in
  compiled code. hgvs_extractor, in Python, takes 14 µs for a 300 bp read with one change
  (Mutalyzer 2: 22 µs), 131 µs with three (Mutalyzer 2: 89 µs, Mutalyzer 3: 267 µs) and 0.7 ms
  for 10 kb with five changes (Mutalyzer 2: 0.75 ms, Mutalyzer 3: 0.93 ms, mutalyzer-algebra:
  19 ms), as it skips the unchanged stretches between changes.

## Results

| tool | version | Python |
|---|---|---|
| hgvs_extractor | 0.1.0 | 3.12.3 |
| hgvs_extractor (pure Python) | 0.1.0 | 3.12.3 |
| countess-variant-caller | 0.1.1 | 3.12.3 |
| Mutalyzer 2 description-extractor | git 155b79d | 3.10.20 |
| Mutalyzer 3 description_extractor | mutalyzer 3.1.1, description-extractor 3.0.0 | 3.12.3 |
| mutalyzer-algebra | 1.5.2 | 3.12.3 |

### HGVS specification examples (25 cases)

| tool | correct and valid HGVS | correct, invalid positions | wrong or error | as expected | changes per case | characters per case | time per case |
|---|---|---|---|---|---|---|---|
| hgvs_extractor | 25 | 0 | 0 | 24 | 1.04 | 7.4 | 30 µs |
| hgvs_extractor (pure Python) | 25 | 0 | 0 | 24 | 1.04 | 7.4 | 29 µs |
| countess-variant-caller | 25 | 0 | 0 | 13 | 1.60 | 9.9 | 8 µs |
| Mutalyzer 2 description-extractor | 25 | 0 | 0 | 19 | 1.16 | 7.6 | 88 µs |
| Mutalyzer 3 description_extractor | 25 | 0 | 0 | 19 | 1.16 | 7.6 | 77 µs |
| mutalyzer-algebra | 25 | 0 | 0 | 14 | 1.16 | 8.8 | 77 µs |

### description-extractor test cases (18 cases)

| tool | correct and valid HGVS | correct, invalid positions | wrong or error | as expected | changes per case | characters per case | time per case |
|---|---|---|---|---|---|---|---|
| hgvs_extractor | 18 | 0 | 0 | 18 | 1.39 | 9.9 | 60 µs |
| hgvs_extractor (pure Python) | 18 | 0 | 0 | 18 | 1.39 | 9.9 | 71 µs |
| countess-variant-caller | 18 | 0 | 0 | 12 | 1.61 | 12.4 | 6 µs |
| Mutalyzer 2 description-extractor | 18 | 0 | 0 | 18 | 1.39 | 9.9 | 56 µs |
| Mutalyzer 3 description_extractor | 18 | 0 | 0 | 17 | 1.39 | 9.8 | 90 µs |
| mutalyzer-algebra | 18 | 0 | 0 | 11 | 1.44 | 11.2 | 60 µs |

### random short pairs (5000 cases)

| tool | correct and valid HGVS | correct, invalid positions | wrong or error | changes per case | characters per case | time per case |
|---|---|---|---|---|---|---|
| hgvs_extractor | 5000 | 0 | 0 | 1.99 | 18.0 | 164 µs |
| hgvs_extractor (pure Python) | 5000 | 0 | 0 | 1.99 | 18.0 | 206 µs |
| countess-variant-caller | 4636 | 217 | 147 | 2.58 | 21.5 | 7 µs |
| Mutalyzer 2 description-extractor | 4796 | 163 | 41 | 1.70 | 16.5 | 38 µs |
| Mutalyzer 3 description_extractor | 4836 | 164 | 0 | 1.70 | 16.5 | 125 µs |
| mutalyzer-algebra | 4787 | 213 | 0 | 1.88 | 18.4 | 112 µs |

Invalid positions and exceptions: countess-variant-caller: 15 overlapping changes, 105 position 0, 97 position after the end; Mutalyzer 2 description-extractor: 64 position 0, 99 position after the end; Mutalyzer 3 description_extractor: 66 position 0, 98 position after the end; mutalyzer-algebra: 85 position 0, 128 position after the end.

### 300 bp, 1 change (1000 cases)

| tool | correct and valid HGVS | correct, invalid positions | wrong or error | changes per case | characters per case | time per case |
|---|---|---|---|---|---|---|
| hgvs_extractor | 1000 | 0 | 0 | 1.00 | 7.2 | 14 µs |
| hgvs_extractor (pure Python) | 1000 | 0 | 0 | 1.00 | 7.2 | 19 µs |
| countess-variant-caller | 1000 | 0 | 0 | 1.00 | 7.2 | 5 µs |
| Mutalyzer 2 description-extractor | 1000 | 0 | 0 | 1.00 | 7.2 | 22 µs |
| Mutalyzer 3 description_extractor | 1000 | 0 | 0 | 1.00 | 7.2 | 70 µs |
| mutalyzer-algebra | 1000 | 0 | 0 | 1.00 | 7.3 | 56 µs |

### 300 bp, 3 changes (1000 cases)

| tool | correct and valid HGVS | correct, invalid positions | wrong or error | changes per case | characters per case | time per case |
|---|---|---|---|---|---|---|
| hgvs_extractor | 1000 | 0 | 0 | 2.98 | 25.0 | 131 µs |
| hgvs_extractor (pure Python) | 1000 | 0 | 0 | 2.98 | 25.0 | 201 µs |
| countess-variant-caller | 989 | 0 | 11 | 3.12 | 25.9 | 12 µs |
| Mutalyzer 2 description-extractor | 1000 | 0 | 0 | 2.96 | 24.9 | 89 µs |
| Mutalyzer 3 description_extractor | 1000 | 0 | 0 | 2.96 | 24.9 | 267 µs |
| mutalyzer-algebra | 1000 | 0 | 0 | 2.96 | 25.7 | 324 µs |

### 10 kb, 5 changes (20 cases)

| tool | correct and valid HGVS | correct, invalid positions | wrong or error | changes per case | characters per case | time per case |
|---|---|---|---|---|---|---|
| hgvs_extractor | 20 | 0 | 0 | 5.00 | 48.6 | 696 µs |
| hgvs_extractor (pure Python) | 20 | 0 | 0 | 5.00 | 48.6 | 953 µs |
| countess-variant-caller | 20 | 0 | 0 | 5.20 | 50.9 | 1,960 µs |
| Mutalyzer 2 description-extractor | 20 | 0 | 0 | 5.00 | 48.6 | 752 µs |
| Mutalyzer 3 description_extractor | 20 | 0 | 0 | 5.00 | 48.4 | 929 µs |
| mutalyzer-algebra | 20 | 0 | 0 | 5.00 | 50.0 | 19 ms |
