#!/usr/bin/env python3
"""Compare variant_extractor with Mutalyzer's description-extractor on random sequence pairs.

https://github.com/mutalyzer/description-extractor (C++ with a SWIG wrapper; it builds with
Python 3.10, SWIG from PyPI, biopython and crossmapper-legacy, see the README of this
repository's tests).  Run with description-extractor importable:

    PYTHONPATH=src:tests:path/to/description-extractor python tools/compare_description_extractor.py

Prints a summary of how the descriptions compare.  With --save FILE it writes the pairs on
which both give the same description to FILE as JSON, for tests/test_description_extractor.py.
"""

import argparse
import collections
import json
import random
import re

from extractor import describe_dna  # pylint: disable=import-error
from hgvs_apply import apply_nucleotide, revcomp

from variant_extractor import extract_variant_string


def mutate(rng: random.Random, seq: str) -> str:
    """1-4 random substitutions, deletions, insertions, duplications or inversions"""
    s = list(seq)
    for _ in range(rng.randint(1, 4)):
        i, n, op = rng.randrange(len(s) + 1), rng.randint(1, 6), rng.choice("sdiDv")
        if op == "s" and i < len(s):
            s[i] = rng.choice("ACGT")
        elif op == "d":
            del s[i : i + n]
        elif op == "i":
            s[i:i] = rng.choices("ACGT", k=n)
        elif op == "D":
            s[i:i] = s[max(0, i - n) : i]
        elif op == "v":
            s[i : i + n] = list(revcomp("".join(s[i : i + n])))
    return "".join(s) or "A"


def split_allele(description: str) -> list[str]:
    """top-level changes of a description-extractor description (composite insertions contain ';')"""
    if description == "=":
        return []
    if not (description.startswith("[") and description.endswith("]")):
        return [description]
    parts, depth, current = [], 0, ""
    for char in description[1:-1]:
        depth += (char == "[") - (char == "]")
        if char == ";" and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += char
    return parts + [current]


def apply_extractor(description: str, ref: str) -> str:
    """apply a description-extractor description, including composite insertions (ins[2_13;GGTC]) and
    insertions at position 0 or after the end, which hgvs_apply does not accept"""

    def inserted(arg: str) -> str:
        if arg.startswith("["):
            return "".join(inserted(x) for x in arg[1:-1].split(";"))
        if m := re.fullmatch(r"(\d+)_(\d+)(inv)?", arg):
            s = ref[int(m.group(1)) - 1 : int(m.group(2))]
            return revcomp(s) if m.group(3) else s
        return arg

    seq = ref
    edits = split_allele(description)
    for edit in sorted(edits, key=lambda e: int(re.match(r"\d+", e).group()), reverse=True):  # type: ignore
        if m := re.fullmatch(r"(\d+)([ACGT])>([ACGT])", edit):
            i = int(m.group(1)) - 1
            seq = seq[:i] + m.group(3) + seq[i + 1 :]
            continue
        m = re.fullmatch(r"(\d+)(?:_(\d+))?(delins|del|dup|inv|ins)(.*)", edit)
        if not m:
            raise ValueError(edit)
        a = int(m.group(1)) - 1
        b = int(m.group(2)) - 1 if m.group(2) else a
        op, arg = m.group(3), m.group(4)
        if op == "del":
            seq = seq[:a] + seq[b + 1 :]
        elif op == "dup":
            seq = seq[: b + 1] + seq[a : b + 1] + seq[b + 1 :]
        elif op == "inv":
            seq = seq[:a] + revcomp(seq[a : b + 1]) + seq[b + 1 :]
        elif op == "delins":
            seq = seq[:a] + inserted(arg) + seq[b + 1 :]
        else:
            seq = seq[: a + 1] + inserted(arg) + seq[a + 1 :]
    return seq


def extractor_problem(description: str, ref: str, alt: str) -> str:
    """why a description-extractor description is not a correct HGVS description, or ''"""
    try:
        if apply_extractor(description, ref) != alt:
            return "does not give the variant sequence"
    except (ValueError, IndexError):
        return "does not give the variant sequence"
    positions = [int(p) for p in re.findall(r"(?<![\[;])\b\d+", re.sub(r"\[.*?\]", "", description))]
    if any(p < 1 or p > len(ref) for p in positions):
        return "position outside the reference (0, or after the end)"
    spans = []
    for edit in split_allele(description):
        nums = [int(x) for x in re.findall(r"\d+", edit.split("ins")[0])]
        spans.append((nums[0], nums[-1], "ins" in edit and "delins" not in edit))
    for (_, e1, ins1), (s2, _, ins2) in zip(spans, spans[1:]):
        if s2 < e1 or (s2 == e1 and not (ins1 or ins2)):
            return "overlapping changes"
    return ""


def compare(cases: int, seed: int):
    rng = random.Random(seed)
    summary: collections.Counter = collections.Counter()
    examples: dict = collections.defaultdict(list)
    agreed = []
    for n in range(cases):
        ref = "".join(rng.choices("ACGT" if n % 3 else "AT", k=rng.randint(5, 60)))
        alt = mutate(rng, ref)
        theirs = str(describe_dna(ref, alt))
        ours = extract_variant_string("g.", ref, alt)[2:]
        assert apply_nucleotide("g." + ours, ref) == alt, (ref, alt, ours)
        if theirs == ours:
            summary["same description"] += 1
            agreed.append((ref, alt, ours))
            continue
        problem = extractor_problem(theirs, ref, alt)
        if problem:
            key = f"different: description-extractor's description is wrong: {problem}"
        elif len(split_allele(theirs)) < len(split_allele(ours)) and "delins" in theirs:
            key = "different: description-extractor joins separate changes in a delins (HGVS: describe them separately)"
        else:
            key = "different: both correct, different alignment"
        summary[key] += 1
        if len(examples[key]) < 3:
            examples[key].append((ref, alt, theirs, ours))
    return summary, examples, agreed


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cases", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--save", help="write the agreed pairs to this JSON file")
    parser.add_argument("--save-limit", type=int, default=500)
    args = parser.parse_args()
    summary, examples, agreed = compare(args.cases, args.seed)
    for key, count in sorted(summary.items(), key=lambda item: -item[1]):
        print(f"{count:7d} {100 * count / args.cases:5.1f}%  {key}")
        for ref, alt, theirs, ours in examples.get(key, []):
            print(f"            {ref} > {alt}")
            print(f"              description-extractor {theirs}")
            print(f"              variant_extractor     {ours}")
    if args.save:
        with open(args.save, "w", encoding="utf-8") as out:  # one pair per line
            out.write("[\n" + ",\n".join(json.dumps(case) for case in agreed[: args.save_limit]) + "\n]\n")


if __name__ == "__main__":
    main()
