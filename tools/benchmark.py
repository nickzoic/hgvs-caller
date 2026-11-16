#!/usr/bin/env python3
"""Benchmark variant extraction tools: HGVS descriptions from a reference and an observed sequence.

Each tool runs in an environment where it is installed; results are saved as JSON and then
scored together, so that every tool's descriptions are checked the same way:

    python tools/benchmark.py run TOOL --out results/TOOL.json     (in the tool's environment)
    python tools/benchmark.py report results/*.json                (anywhere: a Markdown report)

Run from the repository root with PYTHONPATH=src:tests:tools.  Tools (see TOOLS):

* variant_extractor, with and without rapidfuzz (this package)
* countess-variant-caller (pip install countess-variant-caller)
* Mutalyzer 2: https://github.com/mutalyzer/description-extractor, built as described in
  tools/compare_description_extractor.py
* Mutalyzer 3: mutalyzer.description_extractor (pip install mutalyzer "setuptools<80"),
  which uses description-extractor 3
* mutalyzer-algebra: algebra.extractor.extract_sequence (pip install mutalyzer-algebra)
"""

import argparse
import collections
import importlib.metadata
import json
import platform
import random
import re
import statistics
import sys
import time
from typing import Callable, Optional

from compare_description_extractor import mutate, split_allele
from hgvs_apply import revcomp

Extractor = Callable[[str, str], str]  # (reference, observed) -> description without "g."

# --------------------------------------------------------------------------- tools


def _variant_extractor(rapidfuzz: bool) -> tuple[Extractor, str]:
    # pylint: disable=import-outside-toplevel,protected-access
    import variant_extractor.extractor as module
    from variant_extractor import extract_variant_string

    if not rapidfuzz:
        module._edit_distance = None
    elif module._edit_distance is None:
        raise ImportError("rapidfuzz is not installed")
    return (lambda ref, alt: extract_variant_string("g.", ref, alt)[2:]), importlib.metadata.version(
        "variant_extractor"
    )


def _countess() -> tuple[Extractor, str]:
    from countess_variant_caller import find_variant_string  # pylint: disable=import-outside-toplevel,import-error

    return (lambda ref, alt: find_variant_string("g.", ref, alt)[2:]), importlib.metadata.version(
        "countess-variant-caller"
    )


def _mutalyzer2() -> tuple[Extractor, str]:
    import extractor  # pylint: disable=import-outside-toplevel,import-error

    if not hasattr(extractor, "describe"):  # description-extractor 3 has a different interface
        raise ImportError("not the Mutalyzer 2 description-extractor")
    return (lambda ref, alt: str(extractor.describe_dna(ref, alt))), "git 155b79d"


def _mutalyzer3() -> tuple[Extractor, str]:
    # pylint: disable=import-outside-toplevel,import-error
    from mutalyzer.description_extractor import description_extractor

    versions = [importlib.metadata.version(p) for p in ("mutalyzer", "description-extractor")]
    return description_extractor, f"mutalyzer {versions[0]}, description-extractor {versions[1]}"


def _algebra() -> tuple[Extractor, str]:
    from algebra.extractor import extract_sequence, to_hgvs  # pylint: disable=import-outside-toplevel,import-error

    return (lambda ref, alt: to_hgvs(extract_sequence(ref, alt)[0], ref)), importlib.metadata.version(
        "mutalyzer-algebra"
    )


TOOLS: dict[str, Callable[[], tuple[Extractor, str]]] = {
    "variant_extractor": lambda: _variant_extractor(True),
    "variant_extractor (pure Python)": lambda: _variant_extractor(False),
    "countess-variant-caller": _countess,
    "Mutalyzer 2 description-extractor": _mutalyzer2,
    "Mutalyzer 3 description_extractor": _mutalyzer3,
    "mutalyzer-algebra": _algebra,
}

# --------------------------------------------------------------------------- datasets


def _renumber(description: str, by: int) -> str:
    """add `by` to every position in a description"""
    return re.sub(r"\d+", lambda m: str(int(m.group()) + by), description)


def _spec_examples() -> list[dict]:
    """nucleotide examples of the HGVS specification, as DNA, numbered from 1"""
    from test_hgvs_spec_examples import EXAMPLES  # pylint: disable=import-outside-toplevel

    cases = []
    for expected, kind, ref, alt, offset, _, _ in EXAMPLES:
        if kind in ("dna", "rna"):
            body = re.sub(  # RNA bases to DNA, leaving del, dup, inv and ins alone
                r"(delins|del|dup|inv|ins)|([acgun]+)",
                lambda m: m.group(1) or m.group(2).upper().replace("U", "T"),
                expected[2:],
            )
            body = _renumber(body, -offset)
            cases.append({"ref": ref.upper().replace("U", "T"), "alt": alt.upper().replace("U", "T"), "expected": body})
    return cases


def _description_extractor_examples() -> list[dict]:
    from test_description_extractor import EXAMPLES  # pylint: disable=import-outside-toplevel

    return [{"ref": ref, "alt": alt, "expected": expected} for _, ref, alt, expected in EXAMPLES]


def _random_pairs(cases: int, seed: int) -> list[dict]:
    """short random pairs, every third AT-rich, with 1-4 random changes"""
    rng = random.Random(seed)
    result = []
    for n in range(cases):
        ref = "".join(rng.choices("ACGT" if n % 3 else "AT", k=rng.randint(5, 60)))
        result.append({"ref": ref, "alt": mutate(rng, ref)})
    return result


def _reads(length: int, changes: int, cases: int, seed: int) -> list[dict]:
    """reads of one random reference with a few changes far apart"""
    rng = random.Random(seed)
    ref = "".join(rng.choices("ACGT", k=length))
    result = []
    for _ in range(cases):
        seq = list(ref)
        for i in sorted(rng.sample(range(5, length - 5), changes), reverse=True):
            op = rng.random()
            if op < 0.6:
                seq[i] = rng.choice([b for b in "ACGT" if b != seq[i]])
            elif op < 0.8:
                del seq[i : i + rng.randint(1, 3)]
            else:
                seq[i:i] = rng.choices("ACGT", k=rng.randint(1, 3))
        result.append({"ref": ref, "alt": "".join(seq)})
    return result


DATASETS: dict[str, Callable[[], list[dict]]] = {
    "HGVS specification examples": _spec_examples,
    "description-extractor test cases": _description_extractor_examples,
    "random short pairs": lambda: _random_pairs(5000, 1),
    "300 bp, 1 change": lambda: _reads(300, 1, 1000, 2),
    "300 bp, 3 changes": lambda: _reads(300, 3, 1000, 3),
    "10 kb, 5 changes": lambda: _reads(10_000, 5, 20, 4),
}

# --------------------------------------------------------------------------- run


def run(tool: str, out: str, time_limit: float) -> None:
    """run one tool on all datasets and save its descriptions and timings"""
    extract, version = TOOLS[tool]()
    results = {}
    for name, make in DATASETS.items():
        cases = make()
        extract(cases[0]["ref"], cases[0]["alt"])  # warm up
        outputs: list[Optional[str]] = []
        times = []
        started = time.perf_counter()
        for case in cases:
            t0 = time.perf_counter()
            try:
                outputs.append(extract(case["ref"], case["alt"]))
            except Exception as error:  # pylint: disable=broad-exception-caught
                outputs.append(f"<{type(error).__name__}>")
            times.append(time.perf_counter() - t0)
            if time.perf_counter() - started > time_limit:
                break
        results[name] = {"outputs": outputs, "times": times, "complete": len(outputs) == len(cases)}
        print(f"{tool}: {name}: {len(outputs)} cases, {statistics.mean(times) * 1e6:.0f} µs each", file=sys.stderr)
    with open(out, "w", encoding="utf-8") as file:
        json.dump({"tool": tool, "version": version, "python": platform.python_version(), "results": results}, file)


# --------------------------------------------------------------------------- scoring


def _inserted(arg: str, ref: str) -> str:
    """an inserted sequence: literal, a copy of ref (1_5, 1_5inv), a repeat (AT[3]), or [a;b]"""
    if arg.startswith("[") and arg.endswith("]") and not re.fullmatch(r"[ACGTN]+\[\d+\]", arg):
        return "".join(_inserted(part, ref) for part in arg[1:-1].split(";"))
    if m := re.fullmatch(r"(\d+)_(\d+)(inv)?", arg):
        copy = ref[int(m.group(1)) - 1 : int(m.group(2))]
        return revcomp(copy) if m.group(3) else copy
    if m := re.fullmatch(r"([ACGTN]+)\[(\d+)\]", arg):
        return m.group(1) * int(m.group(2))
    if not re.fullmatch(r"[ACGTN]+", arg):
        raise ValueError(f"inserted sequence {arg}")
    return arg


def apply_description(description: str, ref: str) -> tuple[str, str]:
    """(the observed sequence a description gives, why its positions are not valid HGVS or
    ""): positions within the reference, changes in order and not overlapping"""
    if description == "=":
        return ref, ""
    changes = []
    positions = []
    for edit in split_allele(description):
        if m := re.fullmatch(r"(\d+)([ACGTN])>([ACGTN])", edit):
            start = end = int(m.group(1))
            if ref[start - 1] != m.group(2):
                raise ValueError(f"{edit}: reference base is {ref[start - 1]}")
            new = m.group(3)
        elif m := re.fullmatch(r"(\d+)(?:_(\d+))?([ACGTN]+)\[(\d+)\]", edit):  # repeat
            start, end = int(m.group(1)), int(m.group(2) or m.group(1))
            new = m.group(3) * int(m.group(4))
        elif m := re.fullmatch(r"(\d+)(?:_(\d+))?(delins|del|dup|inv|ins)(.*)", edit):
            start, end, op, arg = int(m.group(1)), int(m.group(2) or m.group(1)), m.group(3), m.group(4)
            if op == "ins":
                if end != start + 1:
                    raise ValueError(f"{edit}: not between adjacent positions")
                positions += [start, end]  # both flanks must exist
                changes.append((start + 1, start, _inserted(arg, ref)))  # (insert before start + 1)
                continue
            old = ref[start - 1 : end] if start >= 1 else ""
            new = {"del": "", "dup": old + old, "inv": revcomp(old)}.get(op, "")
            if op == "delins":
                new = _inserted(arg, ref)
        else:
            raise ValueError(f"cannot apply {edit}")
        positions += [start, end]
        changes.append((start, end, new))
    problem = ""
    if min(positions) < 1:
        problem = "position 0"
    elif max(positions) > len(ref):
        problem = "position after the end"
    elif any(a[1] >= b[0] for a, b in zip(changes, changes[1:])):
        problem = "overlapping changes"
    seq = ref
    for start, end, new in sorted(changes, reverse=True):
        seq = seq[: start - 1] + new + seq[end:]
    return seq, problem


def score(cases: list[dict], outputs: list[str]) -> collections.Counter:
    """counts for one tool on one dataset"""
    counts: collections.Counter = collections.Counter()
    for case, output in zip(cases, outputs):
        counts["cases"] += 1
        if output.startswith("<"):
            counts["errors"] += 1
            continue
        try:
            seq, problem = apply_description(output, case["ref"])
        except (ValueError, IndexError):
            seq, problem = "", ""
        counts["correct"] += seq == case["alt"]
        if seq == case["alt"]:
            counts[f"invalid: {problem}" if problem else "correct and valid"] += 1
        counts["changes"] += len(split_allele(output)) if output != "=" else 0
        counts["characters"] += len(output)
        if "expected" in case:
            counts["as expected"] += output == case["expected"]
    return counts


def report(files: list[str]) -> None:
    """a Markdown report of the saved results"""
    runs = [json.load(open(f, encoding="utf-8")) for f in files]  # pylint: disable=consider-using-with
    datasets = {name: make() for name, make in DATASETS.items()}
    print("| tool | version | Python |\n|---|---|---|")
    for r in runs:
        print(f"| {r['tool']} | {r['version']} | {r['python']} |")
    for name, cases in datasets.items():
        print(f"\n### {name} ({len(cases)} cases)\n")
        expected = "expected" in cases[0]
        head = "| tool | correct and valid HGVS | correct, invalid positions | wrong or error |"
        head += " as expected |" if expected else ""
        head += " changes per case | characters per case | time per case |"
        print(head + "\n|" + "---|" * (head.count("|") - 1))
        notes = []
        for r in runs:
            result = r["results"][name]
            c = score(cases, result["outputs"])
            n = c["cases"]
            row = f"| {r['tool']} | {c['correct and valid']} | {c['correct'] - c['correct and valid']} |"
            row += f" {n - c['correct']} |"
            if expected:
                row += f" {c['as expected']} |"
            t = statistics.mean(result["times"]) * 1e6
            timing = f"{t:,.0f} µs" if t < 10_000 else f"{t / 1000:,.0f} ms"
            if not result["complete"]:
                timing += f" (stopped after {n})"
            row += f" {c['changes'] / n:.2f} | {c['characters'] / n:.1f} | {timing} |"
            print(row)
            problems = [f"{v} {k[9:]}" for k, v in sorted(c.items()) if k.startswith("invalid: ")]
            if c["errors"]:
                problems.append(f"{c['errors']} exceptions")
            if problems:
                notes.append(f"{r['tool']}: " + ", ".join(problems))
        if notes:
            print("\nInvalid positions and exceptions: " + "; ".join(notes) + ".")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    run_parser = sub.add_parser("run")
    run_parser.add_argument("tool", choices=TOOLS)
    run_parser.add_argument("--out", required=True)
    run_parser.add_argument("--time-limit", type=float, default=120, help="seconds per dataset")
    report_parser = sub.add_parser("report")
    report_parser.add_argument("files", nargs="+")
    args = parser.parse_args()
    if args.command == "run":
        run(args.tool, args.out, args.time_limit)
    else:
        report(args.files)


if __name__ == "__main__":
    main()
