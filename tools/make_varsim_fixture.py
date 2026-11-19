#!/usr/bin/env python3
"""Generate tests/data/varsim_extract_results.json: our test cases run through VarSim's extraction.

VarSim (https://github.com/liu-sun/VarSim, MIT) is built from a local clone, then its
``extract_with_positions`` is run on each reference / observed pair, with ``g.`` coordinates and
no accession. The pairs come from the benchmark datasets (tools/benchmark.py): the HGVS
specification examples, description-extractor's test cases, and samples of the random short pairs
and the simulated reads. Each entry is {"ref", "alt", "varsim", "source"}.

    git clone --depth 1 https://github.com/liu-sun/VarSim.git /tmp/VarSim
    EMAIL=you@example.org pip install /tmp/VarSim
    VARSIM_SOURCE=/tmp/VarSim PYTHONPATH=src:tests:tools python tools/make_varsim_fixture.py \
        tests/data/varsim_extract_results.json

Run it in an environment with the benchmark's dependencies (see tools/benchmark.py).
"""

import json
import os
import subprocess
import sys

import benchmark

SAMPLES = {  # dataset name: number of cases to take (all of the small ones)
    "HGVS specification examples": None,
    "description-extractor test cases": None,
    "random short pairs": 500,
    "300 bp, 1 change": 300,
    "300 bp, 3 changes": 150,
}


def _varsim_commit(path: str) -> str:
    """the commit of the local VarSim clone, if the path is a git checkout"""
    try:
        return subprocess.run(
            ["git", "-C", path, "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main() -> None:
    output = sys.argv[1]
    extract, version = benchmark.TOOLS["VarSim extract_with_positions"]()
    entries = []
    for name, limit in SAMPLES.items():
        cases = benchmark.DATASETS[name]()
        for case in cases[:limit]:
            entries.append(
                {"ref": case["ref"], "alt": case["alt"], "varsim": extract(case["ref"], case["alt"]), "source": name}
            )
    document = {
        "tool": "VarSim extract_with_positions",
        "varsim_version": version,
        "varsim_commit": _varsim_commit(os.environ.get("VARSIM_SOURCE", ".")),
        "python": sys.version.split()[0],
        "entries": entries,
    }
    with open(output, "w", encoding="utf-8") as file:
        json.dump(document, file, indent=1)
    print(f"{len(entries)} entries written to {output}", file=sys.stderr)


if __name__ == "__main__":
    main()
