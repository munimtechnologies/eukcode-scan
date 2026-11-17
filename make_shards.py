#!/usr/bin/env python3
# eukcode-scan: write screen shards from the pilot census.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Usage: python make_shards.py census_species.csv shards/ --per-shard 30"""

import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from eukcode import shards  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("census")
    ap.add_argument("out_dir")
    ap.add_argument("--per-shard", type=int, default=30)
    ap.add_argument("--max-genome-gb", type=float, default=3.0)
    args = ap.parse_args(argv)
    rows = shards.read_species(args.census)
    plan, skipped = shards.plan(rows, args.per_shard, int(args.max_genome_gb * 1e9))
    paths = shards.write_shards(plan, args.out_dir)
    with open(os.path.join(args.out_dir, "skipped.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["accession", "organism", "group", "length", "skip_reason"],
                           extrasaction="ignore")
        w.writeheader()
        for r in skipped:
            w.writerow(r)
    summary = {"species_in": len(rows), "shards": len(paths), "species_planned": sum(len(p) for p in plan),
               "skipped": len(skipped)}
    with open(os.path.join(args.out_dir, "shards_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
