# eukcode-scan: split the census into screen shards.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Turn census_species.csv (one assembly per taxon, from run_pilot.py) into shard CSVs.

Deterministic: rows are filtered, sorted by (group, organism, accession) and cut into
consecutive shards, so the same census always gives the same shards.
"""

import csv
import os

SHARD_FIELDS = ["shard", "group", "accession", "organism", "tax_id", "level", "length",
                "contig_n50", "ncbi_gencode", "lineage"]


def _int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0


def read_species(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def plan(species, per_shard=30, max_genome_bp=3_000_000_000, min_genome_bp=1_000_000):
    """Return (shards, skipped). shards: list of lists of rows; skipped: rows with 'skip_reason'."""
    if per_shard < 1:
        raise ValueError("per_shard must be >= 1")
    keep, skipped = [], []
    seen = set()
    for row in species:
        acc = (row.get("accession") or "").strip()
        length = _int(row.get("length"))
        reason = ""
        if not acc:
            reason = "no accession"
        elif acc in seen:
            reason = "duplicate accession"
        elif length < min_genome_bp:
            reason = "assembly shorter than %d bp" % min_genome_bp
        elif length > max_genome_bp:
            reason = "assembly longer than %d bp" % max_genome_bp
        if reason:
            r = dict(row)
            r["skip_reason"] = reason
            skipped.append(r)
            continue
        seen.add(acc)
        keep.append(row)
    keep.sort(key=lambda r: (r.get("group", ""), r.get("organism", ""), r.get("accession", "")))
    shards = [keep[i:i + per_shard] for i in range(0, len(keep), per_shard)]
    return shards, skipped


def write_shards(shards, out_dir, prefix="shard"):
    """Write shard_00.csv, shard_01.csv, ... and return their paths."""
    os.makedirs(out_dir, exist_ok=True)
    paths = []
    for n, rows in enumerate(shards):
        name = "%s_%02d" % (prefix, n)
        p = os.path.join(out_dir, name + ".csv")
        with open(p, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=SHARD_FIELDS, extrasaction="ignore")
            w.writeheader()
            for r in rows:
                r2 = dict(r)
                r2["shard"] = name
                w.writerow(r2)
        paths.append(p)
    return paths


def gencode_of(row):
    """NCBI nuclear code id from a census/shard row, or None."""
    v = (row.get("ncbi_gencode") or "").strip()
    return int(v) if v.isdigit() else None
