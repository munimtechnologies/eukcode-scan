# eukcode-scan: FASTA reading and deterministic subsampling.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Read FASTA and draw a reproducible sample of fixed-size chunks.

Large eukaryotic assemblies are cut into chunks (default 100 kb), chunks that are short or
mostly N are dropped, and a seeded shuffle picks chunks until a cap in base pairs is reached.
The same input, cap and seed always give the same sample.
"""

import gzip
import io
import random


def _open_text(path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", errors="replace")
    return open(path, "r", encoding="ascii", errors="replace")


def read_fasta(path_or_handle):
    """Yield (name, sequence) with the sequence upper-cased and whitespace removed."""
    handle = _open_text(path_or_handle) if isinstance(path_or_handle, (str, bytes)) or hasattr(path_or_handle, "__fspath__") else path_or_handle
    name, parts = None, []
    try:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper()
                name = line[1:].split()[0] if len(line) > 1 else "unnamed"
                parts = []
            else:
                parts.append(line)
        if name is not None:
            yield name, "".join(parts).upper()
    finally:
        if handle is not path_or_handle:
            handle.close()


def chunks(records, chunk_bp=100_000, min_chunk_bp=5_000, max_n_frac=0.1):
    """Cut records into chunks. Returns a list of (chunk_name, seq) in input order."""
    out = []
    for idx, (_name, seq) in enumerate(records):
        for start in range(0, len(seq), chunk_bp):
            piece = seq[start:start + chunk_bp]
            if len(piece) < min_chunk_bp:
                continue
            n_count = piece.count("N")
            if n_count / len(piece) > max_n_frac:
                continue
            out.append(("s%d_%d" % (idx, start), piece))
    return out


def sample_chunks(chunk_list, cap_bp, seed=0):
    """Pick chunks in a seeded random order until cap_bp is reached.

    If everything fits under the cap, all chunks are returned in input order.
    """
    total = sum(len(s) for _, s in chunk_list)
    if total <= cap_bp:
        return list(chunk_list)
    order = list(range(len(chunk_list)))
    random.Random(seed).shuffle(order)
    picked, size = [], 0
    for i in order:
        if size >= cap_bp:
            break
        picked.append(i)
        size += len(chunk_list[i][1])
    picked.sort()
    return [chunk_list[i] for i in picked]


def write_fasta(records, path, width=80):
    with open(path, "w", encoding="ascii") as fh:
        for name, seq in records:
            fh.write(">%s\n" % name)
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def gc_fraction(records):
    gc = at = 0
    for _, seq in records:
        gc += seq.count("G") + seq.count("C")
        at += seq.count("A") + seq.count("T")
    return gc / (gc + at) if gc + at else 0.0
