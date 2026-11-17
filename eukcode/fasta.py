# eukcode-scan: FASTA reading and deterministic subsampling.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Read FASTA and draw a reproducible sample of fixed-size chunks.

Large eukaryotic assemblies are cut into chunks (default 100 kb), chunks that are short or
mostly N are dropped, and a seeded shuffle picks chunks until a cap in base pairs is reached.
The same input, cap and seed always give the same sample.

`sample_file` gives exactly the same sample as `sample_chunks(chunks(read_fasta(path)))` but
reads the file twice instead of holding every chunk in memory, so multi-gigabase protist
assemblies fit in a small container. `disjoint_from` draws the pre-registered second sample
(seed 1) from the chunks that the first sample (seed 0) did not use.
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


def _keep(piece, min_chunk_bp, max_n_frac):
    if len(piece) < min_chunk_bp:
        return False
    return piece.count("N") / len(piece) <= max_n_frac


def chunks(records, chunk_bp=100_000, min_chunk_bp=5_000, max_n_frac=0.1):
    """Cut records into chunks. Returns a list of (chunk_name, seq) in input order."""
    out = []
    for idx, (_name, seq) in enumerate(records):
        for start in range(0, len(seq), chunk_bp):
            piece = seq[start:start + chunk_bp]
            if _keep(piece, min_chunk_bp, max_n_frac):
                out.append(("s%d_%d" % (idx, start), piece))
    return out


def chunk_index(records, chunk_bp=100_000, min_chunk_bp=5_000, max_n_frac=0.1):
    """Like chunks(), but returns (chunk_name, record_index, start, length) without sequence."""
    out = []
    for idx, (_name, seq) in enumerate(records):
        for start in range(0, len(seq), chunk_bp):
            piece = seq[start:start + chunk_bp]
            if _keep(piece, min_chunk_bp, max_n_frac):
                out.append(("s%d_%d" % (idx, start), idx, start, len(piece)))
    return out


def choose(lengths, cap_bp, seed=0):
    """Indices picked by the seeded shuffle until cap_bp is reached, in increasing order.

    If everything fits under the cap, every index is returned.
    """
    if sum(lengths) <= cap_bp:
        return list(range(len(lengths)))
    order = list(range(len(lengths)))
    random.Random(seed).shuffle(order)
    picked, size = [], 0
    for i in order:
        if size >= cap_bp:
            break
        picked.append(i)
        size += lengths[i]
    picked.sort()
    return picked


def sample_chunks(chunk_list, cap_bp, seed=0):
    """Pick chunks in a seeded random order until cap_bp is reached.

    If everything fits under the cap, all chunks are returned in input order.
    """
    picked = choose([len(s) for _, s in chunk_list], cap_bp, seed)
    return [chunk_list[i] for i in picked]


def _extract(path, wanted):
    """wanted: list of (name, record_index, start, length). Returns [(name, seq)] in that order."""
    by_rec = {}
    for name, rec, start, length in wanted:
        by_rec.setdefault(rec, []).append((name, start, length))
    got = {}
    for idx, (_n, seq) in enumerate(read_fasta(path)):
        for name, start, length in by_rec.get(idx, ()):
            got[name] = seq[start:start + length]
    return [(name, got[name]) for name, _r, _s, _l in wanted]


def sample_file(path, cap_bp, seed=0, disjoint_from=None, chunk_bp=100_000, min_chunk_bp=5_000,
                max_n_frac=0.1):
    """Sample a FASTA file with low memory. Returns (sample_records, n_chunks_total).

    disjoint_from: a seed whose sample (at the same cap) is excluded first, so the result
    shares no chunk with it.
    """
    idx = chunk_index(read_fasta(path), chunk_bp, min_chunk_bp, max_n_frac)
    pool = list(range(len(idx)))
    if disjoint_from is not None:
        first = set(choose([x[3] for x in idx], cap_bp, disjoint_from))
        pool = [i for i in pool if i not in first]
    sub = [idx[i] for i in pool]
    picked = choose([x[3] for x in sub], cap_bp, seed)
    return _extract(path, [sub[i] for i in picked]), len(idx)


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
