# eukcode-scan: readers for Codetta 2.0 output files.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Parse Codetta's one-line results summary (-s) and its .genetic_code.out file.

Formats follow codetta.py (github.com/kshulgina/codetta, GPL-3.0; not bundled here):
summary header  prefix,profile_db,evalue_threshold,prob_threshold,max_fraction,
                excluded_domains,inferred_gencode,inferred_gencode_best_models
.out table rows '%-10s%-11s%-10s%-9s%-11i%-10i' = codon, inference, std code, diff?,
                N aligned, N used; then log decoding probabilities; then the final code.
"""

import csv
import io

from .codes import normalize_codon


def parse_summary(text):
    """Return a list of dicts, one per summary line."""
    rows = []
    reader = csv.reader(io.StringIO(text))
    header = None
    for row in reader:
        if not row:
            continue
        if row[0] == "prefix":
            header = row
            continue
        if header is None:
            header = ["prefix", "profile_db", "evalue_threshold", "prob_threshold", "max_fraction",
                      "excluded_domains", "inferred_gencode", "inferred_gencode_best_models"]
        rows.append(dict(zip(header, row)))
    for r in rows:
        for key in ("inferred_gencode", "inferred_gencode_best_models"):
            if key in r and len(r[key]) != 64:
                raise ValueError("%s has %d characters, expected 64" % (key, len(r[key])))
    return rows


def parse_inference_out(text):
    """Parse a .genetic_code.out file.

    Returns {'args': {...}, 'codons': {codon: {...}}, 'logp_header': [...],
             'logp': {codon: [21 floats]}, 'final': str or None}.
    """
    args, codons, logp = {}, {}, {}
    logp_header = []
    final = None
    section = "args"
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        s = line.strip()
        if s.startswith("# Analysis arguments"):
            section = "args"
        elif s.startswith("# codon") and "inference" in s:
            section = "table"
        elif s.startswith("# Log decoding probabilities"):
            section = "logp"
            if i + 1 < len(lines) and lines[i + 1].strip().startswith("# codon"):
                logp_header = lines[i + 1].replace("# codon", "").split()
                i += 1
        elif s.startswith("# Final genetic code inference"):
            if i + 1 < len(lines):
                final = lines[i + 1].strip()
            break
        elif s.startswith("#") or not s:
            pass
        else:
            parts = s.split()
            if section == "args" and len(parts) >= 1:
                args[parts[0]] = " ".join(parts[1:])
            elif section == "table" and len(parts) >= 6:
                codon = normalize_codon(parts[0])
                codons[codon] = {
                    "inference": parts[1],
                    "std": parts[2],
                    "diff": parts[3],
                    "n_aligned": int(parts[4]),
                    "n_used": int(parts[5]),
                }
            elif section == "logp" and len(parts) >= 2:
                codon = normalize_codon(parts[0])
                logp[codon] = [float(x) for x in parts[1:]]
        i += 1
    if final is not None and len(final) != 64:
        raise ValueError("final genetic code has %d characters" % len(final))
    return {"args": args, "codons": codons, "logp_header": logp_header, "logp": logp, "final": final}
