# eukcode-scan: a systematic screen for alternative nuclear genetic codes in eukaryotes.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Read and grade KACI output (addendum A8, the fast triage path).

KACI (Melnykov 2026, PNAS 123:e2610659123; lib-kaci, MIT licence) writes one JSON
object per query: a few metadata keys plus one key per codon (TTT..GGG). Each codon
entry holds a probability for each of the 20 amino acids and "?", and "n_", the
number of reference k-mer hits that informed it.

The calling rule below is the one in lib-kaci's examples/2_evaluate_results.py
(infer_aa): a codon with n_ <= 3 is uncalled; otherwise the most probable symbol
(first in AAS order on a tie) is called if its probability is strictly above
0.9999, and a confident "?" stays uncalled. A8 pre-registers this rule unchanged.
"""
from __future__ import annotations

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AAS = "ACDEFGHIKLMNPQRSTVWY?"
MIN_N = 3
HIGH_CONF = 0.9999

# NCBI translation tables in TCAG order, "*" = stop.
TABLES = {
    1: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
}
STANDARD = TABLES[1]

# lib-kaci's evaluate script does not flag AGA/AGG read as lysine; A8 reports
# them separately and does not count them as unexpected.
KNOWN_ARTEFACTS = {("AGA", "K"), ("AGG", "K")}

ORGANELLE_WORDS = ("mitochondri", "apicoplast", "chloroplast", "plastid", "kinetoplast")


def codon_index(codon: str) -> int:
    codon = codon.upper()
    return 16 * BASES.index(codon[0]) + 4 * BASES.index(codon[1]) + BASES.index(codon[2])


def infer_codon(entry: dict, min_n: int = MIN_N, high_conf: float = HIGH_CONF) -> str:
    """One codon's call: an amino acid letter, or "?" when uncalled."""
    if int(entry["n_"]) <= min_n:
        return "?"
    probs = [float(entry[a]) for a in AAS]
    best = probs.index(max(probs))
    if probs[best] > high_conf:
        return AAS[best]
    return "?"


def infer_code(result: dict, min_n: int = MIN_N, high_conf: float = HIGH_CONF) -> str:
    """64-letter inferred code in TCAG order ("?" = uncalled)."""
    out = []
    for codon in CODONS:
        if codon not in result:
            raise ValueError(f"KACI result has no entry for codon {codon}")
        out.append(infer_codon(result[codon], min_n, high_conf))
    return "".join(out)


def compare(inferred: str, expected: str) -> dict:
    """Compare an inferred code with an expected NCBI table."""
    if len(inferred) != 64 or len(expected) != 64:
        raise ValueError("codes must have 64 letters")
    res = {"sense_expected": [], "sense_unexpected": [], "sense_uncalled": [],
           "stop_called": [], "stop_uncalled": [], "artefacts": []}
    for i, codon in enumerate(CODONS):
        exp, got = expected[i], inferred[i]
        if exp == "*":
            if got == "?":
                res["stop_uncalled"].append(codon)
            else:
                res["stop_called"].append([codon, got])
        elif got == "?":
            res["sense_uncalled"].append(codon)
        elif got == exp:
            res["sense_expected"].append(codon)
        elif (codon, got) in KNOWN_ARTEFACTS:
            res["artefacts"].append([codon, exp, got])
        else:
            res["sense_unexpected"].append([codon, exp, got])
    return res


def grade_control(inferred: str, control: dict) -> dict:
    """Grade one A8 control.

    control: {"table": <NCBI table number>, "must": {"CTG": "A", ...}}.
    Pass iff every "must" codon is called as stated, no sense codon is called
    against the control's own table (AGA/AGG->K artefacts excepted), and no stop
    codon of that table is called as an amino acid.
    """
    expected = TABLES[int(control["table"])]
    cmp = compare(inferred, expected)
    must = control.get("must") or {}
    must_ok = {c: inferred[codon_index(c)] == aa for c, aa in must.items()}
    passed = all(must_ok.values()) and not cmp["sense_unexpected"] and not cmp["stop_called"]
    return {
        "pass": passed,
        "must": {c: {"wanted": must[c], "called": inferred[codon_index(c)], "ok": ok}
                 for c, ok in must_ok.items()},
        "n_sense_called": len(cmp["sense_expected"]) + len(cmp["sense_unexpected"]) + len(cmp["artefacts"]),
        "comparison": cmp,
    }


def is_organelle(header: str) -> bool:
    """True for FASTA headers of organelle sequences (excluded before KACI runs)."""
    low = header.lower()
    return any(word in low for word in ORGANELLE_WORDS)


def filter_fasta(in_path: str, out_path: str) -> dict:
    """Copy a FASTA file without its organelle records. Returns counts."""
    kept = dropped = 0
    dropped_names = []
    keep = False
    with open(in_path, "r", encoding="ascii", errors="replace") as src, \
            open(out_path, "w", encoding="ascii", newline="\n") as dst:
        for line in src:
            if line.startswith(">"):
                keep = not is_organelle(line)
                if keep:
                    kept += 1
                else:
                    dropped += 1
                    dropped_names.append(line[1:].strip()[:120])
            if keep:
                dst.write(line if line.endswith("\n") else line + "\n")
    return {"kept": kept, "dropped": dropped, "dropped_names": dropped_names}
