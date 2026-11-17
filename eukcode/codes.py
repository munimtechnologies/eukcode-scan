# eukcode-scan: NCBI genetic code tables and comparison helpers.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Genetic code tables in Codetta / NCBI codon order (bases T, C, A, G; first base slowest).

Only nuclear (and bacterial) tables are listed: the screen is about nuclear codes.
Tables are built from the standard code plus the NCBI-documented changes, so a typo in a
64-character string cannot slip in.
"""

BASES = "TCAG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODON_INDEX = {c: i for i, c in enumerate(CODONS)}
STANDARD = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

# NCBI translation table id -> {codon: amino acid} changes relative to table 1.
# Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
_CHANGES = {
    1: {},
    4: {"TGA": "W"},                                  # Mycoplasma/Spiroplasma (nuclear-type use here)
    6: {"TAA": "Q", "TAG": "Q"},                      # ciliate, dasycladacean, Hexamita
    10: {"TGA": "C"},                                 # Euplotid
    11: {},                                           # bacterial/archaeal/plastid (same codon meanings)
    12: {"CTG": "S"},                                 # alternative yeast
    15: {"TAG": "Q"},                                 # Blepharisma
    25: {"TGA": "G"},                                 # candidate division SR1 / Gracilibacteria
    26: {"CTG": "A"},                                 # Pachysolen tannophilus
    27: {"TAA": "Q", "TAG": "Q", "TGA": "W"},         # karyorelict (TGA also stop)
    28: {"TAA": "Q", "TAG": "Q", "TGA": "W"},         # Condylostoma (all three also stop)
    29: {"TAA": "Y", "TAG": "Y"},                     # Mesodinium
    30: {"TAA": "E", "TAG": "E"},                     # Peritrich
    31: {"TAA": "E", "TAG": "E", "TGA": "W"},         # Blastocrithidia (TAA also stop)
}

# Codons that NCBI documents as having a dual meaning (sense or stop) in a table.
# A call of either the amino acid or "no call" is acceptable there.
DUAL_MEANING = {
    27: {"TGA"},
    28: {"TAA", "TAG", "TGA"},
    31: {"TAA"},
}

NONCALL = "?"


def table(table_id):
    """Return the 64-character code string for an NCBI table id."""
    if table_id not in _CHANGES:
        raise KeyError("translation table %r is not in eukcode.codes" % (table_id,))
    code = list(STANDARD)
    for codon, aa in _CHANGES[table_id].items():
        code[CODON_INDEX[codon]] = aa
    return "".join(code)


def known_tables():
    return sorted(_CHANGES)


def normalize_codon(codon):
    codon = codon.strip().upper().replace("U", "T")
    if codon not in CODON_INDEX:
        raise ValueError("not a codon: %r" % (codon,))
    return codon


def compare(inferred, table_id):
    """Compare a Codetta 64-char inferred code with an NCBI table.

    Returns a list of disagreements, one dict per called codon whose call differs
    from the table. '?' (no call) never counts as a disagreement. Codons with a
    documented dual meaning (sense or stop) carry the sense amino acid in the table,
    so calling that amino acid, or making no call, is not a disagreement.
    """
    if len(inferred) != 64:
        raise ValueError("inferred code must have 64 characters, got %d" % len(inferred))
    expected = table(table_id)
    out = []
    for i, codon in enumerate(CODONS):
        call = inferred[i].upper()
        if call == NONCALL:
            continue
        exp = expected[i]
        if call == exp:
            continue
        kind = "stop_to_sense" if exp == "*" else "sense_to_sense"
        out.append({"codon": codon, "called": call, "expected": exp, "kind": kind})
    return out


def diagnostic_check(inferred, diagnostic):
    """Check named codon calls, e.g. {'CTG': 'S'}. Returns (all_ok, details)."""
    details = {}
    ok = True
    for codon, aa in diagnostic.items():
        c = normalize_codon(codon)
        call = inferred[CODON_INDEX[c]].upper()
        details[c] = {"expected": aa, "called": call, "ok": call == aa}
        ok = ok and call == aa
    return ok, details


def called_count(inferred):
    return sum(1 for ch in inferred if ch != NONCALL)
