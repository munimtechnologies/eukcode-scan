# eukcode-scan: prediction P11, the P10K CUG question in ciliates (PREREGISTRATION.md, A5).
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Do any ciliates read CUG as serine?

The P10K database paper (Gao et al., NAR 2024, 52:D747) assigns NCBI table 12 (the alternative
yeast code, CUG = Ser) to some Spirotrichea. No sense-codon reassignment has ever been reported
in a ciliate. Addendum A5 registers prediction P11 before any eukcode output exists:

  P11: no ciliate in the screen has a Codetta CTG = S call that survives check 2 (pass or
  pass_unplaced). A CTG = S call whose check 2 fails or is inconclusive counts as explained
  (or unresolved), not as a reassignment.

Pure functions only; inputs are the per-assembly JSON records of run_screen.py / run_tsa.py and
the check2/<case>.json results of run_check2.py.
"""

from .codes import CODON_INDEX, NONCALL

CTG = CODON_INDEX["CTG"]
SURVIVES = ("pass", "pass_unplaced")


def is_ciliate(lineage):
    return "Ciliophora" in [x.strip() for x in (lineage or "").split(";")]


def genus_of(organism):
    parts = (organism or "").strip().split()
    return parts[0] if parts else ""


def ctg_rows(records):
    """One row per ciliate record that finished: accession, organism, genus, CTG call."""
    out = []
    for r in records:
        if "error" in r or not is_ciliate(r.get("lineage", "")):
            continue
        inferred = r.get("inferred") or ""
        if len(inferred) != 64:
            continue
        out.append({"accession": r.get("accession", ""), "organism": r.get("organism", ""),
                    "genus": genus_of(r.get("organism", "")), "lineage": r.get("lineage", ""),
                    "data_type": r.get("data_type", "genome"), "ctg_call": inferred[CTG].upper()})
    return out


def call_counts(rows):
    out = {}
    for r in rows:
        k = r["ctg_call"] if r["ctg_call"] in ("L", "S", NONCALL) else "other"
        out[k] = out.get(k, 0) + 1
    return out


def ctg_check2(result):
    """Check 2 outcome for CTG from one run_check2 case result, or None if absent."""
    if not result or "error" in result:
        return None
    entry = (result.get("codons") or {}).get("CTG")
    if not entry:
        return None
    return entry.get("check2", {}).get("outcome")


def grade_p11(rows, check2_by_accession):
    """rows: ctg_rows(); check2_by_accession: {accession: run_check2 case result}.

    Graded only when every CTG = S ciliate has a check 2 result for CTG.
    """
    flagged = [r for r in rows if r["ctg_call"] == "S"]
    survived, explained, pending = [], [], []
    for r in flagged:
        outcome = ctg_check2(check2_by_accession.get(r["accession"]))
        if outcome is None:
            pending.append(r["accession"])
        elif outcome in SURVIVES:
            survived.append({"accession": r["accession"], "organism": r["organism"], "check2": outcome})
        else:
            explained.append({"accession": r["accession"], "organism": r["organism"], "check2": outcome})
    graded = bool(rows) and not pending
    return {
        "rule": "no ciliate has a Codetta CTG = S call that survives check 2 (pass or pass_unplaced)",
        "ciliates_screened": len(rows), "ctg_calls": call_counts(rows),
        "ctg_s": [r["accession"] for r in flagged], "survived_check2": survived,
        "explained_or_unresolved": explained, "pending_check2": pending,
        "graded": graded, "holds": graded and not survived,
    }


def p10k_targets(targets, rows):
    """For each genus P10K labels with table 12: the screened datasets and their CTG calls.

    A CTG call of L contradicts P10K's table 12 for that dataset; '?' means Codetta saw too
    little evidence to call CTG; S agrees with P10K and goes to check 2.
    """
    by_genus = {}
    for r in rows:
        by_genus.setdefault(r["genus"], []).append(r)
    out = []
    for g in targets:
        hits = by_genus.get(g, [])
        calls = sorted(set(h["ctg_call"] for h in hits))
        if not hits:
            status = "not_screened"
        elif calls == ["S"]:
            status = "agrees_with_p10k"
        elif "S" in calls:
            status = "mixed"
        elif "L" in calls:
            status = "contradicts_p10k"
        else:
            status = "no_call"
        out.append({"genus": g, "datasets": [h["accession"] for h in hits], "ctg_calls": calls, "status": status})
    return out
