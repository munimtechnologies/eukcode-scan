# eukcode-scan: grading of the transcriptome control pilot (P5, P6) and of prediction P7.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Pure functions for PREREGISTRATION.md addenda A2 and A3.

- control_list: flatten tsa_controls.json into one dict per control.
- evaluate: compare one Codetta inferred code with the control's registered table.
- grade_p5 / grade_p6: sensitivity and specificity of Codetta on transcriptome samples.
- gate: may the TSA screen start (needs P1, P2 from the genome pilot and P5, P6)?
- grade_p7: the biological prediction on tier-`new` non-ciliate species.

Counting rules (A3.2, fixed before any transcriptome output exists): a control "resolves" when
an INSDC TSA set was found for it. A resolved positive control that ends without a Codetta
result counts as a failure for P5. P6 counts only standard controls that finished, and needs at
least 3 of them.
"""

from .codes import called_count, compare, diagnostic_check

MIN_RESOLVED = 3
P7_THRESHOLD = 0.8


def control_list(cfg):
    """[{'label', 'role', 'taxa', 'table', 'diagnostic'}] from the tsa_controls.json dict."""
    out = []
    for role in ("positive", "standard"):
        for c in cfg.get(role, []):
            out.append({"label": c["label"], "role": role, "taxa": list(c["taxa"]),
                        "table": int(c["table"]), "diagnostic": dict(c.get("diagnostic", {}))})
    labels = [c["label"] for c in out]
    if len(labels) != len(set(labels)):
        raise ValueError("duplicate control label in tsa_controls.json")
    for c in out:
        if c["role"] == "positive" and not c["diagnostic"]:
            raise ValueError("positive control %s has no diagnostic codon" % c["label"])
    return out


def evaluate(control, inferred):
    """Verdict for one finished control: disagreements with its registered table, diagnostics."""
    ok, diag = diagnostic_check(inferred, control.get("diagnostic", {}))
    return {
        "label": control["label"],
        "role": control["role"],
        "table": control["table"],
        "n_called": called_count(inferred),
        "disagreements": compare(inferred, control["table"]),
        "diagnostic_ok": ok,
        "diagnostic": diag,
    }


def grade_p5(status, min_resolved=MIN_RESOLVED):
    """status: [{'label', 'role', 'resolved': bool, 'verdict': dict or None}]."""
    pos = [s for s in status if s["role"] == "positive" and s.get("resolved")]
    passed = [s["label"] for s in pos if s.get("verdict") and s["verdict"]["diagnostic_ok"]]
    unfinished = [s["label"] for s in pos if not s.get("verdict")]
    failed = [s["label"] for s in pos if s.get("verdict") and not s["verdict"]["diagnostic_ok"]]
    testable = len(pos) >= min_resolved
    return {
        "rule": "at least %d positive controls resolve, and all but at most one of them call every "
                "diagnostic codon as expected (unfinished counts as a failure)" % min_resolved,
        "resolved": len(pos), "passed": len(passed), "failed_labels": failed,
        "unfinished_labels": unfinished, "testable": testable,
        "holds": testable and len(passed) >= len(pos) - 1,
    }


def grade_p6(status, min_finished=MIN_RESOLVED):
    std = [s for s in status if s["role"] == "standard" and s.get("resolved")]
    done = [s for s in std if s.get("verdict")]
    bad = sum(len(s["verdict"]["disagreements"]) for s in done)
    testable = len(done) >= min_finished
    return {
        "rule": "at least %d standard-code controls finish, and across all of them at most one called "
                "codon contradicts table 1" % min_finished,
        "resolved": len(std), "finished": len(done), "disagreements": bad,
        "details": {s["label"]: s["verdict"]["disagreements"] for s in done if s["verdict"]["disagreements"]},
        "unfinished_labels": [s["label"] for s in std if not s.get("verdict")],
        "testable": testable,
        "holds": testable and bad <= 1,
    }


def gate(p5, p6, p1_holds=None, p2_holds=None):
    """The TSA screen may start only if P1, P2 (genome pilot) and P5, P6 all hold."""
    reasons = []
    for name, value in (("P1", p1_holds), ("P2", p2_holds)):
        if value is None:
            reasons.append("%s not graded yet (genome pilot has no result)" % name)
        elif not value:
            reasons.append("%s does not hold" % name)
    for name, g in (("P5", p5), ("P6", p6)):
        if not g["testable"]:
            reasons.append("%s cannot be tested (too few controls)" % name)
        elif not g["holds"]:
            reasons.append("%s does not hold" % name)
    return {"tsa_screen_may_start": not reasons, "reasons": reasons}


def is_ciliate(lineage):
    return "Ciliophora" in [x.strip() for x in (lineage or "").split(";")]


def grade_p7(records, threshold=P7_THRESHOLD):
    """records: per-set JSON records from run_tsa.py screen."""
    eligible = [r for r in records
                if r.get("tier") == "new" and "error" not in r and not is_ciliate(r.get("lineage", ""))]
    flagged = [r.get("accession", "") for r in eligible if r.get("candidates")]
    clean = len(eligible) - len(flagged)
    frac = clean / len(eligible) if eligible else None
    return {
        "rule": "in tier new, at least %d%% of non-ciliate species that finish have no candidate" % round(threshold * 100),
        "eligible": len(eligible), "no_candidate": clean,
        "fraction": None if frac is None else round(frac, 4),
        "with_candidates": flagged, "testable": bool(eligible),
        "holds": frac is not None and frac >= threshold,
    }


def report_md(status, p5, p6, g, meta):
    """Markdown summary of the transcriptome control pilot."""
    lines = ["# eukcode-scan transcriptome control pilot", "",
             "eukcode-scan %s, Codetta %s, cap %s Mb, finished %s UTC." % (
                 meta.get("eukcode_version", "?"), meta.get("codetta_commit", "?"), meta.get("cap_mb", "?"),
                 meta.get("finished_utc", "?")), "",
             "| label | role | TSA set | organism | NCBI code | registered table | called | diagnostic | disagreements | note |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for s in status:
        v = s.get("verdict")
        if v:
            diag = ", ".join("%s=%s(%s)" % (k, d["called"], "ok" if d["ok"] else "expected " + d["expected"])
                             for k, d in sorted(v["diagnostic"].items())) or "-"
            dis = ", ".join("%s %s->%s" % (d["codon"], d["expected"], d["called"]) for d in v["disagreements"]) or "none"
            called = str(v["n_called"])
        else:
            diag, dis, called = "-", "-", "-"
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            s["label"], s["role"], s.get("accession") or "-", s.get("organism") or "-",
            s.get("ncbi_gencode") or "-", s.get("table", "-"), called, diag, dis, (s.get("note") or "")[:90]))
    lines += ["", "## Pre-registered predictions", ""]
    for name, p in (("P5", p5), ("P6", p6)):
        verdict = "HOLDS" if p["holds"] else ("CANNOT BE TESTED" if not p["testable"] else "DOES NOT HOLD")
        extra = {k: v for k, v in p.items() if k not in ("rule", "holds", "testable")}
        lines.append("- **%s** (%s): %s. %s" % (name, p["rule"], verdict, extra))
    lines += ["", "**Gate for the TSA screen:** %s%s" % (
        "open" if g["tsa_screen_may_start"] else "closed",
        "" if g["tsa_screen_may_start"] else " (" + "; ".join(g["reasons"]) + ")"), ""]
    return "\n".join(lines) + "\n"
