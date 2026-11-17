# eukcode-scan: grading of the pre-registered pilot predictions.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Pure functions that turn per-control results into pass/fail for P1-P4 (PREREGISTRATION.md)."""

import statistics

from .codes import compare, diagnostic_check


def evaluate_control(control, inferred):
    """Return the per-control verdict dict for one finished Codetta run."""
    disagreements = compare(inferred, control["table"])
    diag_ok, diag = diagnostic_check(inferred, control.get("diagnostic", {}))
    return {
        "label": control["label"],
        "table": control["table"],
        "role": control.get("role", "main"),
        "n_called": sum(1 for ch in inferred if ch != "?"),
        "disagreements": disagreements,
        "diagnostic_ok": diag_ok,
        "diagnostic": diag,
    }


def grade(verdicts, runtimes_min_per_mb):
    """verdicts: list from evaluate_control (only finished runs). Returns dict P1..P4."""
    main = [v for v in verdicts if v["role"] == "main"]
    positives = [v for v in main if v["diagnostic"]]
    standard = [v for v in main if v["table"] == 1]
    power = [v for v in verdicts if v["role"] == "power"]

    n_pos = len(positives)
    n_pos_ok = sum(1 for v in positives if v["diagnostic_ok"])
    p1 = {
        "rule": "positive controls with every diagnostic codon called as expected >= n_resolved - 1",
        "resolved": n_pos, "passed": n_pos_ok,
        "holds": n_pos > 0 and n_pos_ok >= n_pos - 1,
        "failed_labels": [v["label"] for v in positives if not v["diagnostic_ok"]],
    }

    std_bad = sum(len(v["disagreements"]) for v in standard)
    all_bad = sum(len(v["disagreements"]) for v in main)
    p2 = {
        "rule": "0 disagreements in standard-code controls and <= 1 disagreement over all main controls",
        "standard_disagreements": std_bad, "all_disagreements": all_bad,
        "holds": len(standard) > 0 and std_bad == 0 and all_bad <= 1,
        "details": {v["label"]: v["disagreements"] for v in main if v["disagreements"]},
    }

    p3 = {
        "rule": "2 Mb samples still call every diagnostic codon (calb_2mb CTG=S; tthe_2mb TAA=Q, TAG=Q)",
        "resolved": len(power), "passed": sum(1 for v in power if v["diagnostic_ok"]),
        "holds": len(power) > 0 and all(v["diagnostic_ok"] for v in power),
    }

    med = statistics.median(runtimes_min_per_mb) if runtimes_min_per_mb else None
    p4 = {
        "rule": "median Codetta wall time <= 4 minutes per Mb of input at 2 hmmscan threads",
        "median_min_per_mb": med,
        "holds": med is not None and med <= 4.0,
    }
    return {"P1": p1, "P2": p2, "P3": p3, "P4": p4}
