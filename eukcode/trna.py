# eukcode-scan: tRNAscan-SE 2.0 output parsing and the A7 tRNA test (pure functions).
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""The A7 tRNA test of P10K's Euplotes CUG label (PREREGISTRATION_A7.md).

A nuclear CUG -> Ser reassignment, as in the CTG-clade yeasts, is carried by a tRNA whose
anticodon reads CUG (5'-CAG-3') but whose body is serine-type. tRNAscan-SE 2.0 scores every
predicted tRNA against isotype-specific covariance models independently of its anticodon
(--detail: "Isotype CM", "Isotype Score", note "IPD" on disagreement). This module parses
that output and applies the pre-registered reading rules. It has no network or tool calls.
"""

import re

# Anticodons (DNA, 5'->3') that can read the codon CUG: CAG pairs with it directly, and a U34
# tRNA (anticodon TAG) can read all four CUN codons.
CUG_ANTICODONS = ("CAG", "TAG")

# Euplotes macronuclear telomeres are C4A4 repeats (G-strand T4G4).
TELOMERE_MOTIFS = ("CCCCAAAACCCC", "GGGGTTTTGGGG")
TELOMERE_WINDOW = 100


def _int(x, default=0):
    try:
        return int(str(x).strip())
    except (TypeError, ValueError):
        return default


def _float(x):
    try:
        return float(str(x).strip())
    except (TypeError, ValueError):
        return None


def _dna(s):
    return (s or "").strip().upper().replace("U", "T")


def parse_row(fields, detail):
    """One tab-split data row of a tRNAscan-SE 2.0 -o table -> dict, or None if malformed."""
    f = [x.strip() for x in fields]
    if len(f) < 9 or not f[0]:
        return None
    num = _int(f[1], default=-1)
    if num < 0:
        return None
    row = {
        "name": f[0],
        "num": num,
        "id": "%s.trna%d" % (f[0], num),
        "begin": _int(f[2]),
        "end": _int(f[3]),
        "type": f[4],
        "anticodon": _dna(f[5]),
        "intron_begin": _int(f[6]),
        "intron_end": _int(f[7]),
        "inf_score": _float(f[8]),
        "hmm_score": None,
        "str_score": None,
        "isotype_cm": "",
        "isotype_score": None,
        "note": "",
    }
    if detail and len(f) >= 13:
        row["hmm_score"] = _float(f[9])
        row["str_score"] = _float(f[10])
        row["isotype_cm"] = f[11]
        row["isotype_score"] = _float(f[12])
        row["note"] = " ".join(x for x in f[13:] if x)
    else:
        row["note"] = " ".join(x for x in f[9:] if x)
    row["strand"] = "-" if row["begin"] > row["end"] else "+"
    return row


def parse_table(text):
    """Parse a tRNAscan-SE 2.0 tabular (-o) file into a list of row dicts.

    The header ends with a line of dashes. With --detail its columns are: name, tRNA #,
    begin, end, type, anticodon, intron begin, intron end, Inf score, HMM score, 2'Str score,
    isotype CM, isotype score, note. Without a header (--brief), a row of 13 or more fields
    is read as --detail.
    """
    lines = text.splitlines()
    dash = None
    for i, line in enumerate(lines):
        if line.startswith("--------"):
            dash = i
            break
    if dash is None:
        header, body, detail = [], lines, None
    else:
        header, body = lines[:dash], lines[dash + 1:]
        detail = any("Isotype" in h for h in header)
    rows = []
    for line in body:
        if not line.strip():
            continue
        fields = line.split("\t")
        d = detail if detail is not None else len(fields) >= 13
        row = parse_row(fields, d)
        if row is not None:
            rows.append(row)
    return rows


_STRUCT_HEAD = re.compile(r"^(\S+\.trna\d+)\s+\((\d+)-(\d+)\)\s+Length:\s*(\d+)")


def parse_struct(text):
    """Parse a tRNAscan-SE secondary-structure (-f) file into {tRNA id: {seq, str, length}}."""
    out = {}
    cur = None
    for line in text.splitlines():
        m = _STRUCT_HEAD.match(line)
        if m:
            cur = m.group(1)
            out[cur] = {"length": int(m.group(4)), "seq": "", "str": ""}
            continue
        if cur is None:
            continue
        if line.startswith("Seq: "):
            out[cur]["seq"] = line[5:].strip()
        elif line.startswith("Str: "):
            out[cur]["str"] = line[5:].strip()
    return out


def contig_stats(seq, window=TELOMERE_WINDOW):
    """Length, GC fraction (of A/C/G/T) and a C4A4 telomere flag at each end of one contig."""
    s = (seq or "").upper()
    acgt = sum(s.count(b) for b in "ACGT")
    gc = (s.count("G") + s.count("C")) / acgt if acgt else 0.0
    head, tail = s[:window], s[-window:] if s else ""
    return {
        "length": len(s),
        "gc": round(gc, 4),
        "tel5": any(m in head for m in TELOMERE_MOTIFS),
        "tel3": any(m in tail for m in TELOMERE_MOTIFS),
    }


def is_pseudo(row):
    return "pseudo" in (row.get("note") or "").lower()


def body_class(row):
    iso = (row.get("isotype_cm") or "").strip()
    if not iso:
        return "no-isotype"
    if iso == "Leu":
        return "Leu-type"
    if iso == "Ser":
        return "Ser-type"
    return "other:" + iso


def summarize(rows, contigs=None, struct=None):
    """Counts and the CUG-decoding tRNAs of one input (pseudogenes counted apart)."""
    contigs = contigs or {}
    struct = struct or {}
    real = [r for r in rows if not is_pseudo(r)]
    isotypes = {}
    for r in real:
        isotypes[r["type"]] = isotypes.get(r["type"], 0) + 1
    decoders = []
    for r in rows:
        if r["anticodon"] not in CUG_ANTICODONS:
            continue
        d = {
            "id": r["id"],
            "contig": r["name"],
            "begin": r["begin"],
            "end": r["end"],
            "type": r["type"],
            "anticodon": r["anticodon"],
            "isotype_cm": r["isotype_cm"],
            "isotype_score": r["isotype_score"],
            "inf_score": r["inf_score"],
            "note": r["note"],
            "pseudo": is_pseudo(r),
            "body": body_class(r),
        }
        c = contigs.get(r["name"])
        if c:
            d["contig_original"] = c.get("original", "")
            d["contig_length"] = c.get("length")
            d["contig_gc"] = c.get("gc")
            d["telomere"] = bool(c.get("tel5") or c.get("tel3"))
        s = struct.get(r["id"])
        if s:
            d["seq"] = s.get("seq", "")
            d["str"] = s.get("str", "")
            d["tail4"] = s.get("seq", "")[-4:]
        decoders.append(d)
    live = [d for d in decoders if not d["pseudo"]]
    cag = [d for d in live if d["anticodon"] == "CAG"]
    tag = [d for d in live if d["anticodon"] == "TAG"]
    return {
        "trnas": len(real),
        "pseudo": len(rows) - len(real),
        "ipd": sum(1 for r in real if "IPD" in (r.get("note") or "")),
        "isotype_counts": dict(sorted(isotypes.items())),
        "cag": len(cag),
        "cag_leu": sum(1 for d in cag if d["body"] == "Leu-type"),
        "cag_ser": sum(1 for d in cag if d["body"] == "Ser-type"),
        "cag_ser_telomere": sum(1 for d in cag if d["body"] == "Ser-type" and d.get("telomere")),
        "cag_other": sum(1 for d in cag if d["body"] not in ("Leu-type", "Ser-type")),
        "tag": len(tag),
        "tag_leu": sum(1 for d in tag if d["body"] == "Leu-type"),
        "tag_ser": sum(1 for d in tag if d["body"] == "Ser-type"),
        "decoders": decoders,
    }


def grade_positive(summary):
    """P12 control: a known CUG-Ser yeast must show >= 1 Ser-type CAG tRNA."""
    if summary is None:
        return "not run"
    return "pass" if summary["cag_ser"] >= 1 else "fail"


def grade_negative(summary):
    """P13 control: every non-pseudo CAG tRNA of a CUG-Leu ciliate is Leu-type (needs >= 1)."""
    if summary is None:
        return "not run"
    if summary["cag"] == 0:
        return "uninformative"
    return "pass" if summary["cag_leu"] == summary["cag"] else "fail"


def euplotes_outcome(summary):
    """A7.4 outcome class for one Euplotes genome (rules fixed before any output)."""
    if summary is None:
        return "not run"
    if summary["cag_ser_telomere"] >= 1:
        return "candidate"
    if summary["cag_ser"] >= 1:
        return "unattributed-ser-type"
    if summary["cag_leu"] >= 1:
        return "cug-leu-supported"
    if summary["cag"] == 0:
        if summary["tag_leu"] >= 1 and summary["tag_ser"] == 0:
            return "no-cag-trna-leu-u34"
        return "no-cug-decoder-found"
    return "cag-other-isotype"


def grade_a7(results):
    """Grade P12-P15 from {label: {"role": ..., "summary": dict or None}} (A7.3-A7.5)."""
    def summ(role):
        return [(lab, r.get("summary")) for lab, r in sorted(results.items()) if r.get("role") == role]

    pos = summ("positive")
    neg = summ("negative")
    main = summ("main")
    other = summ("euplotes")
    tsa = summ("transcriptome")

    pos_g = {lab: grade_positive(s) for lab, s in pos}
    if pos_g and all(g == "pass" for g in pos_g.values()):
        p12 = "pass"
    elif any(g == "fail" for g in pos_g.values()):
        p12 = "fail"
    else:
        p12 = "not run"
    neg_g = {lab: grade_negative(s) for lab, s in neg}
    vals = list(neg_g.values())
    if vals and all(g == "pass" for g in vals):
        p13 = "pass"
    elif any(g == "fail" for g in vals):
        p13 = "fail"
    elif vals and all(g == "uninformative" for g in vals):
        p13 = "uninformative"
    else:
        p13 = "not run"
    main_o = {lab: euplotes_outcome(s) for lab, s in main}
    other_o = {lab: euplotes_outcome(s) for lab, s in other}
    tsa_o = {lab: euplotes_outcome(s) for lab, s in tsa}

    notes = []
    if p12 != "pass":
        notes.append("P12 did not pass: Ser-type detection is not validated, so a result without "
                     "a Ser-type CAG tRNA does not count against CUG = Ser.")
    if p13 == "fail":
        notes.append("P13 failed: Leu-type calls are not reliable in a ciliate genome.")
    outs = set(main_o.values())
    if "candidate" in outs:
        notes.append("Candidate: a Ser-type CAG tRNA on a telomere-capped Euplotes contig. Checks "
                     "and a person's review come before anything is said; not a discovery.")
    elif p12 == "pass" and outs and outs <= {"cug-leu-supported", "no-cag-trna-leu-u34"}:
        notes.append("No Ser-type CUG-reading tRNA in the TN1 genome while the CUG-Ser controls show "
                     "theirs: P10K's table-12 label for Euplotes is not supported at the tRNA level.")
    if any(o in ("candidate", "unattributed-ser-type") for o in tsa_o.values()) and \
            not any(o in ("candidate", "unattributed-ser-type") for o in outs):
        notes.append("A Ser-type CAG tRNA appears in a TN1 transcriptome but not in the TN1 genome: "
                     "consistent with a contaminant in that culture (A7.5).")
    return {
        "P12": p12, "P12_detail": pos_g,
        "P13": p13, "P13_detail": neg_g,
        "P14": main_o,
        "P15": other_o,
        "transcriptomes": tsa_o,
        "notes": notes,
    }
