#!/usr/bin/env python3
# eukcode-scan A7: tRNA test of P10K's Euplotes CUG label with tRNAscan-SE 2.0 isotype models.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Run tRNAscan-SE 2.0 (--detail) on the A7 inputs and grade P12-P15 (PREREGISTRATION_A7.md).

Usage: python run_trna.py --cases trna_cases.json --out-dir /out [--threads 2]

Each input is downloaded to temporary storage, its contigs are renamed c1, c2, ... (the map
and per-contig length, GC and Euplotes-telomere flags are written next to the results), it is
scanned, and the temporary copy is deleted. Nothing is cached in the shared data volume.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from eukcode import ncbi, trna, tsa  # noqa: E402

TOOLS_BIN = os.path.join(HERE, "trna-tools", "bin")

# Human tRNA-Leu-CAG-1-1 (GtRNAdb hg38: Leu 124.9 bits, Ser 82.0 bits), flanked by filler.
SMOKE_TRNA = ("GTCAGGATGGCCGAGCGGTCTAAGGCGCTGCGTTCAGGTCGCAGTCTCCCCTGGAGGCGTGGGTTCGAATCCCACTCCTGACA")
SMOKE_FLANK = "ATTTAGCATAATCGATTAACTTAGTCAATACTTAAGCTATTAGAT"


def tool_env():
    env = dict(os.environ)
    env["PATH"] = TOOLS_BIN + os.pathsep + env.get("PATH", "")
    return env


def trnascan_path():
    p = os.path.join(TOOLS_BIN, "tRNAscan-SE")
    return p if os.path.exists(p) else (shutil.which("tRNAscan-SE") or "tRNAscan-SE")


def run_trnascan(fasta, prefix, threads, timeout_s):
    """Run tRNAscan-SE 2.0 in eukaryotic mode with isotype detail; return (paths, seconds, rc)."""
    paths = {k: "%s.%s" % (prefix, k) for k in ("tab", "struct", "iso", "stats", "log")}
    for p in paths.values():
        if os.path.exists(p):
            os.remove(p)
    cmd = [trnascan_path(), "-E", "--detail", "-q", "--thread", str(threads),
           "-o", paths["tab"], "-f", paths["struct"], "-s", paths["iso"], "-m", paths["stats"],
           fasta]
    t0 = time.time()
    with open(paths["log"], "w") as log:
        try:
            rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                env=tool_env(), timeout=timeout_s).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return paths, round(time.time() - t0, 1), rc


def read_text(path):
    if not os.path.exists(path):
        return ""
    with open(path, "r", encoding="ascii", errors="replace") as fh:
        return fh.read()


def smoke_test(out_dir, threads):
    """Install check: the human Leu-CAG tRNA must be found with anticodon CAG."""
    d = os.path.join(out_dir, "smoke")
    os.makedirs(d, exist_ok=True)
    fa = os.path.join(d, "smoke.fa")
    with open(fa, "w") as fh:
        fh.write(">smoke\n%s%s%s\n" % (SMOKE_FLANK, SMOKE_TRNA, SMOKE_FLANK))
    paths, secs, rc = run_trnascan(fa, os.path.join(d, "smoke"), threads, 600)
    rows = trna.parse_table(read_text(paths["tab"]))
    found = [r for r in rows if r["anticodon"] == "CAG"]
    result = {"rc": rc, "seconds": secs, "rows": rows,
              "ok": bool(found), "isotype_cm": found[0]["isotype_cm"] if found else None}
    with open(os.path.join(out_dir, "smoke.json"), "w") as fh:
        json.dump(result, fh, indent=1)
    print("smoke test:", json.dumps({k: result[k] for k in ("rc", "seconds", "ok", "isotype_cm")}))
    if not found:
        print(read_text(paths["log"])[-3000:])
    return result


def fetch(case, tmp):
    """Download one input to tmp; return a FASTA or FASTA.gz path."""
    if case["kind"] == "genome":
        zip_path = ncbi.download_genome_zip(case["accession"], tmp)
        fa = os.path.join(tmp, case["label"] + ".raw.fna")
        ncbi.extract_fasta(zip_path, fa)
        os.remove(zip_path)
        return fa
    if case["kind"] == "tsa":
        last = None
        for url in tsa.fasta_urls(case["accession"]):
            path = os.path.join(tmp, case["label"] + ".raw.fa.gz")
            try:
                ncbi.http_download(url, path)
                return path
            except Exception as exc:  # try the next mirror
                last = exc
        raise last
    raise ValueError("unknown kind %r" % case["kind"])


MAP_HEADER = "short\toriginal\tlength\tgc\ttel5\ttel3\n"


def map_line(short, st):
    return "%s\t%s\t%d\t%.4f\t%d\t%d\n" % (short, st["original"], st["length"], st["gc"],
                                           st["tel5"], st["tel3"])


def rename(src, fasta_out):
    """Write contigs as c1, c2, ...; return ({short: stats}, totals)."""
    contigs = {}
    n = gc = acgt = 0
    with open(fasta_out, "w") as fo:
        for header, seq in tsa.iter_records(src):
            n += 1
            short = "c%d" % n
            st = trna.contig_stats(seq)
            st["original"] = header.split()[0] if header.split() else header
            contigs[short] = st
            gc += seq.count("G") + seq.count("C")
            acgt += sum(seq.count(b) for b in "ACGT")
            fo.write(">%s\n" % short)
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i + 80] + "\n")
    tel = sum(1 for c in contigs.values() if c["tel5"] or c["tel3"])
    return contigs, {"contigs": n, "bp": sum(c["length"] for c in contigs.values()),
                     "gc": round(gc / acgt, 4) if acgt else 0.0, "telomere_contigs": tel}


def iso_lines(path, ids):
    keep = []
    want = set(ids)
    for line in read_text(path).splitlines():
        first = line.split("\t")[0].strip() if line else ""
        if first in want or any(line.startswith(i + " ") or line.startswith(i + "\t") for i in want):
            keep.append(line)
    return keep


def run_case(case, out_dir, tmp, threads, timeout_s):
    label = case["label"]
    rec = {"label": label, "role": case["role"], "kind": case["kind"],
           "accession": case["accession"], "organism": case.get("organism", ""), "summary": None}
    t0 = time.time()
    try:
        raw = fetch(case, tmp)
        rec["download_seconds"] = round(time.time() - t0, 1)
        fa = os.path.join(tmp, label + ".fa")
        contigs, totals = rename(raw, fa)
        os.remove(raw)
        rec["input"] = totals
        paths, secs, rc = run_trnascan(fa, os.path.join(out_dir, "raw", label), threads, timeout_s)
        os.remove(fa)
        rec["scan_seconds"], rec["rc"] = secs, rc
        rows = trna.parse_table(read_text(paths["tab"]))
        with open(os.path.join(out_dir, "raw", label + ".trna_contigs.tsv"), "w") as fm:
            fm.write(MAP_HEADER)
            for short in sorted({r["name"] for r in rows}, key=lambda x: int(x[1:]) if x[1:].isdigit() else 0):
                if short in contigs:
                    fm.write(map_line(short, contigs[short]))
        struct = trna.parse_struct(read_text(paths["struct"]))
        summary = trna.summarize(rows, contigs=contigs, struct=struct)
        summary["iso_lines"] = iso_lines(paths["iso"], [d["id"] for d in summary["decoders"]])
        if rc != 0:
            rec["error"] = "tRNAscan-SE exit %s" % rc
        if rc == 0 or rows:
            rec["summary"] = summary
    except Exception as exc:
        rec["error"] = "%s: %s" % (type(exc).__name__, str(exc)[:300])
    finally:
        for name in os.listdir(tmp):
            if name.startswith(label + ".") or name.startswith(case["accession"]):
                try:
                    os.remove(os.path.join(tmp, name))
                except OSError:
                    pass
    rec["seconds"] = round(time.time() - t0, 1)
    s = rec["summary"]
    print("%-10s %-6s trnas=%s cag=%s cag_leu=%s cag_ser=%s tag=%s rc=%s err=%s" % (
        label, case["role"], s and s["trnas"], s and s["cag"], s and s["cag_leu"],
        s and s["cag_ser"], s and s["tag"], rec.get("rc"), rec.get("error", "")), flush=True)
    return rec


def fmt_decoder(d):
    return "| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
        d["id"], d.get("contig_original", ""), d["anticodon"], d["isotype_cm"] or "-",
        d["isotype_score"], d["inf_score"], d["note"] or "", d.get("telomere", ""),
        d.get("tail4", ""))


def write_report(out_dir, records, grade, smoke, cases_file):
    lines = ["# A7 tRNA test of P10K's Euplotes CUG label", "",
             "Pre-registered in PREREGISTRATION_A7.md (eukcode-scan v0.7.0) before any output. "
             "Cases: `%s`. Tool: tRNAscan-SE 2.0.12 (Infernal 1.1.5), eukaryotic mode, --detail. "
             "Code written by AI agents working in MV2's name; no person has checked these "
             "results yet." % os.path.basename(cases_file), "",
             "Install smoke test (human tRNA-Leu-CAG-1-1): found=%s, isotype=%s." % (
                 smoke.get("ok"), smoke.get("isotype_cm")), "",
             "## Grades", "",
             "- **P12** (CUG-Ser yeasts show a Ser-type CAG tRNA): **%s** %s" % (
                 grade["P12"], json.dumps(grade["P12_detail"])),
             "- **P13** (ciliate CUG-Leu control: CAG tRNAs Leu-type): **%s** %s" % (
                 grade["P13"], json.dumps(grade["P13_detail"])),
             "- **P14** (E. focardii TN1 genome outcome): **%s**" % json.dumps(grade["P14"]),
             "- **P15** (other Euplotes genomes): %s" % json.dumps(grade["P15"]),
             "- Transcriptomes (secondary, A7.5): %s" % json.dumps(grade["transcriptomes"]), ""]
    for n in grade["notes"]:
        lines.append("> " + n)
        lines.append("")
    lines += ["## Per input", "",
              "| label | role | accession | contigs | Mb | GC | telomere contigs | tRNAs | pseudo | "
              "IPD | CAG (Leu/Ser/other) | TAG (Leu/Ser) | minutes | error |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in records:
        s, t = r.get("summary"), r.get("input", {})
        lines.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            r["label"], r["role"], r["accession"], t.get("contigs", ""),
            round(t.get("bp", 0) / 1e6, 1) if t else "", t.get("gc", ""),
            t.get("telomere_contigs", ""), s and s["trnas"], s and s["pseudo"], s and s["ipd"],
            s and "%d (%d/%d/%d)" % (s["cag"], s["cag_leu"], s["cag_ser"], s["cag_other"]),
            s and "%d (%d/%d)" % (s["tag"], s["tag_leu"], s["tag_ser"]),
            round(r.get("seconds", 0) / 60, 1), r.get("error", "")))
    lines += ["", "## Every CUG-reading tRNA (anticodon CAG or TAG, pseudogenes included)", ""]
    for r in records:
        s = r.get("summary")
        if not s:
            continue
        lines += ["### %s (%s, %s)" % (r["label"], r.get("organism", ""), r["accession"]), "",
                  "| id | contig | anticodon | isotype CM | isotype score | Inf score | note | "
                  "telomere | 3' end |", "|---|---|---|---|---|---|---|---|---|"]
        lines += [fmt_decoder(d) for d in s["decoders"]] or ["| none | | | | | | | | |"]
        lines.append("")
    with open(os.path.join(out_dir, "TRNA_RESULT.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="trna_cases.json")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--tmp", default="/tmp/eukcode-trna")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--per-case-max-min", type=float, default=90)
    a = ap.parse_args(argv)
    os.makedirs(os.path.join(a.out_dir, "raw"), exist_ok=True)
    os.makedirs(a.tmp, exist_ok=True)
    with open(a.cases) as fh:
        cases = json.load(fh)["cases"]
    smoke = smoke_test(a.out_dir, a.threads)
    if not smoke["ok"]:
        print("tRNAscan-SE install check failed; stopping before any input is read")
        return 2
    records = []
    for case in cases:
        records.append(run_case(case, a.out_dir, a.tmp, a.threads, int(a.per_case_max_min * 60)))
        with open(os.path.join(a.out_dir, "trna_results.json"), "w") as fh:
            json.dump({"cases_file": os.path.basename(a.cases), "smoke": smoke,
                       "records": records}, fh, indent=1)
    grade = trna.grade_a7({r["label"]: r for r in records})
    with open(os.path.join(a.out_dir, "trna_results.json"), "w") as fh:
        json.dump({"cases_file": os.path.basename(a.cases), "smoke": smoke, "grade": grade,
                   "records": records}, fh, indent=1)
    write_report(a.out_dir, records, grade, smoke, a.cases)
    print(json.dumps(grade, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
