#!/usr/bin/env python3
# eukcode-scan check 2: where does a candidate codon's evidence come from?
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Check 2 of PREREGISTRATION.md (addendum A4): contamination and source of the signal.

  python run_check2.py demo --cases check2_demo.json --data-dir /data --out-dir /out
      The pre-registered validation (P8-P10): a Euplotes transcriptome (UGA = Cys) and a
      constructed 50/50 mix of a UAR = Gln ciliate and a standard-code alga.

  python run_check2.py records --records-dir out-t00/species --data-dir /data --out-dir /out
      Every screened set (JSON record from run_screen.py or run_tsa.py) that has candidates:
      redraws the same sample, reruns Codetta keeping its alignment output, traces each
      candidate codon's observations to the input sequences, sends the top 20 sequences to
      blastx (EMBL-EBI, UniProtKB, the species' own taxid excluded) and grades check 2.

Writes CHECK2.md, check2_summary.json and check2/<label>.json.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from eukcode import __version__, codes, codetta_io, ebi_blast, fasta, ncbi, support, tsa  # noqa: E402

_lock = threading.Lock()
LOG_PATH = None


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    with _lock:
        print(line, flush=True)
        if LOG_PATH:
            with open(LOG_PATH, "a", encoding="utf-8") as fh:
                fh.write(line + "\n")


# ------------------------------------------------------------------ inputs

def resolve_taxid(name):
    for url in (tsa.taxid_search_url(name),
                "%s/esearch.fcgi?db=taxonomy&retmode=json&term=%s" % (tsa.EUTILS, urllib.parse.quote(name))):
        try:
            ids = json.loads(ncbi.http_get(url, pause=0.4).decode("utf-8")).get("esearchresult", {}).get("idlist", [])
        except Exception:
            continue
        if len(ids) == 1:
            return int(ids[0])
    return None


def resolve_tsa(taxa, min_bp=1_000_000):
    """First name in `taxa` with a usable TSA set -> its largest set's ENA row, or None."""
    for name in taxa:
        tid = resolve_taxid(name)
        if tid is None:
            log("  %s: taxid not resolved" % name)
            continue
        try:
            rows = tsa.parse_tsv(ncbi.http_get(tsa.ena_search_url(tid), timeout=300, pause=0.5).decode("utf-8"))
        except Exception as exc:
            log("  %s: ENA error %s" % (name, exc))
            continue
        picked, _sk = tsa.choose_sets(rows, min_bp=min_bp)
        if picked:
            best = max(picked, key=lambda r: tsa._int(r.get("base_count")))
            best["resolved_from"] = name
            return best
        log("  %s (taxid %d): no usable TSA set" % (name, tid))
    return None


def load_input(spec, work):
    """Download and sample one input. Returns (records, info)."""
    os.makedirs(work, exist_ok=True)
    cap = int(float(spec.get("cap_mb", 2.0)) * 1e6)
    seed, disjoint = int(spec.get("seed", 0)), spec.get("disjoint_from")
    kind = spec.get("kind", "tsa")
    info = {"kind": kind, "cap_mb": spec.get("cap_mb", 2.0), "seed": seed, "disjoint_from": disjoint,
            "prefix": spec.get("prefix", "")}
    if kind == "genome":
        acc = spec["accession"]
        zp = ncbi.download_genome_zip(acc, work)
        raw = os.path.join(work, "raw.fna")
        ncbi.extract_fasta(zp, raw)
        os.remove(zp)
        recs, n_total = fasta.sample_file(raw, cap, seed=seed, disjoint_from=disjoint)
        os.remove(raw)
        info.update({"accession": acc, "units_total": n_total})
    else:
        if spec.get("accession"):
            row = {"accession": spec["accession"], "scientific_name": spec.get("organism", ""),
                   "tax_id": spec.get("tax_id", "")}
        else:
            row = resolve_tsa(spec["taxa"])
            if row is None:
                raise RuntimeError("no usable TSA set for %s" % spec["taxa"])
        acc = row["accession"]
        gz = os.path.join(work, "set.fsa_nt.gz")
        last = None
        for u in tsa.fasta_urls(acc):
            try:
                ncbi.http_download(u, gz, timeout=900)
                info["source_url"] = u
                last = None
                break
            except Exception as exc:
                last = exc
        if last is not None:
            raise RuntimeError("download failed: %s" % last)
        recs, n_total = tsa.sample_transcripts(gz, cap, seed=seed, disjoint_from=disjoint,
                                               min_len=int(spec.get("min_len", 300)))
        os.remove(gz)
        info.update({"accession": acc, "organism": row.get("scientific_name", ""), "tax_id": row.get("tax_id", ""),
                     "first_public": row.get("first_public", ""), "units_total": n_total,
                     "resolved_from": row.get("resolved_from", "")})
    if not recs:
        raise RuntimeError("empty sample")
    pre = spec.get("prefix", "")
    recs = [("%s|%s" % (pre, n) if pre else n, s) for n, s in recs]
    info.update({"units_used": len(recs), "sample_bp": sum(len(s) for _, s in recs)})
    return recs, info


# ------------------------------------------------------------------ Codetta, keeping the alignment

def ensure_resources(codetta_dir, data_dir):
    """run_pilot.setup_resources under the same file lock the other eukcode jobs use."""
    import fcntl  # Linux only; the job runs in a Linux container
    import run_pilot
    os.makedirs(data_dir, exist_ok=True)
    with open(os.path.join(data_dir, ".eukcode-pfam.lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            part = os.path.join(data_dir, "Pfam-A_enone.tar.gz.part")
            waited = 0
            while os.path.exists(part) and waited < 3600:
                s0 = os.path.getsize(part)
                time.sleep(30)
                waited += 30
                if not os.path.exists(part) or os.path.getsize(part) == s0:
                    break
            return run_pilot.setup_resources(codetta_dir, data_dir)
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


def run_codetta_keep(label, inp, codetta_dir, res, out_dir, work, threads, timeout_s):
    align_dir = os.path.join(work, "align")
    os.makedirs(align_dir, exist_ok=True)
    cdir = os.path.join(out_dir, "codetta")
    os.makedirs(cdir, exist_ok=True)
    out_file = os.path.join(cdir, "%s.genetic_code.out" % label)
    summ_file = os.path.join(cdir, "%s.summary.csv" % label)
    for p in (out_file, summ_file):
        if os.path.exists(p):
            os.remove(p)
    cmd = [sys.executable, os.path.join(codetta_dir, "codetta.py"), inp, "--resource_directory", res,
           "--align_output", os.path.join(align_dir, label), "--inference_output", out_file, "-s", summ_file]
    os.makedirs(os.path.join(out_dir, "logs"), exist_ok=True)
    t0 = time.time()
    with open(os.path.join(out_dir, "logs", "%s.codetta.log" % label), "w", encoding="utf-8") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, cwd=codetta_dir,
                              env=dict(os.environ, HMMER_NCPU=str(threads)), timeout=max(60, timeout_s))
    if proc.returncode != 0:
        raise RuntimeError("codetta exited %d" % proc.returncode)
    with open(summ_file, encoding="utf-8") as fh:
        rows = codetta_io.parse_summary(fh.read())
    if not rows:
        raise RuntimeError("empty Codetta summary")
    found = sorted(glob.glob(os.path.join(align_dir, "*alignment_output*")))
    pieces = sorted(glob.glob(os.path.join(align_dir, "*sequence_pieces*")) +
                    glob.glob(os.path.join(os.path.dirname(inp), "*sequence_pieces*")))
    return {"inferred": rows[-1]["inferred_gencode"], "minutes": round((time.time() - t0) / 60, 2),
            "alignment_files": found, "pieces_files": pieces}


# ------------------------------------------------------------------ one case

def lineage_of(taxid):
    if not taxid:
        return ""
    t = ncbi.taxonomy_lookup([int(taxid)])
    return (t.get(int(taxid)) or {}).get("lineage", "")


def run_case(case, args, res):
    label = case["label"]
    work = os.path.join(args.tmp_dir, label)
    t0 = time.time()
    out = {"label": label, "kind": case.get("kind", "record"), "codons_requested": case["codons"]}
    try:
        recs, infos = [], []
        for i, spec in enumerate(case["inputs"]):
            r, info = load_input(spec, os.path.join(work, "in%d" % i))
            info["first_record"] = len(recs)
            recs.extend(r)
            infos.append(info)
            log("%s: input %d %s %s, %d units, %.2f Mb" % (label, i, info.get("accession"), info.get("organism", ""),
                                                          info["units_used"], info["sample_bp"] / 1e6))
        out["inputs"] = infos
        sp_taxid = case.get("species_tax_id") or infos[0].get("tax_id")
        sp_lineage = case.get("species_lineage") or lineage_of(sp_taxid)
        out.update({"species_tax_id": sp_taxid, "species_lineage": sp_lineage,
                    "species_major_clade": list(support.major_clade(sp_lineage) or [])})
        inp = os.path.join(work, label + ".fna")
        fasta.write_fasta(recs, inp)
        c = run_codetta_keep(label, inp, args.codetta_dir, res, args.out_dir, work, args.threads,
                             args.per_case_max_min * 60)
        out.update({"inferred": c["inferred"], "codetta_minutes": c["minutes"]})
        if not c["alignment_files"]:
            raise RuntimeError("Codetta left no alignment output file")
        names = [n for n, _ in recs]
        seqs = [s for _, s in recs]
        pieces = support.codetta_pieces([len(s) for s in seqs])
        pseqs = support.piece_sequences(seqs, pieces)
        pf = "absent"
        if c["pieces_files"]:
            got = support.read_pieces_file(c["pieces_files"][0])
            pf = "agrees" if (len(got) == len(pseqs) and all(got.get(i, "").upper() == s for i, s in enumerate(pseqs))) else "differs"
        codons = [x.upper() for x in case["codons"]]
        t = support.tally(support.iter_alignment(c["alignment_files"][0]), pieces, pseqs, codons)
        out["mapping"] = {"rows": t["rows"], "mismatches": t["mismatches"], "examples": t["mismatch_examples"],
                          "pieces": len(pieces), "pieces_file": pf}
        log("%s: %s, %d alignment rows, %d mismatches, pieces file %s" % (label, c["inferred"], t["rows"],
                                                                           t["mismatches"], pf))
        idx = codes.CODON_INDEX
        blast_cache = {}
        out["codons"] = {}
        for cd in codons:
            d = t["by_codon"][cd]
            top = support.top_records(d["per_record"], args.top_n)
            jobs = []
            for rec in top:
                if rec not in blast_cache:
                    pr = d["per_record"][rec]
                    s, e = support.region(pr["start"], pr["end"], len(seqs[rec]))
                    jobs.append((rec, s, e))
            with ThreadPoolExecutor(max_workers=args.blast_workers) as ex:
                done = list(ex.map(lambda j: (j, ebi_blast.search(seqs[j[0]][j[1]:j[2]], title="%s_%d" % (label, j[0]),
                                                                    exclude_taxids=[sp_taxid] if sp_taxid else [])), jobs))
            for (rec, s, e), r in done:
                blast_cache[rec] = dict(r, region=[s, e])
            taxids = [blast_cache[r]["hit"]["taxid"] for r in top if blast_cache[r]["hit"] and blast_cache[r]["hit"]["taxid"]]
            lin = ncbi.taxonomy_lookup(taxids) if taxids else {}
            classes, rows = {}, []
            for rec in top:
                b = blast_cache[rec]
                h = b["hit"]
                if not b["status"].startswith("ok"):
                    cls, clade, hl = "error", None, ""
                elif h is None:
                    cls, clade, hl = "no_hit", None, ""
                else:
                    hl = (lin.get(h["taxid"]) or {}).get("lineage", "") if h["taxid"] else ""
                    cls = support.classify_hit(hl, sp_lineage)  # "" (no taxid or lookup failed) -> unknown
                    clade = support.major_clade(hl)
                classes[rec] = (cls, clade)
                pr = d["per_record"][rec]
                rows.append({"record": rec, "name": names[rec], "obs": pr["obs"], "profiles": pr["profiles"],
                             "region": b.get("region"), "class": cls, "hit": h, "hit_lineage": hl,
                             "search_status": b["status"]})
            g = support.grade_check2({r: d["per_record"][r]["obs"] for r in top}, classes, d["total"])
            entry = {"call": c["inferred"][idx[cd]], "total_obs": d["total"], "records_with_obs": len(d["per_record"]),
                     "check2": g, "top": rows}
            if len(infos) > 1:
                bounds = [(i["first_record"], i["prefix"] or str(k)) for k, i in enumerate(infos)]

                def origin(rec, bounds=bounds):
                    lab = bounds[0][1]
                    for first, name in bounds:
                        if rec >= first:
                            lab = name
                    return lab
                entry["origin_shares"] = support.origin_shares(d["per_record"], origin)
            out["codons"][cd] = entry
            log("%s %s: call %s, %d obs in %d records, check 2 %s %s" % (label, cd, entry["call"], d["total"],
                                                                         len(d["per_record"]), g["outcome"],
                                                                         entry.get("origin_shares", "")))
    except subprocess.TimeoutExpired:
        out["error"] = "codetta timeout"
    except Exception as exc:
        out["error"] = str(exc)[:300]
        log("%s: error %s" % (label, out["error"]))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    out["wall_minutes"] = round((time.time() - t0) / 60, 2)
    os.makedirs(os.path.join(args.out_dir, "check2"), exist_ok=True)
    with open(os.path.join(args.out_dir, "check2", label + ".json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1, default=str)
    return out


# ------------------------------------------------------------------ grading and report

def grade_demo(results):
    """P8-P10 (A4.5)."""
    by = {r["label"]: r for r in results}
    finished = [r for r in results if "error" not in r]
    rows = sum(r["mapping"]["rows"] for r in finished)
    bad = sum(r["mapping"]["mismatches"] for r in finished)
    p8 = {"rule": "every alignment row of every finished case maps back to the same codon (0 mismatches, >0 rows)",
          "rows": rows, "mismatches": bad, "testable": rows > 0, "holds": rows > 0 and bad == 0}
    mix = by.get("mix_uar")
    p9 = {"rule": "in the ciliate+alga mix, >= 90% of TAA observations come from ciliate transcripts and check 2 on TAA passes",
          "testable": bool(mix and "error" not in mix and "TAA" in mix.get("codons", {}))}
    if p9["testable"]:
        e = mix["codons"]["TAA"]
        p9.update({"ciliate_share": e.get("origin_shares", {}).get("CIL", 0.0), "check2": e["check2"]["outcome"]})
        p9["holds"] = p9["ciliate_share"] >= 0.9 and e["check2"]["outcome"] == "pass"
    else:
        p9["holds"] = False
    eup = by.get("euplotes_uga")
    p10 = {"rule": "check 2 on the Euplotes UGA = Cys signal passes", "testable": bool(eup and "error" not in eup and "TGA" in eup.get("codons", {}))}
    if p10["testable"]:
        p10["check2"] = eup["codons"]["TGA"]["check2"]["outcome"]
        p10["holds"] = p10["check2"] == "pass"
    else:
        p10["holds"] = False
    return {"P8": p8, "P9": p9, "P10": p10}


def report_md(results, grades, meta):
    L = ["# eukcode-scan check 2", "",
         "eukcode-scan %s, finished %s UTC. blastx at EMBL-EBI against UniProtKB (E <= 1e-5), the species' own taxid excluded." % (
             meta["eukcode_version"], meta["finished_utc"]), "",
         "| case | inputs | call | codon | obs | records | check 2 | eukaryotic | same clade | top foreign clade | coverage | origin shares |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in results:
        ins = "; ".join("%s %s" % (i.get("accession", "?"), i.get("organism", "")) for i in r.get("inputs", []))
        if "error" in r:
            L.append("| %s | %s | - | - | - | - | error: %s | | | | | |" % (r["label"], ins, r["error"][:80]))
            continue
        for cd, e in r["codons"].items():
            g = e["check2"]
            L.append("| %s | %s | %s | %s | %d | %d | %s%s | %.2f | %.2f | %s %.2f | %.2f | %s |" % (
                r["label"], ins, e["call"], cd, e["total_obs"], e["records_with_obs"], g["outcome"],
                " (concentrated)" if g["concentrated"] else "", g["eukaryotic_share"], g["same_clade_share"],
                g["top_foreign_clade"] or "-", g["top_foreign_share"], g["coverage"],
                json.dumps(e.get("origin_shares", "")) if e.get("origin_shares") else "-"))
    if grades:
        L += ["", "## Pre-registered predictions (A4.5)", ""]
        for k in ("P8", "P9", "P10"):
            p = grades[k]
            v = "HOLDS" if p["holds"] else ("CANNOT BE TESTED" if not p["testable"] else "DOES NOT HOLD")
            L.append("- **%s** (%s): %s. %s" % (k, p["rule"], v, {a: b for a, b in p.items() if a not in ("rule", "holds", "testable")}))
    return "\n".join(L) + "\n"


# ------------------------------------------------------------------ CLI

def cases_from_records(records_dir, only=""):
    keep = set(x.strip() for x in only.split(",") if x.strip())
    cases = []
    for path in sorted(glob.glob(os.path.join(records_dir, "*.json"))):
        with open(path, encoding="utf-8") as fh:
            r = json.load(fh)
        if "error" in r or not r.get("candidates"):
            continue
        if keep and r.get("accession") not in keep:
            continue
        kind = "tsa" if r.get("data_type") == "transcriptome" else "genome"
        cases.append({"label": "%s_c2" % r.get("label", r["accession"]), "kind": "record",
                      "inputs": [{"kind": kind, "accession": r["accession"], "organism": r.get("organism", ""),
                                  "tax_id": r.get("tax_id", ""), "cap_mb": r.get("cap_mb", 2.0),
                                  "seed": r.get("seed", 0), "disjoint_from": r.get("disjoint_from")}],
                      "codons": sorted(set(c["codon"] for c in r["candidates"])),
                      "species_tax_id": r.get("tax_id", ""), "species_lineage": r.get("lineage", "")})
    return cases


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["demo", "records"])
    ap.add_argument("--cases", default=os.path.join(HERE, "check2_demo.json"))
    ap.add_argument("--records-dir", default="")
    ap.add_argument("--only", default="")
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--codetta-dir", default=os.path.join(HERE, "codetta"))
    ap.add_argument("--tmp-dir", default="/tmp/eukcode-check2")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--top-n", type=int, default=support.TOP_N)
    ap.add_argument("--blast-workers", type=int, default=4)
    ap.add_argument("--per-case-max-min", type=float, default=120)
    return ap


def main(argv=None):
    global LOG_PATH
    args = build_parser().parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    LOG_PATH = os.path.join(args.out_dir, "log.txt")
    import run_pilot
    run_pilot.LOG_PATH = LOG_PATH
    if args.mode == "demo":
        with open(args.cases, encoding="utf-8") as fh:
            cases = json.load(fh)["cases"]
    else:
        cases = cases_from_records(args.records_dir, args.only)
    hmmscan = os.path.join(args.codetta_dir, "hmmer-3.3.2", "bin", "hmmscan")
    if not os.path.exists(hmmscan):
        log("ERROR: Codetta's HMMER build is missing at %s" % hmmscan)
        return 2
    res = ensure_resources(args.codetta_dir, args.data_dir)
    log("eukcode-scan %s check 2, %d cases" % (__version__, len(cases)))
    results = [run_case(c, args, res) for c in cases]
    grades = grade_demo(results) if args.mode == "demo" else None
    meta = {"eukcode_version": __version__, "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "top_n": args.top_n, "database": "uniprotkb", "exp": "1e-5"}
    with open(os.path.join(args.out_dir, "check2_summary.json"), "w", encoding="utf-8") as fh:
        json.dump({"meta": meta, "predictions": grades,
                   "cases": [{k: v for k, v in r.items() if k != "codons"} |
                             {"codons": {cd: {k: v for k, v in e.items() if k != "top"} for cd, e in r.get("codons", {}).items()}}
                             for r in results]}, fh, indent=1, default=str)
    with open(os.path.join(args.out_dir, "CHECK2.md"), "w", encoding="utf-8") as fh:
        fh.write(report_md(results, grades, meta))
    log("check 2 done: %s" % json.dumps(grades, default=str))
    shutil.rmtree(args.tmp_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
