#!/usr/bin/env python3
# eukcode-scan stage 2: run Codetta on one shard of the protist census.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Screen one shard CSV (from make_shards.py) for alternative nuclear genetic codes.

For each assembly: download from NCBI Datasets into a temporary folder (never into the
persistent data volume), draw the seeded chunk sample with low memory, run Codetta 2.0, compare
with the NCBI-assigned code and classify every difference (eukcode/known.py). Each species
writes its own JSON record as soon as it finishes, so a shard cut short by the time budget
still counts. Only the Pfam database lives in --data-dir.

  python run_screen.py --shard shards/shard_00.csv --data-dir /data --out-dir /out
  python run_screen.py --shard shards/shard_00.csv --only GCA_x.1 --seed 1 --disjoint-from 0 --cap-mb 4
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import run_pilot  # noqa: E402
from eukcode import __version__, codes, fasta, known, ncbi, shards  # noqa: E402

log = run_pilot.log


def label_for(accession, seed, disjoint_from):
    base = accession.replace(".", "_")
    if seed == 0 and disjoint_from is None:
        return base
    return "%s_s%d" % (base, seed)


def screen_one(row, args, res, deadline):
    acc = row["accession"].strip()
    label = label_for(acc, args.seed, args.disjoint_from)
    rec_path = os.path.join(args.out_dir, "species", label + ".json")
    if os.path.exists(rec_path):
        with open(rec_path, encoding="utf-8") as fh:
            old = json.load(fh)
        if "error" not in old:
            return old
    gencode = shards.gencode_of(row)
    rec = {"accession": acc, "label": label, "organism": row.get("organism", ""),
           "group": row.get("group", ""), "tax_id": row.get("tax_id", ""), "ncbi_gencode": gencode,
           "lineage": row.get("lineage", ""), "seed": args.seed, "disjoint_from": args.disjoint_from}
    work = os.path.join(args.tmp_dir, label)
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    try:
        zp = ncbi.download_genome_zip(acc, work)
        raw = os.path.join(work, "raw.fna")
        ncbi.extract_fasta(zp, raw)
        os.remove(zp)
        cap_mb = args.cap_mb
        sample, n_chunks = fasta.sample_file(raw, int(cap_mb * 1e6), seed=args.seed,
                                             disjoint_from=args.disjoint_from)
        gc = fasta.gc_fraction(sample)
        if (gc < 0.25 or gc > 0.70) and args.cap_mb_extreme_gc > cap_mb:
            cap_mb = args.cap_mb_extreme_gc
            sample, n_chunks = fasta.sample_file(raw, int(cap_mb * 1e6), seed=args.seed,
                                                 disjoint_from=args.disjoint_from)
            gc = fasta.gc_fraction(sample)
        os.remove(raw)
        if not sample:
            raise RuntimeError("no usable chunks after filtering")
        inp = os.path.join(work, label + ".fna")
        fasta.write_fasta(sample, inp)
        rec.update({"cap_mb": cap_mb, "chunks_total": n_chunks, "chunks_used": len(sample),
                    "sample_bp": sum(len(s) for _, s in sample), "sample_gc": round(gc, 4),
                    "prep_minutes": round((time.time() - t0) / 60, 2)})
        del sample
        remaining = deadline - time.time()
        if remaining < 300:
            raise RuntimeError("time budget reached before Codetta")
        timeout = min(remaining, args.per_species_max_min * 60)
        r = run_pilot.run_codetta({"label": label}, {"input_path": inp}, args.codetta_dir, res,
                                  args.out_dir, work, args.threads, timeout)
        cands = known.candidates(r["inferred"], gencode, rec["lineage"], rec["organism"])
        keep = set(c["codon"] for c in cands) | {"TAA", "TAG", "TGA"}
        rec.update({"inferred": r["inferred"], "best_models": r["best_models"],
                    "codetta_minutes": r["minutes"], "n_called": codes.called_count(r["inferred"]),
                    "stop_calls": known.stop_calls(r["inferred"]), "candidates": cands,
                    "top_class": known.most_interesting(cands),
                    "codon_detail": {k: v for k, v in r["codon_table"].items() if k in keep}})
        log("%s %s: %s cands=%s (%.1f min)" % (acc, rec["organism"], r["inferred"],
                                               ",".join("%s=%s/%s" % (c["codon"], c["called"], c["class"]) for c in cands) or "none",
                                               r["minutes"]))
    except subprocess.TimeoutExpired:
        rec["error"] = "codetta timeout"
    except Exception as exc:
        rec["error"] = str(exc)[:300]
        log("%s: error %s" % (acc, rec["error"]))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    rec["wall_minutes"] = round((time.time() - t0) / 60, 2)
    os.makedirs(os.path.dirname(rec_path), exist_ok=True)
    with open(rec_path, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=1)
    return rec


ROW_FIELDS = ["accession", "organism", "group", "tax_id", "ncbi_gencode", "seed", "cap_mb", "sample_bp",
              "sample_gc", "codetta_minutes", "n_called", "inferred", "stops", "n_candidates", "top_class",
              "candidates", "error"]
CAND_FIELDS = ["accession", "organism", "group", "ncbi_gencode", "codon", "called", "expected", "kind",
               "class", "entries", "table_fallback", "n_used", "lineage"]


def write_tables(out_dir, recs):
    with open(os.path.join(out_dir, "screen_rows.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in recs:
            row = dict(r)
            row["stops"] = "".join(r.get("stop_calls", {}).get(s, "") for s in ("TAA", "TAG", "TGA"))
            row["n_candidates"] = len(r.get("candidates", []))
            row["candidates"] = " ".join("%s=%s(%s)" % (c["codon"], c["called"], c["class"])
                                         for c in r.get("candidates", []))
            w.writerow(row)
    with open(os.path.join(out_dir, "candidates.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=CAND_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in recs:
            for c in r.get("candidates", []):
                row = dict(c)
                row.update({k: r.get(k, "") for k in ("accession", "organism", "group", "ncbi_gencode", "lineage")})
                row["entries"] = ";".join(c.get("entries", []))
                row["n_used"] = r.get("codon_detail", {}).get(c["codon"], {}).get("n_used", "")
                w.writerow(row)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shard", required=True)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--codetta-dir", default=os.path.join(HERE, "codetta"))
    ap.add_argument("--tmp-dir", default="/tmp/eukcode-screen")
    ap.add_argument("--cap-mb", type=float, default=2.0)
    ap.add_argument("--cap-mb-extreme-gc", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--disjoint-from", type=int, default=None)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--per-species-max-min", type=float, default=90)
    ap.add_argument("--time-budget-min", type=float, default=280)
    ap.add_argument("--only", default="", help="comma-separated accessions")
    args = ap.parse_args(argv)

    t_start = time.time()
    deadline = t_start + args.time_budget_min * 60
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.tmp_dir, exist_ok=True)
    run_pilot.LOG_PATH = os.path.join(args.out_dir, "log.txt")
    rows = shards.read_species(args.shard)
    if args.only:
        keep = set(x.strip() for x in args.only.split(",") if x.strip())
        rows = [r for r in rows if r.get("accession", "").strip() in keep]
    hmmscan = os.path.join(args.codetta_dir, "hmmer-3.3.2", "bin", "hmmscan")
    if not os.path.exists(hmmscan):
        log("ERROR: Codetta's HMMER build is missing at %s" % hmmscan)
        return 2
    commit = subprocess.run(["git", "-C", args.codetta_dir, "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip() or "unknown"
    log("eukcode-scan %s screen, Codetta %s, %d assemblies in %s" % (__version__, commit, len(rows), args.shard))
    res = run_pilot.setup_resources(args.codetta_dir, args.data_dir)

    recs, not_run = {}, []

    def _do(row):
        if deadline - time.time() < 600:
            return row, None
        return row, screen_one(row, args, res, deadline)

    if rows:
        row, rec = _do(rows[0])  # alone first: builds Codetta's emission cache if it is missing
        if rec is None:
            not_run.append(row["accession"])
        else:
            recs[row["accession"]] = rec
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
            futs = [ex.submit(_do, r) for r in rows[1:]]
            for f in as_completed(futs):
                row, rec = f.result()
                if rec is None:
                    not_run.append(row["accession"])
                else:
                    recs[row["accession"]] = rec
    ordered = [recs[r["accession"]] for r in rows if r["accession"] in recs]
    write_tables(args.out_dir, ordered)
    classes = {}
    for r in ordered:
        for c in r.get("candidates", []):
            classes[c["class"]] = classes.get(c["class"], 0) + 1
    summary = {"eukcode_version": __version__, "codetta_commit": commit, "shard": os.path.basename(args.shard),
               "seed": args.seed, "disjoint_from": args.disjoint_from, "cap_mb": args.cap_mb,
               "cap_mb_extreme_gc": args.cap_mb_extreme_gc, "assemblies": len(rows),
               "finished": sum(1 for r in ordered if "error" not in r),
               "errors": sum(1 for r in ordered if "error" in r), "not_run": sorted(not_run),
               "species_with_candidates": sum(1 for r in ordered if r.get("candidates")),
               "candidate_classes": classes,
               "wall_minutes": round((time.time() - t_start) / 60, 1),
               "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    with open(os.path.join(args.out_dir, "screen_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    log("screen done: %s" % json.dumps(summary))
    shutil.rmtree(args.tmp_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
