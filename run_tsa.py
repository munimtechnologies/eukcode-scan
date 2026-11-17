#!/usr/bin/env python3
# eukcode-scan stage 3: the transcriptome (INSDC TSA) tranche.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Census and screen of protist transcriptome assemblies (PREREGISTRATION.md, addendum A2).

  python run_tsa.py census --out-dir /out
      Lists every INSDC TSA set under the screen's protist groups (ENA portal API), keeps one
      set per species, looks up NCBI's nuclear genetic code, checks the download URL, resolves
      the transcriptome controls (tsa_controls.json), and writes:
      tsa_sets.csv, tsa_species.csv, tsa_skipped.csv, tsa_controls_resolved.csv,
      tsa_census_summary.json. No Codetta is run.

  python run_tsa.py screen --census /work/tsa_species.csv --shard-index 0 --per-shard 12 \\
      --data-dir /data --out-dir /out --cap-mb 2
      Runs Codetta on one deterministic shard (eukcode.tsa.plan_shards). Needs Codetta
      (install_codetta.sh) and the Pfam cache; output matches run_screen.py.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import threading
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from eukcode import __version__, ncbi, tsa  # noqa: E402

DEFAULT_GROUPS = ["Sar", "Discoba", "Metamonada", "Amoebozoa", "Haptista", "Cryptophyceae", "Rhodophyta",
                  "Glaucocystophyceae", "Chlorophyta", "Apusozoa", "Ichthyosporea", "Choanoflagellata",
                  "Rotosphaerida", "Ancyromonadida", "CRuMs", "Hemimastigophora", "Provora", "Picozoa",
                  "Telonemia", "Rhodelphea", "Filasterea", "Aphelida", "Breviatea"]

_lock = threading.Lock()
LOG_PATH = None


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    with _lock:
        print(line, flush=True)
        if LOG_PATH:
            with open(LOG_PATH, "a", encoding="utf-8") as fh:
                fh.write(line + "\n")


def resolve_taxid(name):
    """NCBI taxid for a scientific name (exact field first, then any name field)."""
    for url in (tsa.taxid_search_url(name),
                "%s/esearch.fcgi?db=taxonomy&retmode=json&term=%s" % (tsa.EUTILS, urllib.parse.quote(name))):
        try:
            data = json.loads(ncbi.http_get(url, pause=0.4).decode("utf-8"))
        except Exception:
            continue
        ids = data.get("esearchresult", {}).get("idlist", [])
        if len(ids) == 1:
            return int(ids[0])
    return None


def ena_sets(taxid):
    return tsa.parse_tsv(ncbi.http_get(tsa.ena_search_url(taxid), timeout=300, pause=0.5).decode("utf-8"))


def url_ok(url, timeout=60):
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": ncbi.UA})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            n = resp.headers.get("Content-Length")
            return resp.status == 200, int(n) if n and n.isdigit() else None
    except Exception:
        return False, None


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


SPECIES_FIELDS = ["group", "accession", "organism", "tax_id", "length", "first_public", "tier", "n_sets",
                  "study_accession", "ncbi_gencode", "lineage", "url", "url_ok", "gz_bytes"]
SET_FIELDS = ["group"] + tsa.ENA_FIELDS + ["tier", "environmental", "tpa", "mmetsp"]


def census(args):
    groups = [g.strip() for g in args.groups.split(",") if g.strip()] if args.groups else DEFAULT_GROUPS
    per_group, all_sets, chosen_rows, skipped_all = {}, [], [], []
    seen_acc, seen_tax = set(), set()
    for g in groups:
        tid = resolve_taxid(g)
        if tid is None:
            per_group[g] = {"error": "taxid not resolved"}
            log("group %s: taxid not resolved" % g)
            continue
        try:
            rows = ena_sets(tid)
        except Exception as exc:
            per_group[g] = {"taxid": tid, "error": str(exc)[:200]}
            log("group %s: ENA error %s" % (g, exc))
            continue
        rows = [r for r in rows if r.get("accession") not in seen_acc]
        for r in rows:
            seen_acc.add(r.get("accession"))
            f = tsa.flags(r)
            s = dict(r)
            s.update({"group": g, "tier": tsa.set_tier(r), "environmental": f["environmental"],
                      "tpa": f["tpa"], "mmetsp": f["mmetsp"]})
            all_sets.append(s)
        chosen, skipped = tsa.choose_sets(rows, min_bp=args.min_bp)
        chosen = [c for c in chosen if c["tax_id"] not in seen_tax]
        for c in chosen:
            seen_tax.add(c["tax_id"])
        for s in skipped:
            s["group"] = g
        skipped_all.extend(skipped)
        chosen_rows.extend((g, c) for c in chosen)
        per_group[g] = {"taxid": tid, "sets": len(rows), "species": len(chosen), "skipped": len(skipped)}
        log("group %s (taxid %d): %d sets, %d species" % (g, tid, len(rows), len(chosen)))

    tax = ncbi.taxonomy_lookup([int(c["tax_id"]) for _g, c in chosen_rows])
    species = [tsa.species_row(c, tax, g) for g, c in chosen_rows]

    if not args.no_url_check:
        log("checking %d download URLs" % len(species))
        for s in species:
            ok, size, used = False, None, ""
            for u in tsa.fasta_urls(s["accession"]):
                ok, size = url_ok(u)
                if ok:
                    used = u
                    break
            s["url"], s["url_ok"], s["gz_bytes"] = used, ok, size if size is not None else ""

    write_csv(os.path.join(args.out_dir, "tsa_sets.csv"), all_sets, SET_FIELDS)
    write_csv(os.path.join(args.out_dir, "tsa_species.csv"), species, SPECIES_FIELDS)
    write_csv(os.path.join(args.out_dir, "tsa_skipped.csv"), skipped_all, ["group"] + tsa.ENA_FIELDS + ["skip_reason"])

    # controls: resolved by name, independent of the group census
    with open(args.controls, encoding="utf-8") as fh:
        cfg = json.load(fh)
    ctrl_rows = []
    for role in ("positive", "standard"):
        for c in cfg.get(role, []):
            row = {"label": c["label"], "role": role, "table": c["table"],
                   "diagnostic": " ".join("%s=%s" % kv for kv in sorted(c.get("diagnostic", {}).items())),
                   "accession": "", "organism": "", "length": "", "first_public": "", "tier": "",
                   "ncbi_gencode": "", "status": "unresolved"}
            for name in c["taxa"]:
                tid = resolve_taxid(name)
                if tid is None:
                    continue
                try:
                    picked, _sk = tsa.choose_sets(ena_sets(tid), min_bp=args.min_bp)
                except Exception as exc:
                    row["status"] = "ENA error: %s" % str(exc)[:100]
                    continue
                if not picked:
                    row["status"] = "no usable TSA set"
                    continue
                best = max(picked, key=lambda r: tsa._int(r.get("base_count")))
                t = ncbi.taxonomy_lookup([int(best["tax_id"])])
                srow = tsa.species_row(best, t, role)
                row.update({k: srow[k] for k in ("accession", "organism", "length", "first_public", "tier",
                                                   "ncbi_gencode")})
                row["status"] = "resolved"
                break
            ctrl_rows.append(row)
            log("control %s: %s %s" % (c["label"], row["status"], row["accession"]))
    write_csv(os.path.join(args.out_dir, "tsa_controls_resolved.csv"), ctrl_rows,
              ["label", "role", "table", "diagnostic", "accession", "organism", "length", "first_public",
               "tier", "ncbi_gencode", "status"])

    def count(key):
        out = {}
        for s in species:
            out[str(s.get(key))] = out.get(str(s.get(key)), 0) + 1
        return out

    by_group_tier = {}
    for s in species:
        by_group_tier.setdefault(s["group"], {}).setdefault(s["tier"], 0)
        by_group_tier[s["group"]][s["tier"]] += 1
    shards = tsa.plan_shards(species, per_shard=args.per_shard)
    summary = {
        "eukcode_version": __version__,
        "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "groups": per_group, "sets": len(all_sets), "species": len(species),
        "species_by_tier": count("tier"), "species_by_group_and_tier": by_group_tier,
        "ncbi_gencode_counts": count("ncbi_gencode"),
        "url_ok": sum(1 for s in species if s.get("url_ok") is True),
        "url_missing": [s["accession"] for s in species if s.get("url_ok") is False][:200],
        "total_bp_species": sum(tsa._int(s["length"]) for s in species),
        "sample_mb_at_cap_2": round(sum(min(2e6, tsa._int(s["length"])) for s in species) / 1e6, 1),
        "per_shard": args.per_shard, "shards": len(shards),
        "shards_tier_new": sum(1 for sh in shards if sh and sh[0].get("tier") == "new"),
        "controls": {r["label"]: r["status"] for r in ctrl_rows},
        "swart_cutoff": tsa.SWART_CUTOFF,
    }
    with open(os.path.join(args.out_dir, "tsa_census_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    log("census done: %d sets, %d species %s" % (len(all_sets), len(species), json.dumps(summary["species_by_tier"])))
    return 0


def screen(args):
    import run_pilot  # needs Codetta; only the screen mode imports these
    import run_screen
    from eukcode import codes, known, shards as shardmod

    run_pilot.LOG_PATH = LOG_PATH
    rows = shardmod.read_species(args.census)
    if args.only:
        keep = set(x.strip() for x in args.only.split(",") if x.strip())
        rows = [r for r in rows if r.get("accession", "").strip() in keep]
    else:
        plan = tsa.plan_shards(rows, per_shard=args.per_shard)
        if args.shard_index >= len(plan):
            log("shard index %d out of range (%d shards)" % (args.shard_index, len(plan)))
            return 2
        rows = plan[args.shard_index]
    hmmscan = os.path.join(args.codetta_dir, "hmmer-3.3.2", "bin", "hmmscan")
    if not os.path.exists(hmmscan):
        log("ERROR: Codetta's HMMER build is missing at %s" % hmmscan)
        return 2
    commit = subprocess.run(["git", "-C", args.codetta_dir, "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip() or "unknown"
    res = run_pilot.setup_resources(args.codetta_dir, args.data_dir)
    t_start = time.time()
    deadline = t_start + args.time_budget_min * 60
    os.makedirs(args.tmp_dir, exist_ok=True)

    def one(row):
        acc = row["accession"].strip()
        label = run_screen.label_for(acc, args.seed, args.disjoint_from)
        rec_path = os.path.join(args.out_dir, "species", label + ".json")
        if os.path.exists(rec_path):
            with open(rec_path, encoding="utf-8") as fh:
                old = json.load(fh)
            if "error" not in old:
                return old
        gencode = shardmod.gencode_of(row)
        rec = {"accession": acc, "label": label, "organism": row.get("organism", ""), "group": row.get("group", ""),
               "tax_id": row.get("tax_id", ""), "ncbi_gencode": gencode, "lineage": row.get("lineage", ""),
               "tier": row.get("tier", ""), "data_type": "transcriptome", "seed": args.seed,
               "disjoint_from": args.disjoint_from}
        work = os.path.join(args.tmp_dir, label)
        os.makedirs(work, exist_ok=True)
        t0 = time.time()
        try:
            if deadline - time.time() < 600:
                raise RuntimeError("time budget reached before download")
            gz = os.path.join(work, "set.fsa_nt.gz")
            last = None
            for u in ([row["url"]] if row.get("url") else []) + tsa.fasta_urls(acc):
                try:
                    ncbi.http_download(u, gz, timeout=900)
                    rec["source_url"] = u
                    last = None
                    break
                except Exception as exc:
                    last = exc
            if last is not None:
                raise RuntimeError("download failed: %s" % last)
            sample, n_usable = tsa.sample_transcripts(gz, int(args.cap_mb * 1e6), seed=args.seed,
                                                      disjoint_from=args.disjoint_from, min_len=args.min_len)
            os.remove(gz)
            if not sample:
                raise RuntimeError("no usable transcripts")
            inp = os.path.join(work, label + ".fna")
            from eukcode import fasta
            fasta.write_fasta(sample, inp)
            rec.update({"cap_mb": args.cap_mb, "transcripts_usable": n_usable, "transcripts_used": len(sample),
                        "sample_bp": sum(len(s) for _, s in sample), "sample_gc": round(fasta.gc_fraction(sample), 4),
                        "prep_minutes": round((time.time() - t0) / 60, 2)})
            del sample
            remaining = deadline - time.time()
            if remaining < 300:
                raise RuntimeError("time budget reached before Codetta")
            r = run_pilot.run_codetta({"label": label}, {"input_path": inp}, args.codetta_dir, res, args.out_dir,
                                      work, args.threads, min(remaining, args.per_species_max_min * 60))
            cands = known.candidates(r["inferred"], gencode, rec["lineage"], rec["organism"])
            keep = set(c["codon"] for c in cands) | {"TAA", "TAG", "TGA"}
            rec.update({"inferred": r["inferred"], "best_models": r["best_models"], "codetta_minutes": r["minutes"],
                        "n_called": codes.called_count(r["inferred"]), "stop_calls": known.stop_calls(r["inferred"]),
                        "candidates": cands, "top_class": known.most_interesting(cands),
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

    log("eukcode-scan %s TSA screen, Codetta %s, %d sets" % (__version__, commit, len(rows)))
    recs = []
    if rows:
        recs.append(one(rows[0]))  # alone first: builds Codetta's emission cache if missing
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
            recs.extend(ex.map(one, rows[1:]))
    run_screen.write_tables(args.out_dir, recs)
    classes = {}
    for r in recs:
        for c in r.get("candidates", []):
            classes[c["class"]] = classes.get(c["class"], 0) + 1
    summary = {"eukcode_version": __version__, "codetta_commit": commit, "data_type": "transcriptome",
               "shard_index": None if args.only else args.shard_index, "per_shard": args.per_shard,
               "seed": args.seed, "disjoint_from": args.disjoint_from, "cap_mb": args.cap_mb, "sets": len(rows),
               "finished": sum(1 for r in recs if "error" not in r), "errors": sum(1 for r in recs if "error" in r),
               "species_with_candidates": sum(1 for r in recs if r.get("candidates")), "candidate_classes": classes,
               "wall_minutes": round((time.time() - t_start) / 60, 1),
               "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    with open(os.path.join(args.out_dir, "screen_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    log("TSA screen done: %s" % json.dumps(summary))
    shutil.rmtree(args.tmp_dir, ignore_errors=True)
    return 0


def main(argv=None):
    global LOG_PATH
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("census")
    c.add_argument("--out-dir", default="/out")
    c.add_argument("--groups", default="")
    c.add_argument("--controls", default=os.path.join(HERE, "tsa_controls.json"))
    c.add_argument("--min-bp", type=int, default=1_000_000)
    c.add_argument("--per-shard", type=int, default=12)
    c.add_argument("--no-url-check", action="store_true")
    s = sub.add_parser("screen")
    s.add_argument("--census", required=True)
    s.add_argument("--shard-index", type=int, default=0)
    s.add_argument("--per-shard", type=int, default=12)
    s.add_argument("--only", default="")
    s.add_argument("--data-dir", default="/data")
    s.add_argument("--out-dir", default="/out")
    s.add_argument("--codetta-dir", default=os.path.join(HERE, "codetta"))
    s.add_argument("--tmp-dir", default="/tmp/eukcode-tsa")
    s.add_argument("--cap-mb", type=float, default=2.0)
    s.add_argument("--min-len", type=int, default=300)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--disjoint-from", type=int, default=None)
    s.add_argument("--workers", type=int, default=2)
    s.add_argument("--threads", type=int, default=2)
    s.add_argument("--per-species-max-min", type=float, default=90)
    s.add_argument("--time-budget-min", type=float, default=280)
    args = ap.parse_args(argv)
    os.makedirs(args.out_dir, exist_ok=True)
    LOG_PATH = os.path.join(args.out_dir, "log.txt")
    return census(args) if args.mode == "census" else screen(args)


if __name__ == "__main__":
    sys.exit(main())
