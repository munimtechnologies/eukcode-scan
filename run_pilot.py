#!/usr/bin/env python3
# eukcode-scan pilot: validate Codetta on nuclear genomes with known codes, and census the targets.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Pilot job for the eukaryote-wide genetic code screen (see PREREGISTRATION.md).

1. Uses Codetta 2.0 (Shulgina & Eddy; GPL-3.0, cloned at setup, not bundled) with its
   Pfam 35 'enone' profile database, cached under --data-dir.
2. Downloads each control genome from NCBI Datasets, draws a seeded chunk sample of
   cap_mb megabases, runs Codetta, and compares the inferred code with the NCBI table.
3. Grades the pre-registered predictions P1-P4.
4. In parallel, lists every current NCBI assembly in the protist groups that the full
   screen will cover, with each species' NCBI-assigned nuclear genetic code.

v0.6.0: genomes are sampled with fasta.sample_file (two streaming passes, the same sample as
the in-memory path, A1.2), and genome zips go to --tmp-dir and are deleted after extraction,
so the job fits in about 1.5 GB and caches nothing but the Pfam database.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tarfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from eukcode import __version__, codes, codetta_io, fasta, grade, ncbi  # noqa: E402

PFAM_URL = "http://eddylab.org/publications/Shulgina21/Pfam-A_enone.tar.gz"
PFAM_NAME = "Pfam-A_enone.hmm"

_log_lock = threading.Lock()
LOG_PATH = None


def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    with _log_lock:
        print(line, flush=True)
        if LOG_PATH:
            with open(LOG_PATH, "a", encoding="utf-8") as fh:
                fh.write(line + "\n")


def setup_resources(codetta_dir, data_dir):
    """Make a persistent resource directory: Codetta's small lists + the Pfam database."""
    res = os.path.join(data_dir, "codetta-res")
    os.makedirs(res, exist_ok=True)
    src = os.path.join(codetta_dir, "resources")
    for name in os.listdir(src):
        p = os.path.join(src, name)
        if os.path.isfile(p) and not name.startswith(PFAM_NAME):
            shutil.copy2(p, os.path.join(res, name))
    if not os.path.exists(os.path.join(res, PFAM_NAME)):
        tgz = os.path.join(data_dir, "Pfam-A_enone.tar.gz")
        if not (os.path.exists(tgz) and tarfile.is_tarfile(tgz)):
            log("downloading %s" % PFAM_URL)
            ncbi.http_download(PFAM_URL, tgz, timeout=1800)
        log("extracting Pfam (%.0f MB)" % (os.path.getsize(tgz) / 1e6))
        ext = os.path.join(data_dir, "pfam-extract")
        os.makedirs(ext, exist_ok=True)
        with tarfile.open(tgz) as tf:
            try:
                tf.extractall(ext, filter="data")
            except TypeError:  # Python without extraction filters
                tf.extractall(ext)
        moved = 0
        for root, _dirs, files in os.walk(ext):
            for f in files:
                if f.startswith(PFAM_NAME):
                    shutil.move(os.path.join(root, f), os.path.join(res, f))
                    moved += 1
        shutil.rmtree(ext, ignore_errors=True)
        if not os.path.exists(os.path.join(res, PFAM_NAME)):
            raise RuntimeError("Pfam-A_enone.hmm not found in the tarball")
        os.remove(tgz)
        log("Pfam ready (%d files)" % moved)
    if not os.path.exists(os.path.join(res, PFAM_NAME + ".h3m")):
        log("pressing Pfam with hmmpress")
        subprocess.run([os.path.join(codetta_dir, "hmmer-3.3.2", "bin", "hmmpress"), "-f",
                        os.path.join(res, PFAM_NAME)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return res


def ensure_resources_locked(codetta_dir, data_dir):
    """setup_resources under the file lock the other eukcode jobs use on the shared volume."""
    try:
        import fcntl  # Linux only; the job runs in a Linux container
    except ImportError:
        return setup_resources(codetta_dir, data_dir)
    os.makedirs(data_dir, exist_ok=True)
    with open(os.path.join(data_dir, ".eukcode-pfam.lock"), "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            return setup_resources(codetta_dir, data_dir)
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


def prepare_input(ctrl, data_dir, tmp_dir):
    """Resolve, download, sample. Returns a dict describing the input or raises."""
    summary = None
    if ctrl.get("accession"):
        try:
            summary = ncbi.accession_report(ctrl["accession"])
            if summary:
                summary["resolved_from"] = "pinned accession"
        except Exception as exc:
            log("%s: pinned accession lookup failed (%s); resolving taxon" % (ctrl["label"], exc))
            summary = None
    if summary is None:
        summary = ncbi.resolve_taxon(ctrl["taxa"])
    if summary is None:
        raise RuntimeError("no NCBI assembly found for %s" % ctrl["taxa"])
    acc = summary["accession"]
    gdir = os.path.join(tmp_dir, "genomes", ctrl["label"])
    zp = ncbi.download_genome_zip(acc, gdir)
    raw = os.path.join(tmp_dir, "%s.raw.fna" % ctrl["label"])
    ncbi.extract_fasta(zp, raw)
    shutil.rmtree(gdir, ignore_errors=True)
    sample, n_chunks = fasta.sample_file(raw, int(ctrl["cap_mb"] * 1_000_000), seed=0)
    os.remove(raw)
    path = os.path.join(tmp_dir, "%s.fna" % ctrl["label"])
    fasta.write_fasta(sample, path)
    bp = sum(len(s) for _, s in sample)
    info = {
        "assembly": summary,
        "chunks_total": n_chunks,
        "chunks_used": len(sample),
        "sample_bp": bp,
        "sample_gc": round(fasta.gc_fraction(sample), 4),
        "input_path": path,
    }
    log("%s: %s %s, sample %.2f Mb of %.1f Mb, GC %.3f" % (
        ctrl["label"], acc, summary.get("organism"), bp / 1e6, summary.get("length", 0) / 1e6, info["sample_gc"]))
    return info


def run_codetta(ctrl, info, codetta_dir, res, out_dir, tmp_dir, threads, timeout_s):
    label = ctrl["label"]
    align_dir = os.path.join(tmp_dir, "align", label)
    os.makedirs(align_dir, exist_ok=True)
    cdir = os.path.join(out_dir, "codetta")
    os.makedirs(cdir, exist_ok=True)
    out_file = os.path.join(cdir, "%s.genetic_code.out" % label)
    summ_file = os.path.join(cdir, "%s.summary.csv" % label)
    for p in (out_file, summ_file):
        if os.path.exists(p):
            os.remove(p)
    cmd = [sys.executable, os.path.join(codetta_dir, "codetta.py"), info["input_path"],
           "--resource_directory", res,
           "--align_output", os.path.join(align_dir, label),
           "--inference_output", out_file,
           "-s", summ_file]
    env = dict(os.environ, HMMER_NCPU=str(threads))
    logs = os.path.join(out_dir, "logs")
    os.makedirs(logs, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(logs, "%s.codetta.log" % label), "w", encoding="utf-8") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, env=env,
                              cwd=codetta_dir, timeout=max(60, timeout_s))
    minutes = (time.time() - t0) / 60.0
    shutil.rmtree(align_dir, ignore_errors=True)
    if proc.returncode != 0:
        raise RuntimeError("codetta exited %d after %.1f min" % (proc.returncode, minutes))
    with open(summ_file, encoding="utf-8") as fh:
        rows = codetta_io.parse_summary(fh.read())
    if not rows:
        raise RuntimeError("empty Codetta summary")
    inferred = rows[-1]["inferred_gencode"]
    best_models = rows[-1].get("inferred_gencode_best_models", "")
    detail, detail_error = {}, ""
    if os.path.exists(out_file):
        try:
            with open(out_file, encoding="utf-8") as fh:
                detail = codetta_io.parse_inference_out(fh.read())
        except Exception as exc:  # keep the inferred code even if the detail table changes format
            detail_error = str(exc)[:200]
    return {"inferred": inferred, "best_models": best_models, "minutes": round(minutes, 2),
            "codon_table": detail.get("codons", {}), "detail_error": detail_error}


def census(groups, out_dir):
    """List current assemblies for each group; one per taxon with NCBI's nuclear code."""
    rows, per_group = [], {}
    for g in groups:
        try:
            found = list(ncbi.taxon_reports(g, page_size=1000, max_pages=30))
        except Exception as exc:
            per_group[g] = {"error": str(exc)[:200]}
            log("census %s: error %s" % (g, exc))
            continue
        for s in found:
            s["group"] = g
        rows.extend(found)
        per_group[g] = {"assemblies": len(found)}
        log("census %s: %d assemblies" % (g, len(found)))
    seen, uniq = set(), []
    for r in rows:  # a group list can overlap another (e.g. CRuMs inside nothing else, but be safe)
        if r["accession"] in seen:
            continue
        seen.add(r["accession"])
        uniq.append(r)
    fields = ["group", "accession", "organism", "tax_id", "level", "refseq_category", "release_date",
              "assembly_name", "length", "contig_n50", "source", "annotated"]
    with open(os.path.join(out_dir, "census_assemblies.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in uniq:
            w.writerow(r)
    species = ncbi.one_per_taxon(uniq)
    tax = ncbi.taxonomy_lookup([s["tax_id"] for s in species])
    gencode_counts = {}
    for s in species:
        t = tax.get(s["tax_id"], {})
        s["ncbi_gencode"] = t.get("gencode")
        s["lineage"] = t.get("lineage", "")
        gencode_counts[str(s["ncbi_gencode"])] = gencode_counts.get(str(s["ncbi_gencode"]), 0) + 1
    with open(os.path.join(out_dir, "census_species.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields + ["ncbi_gencode", "lineage"], extrasaction="ignore")
        w.writeheader()
        for s in species:
            w.writerow(s)
    for g in per_group:
        if "error" not in per_group[g]:
            per_group[g]["taxa"] = sum(1 for s in species if s.get("group") == g)
    plan = {}
    for cap in (2, 4, 8):
        plan["sample_mb_at_cap_%d" % cap] = round(sum(min(cap * 1e6, s["length"]) for s in species) / 1e6, 1)
    summary = {"groups": per_group, "assemblies": len(uniq), "taxa": len(species),
               "ncbi_gencode_counts": gencode_counts, "plan": plan}
    with open(os.path.join(out_dir, "census_summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=1)
    log("census done: %d assemblies, %d taxa" % (len(uniq), len(species)))
    return summary


def write_results_md(path, controls, results, graded, meta):
    lines = ["# eukcode-scan pilot results", "",
             "Codetta commit %s, eukcode-scan %s, run %s UTC." % (meta["codetta_commit"], __version__, meta["finished_utc"]),
             "", "| label | group | assembly | organism | sample Mb | min | called | diagnostic | disagreements |",
             "|---|---|---|---|---|---|---|---|---|"]
    for c in controls:
        r = results.get(c["label"], {})
        if "error" in r:
            lines.append("| %s | %s | - | - | - | - | - | error: %s | - |" % (c["label"], c["group"], r["error"][:80]))
            continue
        if "verdict" not in r:
            lines.append("| %s | %s | - | - | - | - | - | not run | - |" % (c["label"], c["group"]))
            continue
        v = r["verdict"]
        diag = ", ".join("%s=%s(%s)" % (k, d["called"], "ok" if d["ok"] else "expected " + d["expected"])
                         for k, d in v["diagnostic"].items()) or "-"
        dis = ", ".join("%s %s->%s" % (d["codon"], d["expected"], d["called"]) for d in v["disagreements"]) or "none"
        a = r["input"]["assembly"]
        lines.append("| %s | %s | %s | %s | %.2f | %.1f | %d | %s | %s |" % (
            c["label"], c["group"], a["accession"], a.get("organism", ""), r["input"]["sample_bp"] / 1e6,
            r["codetta"]["minutes"], v["n_called"], diag, dis))
    lines += ["", "## Pre-registered predictions", ""]
    for k in ("P1", "P2", "P3", "P4"):
        g = graded[k]
        lines.append("- **%s** (%s): %s. %s" % (k, g["rule"], "HOLDS" if g["holds"] else "DOES NOT HOLD",
                                                 json.dumps({x: y for x, y in g.items() if x not in ("rule", "holds")})))
    lines += ["", "## Inferred codes", "", "```"]
    for c in controls:
        r = results.get(c["label"], {})
        if "codetta" in r:
            lines.append("%-9s %s" % (c["label"], r["codetta"]["inferred"]))
    lines.append("%-9s %s" % ("standard", codes.STANDARD))
    lines.append("```")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


def main(argv=None):
    global LOG_PATH
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--controls", default=os.path.join(HERE, "controls.json"))
    ap.add_argument("--codetta-dir", default=os.path.join(HERE, "codetta"))
    ap.add_argument("--tmp-dir", default="/tmp/eukcode")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--time-budget-min", type=float, default=230)
    ap.add_argument("--skip-census", action="store_true")
    ap.add_argument("--only", default="", help="comma-separated labels")
    args = ap.parse_args(argv)

    t_start = time.time()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.tmp_dir, exist_ok=True)
    LOG_PATH = os.path.join(args.out_dir, "log.txt")
    with open(args.controls, encoding="utf-8") as fh:
        cfg = json.load(fh)
    controls = cfg["controls"]
    if args.only:
        keep = set(args.only.split(","))
        controls = [c for c in controls if c["label"] in keep]

    commit = subprocess.run(["git", "-C", args.codetta_dir, "rev-parse", "HEAD"],
                            capture_output=True, text=True).stdout.strip() or "unknown"
    hmmscan = os.path.join(args.codetta_dir, "hmmer-3.3.2", "bin", "hmmscan")
    if not os.path.exists(hmmscan):
        log("ERROR: Codetta's HMMER build is missing at %s" % hmmscan)
        return 2
    log("eukcode-scan %s, Codetta %s" % (__version__, commit))

    census_result = {}
    census_thread = None
    if not args.skip_census:
        def _census():
            try:
                census_result.update(census(cfg.get("census_groups", []), args.out_dir))
            except Exception as exc:
                census_result["error"] = str(exc)[:300]
                log("census failed: %s" % exc)
        census_thread = threading.Thread(target=_census, daemon=True)
        census_thread.start()

    res = ensure_resources_locked(args.codetta_dir, args.data_dir)

    results = {}
    for c in controls:
        try:
            results[c["label"]] = {"input": prepare_input(c, args.data_dir, args.tmp_dir)}
        except Exception as exc:
            results[c["label"]] = {"error": "input: %s" % exc}
            log("%s: input failed: %s" % (c["label"], exc))

    ready = [c for c in controls if "input" in results[c["label"]]]
    ready.sort(key=lambda c: results[c["label"]]["input"]["sample_bp"])
    budget_s = args.time_budget_min * 60

    def _do(c):
        remaining = budget_s - (time.time() - t_start)
        if remaining < 300:
            return c, {"skipped": "time budget"}
        try:
            r = run_codetta(c, results[c["label"]]["input"], args.codetta_dir, res, args.out_dir,
                            args.tmp_dir, args.threads, remaining)
            v = grade.evaluate_control(c, r["inferred"])
            log("%s: %s  (%.1f min) diag_ok=%s disagreements=%d" % (
                c["label"], r["inferred"], r["minutes"], v["diagnostic_ok"], len(v["disagreements"])))
            return c, {"codetta": r, "verdict": v}
        except Exception as exc:
            log("%s: codetta failed: %s" % (c["label"], exc))
            return c, {"error": "codetta: %s" % exc}

    if ready:  # the first (smallest) run builds Codetta's emission cache; run it alone
        c, r = _do(ready[0])
        results[c["label"]].update(r)
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as ex:
            futs = [ex.submit(_do, c) for c in ready[1:]]
            for f in as_completed(futs):
                c, r = f.result()
                results[c["label"]].update(r)

    verdicts, rates = [], []
    for c in controls:
        r = results[c["label"]]
        if "verdict" in r:
            verdicts.append(r["verdict"])
            mb = r["input"]["sample_bp"] / 1e6
            if mb >= 1.0:
                rates.append(r["codetta"]["minutes"] / mb)
    graded = grade.grade(verdicts, rates)

    if census_thread is not None:
        census_thread.join(timeout=max(10, budget_s + 1800 - (time.time() - t_start)))

    meta = {"eukcode_version": __version__, "codetta_commit": commit,
            "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "wall_minutes": round((time.time() - t_start) / 60, 1),
            "workers": args.workers, "hmmer_threads": args.threads}
    for r in results.values():
        if "input" in r:
            r["input"].pop("input_path", None)
    out = {"meta": meta, "predictions": graded, "results": results, "census": census_result}
    with open(os.path.join(args.out_dir, "results.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    write_results_md(os.path.join(args.out_dir, "RESULTS.md"), controls, results, graded, meta)
    log("done in %.1f min" % meta["wall_minutes"])
    shutil.rmtree(args.tmp_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
