# eukcode-scan addendum A8: run KACI on the six pre-registered control genomes and grade P16-P21.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Usage (inside the compete container, after bash install_kaci.sh):

    python run_kaci.py --data-dir /data --out-dir /out --threads 2

The KACI reference (Zenodo 10.5281/zenodo.19318166, CC BY 4.0) is cached in the
shared data volume under a lock; genomes are per-job inputs and live in a temporary
folder that is deleted after each genome. Grading follows PREREGISTRATION_A8.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import tarfile
import tempfile
import time
import traceback
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from eukcode.kaci_io import STANDARD, filter_fasta, grade_control, infer_code  # noqa: E402

UA = {"User-Agent": "eukcode-scan/0.8.0 (open research pipeline)"}
DATASETS = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
            "?include_annotation_type=GENOME_FASTA")


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def peak_rss_gb() -> float | None:
    try:
        with open("/proc/self/status", encoding="ascii") as fh:
            for line in fh:
                if line.startswith("VmHWM:"):
                    return round(int(line.split()[1]) / 1024 / 1024, 3)
    except OSError:
        return None
    return None


def write_status(out: Path, status: dict) -> None:
    (out / "status.json").write_text(json.dumps(status, indent=1), encoding="utf-8")


def download(url: str, dest: Path, md5: str | None = None, tries: int = 3) -> None:
    last = None
    for attempt in range(tries):
        try:
            part = dest.with_name(dest.name + ".part")
            digest = hashlib.md5()
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=300) as r, \
                    open(part, "wb") as fh:
                while True:
                    block = r.read(1 << 20)
                    if not block:
                        break
                    fh.write(block)
                    digest.update(block)
            if md5 and digest.hexdigest() != md5:
                raise IOError(f"md5 mismatch for {url}: {digest.hexdigest()} != {md5}")
            part.replace(dest)
            return
        except Exception as error:  # noqa: BLE001 - retried, then raised
            last = error
            log(f"download attempt {attempt + 1} failed: {error}")
            time.sleep(20 * (attempt + 1))
    raise IOError(f"download failed: {url}: {last}")


def ensure_reference(data_dir: Path, ref: dict) -> Path:
    """Download and unpack the KACI reference once per data volume; return the .bin file."""
    import fcntl  # Linux only; kept here so the tests import this file on Windows

    root = data_dir / "kaci-ref"
    root.mkdir(parents=True, exist_ok=True)
    ready = root / "READY.json"
    with open(data_dir / "kaci-ref.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not ready.exists():
            tgz = root / ref["file"]
            log(f"downloading KACI reference ({ref['bytes'] / 1e9:.2f} GB)")
            download(ref["url"], tgz, ref["md5"])
            unpack = root / "ref"
            if unpack.exists():
                shutil.rmtree(unpack)
            unpack.mkdir()
            with tarfile.open(tgz, "r:gz") as tar:
                try:
                    tar.extractall(unpack, filter="data")
                except TypeError:
                    tar.extractall(unpack)
            tgz.unlink()
            files = sorted(str(p.relative_to(unpack)) for p in unpack.rglob("*") if p.is_file())
            ready.write_text(json.dumps({"md5": ref["md5"], "files": files}, indent=1), encoding="utf-8")
        fcntl.flock(lock, fcntl.LOCK_UN)
    bins = sorted((root / "ref").rglob("*.bin"))
    if len(bins) != 1:
        listing = json.loads(ready.read_text(encoding="utf-8"))["files"]
        raise RuntimeError(f"expected one .bin reference file, found {len(bins)}; files: {listing[:50]}")
    return bins[0]


def fetch_genome(acc: str, tmp: Path) -> tuple[Path, dict]:
    """Genome FASTA from NCBI Datasets, organelle records removed."""
    zpath = tmp / f"{acc}.zip"
    download(DATASETS.format(acc=acc), zpath)
    raw = tmp / f"{acc}_raw.fna"
    with zipfile.ZipFile(zpath) as z, open(raw, "wb") as dst:
        members = [m for m in z.namelist() if m.endswith("_genomic.fna")
                   and "cds_from" not in m and "rna_from" not in m]
        if not members:
            raise RuntimeError(f"{acc}: no genomic FASTA in the Datasets package: {z.namelist()}")
        for m in members:
            with z.open(m) as src:
                shutil.copyfileobj(src, dst, 1 << 20)
    zpath.unlink()
    clean = tmp / f"{acc}_genomic.fna"
    counts = filter_fasta(str(raw), str(clean))
    counts["members"] = members
    counts["bytes"] = clean.stat().st_size
    raw.unlink()
    return clean, counts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--only", default="", help="comma-separated control ids (debug)")
    args = ap.parse_args()
    out = Path(args.out_dir)
    (out / "kaci").mkdir(parents=True, exist_ok=True)
    spec = json.loads((HERE / "kaci_controls.json").read_text(encoding="utf-8"))
    controls = spec["controls"]
    if args.only:
        keep = set(args.only.split(","))
        controls = [c for c in controls if c["id"] in keep]
    status = {"stage": "reference", "started": time.time(), "genomes": {}}
    write_status(out, status)
    try:
        ref_bin = ensure_reference(Path(args.data_dir), spec["reference"])
    except Exception as error:  # noqa: BLE001 - reported, not hidden
        status.update(stage="reference-failed", error=f"{type(error).__name__}: {error}")
        write_status(out, status)
        traceback.print_exc()
        return 2
    status["reference_bin"] = str(ref_bin)
    status["reference_bin_bytes"] = ref_bin.stat().st_size

    sys.path.insert(0, str(HERE / "lib-kaci"))
    from src.libkaci_wrapper import LibKaci  # type: ignore  # noqa: E402

    lib = LibKaci(str(HERE / "lib-kaci" / "src" / "libkaci.so"))
    p = spec["params"]
    lib.set_kmer_size(int(p["kmer_size"]))
    lib.set_translation_table(STANDARD)
    lib.set_verbose(False)
    lib.set_report_query_details(False)
    raw_out = out / "kaci" / "raw"
    raw_out.mkdir(parents=True, exist_ok=True)
    lib.set_output_dir(str(raw_out))
    lib.set_max_query_count(int(p["max_query_count"]))
    lib.set_n_threads(int(args.threads))
    status["stage"] = "loading reference"
    write_status(out, status)
    t0 = time.time()
    lib.load_ref_bin(str(ref_bin))
    status["reference_load_s"] = round(time.time() - t0, 1)
    status["peak_rss_gb_after_load"] = peak_rss_gb()
    log(f"reference loaded in {status['reference_load_s']} s, peak RSS {status['peak_rss_gb_after_load']} GB")

    grades = {}
    tmp = Path(tempfile.mkdtemp(prefix="kaci-"))
    for c in controls:
        acc = c["accession"]
        status["stage"] = f"genome {c['id']} {acc}"
        write_status(out, status)
        g = {"id": c["id"], "accession": acc, "organism": c["organism"], "table": c["table"],
             "must": c.get("must") or {}}
        try:
            genome, counts = fetch_genome(acc, tmp)
            g["fasta"] = counts
            lib.set_input_dir(str(tmp))
            before = {q.name for q in raw_out.glob("*.json")}
            t1 = time.time()
            lib.process_query(str(genome))
            g["kaci_s"] = round(time.time() - t1, 1)
            new = sorted(q for q in raw_out.glob("*.json") if q.name not in before)
            if len(new) != 1:
                raise RuntimeError(f"expected one new KACI output, found {[q.name for q in new]}")
            dest = out / "kaci" / f"{acc}.kaci.json"
            new[0].replace(dest)
            result = json.loads(dest.read_text(encoding="utf-8"))
            code = infer_code(result, int(p["min_n"]), float(p["high_conf"]))
            g["inferred_code"] = code
            g.update(grade_control(code, c))
            genome.unlink()
        except Exception as error:  # noqa: BLE001 - one genome's failure is recorded
            g["pass"] = None
            g["error"] = f"{type(error).__name__}: {error}"
            traceback.print_exc()
        g["peak_rss_gb"] = peak_rss_gb()
        grades[c["id"]] = g
        status["genomes"][c["id"]] = {k: g.get(k) for k in ("pass", "error", "kaci_s", "peak_rss_gb")}
        write_status(out, status)
        log(f"{c['id']} {acc}: pass={g.get('pass')} {g.get('error', '')}")
    shutil.rmtree(tmp, ignore_errors=True)
    lib.free_resources()

    verdicts = [g.get("pass") for g in grades.values()]
    complete = len(grades) == 6 and all(v is not None for v in verdicts)
    summary = {
        "addendum": "A8",
        "complete": complete,
        "adopt_kaci_triage": bool(complete and all(verdicts)),
        "passed": [k for k, g in grades.items() if g.get("pass") is True],
        "failed": [k for k, g in grades.items() if g.get("pass") is False],
        "not_graded": [k for k, g in grades.items() if g.get("pass") is None],
        "reference_load_s": status.get("reference_load_s"),
        "peak_rss_gb": peak_rss_gb(),
        "controls": grades,
    }
    (out / "a8_grade.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    lines = ["# A8 KACI control grades", "",
             "| id | organism | table | must | called | unexpected sense | stops called | pass |",
             "|---|---|---|---|---|---|---|---|"]
    for k, g in grades.items():
        cmp = g.get("comparison") or {}
        called = ", ".join(f"{cd}={v['called']}" for cd, v in (g.get("must") or {}).items()
                           if isinstance(v, dict))
        lines.append(f"| {k} | {g['organism']} | {g['table']} | "
                     f"{json.dumps(c_must(controls, k))} | {called} | "
                     f"{cmp.get('sense_unexpected', '')} | {cmp.get('stop_called', '')} | "
                     f"{g.get('pass')} {g.get('error', '')} |")
    lines += ["", f"adopt_kaci_triage: {summary['adopt_kaci_triage']} (complete: {complete})"]
    (out / "A8_RESULT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    status["stage"] = "done"
    write_status(out, status)
    return 0


def c_must(controls: list[dict], cid: str) -> dict:
    for c in controls:
        if c["id"] == cid:
            return c.get("must") or {}
    return {}


if __name__ == "__main__":
    sys.exit(main())
