# eukcode-scan: NCBI Datasets v2 and E-utilities helpers (public endpoints, no key).
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Find genome assemblies by taxon, download their FASTA, and look up NCBI genetic codes.

Network calls are thin wrappers; the parsing functions are pure and unit-tested.
"""

import json
import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "eukcode-scan/0.1 (MV2 of Munim, Inc.; mv2@munimtech.com)"

LEVEL_RANK = {"Complete Genome": 4, "Chromosome": 3, "Scaffold": 2, "Contig": 1}


def http_get(url, timeout=180, retries=4, pause=0.4):
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            time.sleep(pause)
            return data
        except Exception as exc:  # network errors are retried, then raised
            last = exc
            time.sleep(2.0 * (attempt + 1))
    raise last


def http_download(url, path, timeout=600, retries=3):
    tmp = path + ".part"
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(tmp, "wb") as fh:
                while True:
                    block = resp.read(1 << 20)
                    if not block:
                        break
                    fh.write(block)
            os.replace(tmp, path)
            return path
        except Exception as exc:
            last = exc
            time.sleep(3.0 * (attempt + 1))
    raise last


def summarize_report(r):
    """Flatten one Datasets genome report into the fields the screen uses."""
    info = r.get("assembly_info", {}) or {}
    stats = r.get("assembly_stats", {}) or {}
    org = r.get("organism", {}) or {}

    def _int(x):
        try:
            return int(x)
        except (TypeError, ValueError):
            return 0

    return {
        "accession": r.get("accession", ""),
        "organism": org.get("organism_name", ""),
        "tax_id": _int(org.get("tax_id")),
        "level": info.get("assembly_level", ""),
        "refseq_category": info.get("refseq_category", ""),
        "release_date": info.get("release_date", ""),
        "assembly_name": info.get("assembly_name", ""),
        "length": _int(stats.get("total_sequence_length")),
        "contig_n50": _int(stats.get("contig_n50")),
        "source": r.get("source_database", ""),
        "annotated": bool(r.get("annotation_info")),
    }


def best_assembly(summaries):
    """Prefer a reference genome, then a higher assembly level, then contig N50, then GCF."""
    if not summaries:
        return None

    def key(s):
        return (
            1 if s.get("refseq_category") == "reference genome" else 0,
            LEVEL_RANK.get(s.get("level", ""), 0),
            s.get("contig_n50", 0),
            1 if s.get("accession", "").startswith("GCF_") else 0,
            s.get("accession", ""),
        )

    return max(summaries, key=key)


def one_per_taxon(summaries):
    by = {}
    for s in summaries:
        by.setdefault(s["tax_id"], []).append(s)
    return [best_assembly(v) for _, v in sorted(by.items())]


def taxon_reports(taxon, reference_only=False, page_size=1000, max_pages=40):
    """Yield flattened reports for all current assemblies under a taxon (name or id)."""
    token = None
    for _ in range(max_pages):
        q = {"page_size": page_size}
        if reference_only:
            q["filters.reference_only"] = "true"
        if token:
            q["page_token"] = token
        url = "%s/genome/taxon/%s/dataset_report?%s" % (
            API, urllib.parse.quote(str(taxon), safe=""), urllib.parse.urlencode(q))
        data = json.loads(http_get(url).decode("utf-8"))
        for r in data.get("reports", []) or []:
            yield summarize_report(r)
        token = data.get("next_page_token")
        if not token:
            break


def resolve_taxon(taxa):
    """Return the best assembly summary for the first taxon name in `taxa` that has one."""
    if isinstance(taxa, str):
        taxa = [taxa]
    for t in taxa:
        for ref_only in (True, False):
            try:
                found = list(taxon_reports(t, reference_only=ref_only, max_pages=3))
            except Exception:
                found = []
            if found:
                best = best_assembly(found)
                best["resolved_from"] = t
                best["reference_only_query"] = ref_only
                return best
    return None


def accession_report(accession):
    url = "%s/genome/accession/%s/dataset_report" % (API, urllib.parse.quote(accession, safe=""))
    data = json.loads(http_get(url).decode("utf-8"))
    reports = data.get("reports", []) or []
    return summarize_report(reports[0]) if reports else None


def download_genome_zip(accession, dest_dir):
    os.makedirs(dest_dir, exist_ok=True)
    path = os.path.join(dest_dir, "%s.zip" % accession)
    if os.path.exists(path) and zipfile.is_zipfile(path):
        return path
    url = "%s/genome/accession/%s/download?include_annotation_type=GENOME_FASTA" % (
        API, urllib.parse.quote(accession, safe=""))
    http_download(url, path)
    if not zipfile.is_zipfile(path):
        os.remove(path)
        raise RuntimeError("download for %s is not a zip" % accession)
    return path


def fasta_members(zip_path):
    with zipfile.ZipFile(zip_path) as z:
        return [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta", ".fna.gz"))]


def extract_fasta(zip_path, out_path):
    members = fasta_members(zip_path)
    if not members:
        raise RuntimeError("no FASTA in %s" % zip_path)
    with zipfile.ZipFile(zip_path) as z, open(out_path, "wb") as out:
        for m in members:
            with z.open(m) as src:
                while True:
                    block = src.read(1 << 20)
                    if not block:
                        break
                    out.write(block)
    return out_path


def parse_taxonomy_xml(xml_bytes):
    """Parse efetch db=taxonomy XML into {tax_id: {...}} for the top-level Taxon records."""
    root = ET.fromstring(xml_bytes)
    out = {}
    for t in root.findall("./Taxon"):
        tid = t.findtext("TaxId")
        if tid is None:
            continue
        gc = t.findtext("GeneticCode/GCId")
        mgc = t.findtext("MitoGeneticCode/MGCId")
        out[int(tid)] = {
            "name": t.findtext("ScientificName") or "",
            "rank": t.findtext("Rank") or "",
            "gencode": int(gc) if gc and gc.isdigit() else None,
            "mito_gencode": int(mgc) if mgc and mgc.isdigit() else None,
            "lineage": t.findtext("Lineage") or "",
        }
    return out


def taxonomy_lookup(tax_ids, batch=200):
    out = {}
    ids = sorted(set(int(x) for x in tax_ids if x))
    for i in range(0, len(ids), batch):
        part = ids[i:i + batch]
        url = "%s/efetch.fcgi?db=taxonomy&id=%s" % (EUTILS, ",".join(str(x) for x in part))
        try:
            out.update(parse_taxonomy_xml(http_get(url, pause=0.5)))
        except Exception as exc:
            for x in part:
                out.setdefault(x, {"error": str(exc)[:200]})
    return out
