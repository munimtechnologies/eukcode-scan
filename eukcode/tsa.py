# eukcode-scan: transcriptome (INSDC TSA) tranche of the screen.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Find, choose and sample protist transcriptome assemblies (INSDC Transcriptome Shotgun Assembly).

Many protist lineages, including uncultured single-cell isolates, have a transcriptome but
no genome assembly. Transcripts are mostly coding sequence, so a small sample gives Codetta
far more aligned codons per megabase than a genome does.

- The census uses the ENA portal API (`result=tsa_set`, `tax_tree(<taxid>)`).
- Each set is fetched from NCBI GenBank
  (https://ftp.ncbi.nlm.nih.gov/genbank/tsa/<first letter>/tsa.<CODE>.1.fsa_nt.gz), or from
  the ENA mirror if that fails.
- Prior art fixes the novelty tiers. Swart et al. 2016 (Cell 166:691) screened every MMETSP
  transcriptome for stop-codon reassignments, so sets derived from MMETSP are tier
  `mmetsp`. Sets first public before 2016-09-01 are tier `pre2016`. Everything newer is
  tier `new`.
- MMETSP-derived sets are recognised by the word MMETSP in their description or keywords,
  or (addendum A6) by belonging to one of the two INSDC projects that hold the 2021 EBI/ELIXIR
  re-assemblies of MMETSP. Those 406 sets (prefix HB) say only "TSA: <species>" or
  "TPA: <species>", so the word test alone never caught them.
- Trinity-style assemblies hold many isoforms of one gene. Sampling keeps only the longest
  isoform per gene, so no gene is counted twice, and a disjoint resample (check 1) shares no
  gene with the first sample.
"""

import csv
import io
import re
import urllib.parse

from . import fasta

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_TSA = "https://ftp.ncbi.nlm.nih.gov/genbank/tsa"
ENA_TSA = "https://ftp.ebi.ac.uk/pub/databases/ena/tsa/public"
ENA_FIELDS = ["accession", "scientific_name", "tax_id", "base_count", "first_public",
              "environmental_sample", "dataclass", "keywords", "description", "study_accession"]
SWART_CUTOFF = "2016-09-01"
TIER_ORDER = ["new", "pre2016", "mmetsp"]
# INSDC projects titled "Assembly and functional analysis of the Marine Microbial Eukaryotic
# Transcriptome Project data set" (ENA study search, 2026-11-19): 406 TSA sets between them.
MMETSP_STUDIES = frozenset({"PRJEB37117", "PRJEB37159"})

_SET_RE = re.compile(r"^([A-Z]{4,6})(\d{2})\d+(?:\.\d+)?$")
_GENE_RES = [
    re.compile(r"(\S*?c\d+_g\d+)_i\d+"),     # Trinity >= 2.0: TRINITY_DN12_c0_g1_i3, foca_c0_g1_i1
    re.compile(r"(comp\d+_c\d+)_seq\d+"),    # Trinity < 2.0: comp12_c0_seq3
]


def ena_search_url(taxid, fields=ENA_FIELDS):
    q = {"result": "tsa_set", "query": "tax_tree(%d)" % int(taxid), "fields": ",".join(fields),
         "format": "tsv", "limit": "0"}
    return "%s?%s" % (ENA_PORTAL, urllib.parse.urlencode(q))


def taxid_search_url(name):
    term = '"%s"[Scientific Name]' % name
    return "%s/esearch.fcgi?db=taxonomy&retmode=json&term=%s" % (EUTILS, urllib.parse.quote(term))


def parse_tsv(text):
    """ENA portal TSV -> list of dicts (empty list for an empty or header-only reply)."""
    text = text.strip("\n")
    if not text:
        return []
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def set_code(accession):
    """'GABU01000000' -> ('GABU', '01'); raises ValueError for anything else."""
    m = _SET_RE.match((accession or "").strip())
    if not m:
        raise ValueError("not a TSA set accession: %r" % (accession,))
    return m.group(1), m.group(2)


def fasta_urls(accession):
    """Candidate download URLs for a TSA set, NCBI first, then the ENA mirror."""
    code, ver = set_code(accession)
    return ["%s/%s/tsa.%s.1.fsa_nt.gz" % (NCBI_TSA, code[0], code),
            "%s/%s/%s%s.fasta.gz" % (ENA_TSA, code[:3].lower(), code, ver)]


def _int(x):
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return 0


def is_mmetsp(row):
    """True if an ENA tsa_set row is derived from MMETSP (word test or re-assembly project)."""
    text = " ".join([row.get("description") or "", row.get("keywords") or ""]).upper()
    study = (row.get("study_accession") or "").strip().upper()
    return "MMETSP" in text or study in MMETSP_STUDIES


def flags(row):
    """Return {'environmental', 'tpa', 'mmetsp'} booleans for one ENA tsa_set row."""
    kw = [k.strip().upper() for k in (row.get("keywords") or "").replace(";", ",").split(",")]
    return {
        "environmental": (row.get("environmental_sample") or "").strip().lower() == "true",
        "tpa": (row.get("dataclass") or "").strip().upper() == "TPA" or "TPA" in kw,
        "mmetsp": is_mmetsp(row),
    }


def set_tier(row):
    f = flags(row)
    if f["mmetsp"]:
        return "mmetsp"
    if (row.get("first_public") or "9999") < SWART_CUTOFF:
        return "pre2016"
    return "new"


def retier_species(species, sets):
    """Recompute species tiers from a census' per-set rows with the current rule (A6.2).

    species: rows with 'tax_id' (and 'tier'); sets: per-set rows with 'tax_id',
    'description', 'keywords', 'study_accession', 'first_public' and an optional
    'skip_reason'. Only usable sets (no skip_reason) count, as in choose_sets. Returns new
    species dicts with 'tier' replaced and the old value kept as 'tier_before_a6'.
    """
    by = {}
    for s in sets:
        if s.get("skip_reason"):
            continue
        by.setdefault(_int(s.get("tax_id")), []).append(set_tier(s))
    out = []
    for sp in species:
        tiers = by.get(_int(sp.get("tax_id")))
        row = dict(sp)
        row["tier_before_a6"] = sp.get("tier", "")
        if tiers:
            row["tier"] = max(tiers, key=TIER_ORDER.index)
        out.append(row)
    return out


def choose_sets(rows, min_bp=1_000_000):
    """One TSA set per tax_id. Returns (chosen, skipped).

    Environmental (metatranscriptome) sets and sets under min_bp are skipped with a reason.
    Among the rest, a set not derived from MMETSP is preferred, then the largest, then the
    newest, then the accession. The species tier is the "oldest" tier of all its usable sets
    (any MMETSP set makes it 'mmetsp'; any pre-2016 set makes it 'pre2016'), because a species
    screened before stays screened even if a newer set exists.
    """
    by, skipped, seen = {}, [], set()
    for r in rows:
        acc = (r.get("accession") or "").strip()
        reason = ""
        try:
            set_code(acc)
        except ValueError:
            reason = "bad accession"
        if not reason and acc in seen:
            reason = "duplicate accession"
        f = flags(r)
        if not reason and f["environmental"]:
            reason = "environmental sample"
        if not reason and _int(r.get("base_count")) < min_bp:
            reason = "under %d bp" % min_bp
        if not reason and not _int(r.get("tax_id")):
            reason = "no tax_id"
        if reason:
            s = dict(r)
            s["skip_reason"] = reason
            skipped.append(s)
            continue
        seen.add(acc)
        by.setdefault(_int(r["tax_id"]), []).append(r)
    chosen = []
    for tid in sorted(by):
        sets = by[tid]
        best = max(sets, key=lambda r: (0 if flags(r)["mmetsp"] else 1, _int(r.get("base_count")),
                                        r.get("first_public", ""), r.get("accession", "")))
        tiers = [set_tier(r) for r in sets]
        out = dict(best)
        out["tier"] = max(tiers, key=TIER_ORDER.index)
        out["n_sets"] = len(sets)
        out["tax_id"] = str(tid)
        chosen.append(out)
    return chosen, skipped


def gene_key(header):
    """Gene id of a transcript from its FASTA header, or None if no assembler pattern matches."""
    for rx in _GENE_RES:
        m = rx.search(header or "")
        if m:
            return m.group(1)
    return None


def iter_records(path):
    """Yield (full_header_without_'>', sequence) from a FASTA or FASTA.gz file."""
    handle = fasta._open_text(path)
    header, parts = None, []
    try:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts).upper()
                header, parts = line[1:], []
            else:
                parts.append(line)
        if header is not None:
            yield header, "".join(parts).upper()
    finally:
        handle.close()


def transcript_index(path, min_len=300, max_n_frac=0.1):
    """Usable transcripts, longest isoform per gene: [(name, record_index, length)] by index.

    A transcript is usable if it is at least min_len long and at most max_n_frac N. Headers
    without a recognisable gene id are treated as their own gene.
    """
    best = {}
    for idx, (header, seq) in enumerate(iter_records(path)):
        n = len(seq)
        if n < min_len or seq.count("N") / n > max_n_frac:
            continue
        name = header.split()[0] if header.split() else "t%d" % idx
        key = gene_key(header) or ("#%d" % idx)
        cur = best.get(key)
        if cur is None or n > cur[2]:
            best[key] = (name, idx, n)
    return sorted(best.values(), key=lambda x: x[1])


def sample_transcripts(path, cap_bp, seed=0, disjoint_from=None, min_len=300, max_n_frac=0.1):
    """Seeded sample of whole transcripts up to cap_bp. Returns (records, n_usable).

    Uses eukcode.fasta.choose, so the rules match the genome sampler: same input, cap and seed
    give the same sample; disjoint_from excludes the transcripts that seed's sample used.
    """
    idx = transcript_index(path, min_len, max_n_frac)
    pool = list(range(len(idx)))
    if disjoint_from is not None:
        first = set(fasta.choose([x[2] for x in idx], cap_bp, disjoint_from))
        pool = [i for i in pool if i not in first]
    sub = [idx[i] for i in pool]
    picked = [sub[i] for i in fasta.choose([x[2] for x in sub], cap_bp, seed)]
    want = {rec: name for name, rec, _n in picked}
    out = {}
    for rec, (_h, seq) in enumerate(iter_records(path)):
        if rec in want:
            out[rec] = (want[rec], seq)
    return [out[rec] for _name, rec, _n in picked], len(idx)


def species_row(chosen, taxinfo, group):
    """Merge one chosen set with NCBI taxonomy into the row format the screen writes."""
    t = taxinfo.get(_int(chosen.get("tax_id")), {}) or {}
    return {
        "group": group,
        "accession": chosen.get("accession", ""),
        "organism": chosen.get("scientific_name", ""),
        "tax_id": chosen.get("tax_id", ""),
        "length": chosen.get("base_count", ""),
        "first_public": chosen.get("first_public", ""),
        "tier": chosen.get("tier", ""),
        "n_sets": chosen.get("n_sets", 1),
        "study_accession": chosen.get("study_accession", ""),
        "ncbi_gencode": "" if t.get("gencode") is None else str(t.get("gencode")),
        "lineage": t.get("lineage", ""),
    }


def plan_shards(species, per_shard=12):
    """Deterministic shards: tier 'new' first, then 'pre2016', then 'mmetsp'; within a tier by
    (group, organism, accession)."""
    if per_shard < 1:
        raise ValueError("per_shard must be >= 1")
    rows = sorted(species, key=lambda r: (TIER_ORDER.index(r.get("tier", "new")) if r.get("tier", "new") in TIER_ORDER else 9,
                                          r.get("group", ""), r.get("organism", ""), r.get("accession", "")))
    return [rows[i:i + per_shard] for i in range(0, len(rows), per_shard)]
