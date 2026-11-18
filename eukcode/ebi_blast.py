# eukcode-scan: best homology hit of a nucleotide region via EMBL-EBI Job Dispatcher NCBI BLAST+.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""blastx of a few sequences against UniProtKB through the public EBI REST service.

REST: POST {BASE}/run (form fields), GET {BASE}/status/<job>, GET {BASE}/result/<job>/json.
The JSON follows EBI's SSS schema (github.com/ebi-jdispatcher/sss_json_schema): a 'hits' list,
each hit with hit_acc, hit_uni_os, hit_uni_ox (taxid as a string) and hit_hsps. Only a few
dozen searches per candidate are sent, at most a handful at a time (EBI fair use).
"""

import json
import re
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
UA = "eukcode-scan (MV2 of Munim, Inc.; mv2@munimtech.com)"
EMAIL = "mv2@munimtech.com"
_OX = re.compile(r"\bOX=(\d+)")


def run_params(seq, title="", database="uniprotkb", exp="1e-5", alignments=5, exclude_taxids=(),
               email=EMAIL):
    p = {"email": email, "program": "blastx", "stype": "dna", "sequence": seq, "database": database,
         "exp": exp, "alignments": str(alignments), "scores": str(alignments), "title": title[:60]}
    ex = [str(int(t)) for t in exclude_taxids if str(t).strip().isdigit()]
    if ex:
        p["negative_taxids"] = ",".join(ex)
    return p


def best_hit(result):
    """Best hit of an SSS JSON result -> dict, or None if there is no hit."""
    hits = sorted(result.get("hits") or [], key=lambda h: int(h.get("hit_num", 0) or 0))
    if not hits:
        return None
    h = hits[0]
    hsps = sorted(h.get("hit_hsps") or [], key=lambda x: -float(x.get("hsp_bit_score", 0) or 0))
    hsp = hsps[0] if hsps else {}
    ox = str(h.get("hit_uni_ox") or "").strip()
    if not ox.isdigit():
        m = _OX.search(h.get("hit_def", "") or "") or _OX.search(h.get("hit_desc", "") or "")
        ox = m.group(1) if m else ""
    return {"acc": h.get("hit_acc", ""), "id": h.get("hit_id", ""),
            "organism": h.get("hit_uni_os") or h.get("hit_os") or "",
            "taxid": int(ox) if ox else None,
            "desc": (h.get("hit_uni_de") or h.get("hit_desc") or "")[:120],
            "evalue": hsp.get("hsp_expect"), "bits": hsp.get("hsp_bit_score"),
            "identity": hsp.get("hsp_identity"), "align_len": hsp.get("hsp_align_len")}


def _req(url, data=None, timeout=120):
    req = urllib.request.Request(url, data=data, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


def submit(params):
    return _req(BASE + "/run", data=urllib.parse.urlencode(params).encode("ascii")).strip()


def wait(job, max_wait_s=1800, poll_s=15):
    t0 = time.time()
    while time.time() - t0 < max_wait_s:
        status = _req("%s/status/%s" % (BASE, job)).strip()
        if status == "FINISHED":
            return status
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            return status
        time.sleep(poll_s)
    return "TIMEOUT"


def search(seq, title="", exclude_taxids=(), database="uniprotkb", retries=2):
    """Submit, wait and return {'status', 'job', 'hit'} for one nucleotide sequence."""
    last = ""
    for attempt in range(retries + 1):
        try:
            job = submit(run_params(seq, title, database=database, exclude_taxids=exclude_taxids))
            status = wait(job)
            if status != "FINISHED":
                last = status
                continue
            result = json.loads(_req("%s/result/%s/json" % (BASE, job)))
            return {"status": "ok", "job": job, "hit": best_hit(result)}
        except Exception as exc:  # network trouble: retry, then report
            last = str(exc)[:200]
            time.sleep(10 * (attempt + 1))
    return {"status": "error: %s" % last, "job": "", "hit": None}
