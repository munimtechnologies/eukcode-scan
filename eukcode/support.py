# eukcode-scan: check 2, where a candidate codon's evidence comes from.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Trace Codetta's codon observations back to the input sequences (PREREGISTRATION.md, A4).

Codetta 2.0 infers each codon's meaning from the Pfam profile columns that the codon is
aligned to. Its alignment output file has one line per aligned codon:

    codon_index, piece.frame+position, codon, piece, frame, query_position, profile, evalue, hmm_position

(written by `GeneticCode.process_hmmscan_results` in codetta.py). Pieces are the input FASTA
records in file order, after `processing_genome` has cut every record longer than 100,000 nt
into halves (first half in place, second half appended to that record's list, repeated until
no piece is longer). Frames 0-2 read the forward strand from offset 0, 1, 2; frames 3-5 read
the reverse complement of seq, seq[:-1] and seq[:-2].

This module rebuilds that mapping, recomputes every codon from the coordinates as a
self-check, counts each candidate codon's observations per input sequence, and grades check 2
from the best homology hit of the sequences that carry most of the evidence.
"""

import gzip
import io

PIECE_MAX = 100_000
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

TOP_N = 20
MIN_COVER = 0.2
MIN_CLASSIFIED = 3
EUK_CLASSES = ("same_clade", "other_eukaryote")
CLASSIFIED = ("same_clade", "other_eukaryote", "bacteria", "archaea", "virus")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def codetta_pieces(lengths, piece_max=PIECE_MAX):
    """[(record_index, start, length)] in Codetta's piece-number order.

    Records of length 0 produce no piece (Codetta never sees a sequence for them).
    """
    out = []
    for rec, n in enumerate(lengths):
        if n <= 0:
            continue
        segs = [(0, n)]
        while any(l > piece_max for _s, l in segs):
            long_inds = [i for i, (_s, l) in enumerate(segs) if l > piece_max]
            for li in long_inds:
                s, l = segs[li]
                h = l // 2
                segs[li] = (s, h)
                segs.append((s + h, l - h))
        out.extend((rec, s, l) for s, l in segs)
    return out


def codon_at(piece_seq, frame, qpos):
    """(codon, (start, end)) that Codetta reads at 0-based protein position qpos of a frame.

    start/end are forward-strand coordinates on the piece; codon is read on the frame's strand.
    Returns ('', (start, end)) if the position falls outside the piece.
    """
    if frame < 3:
        s = frame + 3 * qpos
        e = s + 3
        if s < 0 or e > len(piece_seq):
            return "", (s, e)
        return piece_seq[s:e].upper(), (s, e)
    lp = len(piece_seq) - (frame - 3)
    e = lp - 3 * qpos
    s = e - 3
    if s < 0 or e > len(piece_seq):
        return "", (s, e)
    return revcomp(piece_seq[s:e]).upper(), (s, e)


def parse_alignment_line(line):
    """One Codetta alignment output line -> dict, or None for anything that does not parse."""
    parts = line.rstrip("\r\n").split(",")
    if len(parts) < 9:
        return None
    try:
        return {"codon_index": int(parts[0]), "codon": parts[2].strip().upper(), "piece": int(parts[3]),
                "frame": int(parts[4]), "qpos": int(parts[5]), "profile": ",".join(parts[6:-2]),
                "evalue": float(parts[-2]), "hmm_pos": int(parts[-1])}
    except ValueError:
        return None


def open_text(path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="ascii", errors="replace")
    return open(path, "r", encoding="ascii", errors="replace")


def iter_alignment(path):
    with open_text(path) as fh:
        for line in fh:
            row = parse_alignment_line(line)
            if row is not None:
                yield row


def read_pieces_file(path):
    """Codetta's <prefix>.sequence_pieces.fna -> {piece_number: sequence}."""
    out, cur = {}, None
    with open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">piece_"):
                cur = int(line[len(">piece_"):])
                out[cur] = []
            elif cur is not None and line:
                out[cur].append(line)
    return {k: "".join(v) for k, v in out.items()}


def piece_sequences(seqs, pieces):
    return [seqs[rec][s:s + l] for rec, s, l in pieces]


def tally(rows, pieces, piece_seqs, codons):
    """Count observations of each codon in `codons` per input record, checking every row.

    rows: parsed alignment rows; pieces: codetta_pieces(); piece_seqs: piece_sequences().
    Returns {'rows', 'mismatches', 'mismatch_examples', 'by_codon': {codon: {'total',
    'per_record': {rec: {'obs', 'profiles', 'start', 'end'}}}}} with record coordinates.
    """
    want = set(codons)
    by = {c: {"total": 0, "per_record": {}} for c in want}
    n = bad = 0
    examples = []
    for r in rows:
        n += 1
        p = r["piece"]
        if p < 0 or p >= len(pieces):
            bad += 1
            if len(examples) < 5:
                examples.append({"row": r, "why": "piece out of range"})
            continue
        got, (s, e) = codon_at(piece_seqs[p], r["frame"], r["qpos"])
        if got != r["codon"]:
            bad += 1
            if len(examples) < 5:
                examples.append({"row": r, "why": "codon %r at piece coordinates" % got})
            continue
        if r["codon"] not in want:
            continue
        rec, off, _l = pieces[p]
        d = by[r["codon"]]
        d["total"] += 1
        pr = d["per_record"].setdefault(rec, {"obs": 0, "profiles": set(), "start": off + s, "end": off + e})
        pr["obs"] += 1
        pr["profiles"].add(r["profile"])
        pr["start"] = min(pr["start"], off + s)
        pr["end"] = max(pr["end"], off + e)
    for d in by.values():
        for pr in d["per_record"].values():
            pr["profiles"] = len(pr["profiles"])
    return {"rows": n, "mismatches": bad, "mismatch_examples": examples, "by_codon": by}


def top_records(per_record, n=TOP_N):
    """Record indices carrying the most observations (ties: lower index first)."""
    return [rec for rec, _ in sorted(per_record.items(), key=lambda kv: (-kv[1]["obs"], kv[0]))[:n]]


def region(start, end, rec_len, pad=300, max_len=3000):
    """Window around the aligned codons to send to the homology search."""
    s, e = max(0, start - pad), min(rec_len, end + pad)
    if e - s > max_len:
        mid = (start + end) // 2
        s = max(0, mid - max_len // 2)
        e = min(rec_len, s + max_len)
        s = max(0, e - max_len)
    return s, e


def lineage_parts(lineage):
    return [x.strip() for x in (lineage or "").split(";") if x.strip() and x.strip() != "cellular organisms"]


def major_clade(lineage, depth=2):
    """The two lineage names below Eukaryota, e.g. ('Sar', 'Alveolata'); None if not eukaryotic."""
    parts = lineage_parts(lineage)
    if not parts or parts[0] != "Eukaryota":
        return None
    return tuple(parts[1:1 + depth])


def classify_hit(hit_lineage, species_lineage):
    """Class of one best hit relative to the screened species (A4.3)."""
    if hit_lineage is None:
        return "no_hit"
    parts = lineage_parts(hit_lineage)
    if not parts:
        return "unknown"
    top = {"Bacteria": "bacteria", "Archaea": "archaea", "Viruses": "virus"}.get(parts[0])
    if top:
        return top
    if parts[0] != "Eukaryota":
        return "unknown"
    sp, h = major_clade(species_lineage), major_clade(hit_lineage)
    if sp and h and h[:len(sp)] == sp:
        return "same_clade"
    return "other_eukaryote"


def grade_check2(obs, classes, total, min_cover=MIN_COVER, min_classified=MIN_CLASSIFIED):
    """Grade check 2 for one candidate codon.

    obs: {record: observations} for the searched records; classes: {record: (class, clade)};
    total: all observations of the codon in the sample.
    """
    cls = {r: classes.get(r, ("no_hit", None)) for r in obs}
    classified = {r: o for r, o in obs.items() if cls[r][0] in CLASSIFIED}
    cobs = sum(classified.values())

    def share(names):
        return sum(o for r, o in classified.items() if cls[r][0] in names) / cobs if cobs else 0.0

    foreign = {}
    for r, o in classified.items():
        if cls[r][0] == "other_eukaryote":
            k = " / ".join(cls[r][1] or ("?",))
            foreign[k] = foreign.get(k, 0) + o
    top_foreign = max(foreign.items(), key=lambda kv: (kv[1], kv[0])) if foreign else ("", 0)
    out = {
        "searched": len(obs), "classified": len(classified),
        "coverage": round(cobs / total, 4) if total else 0.0,
        "eukaryotic_share": round(share(EUK_CLASSES), 4),
        "same_clade_share": round(share(("same_clade",)), 4),
        "outside_clade_share": round(1 - share(("same_clade",)), 4) if cobs else 0.0,
        "prokaryote_or_virus_share": round(share(("bacteria", "archaea", "virus")), 4),
        "top_foreign_clade": top_foreign[0],
        "top_foreign_share": round(top_foreign[1] / cobs, 4) if cobs else 0.0,
        "top1_share_of_total": round(max(obs.values()) / total, 4) if obs and total else 0.0,
    }
    if out["classified"] < min_classified or (cobs / total if total else 0) < min_cover:
        out["outcome"] = "inconclusive"
    elif out["eukaryotic_share"] <= 0.5:
        out["outcome"] = "fail"
    elif out["same_clade_share"] > 0.5:
        out["outcome"] = "pass"
    elif out["top_foreign_share"] > 0.5:
        out["outcome"] = "fail"
    else:
        out["outcome"] = "pass_unplaced"
    out["concentrated"] = out["top1_share_of_total"] > 0.5
    return out


def origin_shares(per_record, origin_of):
    """For a constructed mix: share of a codon's observations by source label."""
    tot = sum(v["obs"] for v in per_record.values())
    out = {}
    for rec, v in per_record.items():
        k = origin_of(rec)
        out[k] = out.get(k, 0) + v["obs"]
    return {k: round(v / tot, 4) for k, v in sorted(out.items())} if tot else {}
