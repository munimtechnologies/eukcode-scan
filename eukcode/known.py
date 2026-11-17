# eukcode-scan: catalogue of documented nuclear codon reassignments, and candidate classification.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
"""Which nuclear reassignments are already documented, and where.

The screen compares every Codetta call with the code NCBI Taxonomy assigns to the species.
A call that differs is a *candidate*. Before any check is run, each candidate is placed in a
pre-registered class (PREREGISTRATION.md, addendum A1) by looking it up here:

  known_taxon      the same codon -> amino acid change is documented for a taxon in this
                   species' lineage (NCBI simply does not record it for the species)
  known_in_clade   documented in the same major clade, but not for this taxon
                   (a possible new independent origin or a new taxon with a known code)
  known_elsewhere  documented only in a different major clade
  unreported_type  this codon -> amino acid change is not documented in any eukaryotic
                   nuclear code in this catalogue

The catalogue was compiled on 2026-11-17 from these reviews and papers (refs therein):
  [McGowan2024] McGowan et al., PLoS Genet 20:e1011512 (2024), Introduction and refs.
  [Panek2017]   Panek et al., BMC Biol 15:8 (2017), doi:10.1186/s12915-017-0353-y.
  [Nenarokova2025] Nenarokova et al./Lukes lab, mBio 16:e00885-25 (2025), Blastocrithidia.
  [NCBI]        NCBI translation tables, https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
It is a starting point, not ground truth: every non-'known_taxon' candidate still gets a
manual literature check, including the 2025 Current Biology review (S0960-9822(25)01302-8).
"""

from .codes import CODONS, NONCALL, compare, known_tables, normalize_codon

# Each entry: id, clade (a name that appears in NCBI lineages), taxa (names that appear in
# NCBI lineages or as the genus), changes (codon -> amino acid), source, note.
CATALOGUE = [
    # --- ciliates (Ciliophora) ---
    {"id": "ciliate-UAR-Q", "clade": "Ciliophora", "taxa": ["Ciliophora"],
     "changes": {"TAA": "Q", "TAG": "Q"}, "source": "McGowan2024; NCBI table 6",
     "note": "the most common ciliate code, many independent lineages (Tetrahymena, Oxytricha, Palmarella)"},
    {"id": "peritrich-UAR-E", "clade": "Ciliophora", "taxa": ["Peritrichia"],
     "changes": {"TAA": "E", "TAG": "E"}, "source": "McGowan2024; NCBI table 30", "note": ""},
    {"id": "mesodinium-UAR-Y", "clade": "Ciliophora", "taxa": ["Mesodinium", "Mesodiniidae"],
     "changes": {"TAA": "Y", "TAG": "Y"}, "source": "McGowan2024; NCBI table 29", "note": ""},
    {"id": "euplotes-UGA-C", "clade": "Ciliophora",
     "taxa": ["Euplotidae", "Euplotes", "Euplotoides", "Moneuplotes", "Euplotopsis"],
     "changes": {"TGA": "C"}, "source": "McGowan2024; NCBI table 10", "note": ""},
    {"id": "heterotrich-UGA-W", "clade": "Ciliophora", "taxa": ["Blepharisma", "Condylostentor"],
     "changes": {"TGA": "W"}, "source": "McGowan2024", "note": "NCBI has no table for this code"},
    {"id": "karyorelict-context", "clade": "Ciliophora",
     "taxa": ["Karyorelictea", "Condylostoma", "Plagiopyla", "Parduczia"],
     "changes": {"TAA": "Q", "TAG": "Q", "TGA": "W"}, "source": "McGowan2024; NCBI tables 27, 28",
     "note": "all three stops reassigned, termination by context"},
    {"id": "ciliate-UAA-K-UAG-E", "clade": "Ciliophora", "taxa": [],
     "changes": {"TAA": "K", "TAG": "E"},
     "source": "McGowan2024 ref 10; PLoS Genet 2023 e1010913",
     "note": "uncultured ciliate; the first code with UAA and UAG meaning different amino acids"},
    {"id": "ciliate-UAG-Q", "clade": "Ciliophora",
     "taxa": ["Protocruzia", "Hartmannula", "Trochilia"],
     "changes": {"TAG": "Q"}, "source": "McGowan2024; NCBI table 15", "note": "UAA kept as stop"},
    {"id": "phyllopharyngea-UAG-L", "clade": "Ciliophora", "taxa": [],
     "changes": {"TAG": "L"}, "source": "McGowan2024 (TARA MAGs ARC_274, ARC_306, SOC_66)",
     "note": "uncultured phyllopharyngean clade"},
    # --- Discoba ---
    {"id": "blastocrithidia", "clade": "Kinetoplastea", "taxa": ["Blastocrithidia"],
     "changes": {"TAA": "E", "TAG": "E", "TGA": "W"},
     "source": "Nenarokova2025; NCBI table 31", "note": "UAA also stop; sister genus Obscuromonas is standard"},
    # --- Metamonada ---
    {"id": "hexamitin-UAR-Q", "clade": "Metamonada", "taxa": ["Hexamitinae", "Hexamita"],
     "changes": {"TAA": "Q", "TAG": "Q"}, "source": "McGowan2024 ref 33; NCBI table 6", "note": ""},
    {"id": "streblomastix-UAR-Q", "clade": "Metamonada", "taxa": ["Streblomastix"],
     "changes": {"TAA": "Q", "TAG": "Q"}, "source": "McGowan2024 ref 34", "note": "oxymonad"},
    {"id": "iotanema-UAG-Q", "clade": "Metamonada", "taxa": ["Iotanema"],
     "changes": {"TAG": "Q"}, "source": "Panek2017", "note": "fornicate; UAA kept as stop"},
    # --- green algae ---
    {"id": "ulvophyte-UAR-Q", "clade": "Chlorophyta",
     "taxa": ["Dasycladales", "Cladophorales", "Acetabularia", "Batophora"],
     "changes": {"TAA": "Q", "TAG": "Q"}, "source": "McGowan2024 ref 32; NCBI table 6",
     "note": "Ulvophyceae"},
    {"id": "scotinosphaera-UAG-S", "clade": "Chlorophyta", "taxa": ["Scotinosphaera", "Scotinosphaerales"],
     "changes": {"TAG": "S"}, "source": "McGowan2024 ref 29", "note": "UAA kept as stop"},
    # --- Rhizaria ---
    {"id": "rhizarian-UAG-L", "clade": "Rhizaria", "taxa": [],
     "changes": {"TAG": "L"}, "source": "Panek2017", "note": "uncultured rhizarian 'exLh'"},
    # --- Opisthokonta: aphelids and yeasts ---
    {"id": "aphelid-UAR-Q", "clade": "Aphelida", "taxa": ["Amoeboaphelidium"],
     "changes": {"TAA": "Q", "TAG": "Q"}, "source": "McGowan2024 ref 35", "note": ""},
    {"id": "yeast-CUG-S", "clade": "Saccharomycotina",
     "taxa": ["Serinales", "Debaryomycetaceae", "Metschnikowiaceae", "Ascoidea", "Saccharomycopsis"],
     "changes": {"CTG": "S"}, "source": "NCBI table 12", "note": "CUG-Ser1 and CUG-Ser2 clades"},
    {"id": "yeast-CUG-A", "clade": "Saccharomycotina", "taxa": ["Pachysolen"],
     "changes": {"CTG": "A"}, "source": "NCBI table 26", "note": ""},
]

CLASS_ORDER = ["known_taxon", "known_in_clade", "known_elsewhere", "unreported_type"]


def lineage_names(lineage, organism=""):
    """Set of taxon names from an NCBI lineage string ('a; b; c') plus the genus and species."""
    names = set()
    for part in (lineage or "").split(";"):
        part = part.strip()
        if part:
            names.add(part)
    organism = (organism or "").strip()
    if organism:
        names.add(organism)
        names.add(organism.split()[0])
    return names


def classify_change(codon, aa, lineage, organism=""):
    """Classify one codon -> amino acid change. Returns {'class': ..., 'entries': [ids]}."""
    codon = normalize_codon(codon)
    aa = aa.upper()
    names = lineage_names(lineage, organism)
    same = [e for e in CATALOGUE if e["changes"].get(codon) == aa]
    hit = [e["id"] for e in same if names.intersection(e["taxa"])]
    if hit:
        return {"class": "known_taxon", "entries": hit}
    hit = [e["id"] for e in same if e["clade"] in names]
    if hit:
        return {"class": "known_in_clade", "entries": hit}
    if same:
        return {"class": "known_elsewhere", "entries": [e["id"] for e in same]}
    return {"class": "unreported_type", "entries": []}


def candidates(inferred, ncbi_table, lineage="", organism=""):
    """All calls that differ from the NCBI-assigned code, each with its pre-registered class.

    If NCBI gives no table, or one eukcode does not model, the standard code is used and every
    candidate carries 'table_fallback': True.
    """
    fallback = ncbi_table not in known_tables()
    table_id = 1 if fallback else ncbi_table
    out = []
    for d in compare(inferred, table_id):
        c = classify_change(d["codon"], d["called"], lineage, organism)
        row = dict(d)
        row["class"] = c["class"]
        row["entries"] = c["entries"]
        row["table_used"] = table_id
        row["table_fallback"] = fallback
        out.append(row)
    return out


def most_interesting(cands):
    """The highest class among candidates (unreported_type > ... > known_taxon), or ''."""
    best = -1
    for c in cands:
        best = max(best, CLASS_ORDER.index(c["class"]))
    return CLASS_ORDER[best] if best >= 0 else ""


def stop_calls(inferred):
    """Calls at the three standard stop codons, e.g. {'TAA': '?', 'TAG': 'Q', 'TGA': '?'}."""
    idx = {c: i for i, c in enumerate(CODONS)}
    return {s: (inferred[idx[s]].upper() if len(inferred) == 64 else NONCALL) for s in ("TAA", "TAG", "TGA")}
