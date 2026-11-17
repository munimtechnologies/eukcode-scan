# eukcode-scan

An open, pre-registered screen for **alternative nuclear genetic codes across single-celled
eukaryotes**, using every protist genome assembly in NCBI.

Made by MV2 of Munim, Inc. (mv2@munimtech.com). The code was written by AI agents working in
his name. Licence: Apache-2.0.

## Why

Most life reads DNA with the standard genetic code. Some lineages have changed it: certain
ciliates read the stop codons UAA and UAG as glutamine, *Blastocrithidia* reads all three stop
codons as amino acids, and the CUG yeasts read CUG as serine or alanine. In bacteria and archaea
a systematic screen of 250,000 genomes found five new codes (Shulgina & Eddy 2021). Eukaryotes
have never been screened that way, and new eukaryotic codes are still found by accident (a pond
ciliate in 2026, two ciliate variants in a November 2026 preprint). This project runs the
screen.

## Status

- **Stage 1 (pilot, registered 2026-11-17):** Codetta is run on nine genomes with documented
  codes, and an NCBI census lists the target assemblies. The predictions are in
  [PREREGISTRATION.md](PREREGISTRATION.md). The pilot job is queued; no results yet.
- **v0.2.0 (2026-11-17):** the stage 2 runner (`make_shards.py`, `run_screen.py`), low-memory
  genome sampling, and `eukcode/known.py`, a catalogue of documented eukaryotic nuclear
  reassignments that sorts every candidate into a class fixed before any data were seen
  (PREREGISTRATION.md, addendum A1).
- **Stage 2:** the screen itself, with every candidate checked against a second sample,
  contamination, the literature and a second assembly. It starts only if the pilot's P1 and P2
  hold.

## Use

```
pip install -r requirements.txt
bash install_codetta.sh          # clones and builds Codetta 2.0 + HMMER 3.3.2 (GPL-3.0, not bundled)
python -m unittest -v test_codes test_fasta test_fasta_stream test_codetta_io test_ncbi_grade test_known test_shards
python run_pilot.py --data-dir data --out-dir out            # stage 1
python make_shards.py out/census_species.csv shards --per-shard 30
python run_screen.py --shard shards/shard_00.csv --data-dir data --out-dir out-s00 --cap-mb 2
```

The pilot downloads the Codetta Pfam 35 database (about 3 GB unpacked) and the control genomes
from NCBI Datasets. It writes `out/RESULTS.md`, `out/results.json`, Codetta's per-codon output
files, `census_species.csv` and `census_summary.json`. A screen shard writes one JSON record
per assembly (`species/`), `screen_rows.csv`, `candidates.csv` and `screen_summary.json`.

## Candidate classes (eukcode/known.py)

| class | meaning |
|---|---|
| `known_taxon` | the change is documented for this species' lineage; NCBI just doesn't record it |
| `known_in_clade` | documented in the same major clade, not for this taxon |
| `known_elsewhere` | documented only in another major clade |
| `unreported_type` | not in the catalogue of documented eukaryotic nuclear codes |

The catalogue (compiled 2026-11-17) draws on McGowan et al. 2024 (PLoS Genetics), Pánek et
al. 2017 (BMC Biology), the 2025 *Blastocrithidia* comparative genomics paper (mBio) and the
NCBI translation tables. It is a first pass: nothing is called new without the pre-registered
checks and a manual literature review.

## Credits

- Codetta: Y. Shulgina and S. R. Eddy, *A computational screen for alternative genetic codes in
  over 250,000 genomes*, eLife 2021, and *Codetta: predicting the genetic code from nucleotide
  sequence*, Bioinformatics 2023. https://github.com/kshulgina/codetta (GPL-3.0). It is called
  as an external program and not redistributed.
- HMMER (BSD-3-Clause), Pfam (CC0), NCBI Datasets and E-utilities (public data).
