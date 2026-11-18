# eukcode-scan

An open, pre-registered screen for **alternative nuclear genetic codes across single-celled
eukaryotes**, using every protist genome assembly in NCBI and, from v0.3.0, every protist
transcriptome assembly in the INSDC Transcriptome Shotgun Assembly (TSA) archive.

Made by MV2 of Munim, Inc. (mv2@munimtech.com). The code was written by AI agents working in
his name. Licence: Apache-2.0.

## Why

Most life reads DNA with the standard genetic code. Some lineages have changed it: certain
ciliates read the stop codons UAA and UAG as glutamine, *Blastocrithidia* reads all three stop
codons as amino acids, and the CUG yeasts read CUG as serine or alanine. In bacteria and archaea
a systematic screen of 250,000 genomes found five new codes (Shulgina & Eddy 2021). In
eukaryotes, new codes are still mostly found one organism at a time. The pond ciliate that reads
UAA as lysine and UAG as glutamate (McGowan et al. 2023) turned up while a sequencing method was
being tested. Earlier screens covered parts of the tree: Swart et al. 2016 (the MMETSP
transcriptomes, stop codons only), Chen et al. 2023 (ciliates) and McGowan et al. 2024 (ciliate
MAGs). The P10K database (Gao et al. 2024) ran Codetta on about 1,760 protist datasets in its
annotation pipeline, but it reports codes for ciliates only. This project reports a code call
for every protist genome and transcriptome in the public INSDC archives, and traces where each
unexpected call's evidence comes from.

**What a success would be:** a genetic code nobody has reported before. The floor result, if
nothing new turns up, is an open map of where the known codes occur across every protist genome
and transcriptome in the public archives.

## Status

- **Stage 1 (pilot, registered 2026-11-17):** Codetta is run on nine genomes with documented
  codes, and an NCBI census lists the target assemblies. The predictions are in
  [PREREGISTRATION.md](PREREGISTRATION.md). The pilot job is queued; no results yet.
- **v0.2.0 (2026-11-17):** the stage 2 runner (`make_shards.py`, `run_screen.py`), low-memory
  genome sampling, and `eukcode/known.py`, a catalogue of documented eukaryotic nuclear
  reassignments that sorts every candidate into a class fixed before any data were seen
  (PREREGISTRATION.md, addendum A1).
- **v0.3.0 (2026-11-17): the transcriptome tranche** (`run_tsa.py`, `eukcode/tsa.py`,
  addendum A2). Many protist lineages, including uncultured single cells, have a transcriptome
  but no genome. The tranche:
  - lists every INSDC TSA set under the target groups through the ENA portal API, keeps one per
    species and fetches it from NCBI GenBank;
  - samples whole transcripts, keeping only the longest isoform of each gene;
  - sorts species into novelty tiers set by prior art. `mmetsp`: Swart et al. 2016 already
    screened all MMETSP transcriptomes for stop-codon changes. `pre2016`: first public before
    that screen. `new`: sets no systematic screen has seen.

  It adds a transcriptome control pilot (P5, P6) and a pre-registered biological prediction
  (P7) before any transcriptome is read. The census of TSA sets is queued.
- **v0.5.0 (2026-11-18)** includes the transcriptome control pilot that was planned as v0.4.0
  and was never released on its own. It adds:
  - **The transcriptome control pilot** (`run_tsa.py controls`, `eukcode/tsa_pilot.py`,
    addendum A3). It resolves 16 controls (8 known alternative codes, 8 standard-code species)
    to their TSA sets and runs Codetta on a 2 Mb sample of each. P5 and P6 are graded with
    counting rules fixed in A3 before any output existed. The job is queued.
  - **Check 2, where a candidate's evidence comes from** (`run_check2.py`, `eukcode/support.py`,
    `eukcode/ebi_blast.py`, addendum A4). A new code called from a contaminated sample would be
    a false discovery: prey inside a predator's cells, or bacteria or yeasts in a culture, can
    lend the sample their codon usage. Check 2:
    - reruns Codetta on the same sample with its alignment output kept;
    - maps every aligned codon back to its transcript or genome chunk, recomputing each codon
      to verify the mapping;
    - sends the 20 sequences carrying most of the candidate codon's evidence to blastx at
      EMBL-EBI (UniProtKB, the species' own taxid excluded);
    - grades the signal as pass, pass_unplaced, fail or inconclusive, with rules fixed in A4.

    Before check 2 is used, it is validated on a known code (*Euplotes*, UGA = Cys) and on a
    constructed 50/50 ciliate + alga mix where the truth is known (P8-P10). The validation job
    is queued.
  - **The P10K CUG question** (`eukcode/cug.py`, addendum A5). The P10K paper assigns NCBI table
    12 (CUG = Ser, the code of the CUG yeasts) to some spirotrich ciliates, without discussion.
    No ciliate has ever been shown to change a sense codon. If the call is real, it would be the
    first nuclear CUG change outside the yeasts. If not, the likely cause is yeast RNA in the
    cultures. Prediction P11 (no ciliate CTG = Ser call survives check 2) was registered before
    any data were seen. A5 also adds P10K to the prior art and narrows the novelty claim.
- **Stage 2:** the screen itself, with every candidate checked against a second sample,
  contamination, the literature and a second assembly. It starts only if the pilot's P1 and P2
  hold (and, for transcriptomes, P5 and P6).

## Use

```
pip install -r requirements.txt
bash install_codetta.sh          # clones and builds Codetta 2.0 + HMMER 3.3.2 (GPL-3.0, not bundled)
python -m unittest -v test_codes test_fasta test_fasta_stream test_codetta_io test_ncbi_grade test_known test_shards test_tsa test_tsa_pilot test_support test_cug
python run_pilot.py --data-dir data --out-dir out            # stage 1
python make_shards.py out/census_species.csv shards --per-shard 30
python run_screen.py --shard shards/shard_00.csv --data-dir data --out-dir out-s00 --cap-mb 2
python run_tsa.py census --out-dir out-tsa                   # transcriptome census (metadata only)
python run_tsa.py controls --data-dir data --out-dir out-tctl --cap-mb 2   # transcriptome control pilot
python run_tsa.py screen --census out-tsa/tsa_species.csv --shard-index 0 --data-dir data --out-dir out-t00
python run_check2.py demo --data-dir data --out-dir out-c2demo             # check 2 validation (P8-P10)
python run_check2.py records --records-dir out-t00/species --data-dir data --out-dir out-c2
```

The pilot downloads the Codetta Pfam 35 database (about 3 GB unpacked) and the control genomes
from NCBI Datasets. It writes `out/RESULTS.md`, `out/results.json`, Codetta's per-codon output
files, `census_species.csv` and `census_summary.json`. A screen shard writes one JSON record
per assembly (`species/`), `screen_rows.csv`, `candidates.csv` and `screen_summary.json`. The
TSA census writes `tsa_sets.csv`, `tsa_species.csv`, `tsa_skipped.csv`,
`tsa_controls_resolved.csv` and `tsa_census_summary.json`. The transcriptome control pilot
writes `TSA_PILOT.md`, `tsa_pilot.json`, `tsa_controls_resolved.csv` and the screen tables.
Check 2 writes `CHECK2.md`, `check2_summary.json` and one `check2/<case>.json` per case, with
every searched sequence, its best hit and its class.

## Candidate classes (eukcode/known.py)

| class | meaning |
|---|---|
| `known_taxon` | the change is documented for this species' lineage; NCBI just doesn't record it |
| `known_in_clade` | documented in the same major clade, not for this taxon |
| `known_elsewhere` | documented only in another major clade |
| `unreported_type` | not in the catalogue of documented eukaryotic nuclear codes |

The catalogue was compiled on 2026-11-17 and re-checked on 2026-11-18 against every nuclear
code listed in McGowan et al. 2024. It draws on:
- McGowan et al. 2024 (PLoS Genetics);
- Pánek et al. 2017 (BMC Biology);
- the 2025 *Blastocrithidia* comparative genomics paper (mBio);
- the NCBI translation tables.

It is a first pass. Nothing is called new without the pre-registered checks and a manual
literature review, including Lukeš et al. 2025 (*Natural and artificial variations of the
standard genetic code*, Current Biology 35:R1104).

## Credits

- Codetta: Y. Shulgina and S. R. Eddy, *A computational screen for alternative genetic codes in
  over 250,000 genomes*, eLife 2021, and *Codetta: predicting the genetic code from nucleotide
  sequence*, Bioinformatics 2023. https://github.com/kshulgina/codetta (GPL-3.0). It is called
  as an external program and not redistributed. Check 2 reads its alignment output format as
  written by `codetta.py` 2.0.
- Prior screens: E. C. Swart et al., *Genetic codes with no dedicated stop codon*, Cell 2016;
  S. M. Heaphy et al., Mol Biol Evol 2016; W. Chen et al., *Stop or not: genome-wide profiling
  of reassigned stop codons in ciliates*, Mol Biol Evol 2023; B. McGowan et al., PLoS Genet
  2024; X. Gao et al., *The P10K database*, Nucleic Acids Res 2024.
- HMMER (BSD-3-Clause), Pfam (CC0), NCBI Datasets and E-utilities, ENA portal API, INSDC TSA
  (public data), EMBL-EBI Job Dispatcher NCBI BLAST+ service and UniProtKB (CC BY 4.0).
