# Pre-registration: a eukaryote-wide screen for alternative nuclear genetic codes

Registered 2026-11-17 by MV2 of Munim, Inc., before any pilot output was seen. The analysis
code was written by AI agents working in his name. Licence: Apache-2.0 (code), CC BY 4.0 (text).

## The question

How many lineages of single-celled eukaryotes (protists) read the genetic code differently from
the standard code, and are there reassignments nobody has reported?

Nuclear code changes in eukaryotes are known mainly from ciliates, a few trypanosomatids
(Blastocrithidia), diplomonads, some green algae and the CUG yeasts. Nearly all of them were
found one organism at a time, often by accident: the 2026 report of a pond ciliate with UAA=Lys
and UAG=Glu came from testing a sequencing method. Systematic screens exist only for bacteria
and archaea (Shulgina & Eddy 2021, eLife 10:e71402, which found five new codes in 250,000
genomes), budding yeasts, and green-algal and animal mitochondria. In 2021 the Codetta authors
named eukaryotes as the next target. As far as a literature search on 2026-11-17 found, no
eukaryote-wide nuclear screen has been published.

**What a success would be:** a reassignment never reported before, supported by a codon-level
signal that survives the checks below. That would be a new genetic code, the kind of result that
gets written up. An honest census of where the known codes occur, with every protist genome in
NCBI checked, is the floor result even if nothing new turns up.

## Method

- Tool: Codetta 2.0 (Shulgina & Eddy, Bioinformatics 2023; GPL-3.0, used unmodified as an
  external program) with its Pfam 35 "enone" profile database and its default settings
  (e-value 1e-10, decoding probability 0.9999, max fraction 0.01, mitochondrial, transposon,
  viral, Sec and Pyl domains excluded).
- Input: each NCBI genome assembly is cut into 100 kb chunks. Chunks that are shorter than 5 kb
  or more than 10% N are dropped, and a seeded shuffle (seed 0) picks chunks up to a cap in Mb.
- Expected code: the nuclear genetic code NCBI Taxonomy assigns to the species (GCId).

## Stage 1: pilot (this registration)

The pilot uses 11 runs on 9 genomes with documented codes: *Mycoplasmoides genitalium* (table
4), *Saccharomyces cerevisiae* (1), *Candida albicans* (12), *Pachysolen tannophilus* (26),
*Plasmodium falciparum* (1), *Tetrahymena thermophila* (6), *Oxytricha trifallax* (6), a
*Euplotes* species (10) and *Blastocrithidia nonstop* (31). Each is run at an 8 Mb cap. Two
power runs repeat *C. albicans* and *T. thermophila* at a 2 Mb cap.

Predictions, graded automatically by `eukcode/grade.py`:

- **P1 (sensitivity):** among the positive controls that resolve to an NCBI assembly, all but at
  most one call every diagnostic codon as expected. The diagnostic codons are Mgen UGA=W, Calb
  CUG=S, Ptan CUG=A, Tthe and Otri UAA=Q and UAG=Q, Euplotes UGA=C, and Bnon UAG=E and UGA=W.
- **P2 (specificity):** zero called codons contradict table 1 in the two standard-code controls,
  and at most one contradiction occurs over all 8 Mb controls. Codons with a documented dual
  meaning (Bnon UAA) are not counted.
- **P3 (power):** both 2 Mb runs still call their diagnostic codons.
- **P4 (cost):** median Codetta wall time is 4 minutes or less per Mb at 2 hmmscan threads.

How each outcome changes stage 2:
- P1 or P2 fails: the screen does not start. The failure is diagnosed first (e.g. stop codons
  in AT-rich genomes, contamination) and reported as is.
- P3 holds: the screen uses a 2 Mb cap (4 Mb for genomes with GC below 25% or above 70%).
- P3 fails: the screen uses a 4-8 Mb cap and fewer genomes per job.

## Stage 2: the screen (to be registered in detail after the pilot, before it runs)

Scope: one assembly per species (reference genome first, then the best assembly level, then
contig N50) in Sar, Discoba, Metamonada, Amoebozoa, Haptista, Cryptophyceae, Rhodophyta,
Glaucocystophyceae, Chlorophyta, Apusozoa, Ichthyosporea, Choanoflagellata and the smaller
orphan lineages. The pilot job lists them in `census_species.csv`.

A **candidate** is any codon where Codetta calls an amino acid that differs from the NCBI table
for that species. Every candidate must pass these checks before anything is claimed:

1. Repeat on a second, disjoint chunk sample (seed 1) and get the same call.
2. Contamination: the codon's supporting alignments must come mostly from contigs whose best
   hits are eukaryotic, i.e. not from bacterial or organellar contamination. The per-contig
   split of the signal is reported.
3. Literature check against the 2025 Current Biology catalogue of known reassignments
   (*Natural and artificial variations of the standard genetic code*) and anything newer.
4. A second assembly of the same species or genus, where one exists, gives the same call.

Each candidate ends in one of three classes: (a) a known code that NCBI's taxonomy does not yet
record, (b) an artefact, or (c) an unreported reassignment. Only (c) would be called new.
Nothing is written to authors or databases until a person has checked the evidence.

## Addendum A1 (2026-11-17, v0.2.0, still before any pilot or screen output was seen)

Added so that the sorting of candidates is fixed before any data exist. Nothing above changes.

**A1.1 Automatic first-pass class.** Each candidate is labelled by `eukcode/known.py`, using a
catalogue of documented eukaryotic nuclear reassignments compiled on 2026-11-17 from McGowan
et al. 2024 (PLoS Genet 20:e1011512) and its references, Pánek et al. 2017 (BMC Biol 15:8),
the 2025 *Blastocrithidia* comparative genomics paper (mBio e00885-25) and the NCBI tables:

- `known_taxon`: the same codon to amino acid change is documented for a taxon in the species'
  NCBI lineage, or for its genus. Expected to end in class (a).
- `known_in_clade`: documented in the same major clade (e.g. Ciliophora, Kinetoplastea,
  Metamonada, Chlorophyta, Rhizaria) but not for this taxon. This could be a new independent
  origin or a new taxon with a known code. Full checks 1-4 apply.
- `known_elsewhere`: documented only in another major clade. Full checks 1-4 apply.
- `unreported_type`: not documented in any eukaryotic nuclear code in the catalogue. Full
  checks 1-4, plus a manual search of the newest literature before anything is said.

This label is only a first pass. Class (c) still requires checks 1-4 and a manual literature
check, including the 2025 Current Biology review.

**A1.2 Data handling.** Genomes are downloaded to temporary storage and deleted after sampling.
Only the Pfam database is cached. Sampling (`fasta.sample_file`) reads each assembly twice
rather than holding it in memory, and gives exactly the in-memory sample. Assemblies under
1 Mb or over 3 Gb are listed as skipped, with the reason, and are not silently dropped.
Check 1 uses `--seed 1 --disjoint-from 0` at the cap the species used in the screen.

**A1.3 Reporting.** Every screened assembly is reported, with its inferred code, stop-codon
calls, candidates and errors (`screen_rows.csv`, `candidates.csv`), whatever the outcome.
