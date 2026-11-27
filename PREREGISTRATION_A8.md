# Addendum A8: a KACI fast triage path, with six control predictions (P16-P21)

Registered 2026-11-27, before any KACI output on any genome exists in this project.
Licence: Apache-2.0. By MV2 of Munim, Inc.; the code and this text were written by AI
agents working in his name, and no human expert has reviewed them yet.

## Why

Stage 2 of eukcode-scan runs Codetta (Shulgina & Eddy 2021, eLife 10:e71402) on every
protist genome and transcriptome. Codetta aligns Pfam profile HMMs to six-frame
translations, which is slow. KACI (K-mer Assisted Code Inference; A. V. Melnykov, PNAS
2026, 123(37):e2610659123, doi:10.1073/pnas.2610659123; code
https://github.com/artmeln/lib-kaci, MIT licence) infers the same per-codon decoding
probabilities from amino-acid k-mers of Pfam families and reports a >100-fold speed-up
at a very slight loss of accuracy against Codetta on more than 200,000 bacterial and
archaeal genomes. Its published screen covers bacteria and archaea only (searched
2026-11-27: bioRxiv 10.64898/2026.04.02.715157, the PNAS paper, the lib-kaci README); we
found no published use of KACI on eukaryote nuclear genomes.

A8 adds KACI as a **triage** step: if it passes the controls below, stage 2 may run KACI
on every target first and send only flagged assemblies to the registered Codetta + check-2
path. KACI alone never makes a claim in this project. Nothing in PREREGISTRATION.md or
addenda A1-A7 changes.

## Pinned tool and data

- lib-kaci commit `cc8601f2121ddb2546b69c5d494652c38294c02a` (2026-07-12), built with
  `make libkaci` (install_kaci.sh).
- Reference: `ref_k11_link20.tar.gz`, Zenodo doi:10.5281/zenodo.19318166 (CC BY 4.0),
  2,169,547,436 bytes, md5 `6149c67b184cfb053bc11d44bf9b1778`. Exactly one `.bin` file
  must be inside; it is loaded with `load_ref_bin`.
- Parameters as in the author's `examples/process_query.py`: k-mer size 11, max query
  count 2, input translation table = standard code. Threads do not change results.
- Genomes: NCBI Datasets `GENOME_FASTA` packages of the accessions below. Records whose
  FASTA header contains mitochondri, apicoplast, chloroplast, plastid or kinetoplast are
  removed first (organelles use other codes; e.g. yeast mitochondria read TGA as Trp).

## Calling rule (the author's, unchanged)

From lib-kaci `examples/2_evaluate_results.py` (`infer_aa`), implemented in
`eukcode/kaci_io.py`: a codon with `n_ <= 3` is uncalled; otherwise the most probable of
`ACDEFGHIKLMNPQRSTVWY?` (first on a tie) is called if its probability is strictly greater
than 0.9999; a confident `?` stays uncalled. As in the author's script, AGA/AGG called as
Lys are reported separately as known artefacts and not counted as unexpected.

## Predictions

Each control is graded against its own NCBI table. A control **passes** iff (a) every
listed codon is called as stated, (b) no sense codon is called as a different amino acid
than its table says (AGA/AGG->K excepted), and (c) no stop codon of its table is called as
an amino acid. Uncalled codons are allowed and reported (they measure sensitivity, not
graded).

| id | accession | organism | table | must be called |
|---|---|---|---|---|
| P16 | GCF_000146045.2 | *Saccharomyces cerevisiae* S288C | 1 | CTG = L |
| P17 | GCF_000002945.2 | *Schizosaccharomyces pombe* 972h- | 1 | none (specificity) |
| P18 | GCF_000002765.6 | *Plasmodium falciparum* 3D7 (80% A+T) | 1 | none (specificity) |
| P19 | GCA_001661245.1 | *Pachysolen tannophilus* NRRL Y-2460 | 26 | CTG = A |
| P20 | GCF_000182965.3 | *Candida albicans* SC5314 | 12 | CTG = S |
| P21 | GCF_000189635.1 | *Tetrahymena thermophila* SB210 | 6 | TAA = Q and TAG = Q |

**Gate:** KACI triage is adopted for stage 2 only if all six pass in one run. If any
fails, KACI is not used for triage and the failure is published as it is, with the
codons involved. If the reference cannot be loaded in the job's memory, or a genome
cannot be fetched, the affected predictions are "not graded" (neither pass nor fail) and
the run is repeated on a host with more memory; they are never dropped.

## What we expect, and why it is not a foregone conclusion

P20 is the hardest: in *C. albicans* most CUG codons sit at positions that are not
conserved, and CUG is read about 97% as Ser and 3% as Leu, so KACI may leave CTG uncalled
(a P20 fail by rule (a)). P18 tests false calls in an extremely A+T-rich genome. P21 tests
whether KACI calls stop-to-sense reassignments in a ciliate, the commonest kind of nuclear
code change in eukaryotes.

## Files

`eukcode/kaci_io.py` (calling rule, grading, organelle filter), `test_kaci.py`,
`kaci_controls.json` (accessions, pins), `install_kaci.sh`, `run_kaci.py` (writes
`a8_grade.json`, `A8_RESULT.md`, `status.json` and raw KACI JSON per genome).
