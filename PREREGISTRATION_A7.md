
## Addendum A7 (2026-11-20, v0.7.0): a tRNA test of P10K's Euplotes CUG label

Registered before any A7 output existed. No eukcode job had produced any output when this was
written: the pilot, transcriptome pilot, census, check 2 validation and Codetta CUG test
(A6.4) were all still queued. A7 adds a test that does not depend on any of them. P1-P11 and
the A6.4 reading rules are unchanged.

**A7.1 Why a tRNA test.** P10K (Gao et al. 2024, NAR 52:D747) labels *Euplotes* with codon
tables "10,12" (A6.3). Table 12 reads CUG as serine. Every known nuclear CUG -> Ser change (the
CTG-clade yeasts) is carried by a tRNA that reads CUG (anticodon 5'-CAG-3') but has a
serine-type body. A codon-usage method such as Codetta (A6.4) asks what amino acids sit at CUG
positions; a tRNA test asks whether the genome even has the adaptor that would put serine
there. The two are independent lines of evidence. A literature check on 2026-11-20 found no
report of a CUG change in *Euplotes*: the 2023 PNAS frameshifting study that produced the second
TN1 transcriptome (Gaydukova et al., PRJNA896607, GKLJ01) translated with table 10 and never
discusses CUG; NCBI Taxonomy gives table 10 for *E. focardii* (taxid 36767).

**A7.2 Method.** tRNAscan-SE 2.0.12 (Chan et al. 2021, NAR 49:9077) with Infernal 1.1.5,
eukaryotic mode, `--detail`. `--detail` reports, for every predicted tRNA, the best-scoring of
the 22 eukaryotic isotype-specific covariance models, which is scored without regard to the
anticodon, and marks disagreement with the anticodon as "IPD". For calibration, GtRNAdb
scores human tRNA-Leu-CAG-1-1 at 124.9 bits on the Leu model and 82.0 on the Ser model. The
job first scans that one tRNA as an install check and stops if it is not found. Each input's
contigs are renamed c1, c2, ...; the map back to INSDC names is kept for every contig that
carries a tRNA, with its length, GC and whether a *Euplotes* C4A4 telomere
(CCCCAAAACCCC or GGGGTTTTGGGG) lies within 100 bp of either end. A "CUG-reading tRNA" has
anticodon CAG, or TAG (U34 can read all CUN). Its body class is its isotype CM: Leu-type,
Ser-type or other. Pseudogenes (note "pseudo") are listed but never counted. Code:
`eukcode/trna.py` (pure, unit-tested in `test_trna.py`) and `run_trna.py`; inputs in
`trna_cases.json`.

**A7.3 Inputs and roles (accessions checked on 2026-11-20).**
- Positive controls, known nuclear CUG -> Ser: *Candida albicans* SC5314 GCA_000182965.3
  (14.3 Mb) and *Debaryomyces hansenii* CBS767 GCA_000006445.2 (12.2 Mb).
- Negative control, a ciliate that reads CUG as Leu: *Tetrahymena thermophila* SB210
  GCA_000189635.1 (103.0 Mb; a 2005 JBC study counted 2 tRNA-Leu-CAG genes).
- Main test: the *Euplotes focardii* strain TN1 macronuclear genome GCA_001880345.2
  (49.3 Mb, 31,248 contigs), the same strain as both TN1 transcriptomes.
- Other *Euplotes* genomes: *E. woodruffi* Iz01 GCA_027382605.1 (72.2 Mb; not the Iz01
  assembly NCBI flags as contaminated) and *E. aediculatus* QDEA GCA_030463445.1 (97.3 Mb).
- Secondary: the TN1 transcriptomes GKLJ01000000 and HBJK01000000. They are poly(A)-selected,
  so few tRNAs are expected; they are read only under A7.5.

**A7.4 Predictions and grading (fixed now).**
- **P12 (method can see a CUG-Ser tRNA).** Each positive control shows at least one
  non-pseudo CAG tRNA whose isotype CM is Ser. Pass only if both do. If P12 fails, no result
  without a Ser-type CAG tRNA is read as evidence against CUG = Ser, and A7 is reported as
  "method not validated for this question".
- **P13 (method does not invent one in a ciliate).** Every non-pseudo CAG tRNA of
  *T. thermophila* is Leu-type, with at least one present. No CAG tRNA = "uninformative".
- **P14 (the main test), outcome classes for E. focardii TN1, in this order:**
  1. `candidate`: a Ser-type CAG tRNA on a contig with a telomere at either end. Goes to
     checks (BLAST of the contig, presence in both transcriptomes, A6.4) and a person's review.
     Never called a discovery on this evidence.
  2. `unattributed-ser-type`: a Ser-type CAG tRNA on a contig without a telomere. Reported;
     the contig is checked for origin before anything is read into it.
  3. `cug-leu-supported`: at least one Leu-type CAG tRNA and no Ser-type one.
  4. `no-cag-trna-leu-u34`: no CAG tRNA, at least one Leu-type TAG tRNA, no Ser-type TAG.
  5. `no-cug-decoder-found` or `cag-other-isotype`: reported as found.
  **Expected: class 3 or 4.** If P12 passes and the TN1 outcome is class 3 or 4, the result
  is stated as: "P10K's table-12 label for *Euplotes* is not supported at the tRNA level: the
  TN1 genome has no serine-type tRNA that reads CUG, while the same method finds the known
  ones in *Candida* and *Debaryomyces*."
- **P15 (descriptive).** The same classes for *E. woodruffi* and *E. aediculatus*. No
  prediction is graded; a Ser-type CAG tRNA in either is reported and checked as in P14.

**A7.5 Contamination reading.** If a TN1 transcriptome has a Ser-type CAG tRNA and the TN1
genome has none, it is reported as consistent with a contaminant (for example a CTG-clade
yeast) in that culture, which together with the A6.4 contamination pattern would explain
P10K's label. That needs the A6.4 result before it is stated as an explanation.

**A7.6 Limits stated in advance.** A tRNA gene count is not a decoding measurement: a
Ser-type tRNA could be unexpressed, and leucine could be inserted at CUG by a tRNA the models
misclass. Draft macronuclear assemblies may miss tRNA genes, so absence is weaker than
presence; that is why P14 class 4 and the U34 route are reported separately. The telomere rule
attributes a contig to *Euplotes* only when the assembly kept its telomeres. Results are
produced by code written by AI agents working in MV2's name and are not checked by a person
before they are first reported.

**A7.7 Compute.** One job, `eukcode-trna-v1`: 2 CPU, 1.5 GB, inputs streamed through temporary
storage and deleted after each scan, nothing written to the shared data volume, at most 90
minutes per input. Inputs run in the order of `trna_cases.json` (controls and the main test
first), so a time-out late in the run still leaves P12-P14 gradable.
