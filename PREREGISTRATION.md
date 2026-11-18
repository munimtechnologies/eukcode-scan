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

## Addendum A2 (2026-11-17, v0.3.0): the transcriptome tranche (stage 3)

Registered before any transcriptome sequence was read and before any pilot output was seen.
Nothing above changes.

**A2.1 Why.** Many protist lineages, especially uncultured cells isolated one at a time, have a
transcriptome but no genome assembly. Transcripts are mostly coding sequence, so a small sample
gives Codetta many more aligned codons per megabase than a genome does.

**A2.2 Prior art, which fixes the novelty tiers.** Swart et al. 2016 (Cell 166:691-702)
screened every transcriptome of the Marine Microbial Eukaryote Transcriptome Sequencing Project
(MMETSP) for stop-codon reassignments and found no new code outside ciliates in 265 non-ciliate
species. Heaphy et al. 2016 (Mol Biol Evol 33:2885) examined the MMETSP ciliates. Each species
therefore gets a tier:
- `mmetsp`: any of its TSA sets mentions MMETSP (for example the 2021 third-party re-assemblies
  of MMETSP, prefix HB).
- `pre2016`: otherwise, if any of its sets was first public before 2016-09-01.
- `new`: otherwise. These are sets no systematic screen could have seen.

Every tier is screened and reported, with `new` first. A candidate in tier `mmetsp` that is a
stop-codon change is compared with Swart et al.'s Data S1 before anything else is said.

**A2.3 Scope.** Every INSDC Transcriptome Shotgun Assembly (TSA) set listed by the ENA portal
API (`result=tsa_set`, `tax_tree(taxid)`) under the stage 2 groups plus Filasterea, Aphelida and
Breviatea. Environmental (metatranscriptome) sets and sets under 1 Mb are skipped with a reason.
For each species, a set not derived from MMETSP is preferred, then the largest. NCBI Taxonomy
gives the expected code. Sequences are downloaded from NCBI GenBank
(`genbank/tsa/<letter>/tsa.<CODE>.1.fsa_nt.gz`), or from the ENA mirror if that fails. The
census (`run_tsa.py census`) only reads metadata and may run before the pilot.

**A2.4 Sampling.** Whole transcripts of at least 300 nt with at most 10% N are used. Only the
longest isoform of each gene is kept (Trinity `..._cN_gN_iN` and `compN_cN_seqN` names). Other
headers count as one gene each. A seeded shuffle (seed 0) then picks transcripts up to a
**fixed 2 Mb cap**. Check 1 uses `--seed 1 --disjoint-from 0`, a sample that shares no gene
with the first one.

**A2.5 Transcriptome pilot and gate.** Before the TSA screen, the controls in
`tsa_controls.json` are run at the 2 Mb cap. Each control is resolved by name to its largest
usable TSA set and is never substituted.
- **P5 (sensitivity):** at least 3 positive controls resolve, and among those, all but at most
  one call every diagnostic codon as expected. The controls are *Euplotes focardii* and *E.
  raikovi* (UGA=C), *Tetrahymena thermophila*, *Paramecium tetraurelia*, *P. caudatum*,
  *Oxytricha trifallax* and *Stylonychia lemnae* (UAA=Q, UAG=Q), and *Candida albicans*
  (CUG=S).
- **P6 (specificity):** at least 3 standard-code controls resolve. Across all of them, at most
  one called codon contradicts table 1. The controls are *Chlamydomonas reinhardtii*,
  *Dictyostelium discoideum*, *Phaeodactylum tricornutum*, *Thalassiosira pseudonana*,
  *Plasmodium falciparum*, *Toxoplasma gondii*, *Gephyrocapsa (Emiliania) huxleyi* and
  *Bigelowiella natans*.
- The TSA screen starts only if the genome pilot's P1 and P2 hold **and** P5 and P6 hold. If
  fewer than 3 controls of a kind resolve, that prediction cannot be tested, the TSA screen
  does not start, and this is reported.

**A2.6 A prediction about the biology.** Swart et al. found no new code in 265 non-ciliate
MMETSP species. **P7:** in tier `new`, at least 80% of the non-ciliate species that finish
without error have no candidate at all. If P7 fails, the excess candidates are reported as a
result, whatever turns out to be their cause (artefact or biology).

**A2.7 Transcriptome-specific checks.** A transcriptome can contain RNA from prey, food
organisms or symbionts (for example the retained cryptophyte nuclei of *Mesodinium*). Check 2
therefore also reports the share of the candidate codon's support coming from transcripts
whose best hit is outside the species' major clade. Check 4 may use a genome assembly of the
same species or genus from stage 2, or a second TSA set. Every set is reported whatever the
outcome (`screen_rows.csv`, `candidates.csv`, one JSON record per set).

## Addendum A3 (2026-11-17, v0.4.0): running and grading the transcriptome control pilot

Registered before any genome or transcriptome output existed: neither the genome pilot nor
the TSA census had run when this was written. Nothing above changes; A3 fixes details A2 left
open.

**A3.1 Order.** The transcriptome control pilot (`run_tsa.py controls`) resolves the controls
itself and does not need the TSA census, so it may run before or alongside the genome pilot.
The gate in A2.5 is unchanged: the TSA screen needs P1, P2, P5 and P6 to hold. The pilot
reports the gate as closed while P1 and P2 are not yet graded.

**A3.2 Counting rules** (`eukcode/tsa_pilot.py`):
- A control *resolves* when an INSDC TSA set of at least 1 Mb, not environmental, is found for
  its name. The largest such set is used. Unresolved controls are listed and not replaced.
- **P5:** a resolved positive control that ends without a Codetta result (download error,
  timeout, time budget) counts as a **failure**. P5 holds if at least 3 positive controls
  resolve and the number that call every diagnostic codon is at least (resolved - 1).
- **P6:** only standard controls that finish are counted, and at least 3 must finish (an
  unfinished control has no calls, so it cannot pass or fail). P6 holds if the total number of
  called codons that contradict table 1 across them is at most 1. A codon with no call ('?')
  never counts as a contradiction.
- Each control is judged against the table registered in `tsa_controls.json`. If NCBI Taxonomy
  gives a different code for the resolved organism, that is reported in a note.
- **P7:** a species is a ciliate if "Ciliophora" is a rank in its NCBI lineage. Only tier `new`
  species that finish without error are counted.

**A3.3 Cost, reported but not predicted.** The pilot reports the median Codetta minutes per Mb
for transcriptome samples, to size the screen shards. There is no pass mark.

**A3.4 Prior art added.** Chen et al. 2023 (*Stop or not: genome-wide profiling of reassigned
stop codons in ciliates*, Mol Biol Evol 40:msad064) profiled stop-codon reassignments across
published ciliate genomes and transcriptomes. It covers ciliates only, so the eukaryote-wide
question stands. Any ciliate candidate is compared with that paper, as well as with Swart et
al. 2016 and McGowan et al. 2024, before anything is said.

## Addendum A4 (2026-11-17, v0.5.0): how check 2 is measured and graded

Registered before any check 2 output, any screen output and any pilot output existed. Nothing
above changes; A4 makes check 2 (and A2.7) exact.

**A4.1 Tracing the evidence.** For one candidate, the species' sample is drawn again with the
same seed and cap (so it is the same sample) and Codetta is rerun with its alignment output
kept. Every line of that output (one aligned codon: piece, frame, protein position, Pfam
profile) is mapped back to an input sequence and position by rebuilding Codetta 2.0's own
piece splitting (records over 100,000 nt are halved until no piece is longer) and its six
frames. As a self-check, the codon at the mapped position is recomputed and must equal the
codon Codetta wrote; any mismatch is reported. The candidate codon's observations are counted
per input sequence (transcript or genome chunk).

**A4.2 Homology search.** The 20 input sequences with the most observations of the candidate
codon (ties: earlier in the sample) are searched with blastx (standard code) at EMBL-EBI
against UniProtKB, E-value at most 1e-5, with the species' own NCBI taxid excluded so that a
contaminated assembly cannot vouch for itself. The query is the stretch covering that
sequence's aligned candidate codons plus 300 nt each side, at most 3,000 nt. The best hit
(first in EBI's order) gives a taxid, and NCBI Taxonomy gives its lineage.

**A4.3 Classes.** A best hit is `same_clade` if it is eukaryotic and its first two lineage
names below Eukaryota (for example Sar / Alveolata, Discoba / Euglenozoa) match the species'
(when the species has only one, that one must match); `other_eukaryote` if eukaryotic
otherwise; `bacteria`, `archaea` or `virus` by its top lineage name; `no_hit` if blastx finds
nothing; `unknown` if the hit has no usable lineage; `error` if the search failed. Only the
first five classes count as *classified*.

**A4.4 Grading** (`eukcode/support.py`, `grade_check2`), with shares weighted by the number of
observations of the candidate codon:
- **inconclusive** if fewer than 3 sequences are classified, or the classified sequences carry
  less than 20% of all observations of the codon;
- otherwise **fail** if the eukaryotic share is 0.5 or less;
- otherwise **pass** if the same-clade share is above 0.5;
- otherwise **fail** if a single other eukaryotic clade carries more than 0.5 (a prey or
  symbiont signal);
- otherwise **pass_unplaced**: eukaryotic, but the hits do not confirm the species' clade, as
  expected for orphan lineages with no close relatives in UniProtKB. This is reported as such
  and is weaker evidence than a pass.
Each result also reports whether one sequence carries more than half of all observations
(`concentrated`), the outside-clade share A2.7 asks for, and the prokaryote or virus share.

**A4.5 Validation, before check 2 is used on any candidate** (`run_check2.py demo`,
`check2_demo.json`). Each input is resolved by name, in the listed order, to its largest usable
INSDC TSA set; the first name that resolves is used and reported.
- `euplotes_uga`: a 2 Mb sample of a *Euplotes* transcriptome (UGA = Cys); codon TGA.
- `mix_uar`: 1 Mb of a UAR = Gln ciliate (*Tetrahymena thermophila* first) plus 1 Mb of a
  standard-code alga or protist (*Phaeodactylum tricornutum* first), one Codetta run on the
  mixed sample; codon TAA. This is a constructed contamination case where the truth is known.

Predictions:
- **P8 (mapping):** in every case that finishes, every alignment line maps back to the codon
  Codetta reported (0 mismatches, more than 0 lines).
- **P9 (attribution):** in `mix_uar`, at least 90% of the TAA observations come from ciliate
  transcripts, and check 2 on TAA (species = the ciliate) passes.
- **P10 (a real reassignment passes):** check 2 on the *Euplotes* TGA signal passes.

If P8 fails, check 2 is not used until the mapping is fixed. If P9 or P10 fails, the failure is
diagnosed and reported, and the rules in A4.2-A4.4 are revised in a new addendum before check
2 is applied to any candidate.

## Addendum A5 (2026-11-18, v0.5.0): P10K prior art, one correction, and prediction P11

Registered before any pilot, census, screen or check 2 output existed (every eukcode job was
still queued). The rules above are unchanged. A5.1 narrows the novelty claim made in "The
question", and A5.2 corrects one date there.

**A5.1 Prior art: the P10K database.** Gao et al. 2024 (*The P10K database: a data portal for
the protist 10 000 genomes project*, Nucleic Acids Res 52:D747, doi:10.1093/nar/gkad992) ran
Codetta and FACIL inside the gene-annotation pipeline of the Protist 10,000 Genomes Project. The
paper says the pipeline was used on the project's 1,101 newly sequenced datasets (mostly
ciliates) and on 659 re-annotated MMETSP transcriptomes. It reports codon tables **for ciliates
only**: tables 1, 4, 6, 10, 12, 27 and 30, no new code type, and disagreements with NCBI
Taxonomy for 31 ciliate species, mostly Heterotrichea. Per-dataset codes appear only on the
portal's gene pages. There is no bulk table and no results for non-ciliates. So eukcode-scan is
not the first project to run Codetta on many protist datasets. The claim is narrowed to this:
it is the first screen to **report** per-assembly code calls across all eukaryotic lineages in
the INSDC archives, with every candidate traced for contamination (check 2). Before anything is
said about them, ciliate candidates are compared with P10K's genus-level codes (its Figure 3) as
well as with the papers in A3.4. Candidates from MMETSP-derived data are compared with the code
P10K gives that dataset.

**A5.2 Correction.** "The question" calls the pond ciliate with UAA = Lys and UAG = Glu a 2026
report. The paper is McGowan et al. 2023 (PLoS Genet 19:e1010913), already in the A1
catalogue. The 2026 items were press coverage of it.

**A5.3 The P10K CUG question (prediction P11).** P10K assigns NCBI table 12 to some Spirotrichea.
Table 12 is the alternative yeast code: CUG = Ser, with UAA, UAG and UGA as stops. Every ciliate
reassignment reported so far is at a stop codon, and no nuclear CUG change is known outside the
budding yeasts. A table 12 call in a ciliate therefore means one of two things:
- **an artefact.** One likely source is RNA or DNA from a CUG-Ser yeast (Serinales) in the
  culture. Spirotrich genomes are AT-rich and use CUG rarely, so even a small contaminant could
  supply most of the CUG evidence. Another is a pipeline that maps the inferred code to the
  nearest NCBI table.
- **a real CUG reassignment.** That would be the first sense-codon change in any ciliate and
  the first nuclear CUG change outside the yeasts.

**P11:** no ciliate in the eukcode screen (genomes and transcriptomes) has a Codetta CTG = S
call that survives check 2 (outcome pass or pass_unplaced). A CTG = S call whose check 2 is
fail or inconclusive counts as explained or unresolved, not as a reassignment. P11 is graded
only after every CTG = S ciliate has a check 2 result for CTG (`eukcode/cug.py`, `grade_p11`).
Our prior expectation is that P11 holds. If it does not, the surviving call goes through
checks 1, 3 and 4 like any candidate. It is reported as a candidate, never as a discovery,
until a person has checked it.

**A5.4 P10K's own targets.** The genera that P10K's Figure 3 labels with table 12 are recorded
in `p10k_table12.json`, together with the date and how they were read. This is done before any
eukcode output for those genera is seen, and the list is not edited afterwards. Every dataset
the screen has for them is reported with its CTG call (`cug.p10k_targets`):
- `agrees_with_p10k` (S), which goes to check 2;
- `contradicts_p10k` (L);
- `no_call` ('?': too few CUG observations to call);
- `mixed`;
- `not_screened`.

A dataset that contradicts P10K is reported as a disagreement with P10K and not as an error by
either side. Codetta at a 2 Mb cap may simply have less evidence than a whole-genome run.
