
## Addendum A6 (2026-11-19, v0.6.0): MMETSP tier fix, P10K targets recorded, and the Euplotes CUG test

Registered before any pilot, census, screen or check 2 output existed: every eukcode job was
still queued when this was written. A6.1 corrects a fact in A2.2. A6.2 fixes the rule that
implements it. The predictions P1-P11 are unchanged.

**A6.1 Correction to A2.2.** A2.2 says MMETSP-derived TSA sets can be recognised because they
"mention MMETSP (for example the 2021 third-party re-assemblies of MMETSP, prefix HB)". That is
wrong for the metadata the census reads. The 2021 EBI/ELIXIR re-assemblies of MMETSP are 406
TSA sets in two INSDC projects, PRJEB37117 and PRJEB37159. Both projects are titled "Assembly
and functional analysis of the Marine Microbial Eukaryotic Transcriptome Project data set". In
the ENA portal their descriptions read only "TSA: <species>" or "TPA: <species>". On
2026-11-19 the ENA portal counted 0 TSA sets with "MMETSP" in the description and 0 with it in
the keywords. Under the v0.3.0-v0.5.0 rule the `mmetsp` tier would have been empty. Every
re-assembled MMETSP set first made public in 2021 would have been tiered `new`, the tier the
project calls unscreened.

**A6.2 The fix.** A TSA set is MMETSP-derived if "MMETSP" appears in its description or
keywords, **or** its study accession is PRJEB37117 or PRJEB37159 (`tsa.is_mmetsp`,
`tsa.MMETSP_STUDIES`). The census job queued under v0.3.0 is not re-run. It writes every
set with its study accession, so species tiers are recomputed from its per-set rows
with `tsa.retier_species` before any screen shard is built or any tier is reported. Both the
old and the new tier are reported. MMETSP data deposited by other submitters under other
projects, without the word MMETSP, would still be missed. The census writes every set to `tsa_sets.csv`
and the unusable ones again, with their reason, to `tsa_skipped.csv`. Before re-tiering, each
set whose accession is in `tsa_skipped.csv` gets that `skip_reason`, so only usable sets count,
as in `choose_sets`.

**A6.3 P10K's targets (A5.4), recorded.** `p10k_table12.json` was filled on 2026-11-19 from
P10K Figure 3, saved at full size on 2026-11-18. Table 12 is labelled only on **Euplotes**
("10,12", with a red arrow for disagreement with NCBI). Table 4 is labelled on
*Condylostoma* and *Blepharisma*. P10K's gene pages give "Genetic Code: 12" for its *Euplotes
focardii* dataset (P10K-NCBI-002372, from MMETSP0205) and "Genetic Code: 10" for its *E.
octocarinatus* genome (P10K-MW-002960). Two INSDC TSA sets exist for *E. focardii*, both from
strain TN1:
- **HBJK01000000** (30.1 Mb), the 2021 EBI re-assembly of the MMETSP reads, probably the same
  sample as P10K's table-12 dataset;
- **GKLJ01000000** (21.3 Mb), an independent sequencing of the same strain, released in 2024.

A4.5's `euplotes_uga` resolves *E. focardii* to its largest usable set, which is HBJK01.
`check2_a6.json` pins that accession, so the input does not depend on the chooser, and adds
CTG to the codons traced. The validation rules P8-P10 are unchanged.

**A6.4 The Euplotes CUG test (exploratory, reported with P11).** `check2_cug.json` runs check 2
for CTG on three more 2 Mb samples: GKLJ01 (seed 0), HBJK01 (seed 1, disjoint from seed 0,
which is check 1 for that set) and *E. octocarinatus* GKLI01000000 (seed 0). It runs before P1,
P2 and P8-P10 are graded, so its results are labelled exploratory until they are. Reading
rules, fixed now:
- **Positive control:** each Euplotes sample should call TGA = C. A sample that does not is
  reported as low-power, and its CTG call is not interpreted.
- **Contamination pattern:** a CTG = S call in HBJK01 whose check 2 fails, with the top
  foreign clade Opisthokonta (fungi), while GKLJ01 calls CTG = L or '?', is read as a
  contaminant in the MMETSP culture. That would explain P10K's table-12 label.
- **Not seen in our data:** if CTG = L or '?' in every sample, P10K's call is reported as not
  reproduced by Codetta on 2 Mb samples. That is a disagreement, not a refutation, since P10K
  used whole transcriptomes.
- **Candidate:** a CTG = S call that passes check 2 in **both** independent TN1 sets (HBJK01
  and GKLJ01) is a candidate sense-codon reassignment. It goes through checks 3 and 4 and a
  person's review before anything is said. It is never called a discovery on this evidence
  alone.

**A6.5 Running the queued jobs in less memory (no change to any result).** None of the eukcode
jobs had started after about 44 hours in the queue. The shared machines mostly have room only
for jobs of about 1-1.5 GB, and these jobs asked for 2-4 GB. They are re-filed at 1.5 GB:
- **Genome pilot:** `run_pilot.py` now samples each genome with `fasta.sample_file` (two
  streaming passes). A1.2 states this gives exactly the in-memory sample: the same chunk names,
  lengths and seeded choice. Genome zips now go to temporary storage and are deleted after
  extraction instead of being kept in the shared data volume. It runs with one Codetta worker
  (2 hmmscan threads) and a 400-minute budget. Controls, caps, seeds and P1-P4 are unchanged.
- **Transcriptome control pilot:** it runs from the unchanged v0.5.0 code, so controls resolve
  exactly as A3 registered, with one worker and a 400-minute budget.
- **Check 2 validation:** it runs as `check2_a6.json` (A6.3).
The TSA census (1 GB) stays queued as it was.
