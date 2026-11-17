# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
import gzip
import os
import tempfile
import unittest

from eukcode import tsa


def ena(acc, tax, bp, date, desc="", env="false", dataclass="TSA", keywords="TSA, Transcriptome Shotgun Assembly"):
    return {"accession": acc, "scientific_name": "Sp %s" % tax, "tax_id": str(tax), "base_count": str(bp),
            "first_public": date, "environmental_sample": env, "dataclass": dataclass, "keywords": keywords,
            "description": desc, "study_accession": "PRJ%s" % tax}


class TestAccessions(unittest.TestCase):
    def test_set_code(self):
        self.assertEqual(tsa.set_code("GABU01000000"), ("GABU", "01"))
        self.assertEqual(tsa.set_code("GKLJ02000000.1"), ("GKLJ", "02"))
        with self.assertRaises(ValueError):
            tsa.set_code("GCA_000001405.1")
        with self.assertRaises(ValueError):
            tsa.set_code("")

    def test_fasta_urls(self):
        self.assertEqual(tsa.fasta_urls("HBRG01000000"), [
            "https://ftp.ncbi.nlm.nih.gov/genbank/tsa/H/tsa.HBRG.1.fsa_nt.gz",
            "https://ftp.ebi.ac.uk/pub/databases/ena/tsa/public/hbr/HBRG01.fasta.gz"])

    def test_search_url(self):
        u = tsa.ena_search_url(5878)
        self.assertIn("result=tsa_set", u)
        self.assertIn("tax_tree%285878%29", u)
        self.assertIn("limit=0", u)
        self.assertIn("format=tsv", u)


class TestTiers(unittest.TestCase):
    def test_parse_tsv(self):
        rows = tsa.parse_tsv("accession\tscientific_name\nGABU01000000\tChilodonella uncinata\n")
        self.assertEqual(rows, [{"accession": "GABU01000000", "scientific_name": "Chilodonella uncinata"}])
        self.assertEqual(tsa.parse_tsv(""), [])

    def test_flags_and_tier(self):
        m = ena("HBRG01000000", 1, 5_000_000, "2021-02-01", desc="Tiarina fusa, MMETSP0472 TPA", dataclass="TPA")
        self.assertEqual(tsa.flags(m), {"environmental": False, "tpa": True, "mmetsp": True})
        self.assertEqual(tsa.set_tier(m), "mmetsp")
        self.assertEqual(tsa.set_tier(ena("GABU01000000", 1, 5_000_000, "2015-08-19")), "pre2016")
        self.assertEqual(tsa.set_tier(ena("GKLJ01000000", 1, 5_000_000, "2024-06-01")), "new")
        self.assertEqual(tsa.set_tier(ena("GKLJ01000000", 1, 5_000_000, "2016-09-01")), "new")
        self.assertTrue(tsa.flags(ena("GAAA01000000", 1, 1, "2020-01-01", env="true"))["environmental"])

    def test_choose_sets(self):
        rows = [
            ena("GAAA01000000", 1, 5_000_000, "2015-01-01"),
            ena("HBAA01000000", 1, 9_000_000, "2021-01-01", desc="MMETSP0001"),
            ena("GBBB01000000", 2, 3_000_000, "2020-01-01"),
            ena("GCCC01000000", 2, 4_000_000, "2022-01-01"),
            ena("GDDD01000000", 3, 500_000, "2022-01-01"),
            ena("GEEE01000000", 4, 9_000_000, "2022-01-01", env="true"),
            ena("GCCC01000000", 2, 4_000_000, "2022-01-01"),
            ena("bad", 5, 9_000_000, "2022-01-01"),
        ]
        chosen, skipped = tsa.choose_sets(rows)
        self.assertEqual([c["accession"] for c in chosen], ["GAAA01000000", "GCCC01000000"])
        self.assertEqual([c["tier"] for c in chosen], ["mmetsp", "new"])
        self.assertEqual([c["n_sets"] for c in chosen], [2, 2])
        self.assertEqual(sorted(s["skip_reason"] for s in skipped),
                         ["bad accession", "duplicate accession", "environmental sample", "under 1000000 bp"])

    def test_species_row_and_shards(self):
        chosen, _ = tsa.choose_sets([ena("GAAA01000000", 5, 2_000_000, "2023-01-01")])
        r = tsa.species_row(chosen[0], {5: {"gencode": 6, "lineage": "a; b"}}, "Sar")
        self.assertEqual((r["ncbi_gencode"], r["lineage"], r["group"], r["tier"]), ("6", "a; b", "Sar", "new"))
        self.assertEqual(tsa.species_row(chosen[0], {}, "Sar")["ncbi_gencode"], "")
        sp = [{"accession": "A", "group": "Sar", "organism": "z", "tier": "mmetsp"},
              {"accession": "B", "group": "Sar", "organism": "y", "tier": "new"},
              {"accession": "C", "group": "Discoba", "organism": "x", "tier": "pre2016"},
              {"accession": "D", "group": "Discoba", "organism": "w", "tier": "new"}]
        plan = tsa.plan_shards(sp, per_shard=2)
        self.assertEqual([[r["accession"] for r in s] for s in plan], [["D", "B"], ["C", "A"]])
        with self.assertRaises(ValueError):
            tsa.plan_shards(sp, per_shard=0)


class TestTranscripts(unittest.TestCase):
    def test_gene_key(self):
        self.assertEqual(tsa.gene_key("GKLJ01000001.1 TSA: Euplotes focardii strain TN1 foca_c0_g1_i1, transcribed RNA sequence"),
                         "foca_c0_g1")
        self.assertEqual(tsa.gene_key("TRINITY_DN12_c0_g1_i3 len=500"), "TRINITY_DN12_c0_g1")
        self.assertEqual(tsa.gene_key("GAAA01000005.1 TSA: Foo bar comp12_c0_seq3 transcribed RNA"), "comp12_c0")
        self.assertIsNone(tsa.gene_key("GABU01000001.1 TSA: Chilodonella uncinata GG3RKAA01BRRRC mRNA sequence"))
        self.assertIsNone(tsa.gene_key("HBRG01000001.1 TPA_asm: Tiarina fusa, MMETSP0472-doi:10.5281-zenodo.249982-Transcript-0"))

    def _write(self, d):
        recs = [
            ("t0 x_c0_g1_i1", "ACGT" * 100),                # 400, gene g1 (shorter isoform)
            ("t1 x_c0_g1_i2", "ACGT" * 200),                # 800, gene g1 (kept)
            ("t2 x_c1_g1_i1", "ACGT" * 50),                 # 200, too short
            ("t3 none", "N" * 100 + "ACGT" * 100),          # 500, 20% N
            ("t4 x_c2_g1_i1", "ACGT" * 150),                # 600
            ("t5 nogene", "ACG" * 100),                     # 300, own gene
        ]
        path = os.path.join(d, "set.fsa_nt.gz")
        with gzip.open(path, "wt", encoding="ascii") as fh:
            for h, s in recs:
                fh.write(">%s\n" % h)
                for i in range(0, len(s), 70):
                    fh.write(s[i:i + 70] + "\n")
        return path

    def test_index_and_sample(self):
        with tempfile.TemporaryDirectory() as d:
            path = self._write(d)
            self.assertEqual(tsa.transcript_index(path), [("t1", 1, 800), ("t4", 4, 600), ("t5", 5, 300)])
            allrec, n = tsa.sample_transcripts(path, 10_000)
            self.assertEqual(n, 3)
            self.assertEqual([x[0] for x in allrec], ["t1", "t4", "t5"])
            self.assertEqual(len(allrec[0][1]), 800)
            a, _ = tsa.sample_transcripts(path, 1000, seed=0)
            b, _ = tsa.sample_transcripts(path, 1000, seed=0)
            self.assertEqual(a, b)
            self.assertGreaterEqual(sum(len(s) for _, s in a), 1000)
            c, _ = tsa.sample_transcripts(path, 1000, seed=1, disjoint_from=0)
            self.assertEqual(set(x[0] for x in a) & set(x[0] for x in c), set())


class TestRunner(unittest.TestCase):
    def test_import_without_codetta(self):
        import run_tsa
        self.assertIn("Sar", run_tsa.DEFAULT_GROUPS)
        self.assertEqual(len(run_tsa.DEFAULT_GROUPS), len(set(run_tsa.DEFAULT_GROUPS)))


if __name__ == "__main__":
    unittest.main()
