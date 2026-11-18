# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
import gzip
import json
import os
import tempfile
import unittest

from eukcode import ebi_blast, support

CIL = "cellular organisms; Eukaryota; Sar; Alveolata; Ciliophora; Spirotrichea"


def line(piece, frame, qpos, codon, profile="P1", evalue="1e-20", hmm=0, ci=0):
    return "%d,%d.%d%06d,%s,%d,%d,%d,%s,%s,%d\n" % (ci, piece, frame, qpos, codon, piece, frame, qpos, profile, evalue, hmm)


class TestPieces(unittest.TestCase):
    def test_codetta_split(self):
        self.assertEqual(support.codetta_pieces([10, 0, 250001]),
                         [(0, 0, 10), (2, 0, 62500), (2, 125000, 62500), (2, 62500, 62500), (2, 187500, 62501)])
        self.assertEqual(support.codetta_pieces([10], piece_max=4), [(0, 0, 2), (0, 5, 2), (0, 2, 3), (0, 7, 3)])
        self.assertEqual(support.piece_sequences(["ABCDEFGHIJ"], support.codetta_pieces([10], piece_max=4)),
                         ["AB", "FG", "CDE", "HIJ"])

    def test_codon_at(self):
        s = "ATGGCCTAAGT"
        self.assertEqual(support.codon_at(s, 0, 0), ("ATG", (0, 3)))
        self.assertEqual(support.codon_at(s, 0, 2), ("TAA", (6, 9)))
        self.assertEqual(support.codon_at(s, 1, 0), ("TGG", (1, 4)))
        self.assertEqual(support.codon_at(s, 2, 1), ("CTA", (5, 8)))
        self.assertEqual(support.codon_at(s, 3, 0), ("ACT", (8, 11)))
        self.assertEqual(support.codon_at(s, 4, 0), ("CTT", (7, 10)))
        self.assertEqual(support.codon_at(s, 5, 1), ("GGC", (3, 6)))
        self.assertEqual(support.codon_at(s, 0, 3), ("", (9, 12)))
        self.assertEqual(support.revcomp("ATGGCCTAAGT"), "ACTTAGGCCAT")


class TestAlignment(unittest.TestCase):
    def test_parse(self):
        r = support.parse_alignment_line("48,3.2000045,TAA,3,2,45,PF00001,1.5e-20,12\n")
        self.assertEqual(r, {"codon_index": 48, "codon": "TAA", "piece": 3, "frame": 2, "qpos": 45,
                             "profile": "PF00001", "evalue": 1.5e-20, "hmm_pos": 12})
        self.assertEqual(support.parse_alignment_line("48,3.2000045,TAA,3,2,45,Foo,bar,1e-11,7")["profile"], "Foo,bar")
        self.assertIsNone(support.parse_alignment_line("x,y"))
        self.assertIsNone(support.parse_alignment_line("a,b,c,d,e,f,g,h,i"))

    def test_files(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.alignment_output.txt.gz")
            with gzip.open(p, "wt", encoding="ascii") as fh:
                fh.write(line(0, 0, 2, "TAA"))
                fh.write("garbage\n")
                fh.write(line(1, 3, 0, "TTA"))
            self.assertEqual([(r["piece"], r["codon"]) for r in support.iter_alignment(p)], [(0, "TAA"), (1, "TTA")])
            q = os.path.join(d, "x.sequence_pieces.fna")
            with open(q, "w", encoding="ascii") as fh:
                fh.write(">piece_0\nACGT\n>piece_1\nGG\n")
            self.assertEqual(support.read_pieces_file(q), {0: "ACGT", 1: "GG"})

    def test_tally(self):
        seqs = ["ATGGCCTAAGT", "TAATAA"]
        pieces = support.codetta_pieces([len(s) for s in seqs])
        self.assertEqual(pieces, [(0, 0, 11), (1, 0, 6)])
        rows = [support.parse_alignment_line(x) for x in (
            line(0, 0, 2, "TAA"), line(1, 0, 0, "TAA"), line(1, 0, 1, "TAA", profile="P2"),
            line(0, 1, 0, "TGG"), line(0, 0, 0, "TGA"), line(5, 0, 0, "TAA"), line(1, 3, 0, "TTA", profile="P3"))]
        t = support.tally(rows, pieces, support.piece_sequences(seqs, pieces), ["TAA", "TTA"])
        self.assertEqual((t["rows"], t["mismatches"], len(t["mismatch_examples"])), (7, 2, 2))
        self.assertEqual(t["by_codon"]["TAA"]["total"], 3)
        self.assertEqual(t["by_codon"]["TAA"]["per_record"],
                         {0: {"obs": 1, "profiles": 1, "start": 6, "end": 9},
                          1: {"obs": 2, "profiles": 2, "start": 0, "end": 6}})
        self.assertEqual(t["by_codon"]["TTA"]["per_record"], {1: {"obs": 1, "profiles": 1, "start": 3, "end": 6}})

    def test_top_and_region(self):
        self.assertEqual(support.top_records({0: {"obs": 1}, 1: {"obs": 2}, 2: {"obs": 2}}, 2), [1, 2])
        self.assertEqual(support.region(1000, 1100, 5000), (700, 1400))
        self.assertEqual(support.region(0, 10, 50), (0, 50))
        self.assertEqual(support.region(0, 10000, 20000), (3500, 6500))


class TestClassify(unittest.TestCase):
    def test_classes(self):
        self.assertEqual(support.major_clade(CIL), ("Sar", "Alveolata"))
        self.assertIsNone(support.major_clade("cellular organisms; Bacteria; Pseudomonadota"))
        c = support.classify_hit
        self.assertEqual(c("cellular organisms; Eukaryota; Sar; Alveolata; Ciliophora; Oligohymenophorea", CIL), "same_clade")
        self.assertEqual(c("cellular organisms; Eukaryota; Sar; Stramenopiles; Ochrophyta", CIL), "other_eukaryote")
        self.assertEqual(c("cellular organisms; Bacteria; Pseudomonadota", CIL), "bacteria")
        self.assertEqual(c("Viruses; Riboviria", CIL), "virus")
        self.assertEqual(c(None, CIL), "no_hit")
        self.assertEqual(c("", CIL), "unknown")
        self.assertEqual(c("other entries; plasmids", CIL), "unknown")
        self.assertEqual(c("Eukaryota; Hemimastigophora; Spironematellidae", "Eukaryota; Hemimastigophora"), "same_clade")


SAME = ("same_clade", ("Sar", "Alveolata"))


class TestGrade(unittest.TestCase):
    def test_pass(self):
        g = support.grade_check2({0: 10, 1: 5, 2: 5, 3: 5}, {0: SAME, 1: SAME, 2: ("bacteria", None), 3: ("no_hit", None)}, 50)
        self.assertEqual((g["outcome"], g["classified"], g["coverage"], g["eukaryotic_share"], g["same_clade_share"],
                          g["prokaryote_or_virus_share"], g["outside_clade_share"], g["top1_share_of_total"], g["concentrated"]),
                         ("pass", 3, 0.4, 0.75, 0.75, 0.25, 0.25, 0.2, False))

    def test_fail_and_unplaced(self):
        stra = ("other_eukaryote", ("Sar", "Stramenopiles"))
        g = support.grade_check2({0: 4, 1: 4, 2: 4}, {0: SAME, 1: stra, 2: stra}, 12)
        self.assertEqual((g["outcome"], g["same_clade_share"], g["top_foreign_clade"], g["top_foreign_share"]),
                         ("fail", 0.3333, "Sar / Stramenopiles", 0.6667))
        g = support.grade_check2({0: 4, 1: 4, 2: 4}, {0: SAME, 1: ("other_eukaryote", ("Discoba", "Euglenozoa")),
                                                      2: ("other_eukaryote", ("Amoebozoa", "Evosea"))}, 12)
        self.assertEqual((g["outcome"], g["top_foreign_clade"], g["top_foreign_share"]), ("pass_unplaced", "Discoba / Euglenozoa", 0.3333))
        g = support.grade_check2({0: 3, 1: 3, 2: 3}, {i: ("bacteria", None) for i in range(3)}, 9)
        self.assertEqual((g["outcome"], g["eukaryotic_share"]), ("fail", 0.0))

    def test_inconclusive_and_concentrated(self):
        self.assertEqual(support.grade_check2({0: 3, 1: 3}, {0: SAME, 1: SAME}, 6)["outcome"], "inconclusive")
        self.assertEqual(support.grade_check2({0: 1, 1: 1, 2: 1}, {0: SAME, 1: SAME, 2: SAME}, 100)["outcome"], "inconclusive")
        g = support.grade_check2({0: 60, 1: 1, 2: 1}, {0: SAME, 1: SAME, 2: SAME}, 100)
        self.assertEqual((g["outcome"], g["concentrated"], g["coverage"]), ("pass", True, 0.62))

    def test_origin(self):
        per = {0: {"obs": 3}, 1: {"obs": 1}, 5: {"obs": 4}}
        self.assertEqual(support.origin_shares(per, lambda r: "A" if r < 5 else "B"), {"A": 0.5, "B": 0.5})
        self.assertEqual(support.origin_shares({}, lambda r: "A"), {})


class TestBlast(unittest.TestCase):
    def test_best_hit(self):
        res = {"hits": [
            {"hit_num": 2, "hit_acc": "X2", "hit_uni_ox": "9606", "hit_hsps": []},
            {"hit_num": 1, "hit_acc": "P1", "hit_id": "ID1", "hit_uni_os": "Euplotes focardii", "hit_uni_ox": "36764",
             "hit_uni_de": "Actin", "hit_hsps": [
                 {"hsp_bit_score": 50.0, "hsp_expect": 1e-8, "hsp_identity": 40.0, "hsp_align_len": 100},
                 {"hsp_bit_score": 80.0, "hsp_expect": 1e-15, "hsp_identity": 60.0, "hsp_align_len": 120}]}]}
        h = ebi_blast.best_hit(res)
        self.assertEqual((h["acc"], h["taxid"], h["organism"], h["desc"], h["evalue"], h["bits"], h["identity"], h["align_len"]),
                         ("P1", 36764, "Euplotes focardii", "Actin", 1e-15, 80.0, 60.0, 120))
        h = ebi_blast.best_hit({"hits": [{"hit_num": 1, "hit_acc": "Q1", "hit_os": "Foo bar",
                                          "hit_def": "SP:Q1 X OS=Foo bar OX=123 GN=x", "hit_hsps": []}]})
        self.assertEqual((h["taxid"], h["organism"], h["evalue"]), (123, "Foo bar", None))
        self.assertIsNone(ebi_blast.best_hit({"hits": []}))
        self.assertIsNone(ebi_blast.best_hit({}))

    def test_params(self):
        p = ebi_blast.run_params("ACGT", "t", exclude_taxids=["5911", "", None])
        self.assertEqual((p["program"], p["stype"], p["database"], p["exp"], p["alignments"], p["negative_taxids"]),
                         ("blastx", "dna", "uniprotkb", "1e-5", "5", "5911"))
        self.assertNotIn("negative_taxids", ebi_blast.run_params("ACGT"))
        self.assertEqual(ebi_blast.run_params("ACGT")["email"], "mv2@munimtech.com")


class TestRunner(unittest.TestCase):
    def test_records_and_parser(self):
        import run_check2
        a = run_check2.build_parser().parse_args(["demo"])
        self.assertEqual((a.mode, a.top_n, a.blast_workers), ("demo", 20, 4))
        with tempfile.TemporaryDirectory() as d:
            recs = {
                "a.json": {"accession": "GAAA01000000", "label": "GAAA01000000", "data_type": "transcriptome",
                           "tax_id": "5911", "lineage": "L", "cap_mb": 2.0, "seed": 0, "disjoint_from": None,
                           "candidates": [{"codon": "TAG"}, {"codon": "TAA"}, {"codon": "TAA"}]},
                "b.json": {"accession": "GBBB01000000", "error": "download failed"},
                "c.json": {"accession": "GCCC01000000", "candidates": []},
                "d.json": {"accession": "GCA_1.1", "label": "GCA_1_1", "candidates": [{"codon": "CTG"}], "cap_mb": 4.0},
            }
            for name, r in recs.items():
                with open(os.path.join(d, name), "w", encoding="utf-8") as fh:
                    json.dump(r, fh)
            cases = run_check2.cases_from_records(d)
            self.assertEqual([c["label"] for c in cases], ["GAAA01000000_c2", "GCA_1_1_c2"])
            self.assertEqual(cases[0]["codons"], ["TAA", "TAG"])
            self.assertEqual((cases[0]["inputs"][0]["kind"], cases[0]["species_tax_id"]), ("tsa", "5911"))
            self.assertEqual((cases[1]["inputs"][0]["kind"], cases[1]["inputs"][0]["cap_mb"]), ("genome", 4.0))
            self.assertEqual([c["label"] for c in run_check2.cases_from_records(d, only="GCA_1.1")], ["GCA_1_1_c2"])

    def test_demo_grades_and_report(self):
        import run_check2
        g = support.grade_check2({0: 10, 1: 5, 2: 5, 3: 5}, {0: SAME, 1: SAME, 2: ("bacteria", None), 3: ("no_hit", None)}, 50)
        eup = {"label": "euplotes_uga", "inputs": [{"accession": "G1", "organism": "Euplotes focardii"}],
               "mapping": {"rows": 100, "mismatches": 0},
               "codons": {"TGA": {"call": "C", "total_obs": 50, "records_with_obs": 30, "check2": g}}}
        mix = {"label": "mix_uar", "inputs": [{"accession": "G2", "organism": "T"}, {"accession": "G3", "organism": "P"}],
               "mapping": {"rows": 200, "mismatches": 0},
               "codons": {"TAA": {"call": "Q", "total_obs": 40, "records_with_obs": 20, "check2": g,
                                  "origin_shares": {"ALG": 0.05, "CIL": 0.95}}}}
        gr = run_check2.grade_demo([eup, mix])
        self.assertEqual((gr["P8"]["rows"], gr["P8"]["holds"], gr["P9"]["ciliate_share"], gr["P9"]["holds"], gr["P10"]["holds"]),
                         (300, True, 0.95, True, True))
        gr2 = run_check2.grade_demo([eup, {"label": "mix_uar", "error": "download failed"}])
        self.assertEqual((gr2["P8"]["rows"], gr2["P9"]["testable"], gr2["P9"]["holds"]), (100, False, False))
        md = run_check2.report_md([eup, mix], gr, {"eukcode_version": "0.5.0", "finished_utc": "x"})
        self.assertIn("| euplotes_uga | G1 Euplotes focardii | C | TGA | 50 | 30 | pass |", md)
        self.assertIn("**P9**", md)
        self.assertIn("HOLDS", md)


if __name__ == "__main__":
    unittest.main()
