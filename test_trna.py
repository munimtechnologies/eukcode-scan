# Tests for eukcode.trna (A7 tRNA test). Hand-made tRNAscan-SE 2.0 output; no tools, no network.
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
import unittest

from eukcode import trna

HEAD_DETAIL = [
    "\t".join(["Sequence", "", "tRNA", "Bounds", "tRNA", "Anti", "Intron Bounds", "Inf", "HMM",
               "2'Str", "Isotype", "Isotype", "      "]),
    "\t".join(["Name    ", "tRNA #", "Begin", "End", "Type", "Codon", "Begin", "End", "Score",
               "Score", "Score", "CM", "Score", "Note"]),
    "\t".join(["--------", "------", "-----", "------", "----", "-----", "-----", "----", "------",
               "-----", "-----", "-------", "-------", "------"]),
]

R_C1 = ["c1      ", "1", "100", "182", "Leu", "CAG", "0", "0", "78.2", "47.50", "30.70", "Leu", "124.9", ""]
R_C2 = ["c2      ", "1", "500", "418", "Leu", "CAG", "0", "0", "60.1", "40.00", "20.10", "Ser", "95.0", "IPD"]
R_C3T = ["c3      ", "2", "10", "90", "Leu", "TAG", "0", "0", "70.0", "45.0", "25.0", "Leu", "110.0", ""]
R_C3P = ["c3      ", "3", "200", "280", "Leu", "CAG", "0", "0", "30.0", "20.0", "10.0", "Leu", "60.0", "pseudo"]
R_C4 = ["c4      ", "1", "5", "80", "Ala", "TGC", "0", "0", "67.1", "46.20", "20.90", "Ala", "77.1", ""]


def table(*rows):
    return "\n".join(HEAD_DETAIL + ["\t".join(r) for r in rows]) + "\n"


def summary_of(*rows, contigs=None):
    return trna.summarize(trna.parse_table(table(*rows)), contigs=contigs)


class ParseTableTests(unittest.TestCase):
    def test_detail_rows(self):
        rows = trna.parse_table(table(R_C1, R_C2, R_C3T, R_C3P, R_C4))
        self.assertEqual(len(rows), 5)
        self.assertEqual(rows[0]["id"], "c1.trna1")
        self.assertEqual(rows[0]["isotype_cm"], "Leu")
        self.assertAlmostEqual(rows[0]["isotype_score"], 124.9)
        self.assertEqual(rows[1]["strand"], "-")
        self.assertEqual(rows[1]["note"], "IPD")
        self.assertEqual(rows[2]["anticodon"], "TAG")
        self.assertEqual(rows[3]["note"], "pseudo")

    def test_no_detail(self):
        text = "\n".join([
            "\t".join(["Sequence", "", "tRNA", "Bounds", "tRNA", "Anti", "Intron Bounds", "Inf", ""]),
            "\t".join(["Name", "tRNA #", "Begin", "End", "Type", "Codon", "Begin", "End", "Score", "Note"]),
            "\t".join(["--------", "------", "-----", "------", "----", "-----", "-----", "----", "------", "------"]),
            "\t".join(["c1", "1", "100", "182", "Leu", "CAG", "0", "0", "78.2", ""]),
        ])
        rows = trna.parse_table(text)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["isotype_cm"], "")
        self.assertEqual(trna.body_class(rows[0]), "no-isotype")

    def test_brief_without_header(self):
        rows = trna.parse_table("\t".join(R_C2) + "\n")
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["isotype_cm"], "Ser")
        self.assertEqual(trna.body_class(rows[0]), "Ser-type")


class StructAndContigTests(unittest.TestCase):
    def test_parse_struct(self):
        text = "\n".join([
            "c1.trna1 (100-182)\tLength: 83 bp",
            "Type: Leu\tAnticodon: CAG at 35-37 (134-136)\tScore: 78.2",
            "HMM Sc=47.50\tSec struct Sc=30.70",
            "         *    |    *    |",
            "Seq: GTCAGGATGGCCGAGCGG",
            "Str: >>>>>>>..>>>>.....",
            "",
            "c2.trna1 (500-418)\tLength: 83 bp",
            "Seq: AAAACCA",
            "Str: .......",
        ])
        s = trna.parse_struct(text)
        self.assertEqual(sorted(s), ["c1.trna1", "c2.trna1"])
        self.assertEqual(s["c1.trna1"]["seq"], "GTCAGGATGGCCGAGCGG")
        self.assertEqual(s["c1.trna1"]["length"], 83)
        self.assertEqual(s["c2.trna1"]["str"], ".......")

    def test_contig_stats(self):
        seq = "CCCCAAAACCCC" + "AT" * 100 + "GGGG"
        st = trna.contig_stats(seq)
        self.assertEqual(st["length"], 216)
        self.assertTrue(st["tel5"])
        self.assertFalse(st["tel3"])
        self.assertAlmostEqual(st["gc"], 0.0556)


class SummaryAndGradeTests(unittest.TestCase):
    def test_summary_counts(self):
        s = summary_of(R_C1, R_C2, R_C3T, R_C3P, R_C4)
        self.assertEqual(s["trnas"], 4)
        self.assertEqual(s["pseudo"], 1)
        self.assertEqual(s["ipd"], 1)
        self.assertEqual(s["isotype_counts"], {"Ala": 1, "Leu": 3})
        self.assertEqual((s["cag"], s["cag_leu"], s["cag_ser"], s["cag_other"]), (2, 1, 1, 0))
        self.assertEqual((s["tag"], s["tag_leu"], s["tag_ser"]), (1, 1, 0))
        self.assertEqual(len(s["decoders"]), 4)
        self.assertEqual(s["cag_ser_telomere"], 0)
        self.assertEqual(trna.euplotes_outcome(s), "unattributed-ser-type")
        self.assertEqual(trna.grade_positive(s), "pass")
        self.assertEqual(trna.grade_negative(s), "fail")

    def test_candidate_needs_telomere(self):
        contigs = {"c2": {"original": "X", "length": 1000, "gc": 0.3, "tel5": True, "tel3": False}}
        s = summary_of(R_C1, R_C2, contigs=contigs)
        self.assertEqual(s["cag_ser_telomere"], 1)
        self.assertEqual(trna.euplotes_outcome(s), "candidate")

    def test_outcome_classes(self):
        self.assertEqual(trna.euplotes_outcome(summary_of(R_C1)), "cug-leu-supported")
        self.assertEqual(trna.grade_negative(summary_of(R_C1)), "pass")
        self.assertEqual(trna.euplotes_outcome(summary_of(R_C3T)), "no-cag-trna-leu-u34")
        self.assertEqual(trna.grade_negative(summary_of(R_C3T)), "uninformative")
        self.assertEqual(trna.euplotes_outcome(summary_of(R_C4)), "no-cug-decoder-found")
        self.assertEqual(trna.grade_positive(summary_of(R_C1)), "fail")
        self.assertEqual(trna.euplotes_outcome(None), "not run")

    def test_grade_a7_not_supported_and_contaminant(self):
        ser, leu = summary_of(R_C2), summary_of(R_C1)
        g = trna.grade_a7({
            "calb": {"role": "positive", "summary": ser},
            "dhan": {"role": "positive", "summary": ser},
            "tthe": {"role": "negative", "summary": leu},
            "efoc": {"role": "main", "summary": leu},
            "hbjk": {"role": "transcriptome", "summary": ser},
        })
        self.assertEqual(g["P12"], "pass")
        self.assertEqual(g["P13"], "pass")
        self.assertEqual(g["P14"], {"efoc": "cug-leu-supported"})
        self.assertEqual(g["transcriptomes"], {"hbjk": "unattributed-ser-type"})
        self.assertEqual(len(g["notes"]), 2)
        self.assertIn("not supported at the tRNA level", g["notes"][0])
        self.assertIn("contaminant", g["notes"][1])

    def test_grade_a7_control_fails(self):
        leu = summary_of(R_C1)
        g = trna.grade_a7({
            "calb": {"role": "positive", "summary": leu},
            "tthe": {"role": "negative", "summary": None},
            "efoc": {"role": "main", "summary": leu},
        })
        self.assertEqual(g["P12"], "fail")
        self.assertEqual(g["P13"], "not run")
        self.assertTrue(g["notes"][0].startswith("P12 did not pass"))


if __name__ == "__main__":
    unittest.main()
