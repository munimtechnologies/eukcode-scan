# eukcode-scan: tests for addendum A8 (KACI fast triage path).
# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
import json
import os
import tempfile
import unittest

from eukcode.kaci_io import (AAS, CODONS, STANDARD, TABLES, codon_index, compare,
                             filter_fasta, grade_control, infer_code, infer_codon,
                             is_organelle)


def entry(aa=None, n=10):
    if aa is None:
        d = {a: 1.0 / 21 for a in AAS}
    else:
        d = {a: 0.0 for a in AAS}
        d[aa] = 1.0
    d["n_"] = n
    return d


def standard_result():
    """A synthetic KACI result that reads every sense codon by the standard code."""
    res = {"version": "0.0.1", "timestamp": 0.0, "referenceFile": "x",
           "inputTranslationTable": STANDARD, "maxQueryCount": "2"}
    for i, codon in enumerate(CODONS):
        aa = STANDARD[i]
        res[codon] = entry(None, 0) if aa == "*" else entry(aa, 10)
    return res


class TestTables(unittest.TestCase):
    def test_tables(self):
        for t in TABLES.values():
            self.assertEqual(len(t), 64)
        self.assertEqual(len(CODONS), 64)
        self.assertEqual(codon_index("CTG"), 19)
        self.assertEqual(STANDARD[codon_index("CTG")], "L")
        self.assertEqual(TABLES[12][19], "S")
        self.assertEqual(TABLES[26][19], "A")
        self.assertEqual(TABLES[6][codon_index("TAA")], "Q")
        self.assertEqual(TABLES[6][codon_index("TAG")], "Q")
        self.assertEqual(TABLES[6][codon_index("TGA")], "*")
        self.assertEqual(STANDARD[codon_index("TGG")], "W")
        self.assertEqual(STANDARD[codon_index("ATG")], "M")
        # tables differ from the standard code only where stated
        diffs = {k: [CODONS[i] for i in range(64) if t[i] != STANDARD[i]] for k, t in TABLES.items()}
        self.assertEqual(diffs, {1: [], 6: ["TAA", "TAG"], 12: ["CTG"], 26: ["CTG"]})


class TestInfer(unittest.TestCase):
    def test_low_n_uncalled(self):
        self.assertEqual(infer_codon(entry("L", 3)), "?")
        self.assertEqual(infer_codon(entry("L", 4)), "L")

    def test_string_n(self):
        e = entry("L", 4)
        e["n_"] = "75"
        self.assertEqual(infer_codon(e), "L")

    def test_real_pachysolen_ctg_single_family(self):
        # CTG entry from lib-kaci test/GCA_001661245.1_Pacta1_2_genomic_expected.json
        # (one Pfam family only): A = 0.9936839919 is below 0.9999, so uncalled.
        e = {"A": 0.9936839919, "C": 0.0003425798, "D": 0.0003284289, "E": 0.0001769189,
             "F": 0.0000743049, "G": 0.0000062787, "H": 0.0013802452, "I": 0.0000778921,
             "K": 0.0002486214, "L": 0.0001489698, "M": 0.0006861751, "N": 0.0000949361,
             "P": 0.0000860523, "Q": 0.0000947116, "R": 0.0001735533, "S": 0.0000062384,
             "T": 0.0000102028, "V": 0.0001424601, "W": 0.0001071249, "Y": 0.0000257753,
             "?": 0.0021045385, "n_": 6}
        self.assertEqual(infer_codon(e), "?")
        e["A"], e["?"] = 0.99995, 0.0
        self.assertEqual(infer_codon(e), "A")

    def test_threshold_is_strict(self):
        e = entry("A", 10)
        e["A"], e["?"] = 0.9999, 0.0001
        self.assertEqual(infer_codon(e), "?")

    def test_flat_and_confident_unknown(self):
        self.assertEqual(infer_codon(entry(None, 50)), "?")
        self.assertEqual(infer_codon(entry("?", 50)), "?")

    def test_infer_code_standard(self):
        code = infer_code(standard_result())
        self.assertEqual(code, STANDARD.replace("*", "?"))

    def test_missing_codon(self):
        res = standard_result()
        del res["GGG"]
        with self.assertRaises(ValueError):
            infer_code(res)


class TestGrade(unittest.TestCase):
    def test_standard_passes(self):
        code = infer_code(standard_result())
        g = grade_control(code, {"table": 1, "must": {"CTG": "L"}})
        self.assertTrue(g["pass"])
        self.assertEqual(g["n_sense_called"], 61)
        self.assertEqual(sorted(g["comparison"]["stop_uncalled"]), ["TAA", "TAG", "TGA"])

    def test_ctg_alanine(self):
        res = standard_result()
        res["CTG"] = entry("A", 20)
        code = infer_code(res)
        self.assertFalse(grade_control(code, {"table": 1})["pass"])
        cmp = compare(code, STANDARD)
        self.assertEqual(cmp["sense_unexpected"], [["CTG", "L", "A"]])
        g26 = grade_control(code, {"table": 26, "must": {"CTG": "A"}})
        self.assertTrue(g26["pass"])
        self.assertFalse(grade_control(code, {"table": 12, "must": {"CTG": "S"}})["pass"])

    def test_must_uncalled_fails(self):
        res = standard_result()
        res["CTG"] = entry("S", 2)  # too few hits: uncalled
        code = infer_code(res)
        g = grade_control(code, {"table": 12, "must": {"CTG": "S"}})
        self.assertFalse(g["pass"])
        self.assertEqual(g["must"]["CTG"]["called"], "?")

    def test_stop_called_fails(self):
        res = standard_result()
        res["TGA"] = entry("W", 12)
        code = infer_code(res)
        g = grade_control(code, {"table": 1})
        self.assertFalse(g["pass"])
        self.assertEqual(g["comparison"]["stop_called"], [["TGA", "W"]])

    def test_ciliate_table6(self):
        res = standard_result()
        res["TAA"] = entry("Q", 30)
        res["TAG"] = entry("Q", 30)
        code = infer_code(res)
        self.assertTrue(grade_control(code, {"table": 6, "must": {"TAA": "Q", "TAG": "Q"}})["pass"])
        self.assertFalse(grade_control(code, {"table": 1})["pass"])

    def test_aga_lysine_artefact(self):
        res = standard_result()
        res["AGA"] = entry("K", 30)
        code = infer_code(res)
        g = grade_control(code, {"table": 1})
        self.assertTrue(g["pass"])
        self.assertEqual(g["comparison"]["artefacts"], [["AGA", "R", "K"]])


class TestOrganelle(unittest.TestCase):
    def test_headers(self):
        self.assertTrue(is_organelle(">NC_001224.1 Saccharomyces cerevisiae S288C mitochondrion, complete genome"))
        self.assertTrue(is_organelle(">X Plasmodium falciparum 3D7 apicoplast, complete sequence"))
        self.assertTrue(is_organelle(">Y Chloroplast genome"))
        self.assertFalse(is_organelle(">NC_001133.9 Saccharomyces cerevisiae S288C chromosome I"))

    def test_filter_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.fna")
            dst = os.path.join(d, "out.fna")
            with open(src, "w", encoding="ascii", newline="\n") as fh:
                fh.write(">chr1 nuclear\nACGT\nACGT\n>chrM mitochondrion\nTTTT\n>chr2 nuclear\nGGCC")
            counts = filter_fasta(src, dst)
            with open(dst, encoding="ascii") as fh:
                text = fh.read()
        self.assertEqual(counts["kept"], 2)
        self.assertEqual(counts["dropped"], 1)
        self.assertEqual(text, ">chr1 nuclear\nACGT\nACGT\n>chr2 nuclear\nGGCC\n")


class TestControlsFile(unittest.TestCase):
    def test_controls_json(self):
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "kaci_controls.json"), encoding="utf-8") as fh:
            spec = json.load(fh)
        ids = [c["id"] for c in spec["controls"]]
        self.assertEqual(ids, ["P16", "P17", "P18", "P19", "P20", "P21"])
        for c in spec["controls"]:
            self.assertIn(int(c["table"]), TABLES)
            for codon, aa in (c.get("must") or {}).items():
                self.assertEqual(TABLES[int(c["table"])][codon_index(codon)], aa)
        self.assertEqual(spec["kaci_commit"], "cc8601f2121ddb2546b69c5d494652c38294c02a")
        self.assertEqual(spec["reference"]["md5"], "6149c67b184cfb053bc11d44bf9b1778")


if __name__ == "__main__":
    unittest.main()
