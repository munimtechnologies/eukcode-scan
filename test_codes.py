# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import unittest

from eukcode import codes


class TestTables(unittest.TestCase):
    def test_codon_order(self):
        self.assertEqual(len(codes.CODONS), 64)
        self.assertEqual(codes.CODONS[0], "TTT")
        self.assertEqual(codes.CODONS[10], "TAA")
        self.assertEqual(codes.CODONS[11], "TAG")
        self.assertEqual(codes.CODONS[14], "TGA")
        self.assertEqual(codes.CODONS[19], "CTG")
        self.assertEqual(codes.CODONS[63], "GGG")

    def test_standard_spot_checks(self):
        s = codes.table(1)
        self.assertEqual(s, codes.STANDARD)
        self.assertEqual(s[codes.CODON_INDEX["ATG"]], "M")
        self.assertEqual(s[codes.CODON_INDEX["TGG"]], "W")
        self.assertEqual(s[codes.CODON_INDEX["AGA"]], "R")
        self.assertEqual(s.count("*"), 3)

    def test_ncbi_strings(self):
        # Verbatim NCBI 'AAs' lines.
        ncbi = {
            4: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            6: "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            10: "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            12: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            26: "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            31: "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        }
        for tid, s in ncbi.items():
            self.assertEqual(codes.table(tid), s, tid)

    def test_unknown_table(self):
        with self.assertRaises(KeyError):
            codes.table(2)


class TestCompare(unittest.TestCase):
    def test_no_calls_never_disagree(self):
        self.assertEqual(codes.compare("?" * 64, 1), [])

    def test_identical(self):
        inferred = codes.STANDARD.replace("*", "?")
        self.assertEqual(codes.compare(inferred, 1), [])

    def test_cug_serine_against_standard(self):
        inferred = list(codes.STANDARD.replace("*", "?"))
        inferred[codes.CODON_INDEX["CTG"]] = "S"
        d = codes.compare("".join(inferred), 1)
        self.assertEqual(len(d), 1)
        self.assertEqual(d[0]["codon"], "CTG")
        self.assertEqual(d[0]["kind"], "sense_to_sense")
        self.assertEqual(codes.compare("".join(inferred), 12), [])

    def test_stop_reassignment_kind(self):
        inferred = list("?" * 64)
        inferred[codes.CODON_INDEX["TAA"]] = "Q"
        d = codes.compare("".join(inferred), 1)
        self.assertEqual(d[0]["kind"], "stop_to_sense")
        self.assertEqual(codes.compare("".join(inferred), 6), [])

    def test_dual_meaning_codon_ok(self):
        inferred = list("?" * 64)
        inferred[codes.CODON_INDEX["TAA"]] = "E"
        self.assertEqual(codes.compare("".join(inferred), 31), [])

    def test_length_check(self):
        with self.assertRaises(ValueError):
            codes.compare("?" * 63, 1)

    def test_diagnostic(self):
        inferred = list("?" * 64)
        inferred[codes.CODON_INDEX["TGA"]] = "C"
        ok, det = codes.diagnostic_check("".join(inferred), {"UGA": "C"})
        self.assertTrue(ok)
        self.assertTrue(det["TGA"]["ok"])
        ok, _ = codes.diagnostic_check("".join(inferred), {"TGA": "W"})
        self.assertFalse(ok)


if __name__ == "__main__":
    unittest.main()
