# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import unittest

from eukcode import codes, known


def with_changes(base_table, changes):
    code = list(codes.table(base_table))
    for codon, aa in changes.items():
        code[codes.CODON_INDEX[codon]] = aa
    return "".join(code)


TTHE_LINEAGE = ("cellular organisms; Eukaryota; Sar; Alveolata; Ciliophora; Intramacronucleata; "
                "Oligohymenophorea; Hymenostomatida; Tetrahymenina; Tetrahymenidae; Tetrahymena")
DINO_LINEAGE = "cellular organisms; Eukaryota; Sar; Alveolata; Dinophyceae; Suessiales; Symbiodiniaceae"
EUPL_LINEAGE = ("cellular organisms; Eukaryota; Sar; Alveolata; Ciliophora; Intramacronucleata; "
                "Spirotrichea; Hypotrichia; Euplotida; Euplotidae")
TRYP_LINEAGE = "cellular organisms; Eukaryota; Discoba; Euglenozoa; Kinetoplastea; Metakinetoplastina; Trypanosomatida"


class TestCatalogue(unittest.TestCase):
    def test_entries_are_valid(self):
        ids = set()
        for e in known.CATALOGUE:
            self.assertNotIn(e["id"], ids)
            ids.add(e["id"])
            self.assertTrue(e["clade"])
            self.assertTrue(e["source"])
            for codon, aa in e["changes"].items():
                self.assertEqual(codes.normalize_codon(codon), codon)
                self.assertIn(aa, "ACDEFGHIKLMNPQRSTVWY")

    def test_ncbi_tables_agree_with_catalogue(self):
        # every change in the eukaryotic NCBI tables eukcode models is in the catalogue
        cat = set()
        for e in known.CATALOGUE:
            for codon, aa in e["changes"].items():
                cat.add((codon, aa))
        for t in (6, 10, 12, 15, 26, 27, 28, 29, 30, 31):
            for d in codes.compare(codes.table(t), 1):
                self.assertIn((d["codon"], d["called"]), cat, "table %d %s" % (t, d))


class TestLineage(unittest.TestCase):
    def test_names(self):
        n = known.lineage_names("a; b ;c", "Tetrahymena thermophila SB210")
        self.assertEqual(n, {"a", "b", "c", "Tetrahymena thermophila SB210", "Tetrahymena"})
        self.assertEqual(known.lineage_names("", ""), set())


class TestClassify(unittest.TestCase):
    def test_known_taxon(self):
        c = known.classify_change("TGA", "C", EUPL_LINEAGE, "Euplotes vannus")
        self.assertEqual(c["class"], "known_taxon")
        self.assertIn("euplotes-UGA-C", c["entries"])
        c = known.classify_change("TAG", "E", TRYP_LINEAGE, "Blastocrithidia nonstop")
        self.assertEqual(c["class"], "known_taxon")

    def test_known_in_clade(self):
        # UGA=Cys in a ciliate outside the Euplotidae
        c = known.classify_change("UGA", "C", TTHE_LINEAGE, "Tetrahymena thermophila")
        self.assertEqual(c["class"], "known_in_clade")
        # UAG=Glu in a trypanosomatid that is not Blastocrithidia
        c = known.classify_change("TAG", "E", TRYP_LINEAGE, "Obscuromonas modryi")
        self.assertEqual(c["class"], "known_in_clade")

    def test_known_elsewhere_and_unreported(self):
        c = known.classify_change("TAG", "Q", DINO_LINEAGE, "Symbiodinium sp.")
        self.assertEqual(c["class"], "known_elsewhere")
        c = known.classify_change("CGG", "W", DINO_LINEAGE, "Symbiodinium sp.")
        self.assertEqual(c["class"], "unreported_type")
        self.assertEqual(c["entries"], [])


class TestCandidates(unittest.TestCase):
    def test_no_candidates_when_code_matches_ncbi(self):
        inferred = codes.table(6)
        self.assertEqual(known.candidates(inferred, 6, TTHE_LINEAGE, "Tetrahymena thermophila"), [])

    def test_noncalls_are_not_candidates(self):
        inferred = "?" * 64
        self.assertEqual(known.candidates(inferred, 1, DINO_LINEAGE, "x"), [])

    def test_ciliate_recorded_as_standard(self):
        inferred = codes.table(6)
        cands = known.candidates(inferred, 1, TTHE_LINEAGE, "Tetrahymena thermophila")
        self.assertEqual(sorted(c["codon"] for c in cands), ["TAA", "TAG"])
        self.assertTrue(all(c["class"] == "known_taxon" for c in cands))
        self.assertEqual(known.most_interesting(cands), "known_taxon")
        self.assertTrue(all(c["kind"] == "stop_to_sense" for c in cands))

    def test_mixed_classes_and_ordering(self):
        inferred = with_changes(1, {"TAG": "Q", "AGG": "M"})
        cands = known.candidates(inferred, 1, DINO_LINEAGE, "Symbiodinium sp.")
        by = {c["codon"]: c for c in cands}
        self.assertEqual(by["TAG"]["class"], "known_elsewhere")
        self.assertEqual(by["AGG"]["class"], "unreported_type")
        self.assertEqual(by["AGG"]["kind"], "sense_to_sense")
        self.assertEqual(known.most_interesting(cands), "unreported_type")
        self.assertEqual(known.most_interesting([]), "")

    def test_table_fallback(self):
        inferred = with_changes(1, {"TGA": "W"})
        cands = known.candidates(inferred, None, DINO_LINEAGE, "x")
        self.assertEqual(len(cands), 1)
        self.assertTrue(cands[0]["table_fallback"])
        self.assertEqual(cands[0]["table_used"], 1)
        cands = known.candidates(inferred, 99, DINO_LINEAGE, "x")
        self.assertTrue(cands[0]["table_fallback"])

    def test_stop_calls(self):
        inferred = with_changes(1, {"TAA": "Q"})
        inferred = inferred.replace("*", "?")
        self.assertEqual(known.stop_calls(inferred), {"TAA": "Q", "TAG": "?", "TGA": "?"})


if __name__ == "__main__":
    unittest.main()
