# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
import unittest

from eukcode import codes, cug

CIL = "cellular organisms; Eukaryota; Sar; Alveolata; Ciliophora; Intramacronucleata; Spirotrichea"
DIA = "cellular organisms; Eukaryota; Sar; Stramenopiles; Ochrophyta; Bacillariophyta"


def inf(table_id=6, ctg=None):
    code = list(codes.table(table_id).replace("*", "?"))
    if ctg is not None:
        code[codes.CODON_INDEX["CTG"]] = ctg
    return "".join(code)


def rec(acc, organism, lineage=CIL, ctg=None, **extra):
    r = {"accession": acc, "organism": organism, "lineage": lineage, "inferred": inf(ctg=ctg)}
    r.update(extra)
    return r


def c2(outcome):
    return {"codons": {"CTG": {"call": "S", "check2": {"outcome": outcome}}}}


class TestRows(unittest.TestCase):
    def test_rows_and_counts(self):
        recs = [rec("A", "Oxytricha trifallax"), rec("B", "Strombidium sp. X", ctg="S"),
                rec("C", "Halteria grandinella", ctg="?"), rec("D", "Phaeodactylum tricornutum", lineage=DIA, ctg="S"),
                rec("E", "Euplotes sp.", error="download failed"), {"accession": "F", "lineage": CIL, "inferred": "short"},
                rec("G", "Stylonychia lemnae", ctg="A")]
        rows = cug.ctg_rows(recs)
        self.assertEqual([r["accession"] for r in rows], ["A", "B", "C", "G"])
        self.assertEqual([r["ctg_call"] for r in rows], ["L", "S", "?", "A"])
        self.assertEqual(rows[1]["genus"], "Strombidium")
        self.assertEqual(cug.call_counts(rows), {"L": 1, "S": 1, "?": 1, "other": 1})
        self.assertFalse(cug.is_ciliate("Eukaryota; Ciliophorales"))
        self.assertEqual(cug.genus_of(""), "")


class TestGrade(unittest.TestCase):
    def setUp(self):
        self.rows = cug.ctg_rows([rec("A", "Oxytricha trifallax"), rec("B", "Strombidium sp.", ctg="S"),
                                  rec("C", "Strombidium sp. 2", ctg="S")])

    def test_pending_then_explained(self):
        g = cug.grade_p11(self.rows, {"B": c2("fail")})
        self.assertEqual((g["graded"], g["holds"], g["pending_check2"]), (False, False, ["C"]))
        g = cug.grade_p11(self.rows, {"B": c2("fail"), "C": c2("inconclusive")})
        self.assertEqual((g["graded"], g["holds"], g["ctg_s"]), (True, True, ["B", "C"]))
        self.assertEqual([e["check2"] for e in g["explained_or_unresolved"]], ["fail", "inconclusive"])

    def test_survives(self):
        g = cug.grade_p11(self.rows, {"B": c2("fail"), "C": c2("pass_unplaced")})
        self.assertEqual((g["graded"], g["holds"]), (True, False))
        self.assertEqual(g["survived_check2"], [{"accession": "C", "organism": "Strombidium sp. 2", "check2": "pass_unplaced"}])

    def test_errors_and_empty(self):
        g = cug.grade_p11(self.rows, {"B": {"error": "codetta timeout"}, "C": {"codons": {}}})
        self.assertEqual(g["pending_check2"], ["B", "C"])
        self.assertEqual(cug.grade_p11([], {})["graded"], False)
        g = cug.grade_p11(cug.ctg_rows([rec("A", "Oxytricha trifallax")]), {})
        self.assertEqual((g["graded"], g["holds"], g["ctg_calls"]), (True, True, {"L": 1}))


class TestTargets(unittest.TestCase):
    def test_status(self):
        rows = cug.ctg_rows([rec("A", "Strombidium sp.", ctg="S"), rec("B", "Halteria grandinella"),
                             rec("C", "Pseudokeronopsis sp.", ctg="?"), rec("D", "Urostyla grandis", ctg="S"),
                             rec("E", "Urostyla sp.")])
        t = {x["genus"]: x for x in cug.p10k_targets(["Strombidium", "Halteria", "Pseudokeronopsis", "Urostyla", "Favella"], rows)}
        self.assertEqual(t["Strombidium"]["status"], "agrees_with_p10k")
        self.assertEqual(t["Halteria"]["status"], "contradicts_p10k")
        self.assertEqual(t["Pseudokeronopsis"]["status"], "no_call")
        self.assertEqual((t["Urostyla"]["status"], t["Urostyla"]["datasets"], t["Urostyla"]["ctg_calls"]),
                         ("mixed", ["D", "E"], ["L", "S"]))
        self.assertEqual(t["Favella"]["status"], "not_screened")


if __name__ == "__main__":
    unittest.main()
