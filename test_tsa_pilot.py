# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Code written by AI agents working in MV2's name.
import json
import os
import unittest

from eukcode import codes, tsa_pilot

HERE = os.path.dirname(os.path.abspath(__file__))


def inf(table_id, **calls):
    """A Codetta-like inferred code: the NCBI table with stops as '?', plus codon overrides."""
    code = list(codes.table(table_id).replace("*", "?"))
    for codon, aa in calls.items():
        code[codes.CODON_INDEX[codon]] = aa
    return "".join(code)


def st(label, role, resolved=True, verdict=None):
    return {"label": label, "role": role, "resolved": resolved, "verdict": verdict}


class TestControls(unittest.TestCase):
    def test_real_control_file(self):
        with open(os.path.join(HERE, "tsa_controls.json"), encoding="utf-8") as fh:
            ctrls = tsa_pilot.control_list(json.load(fh))
        self.assertEqual(len(ctrls), 16)
        self.assertEqual(sum(1 for c in ctrls if c["role"] == "positive"), 8)
        self.assertEqual(sum(1 for c in ctrls if c["role"] == "standard"), 8)
        tthe = [c for c in ctrls if c["label"] == "t_tthe"][0]
        self.assertEqual((tthe["table"], tthe["diagnostic"]), (6, {"TAA": "Q", "TAG": "Q"}))
        self.assertTrue(all(c["table"] == 1 for c in ctrls if c["role"] == "standard"))
        for c in ctrls:  # every registered table is one eukcode models
            codes.table(c["table"])

    def test_bad_control_files(self):
        with self.assertRaises(ValueError):
            tsa_pilot.control_list({"positive": [{"label": "a", "taxa": ["x"], "table": 6, "diagnostic": {"TAA": "Q"}}],
                                    "standard": [{"label": "a", "taxa": ["y"], "table": 1}]})
        with self.assertRaises(ValueError):
            tsa_pilot.control_list({"positive": [{"label": "a", "taxa": ["x"], "table": 6}]})


class TestEvaluate(unittest.TestCase):
    tthe = {"label": "t_tthe", "role": "positive", "table": 6, "diagnostic": {"TAA": "Q", "TAG": "Q"}}
    crei = {"label": "t_crei", "role": "standard", "table": 1, "diagnostic": {}}

    def test_positive(self):
        v = tsa_pilot.evaluate(self.tthe, inf(6))
        self.assertTrue(v["diagnostic_ok"])
        self.assertEqual(v["disagreements"], [])
        self.assertEqual(v["n_called"], 63)  # only TGA (stop in table 6) is '?'
        v = tsa_pilot.evaluate(self.tthe, inf(6, TAA="?"))
        self.assertFalse(v["diagnostic_ok"])
        self.assertEqual(v["diagnostic"]["TAA"], {"expected": "Q", "called": "?", "ok": False})

    def test_standard(self):
        self.assertEqual(tsa_pilot.evaluate(self.crei, inf(1))["disagreements"], [])
        d = tsa_pilot.evaluate(self.crei, inf(1, TGG="C"))["disagreements"]
        self.assertEqual(d, [{"codon": "TGG", "called": "C", "expected": "W", "kind": "sense_to_sense"}])
        d = tsa_pilot.evaluate(self.crei, inf(1, TGA="W"))["disagreements"]
        self.assertEqual(d[0]["kind"], "stop_to_sense")


OK = {"diagnostic_ok": True, "disagreements": []}
BAD = {"diagnostic_ok": False, "disagreements": []}


def std_v(n_bad):
    return {"diagnostic_ok": True, "disagreements": [{"codon": "TGG", "called": "C", "expected": "W"}] * n_bad}


class TestGrades(unittest.TestCase):
    def test_p5(self):
        s = [st("a", "positive", verdict=OK), st("b", "positive", verdict=OK), st("c", "positive", verdict=OK),
             st("d", "positive", verdict=BAD), st("e", "positive", resolved=False), st("x", "standard", verdict=OK)]
        g = tsa_pilot.grade_p5(s)
        self.assertEqual((g["resolved"], g["passed"], g["failed_labels"], g["testable"], g["holds"]),
                         (4, 3, ["d"], True, True))
        g = tsa_pilot.grade_p5(s + [st("f", "positive", verdict=None)])  # unfinished counts as a failure
        self.assertEqual((g["resolved"], g["passed"], g["unfinished_labels"], g["holds"]), (5, 3, ["f"], False))
        g = tsa_pilot.grade_p5([st("a", "positive", verdict=OK), st("b", "positive", verdict=OK)])
        self.assertEqual((g["testable"], g["holds"]), (False, False))

    def test_p6(self):
        s = [st("a", "standard", verdict=std_v(0)), st("b", "standard", verdict=std_v(1)),
             st("c", "standard", verdict=std_v(0)), st("d", "standard", resolved=False)]
        g = tsa_pilot.grade_p6(s)
        self.assertEqual((g["finished"], g["disagreements"], g["holds"]), (3, 1, True))
        self.assertEqual(list(g["details"]), ["b"])
        g = tsa_pilot.grade_p6(s + [st("e", "standard", verdict=std_v(1))])
        self.assertEqual((g["disagreements"], g["holds"]), (2, False))
        g = tsa_pilot.grade_p6([st("a", "standard", verdict=std_v(0)), st("b", "standard", verdict=std_v(0)),
                                st("c", "standard", verdict=None)])
        self.assertEqual((g["finished"], g["unfinished_labels"], g["testable"], g["holds"]), (2, ["c"], False, False))

    def test_gate(self):
        yes = {"testable": True, "holds": True}
        no = {"testable": True, "holds": False}
        untestable = {"testable": False, "holds": False}
        self.assertEqual(tsa_pilot.gate(yes, yes, True, True), {"tsa_screen_may_start": True, "reasons": []})
        g = tsa_pilot.gate(yes, yes)
        self.assertFalse(g["tsa_screen_may_start"])
        self.assertEqual(len(g["reasons"]), 2)
        g = tsa_pilot.gate(no, untestable, True, False)
        self.assertEqual(g["reasons"], ["P2 does not hold", "P5 does not hold", "P6 cannot be tested (too few controls)"])

    def test_p7(self):
        lin = "Eukaryota; Sar; Stramenopiles"
        cil = "Eukaryota; Sar; Alveolata; Ciliophora; Spirotrichea"
        recs = [{"accession": "A%d" % i, "tier": "new", "lineage": lin, "candidates": []} for i in range(4)]
        recs.append({"accession": "B", "tier": "new", "lineage": lin, "candidates": [{"codon": "TGA"}]})
        recs.append({"accession": "C", "tier": "new", "lineage": cil, "candidates": [{"codon": "TAA"}]})
        recs.append({"accession": "D", "tier": "new", "lineage": lin, "error": "download failed"})
        recs.append({"accession": "E", "tier": "pre2016", "lineage": lin, "candidates": [{"codon": "TGA"}]})
        g = tsa_pilot.grade_p7(recs)
        self.assertEqual((g["eligible"], g["no_candidate"], g["fraction"], g["with_candidates"], g["holds"]),
                         (5, 4, 0.8, ["B"], True))
        recs[0]["candidates"] = [{"codon": "CTG"}]
        g = tsa_pilot.grade_p7(recs)
        self.assertEqual((g["fraction"], g["holds"]), (0.6, False))
        self.assertEqual(tsa_pilot.grade_p7([])["testable"], False)
        self.assertTrue(tsa_pilot.is_ciliate(cil))
        self.assertFalse(tsa_pilot.is_ciliate("Eukaryota; Ciliophorales"))

    def test_report(self):
        s = [st("t_tthe", "positive", verdict=tsa_pilot.evaluate(TestEvaluate.tthe, inf(6)))]
        s[0].update({"accession": "GAAA01000000", "organism": "Tetrahymena thermophila", "ncbi_gencode": "6", "table": 6})
        p5, p6 = tsa_pilot.grade_p5(s), tsa_pilot.grade_p6(s)
        md = tsa_pilot.report_md(s, p5, p6, tsa_pilot.gate(p5, p6), {"eukcode_version": "0.4.0"})
        self.assertIn("| t_tthe | positive | GAAA01000000 | Tetrahymena thermophila | 6 | 6 | 63 |", md)
        self.assertIn("TAA=Q(ok), TAG=Q(ok)", md)
        self.assertIn("CANNOT BE TESTED", md)
        self.assertIn("**Gate for the TSA screen:** closed", md)


class TestParser(unittest.TestCase):
    def test_controls_mode(self):
        import run_tsa
        a = run_tsa.build_parser().parse_args(["controls", "--p1-holds", "yes", "--out-dir", "o"])
        self.assertEqual((a.mode, a.p1_holds, a.p2_holds, a.cap_mb, a.seed, a.min_bp), ("controls", True, None, 2.0, 0, 1_000_000))
        a = run_tsa.build_parser().parse_args(["screen", "--census", "c.csv", "--shard-index", "3"])
        self.assertEqual((a.mode, a.shard_index, a.per_shard, a.cap_mb), ("screen", 3, 12, 2.0))
        with self.assertRaises(SystemExit):
            run_tsa.build_parser().parse_args(["screen"])


if __name__ == "__main__":
    unittest.main()
