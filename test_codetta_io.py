# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import unittest

from eukcode import codes, codetta_io

STD_Q = codes.STANDARD.replace("*", "?")


def fake_out(final):
    """Build text in the exact layout codetta.py writes (format strings copied from it)."""
    lines = ["# Analysis arguments",
             "alignment_prefix   /tmp/x/calb",
             "profile_database   Pfam-A_enone.hmm",
             "results_summary    /out/calb.summary.csv",
             "evalue_threshold   1e-10",
             "prob_threshold     0.9999",
             "max_fraction       0.01",
             "excluded_pfams     mtvuy",
             "#\n# Codon inferences                      Consensus columns",
             "# codon   inference  std code  diff?    N aligned  N used"]
    text = "\n".join(lines) + "\n"
    for i, c in enumerate(codes.CODONS):
        inf = final[i]
        std = codes.STANDARD[i]
        diff = "?" if inf == "?" else ("N" if inf == std else "Y")
        text += "%-10s%-11s%-10s%-9s%-11i%-10i\n" % (c, inf, std, diff, 1000 + i, 900 + i)
    text += "#\n# Log decoding probabilities\n# codon      "
    text += "".join(["%-13s"] * 21) % tuple("logP(%s)" % a for a in "ACDEFGHIKLMNPQRSTVWY?") + "\n"
    for c in codes.CODONS:
        text += "%-6s    %s\n" % (c, " ".join(["%12.4f"] * 21) % tuple([-1.5] * 21))
    text += "#\n# Final genetic code inference\n%s" % final
    return text


class TestSummary(unittest.TestCase):
    def test_with_header(self):
        best = STD_Q.lower()
        text = ("prefix,profile_db,evalue_threshold,prob_threshold,max_fraction,"
                "excluded_domains,inferred_gencode,inferred_gencode_best_models\n"
                "/tmp/x/calb,Pfam-A_enone.hmm,1e-10,0.9999,0.01,mtvuy,%s,%s\n" % (STD_Q, best))
        rows = codetta_io.parse_summary(text)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["inferred_gencode"], STD_Q)
        self.assertEqual(rows[0]["excluded_domains"], "mtvuy")

    def test_bad_length(self):
        text = "prefix,profile_db,evalue_threshold,prob_threshold,max_fraction,excluded_domains,inferred_gencode,inferred_gencode_best_models\na,b,1,1,1,x,FFL,FFL\n"
        with self.assertRaises(ValueError):
            codetta_io.parse_summary(text)


class TestOut(unittest.TestCase):
    def test_parse(self):
        final = list(STD_Q)
        final[codes.CODON_INDEX["CTG"]] = "S"
        final = "".join(final)
        p = codetta_io.parse_inference_out(fake_out(final))
        self.assertEqual(p["final"], final)
        self.assertEqual(p["args"]["prob_threshold"], "0.9999")
        self.assertEqual(len(p["codons"]), 64)
        self.assertEqual(p["codons"]["CTG"]["inference"], "S")
        self.assertEqual(p["codons"]["CTG"]["diff"], "Y")
        self.assertEqual(p["codons"]["TTT"]["n_aligned"], 1000)
        self.assertEqual(p["codons"]["GGG"]["n_used"], 963)
        self.assertEqual(len(p["logp"]), 64)
        self.assertEqual(len(p["logp"]["ATG"]), 21)
        self.assertEqual(p["logp_header"][0], "logP(A)")

    def test_rna_codons_accepted(self):
        text = fake_out(STD_Q).replace("TTT       ", "UUU       ", 1)
        p = codetta_io.parse_inference_out(text)
        self.assertIn("TTT", p["codons"])


if __name__ == "__main__":
    unittest.main()
