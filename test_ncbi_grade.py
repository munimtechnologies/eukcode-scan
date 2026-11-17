# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import io
import os
import tempfile
import unittest
import zipfile

from eukcode import codes, grade, ncbi

REPORT = {
    "accession": "GCF_000182965.3",
    "source_database": "SOURCE_DATABASE_REFSEQ",
    "organism": {"organism_name": "Candida albicans SC5314", "tax_id": 237561},
    "assembly_info": {"assembly_level": "Chromosome", "refseq_category": "reference genome",
                      "release_date": "2016-05-23", "assembly_name": "ASM18296v3"},
    "assembly_stats": {"total_sequence_length": "14282666", "contig_n50": 2231883},
    "annotation_info": {"name": "x"},
}

TAX_XML = b"""<?xml version="1.0" ?>
<TaxaSet><Taxon>
  <TaxId>5911</TaxId><ScientificName>Tetrahymena thermophila</ScientificName>
  <Rank>species</Rank>
  <GeneticCode><GCId>6</GCId><GCName>Ciliate Nuclear</GCName></GeneticCode>
  <MitoGeneticCode><MGCId>4</MGCId></MitoGeneticCode>
  <Lineage>cellular organisms; Eukaryota; Sar; Alveolata; Ciliophora</Lineage>
  <LineageEx><Taxon><TaxId>2759</TaxId><ScientificName>Eukaryota</ScientificName><Rank>superkingdom</Rank></Taxon></LineageEx>
</Taxon></TaxaSet>"""


class TestNcbi(unittest.TestCase):
    def test_summarize(self):
        s = ncbi.summarize_report(REPORT)
        self.assertEqual(s["accession"], "GCF_000182965.3")
        self.assertEqual(s["tax_id"], 237561)
        self.assertEqual(s["length"], 14282666)
        self.assertEqual(s["level"], "Chromosome")
        self.assertTrue(s["annotated"])

    def test_summarize_missing_fields(self):
        s = ncbi.summarize_report({"accession": "GCA_1"})
        self.assertEqual(s["length"], 0)
        self.assertEqual(s["tax_id"], 0)

    def test_best_assembly(self):
        a = {"accession": "GCA_2", "refseq_category": "", "level": "Contig", "contig_n50": 10**6, "tax_id": 1}
        b = {"accession": "GCA_3", "refseq_category": "", "level": "Scaffold", "contig_n50": 10, "tax_id": 1}
        c = {"accession": "GCA_4", "refseq_category": "reference genome", "level": "Contig", "contig_n50": 1, "tax_id": 2}
        self.assertEqual(ncbi.best_assembly([a, b])["accession"], "GCA_3")
        self.assertEqual(ncbi.best_assembly([a, b, c])["accession"], "GCA_4")
        self.assertIsNone(ncbi.best_assembly([]))
        self.assertEqual([x["accession"] for x in ncbi.one_per_taxon([a, b, c])], ["GCA_3", "GCA_4"])

    def test_taxonomy_xml_top_level_only(self):
        t = ncbi.parse_taxonomy_xml(TAX_XML)
        self.assertEqual(list(t), [5911])
        self.assertEqual(t[5911]["gencode"], 6)
        self.assertEqual(t[5911]["mito_gencode"], 4)
        self.assertIn("Ciliophora", t[5911]["lineage"])

    def test_extract_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            zp = os.path.join(d, "g.zip")
            with zipfile.ZipFile(zp, "w") as z:
                z.writestr("ncbi_dataset/data/GCA_1/GCA_1_genomic.fna", ">c1\nACGT\n")
                z.writestr("ncbi_dataset/data/assembly_data_report.jsonl", "{}")
            self.assertEqual(ncbi.fasta_members(zp), ["ncbi_dataset/data/GCA_1/GCA_1_genomic.fna"])
            out = ncbi.extract_fasta(zp, os.path.join(d, "o.fna"))
            with open(out) as fh:
                self.assertEqual(fh.read(), ">c1\nACGT\n")


class TestGrade(unittest.TestCase):
    def _code(self, changes):
        s = list(codes.STANDARD.replace("*", "?"))
        for k, v in changes.items():
            s[codes.CODON_INDEX[k]] = v
        return "".join(s)

    def test_grade_all_hold(self):
        ctrls = [
            ({"label": "scer", "table": 1, "diagnostic": {}}, {}),
            ({"label": "calb", "table": 12, "diagnostic": {"CTG": "S"}}, {"CTG": "S"}),
            ({"label": "tthe", "table": 6, "diagnostic": {"TAA": "Q", "TAG": "Q"}}, {"TAA": "Q", "TAG": "Q"}),
            ({"label": "calb_2mb", "table": 12, "diagnostic": {"CTG": "S"}, "role": "power"}, {"CTG": "S"}),
        ]
        verdicts = [grade.evaluate_control(c, self._code(ch)) for c, ch in ctrls]
        g = grade.grade(verdicts, [2.0, 3.0, 5.0])
        self.assertTrue(g["P1"]["holds"])
        self.assertTrue(g["P2"]["holds"])
        self.assertTrue(g["P3"]["holds"])
        self.assertTrue(g["P4"]["holds"])
        self.assertEqual(g["P4"]["median_min_per_mb"], 3.0)

    def test_grade_failures(self):
        ctrls = [
            ({"label": "scer", "table": 1, "diagnostic": {}}, {"CTG": "S"}),
            ({"label": "calb", "table": 12, "diagnostic": {"CTG": "S"}}, {}),
            ({"label": "tthe", "table": 6, "diagnostic": {"TAA": "Q", "TAG": "Q"}}, {"TAA": "Q"}),
        ]
        verdicts = [grade.evaluate_control(c, self._code(ch)) for c, ch in ctrls]
        g = grade.grade(verdicts, [])
        self.assertFalse(g["P1"]["holds"])  # 0 of 2 positives pass
        self.assertFalse(g["P2"]["holds"])  # scer has CTG=S
        self.assertFalse(g["P3"]["holds"])  # no power runs
        self.assertFalse(g["P4"]["holds"])


if __name__ == "__main__":
    unittest.main()
