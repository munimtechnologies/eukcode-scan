# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import io
import os
import tempfile
import unittest

from eukcode import fasta


class TestFasta(unittest.TestCase):
    def test_read(self):
        h = io.StringIO(">a desc\nacgt\nNNAC\n\n>b\nGG\n")
        recs = list(fasta.read_fasta(h))
        self.assertEqual(recs, [("a", "ACGTNNAC"), ("b", "GG")])

    def test_read_file_and_write_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "x.fna")
            fasta.write_fasta([("s1", "A" * 170), ("s2", "C" * 5)], p, width=80)
            self.assertEqual(list(fasta.read_fasta(p)), [("s1", "A" * 170), ("s2", "C" * 5)])

    def test_chunks_filters(self):
        recs = [("c1", "A" * 250), ("c2", "N" * 120), ("c3", "ACGT" * 10)]
        ch = fasta.chunks(recs, chunk_bp=100, min_chunk_bp=60, max_n_frac=0.1)
        # c1 -> 100, 100, 50 (dropped); c2 all N (dropped); c3 40 (dropped)
        self.assertEqual([n for n, _ in ch], ["s0_0", "s0_100"])

    def test_sample_all_when_under_cap(self):
        ch = [("a", "A" * 10), ("b", "C" * 10)]
        self.assertEqual(fasta.sample_chunks(ch, 100), ch)

    def test_sample_deterministic_and_capped(self):
        ch = [("c%d" % i, "A" * 100) for i in range(50)]
        s1 = fasta.sample_chunks(ch, 1000, seed=0)
        s2 = fasta.sample_chunks(ch, 1000, seed=0)
        s3 = fasta.sample_chunks(ch, 1000, seed=1)
        self.assertEqual(s1, s2)
        self.assertNotEqual(s1, s3)
        self.assertEqual(sum(len(x) for _, x in s1), 1000)
        names = [n for n, _ in s1]
        self.assertEqual(names, sorted(names, key=lambda n: int(n[1:])))

    def test_gc(self):
        self.assertAlmostEqual(fasta.gc_fraction([("a", "GGCCAATT")]), 0.5)
        self.assertEqual(fasta.gc_fraction([("a", "NNNN")]), 0.0)


if __name__ == "__main__":
    unittest.main()
