# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import gzip
import os
import random
import tempfile
import unittest

from eukcode import fasta


def random_records(seed=7):
    rng = random.Random(seed)
    recs = []
    for i, n in enumerate([1250, 480, 900, 1300, 160, 1010]):  # 50 usable chunks, ~4.9 kb
        seq = "".join(rng.choice("ACGT") for _ in range(n))
        if i == 2:  # a block of Ns that kills one chunk
            seq = seq[:100] + "N" * 100 + seq[200:]
        recs.append(("rec%d" % i, seq))
    return recs


class TestStreamingSample(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.recs = random_records()
        self.path = os.path.join(self.tmp.name, "g.fna")
        fasta.write_fasta(self.recs, self.path, width=60)
        self.kw = dict(chunk_bp=100, min_chunk_bp=60, max_n_frac=0.1)

    def tearDown(self):
        self.tmp.cleanup()

    def reference(self, cap, seed):
        ch = fasta.chunks(fasta.read_fasta(self.path), **self.kw)
        return fasta.sample_chunks(ch, cap, seed=seed), len(ch)

    def test_index_matches_chunks(self):
        ch = fasta.chunks(fasta.read_fasta(self.path), **self.kw)
        idx = fasta.chunk_index(fasta.read_fasta(self.path), **self.kw)
        self.assertEqual([c[0] for c in ch], [x[0] for x in idx])
        self.assertEqual([len(c[1]) for c in ch], [x[3] for x in idx])

    def test_same_sample_as_in_memory(self):
        for cap in (500, 1500, 10 ** 9):
            for seed in (0, 1, 5):
                got, n = fasta.sample_file(self.path, cap, seed=seed, **self.kw)
                ref, n_ref = self.reference(cap, seed)
                self.assertEqual(got, ref, "cap %d seed %d" % (cap, seed))
                self.assertEqual(n, n_ref)

    def test_gzip_input(self):
        gz = self.path + ".gz"
        with open(self.path, "rb") as src, gzip.open(gz, "wb") as dst:
            dst.write(src.read())
        a, _ = fasta.sample_file(self.path, 700, seed=0, **self.kw)
        b, _ = fasta.sample_file(gz, 700, seed=0, **self.kw)
        self.assertEqual(a, b)

    def test_disjoint_second_sample(self):
        first, n = fasta.sample_file(self.path, 800, seed=0, **self.kw)
        second, _ = fasta.sample_file(self.path, 800, seed=1, disjoint_from=0, **self.kw)
        names1 = set(x[0] for x in first)
        names2 = set(x[0] for x in second)
        self.assertTrue(names1)
        self.assertTrue(names2)
        self.assertFalse(names1 & names2)
        self.assertGreaterEqual(sum(len(s) for _, s in second), 800)
        # every sampled chunk is the true slice of its record
        recs = dict((str(i), s) for i, (_n, s) in enumerate(fasta.read_fasta(self.path)))
        for name, seq in first + second:
            rec, start = name[1:].split("_")
            self.assertEqual(recs[rec][int(start):int(start) + len(seq)], seq)

    def test_disjoint_when_all_fit(self):
        second, n = fasta.sample_file(self.path, 10 ** 9, seed=1, disjoint_from=0, **self.kw)
        self.assertEqual(second, [])
        self.assertGreater(n, 0)

    def test_choose(self):
        self.assertEqual(fasta.choose([10, 10], 100), [0, 1])
        p = fasta.choose([10] * 20, 50, seed=3)
        self.assertEqual(len(p), 5)
        self.assertEqual(p, sorted(p))
        self.assertEqual(p, fasta.choose([10] * 20, 50, seed=3))


if __name__ == "__main__":
    unittest.main()
