# Copyright 2026 MV2 of Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
import csv
import os
import tempfile
import unittest

from eukcode import shards


def row(acc, org, group, length, gencode="1"):
    return {"accession": acc, "organism": org, "group": group, "tax_id": "1", "length": str(length),
            "ncbi_gencode": gencode, "lineage": "x; y"}


class TestShards(unittest.TestCase):
    def setUp(self):
        self.rows = [
            row("GCA_3.1", "Zeta z", "Sar", 50_000_000),
            row("GCA_1.1", "Alpha a", "Sar", 20_000_000, "6"),
            row("GCA_2.1", "Beta b", "Discoba", 30_000_000, ""),
            row("GCA_4.1", "Huge h", "Sar", 5_000_000_000),
            row("GCA_5.1", "Tiny t", "Sar", 500_000),
            row("GCA_1.1", "Alpha a", "Sar", 20_000_000, "6"),
            row("", "Nobody", "Sar", 10_000_000),
        ]

    def test_plan(self):
        plan, skipped = shards.plan(self.rows, per_shard=2)
        flat = [r["accession"] for s in plan for r in s]
        self.assertEqual(flat, ["GCA_2.1", "GCA_1.1", "GCA_3.1"])
        self.assertEqual([len(s) for s in plan], [2, 1])
        reasons = sorted(r["skip_reason"].split()[0] for r in skipped)
        self.assertEqual(reasons, ["assembly", "assembly", "duplicate", "no"])

    def test_plan_rejects_bad_size(self):
        with self.assertRaises(ValueError):
            shards.plan(self.rows, per_shard=0)

    def test_write_and_read_back(self):
        plan, _ = shards.plan(self.rows, per_shard=2)
        with tempfile.TemporaryDirectory() as d:
            paths = shards.write_shards(plan, d)
            self.assertEqual([os.path.basename(p) for p in paths], ["shard_00.csv", "shard_01.csv"])
            back = shards.read_species(paths[0])
            self.assertEqual([r["accession"] for r in back], ["GCA_2.1", "GCA_1.1"])
            self.assertEqual(back[0]["shard"], "shard_00")
            self.assertIsNone(shards.gencode_of(back[0]))
            self.assertEqual(shards.gencode_of(back[1]), 6)
            with open(paths[1], newline="", encoding="utf-8") as fh:
                header = next(csv.reader(fh))
            self.assertEqual(header, shards.SHARD_FIELDS)


if __name__ == "__main__":
    unittest.main()
