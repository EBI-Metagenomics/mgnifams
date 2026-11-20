#!/usr/bin/env python3
"""Tests for classify_multi_hits.py: model spans from RF, row names, copy chaining, protein classes, one full MSA.
Run: python3 bin/test_classify_multi_hits.py"""

import contextlib
import gzip
import io
import tempfile
import unittest
from pathlib import Path

from classify_multi_hits import chain, classify, ends, is_full, main, model_span, region, summarise

M = 100

# Model of 10 positions between two leading and two trailing insert columns.
RF = "..xxxxxxxxxx.."
ROWS = [
    ("1/1-10", "..ABCDEFGHIJ.."),  # single, full
    ("2/1-10", "..ABCDEFGHIJ.."),  # tandem: two full copies
    ("2/21-30", "..ABCDEFGHIJ.."),
    ("3/1-5", "..ABCDE-----.."),  # split: model 1-5 then 6-10, one full copy
    ("3/9-13", "..-----FGHIJ.."),
    ("4/2-7", "..-BCDEFG---.."),  # single, partial (jagged both ends)
    ("5/1-10", "..ABCDEFGHIJ.."),  # mixed: full copy + partial copy
    ("5/30-33", "..ABCD------.."),
    ("6/1-6", "..ABCDEF----.."),  # partial repeat: model 1-6 twice
    ("6/20-25", "..ABCDEF----.."),
    ("7", "abCDEFGHIJKLmn"),  # whole protein, inserts not counted: single, full
]


def sto(rows, rf=RF):
    width = max(len(name) for name, _ in rows) + 1
    body = "".join(f"{name:<{width}}{aligned}\n" for name, aligned in rows)
    return f"# STOCKHOLM 1.0\n#=GF ID fam\n\n{body}{'#=GC RF':<{width}}{rf}\n//\n"


class Spans(unittest.TestCase):
    def test_match_columns_only(self):
        self.assertEqual(model_span("..ABCDEFGHIJ..", RF), (1, 10))
        self.assertEqual(model_span("ab-BCDEFGHI-mn", RF), (2, 9))
        self.assertEqual(model_span("..-----FGHIJ..", RF), (6, 10))

    def test_inserts_only_row_has_no_span(self):
        self.assertIsNone(model_span("ab----------mn", RF))

    def test_region_names(self):
        self.assertEqual(region("1218837520/290-317", "x"), ("1218837520", 290, 317))
        self.assertEqual(region("7", "abCD-EF..g"), ("7", 1, 7))

    def test_ends_follow_pfam_with_5_percent(self):
        self.assertEqual(ends(1, M, M), ("curved", "curved"))
        self.assertEqual(ends(6, M - 5, M), ("curved", "curved"))  # 5 of 100 positions missing at each end
        self.assertEqual(ends(7, M, M), ("jagged", "curved"))
        self.assertEqual(ends(1, M - 6, M), ("curved", "jagged"))
        self.assertEqual(ends(3, 275, 278), ("curved", "curved"))  # MGnifam 25's widest rows


class Classes(unittest.TestCase):
    def run_case(self, rows):
        copies = chain(rows, M)
        return classify(copies, M), copies

    def test_single(self):
        self.assertEqual(self.run_case([(1, 100, 1, 100)])[0], "single")

    def test_tandem(self):
        cls, copies = self.run_case([(120, 215, 3, 98), (1, 100, 1, 100)])
        self.assertEqual((cls, len(copies)), ("tandem", 2))

    def test_split_is_one_copy(self):
        cls, copies = self.run_case([(1, 50, 1, 50), (90, 140, 51, 100)])
        self.assertEqual((cls, len(copies), is_full(copies[0], M)), ("split", 1, True))

    def test_split_allows_small_overlap(self):
        # overlap 6 <= 20% of the shorter row (51)
        self.assertEqual(self.run_case([(1, 55, 1, 55), (80, 130, 50, 100)])[0], "split")

    def test_split_partial(self):
        cls, copies = self.run_case([(1, 30, 1, 30), (50, 80, 31, 60)])
        self.assertEqual((cls, is_full(copies[0], M)), ("split", False))

    def test_split_up_to_half_overlap(self):
        # 4542 MGYP 1470485128 (model 106): overlap 17 of 58 (29%), adjacent on the protein
        copies = chain([(2, 164, 1, 65), (160, 261, 49, 106)], 106)
        self.assertEqual((classify(copies, 106), len(copies), is_full(copies[0], 106)), ("split", 1, True))

    def test_split_full_by_ends_not_covered_positions(self):
        # 14909 MGYP 935219005 (model 223): model 70-130 has no residues; ends 2 and 223
        copies = chain([(5, 74, 2, 69), (55, 150, 131, 223)], 223)
        self.assertEqual((classify(copies, 223), len(copies), is_full(copies[0], 223)), ("split", 1, True))

    def test_large_overlap_is_two_copies(self):
        # overlap 36 > 50% of the shorter row (60)
        cls, copies = self.run_case([(1, 60, 1, 60), (80, 140, 25, 100)])
        self.assertEqual((cls, len(copies)), ("partial_repeat", 2))

    def test_mixed(self):
        self.assertEqual(self.run_case([(1, 100, 1, 100), (130, 170, 1, 40)])[0], "mixed")

    def test_partial_repeat(self):
        self.assertEqual(self.run_case([(1, 50, 10, 60), (80, 135, 15, 70)])[0], "partial_repeat")

    def test_model_order_must_follow_protein_order(self):
        # C-terminal half of the model first on the protein: not one split copy
        cls, copies = self.run_case([(1, 50, 51, 100), (60, 110, 1, 50)])
        self.assertEqual((cls, len(copies)), ("partial_repeat", 2))

    def test_tandem_with_one_split_copy(self):
        cls, copies = self.run_case([(1, 100, 1, 100), (150, 200, 1, 50), (230, 280, 51, 100)])
        self.assertEqual((cls, len(copies)), ("tandem", 2))


class File(unittest.TestCase):
    def write(self, folder, text, name="fam.sto.gz"):
        path = Path(folder) / name
        with gzip.open(path, "wt") as handle:
            handle.write(text)
        return path

    def test_summary_and_rows(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.write(folder, sto(ROWS))
            rows = io.StringIO()
            got = summarise(path, rows)
        want = dict(
            family="fam",
            model_length=10,
            regions=11,
            mgyps=7,
            full_length_regions=5,
            single_full=2,
            single_partial=1,
            multi_hit_mgyps=4,
            split=1,
            split_full=1,
            tandem=1,
            mixed=1,
            partial_repeat=1,
            copies=10,
            full_copies=6,
        )
        self.assertEqual(got, want)
        lines = rows.getvalue().splitlines()
        self.assertEqual(len(lines), 11)
        self.assertIn("fam\t3\t9\t13\t6\t10\tjagged\tcurved\t1\tTrue\tsplit", lines)
        self.assertIn("fam\t4\t2\t7\t2\t7\tjagged\tjagged\t1\tFalse\tsingle", lines)

    def test_main_writes_header_and_totals(self):
        with tempfile.TemporaryDirectory() as folder:
            self.write(folder, sto(ROWS))
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                main([folder])
        text = out.getvalue().splitlines()
        self.assertTrue(text[0].startswith("family\tmodel_length\tregions"))
        self.assertTrue(text[1].startswith("fam\t10\t11\t7\t5"))
        self.assertTrue(text[-1].startswith("# total families=1 regions=11 mgyps=7"))

    def test_repeated_row_name_fails(self):
        # a second Stockholm block repeats names; only single-block (pfam format) files are read
        with tempfile.TemporaryDirectory() as folder:
            path = self.write(folder, sto(ROWS[:1] + ROWS[:1]))
            with self.assertRaises(ValueError):
                summarise(path)

    def test_missing_rf_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            path = self.write(folder, sto(ROWS).replace("#=GC RF", "#=GC XX"))
            with self.assertRaises(ValueError):
                summarise(path)


if __name__ == "__main__":
    unittest.main()
