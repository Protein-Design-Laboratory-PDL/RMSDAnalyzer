import unittest

from rmsd_cal.sequence import GlobalSequenceAligner


class GlobalSequenceAlignerTest(unittest.TestCase):
    def test_aligns_same_sequence_with_numbering_independent_indices(self):
        alignment = GlobalSequenceAligner().align(
            ["ALA", "GLY", "SER"],
            ["ALA", "GLY", "SER"],
        )

        self.assertEqual(alignment.index_pairs, ((0, 0), (1, 1), (2, 2)))
        self.assertEqual(alignment.identity, 1.0)

    def test_records_gap_without_pairing_inserted_residue(self):
        alignment = GlobalSequenceAligner().align(
            ["ALA", "GLY", "SER"],
            ["ALA", "VAL", "GLY", "SER"],
        )

        self.assertEqual(alignment.index_pairs, ((0, 0), (1, 2), (2, 3)))


if __name__ == "__main__":
    unittest.main()
