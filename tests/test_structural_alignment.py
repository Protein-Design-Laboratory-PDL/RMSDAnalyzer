import unittest

import numpy as np

from rmsd_cal.models import ChainMatch, Residue, ResidueId, ResiduePair
from rmsd_cal.structural_alignment import StructuralAligner


def residue_pair(number, first_point, second_point):
    return ResiduePair(
        Residue(ResidueId("A", number), "ALA", {"CA": np.array(first_point, dtype=float)}),
        Residue(ResidueId("X", number), "ALA", {"CA": np.array(second_point, dtype=float)}),
    )


class StructuralAlignerTest(unittest.TestCase):
    def test_iterative_fit_rejects_large_outlier(self):
        reference = [
            (0, 0, 0),
            (2, 0, 0),
            (0, 2, 0),
            (0, 0, 2),
            (2, 2, 2),
            (1, 1, 3),
        ]
        translation = np.array((10, -4, 7))
        mobile = [tuple(np.array(point) + translation) for point in reference]
        mobile[-1] = tuple(np.array(mobile[-1]) + np.array((20, 20, 20)))
        pairs = tuple(
            residue_pair(index, first_point, second_point)
            for index, (first_point, second_point) in enumerate(zip(reference, mobile), start=1)
        )
        match = ChainMatch("A", "X", 1.0, 1.0, pairs, (), ())

        result = StructuralAligner(cutoff=8.0, max_cycles=5).align((match,), prealigned=False)

        self.assertEqual(result.fit_pair_count, 5)
        self.assertEqual(result.cycles, 1)
        self.assertAlmostEqual(result.final_rmsd, 0.0, places=6)

    def test_prealigned_mode_does_not_mark_pairs_as_fitting_atoms(self):
        points = ((0, 0, 0), (1, 0, 0), (0, 1, 0))
        pairs = tuple(
            residue_pair(index, point, point)
            for index, point in enumerate(points, start=1)
        )
        match = ChainMatch("A", "A", 1.0, 1.0, pairs, (), ())

        result = StructuralAligner().align((match,), prealigned=True)

        self.assertEqual(result.final_rmsd, 0.0)
        self.assertEqual(result.fit_residue_keys, frozenset())


if __name__ == "__main__":
    unittest.main()
