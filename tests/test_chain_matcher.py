from pathlib import Path
import unittest

import numpy as np

from rmsd_cal.chain_matcher import ChainMatcher
from rmsd_cal.models import Chain, Residue, ResidueId, Structure
from rmsd_cal.sequence import GlobalSequenceAligner


def chain(chain_id, centers):
    names = ("ALA", "GLY", "SER")
    residues = []
    for index, (name, center) in enumerate(zip(names, centers), start=1):
        residues.append(
            Residue(ResidueId(chain_id, index), name, {"CA": np.array(center, dtype=float)})
        )
    return Chain(chain_id, residues)


class ChainMatcherTest(unittest.TestCase):
    def test_sequence_tie_is_resolved_by_combined_structure_rmsd(self):
        first = Structure(
            Path("first.pdb"),
            1,
            {
                "A": chain("A", [(0, 0, 0), (1, 0, 0), (0, 2, 0)]),
                "B": chain("B", [(8, 1, 1), (9, 1, 1), (8, 4, 2)]),
            },
        )
        translation = np.array((20, -7, 3))
        second = Structure(
            Path("second.pdb"),
            1,
            {
                "X": chain("X", [tuple(np.array(point) + translation) for point in [(8, 1, 1), (9, 1, 1), (8, 4, 2)]]),
                "Y": chain("Y", [tuple(np.array(point) + translation) for point in [(0, 0, 0), (1, 0, 0), (0, 2, 0)]]),
            },
        )

        matches = ChainMatcher(GlobalSequenceAligner()).match(first, second)

        self.assertEqual(
            [(match.first_chain_id, match.second_chain_id) for match in matches],
            [("A", "Y"), ("B", "X")],
        )


if __name__ == "__main__":
    unittest.main()
