from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from rmsd_cal import AnalysisConfig, RmsdAnalysisService, TwoFilesInput
from rmsd_cal.report import CsvResultWriter
from tests.test_pdb_parser import atom_line


BACKBONE_OFFSETS = {
    "N": (-0.4, 0.0, 0.0),
    "CA": (0.0, 0.0, 0.0),
    "C": (0.4, 0.0, 0.0),
    "O": (0.6, 0.2, 0.0),
}


def structure_text(chain, start_number, translation):
    residues = (("ALA", (0, 0, 0)), ("GLY", (2, 1, 0)), ("SER", (1, 3, 1)))
    lines = []
    serial = 1
    for index, (name, center) in enumerate(residues):
        for atom_name, offset in BACKBONE_OFFSETS.items():
            point = tuple(center[axis] + offset[axis] + translation[axis] for axis in range(3))
            lines.append(
                atom_line(serial, atom_name, name, chain, start_number + index, "", *point)
            )
            serial += 1
    return "\n".join(lines)


class AnalysisTest(unittest.TestCase):
    def test_aligns_different_numbering_and_writes_both_metrics(self):
        with TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "first.pdb"
            second = root / "second.pdb"
            output = root / "result.csv"
            first.write_text(structure_text("A", 1, (0, 0, 0)), encoding="utf-8")
            second.write_text(structure_text("X", 101, (10, -3, 4)), encoding="utf-8")
            config = AnalysisConfig(
                input=TwoFilesInput(first, second),
                alignment_mode="align",
                output_csv=output,
            )

            result = RmsdAnalysisService().analyze(config)
            CsvResultWriter().write(result, output)

            self.assertEqual(result.chain_matches[0].first_chain_id, "A")
            self.assertEqual(result.chain_matches[0].second_chain_id, "X")
            self.assertEqual(result.residues[0].first_id.number, 1)
            self.assertEqual(result.residues[0].second_id.number, 101)
            self.assertAlmostEqual(result.alignment.final_rmsd, 0.0, places=6)
            self.assertAlmostEqual(result.residues[0].ca_rmsd, 0.0, places=6)
            self.assertAlmostEqual(result.residues[0].backbone_rmsd, 0.0, places=6)
            csv_content = output.read_text(encoding="utf-8")
            self.assertIn("ca_rmsd_angstrom", csv_content)
            self.assertIn("backbone_rmsd_angstrom", csv_content)


if __name__ == "__main__":
    unittest.main()
