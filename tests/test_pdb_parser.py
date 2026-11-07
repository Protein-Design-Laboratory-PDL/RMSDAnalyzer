from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from rmsd_cal.pdb_parser import PdbParser


def atom_line(serial, atom, residue, chain, number, insertion, x, y, z, altloc=""):
    return (
        f"ATOM  {serial:5d} {atom:>4}{altloc:1}{residue:>3} {chain:1}"
        f"{number:4d}{insertion:1}   {x:8.3f}{y:8.3f}{z:8.3f}"
        "  1.00 20.00           C"
    )


class PdbParserTest(unittest.TestCase):
    def test_reads_requested_model_and_insertion_code(self):
        content = "\n".join(
            [
                "MODEL        1",
                atom_line(1, "CA", "ALA", "A", 42, "", 0, 0, 0),
                "ENDMDL",
                "MODEL        2",
                atom_line(2, "CA", "GLY", "B", 42, "A", 1, 2, 3),
                "ENDMDL",
            ]
        )
        with TemporaryDirectory() as directory:
            path = Path(directory) / "models.pdb"
            path.write_text(content, encoding="utf-8")
            structure = PdbParser().parse(path, 2)

        residue = structure.chains["B"].residues[0]
        self.assertEqual(residue.id.number, 42)
        self.assertEqual(residue.id.insertion_code, "A")
        self.assertEqual(residue.id.label(), "B:42A")

    def test_prefers_blank_alternate_location_over_a(self):
        content = "\n".join(
            [
                atom_line(1, "CA", "ALA", "A", 1, "", 5, 5, 5, "A"),
                atom_line(2, "CA", "ALA", "A", 1, "", 1, 2, 3, ""),
            ]
        )
        with TemporaryDirectory() as directory:
            path = Path(directory) / "altloc.pdb"
            path.write_text(content, encoding="utf-8")
            structure = PdbParser().parse(path, 1)

        self.assertEqual(structure.chains["A"].residues[0].atoms["CA"].tolist(), [1, 2, 3])


if __name__ == "__main__":
    unittest.main()
