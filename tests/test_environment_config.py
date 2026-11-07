from pathlib import Path
import unittest

from rmsd_cal.environment_config import EnvironmentConfigLoader
from rmsd_cal.errors import AnalysisError
from rmsd_cal.input_source import TwoFilesInput, TwoModelsInput


class EnvironmentConfigLoaderTest(unittest.TestCase):
    def setUp(self):
        self.base_values = {
            "INPUT_MODE": "two_models",
            "PDB_FILE_1": "data/models.pdb",
            "MODEL_1": "1",
            "MODEL_2": "2",
            "ALIGNMENT_MODE": "align",
            "TARGET_RESIDUES_CSV": "A:42,A:42A",
            "TOP_COUNT": "10",
            "OUTPUT_CSV": "output/result.csv",
        }

    def test_loads_two_model_settings(self):
        config = EnvironmentConfigLoader().load(self.base_values)

        self.assertIsInstance(config.input, TwoModelsInput)
        self.assertEqual(config.target_residues, ("A:42", "A:42A"))
        self.assertEqual(config.output_csv, Path("output/result.csv"))

    def test_loads_two_file_settings(self):
        values = dict(self.base_values)
        values.update({"INPUT_MODE": "two_files", "PDB_FILE_2": "data/second.pdb"})

        config = EnvironmentConfigLoader().load(values)

        self.assertIsInstance(config.input, TwoFilesInput)
        self.assertEqual(config.input.second_file, Path("data/second.pdb"))

    def test_rejects_missing_second_file(self):
        values = dict(self.base_values)
        values["INPUT_MODE"] = "two_files"

        with self.assertRaisesRegex(AnalysisError, "PDB_FILE_2"):
            EnvironmentConfigLoader().load(values)


if __name__ == "__main__":
    unittest.main()
