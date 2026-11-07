"""Translate settings exported by run_analysis.sh into application config."""

from __future__ import annotations

from os import environ
from pathlib import Path
from typing import Mapping

from .analysis import AnalysisConfig
from .errors import AnalysisError
from .input_source import TwoFilesInput, TwoModelsInput


class EnvironmentConfigLoader:
    def load(self, values: Mapping[str, str] = environ) -> AnalysisConfig:
        input_mode = self._required(values, "INPUT_MODE")
        first_file = Path(self._required(values, "PDB_FILE_1"))
        first_model = self._positive_integer(values, "MODEL_1")
        second_model = self._positive_integer(values, "MODEL_2")

        if input_mode == "two_models":
            structure_input = TwoModelsInput(first_file, first_model, second_model)
        elif input_mode == "two_files":
            structure_input = TwoFilesInput(
                first_file,
                Path(self._required(values, "PDB_FILE_2")),
                first_model,
                second_model,
            )
        else:
            raise AnalysisError(
                "INPUT_MODE must be 'two_models' or 'two_files', "
                f"got {input_mode!r}."
            )

        targets = tuple(
            target.strip()
            for target in values.get("TARGET_RESIDUES_CSV", "").split(",")
            if target.strip()
        )
        return AnalysisConfig(
            input=structure_input,
            alignment_mode=self._required(values, "ALIGNMENT_MODE"),
            target_residues=targets,
            top_count=self._positive_integer(values, "TOP_COUNT"),
            output_csv=Path(self._required(values, "OUTPUT_CSV")),
        )

    @staticmethod
    def _required(values: Mapping[str, str], name: str) -> str:
        value = values.get(name, "").strip()
        if not value:
            raise AnalysisError(f"Required setting {name} is missing.")
        return value

    def _positive_integer(self, values: Mapping[str, str], name: str) -> int:
        raw_value = self._required(values, name)
        try:
            value = int(raw_value)
        except ValueError as exc:
            raise AnalysisError(f"{name} must be an integer, got {raw_value!r}.") from exc
        if value < 1:
            raise AnalysisError(f"{name} must be at least 1, got {value}.")
        return value
