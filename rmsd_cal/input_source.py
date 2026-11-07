"""Input definitions for two PDB files or two models in one PDB file."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, Tuple

from .models import Structure
from .pdb_parser import PdbParser


class StructureInput(Protocol):
    def load(self, parser: PdbParser) -> Tuple[Structure, Structure]: ...


@dataclass(frozen=True)
class TwoFilesInput:
    first_file: Path
    second_file: Path
    first_model: int = 1
    second_model: int = 1

    def load(self, parser: PdbParser) -> Tuple[Structure, Structure]:
        return (
            parser.parse(self.first_file, self.first_model),
            parser.parse(self.second_file, self.second_model),
        )


@dataclass(frozen=True)
class TwoModelsInput:
    file: Path
    first_model: int = 1
    second_model: int = 2

    def load(self, parser: PdbParser) -> Tuple[Structure, Structure]:
        return (
            parser.parse(self.file, self.first_model),
            parser.parse(self.file, self.second_model),
        )
