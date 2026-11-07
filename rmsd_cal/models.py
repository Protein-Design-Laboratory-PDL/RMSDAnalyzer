"""Domain models used by the RMSD analysis pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class ResidueId:
    chain_id: str
    number: int
    insertion_code: str = ""

    def label(self) -> str:
        return f"{self.chain_id}:{self.number}{self.insertion_code}"


@dataclass
class Residue:
    id: ResidueId
    name: str
    atoms: Dict[str, np.ndarray] = field(default_factory=dict)


@dataclass
class Chain:
    id: str
    residues: List[Residue] = field(default_factory=list)


@dataclass
class Structure:
    source: Path
    model_number: int
    chains: Dict[str, Chain]


@dataclass(frozen=True)
class ResiduePair:
    first: Residue
    second: Residue


@dataclass(frozen=True)
class ChainMatch:
    first_chain_id: str
    second_chain_id: str
    sequence_score: float
    sequence_identity: float
    residue_pairs: Tuple[ResiduePair, ...]
    unmatched_first: Tuple[ResidueId, ...]
    unmatched_second: Tuple[ResidueId, ...]


@dataclass(frozen=True)
class AlignmentResult:
    rotation: np.ndarray
    translation: np.ndarray
    initial_rmsd: float
    final_rmsd: float
    fit_pair_count: int
    cycles: int
    fit_residue_keys: frozenset[Tuple[ResidueId, ResidueId]]


@dataclass(frozen=True)
class ResidueResult:
    first_id: ResidueId
    first_name: str
    second_id: ResidueId
    second_name: str
    ca_rmsd: Optional[float]
    backbone_rmsd: Optional[float]
    backbone_atom_count: int
    used_for_fitting: bool


@dataclass(frozen=True)
class AnalysisResult:
    first_source: Path
    second_source: Path
    first_model: int
    second_model: int
    alignment_mode: str
    chain_matches: Tuple[ChainMatch, ...]
    unmatched_first_chains: Tuple[str, ...]
    unmatched_second_chains: Tuple[str, ...]
    alignment: AlignmentResult
    residues: Tuple[ResidueResult, ...]
