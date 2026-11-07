"""Application service coordinating structure comparison."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np

from .chain_matcher import ChainMatcher
from .errors import AnalysisError
from .geometry import rmsd, transform
from .input_source import StructureInput
from .models import AnalysisResult, ResidueResult
from .pdb_parser import PdbParser
from .sequence import GlobalSequenceAligner
from .structural_alignment import StructuralAligner


BACKBONE_ATOMS = ("N", "CA", "C", "O")
ALIGNMENT_MODES = ("align", "prealigned")


@dataclass(frozen=True)
class AnalysisConfig:
    input: StructureInput
    alignment_mode: str = "align"
    target_residues: Tuple[str, ...] = ()
    top_count: int = 10
    output_csv: Path = Path("output/residue_rmsd.csv")

    def validate(self) -> None:
        if self.alignment_mode not in ALIGNMENT_MODES:
            raise AnalysisError(
                f"alignment_mode must be one of {ALIGNMENT_MODES}, got {self.alignment_mode!r}"
            )
        if self.top_count < 1:
            raise AnalysisError("top_count must be at least 1.")


class RmsdAnalysisService:
    def __init__(
        self,
        parser: PdbParser | None = None,
        chain_matcher: ChainMatcher | None = None,
        structural_aligner: StructuralAligner | None = None,
    ):
        self.parser = parser or PdbParser()
        self.chain_matcher = chain_matcher or ChainMatcher(GlobalSequenceAligner())
        self.structural_aligner = structural_aligner or StructuralAligner()

    def analyze(self, config: AnalysisConfig) -> AnalysisResult:
        config.validate()
        first, second = config.input.load(self.parser)
        chain_matches = self.chain_matcher.match(first, second)
        alignment = self.structural_aligner.align(
            chain_matches, prealigned=config.alignment_mode == "prealigned"
        )
        residue_results = self._calculate_residue_results(chain_matches, alignment)
        matched_first_chains = {match.first_chain_id for match in chain_matches}
        matched_second_chains = {match.second_chain_id for match in chain_matches}
        return AnalysisResult(
            first.source,
            second.source,
            first.model_number,
            second.model_number,
            config.alignment_mode,
            chain_matches,
            tuple(chain for chain in first.chains if chain not in matched_first_chains),
            tuple(chain for chain in second.chains if chain not in matched_second_chains),
            alignment,
            tuple(residue_results),
        )

    @staticmethod
    def _calculate_residue_results(chain_matches, alignment) -> List[ResidueResult]:
        results: List[ResidueResult] = []
        for chain_match in chain_matches:
            for pair in chain_match.residue_pairs:
                first_ca = pair.first.atoms.get("CA")
                second_ca = pair.second.atoms.get("CA")
                ca_value = None
                if first_ca is not None and second_ca is not None:
                    transformed_ca = transform(
                        second_ca.reshape(1, 3), alignment.rotation, alignment.translation
                    )[0]
                    ca_value = float(np.linalg.norm(first_ca - transformed_ca))

                common_backbone = [
                    atom_name
                    for atom_name in BACKBONE_ATOMS
                    if atom_name in pair.first.atoms and atom_name in pair.second.atoms
                ]
                backbone_value = None
                if common_backbone:
                    first_backbone = np.vstack(
                        [pair.first.atoms[name] for name in common_backbone]
                    )
                    second_backbone = np.vstack(
                        [pair.second.atoms[name] for name in common_backbone]
                    )
                    transformed_backbone = transform(
                        second_backbone, alignment.rotation, alignment.translation
                    )
                    backbone_value = rmsd(first_backbone, transformed_backbone)

                fit_key = (pair.first.id, pair.second.id)
                results.append(
                    ResidueResult(
                        pair.first.id,
                        pair.first.name,
                        pair.second.id,
                        pair.second.name,
                        ca_value,
                        backbone_value,
                        len(common_backbone),
                        fit_key in alignment.fit_residue_keys,
                    )
                )
        return results
