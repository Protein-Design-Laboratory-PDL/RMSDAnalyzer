"""Automatic chain and residue correspondence by sequence and structure."""

from __future__ import annotations

from itertools import permutations
from math import perm
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from .errors import AnalysisError
from .geometry import kabsch_fit, rmsd, transform
from .models import Chain, ChainMatch, ResiduePair, Structure
from .sequence import GlobalSequenceAligner, SequenceAlignment


MAX_EXHAUSTIVE_ASSIGNMENTS = 100_000


class ChainMatcher:
    """Find a one-to-one chain mapping with a structural tie-breaker.

    Sequence score is the primary criterion. If multiple mappings have the
    same best score, the mapping with the lowest combined C-alpha fitting RMSD
    is selected.
    """

    def __init__(self, sequence_aligner: GlobalSequenceAligner):
        self.sequence_aligner = sequence_aligner

    def match(self, first: Structure, second: Structure) -> Tuple[ChainMatch, ...]:
        first_chains = list(first.chains.values())
        second_chains = list(second.chains.values())
        pair_alignments = self._align_all_pairs(first_chains, second_chains)
        assignment_count = perm(max(len(first_chains), len(second_chains)), min(len(first_chains), len(second_chains)))
        if assignment_count > MAX_EXHAUSTIVE_ASSIGNMENTS:
            raise AnalysisError(
                "Automatic chain mapping would require evaluating "
                f"{assignment_count:,} assignments. Reduce the input to the "
                "biological assembly of interest."
            )

        assignments = self._assignments(len(first_chains), len(second_chains))
        best_score = -np.inf
        best_structural_rmsd = np.inf
        best_assignment: Tuple[Tuple[int, int], ...] | None = None
        for assignment in assignments:
            score = sum(pair_alignments[pair].score for pair in assignment)
            if score < best_score:
                continue
            structural_rmsd = self._combined_ca_rmsd(
                first_chains, second_chains, assignment, pair_alignments
            )
            if score > best_score or structural_rmsd < best_structural_rmsd:
                best_score = score
                best_structural_rmsd = structural_rmsd
                best_assignment = assignment

        if best_assignment is None:
            raise AnalysisError("No chain mapping could be determined.")

        matches: List[ChainMatch] = []
        for first_index, second_index in best_assignment:
            first_chain = first_chains[first_index]
            second_chain = second_chains[second_index]
            alignment = pair_alignments[(first_index, second_index)]
            residue_pairs = tuple(
                ResiduePair(first_chain.residues[a], second_chain.residues[b])
                for a, b in alignment.index_pairs
            )
            matched_first_indices = {a for a, _ in alignment.index_pairs}
            matched_second_indices = {b for _, b in alignment.index_pairs}
            matches.append(
                ChainMatch(
                    first_chain.id,
                    second_chain.id,
                    alignment.score,
                    alignment.identity,
                    residue_pairs,
                    tuple(
                        residue.id for index, residue in enumerate(first_chain.residues)
                        if index not in matched_first_indices
                    ),
                    tuple(
                        residue.id for index, residue in enumerate(second_chain.residues)
                        if index not in matched_second_indices
                    ),
                )
            )
        return tuple(sorted(matches, key=lambda match: match.first_chain_id))

    def _align_all_pairs(
        self, first_chains: Sequence[Chain], second_chains: Sequence[Chain]
    ) -> Dict[Tuple[int, int], SequenceAlignment]:
        return {
            (first_index, second_index): self.sequence_aligner.align(
                [residue.name for residue in first_chain.residues],
                [residue.name for residue in second_chain.residues],
            )
            for first_index, first_chain in enumerate(first_chains)
            for second_index, second_chain in enumerate(second_chains)
        }

    @staticmethod
    def _assignments(first_count: int, second_count: int) -> Iterable[Tuple[Tuple[int, int], ...]]:
        if first_count <= second_count:
            for second_indices in permutations(range(second_count), first_count):
                yield tuple(enumerate(second_indices))
        else:
            for first_indices in permutations(range(first_count), second_count):
                yield tuple((first_index, second_index) for second_index, first_index in enumerate(first_indices))

    @staticmethod
    def _combined_ca_rmsd(
        first_chains: Sequence[Chain],
        second_chains: Sequence[Chain],
        assignment: Tuple[Tuple[int, int], ...],
        pair_alignments: Dict[Tuple[int, int], SequenceAlignment],
    ) -> float:
        first_coordinates: List[np.ndarray] = []
        second_coordinates: List[np.ndarray] = []
        for first_index, second_index in assignment:
            first_chain = first_chains[first_index]
            second_chain = second_chains[second_index]
            for first_residue_index, second_residue_index in pair_alignments[(first_index, second_index)].index_pairs:
                first_ca = first_chain.residues[first_residue_index].atoms.get("CA")
                second_ca = second_chain.residues[second_residue_index].atoms.get("CA")
                if first_ca is not None and second_ca is not None:
                    first_coordinates.append(first_ca)
                    second_coordinates.append(second_ca)
        if len(first_coordinates) < 3:
            return np.inf
        first_array = np.vstack(first_coordinates)
        second_array = np.vstack(second_coordinates)
        rotation, translation = kabsch_fit(first_array, second_array)
        return rmsd(first_array, transform(second_array, rotation, translation))
