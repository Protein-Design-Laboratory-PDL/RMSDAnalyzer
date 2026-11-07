"""PyMOL-align-inspired iterative C-alpha structural fitting."""

from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

from .errors import AnalysisError
from .geometry import kabsch_fit, rmsd, transform
from .models import AlignmentResult, ChainMatch, ResidueId, ResiduePair


class StructuralAligner:
    def __init__(self, cutoff: float = 2.0, max_cycles: int = 5):
        if cutoff <= 0 or max_cycles < 0:
            raise ValueError("Alignment cutoff must be positive and cycles non-negative.")
        self.cutoff = cutoff
        self.max_cycles = max_cycles

    def align(self, matches: Sequence[ChainMatch], prealigned: bool) -> AlignmentResult:
        ca_pairs = [
            pair
            for match in matches
            for pair in match.residue_pairs
            if "CA" in pair.first.atoms and "CA" in pair.second.atoms
        ]
        if not ca_pairs:
            raise AnalysisError("No matched C-alpha atoms were found.")
        reference = np.vstack([pair.first.atoms["CA"] for pair in ca_pairs])
        mobile = np.vstack([pair.second.atoms["CA"] for pair in ca_pairs])

        if prealigned:
            rotation = np.eye(3)
            translation = np.zeros(3)
            current_rmsd = rmsd(reference, mobile)
            return AlignmentResult(
                rotation, translation, current_rmsd, current_rmsd,
                len(ca_pairs), 0, frozenset(),
            )

        if len(ca_pairs) < 3:
            raise AnalysisError("At least three matched C-alpha atoms are required to align structures.")
        active = np.arange(len(ca_pairs))
        rotation, translation = kabsch_fit(reference, mobile)
        initial_rmsd = rmsd(reference, transform(mobile, rotation, translation))
        cycles = 0

        for _ in range(self.max_cycles):
            transformed = transform(mobile, rotation, translation)
            distances = np.linalg.norm(reference - transformed, axis=1)
            retained = active[distances[active] <= self.cutoff]
            if len(retained) < 3:
                raise AnalysisError(
                    "Outlier rejection left fewer than three C-alpha pairs; "
                    "the structures cannot be fitted reliably."
                )
            if np.array_equal(retained, active):
                break
            active = retained
            rotation, translation = kabsch_fit(reference[active], mobile[active])
            cycles += 1

        final_mobile = transform(mobile[active], rotation, translation)
        final_rmsd = rmsd(reference[active], final_mobile)
        return AlignmentResult(
            rotation,
            translation,
            initial_rmsd,
            final_rmsd,
            len(active),
            cycles,
            self._fit_keys(ca_pairs, active),
        )

    @staticmethod
    def _fit_keys(
        pairs: Sequence[ResiduePair], indices: np.ndarray
    ) -> frozenset[Tuple[ResidueId, ResidueId]]:
        return frozenset((pairs[index].first.id, pairs[index].second.id) for index in indices)
