"""Rigid-body fitting and RMSD primitives."""

from __future__ import annotations

from typing import Tuple

import numpy as np

from .errors import AnalysisError


def kabsch_fit(reference: np.ndarray, mobile: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return a rotation and translation mapping ``mobile`` to ``reference``."""
    if reference.shape != mobile.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise AnalysisError("Kabsch fitting requires equally shaped N x 3 coordinates.")
    if len(reference) < 3:
        raise AnalysisError("At least three matched C-alpha atoms are required for fitting.")

    reference_center = reference.mean(axis=0)
    mobile_center = mobile.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    left, _, right_transposed = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right_transposed))
    rotation = left @ correction @ right_transposed
    translation = reference_center - mobile_center @ rotation
    return rotation, translation


def transform(coordinates: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return coordinates @ rotation + translation


def rmsd(first: np.ndarray, second: np.ndarray) -> float:
    if first.shape != second.shape or not len(first):
        raise AnalysisError("RMSD requires non-empty, equally shaped coordinates.")
    squared_distances = np.sum((first - second) ** 2, axis=1)
    return float(np.sqrt(np.mean(squared_distances)))
