"""Global protein sequence alignment with BLOSUM62 and affine gap costs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np


AMINO_ACID_CODES = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "X", "PYL": "X", "ASX": "B", "GLX": "Z", "UNK": "X",
}

_BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1
"""


def _load_blosum62() -> dict[Tuple[str, str], int]:
    rows = [line.split() for line in _BLOSUM62_TEXT.strip().splitlines()]
    columns = rows[0]
    return {
        (row[0], column): int(score)
        for row in rows[1:]
        for column, score in zip(columns, row[1:])
    }


BLOSUM62 = _load_blosum62()


@dataclass(frozen=True)
class SequenceAlignment:
    score: float
    index_pairs: Tuple[Tuple[int, int], ...]
    identity: float


class GlobalSequenceAligner:
    """Needleman-Wunsch alignment using PyMOL align's default penalties."""

    def __init__(self, gap_open: float = -10.0, gap_extend: float = -0.5):
        self.gap_open = gap_open
        self.gap_extend = gap_extend

    def align(self, first: Sequence[str], second: Sequence[str]) -> SequenceAlignment:
        first_seq = [AMINO_ACID_CODES.get(name.upper(), "X") for name in first]
        second_seq = [AMINO_ACID_CODES.get(name.upper(), "X") for name in second]
        n, m = len(first_seq), len(second_seq)
        negative_infinity = -np.inf

        match = np.full((n + 1, m + 1), negative_infinity)
        gap_second = np.full_like(match, negative_infinity)
        gap_first = np.full_like(match, negative_infinity)
        trace_match = np.full((n + 1, m + 1), -1, dtype=np.int8)
        trace_gap_second = np.full_like(trace_match, -1)
        trace_gap_first = np.full_like(trace_match, -1)
        match[0, 0] = 0.0

        for i in range(1, n + 1):
            gap_second[i, 0] = self.gap_open + (i - 1) * self.gap_extend
            trace_gap_second[i, 0] = 0 if i == 1 else 1
        for j in range(1, m + 1):
            gap_first[0, j] = self.gap_open + (j - 1) * self.gap_extend
            trace_gap_first[0, j] = 0 if j == 1 else 2

        for i in range(1, n + 1):
            for j in range(1, m + 1):
                previous = (match[i - 1, j - 1], gap_second[i - 1, j - 1], gap_first[i - 1, j - 1])
                trace_match[i, j] = int(np.argmax(previous))
                match[i, j] = max(previous) + BLOSUM62[(first_seq[i - 1], second_seq[j - 1])]

                gap_second_options = (
                    match[i - 1, j] + self.gap_open,
                    gap_second[i - 1, j] + self.gap_extend,
                    gap_first[i - 1, j] + self.gap_open,
                )
                trace_gap_second[i, j] = int(np.argmax(gap_second_options))
                gap_second[i, j] = max(gap_second_options)

                gap_first_options = (
                    match[i, j - 1] + self.gap_open,
                    gap_second[i, j - 1] + self.gap_open,
                    gap_first[i, j - 1] + self.gap_extend,
                )
                trace_gap_first[i, j] = int(np.argmax(gap_first_options))
                gap_first[i, j] = max(gap_first_options)

        state = int(np.argmax((match[n, m], gap_second[n, m], gap_first[n, m])))
        score = max(match[n, m], gap_second[n, m], gap_first[n, m])
        pairs: List[Tuple[int, int]] = []
        i, j = n, m
        while i > 0 or j > 0:
            if state == 0:
                pairs.append((i - 1, j - 1))
                state = int(trace_match[i, j])
                i -= 1
                j -= 1
            elif state == 1:
                state = int(trace_gap_second[i, j])
                i -= 1
            else:
                state = int(trace_gap_first[i, j])
                j -= 1

        pairs.reverse()
        identical = sum(first_seq[a] == second_seq[b] for a, b in pairs)
        identity = identical / len(pairs) if pairs else 0.0
        return SequenceAlignment(float(score), tuple(pairs), identity)
