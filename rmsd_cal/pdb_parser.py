"""Minimal, deterministic parser for protein ATOM records in PDB files."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

from .errors import PdbFormatError
from .models import Chain, Residue, ResidueId, Structure


class PdbParser:
    """Parse one model from a PDB file.

    Only protein ``ATOM`` records are retained. For alternate locations, a
    blank location is preferred over ``A``; other locations are ignored.
    """

    def parse(self, path: Path, model_number: int) -> Structure:
        if model_number < 1:
            raise PdbFormatError("Model numbers must be positive integers.")
        if not path.is_file():
            raise PdbFormatError(f"PDB file does not exist: {path}")

        chains: Dict[str, Chain] = OrderedDict()
        residue_lookup: Dict[Tuple[str, int, str], Residue] = {}
        atom_altloc: Dict[Tuple[str, int, str, str], str] = {}
        saw_model_record = False
        active_model = 1
        selected_model_seen = False

        try:
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        except OSError as exc:
            raise PdbFormatError(f"Could not read PDB file: {path}") from exc

        for line_number, line in enumerate(lines, start=1):
            record = line[0:6].strip() if len(line) >= 6 else ""
            if record == "MODEL":
                saw_model_record = True
                try:
                    active_model = int(line[10:14].strip())
                except ValueError as exc:
                    raise PdbFormatError(
                        f"Invalid MODEL record at {path}:{line_number}"
                    ) from exc
                if active_model == model_number:
                    selected_model_seen = True
                continue
            if record == "ENDMDL":
                if saw_model_record and active_model == model_number:
                    break
                continue
            if record != "ATOM":
                continue
            if active_model != model_number:
                continue
            if not saw_model_record and model_number == 1:
                selected_model_seen = True
            if len(line) < 54:
                raise PdbFormatError(f"Short ATOM record at {path}:{line_number}")

            altloc = line[16].strip()
            if altloc not in ("", "A"):
                continue
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip().upper()
            chain_id = line[21].strip() or "_"
            insertion_code = line[26].strip()
            try:
                residue_number = int(line[22:26].strip())
                coordinates = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
            except ValueError as exc:
                raise PdbFormatError(
                    f"Invalid coordinate or residue number at {path}:{line_number}"
                ) from exc

            residue_key = (chain_id, residue_number, insertion_code)
            residue = residue_lookup.get(residue_key)
            if residue is None:
                residue = Residue(
                    id=ResidueId(chain_id, residue_number, insertion_code),
                    name=residue_name,
                )
                residue_lookup[residue_key] = residue
                chains.setdefault(chain_id, Chain(chain_id)).residues.append(residue)
            elif residue.name != residue_name:
                raise PdbFormatError(
                    f"Conflicting residue names for {residue.id.label()} in {path}"
                )

            atom_key = (*residue_key, atom_name)
            previous_altloc = atom_altloc.get(atom_key)
            if previous_altloc == "":
                continue
            if previous_altloc == "A" and altloc == "":
                residue.atoms[atom_name] = coordinates
                atom_altloc[atom_key] = altloc
            elif previous_altloc is None:
                residue.atoms[atom_name] = coordinates
                atom_altloc[atom_key] = altloc

        if not selected_model_seen:
            raise PdbFormatError(f"Model {model_number} was not found in {path}")
        if not chains:
            raise PdbFormatError(
                f"Model {model_number} in {path} contains no usable ATOM records"
            )
        return Structure(path, model_number, chains)
