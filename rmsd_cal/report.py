"""Terminal and CSV reporting for RMSD analysis results."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, TextIO

from .models import AnalysisResult, ResidueResult


CSV_FIELDS = (
    "chain_1",
    "residue_number_1",
    "insertion_code_1",
    "residue_name_1",
    "chain_2",
    "residue_number_2",
    "insertion_code_2",
    "residue_name_2",
    "ca_rmsd_angstrom",
    "backbone_rmsd_angstrom",
    "backbone_atom_count",
    "used_for_fitting",
)


class CsvResultWriter:
    def write(self, result: AnalysisResult, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=CSV_FIELDS)
            writer.writeheader()
            for residue in result.residues:
                writer.writerow(self._row(residue))

    @staticmethod
    def _row(residue: ResidueResult) -> dict[str, object]:
        return {
            "chain_1": residue.first_id.chain_id,
            "residue_number_1": residue.first_id.number,
            "insertion_code_1": residue.first_id.insertion_code,
            "residue_name_1": residue.first_name,
            "chain_2": residue.second_id.chain_id,
            "residue_number_2": residue.second_id.number,
            "insertion_code_2": residue.second_id.insertion_code,
            "residue_name_2": residue.second_name,
            "ca_rmsd_angstrom": _format_optional(residue.ca_rmsd, 6),
            "backbone_rmsd_angstrom": _format_optional(residue.backbone_rmsd, 6),
            "backbone_atom_count": residue.backbone_atom_count,
            "used_for_fitting": str(residue.used_for_fitting).lower(),
        }


class TerminalReporter:
    def print(
        self,
        result: AnalysisResult,
        top_count: int,
        target_residues: Sequence[str],
        output_csv: Path,
        stream: TextIO,
    ) -> None:
        print("RMSD analysis", file=stream)
        print(
            f"Structures: {result.first_source} (model {result.first_model}) vs "
            f"{result.second_source} (model {result.second_model})",
            file=stream,
        )
        print("\nChain mapping:", file=stream)
        for match in result.chain_matches:
            print(
                f"  {match.first_chain_id} -> {match.second_chain_id}  "
                f"identity={match.sequence_identity * 100:.1f}%  "
                f"aligned_residues={len(match.residue_pairs)}",
                file=stream,
            )
        self._print_warnings(result, stream)

        alignment = result.alignment
        print("\nC-alpha alignment:", file=stream)
        print(f"  Mode: {result.alignment_mode}", file=stream)
        if result.alignment_mode == "prealigned":
            print(f"  Current-coordinate RMSD: {alignment.final_rmsd:.4f} A", file=stream)
            print(f"  Matched C-alpha pairs: {alignment.fit_pair_count}", file=stream)
        else:
            print(f"  RMSD before refinement: {alignment.initial_rmsd:.4f} A", file=stream)
            print(f"  RMSD after refinement:  {alignment.final_rmsd:.4f} A", file=stream)
            print(f"  Retained C-alpha pairs: {alignment.fit_pair_count}", file=stream)
            print(f"  Refinement cycles: {alignment.cycles}", file=stream)

        ranked = sorted(
            (residue for residue in result.residues if residue.backbone_rmsd is not None),
            key=lambda residue: (
                residue.backbone_rmsd,
                float("inf") if residue.ca_rmsd is None else residue.ca_rmsd,
            ),
        )[:top_count]
        print(f"\nLowest backbone RMSD residues (top {len(ranked)}):", file=stream)
        self._print_residue_table(ranked, stream)

        if target_residues:
            print("\nSelected residues (structure 1 numbering):", file=stream)
            for target in target_residues:
                matches = [
                    residue for residue in result.residues
                    if residue.first_id.label().upper() == target.strip().upper()
                ]
                if matches:
                    self._print_residue_table(matches, stream, include_header=False)
                else:
                    print(f"  {target}: not found among aligned residues", file=stream)
        print(f"\nCSV written to: {output_csv}", file=stream)

    @staticmethod
    def _print_warnings(result: AnalysisResult, stream: TextIO) -> None:
        if result.unmatched_first_chains:
            print(f"  Warning: unmatched chains in structure 1: {', '.join(result.unmatched_first_chains)}", file=stream)
        if result.unmatched_second_chains:
            print(f"  Warning: unmatched chains in structure 2: {', '.join(result.unmatched_second_chains)}", file=stream)
        for match in result.chain_matches:
            if match.unmatched_first or match.unmatched_second:
                print(
                    f"  Warning: {match.first_chain_id}->{match.second_chain_id} has "
                    f"{len(match.unmatched_first)} and {len(match.unmatched_second)} "
                    "unaligned residues, respectively",
                    file=stream,
                )

    @staticmethod
    def _print_residue_table(
        residues: Iterable[ResidueResult], stream: TextIO, include_header: bool = True
    ) -> None:
        if include_header:
            print("  RESIDUE  MATCHED_TO  NAME  CA_RMSD(A)  BACKBONE_RMSD(A)  BB_ATOMS", file=stream)
        for residue in residues:
            print(
                f"  {residue.first_id.label():<8} "
                f"{residue.second_id.label():<11} "
                f"{residue.first_name:<5} "
                f"{_format_optional(residue.ca_rmsd, 4):>10}  "
                f"{_format_optional(residue.backbone_rmsd, 4):>16}  "
                f"{residue.backbone_atom_count:>8}",
                file=stream,
            )


def _format_optional(value: Optional[float], digits: int) -> str:
    return "" if value is None else f"{value:.{digits}f}"
