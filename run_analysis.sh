#!/usr/bin/env bash
set -euo pipefail

# =============================================================================
# USER SETTINGS
# =============================================================================

# Choose one: "two_models" or "two_files"
INPUT_MODE="two_models"

# -----------------------------------------------------------------------------
# Input option 1: compare two MODEL records in one PDB file
# Used when INPUT_MODE="two_models"
# -----------------------------------------------------------------------------
PDB_FILE_1="data/CXCR1_align.pdb"
MODEL_1=1
MODEL_2=2

# -----------------------------------------------------------------------------
# Input option 2: compare two PDB files
# Uncomment PDB_FILE_2 and set INPUT_MODE="two_files"
# MODEL_1 and MODEL_2 above select the model from each file.
# -----------------------------------------------------------------------------
# PDB_FILE_1="data/structure_1.pdb"
# PDB_FILE_2="data/structure_2.pdb"

# "align": sequence alignment plus iterative C-alpha structural fitting
# "prealigned": preserve coordinates and calculate RMSD without fitting
ALIGNMENT_MODE="prealigned"

# Residues use structure 1 numbering: "CHAIN:NUMBER" or "CHAIN:NUMBER+ICODE"
TARGET_RESIDUES=(
  # "A:42"
  # "A:42A"
)

TOP_COUNT=10
OUTPUT_CSV="output/residue_rmsd.csv"

# =============================================================================
# END USER SETTINGS
# =============================================================================

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ "$INPUT_MODE" == "two_files" && -z "${PDB_FILE_2:-}" ]]; then
  echo "Error: PDB_FILE_2 is required when INPUT_MODE=\"two_files\"." >&2
  exit 1
fi

TARGET_RESIDUES_CSV=""
if (( ${#TARGET_RESIDUES[@]} > 0 )); then
  printf -v TARGET_RESIDUES_CSV '%s,' "${TARGET_RESIDUES[@]}"
  TARGET_RESIDUES_CSV="${TARGET_RESIDUES_CSV%,}"
fi

export INPUT_MODE PDB_FILE_1 MODEL_1 MODEL_2 ALIGNMENT_MODE
export TARGET_RESIDUES_CSV TOP_COUNT OUTPUT_CSV
if [[ -n "${PDB_FILE_2:-}" ]]; then
  export PDB_FILE_2
fi

python3 run_analysis.py
