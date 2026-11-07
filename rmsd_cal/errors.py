"""User-facing errors for invalid structures and analysis settings."""


class RmsdCalError(Exception):
    """Base class for expected analysis failures."""


class PdbFormatError(RmsdCalError):
    """Raised when a PDB file cannot provide the requested structure."""


class AnalysisError(RmsdCalError):
    """Raised when structures cannot be compared safely."""
