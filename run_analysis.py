"""Python entry point invoked by run_analysis.sh."""

import sys

from rmsd_cal import RmsdAnalysisService
from rmsd_cal.environment_config import EnvironmentConfigLoader
from rmsd_cal.errors import RmsdCalError
from rmsd_cal.report import CsvResultWriter, TerminalReporter


def main() -> int:
    try:
        config = EnvironmentConfigLoader().load()
        result = RmsdAnalysisService().analyze(config)
        CsvResultWriter().write(result, config.output_csv)
        TerminalReporter().print(
            result,
            config.top_count,
            config.target_residues,
            config.output_csv,
            sys.stdout,
        )
    except RmsdCalError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
