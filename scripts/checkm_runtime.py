"""Runtime checks and configuration helpers for the CheckM1 dependency."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


def _candidate_database_files(root: Path) -> list[Path]:
    return [
        root / "pfam" / "Pfam-A.hmm.dat",
        root / "Pfam-A.hmm.dat",
    ]


def resolve_checkm1_data_path(configured: str | None = None) -> Path | None:
    """Return a configured CheckM1 database root, if one is available."""
    placeholders = {
        "",
        "/path/to/checkm1_database",
        "/path/to/checkm1_data",
    }
    values = [
        configured,
        os.environ.get("CHECKM_DATA_PATH"),
        os.environ.get("CHECKM_DATA_ROOT"),
    ]
    for value in values:
        if value and str(value).strip() not in placeholders:
            return Path(os.path.expanduser(str(value))).resolve()
    return None


def ensure_checkm1_ready(configured: str | None = None) -> Path:
    """Validate CheckM1 data and configure CheckM's data root when possible.

    ArchLink uses CheckM1 through the bundled UniItem code to select the best
    Leiden result. CheckM2 is a separate downstream quality assessment.
    """
    data_root = resolve_checkm1_data_path(configured)
    if data_root is None:
        raise RuntimeError(
            "CheckM1 database path is not configured. Set "
            "common.path.checkm1_data_path in the YAML file or export "
            "CHECKM_DATA_PATH before running ArchLink."
        )

    expected = next((path for path in _candidate_database_files(data_root) if path.is_file()), None)
    if expected is None:
        checked = ", ".join(str(path) for path in _candidate_database_files(data_root))
        raise FileNotFoundError(
            "CheckM1 database is incomplete: Pfam-A.hmm.dat was not found. "
            f"Checked: {checked}. Download/initialize CheckM1 at {data_root}."
        )

    os.environ["CHECKM_DATA_PATH"] = str(data_root)
    os.environ["CHECKM_DATA_ROOT"] = str(data_root)

    checkm = shutil.which("checkm")
    if checkm is None:
        raise RuntimeError(
            "The CheckM1 executable `checkm` was not found on PATH. "
            "Activate the ArchLink environment before running."
        )

    # CheckM1 stores its data root in its own configuration. Synchronize it
    # with the explicit ArchLink setting so it cannot silently fall back to
    # ~/.checkm.
    result = subprocess.run(
        [checkm, "data", "setRoot", str(data_root)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Unable to configure CheckM1 data root. "
            f"Command: {checkm} data setRoot {data_root}\n"
            f"stderr: {result.stderr.strip()}"
        )
    return data_root
