"""Validation helpers for detecting incomplete ArchLink runs."""

from __future__ import annotations

import csv
from pathlib import Path


def _nonempty_files(directory: Path, suffixes: tuple[str, ...]) -> list[Path]:
    if not directory.is_dir():
        return []
    return [
        path
        for path in directory.iterdir()
        if path.is_file() and path.suffix.lower() in suffixes and path.stat().st_size > 0
    ]


def _after_link_ids(directory: Path) -> set[str]:
    ids = set()
    for path in _nonempty_files(directory, (".fa", ".fna", ".fasta")):
        if path.stem.endswith(("_c", "_g")):
            ids.add(path.stem[:-2])
    return ids


def _selected_before_ids(report: Path, contamination_max: float) -> set[str]:
    if not report.is_file() or report.stat().st_size == 0:
        raise FileNotFoundError(
            f"CheckM2 quality report is missing or empty: {report}"
        )

    with report.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"Name", "Contamination"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            found = ", ".join(reader.fieldnames or [])
            raise ValueError(
                "CheckM2 quality report must contain Name and Contamination "
                f"columns; found: {found}"
            )

        selected = set()
        for row in reader:
            try:
                contamination = float(row["Contamination"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid CheckM2 contamination value for bin {row.get('Name')!r}: "
                    f"{row.get('Contamination')!r}"
                ) from exc
            if contamination < contamination_max:
                selected.add(str(row["Name"]))
    return selected


def validate_archlink_output(output_path: str | Path) -> dict[str, int]:
    """Validate the before-link and after-link MAG collections.

    The check is intentionally independent of CheckM2 table formatting: it
    verifies that the files consumed by downstream QC actually exist.
    """
    root = Path(output_path).expanduser().resolve()
    before_dir = root / "binning" / "bins"
    after_dir = root / "linking" / "connect"
    selected_ids = _selected_before_ids(
        root / "binning" / "checkm2_bins" / "quality_report.tsv",
        contamination_max=10.0,
    )
    if not selected_ids:
        raise RuntimeError(
            "CheckM2 selected zero bins with contamination < 10%; "
            "the ArchLink run cannot be marked successful."
        )

    before_files = {
        path.stem: path
        for path in _nonempty_files(before_dir, (".fa", ".fna", ".fasta"))
    }
    before_ids = selected_ids.intersection(before_files)
    missing_before = sorted(selected_ids.difference(before_ids))
    if missing_before:
        preview = ", ".join(missing_before[:10])
        suffix = " ..." if len(missing_before) > 10 else ""
        raise RuntimeError(
            "CheckM2-qualified before-link bins are missing from binning/bins: "
            f"{len(missing_before)} ({preview}{suffix})."
        )

    after_ids = _after_link_ids(after_dir)
    missing_after = sorted(before_ids.difference(after_ids))
    extra_after = sorted(after_ids.difference(before_ids))
    if missing_after or extra_after:
        details = []
        if missing_after:
            preview = ", ".join(missing_after[:10])
            suffix = " ..." if len(missing_after) > 10 else ""
            details.append(
                f"missing={len(missing_after)} ({preview}{suffix})"
            )
        if extra_after:
            preview = ", ".join(extra_after[:10])
            suffix = " ..." if len(extra_after) > 10 else ""
            details.append(f"extra={len(extra_after)} ({preview}{suffix})")
        raise RuntimeError(
            "ArchLink produced non-matching paired MAG output: "
            + "; ".join(details)
        )
    return {"before_link": len(before_ids), "after_link": len(before_ids)}


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Validate ArchLink MAG outputs.")
    parser.add_argument("--output-path", required=True)
    args = parser.parse_args()
    counts = validate_archlink_output(args.output_path)
    print(
        f"ArchLink output is valid: before_link={counts['before_link']}, "
        f"after_link={counts['after_link']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
