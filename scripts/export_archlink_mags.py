#!/usr/bin/env python3
"""Export paired ArchLink MAG sets for downstream quality control."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd


def _fasta_files(directory: Path) -> dict[str, Path]:
    if not directory.is_dir():
        return {}
    files = {}
    for path in directory.iterdir():
        if (
            path.is_file()
            and path.suffix.lower() in {".fa", ".fna", ".fasta"}
            and path.stat().st_size > 0
        ):
            files[path.stem] = path
    return files


def _read_high_purity_ids(report: Path, contamination_max: float) -> set[str]:
    if not report.is_file():
        raise FileNotFoundError(f"CheckM2 quality report not found: {report}")
    table = pd.read_csv(report, sep="\t")
    required = {"Name", "Contamination"}
    missing = required.difference(table.columns)
    if missing:
        raise ValueError(
            f"CheckM2 report is missing columns: {', '.join(sorted(missing))}"
        )
    return {
        str(row["Name"])
        for _, row in table.iterrows()
        if float(row["Contamination"]) < contamination_max
    }


def _copy_selected(files: dict[str, Path], selected: set[str], destination: Path) -> int:
    if destination.is_dir():
        for old_path in destination.iterdir():
            if old_path.is_dir():
                shutil.rmtree(old_path)
            else:
                old_path.unlink()
    destination.mkdir(parents=True, exist_ok=True)
    copied = 0
    for bin_id in sorted(selected):
        source = files.get(bin_id)
        if source is None:
            continue
        shutil.copy2(source, destination / f"{bin_id}.fa")
        copied += 1
    return copied


def export_mags(
    output_path: str | Path,
    destination: str | Path | None = None,
    contamination_max: float = 10.0,
) -> dict[str, int]:
    root = Path(output_path).expanduser().resolve()
    export_root = (
        Path(destination).expanduser().resolve()
        if destination
        else root / "mags" / "export"
    )
    selected = _read_high_purity_ids(
        root / "binning" / "checkm2_bins" / "quality_report.tsv",
        contamination_max,
    )
    before_files = _fasta_files(root / "binning" / "bins")
    after_raw = _fasta_files(root / "linking" / "connect")

    # Linking writes <bin>_c.fa for scaffolded output and <bin>_g.fa as a
    # fallback. Prefer scaffolded files and retain only matching bin IDs.
    after_files: dict[str, Path] = {}
    for stem, path in after_raw.items():
        if stem.endswith("_c"):
            after_files[stem[:-2]] = path
    for stem, path in after_raw.items():
        if stem.endswith("_g") and stem[:-2] not in after_files:
            after_files[stem[:-2]] = path

    before_count = _copy_selected(
        before_files, selected, export_root / "before_link"
    )
    after_count = _copy_selected(
        after_files, selected, export_root / "after_link"
    )
    missing_after = sorted(selected.difference(after_files))
    missing_before = sorted(selected.difference(before_files))

    export_root.mkdir(parents=True, exist_ok=True)
    with open(export_root / "EXPORT.txt", "w", encoding="utf-8") as handle:
        handle.write("before_link source: {output}/binning/bins\n")
        handle.write(
            "before_link filter: CheckM2 Contamination < "
            f"{contamination_max:g}\n"
        )
        handle.write("after_link source: {output}/linking/connect\n")
        handle.write("after_link matching: same bin IDs as before_link\n")
        handle.write(f"before_link_count: {before_count}\n")
        handle.write(f"after_link_count: {after_count}\n")
        if missing_after:
            handle.write("missing_after_link_ids: " + ",".join(missing_after) + "\n")
        if missing_before:
            handle.write("missing_before_link_ids: " + ",".join(missing_before) + "\n")

    if before_count == 0:
        raise RuntimeError(
            "No before_link MAGs were exported. Check the CheckM2 report and "
            f"contamination threshold < {contamination_max:g}."
        )
    if missing_before:
        raise RuntimeError(
            "CheckM2-qualified bins are missing from binning/bins: "
            + ", ".join(missing_before)
        )
    if missing_after:
        raise RuntimeError(
            "CheckM2-qualified bin IDs are missing from linking/connect: "
            + ", ".join(missing_after)
        )
    if after_count == 0:
        raise RuntimeError(
            "No after_link MAGs were exported. Check linking/connect output."
        )
    return {"before_link": before_count, "after_link": after_count}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-path", required=True)
    parser.add_argument("--destination", default=None)
    parser.add_argument("--contamination-max", type=float, default=10.0)
    args = parser.parse_args()
    counts = export_mags(
        args.output_path,
        args.destination,
        args.contamination_max,
    )
    print(
        f"Exported before_link={counts['before_link']} and "
        f"after_link={counts['after_link']} MAGs."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
