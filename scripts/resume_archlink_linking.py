#!/usr/bin/env python3
"""Resume ArchLink after contrastive learning and primary clustering.

This is a real Python module rather than a shell here-document. That matters
when downstream stages use multiprocessing with the ``spawn`` start method.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Make the repository root importable when this file is launched by path.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from archlink import Args
from binning02 import binning_main
from connect04 import connect_main
from generate_cos03 import generate_cos_main
from scripts.checkm_runtime import ensure_checkm1_ready
from scripts.clustering_config import normalize_clustering_mode
from scripts.output_validation import validate_archlink_output


def _has_files(directory: Path, suffixes: tuple[str, ...]) -> bool:
    return directory.is_dir() and any(
        path.is_file() and path.suffix.lower() in suffixes and path.stat().st_size > 0
        for path in directory.iterdir()
    )


def _has_complete_linking_outputs(input_dir: Path, output_dir: Path) -> bool:
    """Return whether every linking input bin has a corresponding output FASTA."""
    if not input_dir.is_dir() or not output_dir.is_dir():
        return False

    input_ids = {
        path.name
        for path in input_dir.iterdir()
        if path.is_dir()
    }
    if not input_ids:
        return False

    output_ids = set()
    for path in output_dir.iterdir():
        if not path.is_file() or path.suffix.lower() not in {".fa", ".fna", ".fasta"}:
            continue
        if path.stem.endswith("_c"):
            output_ids.add(path.stem[:-2])
        elif path.stem.endswith("_g"):
            output_ids.add(path.stem[:-2])
    return input_ids == output_ids


def _configure_logger(output_path: Path) -> logging.Logger:
    output_path.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("ArchLinkResume")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)
    file_handler = logging.FileHandler(output_path / "ArchLink.resume.log")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger


def _load_args(config_path: Path, clustering_mode: str | None):
    args = Args.from_yaml(str(config_path))
    args.clustering_mode = normalize_clustering_mode(
        clustering_mode or getattr(args, "clustering_mode", "full")
    )
    args.checkm1_data_path = ensure_checkm1_ready(
        getattr(args, "checkm1_data_path", None)
    )
    return args


def resume(args, logger, start_at: str = "auto") -> dict[str, int]:
    output_path = Path(args.output_path).expanduser().resolve()
    estimate = output_path / "cluster_res" / "estimate_res.txt"
    if not estimate.is_file() or estimate.stat().st_size == 0:
        raise FileNotFoundError(
            f"Primary clustering summary is missing: {estimate}. "
            "Run the full ArchLink pipeline before using this resume script."
        )

    checkm2_report = output_path / "binning" / "checkm2_bins" / "quality_report.tsv"
    binning_bins_dir = output_path / "binning" / "bins"
    cosine_file = output_path / "linking" / "cosine" / "cosine_model_features_softmax.pkl"
    linking_bins_dir = output_path / "linking" / "bins_0.9"
    connect_dir = output_path / "linking" / "connect"

    run_binning = start_at == "binning" or (
        start_at == "auto"
        and not (
            checkm2_report.is_file()
            and checkm2_report.stat().st_size > 0
            and _has_files(binning_bins_dir, (".fa", ".fna", ".fasta"))
        )
    )
    if start_at == "binning" or run_binning:
        logger.info("Resume stage: secondary binning -> CheckM2")
        binning_main.binning_init(args, logger)
    else:
        logger.info("Secondary binning output is present; skipping binning stage.")

    if not _has_files(binning_bins_dir, (".fa", ".fna", ".fasta")):
        raise RuntimeError(
            f"Secondary bin FASTA files are missing or empty: {binning_bins_dir}. "
            "Check the CheckM2 stage and its quality report."
        )
    graph_file = output_path / "bam.graph"
    if not graph_file.is_file() or graph_file.stat().st_size == 0:
        raise FileNotFoundError(f"Required BAM graph is missing or empty: {graph_file}")

    run_cosine = start_at in ("auto", "binning", "cosine") and not (
        cosine_file.is_file()
        and cosine_file.stat().st_size > 0
        and _has_files(linking_bins_dir, (".fasta", ".fa", ".fna"))
    )
    if start_at == "cosine" or run_cosine:
        logger.info("Resume stage: context-aware cosine features")
        generate_cos_main.main(args)
    else:
        logger.info("Cosine feature output is present; skipping cosine stage.")

    if not cosine_file.is_file() or cosine_file.stat().st_size == 0:
        raise RuntimeError(f"Cosine feature output is missing or empty: {cosine_file}")
    if not _has_files(linking_bins_dir, (".fasta", ".fa", ".fna")):
        raise RuntimeError(
            f"Linking input bins are missing or empty: {linking_bins_dir}"
        )

    run_connect = start_at in ("auto", "binning", "cosine", "connect") and not (
        _has_complete_linking_outputs(linking_bins_dir, connect_dir)
    )
    if start_at == "connect" or run_connect:
        logger.info("Resume stage: RF linking and scaffolding")
        connect_main.connect_main(args)
    else:
        logger.info("Connected MAG output is present; skipping connect stage.")

    counts = validate_archlink_output(output_path)
    logger.info(
        "Resume completed: before_link=%d, after_link=%d",
        counts["before_link"],
        counts["after_link"],
    )
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Resume ArchLink from secondary binning or a later linking stage."
    )
    parser.add_argument("--config", default="configuration.yaml")
    parser.add_argument(
        "--start-at",
        choices=("auto", "binning", "cosine", "connect"),
        default="auto",
        help="Stage to run. auto skips complete stages when their outputs are valid.",
    )
    parser.add_argument(
        "--clustering-mode",
        choices=("full", "fast"),
        default=None,
        help="Override the mode in the configuration file.",
    )
    cli_args = parser.parse_args()
    config_path = Path(cli_args.config).expanduser().resolve()
    args = _load_args(config_path, cli_args.clustering_mode)
    logger = _configure_logger(Path(args.output_path))
    logger.info("Configuration: %s", config_path)
    logger.info("Output path: %s", args.output_path)
    resume(args, logger, cli_args.start_at)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
