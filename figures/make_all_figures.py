#!/usr/bin/env python3
"""Regenerate all data-driven ArchLink figure panels from Source_Data.xlsx."""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
import numpy
import pandas

import fig2_binning_benchmark
import fig3_link_validation
import fig4_taxonomy_metabolism
import fig5_bgc_analysis
import fig6_quantitative
from common import DEFAULT_DPI, DEFAULT_OUTPUT, DEFAULT_SEED, DEFAULT_SOURCE_DATA, sha256


GENERATORS = {
    "2": fig2_binning_benchmark.generate,
    "3": fig3_link_validation.generate,
    "4": fig4_taxonomy_metabolism.generate,
    "5": fig5_bgc_analysis.generate,
    "6": fig6_quantitative.generate,
}
SHEETS = {
    "2": ["Fig2a_ARI_ACC", "Fig2b_CAMI_bins", "Fig2c_real_MAGs"],
    "3": ["Fig3a_precision_TP", "Fig3b_longread_main", "Fig3c_N50_raw", "Fig3c_N50_summary"],
    "4": ["Fig4a_tax_rank_counts", "Fig4b_align_frac_raw", "Fig4c_marker_summary", "Fig4d_DRAM_raw"],
    "5": ["Fig5a_BGC_counts", "Fig5b_BGC_lengths", "Fig5c_BGC_class_resolved"],
    "6": ["Fig6a_scaffold_summary", "Fig6d_PAE_matrix", "Fig6e_depth_profile"],
}


def git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-data", type=Path, default=DEFAULT_SOURCE_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figures", nargs="+", choices=tuple(GENERATORS), default=tuple(GENERATORS))
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--formats", nargs="+", choices=("pdf", "png", "svg"), default=("pdf",))
    args = parser.parse_args()
    source_data = args.source_data.resolve()
    output = args.output.resolve()
    generated = []
    for figure in args.figures:
        generated.extend(
            GENERATORS[figure](source_data, output, args.dpi, args.seed, args.formats)
        )

    root = Path(__file__).resolve().parents[1]
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository_commit": git_commit(root),
        "source_data": str(source_data),
        "source_data_sha256": sha256(source_data),
        "figures": list(args.figures),
        "sheets": {figure: SHEETS[figure] for figure in args.figures},
        "random_seed": args.seed,
        "dpi": args.dpi,
        "formats": list(args.formats),
        "python": platform.python_version(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "matplotlib": matplotlib.__version__,
        "openpyxl": importlib.metadata.version("openpyxl"),
        "outputs": [str(path) for path in generated],
    }
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "figure_build_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    for path in generated:
        print(path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
