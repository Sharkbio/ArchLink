#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path


REQUIRED_FILES = [
    "README.md",
    "LICENSE",
    "environment.yml",
    "configuration.yaml",
    "archlink.py",
    "save_models/bacteria_transformer2.pth",
    "save_models/best_random_forest_model_focus0_D_B2.pkl",
    "save_models/feature_columns_focus0_D_B2.pkl",
    "save_models/best_random_forest_model_gas_connect_COMB_A_weight1_A_weight23.pkl",
    "save_models/feature_columns_gas_connect_COMB_A_weight1_A_weight23.pkl",
    "save_models/best_random_forest_model_gas_connect_COMB_C1_cosine_C2_cosine3.pkl",
    "save_models/feature_columns_gas_connect_COMB_C1_cosine_C2_cosine3.pkl",
    "save_models/generateG13",
    "save_models/matching",
    "FragGeneScan-master/run_FragGeneScan.pl",
    "FragGeneScan-master/Makefile",
    "scripts/build_fraggenescan.sh",
    "scripts/clustering_config.py",
    "scripts/checkm_runtime.py",
    "scripts/output_validation.py",
    "scripts/export_archlink_mags.py",
    "scripts/resume_archlink_linking.py",
    "scripts/resume_archlink_linking.sh",
    "docs/TROUBLESHOOTING.md",
    "contrastive_learning/train_CLmodel.py",
    "connect04/random_forest_predict_dir1.py",
    "binning02/random_forest_predict_bd0.py",
    "example/README.md",
    "example/config.minimal.yaml",
    "example/config.toy.yaml",
    "example/create_toy_data.py",
    "example/run_toy.py",
    "example/run_toy.sh",
    "example/run_toy.ps1",
    "example/toy_bam.py",
    "example/data/toy_contigs.fasta",
    "example/data/toy_reads.sorted.bam",
    "example/data/toy_reads.sorted.bam.bai",
    "example/expected/expected_bins.tsv",
    "example/expected/expected_candidates.tsv",
    "example/expected/expected_links.tsv",
    "example/expected/expected_unconnected_ends.tsv",
    "example/expected/expected_scaffolds.fasta",
    "scripts/check_toy_output.py",
    "figures/README.md",
    "figures/make_all_figures.py",
    "figures/fig2_binning_benchmark.py",
    "figures/fig3_link_validation.py",
    "figures/fig4_taxonomy_metabolism.py",
    "figures/fig5_bgc_analysis.py",
    "figures/fig6_quantitative.py",
    "figures/requirements.txt",
    "source_data/Source_Data.xlsx",
    "source_data/Figure_Panel_Data_Source_Map_v3.csv",
    "source_data/SHA256SUMS.txt",
    "save_models/SHA256SUMS.txt",
    "RELEASE_NOTES_v1.0.0-nbt-submission.md",
    "benchmarks/README.md",
]


OPTIONAL_PATHS = []


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check whether the ArchLink repository contains the core files expected for public release."
    )
    parser.add_argument(
        "--root",
        default=Path(__file__).resolve().parents[1],
        type=Path,
        help="Repository root directory. Defaults to the parent of scripts/.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    missing = []

    print(f"Auditing repository root: {root}")

    for rel_path in REQUIRED_FILES:
        candidate = root / rel_path
        if candidate.exists():
            print(f"[OK] {rel_path}")
        else:
            print(f"[MISSING] {rel_path}")
            missing.append(rel_path)

    for rel_path in OPTIONAL_PATHS:
        candidate = root / rel_path
        if candidate.exists():
            print(f"[INFO] Optional path present: {rel_path}")
        else:
            print(f"[INFO] Optional path not present: {rel_path}")

    if missing:
        print("\nRepository audit failed.")
        print("Missing required paths:")
        for rel_path in missing:
            print(f"  - {rel_path}")
        return 1

    print("\nRepository audit passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
