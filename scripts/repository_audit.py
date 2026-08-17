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
    "contrastive_learning/train_CLmodel.py",
    "connect04/random_forest_predict_dir1.py",
    "binning02/random_forest_predict_bd0.py",
    "example/README.md",
    "example/config.minimal.yaml",
    "benchmarks/README.md",
]


OPTIONAL_PATHS = [
    "figures",
    "source_data",
]


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
