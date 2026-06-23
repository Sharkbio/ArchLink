#!/usr/bin/env python3

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path


def reverse_direction(direction: str) -> str:
    return "+" if direction == "-" else "-"


def parse_ground_truth(path: Path):
    directional = set()
    agnostic = set()

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) != 2:
                continue

            contig1_full, contig2_full = parts
            contig1_id, contig1_dir = contig1_full[:-1], contig1_full[-1]
            contig2_id, contig2_dir = contig2_full[:-1], contig2_full[-1]

            pair = ((contig1_id, contig1_dir), (contig2_id, contig2_dir))
            reverse_pair = (
                (contig2_id, reverse_direction(contig2_dir)),
                (contig1_id, reverse_direction(contig1_dir)),
            )
            directional.add(pair)
            directional.add(reverse_pair)
            agnostic.add(tuple(sorted((contig1_id, contig2_id))))

    return directional, agnostic


def evaluate_prediction_dir(prediction_dir: Path) -> dict:
    results = defaultdict(dict)
    bin_dirs = [path for path in prediction_dir.iterdir() if path.is_dir()]

    for bin_dir in sorted(bin_dirs):
        truth_file = bin_dir / "ground_truth.txt"
        prediction_file = bin_dir / "connect_dir.r"

        if not truth_file.exists() or not prediction_file.exists():
            continue

        truth_directional, truth_agnostic = parse_ground_truth(truth_file)
        if not truth_agnostic:
            continue

        total_predicted = 0
        correct_directional = 0
        correct_agnostic = 0

        with prediction_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                parts = line.strip().split()
                if len(parts) < 2:
                    continue

                for index in range(len(parts) - 1):
                    total_predicted += 1
                    contig_a = parts[index]
                    contig_b = parts[index + 1]
                    contig_a_id, contig_a_dir = contig_a[:-1], contig_a[-1]
                    contig_b_id, contig_b_dir = contig_b[:-1], contig_b[-1]

                    directional_pair = (
                        (contig_a_id, contig_a_dir),
                        (contig_b_id, contig_b_dir),
                    )
                    agnostic_pair = tuple(sorted((contig_a_id, contig_b_id)))

                    if directional_pair in truth_directional:
                        correct_directional += 1
                    if agnostic_pair in truth_agnostic:
                        correct_agnostic += 1

        truth_directional_total = len(truth_directional) / 2
        truth_agnostic_total = len(truth_agnostic)

        results[bin_dir.name]["directional"] = {
            "correct": correct_directional,
            "predicted_total": total_predicted,
            "ground_truth_total": truth_directional_total,
            "precision": correct_directional / total_predicted if total_predicted else 0.0,
            "recall": correct_directional / truth_directional_total if truth_directional_total else 0.0,
        }
        results[bin_dir.name]["agnostic"] = {
            "correct": correct_agnostic,
            "predicted_total": total_predicted,
            "ground_truth_total": truth_agnostic_total,
            "precision": correct_agnostic / total_predicted if total_predicted else 0.0,
            "recall": correct_agnostic / truth_agnostic_total if truth_agnostic_total else 0.0,
        }

    return results


def print_summary(results: dict) -> None:
    if not results:
        print("No evaluable benchmark bins were found.")
        return

    total_correct_directional = sum(item["directional"]["correct"] for item in results.values())
    total_correct_agnostic = sum(item["agnostic"]["correct"] for item in results.values())
    total_predicted = sum(item["directional"]["predicted_total"] for item in results.values())
    total_truth_directional = sum(item["directional"]["ground_truth_total"] for item in results.values())
    total_truth_agnostic = sum(item["agnostic"]["ground_truth_total"] for item in results.values())

    print("Per-bin summary:")
    for bin_name, item in sorted(results.items()):
        directional = item["directional"]
        agnostic = item["agnostic"]
        print(
            f"  {bin_name}: directional P={directional['precision']:.3f}, R={directional['recall']:.3f}; "
            f"agnostic P={agnostic['precision']:.3f}, R={agnostic['recall']:.3f}"
        )

    print("\nOverall summary:")
    print(
        f"  directional P={total_correct_directional / total_predicted if total_predicted else 0.0:.3f}, "
        f"R={total_correct_directional / total_truth_directional if total_truth_directional else 0.0:.3f}"
    )
    print(
        f"  agnostic P={total_correct_agnostic / total_predicted if total_predicted else 0.0:.3f}, "
        f"R={total_correct_agnostic / total_truth_agnostic if total_truth_agnostic else 0.0:.3f}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate ArchLink linking predictions against per-bin ground truth."
    )
    parser.add_argument("prediction_dir", type=Path, help="Directory containing benchmark bin subdirectories.")
    args = parser.parse_args()

    results = evaluate_prediction_dir(args.prediction_dir.resolve())
    print_summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
