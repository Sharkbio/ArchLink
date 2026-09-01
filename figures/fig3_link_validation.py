#!/usr/bin/env python3
"""Generate the quantitative panels for ArchLink Figure 3."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from common import (
    TOOL_COLORS,
    add_common_arguments,
    configure,
    load_sheet,
    ordered,
    panel_label,
    save_figure,
    style_axis,
)


TOOLS = [
    "MetaCarvel (Default)",
    "MetaCarvel (b=1)",
    "ArchLink (T Only)",
    "ArchLink (P Only)",
    "ArchLink (P+T Full)",
]
DATASETS = ["Airways", "Oral", "Skin"]


def panel_a(source_data: Path) -> plt.Figure:
    frame = load_sheet(source_data, "Fig3a_precision_TP", ["Dataset", "Tool", "TP", "FP", "Precision"])
    grouped = frame.groupby("Tool", sort=False).agg(
        precision=("Precision", "mean"),
        precision_sd=("Precision", "std"),
        tp=("TP", "mean"),
        tp_sd=("TP", "std"),
    )
    tools = [tool for tool in TOOLS if tool in grouped.index]
    x = np.arange(len(tools))
    colors = [TOOL_COLORS[tool] for tool in tools]
    figure, axes = plt.subplots(1, 2, figsize=(8.4, 3.2))
    axes[0].bar(x, grouped.loc[tools, "precision"], color=colors, edgecolor="#333333", linewidth=0.6)
    axes[0].errorbar(x, grouped.loc[tools, "precision"], yerr=grouped.loc[tools, "precision_sd"], fmt="none", color="#333333", capsize=2)
    axes[0].set_ylabel("Mean precision across datasets")
    axes[0].set_ylim(0, 1.08)
    axes[1].bar(x, grouped.loc[tools, "tp"], color=colors, edgecolor="#333333", linewidth=0.6)
    axes[1].errorbar(x, grouped.loc[tools, "tp"], yerr=grouped.loc[tools, "tp_sd"], fmt="none", color="#333333", capsize=2)
    axes[1].set_ylabel("Mean true-positive joins")
    for axis in axes:
        axis.set_xticks(x, [tool.replace("ArchLink ", "") for tool in tools], rotation=25, ha="right")
        style_axis(axis)
    panel_label(axes[0], "a")
    figure.suptitle("Scaffolding precision and recovered true-positive joins", y=1.02, fontsize=10)
    figure.tight_layout()
    return figure


def panel_b(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig3b_longread_main",
        ["Tool", "Validation_Rate", "Depth_Consistency", "Structural_Purity", "Anti_Chimeric_Safety"],
    )
    metrics = [
        ("Validation_Rate", "Validation rate"),
        ("Depth_Consistency", "Depth consistency"),
        ("Structural_Purity", "Structural purity"),
        ("Anti_Chimeric_Safety", "Anti-chimeric safety"),
    ]
    tools = ["MetaCarvel (Default)", "MetaCarvel (b=1)", "ArchLink"]
    figure, axes = plt.subplots(1, 4, figsize=(9.2, 2.7), sharey=True)
    for axis, (column, title) in zip(axes, metrics):
        for y, tool in enumerate(tools):
            value = float(frame.loc[frame["Tool"] == tool, column].iloc[0])
            axis.hlines(y, 0, value, color=TOOL_COLORS[tool], linewidth=1.2)
            axis.plot(value, y, "o", color=TOOL_COLORS[tool], markersize=5)
            axis.text(value + 0.7, y, f"{value:.1f}", va="center", fontsize=6.5)
        axis.set_title(title)
        axis.set_xlabel("Score (%)")
        axis.set_xlim(45 if column != "Structural_Purity" else 88, 102)
        axis.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    axes[0].set_yticks(range(len(tools)), tools)
    panel_label(axes[0], "b")
    figure.tight_layout()
    return figure


def panel_c(source_data: Path, seed: int) -> plt.Figure:
    raw = load_sheet(source_data, "Fig3c_N50_raw", ["Dataset", "Tool", "Gain_bp"])
    summary = load_sheet(
        source_data,
        "Fig3c_N50_summary",
        ["Dataset", "Tool", "Mean_Gain_bp", "CI_Lower_bp", "CI_Upper_bp", "Effect_Size_r_approx"],
    )
    tools = ["ArchLink", "MetaCarvel (b=1)"]
    offsets = {"ArchLink": 0.14, "MetaCarvel (b=1)": -0.14}
    rng = np.random.default_rng(seed)
    figure, axis = plt.subplots(figsize=(9.2, 3.5))
    for dataset_index, dataset in enumerate(DATASETS):
        for tool in tools:
            values = raw[(raw["Dataset"] == dataset) & (raw["Tool"] == tool)]["Gain_bp"].astype(float).to_numpy()
            y = dataset_index + offsets[tool] + rng.normal(0, 0.035, len(values))
            axis.scatter(values, y, s=7, alpha=0.20, color=TOOL_COLORS[tool], edgecolors="none")
            row = summary[(summary["Dataset"] == dataset) & (summary["Tool"] == tool)].iloc[0]
            mean = float(row["Mean_Gain_bp"])
            low = float(row["CI_Lower_bp"])
            high = float(row["CI_Upper_bp"])
            axis.errorbar(
                mean,
                dataset_index + offsets[tool],
                xerr=[[mean - low], [high - mean]],
                fmt="o",
                color=TOOL_COLORS[tool],
                capsize=2,
                markersize=4.5,
                linewidth=1.2,
                label=tool if dataset_index == 0 else None,
            )
            axis.text(
                high,
                dataset_index + offsets[tool] + (0.08 if tool == "ArchLink" else -0.08),
                f"r={float(row['Effect_Size_r_approx']):.2f}",
                color=TOOL_COLORS[tool],
                fontsize=6,
                ha="right",
            )
    axis.axvline(0, color="#777777", linestyle="--", linewidth=0.8)
    axis.set_xscale("symlog", linthresh=50)
    axis.set_yticks(range(len(DATASETS)), DATASETS)
    axis.set_xlabel("N50 gain (bp; symmetric log scale)")
    axis.set_ylabel("Dataset")
    axis.legend(frameon=False, loc="lower right")
    style_axis(axis)
    panel_label(axis, "c")
    figure.tight_layout()
    return figure


def generate(
    source_data: Path,
    output: Path,
    dpi: int,
    seed: int,
    formats: Sequence[str],
) -> list[Path]:
    configure(seed)
    paths = []
    paths += save_figure(panel_a(source_data), output, "Fig3a", formats, dpi)
    paths += save_figure(panel_b(source_data), output, "Fig3b", formats, dpi)
    paths += save_figure(panel_c(source_data, seed), output, "Fig3c", formats, dpi)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(parser)
    args = parser.parse_args()
    for path in generate(args.source_data, args.output, args.dpi, args.seed, args.formats):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
