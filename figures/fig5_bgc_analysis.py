#!/usr/bin/env python3
"""Generate the quantitative panels for ArchLink Figure 5."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from common import (
    MAIN_REAL_DATASET_ORDER,
    TOOL_COLORS,
    add_common_arguments,
    configure,
    load_sheet,
    panel_label,
    save_figure,
    style_axis,
)


def panel_a(source_data: Path) -> plt.Figure:
    frame = load_sheet(source_data, "Fig5a_BGC_counts", ["Dataset", "Tool", "Category", "Count"])
    tools = ["MetaBAT2", "VAMB", "COMEBin", "ArchLink"]
    category_colors = {
        "Total": "#888888",
        "Bin only": "#A8C8D8",
        "Core extension": "#5B9CC5",
        "Novel BGC recovery": "#0B559F",
    }
    figure, axes = plt.subplots(1, 3, figsize=(9.2, 3.1), sharey=False)
    for axis, dataset in zip(axes, MAIN_REAL_DATASET_ORDER):
        subset = frame[frame["Dataset"] == dataset]
        x = np.arange(len(tools))
        bottom = np.zeros(len(tools))
        for category in ("Total", "Bin only", "Core extension", "Novel BGC recovery"):
            values = []
            for tool in tools:
                row = subset[(subset["Tool"] == tool) & (subset["Category"] == category)]
                values.append(float(row["Count"].iloc[0]) if not row.empty else 0.0)
            axis.bar(x, values, bottom=bottom, color=category_colors[category], edgecolor="#333333", linewidth=0.5, label=category)
            bottom += np.asarray(values)
        axis.set_title(dataset)
        axis.set_xticks(x, tools, rotation=25, ha="right")
        axis.set_ylabel("Number of BGCs")
        style_axis(axis)
    axes[-1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    panel_label(axes[0], "a")
    figure.tight_layout()
    return figure


def panel_b(source_data: Path) -> plt.Figure:
    frame = load_sheet(source_data, "Fig5b_BGC_lengths", ["Dataset_Label", "State", "Total_Length_bp"])
    states = ["ArchLink_bin_only", "ArchLink_scaffolded"]
    colors = ["#E4A07A", "#6E9CBF"]
    figure, axis = plt.subplots(figsize=(7.4, 3.8))
    positions = []
    labels = []
    datasets = ["Waste Water", "Deep Marine", "VD Stool"]
    for dataset_index, dataset in enumerate(datasets):
        for state_index, state in enumerate(states):
            values = frame[(frame["Dataset_Label"] == dataset) & (frame["State"] == state)]["Total_Length_bp"].astype(float).to_numpy()
            position = dataset_index * 3 + state_index
            parts = axis.violinplot(values, positions=[position], widths=0.8, showmeans=False, showmedians=True, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor(colors[state_index])
                body.set_edgecolor("#333333")
                body.set_alpha(0.75)
            parts["cmedians"].set_color("white")
            parts["cmedians"].set_linewidth(1.4)
            positions.append(position)
            labels.append("Bin only" if state_index == 0 else "Scaffolded")
    axis.set_yscale("log")
    axis.set_xticks(positions, labels, rotation=25, ha="right")
    for dataset_index, dataset in enumerate(datasets):
        axis.text(dataset_index * 3 + 0.5, axis.get_ylim()[1] * 0.90, dataset, ha="center", va="top", fontsize=8)
    axis.set_ylabel("Total BGC length (bp; log scale)")
    style_axis(axis)
    panel_label(axis, "b")
    figure.tight_layout()
    return figure


def panel_c(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig5c_BGC_class_resolved",
        ["Dataset_Label", "BGC_Class", "Median_Length_Gain_kb", "Complete_BGCs_Increased", "Total_Class_Count"],
    )
    class_order = (
        frame.groupby("BGC_Class")["Total_Class_Count"].sum().sort_values().index.tolist()
    )
    datasets = ["Deep Marine", "VD Stool", "Waste Water"]
    x_map = {dataset: index for index, dataset in enumerate(datasets)}
    y_map = {bgc_class: index for index, bgc_class in enumerate(class_order)}
    x = frame["Dataset_Label"].map(x_map).to_numpy()
    y = frame["BGC_Class"].map(y_map).to_numpy()
    size = np.sqrt(frame["Total_Class_Count"].astype(float).to_numpy()) * 14
    color = frame["Median_Length_Gain_kb"].astype(float).to_numpy()
    figure, axis = plt.subplots(figsize=(6.2, 4.6))
    scatter = axis.scatter(x, y, s=size, c=color, cmap="viridis", edgecolor="#333333", linewidth=0.4)
    for row, xi, yi in zip(frame.itertuples(), x, y):
        increased = int(row.Complete_BGCs_Increased)
        if increased:
            axis.text(xi, yi, str(increased), ha="center", va="center", fontsize=5.5, color="white" if row.Median_Length_Gain_kb > 8 else "black")
    axis.set_xticks(range(len(datasets)), datasets, rotation=20, ha="right")
    axis.set_yticks(range(len(class_order)), class_order, fontsize=6.5)
    axis.set_xlim(-0.6, len(datasets) - 0.4)
    axis.set_xlabel("Number inside bubble: boundary-complete BGCs increased")
    figure.colorbar(scatter, ax=axis, label="Median positive length gain (kb)", fraction=0.05)
    panel_label(axis, "c")
    figure.tight_layout()
    return figure


def generate(source_data: Path, output: Path, dpi: int, seed: int, formats: Sequence[str]) -> list[Path]:
    configure(seed)
    paths = []
    paths += save_figure(panel_a(source_data), output, "Fig5a", formats, dpi)
    paths += save_figure(panel_b(source_data), output, "Fig5b", formats, dpi)
    paths += save_figure(panel_c(source_data), output, "Fig5c", formats, dpi)
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
