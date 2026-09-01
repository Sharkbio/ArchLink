#!/usr/bin/env python3
"""Generate quantitative summaries for ArchLink Figure 4."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    MAIN_REAL_DATASET_ORDER,
    add_common_arguments,
    configure,
    load_sheet,
    panel_label,
    save_figure,
    style_axis,
)


RANKS = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def panel_a(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig4a_tax_rank_counts",
        [
            "Dataset",
            "Rank",
            "Shared_taxa_count",
            "ArchLink_bin_only_count",
            "ArchLink_scaffolded_only_count",
        ],
    )
    figure, axes = plt.subplots(1, 3, figsize=(9.2, 3.1), sharey=True)
    for axis, dataset in zip(axes, ["Deep Marine", "Waste Water", "VD Stool"]):
        subset = frame[frame["Dataset"] == dataset].set_index("Rank").reindex(RANKS)
        pre = subset["Shared_taxa_count"] + subset["ArchLink_bin_only_count"]
        post = subset["Shared_taxa_count"] + subset["ArchLink_scaffolded_only_count"]
        y = np.arange(len(RANKS))
        axis.plot(pre, y, "o-", color="#A8B6BE", label="Unscaffolded")
        axis.plot(post, y, "o-", color="#1F78B4", label="Scaffolded")
        for yi, left, right in zip(y, pre, post):
            axis.plot([left, right], [yi, yi], color="#D7D7D7", linewidth=0.7, zorder=0)
        axis.set_title(dataset)
        axis.set_xlabel("Taxa count")
        axis.set_yticks(y, RANKS)
        axis.invert_yaxis()
        axis.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    axes[-1].legend(frameon=False, loc="lower right")
    panel_label(axes[0], "a")
    figure.suptitle("Rank-count summary underlying the taxonomic transition panel", y=1.03, fontsize=10)
    figure.tight_layout()
    return figure


def panel_b(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig4b_align_frac_raw",
        ["Dataset", "Group_ID", "State", "ANI", "Aligned_Fraction"],
    )
    figure, axes = plt.subplots(1, 3, figsize=(9.2, 3.6), sharex=True)
    for axis, dataset in zip(axes, ["Waste Water", "Deep Marine", "VD Stool"]):
        subset = frame[frame["Dataset"] == dataset]
        groups = list(dict.fromkeys(subset["Group_ID"].tolist()))
        for y, group in enumerate(groups):
            rows = subset[subset["Group_ID"] == group].set_index("State")
            if not {"ArchLink_bin_only", "ArchLink_scaffolded"}.issubset(rows.index):
                continue
            pre = float(rows.loc["ArchLink_bin_only", "Aligned_Fraction"])
            post = float(rows.loc["ArchLink_scaffolded", "Aligned_Fraction"])
            axis.plot([pre, post], [y, y], color="#B0B0B0", linewidth=1)
            axis.plot(pre, y, "o", color="#A8B6BE", markersize=4)
            axis.plot(post, y, "o", color="#1F78B4", markersize=4)
        axis.set_title(dataset)
        axis.set_yticks(range(len(groups)), [f"group{int(group)}" for group in groups], fontsize=6)
        axis.set_xlabel("Aligned reference fraction (%)")
        axis.grid(axis="x", color="#E6E6E6", linewidth=0.6)
    axes[0].plot([], [], "o", color="#A8B6BE", label="Unscaffolded")
    axes[0].plot([], [], "o", color="#1F78B4", label="Scaffolded")
    axes[0].legend(frameon=False, loc="lower right")
    panel_label(axes[0], "b")
    figure.tight_layout()
    return figure


def panel_c(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig4c_marker_summary",
        ["Dataset_Label", "Delta_SCG", "Delta_Duplicated_Markers"],
    )
    frame = frame[frame["Dataset_Label"].isin(MAIN_REAL_DATASET_ORDER)].set_index("Dataset_Label").reindex(MAIN_REAL_DATASET_ORDER)
    x = np.arange(len(frame))
    width = 0.34
    figure, axis = plt.subplots(figsize=(5.2, 3.2))
    axis.bar(x - width / 2, frame["Delta_SCG"], width, color="#A8C6D6", label="SCG net change")
    axis.bar(x + width / 2, frame["Delta_Duplicated_Markers"], width, color="#E7B7A8", label="Duplicated-marker net change")
    axis.axhline(0, color="#444444", linewidth=0.8)
    axis.set_xticks(x, frame.index)
    axis.set_ylabel("Net marker count change")
    axis.legend(frameon=False)
    style_axis(axis)
    panel_label(axis, "c")
    figure.tight_layout()
    return figure


def metabolic_categories(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = frame[frame["Value_Type"] == "numeric_fraction"].copy()
    numeric["Value_Post"] = pd.to_numeric(numeric["Value_Post"], errors="coerce")
    categories = {
        "ETC complexes": ("Complex I:", "Complex II:", "Complex III:", "Complex IV", "Complex V:"),
        "Carbon fixation": (
            "Reductive pentose phosphate cycle",
            "Reductive citrate cycle",
            "3-Hydroxypropionate",
            "Hydroxypropionate-hydroxybutylate",
            "Dicarboxylate-hydroxybutyrate",
            "Acetyl-CoA pathway",
        ),
        "Central carbon": (
            "Glycolysis",
            "Pentose phosphate pathway",
            "Citrate cycle",
            "SCFA and alcohol conversions",
        ),
    }
    rows = []
    for (dataset, group), subset in numeric.groupby(["Dataset_Label", "Group_ID"], sort=False):
        row = {"Dataset_Label": dataset, "Group_ID": group}
        for label, prefixes in categories.items():
            selected = subset[subset["Module_Name"].astype(str).str.startswith(prefixes)]
            row[label] = selected["Value_Post"].mean() if not selected.empty else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def panel_d(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig4d_DRAM_raw",
        ["Dataset_Label", "Group_ID", "Module_Name", "Value_Type", "Value_Post"],
    )
    summary = metabolic_categories(frame)
    order = {name: index for index, name in enumerate(MAIN_REAL_DATASET_ORDER)}
    summary["dataset_order"] = summary["Dataset_Label"].map(order)
    summary = summary.sort_values(["dataset_order", "Group_ID"])
    columns = ["ETC complexes", "Carbon fixation", "Central carbon"]
    matrix = summary[columns].to_numpy(dtype=float)
    labels = [f"{row.Dataset_Label[:2]} group{int(row.Group_ID)}" for row in summary.itertuples()]
    figure, axis = plt.subplots(figsize=(5.4, 5.0))
    image = axis.imshow(matrix, aspect="auto", cmap="YlGn", vmin=0, vmax=1)
    axis.set_xticks(range(len(columns)), columns, rotation=35, ha="right")
    axis.set_yticks(range(len(labels)), labels, fontsize=6)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            if np.isfinite(matrix[row, column]):
                axis.text(column, row, f"{matrix[row, column]:.2f}", ha="center", va="center", fontsize=5.5)
    figure.colorbar(image, ax=axis, label="Mean post-scaffolding module completeness", fraction=0.04)
    panel_label(axis, "d")
    figure.tight_layout()
    return figure


def generate(source_data: Path, output: Path, dpi: int, seed: int, formats: Sequence[str]) -> list[Path]:
    configure(seed)
    paths = []
    paths += save_figure(panel_a(source_data), output, "Fig4a", formats, dpi)
    paths += save_figure(panel_b(source_data), output, "Fig4b", formats, dpi)
    paths += save_figure(panel_c(source_data), output, "Fig4c", formats, dpi)
    paths += save_figure(panel_d(source_data), output, "Fig4d", formats, dpi)
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
