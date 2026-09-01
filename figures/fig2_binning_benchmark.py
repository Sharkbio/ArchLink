#!/usr/bin/env python3
"""Generate the quantitative panels for ArchLink Figure 2."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt

from common import (
    TOOL_COLORS,
    TOOL_ORDER,
    add_common_arguments,
    configure,
    load_sheet,
    ordered,
    panel_label,
    save_figure,
    style_axis,
)


SIMULATED = ["Airways", "Skin", "Oral"]
REAL = ["Waste Water", "Deep Marine", "VD Stool"]


def panel_a(source_data: Path) -> plt.Figure:
    frame = load_sheet(source_data, "Fig2a_ARI_ACC", ["Dataset", "Metric", "Tool", "Score"])
    metrics = ["ARI (seq)", "Accuracy"]
    figure, axes = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=False)
    for axis, metric in zip(axes, metrics):
        subset = frame[frame["Metric"] == metric]
        for tool in ordered(subset["Tool"], TOOL_ORDER):
            values = (
                subset[subset["Tool"] == tool]
                .set_index("Dataset")
                .reindex(SIMULATED)["Score"]
            )
            axis.plot(
                SIMULATED,
                values,
                marker="o",
                linewidth=1.5,
                markersize=4,
                color=TOOL_COLORS.get(tool, "#333333"),
                label=tool,
            )
        axis.set_title(metric)
        axis.set_ylabel("Score")
        axis.set_ylim(bottom=0)
        style_axis(axis)
    axes[-1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    panel_label(axes[0], "a")
    figure.suptitle("Assembly-level binning quality", y=1.03, fontsize=10)
    figure.tight_layout()
    return figure


def recovery_panel(
    source_data: Path,
    sheet: str,
    datasets: list[str],
    count_column: str,
    precision: float,
    minimum_completeness: float,
    label: str,
    title: str,
) -> plt.Figure:
    frame = load_sheet(
        source_data,
        sheet,
        ["Dataset", "Tool", "Precision", "Completeness_Threshold", count_column],
    )
    frame = frame[
        (frame["Precision"].round(4) == precision)
        & (frame["Completeness_Threshold"] >= minimum_completeness)
    ]
    figure, axes = plt.subplots(1, 3, figsize=(9.2, 2.85), sharey=False)
    for axis, dataset in zip(axes, datasets):
        subset = frame[frame["Dataset"] == dataset]
        for tool in ordered(subset["Tool"], TOOL_ORDER):
            values = subset[subset["Tool"] == tool].sort_values("Completeness_Threshold")
            axis.plot(
                values["Completeness_Threshold"],
                values[count_column],
                marker="o",
                linewidth=1.4,
                markersize=3.5,
                color=TOOL_COLORS.get(tool, "#333333"),
                label=tool,
            )
        axis.set_title(dataset)
        axis.set_xlabel("Completeness threshold")
        axis.set_ylabel("Recovered bins" if count_column == "Recovered_Bin_Count" else "Recovered MAGs")
        style_axis(axis)
    axes[-1].legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    panel_label(axes[0], label)
    figure.suptitle(title, y=1.03, fontsize=10)
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
    paths += save_figure(panel_a(source_data), output, "Fig2a", formats, dpi)
    paths += save_figure(
        recovery_panel(
            source_data,
            "Fig2b_CAMI_bins",
            SIMULATED,
            "Recovered_Bin_Count",
            precision=0.90,
            minimum_completeness=0.10,
            label="b",
            title="CAMI II recovered bins (precision >= 90%)",
        ),
        output,
        "Fig2b",
        formats,
        dpi,
    )
    paths += save_figure(
        recovery_panel(
            source_data,
            "Fig2c_real_MAGs",
            REAL,
            "Recovered_MAG_Count",
            precision=0.95,
            minimum_completeness=0.50,
            label="c",
            title="Real-metagenome recovered MAGs (precision >= 95%)",
        ),
        output,
        "Fig2c",
        formats,
        dpi,
    )
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
