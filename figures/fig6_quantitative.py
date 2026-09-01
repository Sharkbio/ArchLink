#!/usr/bin/env python3
"""Generate the quantitative/data-driven panels available for Figure 6."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from common import add_common_arguments, configure, load_sheet, panel_label, save_figure, style_axis


def panel_a(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig6a_scaffold_summary",
        ["Record_Type", "Feature_ID", "Source_Label", "Contig_Length_bp", "Strand", "Pre_Length_kb", "Post_Length_kb", "Junction_Identity_percent"],
    )
    contigs = frame[frame["Record_Type"] == "Contig"].copy()
    junctions = frame[frame["Record_Type"] == "Junction"].copy()
    bgcs = frame[frame["Record_Type"] == "BGC"].copy()
    figure, axis = plt.subplots(figsize=(9.2, 2.2))
    x = 0.0
    gap_width = 12000.0
    colors = ["#6E6E6E", "#8A8A8A"]
    for index, row in enumerate(contigs.itertuples()):
        length = float(row.Contig_Length_bp)
        axis.barh(0, length, left=x, height=0.32, color=colors[index % 2], edgecolor="#333333")
        axis.text(x + length / 2, 0, str(row.Source_Label), ha="center", va="center", color="white", fontsize=6)
        axis.text(x + length / 2, -0.28, f"{length/1000:.1f} kb ({row.Strand})", ha="center", fontsize=6)
        x += length
        if index < len(junctions):
            junction = junctions.iloc[index]
            axis.plot([x, x + gap_width], [0, 0], color="#D95F02", linewidth=2)
            axis.text(x + gap_width / 2, 0.2, f"{float(junction['Junction_Identity_percent']):.1f}%", ha="center", fontsize=6, color="#D95F02")
            x += gap_width
    for offset, row in enumerate(bgcs.itertuples()):
        axis.text(
            0.02,
            0.92 - offset * 0.13,
            f"{row.Source_Label}: {float(row.Pre_Length_kb):.1f} -> {float(row.Post_Length_kb):.1f} kb",
            transform=axis.transAxes,
            color="#1B9E77" if offset == 0 else "#D95F02",
            fontsize=7,
        )
    axis.set_xlim(0, x)
    axis.set_ylim(-0.5, 0.6)
    axis.set_yticks([])
    axis.set_xlabel("Scaffold coordinate (bp; gap widths schematic)")
    panel_label(axis, "a")
    figure.tight_layout()
    return figure


def panel_d(source_data: Path) -> plt.Figure:
    frame = load_sheet(source_data, "Fig6d_PAE_matrix", ["Residue_ID"])
    matrix = frame.drop(columns=["Residue_ID"]).to_numpy(dtype=float)
    figure, axis = plt.subplots(figsize=(5.0, 4.4))
    image = axis.imshow(matrix, cmap="viridis", vmin=0, vmax=30, interpolation="nearest", rasterized=True)
    axis.set_xlabel("Scored residue")
    axis.set_ylabel("Aligned residue")
    figure.colorbar(image, ax=axis, label="Predicted aligned error (A)", fraction=0.05)
    panel_label(axis, "d")
    figure.tight_layout()
    return figure


def panel_e(source_data: Path) -> plt.Figure:
    frame = load_sheet(
        source_data,
        "Fig6e_depth_profile",
        ["Position_bp", "Depth", "Clipped_Read_Count", "Target_Gene", "Gene_Start_bp", "Gene_End_bp"],
    )
    position_kb = frame["Position_bp"].astype(float) / 1000.0
    figure, axis = plt.subplots(figsize=(8.2, 3.0))
    axis.plot(position_kb, frame["Depth"], color="#1F78B4", linewidth=1, label="Coverage")
    axis.set_xlabel("Genomic position (kb)")
    axis.set_ylabel("Coverage (x)", color="#1F78B4")
    axis.tick_params(axis="y", colors="#1F78B4")
    second = axis.twinx()
    second.bar(position_kb, frame["Clipped_Read_Count"], width=0.002, color="#E67E22", alpha=0.7, label="Clipped reads")
    second.set_ylabel("Clipped reads", color="#E67E22")
    second.tick_params(axis="y", colors="#E67E22")
    start = float(frame["Gene_Start_bp"].iloc[0]) / 1000.0
    end = float(frame["Gene_End_bp"].iloc[0]) / 1000.0
    axis.axvspan(start, end, color="#F7E36D", alpha=0.35, label=str(frame["Target_Gene"].iloc[0]))
    axis.legend(frameon=False, loc="upper left")
    style_axis(axis)
    panel_label(axis, "e")
    figure.tight_layout()
    return figure


def generate(source_data: Path, output: Path, dpi: int, seed: int, formats: Sequence[str]) -> list[Path]:
    configure(seed)
    paths = []
    paths += save_figure(panel_a(source_data), output, "Fig6a", formats, dpi)
    paths += save_figure(panel_d(source_data), output, "Fig6d", formats, dpi)
    paths += save_figure(panel_e(source_data), output, "Fig6e", formats, dpi)
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
