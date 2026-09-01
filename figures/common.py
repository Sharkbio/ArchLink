"""Shared deterministic plotting utilities for ArchLink manuscript figures."""

from __future__ import annotations

import argparse
import hashlib
import random
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE_DATA = ROOT / "source_data" / "Source_Data.xlsx"
DEFAULT_OUTPUT = Path(__file__).resolve().parent / "output"
DEFAULT_SEED = 20260831
DEFAULT_DPI = 300

DATASET_ORDER = ["Airways", "Skin", "Oral", "Waste Water", "Deep Marine", "VD Stool"]
MAIN_REAL_DATASET_ORDER = ["Waste Water", "Deep Marine", "VD Stool"]
TOOL_ORDER = ["MetaBAT2", "VAMB", "TaxVAMB", "COMEBin", "ArchLink"]
TOOL_COLORS = {
    "MetaBAT2": "#9E9E9E",
    "VAMB": "#F28E7F",
    "TaxVAMB": "#F6B2A6",
    "COMEBin": "#78C2B3",
    "ArchLink": "#0B559F",
    "ArchLink (P Only)": "#75AADB",
    "ArchLink (T Only)": "#B7D4DF",
    "ArchLink (P+T Full)": "#0B559F",
    "MetaCarvel (Default)": "#9E9E9E",
    "MetaCarvel (b=1)": "#555555",
}


def configure(seed: int = DEFAULT_SEED) -> None:
    random.seed(seed)
    np.random.seed(seed)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.bbox": "tight",
        }
    )


def load_sheet(source_data: Path, sheet: str, required: Iterable[str]) -> pd.DataFrame:
    if not source_data.is_file():
        raise FileNotFoundError(
            f"Source Data workbook not found: {source_data}. "
            "Pass --source-data or use source_data/Source_Data.xlsx."
        )
    frame = pd.read_excel(source_data, sheet_name=sheet, header=2, engine="openpyxl")
    frame = frame.dropna(how="all")
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Sheet {sheet!r} is missing columns: {missing}")
    return frame


def ordered(values: Iterable[str], preferred: Sequence[str]) -> list[str]:
    present = list(dict.fromkeys(str(value) for value in values))
    return [value for value in preferred if value in present] + sorted(
        value for value in present if value not in preferred
    )


def style_axis(axis: plt.Axes) -> None:
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.6, zorder=0)
    axis.tick_params(length=3, width=0.7)


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.16,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )


def save_figure(
    figure: plt.Figure,
    output_dir: Path,
    basename: str,
    formats: Sequence[str],
    dpi: int,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for extension in formats:
        path = output_dir / f"{basename}.{extension}"
        figure.savefig(path, dpi=dpi)
        paths.append(path)
    plt.close(figure)
    return paths


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--source-data", type=Path, default=DEFAULT_SOURCE_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--formats",
        nargs="+",
        choices=("pdf", "png", "svg"),
        default=("pdf",),
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
