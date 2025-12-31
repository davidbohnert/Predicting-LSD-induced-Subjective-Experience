#!/usr/bin/env python3
"""
Create lower-triangular network-by-network edge-count matrix plots from CPM edge masks.

Expected inputs
---------------
- Edge-mask CSVs (no header), e.g.:
    LSD_GDE_pos.csv
    LSD_GDE_difference_pos.csv
    LSD_AED_pos.csv
    LSD_AED_neg.csv
    LSD_MEQ30_difference_pos.csv

- Atlas label file with at least two columns, where column 2 contains the network label
  (e.g., Schaefer416_8networks.txt as whitespace-delimited text).

Outputs
-------
- PNG files in: <output_dir>/matrix_plots/
"""

from __future__ import annotations

from pathlib import Path
import logging

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------- Configuration -------------------------------

# Edit these paths to match your project layout (or pass via CLI if you add argparse).
DEFAULT_OUTPUT_DIR = Path("results/visualization/plots")
DEFAULT_EDGES_DIR = Path("results/edges")
DEFAULT_ATLAS_FILE = Path("data/Atlas/Schaefer416_8networks.txt")

SIZE_TITLE = 25
SIZE_SUBTITLE = SIZE_TITLE - 6
SIZE_AXIS = 18

# Specify exactly the plots you want to generate (run once).
PLOTS_TO_GENERATE = [
    ("GDE", "pos"),
    ("GDE_difference", "pos"),
    ("AED", "pos"),
    ("AED", "neg"),
    ("MEQ30_difference", "pos"),
]


# ------------------------------- Plot helpers -------------------------------

def set_title_and_subtitle(
    ax: plt.Axes,
    title: str,
    subtitle: str,
    title_size: int,
    subtitle_size: int,
    pad: int = 35,
    subtitle_offset: float = 1.02,
) -> None:
    ax.set_title(title, fontsize=title_size, fontweight="bold", pad=pad)
    ax.text(
        0.5,
        subtitle_offset,
        subtitle,
        transform=ax.transAxes,
        ha="center",
        va="bottom",
        fontsize=subtitle_size,
        fontweight="normal",
    )


def load_atlas_networks(atlas_file: Path) -> tuple[list[str], dict[str, np.ndarray]]:
    """
    Returns:
      network_names: list of unique network labels in order of appearance
      network_indices: mapping network label -> node indices (np.ndarray)
    """
    atlas_df = pd.read_csv(atlas_file, sep=r"\s+", header=None)

    # Keep order of appearance (drop duplicates without sorting)
    network_names = atlas_df[1].drop_duplicates().tolist()
    network_indices = {net: np.where(atlas_df[1].values == net)[0] for net in network_names}

    return network_names, network_indices


def compute_edge_counts(
    connectivity_matrix: np.ndarray,
    network_names: list[str],
    network_indices: dict[str, np.ndarray],
) -> np.ndarray:
    """
    Computes the network-by-network sum of edges, and sets the upper triangle to NaN.
    """
    n = len(network_names)
    edge_counts = np.zeros((n, n), dtype=float)

    for i, net_i in enumerate(network_names):
        idx_i = network_indices[net_i]
        for j, net_j in enumerate(network_names):
            idx_j = network_indices[net_j]
            sub = connectivity_matrix[np.ix_(idx_i, idx_j)]
            edge_counts[i, j] = np.sum(sub)

    # Exclude upper triangle
    iu = np.triu_indices_from(edge_counts, k=1)
    edge_counts[iu] = np.nan

    return edge_counts


def choose_cmap(outcome: str) -> str:
    """
    Match the original behavior:
      GDE -> Reds, MEQ30 -> Greens, AED -> Blues, else gray.
    """
    outcome_lower = outcome.lower()
    if "gde" in outcome_lower:
        return "Reds"
    if "meq" in outcome_lower:
        return "Greens"
    if "aed" in outcome_lower:
        return "Blues"
    return "gray"


def format_titles(outcome: str, model: str) -> tuple[str, str]:
    """
    Returns (title, subtitle) consistent with the original intent.
    """
    outcome_is_diff = "difference" in outcome.lower()
    base_outcome = outcome.replace("_difference", "")

    title = f"Edges predictive of {base_outcome}"
    model_label = "Positive Model" if model.lower() == "pos" else "Negative Model"
    connectome_label = "LSD–Placebo Connectome" if outcome_is_diff else "LSD Connectome"
    subtitle = f"{model_label}, {connectome_label}"
    return title, subtitle


def plot_edge_count_matrix(
    edge_counts: np.ndarray,
    network_names: list[str],
    cmap: str,
    title: str,
    subtitle: str,
    out_png: Path,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 8))
    cax = ax.matshow(edge_counts, cmap=cmap, interpolation="none")

    # Colorbar styling (as in original)
    cbar = plt.colorbar(cax, fraction=0.046, pad=0.04)
    cbar.set_label("Number of Edges", fontsize=SIZE_AXIS + 2, fontweight="bold", labelpad=15)
    cbar.ax.tick_params(labelsize=12)
    for t in cbar.ax.get_yticklabels():
        t.set_fontweight("bold")

    # Gridlines for lower triangle only (including diagonal)
    n = len(network_names)
    for i in range(n):
        for j in range(n):
            if i >= j:
                rect = plt.Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="black",
                    linewidth=1,
                )
                ax.add_patch(rect)

    # Axis labels / ticks
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(network_names, rotation=90, fontsize=SIZE_AXIS, fontweight="bold")
    ax.set_yticklabels(network_names, fontsize=SIZE_AXIS, fontweight="bold")
    ax.xaxis.tick_bottom()

    ax.set_xlabel("Networks", fontsize=SIZE_AXIS + 2, fontweight="bold")
    ax.set_ylabel("Networks", fontsize=SIZE_AXIS + 2, fontweight="bold")

    set_title_and_subtitle(ax, title, subtitle, SIZE_TITLE, SIZE_SUBTITLE)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------- Main ------------------------------------

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    output_dir = DEFAULT_OUTPUT_DIR
    edges_dir = DEFAULT_EDGES_DIR
    atlas_file = DEFAULT_ATLAS_FILE

    network_names, network_indices = load_atlas_networks(atlas_file)

    for outcome, model in PLOTS_TO_GENERATE:
        # Build CSV name exactly as requested (good->GDE, AIA->AED, no MEQ30 without difference)
        csv_file = edges_dir / f"LSD_{outcome}_{model}.csv"
        if not csv_file.exists():
            logging.warning(f"Missing input file, skipping: {csv_file}")
            continue

        connectivity_matrix = pd.read_csv(csv_file, header=None).values
        edge_counts = compute_edge_counts(connectivity_matrix, network_names, network_indices)

        cmap = choose_cmap(outcome)
        title, subtitle = format_titles(outcome, model)

        out_png = output_dir / "matrix_plots" / f"LSD_{outcome}_{model}.png"
        plot_edge_count_matrix(edge_counts, network_names, cmap, title, subtitle, out_png)

        logging.info(f"Saved: {out_png}")


if __name__ == "__main__":
    main()
