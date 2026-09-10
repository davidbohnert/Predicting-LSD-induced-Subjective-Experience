#!/usr/bin/env python3
"""Render network concentration using a shared linear white-to-red scale.

Cells show observed/expected edge concentration and unique edge counts.
Inference is reported in captions and tables, without significance markers.
The default output is three paper figures. --diagnostics adds supporting
figures; --scale log2 --suffix _log2 requests an alternative diagnostic scale.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


OUTCOMES = ("GDE", "MEQ30", "VRS")
NETWORKS = ("VIS", "SMN", "DAN", "SAL", "LIM", "FPN", "DMN", "SBC")
PNG_DPI = 300

# Figures are laid out at their true reproduction width so nothing is downscaled at
# submission. 140 mm is a 1.5-column figure; 180 mm is full double-column width.
MM_PER_INCH = 25.4
SINGLE_WIDTH_IN = 140.0 / MM_PER_INCH
DOUBLE_WIDTH_IN = 180.0 / MM_PER_INCH

# Shared limits preserve comparability across outcomes.
LINEAR_VMIN, LINEAR_VMAX = 0.0, 5.5
LOG2_VMIN, LOG2_VMAX = -3.7, 2.5

# Above this fraction of the ramp the fill is dark enough to need light text.
LIGHT_TEXT_FRACTION = 0.60

ZERO_EDGE_COLOUR = "#BDBDBD"
INK = "#1A1A1A"
MUTED_INK = "#5F5F5F"
# Scale names shown as the centred plot title. Kept short: the connectome, the
# covariates and the definition of the ratio all belong in the caption.
BEHAVIOUR_TITLES = {
    "GDE": "Good Drug Effect (GDE)",
    "MEQ30": "Mystical Experience (MEQ30)",
    "VRS": "Visionary Restructuralization (VRS)",
}


def bh_fdr(values: np.ndarray) -> np.ndarray:
    """Return Benjamini-Hochberg q-values in the original row order."""
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def prepare_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Validate whole-mask inference and support provenance-only legacy redraws."""
    required = {"Outcome", "TotalMaskPermutationP"}
    missing = required.difference(summary.columns)
    if missing:
        raise ValueError(f"Missing outcome-summary columns: {sorted(missing)}")
    summary = summary.copy()
    expected_q = bh_fdr(summary["TotalMaskPermutationP"].to_numpy(dtype=float))
    if "TotalMaskFDRQ" in summary.columns:
        if not np.allclose(
            summary["TotalMaskFDRQ"].to_numpy(dtype=float),
            expected_q,
            rtol=0,
            atol=1e-12,
        ):
            raise ValueError("Whole-mask q-values are inconsistent with BH correction.")
    else:
        # Earlier aggregate tables omitted the derived whole-mask q-value.
        summary["TotalMaskFDRQ"] = expected_q
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--diagnostics", action="store_true")
    parser.add_argument(
        "--scale",
        choices=("linear", "log2"),
        default="linear",
        help="Mapping from the observed / expected ratio to colour.",
    )
    parser.add_argument(
        "--suffix",
        default="",
        help=(
            "Appended to enrichment figure stems. When set, only the enrichment "
            "figures are rendered, so a scale-comparison pass does not redraw the "
            "raw-count and null-count figures."
        ),
    )
    return parser.parse_args()


def apply_publication_style() -> None:
    """One rcParams pass so every figure in the set reads as a single system."""
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            # TrueType rather than Type 3, which many journals reject outright.
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            # Keep SVG text as text so production can restyle it.
            "svg.fonttype": "none",
            "axes.titlesize": 9,
            "axes.labelsize": 7,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "legend.fontsize": 7,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def enrichment_cmap() -> LinearSegmentedColormap:
    """Single-hue white-to-red ramp: legible in greyscale and CVD-safe."""
    cmap = LinearSegmentedColormap.from_list(
        "enrichment_sequential",
        ["#FFFFFF", "#FDDCCB", "#F9A582", "#E86A4A", "#C42D22", "#8B0F13"],
        N=256,
    )
    cmap.set_under(ZERO_EDGE_COLOUR)
    cmap.set_bad("#FFFFFF")
    return cmap


def enrichment_norm(scale: str) -> Normalize:
    if scale == "linear":
        return Normalize(vmin=LINEAR_VMIN, vmax=LINEAR_VMAX)
    return Normalize(vmin=LOG2_VMIN, vmax=LOG2_VMAX)


def enrichment_display(enrich: np.ndarray, counts: np.ndarray, scale: str) -> np.ndarray:
    """Map the observed / expected ratio onto the plotted scale, zeros under vmin."""
    display = np.full_like(enrich, np.nan)
    positive = np.isfinite(enrich) & (enrich > 0)
    if scale == "linear":
        display[positive] = np.clip(enrich[positive], LINEAR_VMIN, LINEAR_VMAX)
        sentinel = LINEAR_VMIN - 1.0
    else:
        display[positive] = np.clip(np.log2(enrich[positive]), LOG2_VMIN, LOG2_VMAX)
        sentinel = LOG2_VMIN - 1.0
    zero_cells = np.isfinite(enrich) & (counts == 0)
    display[zero_cells] = sentinel
    return display


def matrix_from_rows(frame: pd.DataFrame, column: str) -> np.ndarray:
    matrix = np.full((len(NETWORKS), len(NETWORKS)), np.nan, dtype=float)
    network_index = {name: index for index, name in enumerate(NETWORKS)}
    for row in frame.itertuples(index=False):
        i = network_index[str(row.Network2)]
        j = network_index[str(row.Network1)]
        if i < j:
            i, j = j, i
        matrix[i, j] = float(getattr(row, column))
    if np.count_nonzero(np.isfinite(matrix)) != 36:
        raise ValueError("Every outcome must map to exactly 36 lower-triangle cells.")
    return matrix


def setup_matrix_axis(
    ax: plt.Axes,
    title: str,
    tick_size: float,
    title_size: float = 11.0,
    title_weight: str = "bold",
    subtitle: str = "",
    subtitle_size: float = 7.0,
    title_pad: float = 18.0,
    panel_letter: str | None = None,
    letter_size: float = 9.0,
) -> None:
    n = len(NETWORKS)
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(NETWORKS, fontsize=tick_size)
    ax.set_yticklabels(NETWORKS, fontsize=tick_size)
    ax.xaxis.tick_bottom()
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(n - 0.5, -0.5)
    # Tick labels are already the network names; an axis label would only repeat them.
    if title:
        ax.set_title(
            title, loc="center", pad=title_pad, fontsize=title_size, fontweight=title_weight
        )
    if subtitle:
        ax.text(
            0.5,
            1.012,
            subtitle,
            transform=ax.transAxes,
            ha="center",
            va="bottom",
            fontsize=subtitle_size,
            color=MUTED_INK,
        )
    # Panel letters sit top-left, the descriptive title stays centred above them.
    if panel_letter:
        ax.set_title(panel_letter, loc="left", fontweight="bold", fontsize=letter_size)
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_cell_borders(ax: plt.Axes, linewidth: float = 0.5) -> None:
    for i in range(len(NETWORKS)):
        for j in range(i + 1):
            ax.add_patch(
                Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="#8C8C8C",
                    linewidth=linewidth,
                    zorder=3,
                )
            )


def annotate_cells(
    ax: plt.Axes,
    enrich: np.ndarray,
    counts: np.ndarray,
    display: np.ndarray,
    norm: Normalize,
    value_size: float,
    count_size: float | None,
) -> None:
    """The ratio leads; the edge count sits beneath it at a lighter weight."""
    for i in range(len(NETWORKS)):
        for j in range(i + 1):
            count = int(counts[i, j])
            value = enrich[i, j]
            if count == 0:
                value_colour, count_colour = INK, MUTED_INK
            else:
                fraction = float(norm(display[i, j]))
                dark_fill = fraction >= LIGHT_TEXT_FRACTION
                value_colour = "#FFFFFF" if dark_fill else INK
                count_colour = "#E8E8E8" if dark_fill else MUTED_INK
            if count_size is None:
                # Narrow cells: the colourbar already carries the unit, so drop the
                # "x" to keep the glyphs at a legible size inside a ~6 mm cell.
                ax.text(
                    j,
                    i,
                    f"{value:.2f}",
                    ha="center",
                    va="center",
                    fontsize=value_size,
                    color=value_colour,
                    zorder=6,
                )
                continue
            ax.text(
                j,
                i - 0.14,
                f"{value:.2f}×",
                ha="center",
                va="center",
                fontsize=value_size,
                color=value_colour,
                zorder=6,
            )
            ax.text(
                j,
                i + 0.21,
                f"({count})",
                ha="center",
                va="center",
                fontsize=count_size,
                color=count_colour,
                zorder=6,
            )


def add_enrichment_colorbar(
    fig: plt.Figure,
    image,
    axes,
    scale: str,
    label_size: float,
    tick_size: float,
    fraction: float = 0.046,
    # aspect is length/width, so this sets bar length at a fixed width. 18.5 puts the
    # bar at about two-thirds of the matrix height: the conventional 0.6-1.0 range,
    # consistent with the composite, and not so tall that it is balanced against the
    # empty upper-right of a lower-triangular matrix.
    aspect: float = 18.5,
):
    colorbar = fig.colorbar(image, ax=axes, fraction=fraction, aspect=aspect, pad=0.03)
    if scale == "linear":
        ticks = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        labels = [f"{tick:g}×" for tick in ticks]
    else:
        raw = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0])
        ticks = np.log2(raw)
        labels = [f"{value:g}×" for value in raw]
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels(labels)
    colorbar.ax.tick_params(labelsize=tick_size, length=2, width=0.5, pad=2)
    colorbar.set_label("Observed / expected edges", fontsize=label_size)
    colorbar.outline.set_linewidth(0.5)
    colorbar.outline.set_edgecolor("#8C8C8C")
    return colorbar


def legend_text(with_counts: bool) -> str:
    cell_line = (
        "Cell: observed / expected edges; consensus-edge count in brackets"
        if with_counts
        else "Cell: observed / expected edges"
    )
    meaning = (
        "1× = the edges this pair would hold if edges were spread evenly "
        "across all connections"
    )
    return f"{cell_line}\n{meaning}"


def save_figure(fig: plt.Figure, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        stem.with_suffix(".png"),
        dpi=PNG_DPI,
        metadata={"Software": "CPM network-interaction pipeline", "dpi": str(PNG_DPI)},
    )
    fig.savefig(stem.with_suffix(".svg"))
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def save_figure_tight(fig: plt.Figure, stem: Path) -> None:
    """Save a tightly cropped matrix for direct placement in the manuscript."""
    stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".svg", ".pdf"):
        fig.savefig(
            stem.with_suffix(suffix),
            dpi=PNG_DPI,
            bbox_inches="tight",
            pad_inches=0.02,
        )
    plt.close(fig)


def render_manuscript_enrichment(
    frame: pd.DataFrame, outcome: str, figure_dir: Path
) -> None:
    """Render the exact stripped layout embedded as Figures 1–3."""
    enrich = matrix_from_rows(frame, "ObservedExpectedRatio")
    counts = matrix_from_rows(frame, "ObservedConsensusEdges")
    display = enrichment_display(enrich, counts, "linear")
    norm = enrichment_norm("linear")

    fig, ax = plt.subplots(
        figsize=(SINGLE_WIDTH_IN, 5.5), constrained_layout=True
    )
    image = ax.imshow(
        display, cmap=enrichment_cmap(), norm=norm, interpolation="none"
    )
    setup_matrix_axis(ax, "", tick_size=7.5, subtitle="", title_pad=6.0)
    draw_cell_borders(ax)
    annotate_cells(
        ax, enrich, counts, display, norm, value_size=7.4, count_size=6.0
    )
    add_enrichment_colorbar(
        fig, image, ax, "linear", label_size=7.0, tick_size=6.5
    )
    ax.set_xlabel(
        f"{int(np.nansum(counts)):,} consensus edges",
        fontsize=7.0,
        color=MUTED_INK,
        labelpad=7,
    )
    save_figure_tight(fig, figure_dir / f"network_enrichment_{outcome}")


def render_enrichment(
    frame: pd.DataFrame, outcome: str, figure_dir: Path, scale: str, suffix: str
) -> None:
    enrich = matrix_from_rows(frame, "ObservedExpectedRatio")
    counts = matrix_from_rows(frame, "ObservedConsensusEdges")
    display = enrichment_display(enrich, counts, scale)
    norm = enrichment_norm(scale)

    fig, ax = plt.subplots(
        figsize=(SINGLE_WIDTH_IN, 5.5), constrained_layout=True
    )
    image = ax.imshow(
        display, cmap=enrichment_cmap(), norm=norm, interpolation="none"
    )
    total_edges = int(np.nansum(counts))
    setup_matrix_axis(
        ax,
        BEHAVIOUR_TITLES[outcome],
        tick_size=7.5,
        title_size=11.0,
        subtitle=f"{total_edges:,} consensus edges",
        subtitle_size=7.0,
        title_pad=18.0,
    )
    draw_cell_borders(ax)
    annotate_cells(ax, enrich, counts, display, norm, value_size=7.4, count_size=6.0)
    add_enrichment_colorbar(fig, image, ax, scale, label_size=7.0, tick_size=6.5)
    ax.set_xlabel(
        legend_text(with_counts=True),
        fontsize=6.3,
        color="#333333",
        labelpad=7,
    )
    save_figure(fig, figure_dir / f"network_enrichment_{outcome}{suffix}")


def render_enrichment_panel(
    pair_results: pd.DataFrame, figure_dir: Path, scale: str, suffix: str
) -> None:
    """Three outcomes at double-column width under one shared colourbar.

    At 180 mm across three panels a cell is about 6 mm, which will not hold two
    numbers at a legible size. Colour already carries the ratio, so the panel
    prints the enrichment value alone; consensus-edge counts stay on the standalone
    panels, in Table S, and in the caption.
    """
    norm = enrichment_norm(scale)
    cmap = enrichment_cmap()
    fig, axes = plt.subplots(
        1, 3, figsize=(DOUBLE_WIDTH_IN, 2.5), constrained_layout=True
    )
    image = None
    for ax, outcome, letter in zip(axes, OUTCOMES, "ABC", strict=True):
        frame = pair_results[pair_results["Outcome"] == outcome]
        enrich = matrix_from_rows(frame, "ObservedExpectedRatio")
        counts = matrix_from_rows(frame, "ObservedConsensusEdges")
        display = enrichment_display(enrich, counts, scale)
        image = ax.imshow(display, cmap=cmap, norm=norm, interpolation="none")
        setup_matrix_axis(
            ax,
            BEHAVIOUR_TITLES[outcome],
            tick_size=6.0,
            title_size=7.5,
            subtitle=f"{int(np.nansum(counts)):,} consensus edges",
            subtitle_size=5.8,
            title_pad=14.0,
            panel_letter=letter,
            letter_size=9.0,
        )
        draw_cell_borders(ax, linewidth=0.4)
        annotate_cells(
            ax, enrich, counts, display, norm, value_size=6.0, count_size=None
        )
    add_enrichment_colorbar(
        fig,
        image,
        list(axes),
        scale,
        label_size=6.5,
        tick_size=6.0,
        fraction=0.018,
        aspect=16.0,
    )
    fig.supxlabel(
        legend_text(with_counts=False), fontsize=6.0, color="#333333"
    )
    save_figure(fig, figure_dir / f"network_enrichment_panel{suffix}")


def render_raw_counts(frame: pd.DataFrame, outcome: str, figure_dir: Path) -> None:
    counts = matrix_from_rows(frame, "ObservedConsensusEdges")
    fig, ax = plt.subplots(
        figsize=(SINGLE_WIDTH_IN, 5.2), constrained_layout=True
    )
    image = ax.imshow(counts, cmap="Blues", interpolation="none", vmin=0)
    setup_matrix_axis(
        ax,
        BEHAVIOUR_TITLES[outcome],
        tick_size=7.5,
        title_size=11.0,
        subtitle=f"{int(np.nansum(counts)):,} consensus edges",
        subtitle_size=7.0,
        title_pad=18.0,
    )
    draw_cell_borders(ax)
    finite_counts = counts[np.isfinite(counts)]
    threshold = 0.58 * float(finite_counts.max()) if finite_counts.size else np.inf
    for i in range(len(NETWORKS)):
        for j in range(i + 1):
            count = int(counts[i, j])
            ax.text(
                j,
                i,
                str(count),
                ha="center",
                va="center",
                fontsize=7.0,
                color="#FFFFFF" if count >= threshold else INK,
                zorder=5,
            )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, aspect=18.5, pad=0.03)
    colorbar.set_label("Consensus edges", fontsize=7.0)
    colorbar.ax.tick_params(labelsize=6.5, length=2, width=0.5, pad=2)
    colorbar.outline.set_linewidth(0.5)
    colorbar.outline.set_edgecolor("#8C8C8C")
    ax.set_xlabel(
        "Each undirected edge counted once", fontsize=6.3, color="#333333", labelpad=7
    )
    save_figure(fig, figure_dir / f"raw_consensus_edge_counts_{outcome}")


def render_null_distributions(
    null_counts: pd.DataFrame, summary: pd.DataFrame, figure_dir: Path
) -> None:
    fig, axes = plt.subplots(
        1, 3, figsize=(DOUBLE_WIDTH_IN, 3.4), constrained_layout=True
    )
    summary = summary.set_index("Outcome")
    for ax, outcome, letter in zip(axes, OUTCOMES, "ABC", strict=True):
        values = null_counts.loc[
            null_counts["Outcome"] == outcome, "PositiveConsensusEdgeCount"
        ].to_numpy(dtype=float)
        if len(values) == 0:
            raise ValueError(f"No null edge counts found for {outcome}.")
        observed = float(summary.loc[outcome, "ObservedPositiveConsensusEdges"])
        transformed = np.log10(values + 1.0)
        observed_transformed = np.log10(observed + 1.0)
        bins = np.linspace(0.0, max(float(transformed.max()), observed_transformed), 51)
        ax.hist(transformed, bins=bins, color="#7A9EBA", edgecolor="white", linewidth=0.3)
        total_p = float(summary.loc[outcome, "TotalMaskPermutationP"])
        ax.axvline(
            observed_transformed,
            color="#B2182B",
            linewidth=1.4,
            label=f"Observed: {observed:g}\np = {total_p:.4f}",
        )
        ax.set_title(f"{letter}   {outcome}", loc="left", pad=5)
        ax.set_xlabel("Consensus edges (log scale)", fontsize=7.0)
        ax.set_ylabel("Permutations", fontsize=7.0)
        tick_values = np.array([0, 1, 10, 100, 1000, 4000], dtype=float)
        tick_values = tick_values[tick_values <= max(values.max(), observed)]
        ax.set_xticks(np.log10(tick_values + 1.0))
        ax.set_xticklabels([f"{value:g}" for value in tick_values], fontsize=6.5)
        ax.tick_params(axis="y", labelsize=6.5)
        ax.legend(frameon=False, fontsize=6.2, handlelength=1.2)
        ax.spines[["top", "right"]].set_visible(False)
    save_figure(fig, figure_dir / "null_positive_consensus_edge_distributions")


def write_captions(
    pair_results: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: Path,
    scale: str,
    null_counts: pd.DataFrame,
) -> None:
    """Draft captions. All inference lives here, not on the figures."""
    scale_word = "linear" if scale == "linear" else "log2"
    counts = null_counts.groupby("Outcome").size()
    if counts.nunique() != 1:
        raise ValueError("Outcomes must have equal permutation counts.")
    permutations = int(counts.iloc[0])
    significant_pairs = int((pair_results["CountFDRQ"] < 0.05).sum())
    lines = ["# Figure captions", ""]
    lines.append(
        "Draft captions for the network-interaction figures. Colour encodes "
        "concentration, never direction: every edge shown is a positively "
        "correlated consensus edge. The figures carry no significance markers; "
        "the permutation results are reported here and in the supplementary table."
    )
    lines.append("")

    def top_pairs(frame: pd.DataFrame, how_many: int = 3) -> str:
        ranked = frame.sort_values("ObservedExpectedRatio", ascending=False).head(how_many)
        parts = [
            f"{row.NetworkPair} ({row.ObservedExpectedRatio:.2f}x, "
            f"p = {row.CountPermutationP:.3f})"
            for row in ranked.itertuples(index=False)
        ]
        return ", ".join(parts)

    lines.append("## Network-pair concentration (composite, panels A-C)")
    lines.append("")
    lines.append(
        "**Concentration of positive CPM consensus edges across network pairs.** "
        "Consensus edges were selected in at least 8 of 10 folds at p < .01 from the "
        "covariate-adjusted LSD-placebo difference connectome (covariates: study, "
        "sex, age, mean framewise displacement) for (A) GDE, (B) MEQ30 and (C) VRS. "
        "Each undirected edge is counted once, and within-network cells use "
        "n(n-1)/2 possible edges. Each cell gives the observed / expected ratio: the "
        "consensus edges in a network pair as a share of the edges possible there, "
        "divided by the same share across all 86,320 possible edges. A value of 1 "
        "means the pair holds exactly as many edges as an even spread would give it. "
        f"Colour is mapped on a {scale_word} scale with limits shared across all "
        "three panels so the outcomes remain comparable; grey cells contain no "
        "consensus edges. Network-pair concentration was tested by comparing the raw "
        f"consensus-edge count in each pair against {permutations:,} Freedman-Lane behavioural "
        "permutations, with Benjamini-Hochberg correction across the 36 pairs within "
        f"each outcome. {significant_pairs} outcome/network-pair tests survived correction. "
        "Concentration was greatest in "
    )
    for outcome in OUTCOMES:
        frame = pair_results[pair_results["Outcome"] == outcome]
        min_q = float(frame["CountFDRQ"].min())
        lines[-1] += (
            f"{outcome}: {top_pairs(frame)}, minimum q = {min_q:.3f}"
            + ("; " if outcome != OUTCOMES[-1] else ". ")
        )
    lines[-1] += (
        "These concentrations are descriptive effect sizes. Full per-pair statistics "
        "are given in the supplementary table."
    )
    lines.append("")

    for outcome in OUTCOMES:
        frame = pair_results[pair_results["Outcome"] == outcome]
        total = int(frame["ObservedConsensusEdges"].sum())
        empty = int((frame["ObservedConsensusEdges"] == 0).sum())
        min_q = float(frame["CountFDRQ"].min())
        lines.append(f"## Network-pair concentration, {outcome} (standalone)")
        lines.append("")
        lines.append(
            f"**{outcome}: concentration of positive consensus edges across network "
            f"pairs.** {total:,} consensus edges distributed over 36 network pairs. "
            "Each cell gives the observed / expected ratio with the consensus-edge "
            f"count beneath it in brackets; {empty} pair(s) contain no edges and are "
            f"shown in grey. Concentration was greatest in {top_pairs(frame)}. "
            f"{int((frame['CountFDRQ'] < 0.05).sum())} pairs survived FDR correction "
            f"across the 36 tests (minimum q = "
            f"{min_q:.3f}). Colour limits are shared with the other outcome panels. "
            "Details as in the composite figure caption."
        )
        lines.append("")

    lines.append("## Consensus-edge counts under permutation")
    lines.append("")
    lines.append(
        "**Total positive consensus-edge counts under the Freedman-Lane null.** "
        "Distribution of the number of edges retained in at least 8 of 10 folds "
        f"across {permutations:,} behavioural permutations, for (A) GDE, (B) MEQ30 and (C) VRS. "
        "The red line marks the observed count, with its one-sided permutation "
        "p-value. The horizontal axis is log-scaled with zero retained. "
        "Benjamini-Hochberg correction was applied across the three whole-mask "
        "tests: "
    )
    summary_by_outcome = summary.set_index("Outcome")
    mask_results = []
    for outcome in OUTCOMES:
        row = summary_by_outcome.loc[outcome]
        mask_results.append(
            f"{outcome}, p = {row.TotalMaskPermutationP:.4f}, "
            f"q = {row.TotalMaskFDRQ:.4f}"
        )
    lines[-1] += "; ".join(mask_results) + "."
    lines.append("")
    lines.append("## Raw consensus-edge counts (supplement)")
    lines.append("")
    lines.append(
        "**Protocol-style consensus-edge counts per network pair.** Unnormalised "
        "counts, with each undirected edge counted once so that within-network cells "
        "are not double-counted. These counts are not comparable across network "
        "pairs because pairs differ in the number of possible edges; the "
        "observed / expected figures provide the size-normalised view."
    )
    lines.append("")
    (output_dir / "figure_captions.md").write_text("\n".join(lines), encoding="utf-8")


def validate_inputs(pair_results: pd.DataFrame, null_counts: pd.DataFrame) -> None:
    required_pair = {
        "Outcome",
        "Network1",
        "Network2",
        "ObservedConsensusEdges",
        "ObservedExpectedRatio",
        "CountPermutationP",
        "CountFDRQ",
    }
    missing = required_pair.difference(pair_results.columns)
    if missing:
        raise ValueError(f"Missing network-pair columns: {sorted(missing)}")
    if set(pair_results["Outcome"].unique()) != set(OUTCOMES):
        raise ValueError("Network-pair outcomes do not match GDE, MEQ30, and VRS.")
    per_outcome = pair_results.groupby("Outcome").size()
    if not (per_outcome == 36).all():
        raise ValueError("Every outcome must contain exactly 36 pair rows.")
    if set(null_counts["Outcome"].unique()) != set(OUTCOMES):
        raise ValueError("Null-count outcomes do not match GDE, MEQ30, and VRS.")


def main() -> None:
    args = parse_args()
    apply_publication_style()
    output_dir = args.output_dir.resolve()
    pair_results = pd.read_csv(output_dir / "network_pair_results.csv")
    null_counts = pd.read_csv(output_dir / "null_consensus_edge_counts.csv")
    summary = pd.read_csv(output_dir / "outcome_summary.csv")
    summary = prepare_summary(summary)
    validate_inputs(pair_results, null_counts)
    # figures/ holds only the three paper figures; everything else sits in other/.
    figure_dir = output_dir / "figures"
    other_dir = figure_dir / "other"
    is_paper_scale = args.scale == "linear" and not args.suffix
    enrichment_dir = figure_dir if is_paper_scale else other_dir

    for outcome in OUTCOMES:
        frame = pair_results[pair_results["Outcome"] == outcome].copy()
        if is_paper_scale:
            render_manuscript_enrichment(frame, outcome, figure_dir)
            if args.diagnostics:
                render_enrichment(frame, outcome, other_dir, "linear", "_standalone")
        else:
            render_enrichment(frame, outcome, enrichment_dir, args.scale, args.suffix)
    if args.diagnostics or not is_paper_scale:
        render_enrichment_panel(pair_results, other_dir, args.scale, args.suffix)

    if is_paper_scale and args.diagnostics:
        for outcome in OUTCOMES:
            frame = pair_results[pair_results["Outcome"] == outcome].copy()
            render_raw_counts(frame, outcome, other_dir)
        render_null_distributions(null_counts, summary, other_dir)
    if is_paper_scale:
        write_captions(pair_results, summary, output_dir, args.scale, null_counts)

    print(
        f"Enrichment figures: {enrichment_dir}\n"
        f"Supporting figures: {other_dir}\n"
        f"(scale={args.scale}, suffix='{args.suffix}')"
    )


if __name__ == "__main__":
    main()
