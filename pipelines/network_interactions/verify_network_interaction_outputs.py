#!/usr/bin/env python3
"""Audit final CPM network-interaction tables and figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


OUTCOMES = ("GDE", "MEQ30", "VRS")
TOTAL_POSSIBLE = 86320


def bh_fdr(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p, kind="stable")
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.minimum(adjusted, 1.0)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir", type=Path, required=True
    )
    parser.add_argument("--expected-permutations", type=int, default=10000)
    parser.add_argument("--require-workbook", action="store_true")
    parser.add_argument("--diagnostics", action="store_true")
    return parser.parse_args()


def verify_png(path: Path) -> None:
    """Figures are laid out at true reproduction size, so check millimetres, not pixels.

    A correctly sized 180 x 64 mm composite is only 750 px tall; the meaningful
    journal criteria are the resolution and the printed width.
    """
    with Image.open(path) as image:
        dpi = image.info.get("dpi")
        if dpi is None or any(abs(float(value) - 300.0) > 0.02 for value in dpi):
            raise AssertionError(f"{path.name} is not tagged as 300 dpi: {dpi}")
        width_mm = image.width / float(dpi[0]) * 25.4
        height_mm = image.height / float(dpi[1]) * 25.4
        if width_mm < 80.0 or height_mm < 40.0:
            raise AssertionError(
                f"{path.name} is unexpectedly small: "
                f"{width_mm:.0f} x {height_mm:.0f} mm"
            )
        if width_mm > 190.0:
            raise AssertionError(
                f"{path.name} is wider than a double-column page: {width_mm:.0f} mm"
            )


def main() -> None:
    args = parse_args()
    output = args.output_dir.resolve()
    pairs = pd.read_csv(output / "network_pair_results.csv")
    summary = pd.read_csv(output / "outcome_summary.csv")
    null_counts = pd.read_csv(output / "null_consensus_edge_counts.csv")

    forbidden_ratio_columns = {
        "EmpiricalPermutationP",
        "FDRQ",
        "FDRSignificant",
        "ValidNullMasks",
        "ZeroEdgeNullMasksExcluded",
        "FoldEnrichment",
        "Log2FoldEnrichment",
    }
    unexpected = forbidden_ratio_columns.intersection(pairs.columns)
    if unexpected:
        raise AssertionError(f"Deprecated column names remain: {sorted(unexpected)}")

    if len(pairs) != 108:
        raise AssertionError(f"Expected 108 network-pair rows, found {len(pairs)}.")
    if len(summary) != 3:
        raise AssertionError(f"Expected 3 summary rows, found {len(summary)}.")
    if len(null_counts) != 3 * args.expected_permutations:
        raise AssertionError("Unexpected number of permutation-level edge counts.")

    required_summary = {"TotalMaskFDRQ", "TotalMaskFDRSignificant"}
    missing_summary = required_summary.difference(summary.columns)
    if missing_summary:
        raise AssertionError(
            f"Missing whole-mask FDR columns: {sorted(missing_summary)}"
        )
    expected_mask_q = bh_fdr(summary["TotalMaskPermutationP"].to_numpy(dtype=float))
    if not np.allclose(
        expected_mask_q,
        summary["TotalMaskFDRQ"].to_numpy(dtype=float),
        rtol=0,
        atol=1e-12,
    ):
        raise AssertionError("Whole-mask BH q-values are inconsistent.")
    if not np.array_equal(
        expected_mask_q < 0.05,
        summary["TotalMaskFDRSignificant"].to_numpy(dtype=bool),
    ):
        raise AssertionError("Whole-mask FDR significance flags are inconsistent.")

    for outcome in OUTCOMES:
        frame = pairs[pairs["Outcome"] == outcome].sort_values("PairIndex")
        if len(frame) != 36 or frame["PairIndex"].nunique() != 36:
            raise AssertionError(f"{outcome} does not contain 36 unique pairs.")
        if int(frame["PossibleEdges"].sum()) != TOTAL_POSSIBLE:
            raise AssertionError(f"{outcome} possible edges do not sum to 86,320.")
        observed = int(frame["ObservedConsensusEdges"].sum())
        summary_row = summary.loc[summary["Outcome"] == outcome].iloc[0]
        if observed != int(summary_row["ObservedPositiveConsensusEdges"]):
            raise AssertionError(f"{outcome} pair counts do not match its summary.")
        expected_density = (
            frame["ObservedConsensusEdges"].to_numpy(dtype=float)
            / frame["PossibleEdges"].to_numpy(dtype=float)
        )
        expected_ratio = expected_density / (observed / TOTAL_POSSIBLE) if observed else np.zeros_like(expected_density)
        if not np.allclose(
            expected_density,
            frame["NetworkPairDensity"].to_numpy(dtype=float),
            rtol=0,
            atol=1e-12,
        ):
            raise AssertionError(f"{outcome} pair densities are inconsistent.")
        if not np.allclose(
            expected_ratio,
            frame["ObservedExpectedRatio"].to_numpy(dtype=float),
            rtol=0,
            atol=1e-12,
        ):
            raise AssertionError(
                f"{outcome} observed/expected ratios are inconsistent."
            )
        expected_q = bh_fdr(frame["CountPermutationP"].to_numpy(dtype=float))
        if not np.allclose(
            expected_q,
            frame["CountFDRQ"].to_numpy(dtype=float),
            rtol=0,
            atol=1e-12,
        ):
            raise AssertionError(f"{outcome} BH q-values are inconsistent.")
        expected_sparse = frame["ObservedConsensusEdges"].to_numpy() < 5
        if not np.array_equal(expected_sparse, frame["SparseObservedCount"].to_numpy()):
            raise AssertionError(f"{outcome} sparse-cell flags are inconsistent.")
        outcome_null = null_counts[null_counts["Outcome"] == outcome]
        if len(outcome_null) != args.expected_permutations:
            raise AssertionError(f"{outcome} has the wrong null-count length.")
        null_values = outcome_null["PositiveConsensusEdgeCount"].to_numpy(dtype=float)
        expected_total_p = (1 + np.count_nonzero(null_values >= observed)) / (
            args.expected_permutations + 1
        )
        if not np.isclose(
            expected_total_p,
            float(summary_row["TotalMaskPermutationP"]),
            rtol=0,
            atol=1e-12,
        ):
            raise AssertionError(f"{outcome} total-mask p-value is inconsistent.")
        if int(summary_row["ZeroEdgeNullMasks"]) != int(np.count_nonzero(null_values == 0)):
            raise AssertionError(f"{outcome} zero-edge null-mask count is inconsistent.")
    # figures/ holds only the paper figures; supporting plots live in figures/other/.
    figure_dir = output / "figures"
    other_dir = figure_dir / "other"
    expected_stems = [(figure_dir, f"network_enrichment_{outcome}") for outcome in OUTCOMES]
    if args.diagnostics:
        expected_stems += [
            *((other_dir, f"network_enrichment_{outcome}_standalone") for outcome in OUTCOMES),
            *((other_dir, f"raw_consensus_edge_counts_{outcome}") for outcome in OUTCOMES),
            (other_dir, "null_positive_consensus_edge_distributions"),
            (other_dir, "network_enrichment_panel"),
        ]
    for directory, stem in expected_stems:
        png = directory / f"{stem}.png"
        svg = directory / f"{stem}.svg"
        pdf = directory / f"{stem}.pdf"
        if not png.is_file() or not svg.is_file() or not pdf.is_file():
            raise AssertionError(f"Missing PNG/SVG/PDF set for {stem}.")
        verify_png(png)

    stray = sorted(
        path.name
        for path in figure_dir.glob("*")
        if path.is_file()
        # Finder drops .DS_Store into any browsed folder; it is not an artifact.
        and not path.name.startswith(".")
        and not path.name.startswith("network_enrichment_")
    )
    if stray:
        raise AssertionError(
            f"figures/ must hold only the paper enrichment figures; found: {stray}"
        )
    if not (output / "figure_captions.md").is_file():
        raise AssertionError("figure_captions.md is missing.")

    if args.require_workbook and not (output / "network_interaction_results.xlsx").is_file():
        raise AssertionError("The final XLSX workbook is missing.")

    print("PASS: network-interaction outputs satisfy all audit checks.")


if __name__ == "__main__":
    main()
