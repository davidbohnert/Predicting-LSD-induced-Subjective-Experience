#!/usr/bin/env python3
"""Create paper-facing CSV and Markdown versions of updated Table S16."""

from __future__ import annotations

import csv
import argparse
from pathlib import Path

EM_DASH = "—"


def p4(value):
    if value in (None, ""):
        return EM_DASH
    return f"{float(value):.4f}".lstrip("0")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.input.open(newline="", encoding="utf-8-sig") as handle:
        source = list(csv.DictReader(handle))
    if len(source) != 6:
        raise RuntimeError("canonical analysis must contain six rows")
    rows = []
    for row in source:
        primary = float(row["EdgeThreshold"]) == 0.01
        rows.append({
            "Feature threshold": f"p < {float(row['EdgeThreshold']):.2f}".replace("0.", "."),
            "Outcome pair": f"{row['OutcomeA']}–{row['OutcomeB']}",
            "Mask edges, K1/K2": f"{int(row['MaskSizeA']):,} / {int(row['MaskSizeB']):,}",
            "Shared edges (expected)": (
                f"{int(row['SharedEdges']):,} ({float(row['ExpectedSharedEdges']):.2f})"),
            "Observed/expected": f"{float(row['ObservedExpectedRatio']):.2f}",
            "Dice": f"{float(row['Dice']):.3f}".lstrip("0"),
            "Permutation p": p4(row["PermutationP"] if primary else None),
            "BH q": p4(row["BH_Q_Primary"] if primary else None),
        })
    fields = list(rows[0])
    csv_path = args.output / "edge_overlap_joint_native.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    md_path = args.output / "edge_overlap_joint_native.md"
    with md_path.open("w", encoding="utf-8") as handle:
        handle.write("# Table S16. Overlap between positive CPM consensus masks\n\n")
        handle.write("| " + " | ".join(fields) + " |\n")
        handle.write("| " + " | ".join(["---"] * len(fields)) + " |\n")
        for row in rows:
            handle.write("| " + " | ".join(row[field] for field in fields) + " |\n")
        handle.write("\nPermutation p values use 10,000 synchronized Freedman–Lane "
                     "permutations of the complete CPM procedure with native masks; "
                     "BH correction covers the three p < .01 comparisons. The p < .05 "
                     "block is descriptive. Dice is an effect size.\n")
    print(csv_path)
    print(md_path)


if __name__ == "__main__":
    main()
