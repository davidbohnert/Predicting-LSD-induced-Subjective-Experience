"""Render prediction accuracy across 100 partitions of the primary model.

Panels show all six outcomes. Bars are grouped by each partition's uncorrected
permutation p-value; the paper's partition is marked. Read the consolidated CSV
or the earlier aggregate workbook. See README.md for the command.
"""

import argparse
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# --- what this figure reports -------------------------------------------------

MODEL_FORM = "LSD-PCB difference"
COVARIATES = "With covariates"
NETWORK = "Positive"

# Panel order follows Table 1 of Results_Tables_final, which orders the six
# dimensions by effect size.
SCALES = [
    ("GDE", "Good Drug Effect"),
    ("MEQ30", "Mystical Experience"),
    ("VRS", "Visionary Restructuralization"),
    ("OBN", "Oceanic Boundlessness"),
    ("BDE", "Bad Drug Effect"),
    ("AED", "Anxious Ego Dissolution"),
]

# --- appearance ---------------------------------------------------------------

# Pure white: this figure is embedded in a white document page, where an
# off-white surface reads as a visible grey box. Validated for the ordinal
# ramp against #ffffff (light end 2.11:1).
SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#dedcd6"

# Ordinal ramp: one hue, light->dark, darkest = most significant.
# Steps 250/400/600 of the reference blue ramp; validated for ordinal use on
# this surface (monotone lightness, adjacent dL >= 0.06, light end 2.06:1).
C_NS = "#86b6ef"
C_P10 = "#3987e5"
C_P05 = "#184f95"

ACCENT_LINE = "#d95926"  # seed-123 rule, 3.78:1 on surface

X_MIN, X_MAX, BIN_W = -0.30, 0.35, 0.01
Y_MAX, Y_TICKS = 24, [0, 8, 16, 24]

FIG_W_IN, FIG_H_IN = 7.09, 5.15  # 180 mm full journal page width


def load_partitions(workbook_path):
    """Return {scale: (seeds, r, p)} for the primary specification."""
    book = openpyxl.load_workbook(workbook_path, data_only=True)

    rows = list(book["Partition data"].iter_rows(values_only=True))
    header = next(i for i, row in enumerate(rows) if row[0] == "Outcome")
    col = {name: i for i, name in enumerate(rows[header]) if name}

    data = {}
    for row in rows[header + 1:]:
        if not row[col["Outcome"]]:
            continue
        if row[col["Model form"]] != MODEL_FORM:
            continue
        if row[col["Covariates"]] != COVARIATES:
            continue
        data.setdefault(row[col["Outcome"]], []).append((
            row[col["CV seed"]],
            row[col[f"{NETWORK} r"]],
            row[col[f"{NETWORK} p"]],
        ))

    out = {}
    for scale, _ in SCALES:
        assert scale in data, f"{scale} missing from the workbook"
        entries = sorted(data[scale])
        assert len(entries) == 100, f"{scale}: expected 100 partitions, got {len(entries)}"
        seeds = np.array([e[0] for e in entries])
        assert 123 in seeds, f"{scale}: seed 123 (the manuscript partition) is absent"
        out[scale] = (
            seeds,
            np.array([e[1] for e in entries], dtype=float),
            np.array([e[2] for e in entries], dtype=float),
        )
    return out, book


def load_partition_csv(csv_path, expected_partitions=100):
    """Return publication-format partition rows from the consolidated CSV."""
    frame = pd.read_csv(csv_path)
    required = {"Outcome", "CVSeed", "PositiveR", "PositiveP"}
    missing = required.difference(frame.columns)
    assert not missing, f"missing CSV columns: {sorted(missing)}"
    data = {}
    for scale, _ in SCALES:
        rows = frame.loc[frame["Outcome"] == scale].sort_values("CVSeed")
        assert len(rows) == expected_partitions, (
            f"{scale}: expected {expected_partitions} partitions, got {len(rows)}"
        )
        seeds = rows["CVSeed"].to_numpy(dtype=int)
        assert 123 in seeds, f"{scale}: seed 123 is absent"
        data[scale] = (
            seeds,
            rows["PositiveR"].to_numpy(dtype=float),
            rows["PositiveP"].to_numpy(dtype=float),
        )
    return data


def check_against_summary(book, data):
    """Assert our recomputed percentages match the workbook's Summary sheet."""
    rows = list(book["Summary"].iter_rows(values_only=True))
    header = next(i for i, row in enumerate(rows) if row[0] == "Outcome")
    col = {name: i for i, name in enumerate(rows[header]) if name}

    seen = 0
    for row in rows[header + 1:]:
        if not row[col["Outcome"]]:
            continue
        if (row[col["Model form"]], row[col["Covariates"]], row[col["Network"]]) != (
            MODEL_FORM, COVARIATES, NETWORK
        ):
            continue
        scale = row[col["Outcome"]]
        if scale not in data:
            continue
        _, _, p = data[scale]
        for thr, label in ((0.05, "Partitions with raw p<0.05 (%)"),
                           (0.10, "Partitions with raw p<0.10 (%)")):
            ours = float((p < thr).sum())
            theirs = float(row[col[label]])
            assert abs(ours - theirs) < 1e-9, (
                f"{scale} p<{thr}: figure says {ours}%, Summary sheet says {theirs}%"
            )
        seen += 1

    assert seen == len(SCALES), f"matched {seen} Summary rows, expected {len(SCALES)}"


def draw(data, out_dir):
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "pdf.fonttype": 42,
    })

    edges = np.arange(X_MIN, X_MAX + BIN_W / 2, BIN_W)
    centers = edges[:-1] + BIN_W / 2

    fig, axes = plt.subplots(
        2, 3, figsize=(FIG_W_IN, FIG_H_IN), sharex=True, sharey=True,
        facecolor=SURFACE,
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.915, bottom=0.155,
                        wspace=0.16, hspace=0.78)

    for idx, ((scale, full_name), ax) in enumerate(zip(SCALES, axes.ravel())):
        seeds, r, p = data[scale]
        ax.set_facecolor(SURFACE)

        # Stack each bin by the partitions' own p-values.
        bands = [
            (r[p >= 0.10], C_NS),
            (r[(p >= 0.05) & (p < 0.10)], C_P10),
            (r[p < 0.05], C_P05),
        ]
        bottom = np.zeros(len(centers))
        for values, color in bands:
            counts, _ = np.histogram(values, bins=edges)
            ax.bar(centers, counts, bottom=bottom, width=BIN_W,
                   color=color, edgecolor=SURFACE, linewidth=0.35, zorder=3)
            bottom += counts

        assert bottom.sum() == len(seeds), (
            f"{scale}: {bottom.sum():.0f} partitions binned, expected {len(seeds)}"
        )
        assert bottom.max() <= Y_MAX, f"{scale}: peak bin {bottom.max():.0f} exceeds y-limit {Y_MAX}"

        # r = 0 reference, then the manuscript's partition.
        ax.axvline(0, color=GRID, linewidth=0.7, zorder=1)

        # The rule is identified once in the shared legend, not labelled in
        # every panel: an in-panel label collides with the tall bars it sits
        # next to, since the manuscript partition tends to fall near the mode.
        seed_r = float(r[seeds == 123][0])
        ax.axvline(seed_r, color=ACCENT_LINE, linewidth=1.0,
                   linestyle=(0, (3, 2)), zorder=4)

        # Percentages sit above the axes, clear of the distribution.
        n_p05 = int((p < 0.05).sum())
        n_p10 = int((p < 0.10).sum())
        ax.set_title(f"{scale} — {full_name}", fontsize=8.5, color=INK,
                     pad=17, loc="left")
        ax.text(0.0, 1.025, f"p < .05: {n_p05}%     p < .10: {n_p10}%",
                transform=ax.transAxes, ha="left", va="bottom",
                fontsize=6.8, color=INK_SECONDARY)

        ax.set_xlim(X_MIN, X_MAX)
        ax.set_ylim(0, Y_MAX)
        ax.set_yticks(Y_TICKS)
        ax.set_xticks(np.arange(-0.3, 0.31, 0.1))
        ax.tick_params(labelsize=7.5, colors=INK_SECONDARY)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(INK_SECONDARY)

        # Both rows carry their own tick labels and axis title, so neither row
        # has to be read against the other row's axis.
        ax.tick_params(labelbottom=True)
        ax.set_xlabel("Prediction accuracy (r)", fontsize=8, color=INK)
        if idx % 3 == 0:
            total_partitions = len(data[SCALES[0][0]][0])
            ax.set_ylabel(f"Partitions (of {total_partitions})", fontsize=8, color=INK)

    handles = [
        Patch(facecolor=C_P05, edgecolor=SURFACE, linewidth=0.35, label="p < .05"),
        Patch(facecolor=C_P10, edgecolor=SURFACE, linewidth=0.35, label=".05 ≤ p < .10"),
        Patch(facecolor=C_NS, edgecolor=SURFACE, linewidth=0.35, label="p ≥ .10"),
        Line2D([0], [0], color=ACCENT_LINE, linewidth=1.0, linestyle=(0, (3, 2)),
               label="Partition used in the manuscript"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=7.5, labelcolor=INK, handlelength=1.6, handleheight=1.0,
               columnspacing=1.8, bbox_to_anchor=(0.5, 0.012))

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "partition_sensitivity_primary"
    fig.savefig(f"{stem}.png", dpi=600, facecolor=SURFACE)
    fig.savefig(f"{stem}.pdf", facecolor=SURFACE)
    plt.close(fig)
    return stem


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path,
                        help="Partition CSV or canonical XLSX workbook")
    parser.add_argument("--output", required=True, type=Path,
                        help="Directory for PNG and PDF output")
    parser.add_argument("--expected-partitions", type=int, default=100,
                        help="Expected rows per outcome; use 1 for a smoke test")
    args = parser.parse_args()
    if args.input.suffix.lower() == ".csv":
        data = load_partition_csv(args.input, args.expected_partitions)
        verification = "verified consolidated CSV"
    else:
        data, book = load_partitions(args.input)
        check_against_summary(book, data)
        verification = "verified against the Summary sheet"

    stem = draw(data, args.output)

    print(f"{MODEL_FORM}, {NETWORK.lower()} network, {COVARIATES.lower()}\n")
    for scale, _ in SCALES:
        seeds, r, p = data[scale]
        seed_r = float(r[seeds == 123][0])
        print(f"  {scale:<6} median r = {np.median(r):+.3f}   "
              f"range [{r.min():+.3f}, {r.max():+.3f}]   "
              f"seed 123 r = {seed_r:+.3f}   "
              f"p<.05 {int((p < 0.05).sum()):>3}%   p<.10 {int((p < 0.10).sum()):>3}%")
    print(f"\n{verification}\nwrote {stem}.png and {stem}.pdf")


if __name__ == "__main__":
    main()
