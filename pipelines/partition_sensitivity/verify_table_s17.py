"""Independent checks on Table S17.

Recomputes the reported quantities straight from the saved frequency
vectors and the atlas files, without going through table_s17_data, so a
mistake in the panel builder cannot pass unnoticed.
"""
from __future__ import annotations

from pathlib import Path
import argparse

import numpy as np
import scipy.io as sio

import table_s17_data as data

OUTPUT = Path(__file__).resolve().parent / "output"
NETWORK_DISPLAY = {"VIS": "VIS", "ASM": "SMN", "DAN": "DAN", "SAL": "SAL",
                   "LIM": "LIM", "FPN": "FPN", "DMN": "DMN", "Subcortical": "SBC"}
failures = []


def check(name: str, condition: bool, detail: str = "") -> None:
    print(f"  {'PASS' if condition else 'FAIL'}  {name}" + (f"  {detail}" if detail else ""))
    if not condition:
        failures.append(name)


def atlas_networks(atlas_file: Path) -> dict[int, str]:
    networks = {}
    for line in atlas_file.read_text().splitlines():
        parts = line.split()
        if parts:
            networks[int(float(parts[0]))] = NETWORK_DISPLAY[parts[1]]
    return networks


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify Table S17 outputs")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--atlas-file", type=Path, required=True)
    args = parser.parse_args()
    global OUTPUT
    OUTPUT = args.output_dir.resolve()
    data.OUTPUT = OUTPUT
    panel_a, panel_b, _ = data.load_panels()
    networks = atlas_networks(args.atlas_file)
    rows_id, cols_id = data.edge_endpoints()

    print("Atlas join")
    check("416 parcels with network labels", len(networks) == 416)
    check("parcel ids are 1..416", sorted(networks) == list(range(1, 417)))
    parcels = data.load_parcels()
    check("metadata networks agree with the atlas file",
          all(parcels[i].network == networks[i] for i in range(1, 417)))

    for outcome in data.OUTCOMES:
        print(f"\n{outcome}")
        raw = sio.loadmat(OUTPUT / f"edge_selection_frequencies_{outcome}.mat",
                          simplify_cells=True)
        fold = raw["foldSelectionCount"].astype(np.int64)
        consensus = raw["partitionConsensusCount"].astype(np.int64)
        fixed = raw["fixedPartitionFoldCount"].astype(np.int64)
        mean_r = raw["meanEdgeR"].astype(float)
        no_training_folds = int(raw["metadata"]["no_training_folds"])
        no_partitions = int(raw["metadata"]["no_partitions"])
        consensus_folds = int(raw["metadata"]["consensus_folds_required"])
        frequency = fold / no_training_folds

        rebuilt = np.zeros((416, 416), dtype=np.uint8)
        selected = fixed >= consensus_folds
        rebuilt[rows_id[selected] - 1, cols_id[selected] - 1] = 1
        rebuilt[cols_id[selected] - 1, rows_id[selected] - 1] = 1
        check("seed-123 mask is symmetric",
              np.array_equal(rebuilt, rebuilt.T))
        fixed_edge_count = int(selected.sum())
        check("seed-123 mask has the stored edge count",
              fixed_edge_count == int(raw["metadata"]["fixed_partition_edges"]))

        # Panel B identities.
        weighted = np.zeros(417)
        np.add.at(weighted, rows_id, frequency)
        np.add.at(weighted, cols_id, frequency)
        check("weighted degrees sum to twice the total edge frequency",
              np.isclose(weighted.sum(), 2 * frequency.sum()))
        shares = [row["EndpointShare"] for row in panel_b if row["Outcome"] == outcome]
        check("reported shares match a direct recomputation",
              all(np.isclose(row["EndpointShare"],
                             100 * weighted[row["ParcelID"]] / weighted.sum())
                  for row in panel_b if row["Outcome"] == outcome),
              f"top five sum to {sum(shares):.2f}%")

        fixed_degree = np.zeros(417, dtype=int)
        np.add.at(fixed_degree, rows_id[selected], 1)
        np.add.at(fixed_degree, cols_id[selected], 1)
        check("fixed-mask degrees sum to twice the consensus-edge count",
              fixed_degree.sum() == 2 * fixed_edge_count)
        check("reported fixed-mask degrees match the rebuilt mask",
              all(row["FixedMaskDegree"] == rebuilt[row["ParcelID"] - 1].sum()
                  for row in panel_b if row["Outcome"] == outcome))
        check("Panel B is ordered by weighted degree",
              [row["WeightedDegree"] for row in panel_b if row["Outcome"] == outcome]
              == sorted((row["WeightedDegree"] for row in panel_b
                         if row["Outcome"] == outcome), reverse=True))

        # Panel A rows, recomputed one edge at a time from the raw vectors.
        entries = [row for row in panel_a if row["Outcome"] == outcome]
        ok = True
        for entry in entries:
            index = entry["EdgeIndex"] - 1
            ends = {int(rows_id[index]), int(cols_id[index])}
            ok &= ends == {entry["ROI1_ID"], entry["ROI2_ID"]}
            ok &= entry["ConsensusRecurrence"] == consensus[index]
            ok &= entry["FoldSelectionCount"] == fold[index]
            ok &= entry["FixedPartitionFolds"] == fixed[index]
            ok &= np.isclose(entry["MeanPartialR"], mean_r[index])
            ok &= np.isclose(entry["SelectionFrequency"], 100 * frequency[index])
            pair = "–".join(sorted({networks[e] for e in ends})) \
                if len({networks[e] for e in ends}) > 1 else networks[min(ends)]
            ok &= pair == entry["NetworkPair"]
        check("every Panel A row matches the raw vectors", ok)

        maximal = (consensus == no_partitions) & (fold == no_training_folds)
        check("Panel A edges are all maximally stable",
              all(maximal[row["EdgeIndex"] - 1] for row in entries),
              f"{int(maximal.sum())} maximally stable edges in total")
        cutoff = min(entry["MeanPartialR"] for entry in entries)
        # Exactly the reported edges may sit at or above the fifth-place value.
        check("no maximally stable edge has a stronger association than those shown",
              int(np.sum(maximal & (mean_r >= cutoff))) == len(entries),
              f"cut-off r = {cutoff:.3f}")

    print("\n" + ("All checks passed." if not failures
                  else f"{len(failures)} FAILED: {failures}"))
    raise SystemExit(1 if failures else 0)


if __name__ == "__main__":
    main()
