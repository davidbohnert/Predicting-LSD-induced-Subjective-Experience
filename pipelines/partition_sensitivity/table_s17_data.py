"""Build descriptive edge and node rankings from saved selection frequencies.

Edges are ordered by consensus recurrence, fold-selection frequency, and mean
edge-behavior partial correlation. Nodes are ordered by stability-weighted degree.
Both rankings retain ties at fifth place. No feature selection is repeated here.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io as sio

OUTPUT = Path(__file__).resolve().parent / "output"
OUTCOMES = ["GDE", "MEQ30", "VRS"]
OUTCOME_LABELS = {
    "GDE": "VAS Good Drug Effect (GDE)",
    "MEQ30": "MEQ30 Total Score",
    "VRS": "5D-ASC Visionary Restructuralization (VRS)",
}
NO_NODES = 416
TOP_N = 5
MINUS = "−"


@dataclass(frozen=True)
class Parcel:
    parcel_id: int
    label: str
    network: str
    schaefer: str
    x: int
    y: int
    z: int
    anatomy: str
    anatomy_overlap: float

    @property
    def coordinates(self) -> str:
        return "({}, {}, {})".format(*(str(v).replace("-", MINUS)
                                       for v in (self.x, self.y, self.z)))

    @property
    def identity(self) -> str:
        """Parcel number and the Schaefer/AAL3 name it carries in the atlas.

        Cortical parcels use the network abbreviation of the figures rather
        than the spelled-out Schaefer family, so the cell stays narrow while
        remaining a unique key back into the atlas: "L SAL, ParOper 4" is
        7Networks_LH_SalVentAttn_ParOper_4.
        """
        if not self.schaefer:
            return f"{self.parcel_id}  {self.label}"
        parts = self.schaefer.split("_")
        side = "L" if parts[1] == "LH" else "R"
        index = parts[-1]
        region = parts[3] if len(parts) == 5 else ""
        suffix = f", {region} {index}" if region else f" {index}"
        return f"{self.parcel_id}  {side} {self.network}{suffix}"

    @property
    def location(self) -> str:
        anatomy = self.anatomy or "unlabeled"
        return f"{anatomy}  {self.coordinates}"


def load_parcels() -> dict[int, Parcel]:
    parcels = {}
    with (OUTPUT / "atlas_parcel_metadata.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            # The AAL descriptor already carries a hemisphere suffix; the
            # label ahead of it does too, so drop the duplicate.
            anatomy = row["AnatomicalRegion"]
            if anatomy.endswith((" L", " R")):
                anatomy = anatomy[:-2]
            parcel_id = int(row["ParcelID"])
            parcels[parcel_id] = Parcel(
                parcel_id=parcel_id,
                label=row["DisplayLabel"],
                network=row["Network"],
                schaefer=row["SchaeferName"],
                x=int(row["X"]), y=int(row["Y"]), z=int(row["Z"]),
                anatomy=anatomy,
                anatomy_overlap=float(row["AnatomicalOverlapPercent"]),
            )
    assert sorted(parcels) == list(range(1, NO_NODES + 1)), "expected parcels 1-416"
    return parcels


def edge_endpoints() -> tuple[np.ndarray, np.ndarray]:
    """Row and column parcel IDs for the 86,320 unique edges, in the
    column-major upper-triangle order the MATLAB layout uses."""
    linear = np.flatnonzero(np.triu(np.ones((NO_NODES, NO_NODES), bool), 1)
                            .ravel(order="F"))
    return linear % NO_NODES + 1, linear // NO_NODES + 1


def load_stability(outcome: str) -> dict:
    data = sio.loadmat(OUTPUT / f"edge_selection_frequencies_{outcome}.mat",
                       simplify_cells=True)
    fold = data["foldSelectionCount"].astype(np.int64)
    consensus = data["partitionConsensusCount"].astype(np.int64)
    fixed = data["fixedPartitionFoldCount"].astype(np.int64)
    mean_r = data["meanEdgeR"].astype(float)
    metadata = data["metadata"]
    no_folds = int(metadata["no_training_folds"])
    assert fold.size == consensus.size == fixed.size == mean_r.size == 86320
    assert fold.max() <= no_folds and consensus.max() <= int(metadata["no_partitions"])
    assert np.all(np.isfinite(mean_r)) and np.abs(mean_r).max() <= 1
    return {
        "fold": fold, "consensus": consensus, "fixed": fixed, "mean_r": mean_r,
        "no_training_folds": no_folds,
        "no_partitions": int(metadata["no_partitions"]),
        "no_folds": int(metadata["no_folds"]),
        "consensus_folds_required": int(metadata["consensus_folds_required"]),
        "n": int(metadata["no_subjects"]),
        "correlation": str(metadata["correlation_type"]),
        "fixed_edges": int(metadata["fixed_partition_edges"]),
    }


def take_with_ties(order: np.ndarray, keys: np.ndarray, limit: int) -> np.ndarray:
    """First LIMIT entries of ORDER, extended through every entry sharing the
    key vector of the last one, so a tie at fifth place is never broken by
    array order."""
    if order.size <= limit:
        return order
    boundary = keys[order[limit - 1]]
    keep = limit
    while keep < order.size and np.array_equal(keys[order[keep]], boundary):
        keep += 1
    return order[:keep]


def panel_a_rows(parcels: dict[int, Parcel]) -> list[dict]:
    rows_id, cols_id = edge_endpoints()
    out = []
    for outcome in OUTCOMES:
        stability = load_stability(outcome)
        consensus, fold, fixed = (stability["consensus"], stability["fold"],
                                  stability["fixed"])
        mean_r = stability["mean_r"]
        # Consensus recurrence first, then overall selection frequency, then
        # the strength of the association -- see the module docstring.
        rounded_r = np.round(mean_r, 6)
        order = np.lexsort((-rounded_r, -fold, -consensus))
        keys = np.stack([consensus, fold, rounded_r], axis=1)
        selected = take_with_ties(order, keys, TOP_N)

        rank, previous = 0, None
        for position, index in enumerate(selected, start=1):
            key = (consensus[index], fold[index], rounded_r[index])
            if key != previous:
                rank, previous = position, key
            first, second = parcels[rows_id[index]], parcels[cols_id[index]]
            pair = "–".join(sorted({first.network, second.network})) \
                if first.network != second.network else first.network
            out.append({
                "Outcome": outcome,
                "Rank": rank,
                "ROI1_ID": first.parcel_id, "ROI1": first,
                "ROI2_ID": second.parcel_id, "ROI2": second,
                "NetworkPair": pair,
                "FixedPartitionFolds": int(fixed[index]),
                "NoFolds": stability["no_folds"],
                "ConsensusRecurrence": int(consensus[index]),
                "NoPartitions": stability["no_partitions"],
                "SelectionFrequency": 100.0 * fold[index] / stability["no_training_folds"],
                "FoldSelectionCount": int(fold[index]),
                "MeanPartialR": float(mean_r[index]),
                "MaximallyStableEdges": int(np.sum((consensus == stability["no_partitions"])
                                                   & (fold == stability["no_training_folds"]))),
                "EdgeIndex": int(index) + 1,
            })
    return out


def panel_b_rows(parcels: dict[int, Parcel]) -> list[dict]:
    rows_id, cols_id = edge_endpoints()
    out = []
    for outcome in OUTCOMES:
        stability = load_stability(outcome)
        frequency = stability["fold"] / stability["no_training_folds"]
        weighted = np.zeros(NO_NODES + 1)
        np.add.at(weighted, rows_id, frequency)
        np.add.at(weighted, cols_id, frequency)

        in_fixed = stability["fixed"] >= stability["consensus_folds_required"]
        fixed_degree = np.zeros(NO_NODES + 1, dtype=int)
        np.add.at(fixed_degree, rows_id[in_fixed], 1)
        np.add.at(fixed_degree, cols_id[in_fixed], 1)

        total = weighted.sum()
        assert np.isclose(total, 2 * frequency.sum()), \
            "weighted degrees must sum to twice the total edge frequency"
        assert fixed_degree.sum() == 2 * stability["fixed_edges"], \
            "fixed-mask degrees must sum to twice the consensus-edge count"

        nodes = np.arange(1, NO_NODES + 1)
        order = nodes[np.argsort(-weighted[1:], kind="stable")]
        # Weighted degrees are continuous; ties are possible only at zero.
        selected = take_with_ties(order - 1, np.round(weighted[1:], 9), TOP_N) + 1

        rank, previous = 0, None
        for position, node in enumerate(selected, start=1):
            key = round(float(weighted[node]), 9)
            if key != previous:
                rank, previous = position, key
            out.append({
                "Outcome": outcome,
                "Rank": rank,
                "ParcelID": int(node),
                "Parcel": parcels[int(node)],
                "Network": parcels[int(node)].network,
                "FixedMaskDegree": int(fixed_degree[node]),
                "WeightedDegree": float(weighted[node]),
                "EndpointShare": 100.0 * weighted[node] / total,
            })
    return out


def load_panels() -> tuple[list[dict], list[dict], dict[str, dict]]:
    parcels = load_parcels()
    stability = {outcome: load_stability(outcome) for outcome in OUTCOMES}
    return panel_a_rows(parcels), panel_b_rows(parcels), stability
