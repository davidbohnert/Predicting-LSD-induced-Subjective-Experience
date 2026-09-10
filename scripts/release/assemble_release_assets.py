#!/usr/bin/env python3
"""Assemble only nonidentifying, verified derived outputs for a release."""
from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path


ASSETS = {
    "internal_validation/full_internal_grid/internal_validation_results.csv": "source_data/internal_full_internal_grid.csv",
    "internal_validation/fold_sensitivity/internal_validation_results.csv": "source_data/internal_fold_sensitivity.csv",
    "internal_validation/threshold_sensitivity/internal_validation_results.csv": "source_data/internal_threshold_sensitivity.csv",
    "internal_validation/alternative_correlation/internal_validation_results.csv": "source_data/internal_alternative_correlation.csv",
    "internal_validation/gsr/internal_validation_results.csv": "source_data/internal_gsr.csv",
    "internal_validation/matched_sample/internal_validation_results.csv": "source_data/internal_matched_sample.csv",
    "internal_validation/study_loso/internal_validation_results.csv": "source_data/internal_study_loso.csv",
    "internal_validation/primary/table_1_internal_validation.csv": "tables/table_1.csv",
    "internal_validation/primary/internal_validation_results.csv": "source_data/internal_primary.csv",
    "cross_drug_loso/table_2_cross_drug.csv": "tables/table_2.csv",
    "cross_drug_loso/cross_drug_loso_results.csv": "source_data/cross_drug_all_thresholds.csv",
    "partition_sensitivity/partition_sensitivity_results.csv": "source_data/figure_S1.csv",
    "learning_curves/output/primary_learning_curve_summary.csv": "source_data/figure_S2.csv",
    "parcellation_sensitivity/table_S13_internal.csv": "tables/table_S13_internal.csv",
    "parcellation_sensitivity/table_S13_cross_drug.csv": "tables/table_S13_cross_drug.csv",
    "network_interactions/outcome_summary.csv": "tables/table_S15_outcomes.csv",
    "network_interactions/network_pair_results.csv": "tables/table_S15_network_pairs.csv",
    "network_overlap/joint_edge_overlap_results.csv": "tables/table_S16.csv",
    "partition_sensitivity/stability/table_S17_panel_A_top_edges.csv": "tables/table_S17_edges.csv",
    "partition_sensitivity/stability/table_S17_panel_B_top_nodes.csv": "tables/table_S17_nodes.csv",
    "partition_sensitivity/figures/partition_sensitivity_primary.png": "figures/figure_S1.png",
    "partition_sensitivity/figures/partition_sensitivity_primary.pdf": "figures/figure_S1.pdf",
    "learning_curves/output/figure_s2_learning_curves.png": "figures/figure_S2.png",
    "learning_curves/output/figure_s2_learning_curves.pdf": "figures/figure_S2.pdf",
    "network_interactions/figures/network_enrichment_GDE.png": "figures/network_interaction_GDE.png",
    "network_interactions/figures/network_enrichment_GDE.svg": "figures/network_interaction_GDE.svg",
    "network_interactions/figures/network_enrichment_GDE.pdf": "figures/network_interaction_GDE.pdf",
    "network_interactions/figures/network_enrichment_MEQ30.png": "figures/network_interaction_MEQ30.png",
    "network_interactions/figures/network_enrichment_MEQ30.svg": "figures/network_interaction_MEQ30.svg",
    "network_interactions/figures/network_enrichment_MEQ30.pdf": "figures/network_interaction_MEQ30.pdf",
    "network_interactions/figures/network_enrichment_VRS.png": "figures/network_interaction_VRS.png",
    "network_interactions/figures/network_enrichment_VRS.svg": "figures/network_interaction_VRS.svg",
    "network_interactions/figures/network_enrichment_VRS.pdf": "figures/network_interaction_VRS.pdf",
}

# Older completed runs used the previous supplement numbering.
LEGACY_SOURCES = {
    "parcellation_sensitivity/table_S13_internal.csv":
        "parcellation_sensitivity/table_S12_internal.csv",
    "parcellation_sensitivity/table_S13_cross_drug.csv":
        "parcellation_sensitivity/table_S12_cross_drug.csv",
}


def resolve_source(root: Path, relative_source: str) -> Path:
    """Prefer current filenames; accept legacy names without changing results."""
    source = root / relative_source
    legacy_name = LEGACY_SOURCES.get(relative_source)
    if legacy_name is None:
        return source
    legacy = root / legacy_name
    if source.is_file() and legacy.is_file():
        if source.read_bytes() != legacy.read_bytes():
            raise SystemExit(
                f"Conflicting current and legacy table files: {source}, {legacy}. "
                "Stage only the intended final table before assembly."
            )
    return source if source.is_file() else legacy


REPOSITORY_ASSETS = {
    "metadata/atlas/Schaefer416_8networks.txt":
        "metadata/atlas/Schaefer416_8networks.txt",
    "metadata/atlas/atlas_parcel_metadata.csv":
        "metadata/atlas/atlas_parcel_metadata.csv",
    "metadata/atlas/Schaefer439_structure_membership.csv":
        "metadata/atlas/Schaefer439_structure_membership.csv",
    "metadata/atlas/SOURCES.md": "metadata/atlas/SOURCES.md",
    "metadata/cross_drug_row_map.csv": "metadata/cross_drug_row_map.csv",
}


def run_directory(root: Path, relative_source: str) -> Path:
    """Locate the owning analysis record, not a figure/table subdirectory."""
    source = root / relative_source
    for directory in (source.parent, *source.parent.parents):
        if directory == root.parent:
            break
        if list(directory.glob("run_record_*.txt")):
            return directory
    return source.parent


def completed_record(directory: Path) -> bool:
    records = sorted(directory.glob("run_record_*.txt"))
    return bool(records) and "status=COMPLETE" in records[-1].read_text().splitlines()


def copy_asset(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def write_manifest(directory: Path) -> None:
    lines = []
    for path in sorted(item for item in directory.rglob("*") if item.is_file()):
        if path.name == "MANIFEST.sha256":
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {path.relative_to(directory).as_posix()}")
    (directory / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    args = parser.parse_args()
    root = args.analysis_root.resolve()
    destination = args.destination.resolve()
    repository = Path(__file__).resolve().parents[2]

    if destination.exists() and any(destination.iterdir()):
        raise SystemExit(
            f"Destination must be absent or empty: {destination}"
        )

    resolved = {source: resolve_source(root, source) for source in ASSETS}
    missing = [source for source, path in resolved.items() if not path.is_file()]
    if missing:
        raise SystemExit("Missing release inputs:\n  " + "\n  ".join(missing))
    run_directories = {
        run_directory(root, source) for source in ASSETS
        if not source.startswith("parcellation_sensitivity/")
    }
    # The S13 builder combines six separately recorded analyses.
    for atlas in ("original", "no_cerebellum", "with_cerebellum"):
        for pipeline in ("internal_validation", "cross_drug_loso"):
            run_directories.add(root / "parcellation_sensitivity" / atlas / pipeline)
    incomplete = [str(path) for path in sorted(run_directories)
                  if not completed_record(path)
                  or not (path / "checkpoint_signature.mat").is_file()]
    if incomplete:
        raise SystemExit("Outputs lack a COMPLETE run record or compatibility signature:\n  " + "\n  ".join(incomplete))

    mask_source = root / "internal_validation/primary/consensus_masks"
    masks = [mask_source / f"positive_consensus_{outcome}.mat"
             for outcome in ("BDE", "GDE", "OBN", "AED", "VRS", "MEQ30")]
    missing_masks = [str(path) for path in masks if not path.is_file()]
    if missing_masks:
        raise SystemExit("Missing aggregate masks:\n  " + "\n  ".join(missing_masks))

    for source, relative_destination in ASSETS.items():
        copy_asset(resolved[source], destination / relative_destination)

    for source, relative_destination in REPOSITORY_ASSETS.items():
        source_path = repository / source
        if not source_path.is_file():
            raise SystemExit(f"Missing repository release input: {source}")
        copy_asset(source_path, destination / relative_destination)

    for source in masks:
        copy_asset(source, destination / "masks" / source.name)

    write_manifest(destination)
    print(f"Release assets assembled in {destination}")


if __name__ == "__main__":
    main()
