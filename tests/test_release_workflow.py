"""Release completeness and result-dependent caption checks."""
from __future__ import annotations

import contextlib
import importlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/release"))
sys.path.insert(0, str(ROOT / "pipelines/network_interactions"))
import assemble_release_assets as release
import render_network_interaction_figures as figures
import verify_network_interaction_outputs as verifier


def make_network_inputs(output: Path):
    """Synthetic 416-node network counts; no participant data."""
    sizes = [50] * 7 + [66]
    all_pairs, summaries, nulls = [], [], []
    for outcome in figures.OUTCOMES:
        rows = []
        for i, first in enumerate(figures.NETWORKS):
            for j in range(i, len(figures.NETWORKS)):
                count = (i + j) % 5
                possible = sizes[i] * sizes[j] if i != j else sizes[i] * (sizes[i] - 1) // 2
                rows.append(dict(Outcome=outcome, PairIndex=len(rows)+1,
                                 Network1=first, Network2=figures.NETWORKS[j],
                                 NetworkPair=f"{first}–{figures.NETWORKS[j]}",
                                 ObservedConsensusEdges=count, PossibleEdges=possible,
                                 NetworkPairDensity=count/possible, CountPermutationP=.5,
                                 CountFDRQ=.5, SparseObservedCount=count < 5))
        total = sum(row["ObservedConsensusEdges"] for row in rows)
        for row in rows:
            row["ObservedExpectedRatio"] = row["NetworkPairDensity"] / (total/86320)
        all_pairs.extend(rows)
        values = np.arange(total - 5, total + 5)
        p = (1 + np.count_nonzero(values >= total)) / 11
        summaries.append(dict(Outcome=outcome, ObservedPositiveConsensusEdges=total,
                              TotalMaskPermutationP=p, TotalMaskFDRQ=p,
                              TotalMaskFDRSignificant=False, ZeroEdgeNullMasks=0))
        nulls.extend(dict(Outcome=outcome, Iteration=i+1, PositiveConsensusEdgeCount=int(value))
                     for i, value in enumerate(values))
    pd.DataFrame(all_pairs).to_csv(output / "network_pair_results.csv", index=False)
    pd.DataFrame(summaries).to_csv(output / "outcome_summary.csv", index=False)
    pd.DataFrame(nulls).to_csv(output / "null_consensus_edge_counts.csv", index=False)


class ReleaseWorkflowTests(unittest.TestCase):
    def test_figure_modes_and_verification(self):
        with tempfile.TemporaryDirectory() as name:
            output = Path(name)
            make_network_inputs(output)
            args = ["render", "--output-dir", str(output)]
            with patch.object(sys, "argv", args), contextlib.redirect_stdout(io.StringIO()):
                figures.main()
            self.assertEqual(len(list((output / "figures").glob("*.png"))), 3)
            self.assertFalse((output / "figures/other").exists())
            check = ["verify", "--output-dir", str(output), "--expected-permutations", "10"]
            with patch.object(sys, "argv", check), contextlib.redirect_stdout(io.StringIO()):
                verifier.main()
            with patch.object(sys, "argv", args + ["--diagnostics"]), contextlib.redirect_stdout(io.StringIO()):
                figures.main()
            with patch.object(sys, "argv", check + ["--diagnostics"]), contextlib.redirect_stdout(io.StringIO()):
                verifier.main()
    def test_supplementary_coverage_and_complete_assembly(self):
        for preset in ("full_internal_grid", "fold_sensitivity",
                       "threshold_sensitivity", "alternative_correlation",
                       "gsr", "matched_sample", "study_loso"):
            self.assertIn(f"internal_validation/{preset}/internal_validation_results.csv",
                          release.ASSETS)
        with tempfile.TemporaryDirectory() as name:
            root = Path(name) / "analyses"
            destination = Path(name) / "release"
            for relative in release.ASSETS:
                source = root / relative
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_bytes(b"synthetic aggregate")
            owners = {root / "internal_validation" / preset for preset in
                      ("primary", "full_internal_grid", "fold_sensitivity",
                       "threshold_sensitivity", "alternative_correlation", "gsr",
                       "matched_sample", "study_loso")}
            owners.update(root / relative for relative in
                          ("cross_drug_loso", "partition_sensitivity",
                           "partition_sensitivity/stability", "learning_curves/output",
                           "network_interactions", "network_overlap"))
            owners.update(root / "parcellation_sensitivity" / atlas / pipeline
                          for atlas in ("original", "no_cerebellum", "with_cerebellum")
                          for pipeline in ("internal_validation", "cross_drug_loso"))
            for owner in owners:
                owner.mkdir(parents=True, exist_ok=True)
                (owner / "run_record_001.txt").write_text("status=COMPLETE\n")
                (owner / "checkpoint_signature.mat").write_bytes(b"synthetic")
            masks = root / "internal_validation/primary/consensus_masks"
            masks.mkdir()
            for outcome in ("BDE", "GDE", "OBN", "AED", "VRS", "MEQ30"):
                (masks / f"positive_consensus_{outcome}.mat").write_bytes(b"synthetic mask")
            (masks / "unexpected_participant_data.mat").write_bytes(b"not allowed")
            argv = ["assemble", "--analysis-root", str(root),
                    "--destination", str(destination)]
            # An incomplete contributing run must block assembly.
            record = root / "internal_validation/gsr/run_record_001.txt"
            record.write_text("status=RUNNING\n")
            with patch.object(sys, "argv", argv), self.assertRaises(SystemExit):
                release.main()
            self.assertFalse(destination.exists())
            record.write_text("status=COMPLETE\n")
            with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()):
                release.main()
            self.assertEqual(len(list((destination / "masks").glob("*.mat"))), 6)
            self.assertFalse((destination / "masks/unexpected_participant_data.mat").exists())
            self.assertTrue((destination / "source_data/internal_gsr.csv").is_file())

    def test_python_dependencies(self):
        for module in ("numpy", "scipy.io", "pandas", "matplotlib", "PIL", "openpyxl"):
            importlib.import_module(module)

    def test_captions_follow_results(self):
        rows = []
        for outcome in figures.OUTCOMES:
            rows.append(dict(Outcome=outcome, NetworkPair="VIS–SMN",
                             ObservedConsensusEdges=5, ObservedExpectedRatio=2.,
                             CountPermutationP=.001, CountFDRQ=.01))
        pairs = pd.DataFrame(rows)
        summary = pd.DataFrame(dict(Outcome=figures.OUTCOMES,
                                    TotalMaskPermutationP=[.1,.2,.3],
                                    TotalMaskFDRQ=[.3,.3,.3]))
        null = pd.DataFrame({"Outcome": np.repeat(figures.OUTCOMES, 7)})
        with tempfile.TemporaryDirectory() as name:
            output = Path(name)
            figures.write_captions(pairs, summary, output, "linear", null)
            caption = (output / "figure_captions.md").read_text()
            self.assertIn("7 Freedman-Lane", caption)
            self.assertIn("3 outcome/network-pair tests survived", caption)
            self.assertNotIn("No network pair survived", caption)
            self.assertNotIn("10,000", caption)


if __name__ == "__main__":
    unittest.main()
