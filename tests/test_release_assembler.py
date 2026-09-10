"""Safety checks for the explicit release-asset allowlist."""
from __future__ import annotations

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts" / "release"))

from assemble_release_assets import (  # noqa: E402
    ASSETS,
    REPOSITORY_ASSETS,
    completed_record,
    resolve_source,
    write_manifest,
)


class ReleaseAssemblerTests(unittest.TestCase):
    def test_allowlist_excludes_participant_and_checkpoint_outputs(self) -> None:
        forbidden = ("prediction", "checkpoint", "audit", "permutation_array")
        paths = "\n".join((*ASSETS.keys(), *ASSETS.values())).lower()
        for token in forbidden:
            self.assertNotIn(token, paths)

    def test_current_and_legacy_parcellation_sources(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            root = Path(directory_name)
            current = root / "parcellation_sensitivity/table_S13_internal.csv"
            legacy = current.with_name("table_S12_internal.csv")
            legacy.parent.mkdir()
            legacy.write_bytes(b"Outcome,R\nGDE,0.288\n")
            relative = current.relative_to(root).as_posix()
            self.assertEqual(resolve_source(root, relative), legacy)
            self.assertEqual(ASSETS[relative], "tables/table_S13_internal.csv")
            current.write_bytes(legacy.read_bytes())
            self.assertEqual(resolve_source(root, relative), current)
            current.write_bytes(b"different result")
            with self.assertRaises(SystemExit):
                resolve_source(root, relative)
        self.assertEqual(ASSETS["network_interactions/outcome_summary.csv"],
                         "tables/table_S15_outcomes.csv")
        self.assertEqual(ASSETS["network_overlap/joint_edge_overlap_results.csv"],
                         "tables/table_S16.csv")

    def test_repository_assets_exist(self) -> None:
        for source in REPOSITORY_ASSETS:
            self.assertTrue((REPOSITORY / source).is_file(), source)

    def test_complete_record_and_checksum_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory_name:
            directory = Path(directory_name)
            record = directory / "run_record_001.txt"
            record.write_text("status=COMPLETE\n")
            payload = directory / "table.csv"
            payload.write_bytes(b"a,b\n1,2\n")

            self.assertTrue(completed_record(directory))
            write_manifest(directory)

            expected = hashlib.sha256(payload.read_bytes()).hexdigest()
            manifest = (directory / "MANIFEST.sha256").read_text()
            self.assertIn(f"{expected}  table.csv", manifest)
            self.assertNotIn("MANIFEST.sha256", manifest)


if __name__ == "__main__":
    unittest.main()
