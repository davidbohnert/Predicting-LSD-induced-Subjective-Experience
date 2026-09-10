# Aggregate release outputs

Assemble complete analysis outputs after their results have been finalized:

```bash
.venv/bin/python scripts/release/assemble_release_assets.py \
  --analysis-root /path/to/outputs \
  --destination release_assets/assembled
```

The destination must be absent or empty. The assembler checks completion records,
includes aggregate numerical tables, figures, masks and atlas metadata, and
writes a SHA-256 manifest. It excludes participant-level predictions, behaviors,
connectomes, checkpoints, and logs.

Supplementary CSVs use analysis names; the paper correspondence is listed in the
root README. MSE for Table S6 is included in the primary internal source table;
threshold-specific cross-drug results for S10 are included in the cross-drug
source table. No participant data or full analysis rerun is needed to assemble
already completed outputs.

Run the command from the repository root after setting up Python as described
in the main README. Parcellation tables are exported as S13, network-interaction
tables as S15, and network overlap as S16. The assembler accepts legacy
`table_S12_internal.csv` / `table_S12_cross_drug.csv` parcellation sources when
S13 filenames are absent. If both versions exist with different contents,
assembly stops so the intended final source can be selected explicitly.
