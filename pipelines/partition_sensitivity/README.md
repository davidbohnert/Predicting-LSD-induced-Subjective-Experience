# Partition sensitivity and stability — Figure S1 / Table S17

`partition_sensitivity_main.m` runs the primary positive, adjusted difference
model for all six outcomes over CV seeds 123–222. Each partition uses 1,000
iterations by default (one observed fit and 999 shuffled fits). `render_partition_sensitivity_figure.py` accepts the
result CSV (or the archived canonical workbook) and writes the PNG/PDF figure.

```bash
.venv/bin/python pipelines/partition_sensitivity/render_partition_sensitivity_figure.py \
  --input /path/to/outputs/partition_sensitivity/partition_sensitivity_results.csv \
  --output /path/to/outputs/partition_sensitivity/figures
```

`edge_node_stability_main.m` uses the same partitions to record edge-selection
recurrence for GDE, MEQ30, and VRS. It checkpoints one MAT file per outcome.
Afterward, build and verify the flat Table S17 source files with:

```bash
.venv/bin/python pipelines/partition_sensitivity/build_table_s17_workbook.py \
  --output-dir /path/to/outputs/partition_sensitivity/stability
.venv/bin/python pipelines/partition_sensitivity/verify_table_s17.py \
  --output-dir /path/to/outputs/partition_sensitivity/stability \
  --atlas-file /path/to/codebase_revised/metadata/atlas/Schaefer416_8networks.txt
```
