# Network overlap — Table S16

`joint_native_edge_overlap_main.m` uses synchronized Freedman–Lane residual
permutations to compare native GDE, MEQ30, and VRS masks. It preserves the
outcome-availability blocks and tests the mask-size-standardized overlap score
at the primary `.01` threshold; `.05` is retained as a sensitivity analysis.

`build_joint_overlap_table.py --input RESULTS.csv --output TABLE_DIR` produces
compact CSV and Markdown tables from the numerical output.

The default 10,000 permutations are null draws in addition to the observed
fit. Run Python commands from the repository root, for example:

```bash
.venv/bin/python pipelines/network_overlap/build_joint_overlap_table.py \
  --input /path/to/outputs/network_overlap/joint_edge_overlap_results.csv \
  --output /path/to/tables
```
