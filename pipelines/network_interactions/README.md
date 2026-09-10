# Network interactions — main figures and Table S15

`network_interaction_freedman_lane_main.m` rebuilds the GDE, MEQ30, and VRS
positive consensus masks and tests their counts across 36 network pairs using
10,000 Freedman–Lane permutations. `build_network_interaction_workbook.m` and
the verification scripts operate on its saved outputs.

Table S15 uses two separate multiple-comparison families: BH correction across
the 36 network-pair tests within each outcome, and BH correction across the
three whole-mask tests (GDE, MEQ30 and VRS) in Panel A.

The MATLAB entry point automatically builds and verifies the workbook, renders
the three paper figures, and runs the independent Python audit. A
run is marked `COMPLETE` only after all of these steps pass. The canonical
figures are written to `OUTPUT/figures/`; optional standalone and diagnostic
renders are written to `OUTPUT/figures/other/`.

```matlab
network_interaction_freedman_lane_main(cfg, 10000, ...
    fullfile(cfg.output_root, 'network_interactions'), 5)
```

To redraw saved results without rerunning permutations:

```bash
.venv/bin/python pipelines/network_interactions/render_network_interaction_figures.py \
  --output-dir /path/to/outputs/network_interactions
.venv/bin/python pipelines/network_interactions/verify_network_interaction_outputs.py \
  --output-dir /path/to/outputs/network_interactions \
  --expected-permutations 10000 --require-workbook
```

The figures use the shared linear 0–5.5 white-to-red scale, count each
undirected edge once, and display observed/expected concentration without
significance markers.

The default renderer writes only the three paper figures. Add `--diagnostics`
to both rendering and verification commands to include supporting figures.
The diagnostic log2 option remains available; manuscript panels use the shared
linear 0–5.5 scale.

Here the permutation count excludes the observed ordering: 10,000 means
10,000 null draws plus the observed fit. MATLAB workbook and statistics
helpers default to `<repository>/outputs/network_interactions`; pass the
actual output directory explicitly when using a custom configuration.
