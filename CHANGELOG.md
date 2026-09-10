# Changes during manuscript revision

## v2.0.0 — Revised manuscript version (unreleased)

### Analysis changes

- Designate the adjusted positive LSD–placebo difference model as primary;
  retain other model configurations as supplementary sensitivity analyses.
- Add VRS to internal validation and the selected cross-drug and network analyses.
- Encode study categorically and omit redundant indicators in restricted samples
  and study-held-out training folds.
- Replace consensus-mask cross-drug testing with subject-wise LOSO and
  within-cohort subject-profile permutations.
- Report participant-bootstrap intervals and specified multiple-comparison corrections.
- Count unique undirected edges in network figures, normalize by network-pair
  capacity, and use a shared colour scale. Add Freedman–Lane network-count and
  synchronized mask-overlap tests.
- Add partition sensitivity, learning curves, alternative parcellations,
  matched-sample and study-held-out validation, and edge/node stability summaries.

### File map from the original release

| Earlier entry point | Revised successor |
|---|---|
| `scripts/cross_validation/cross_validation_main.m` | [`pipelines/internal_validation/internal_validation_main.m`](pipelines/internal_validation/internal_validation_main.m) |
| `scripts/external_validation/external_validation_main.m` | [`pipelines/cross_drug_loso/cross_drug_loso_main.m`](pipelines/cross_drug_loso/cross_drug_loso_main.m) |
| `scripts/visualization/visualize_networks.py` | [`pipelines/network_interactions/render_network_interaction_figures.py`](pipelines/network_interactions/render_network_interaction_figures.py) |

These are workflow successors, not drop-in APIs:
see the pipeline READMEs for configuration and input requirements. Network
figures now render saved network-interaction results.

### Computation and organization

- Cache fixed quantities within training folds, reuse edge associations across
  thresholds, and store unique undirected edges where applicable.
- Share the internal-CV engine across primary, partition, and learning-curve runs.
- Provide named entry points, shared configuration, compatible checkpoints,
  and compact run records.
- Align supplementary table labels with the revised supplement (S5–S17).
  Parcellation exports now use S13; `build_table_s12` remains a compatibility
  wrapper, and the optional assembler accepts legacy S12 source filenames.
- Remove command-window clearing from entry functions and correct helper
  output defaults to `outputs/network_interactions`.
- Document MATLAB R2022b as the minimum, Python environment selection, and the
  observed-ordering convention for iteration counts.
- Update joint copyright attribution and clarify original and revised
  software contributions.

Computational changes preserve the corresponding model specification; they do
not imply equivalence between scientifically revised analyses and the earlier
release. Outcome-specific correlation methods, difference-behavior fallback,
parallel permutation testing, and consensus masks already existed publicly.

## Earlier public versions

- Original release: tag `v1.0.0`, commit `3fe0c02f610cd779ce774e85695e2a5165a70563`.
- Public branch before this revision: commit
  `4035d87b8a7637b0d19fc5af2fc8f59f573c49c8` (18 March 2026),
  including fixes and documentation changes made after `v1.0.0`.

Both remain identifiable through Git history. Reference implementations used
by numerical tests identify their provenance in their function headers.
