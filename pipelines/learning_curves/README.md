# Learning curves — Figure S2

`learning_curve_main.m` creates deterministic nested, study-stratified samples
and repeats the primary model 100 times at each sample size. Each fit uses
1,000 iterations (one observed and 999 shuffled). The entry point saves the sample plan, numerical summary, and Figure S2. `render_learning_curve_figure(runRoot)` redraws Figure S2 from the
saved summary without rerunning CPM.

```matlab
runRoot = fullfile(cfg.output_root, 'learning_curves');
learning_curve_main(cfg, output_path=runRoot);
% To redraw later:
render_learning_curve_figure(runRoot);
```

This produces the release paths
`learning_curves/output/primary_learning_curve_summary.csv` and
`learning_curves/output/figure_s2_learning_curves.{png,pdf}`.

Optional supporting curves can be generated with
`summarize_learning_curves(runRoot, true)`.
