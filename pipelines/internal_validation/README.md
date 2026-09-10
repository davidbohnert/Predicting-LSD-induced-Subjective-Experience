# Internal validation

`internal_validation_main(cfg, preset)` runs the chosen specification.
Outputs are a flat CSV, a MAT file with numerical details, and a text summary.
The primary preset additionally reports bootstrap intervals and BH q-values
across the six outcomes.

| Preset | Paper output | Specification |
|---|---|---|
| `primary` | Table 1; S6 | Adjusted positive difference model; MSE included |
| `full_internal_grid` | S5 | LSD/difference, adjusted/unadjusted, all three feature sets |
| `fold_sensitivity` | S7 | 5-fold, 10-fold, and leave-one-subject-out |
| `threshold_sensitivity` | S8 | Edge p thresholds .005, .01, .05 |
| `alternative_correlation` | S9 | Alternative Pearson/Spearman selection |
| `gsr` | S14 | Difference connectomes with global signal regression |
| `matched_sample` | S11 | GDE, VRS, MEQ30 in LSD rows 20–67 |
| `study_loso` | S12 | Hold out each source study |

```matlab
internal_validation_main(cfg, "full_internal_grid");
```

Each completed outcome/connectome/adjustment/fold specification is checkpointed.
Call with the same settings to resume. Changed numerical settings or inputs
require a fresh output directory. Bootstrap resample count may be changed to
rebuild intervals from saved predictions without repeating CPM.

Study covariates are selected for the represented sample, or separately within
study-held-out training sets. Positive/negative consensus counts and MSE remain
available in the numerical outputs even when not displayed in a paper table.

The iteration count includes one observed ordering; 10,000 iterations means
9,999 shuffled fits plus the observed fit.
