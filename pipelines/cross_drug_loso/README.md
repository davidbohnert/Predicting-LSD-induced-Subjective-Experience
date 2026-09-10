# Cross-drug subject-wise LOSO

`cross_drug_loso_main.m` removes a participant's LSD observation before every
supervised step and then predicts both of that participant's other-drug
sessions. Permutations move the LSD and other-drug outcomes as one profile
within cohort. The entry point accepts outcome, threshold, connectome-form,
parcellation-directory, iteration, worker, and output overrides.

The row-index mapping in `metadata/cross_drug_row_map.csv` contains no subject
identifiers. Correlations over pooled sessions are evaluated by permuting whole
participant profiles rather than treating sessions as independent.
The complete result CSV retains raw permutation p-values; the separate
`table_2_cross_drug.csv` applies Benjamini–Hochberg correction across GDE,
MEQ30, and VRS separately within each pharmacological class.

```matlab
cross_drug_loso_main(cfg); % primary .01 and sensitivity .05/.005 thresholds
```

Individual-drug files suffice; pooled files, when supplied, are checked for
consistent order. Completed cells resume only with compatible settings and inputs.

Table S10 reports threshold sensitivity. The iteration count includes the
observed ordering: 10,000 iterations means 9,999 subject-profile permutations
plus the observed fit.
