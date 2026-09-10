# Parcellation sensitivity

Each atlas run performs internal validation followed by subject-wise cross-drug
LOSO. The 439-node condition also reports consensus edges involving cerebellum.

```matlab
parcellation_sensitivity_main(cfg, "original");
parcellation_sensitivity_main(cfg, "no_cerebellum");
parcellation_sensitivity_main(cfg, "with_cerebellum");
build_table_s13(fullfile(cfg.output_root, 'parcellation_sensitivity'));
```

Run one atlas at a time. The table builder reads the completed numerical outputs
without repeating CPM. Current paper correspondence: Table S13.
The builder writes `table_S13_internal.csv` and `table_S13_cross_drug.csv`.
`build_table_s12` remains a compatibility alias; it writes the new filenames
and leaves previously saved S12 files unchanged.
